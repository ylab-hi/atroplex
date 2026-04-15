/*
 * Tests for the full builder pipeline: builder::build_from_samples,
 * remove_tombstones(), build_counters aggregation, and build_summary::collect /
 * write_summary.
 *
 * These exercise the orchestration layer that sits on top of build_gff /
 * build_bam, which the other test suites bypass by calling build_gff::build
 * directly. Covers both the default tombstone-kept-in-tree path and the
 * --prune-tombstones physical-removal path.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "build_summary.hpp"
#include "builder.hpp"
#include "genomic_feature.hpp"
#include "sample_info.hpp"

namespace fs = std::filesystem;

class BuilderPipelineTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        // Per-test unique directory — avoids any cross-test aliasing when
        // all tests in the binary share the same gtest random_seed().
        // Using the test name guarantees uniqueness across the suite and
        // makes failing fixtures easy to inspect.
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_builder_test_";
        dir_name += info ? info->name() : "unknown";
        tmp_dir = fs::temp_directory_path() / dir_name;

        std::error_code ec;
        fs::remove_all(tmp_dir, ec);             // start clean
        fs::create_directories(tmp_dir, ec);
        ASSERT_FALSE(ec) << "Failed to create tmp_dir " << tmp_dir << ": " << ec.message();
        ASSERT_TRUE(fs::exists(tmp_dir)) << "tmp_dir missing after create: " << tmp_dir;
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
    }

    fs::path tmp_dir;

    // Helper: write `contents` to tmp_dir/`name`, flush+close explicitly,
    // and verify the file exists with nonzero size before returning.
    // Any failure becomes an immediate test assertion (fast-fail), so we
    // never hand a half-written path to the builder.
    fs::path write_gtf(const std::string& name, const std::string& contents) {
        fs::path p = tmp_dir / name;
        {
            std::ofstream out(p, std::ios::binary | std::ios::trunc);
            if (!out.is_open()) {
                ADD_FAILURE() << "Failed to open " << p << " for writing";
                return p;
            }
            out << contents;
            out.flush();
            out.close();
            if (!out) {
                ADD_FAILURE() << "Write/flush/close failed for " << p;
                return p;
            }
        }
        if (!fs::exists(p)) {
            ADD_FAILURE() << "File does not exist after write: " << p;
        } else if (fs::file_size(p) == 0) {
            ADD_FAILURE() << "File is empty after write: " << p;
        }
        return p;
    }

    // Write a GTF containing a 3-exon ISM (G1: ISM_FIRST)
    fs::path write_ism_gtf() {
        return write_gtf("ism.gtf",
            "chr1\tTEST\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\ttranscript\t2000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\texon\t2000\t2300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"1\";\n"
            "chr1\tTEST\texon\t3500\t3800\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"2\";\n"
            "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"3\";\n");
    }

    // Write a GTF containing the 4-exon parent that absorbs the ISM above
    fs::path write_parent_gtf() {
        return write_gtf("parent.gtf",
            "chr1\tTEST\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\ttranscript\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"1\";\n"
            "chr1\tTEST\texon\t2000\t2300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"2\";\n"
            "chr1\tTEST\texon\t3500\t3800\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"3\";\n"
            "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"4\";\n");
    }

    // Walk the live grove (B+ tree leaves) and collect every segment key.
    // Mirrors analysis_report::collect's traversal pattern.
    std::vector<const segment_feature*> walk_live_segments(grove_type& grove) const {
        std::vector<const segment_feature*> out;
        for (auto& [seqid, root] : grove.get_root_nodes()) {
            if (!root) continue;
            auto* node = root;
            while (!node->get_is_leaf()) {
                auto& children = node->get_children();
                if (children.empty()) break;
                node = children[0];
            }
            while (node) {
                for (auto* key : node->get_keys()) {
                    auto& feature = key->get_data();
                    if (is_segment(feature)) {
                        out.push_back(&get_segment(feature));
                    }
                }
                node = node->get_next();
            }
        }
        return out;
    }
};

// ── Tier 1.1: Full pipeline counter population ─────────────────────────────
//
// Build cross-file: ISM arrives first (as sample), 4-exon parent arrives
// second (as sample). Reverse absorption tombstones the ISM; the sweep
// physically removes it. Verify every counter surfaces the right value.
TEST_F(BuilderPipelineTest, BuilderFullPipeline_CountersPopulated) {
    auto ism_path = write_ism_gtf();
    auto parent_path = write_parent_gtf();

    std::vector<sample_info> samples;
    samples.emplace_back("ism_sample", ism_path);
    samples.back().type = "sample";
    samples.emplace_back("parent_sample", parent_path);
    samples.back().type = "sample";

    grove_type grove(3);
    auto summary = builder::build_from_samples(
        grove, samples,
        /*threads=*/1, /*min_expression=*/-1.0f, /*absorb=*/true,
        /*min_replicates=*/0, /*fuzzy_tolerance=*/5);

    // 2 transcripts read from the two files
    EXPECT_EQ(summary.counters.input_transcripts, 2u);

    // Neither was forward-merged: ISM created a segment (no parent existed yet),
    // parent created a new segment (different structure) and reverse-absorbed
    // the ISM. So merged_transcripts stays 0.
    EXPECT_EQ(summary.counters.merged_transcripts, 0u);

    // Exactly one tombstoned segment (the ISM), removed by the sweep
    EXPECT_EQ(summary.counters.absorbed_segments, 1u);

    // No filter drops
    EXPECT_EQ(summary.counters.discarded_transcripts, 0u);

    // No replicates configured
    EXPECT_EQ(summary.counters.replicates_merged, 0u);

    // After sweep: only the parent segment remains
    EXPECT_EQ(summary.total_segments, 1u);
    EXPECT_EQ(summary.total_genes, 1u);
}

// ── Tier 1.2: default (no --prune-tombstones) counts but keeps in tree ─────
//
// Default behavior: tombstones are counted and pruned from segment_caches
// and gene_indices so downstream stats never see them, but the segment
// keys remain physically in the B+ tree. All grove consumers defensively
// filter `seg.absorbed`, so this is correctness-safe — the tradeoff is
// a bit of memory / .ggx size. Opt into the slow path with
// --prune-tombstones when producing a distributable index.
TEST_F(BuilderPipelineTest, RemoveTombstones_DefaultKeepsInTree) {
    auto ism_path = write_ism_gtf();
    auto parent_path = write_parent_gtf();

    std::vector<sample_info> samples;
    samples.emplace_back("ism_sample", ism_path);
    samples.back().type = "sample";
    samples.emplace_back("parent_sample", parent_path);
    samples.back().type = "sample";

    grove_type grove(3);
    // prune_tombstones = false (default)
    auto summary = builder::build_from_samples(
        grove, samples,
        /*threads=*/1, /*min_expression=*/-1.0f, /*absorb=*/true,
        /*min_replicates=*/0, /*fuzzy_tolerance=*/5,
        /*prune_tombstones=*/false);

    ASSERT_EQ(summary.counters.absorbed_segments, 1u)
        << "Fixture must produce exactly one tombstone";
    EXPECT_EQ(summary.total_segments, 1u)
        << "total_segments should subtract tombstones";

    auto live = walk_live_segments(grove);
    EXPECT_EQ(live.size(), 2u)
        << "Default: both segments still present in B+ tree "
           "(opt into --prune-tombstones for physical removal)";

    size_t absorbed_in_tree = 0;
    size_t live_in_tree = 0;
    for (const auto* seg : live) {
        if (seg->absorbed) ++absorbed_in_tree;
        else               ++live_in_tree;
    }
    EXPECT_EQ(absorbed_in_tree, 1u);
    EXPECT_EQ(live_in_tree, 1u);
}

// ── Tier 1.2b: --prune-tombstones physically removes the key ───────────────
//
// When the flag is set, remove_tombstones calls grove.remove_key() per
// tombstone and runs a remove_edges_if pass to drop orphan EXON_TO_EXON
// chain edges. The tombstoned segment should no longer surface in a tree
// walk, and the edge count should drop.
TEST_F(BuilderPipelineTest, RemoveTombstones_PruneFlagPhysicallyRemoves) {
    auto ism_path = write_ism_gtf();
    auto parent_path = write_parent_gtf();

    std::vector<sample_info> samples;
    samples.emplace_back("ism_sample", ism_path);
    samples.back().type = "sample";
    samples.emplace_back("parent_sample", parent_path);
    samples.back().type = "sample";

    grove_type grove(3);
    auto summary = builder::build_from_samples(
        grove, samples,
        /*threads=*/1, /*min_expression=*/-1.0f, /*absorb=*/true,
        /*min_replicates=*/0, /*fuzzy_tolerance=*/5,
        /*prune_tombstones=*/true);

    ASSERT_EQ(summary.counters.absorbed_segments, 1u);

    auto live = walk_live_segments(grove);
    EXPECT_EQ(live.size(), 1u)
        << "Expected exactly one live segment (parent) after pruning";
    for (const auto* seg : live) {
        EXPECT_FALSE(seg->absorbed)
            << "Tree traversal should not surface any absorbed segment "
               "when --prune-tombstones is set";
    }
}

// ── Tier 1.3: --prune-tombstones drops orphan EXON_TO_EXON edges ───────────
//
// Under --prune-tombstones, remove_tombstones calls remove_edges_if to
// strip chain edges whose metadata.id matches a tombstoned
// segment_index. Verified by comparing the edge count of a pruned grove
// against a non-pruned baseline built from the same files.
TEST_F(BuilderPipelineTest, RemoveTombstones_PruneFlagDropsOrphanEdges) {
    auto ism_path = write_ism_gtf();
    auto parent_path = write_parent_gtf();

    // Build A: --prune-tombstones ON
    size_t pruned_edge_count = 0;
    {
        std::vector<sample_info> samples;
        samples.emplace_back("ism_sample", ism_path);
        samples.back().type = "sample";
        samples.emplace_back("parent_sample", parent_path);
        samples.back().type = "sample";

        grove_type grove(3);
        auto summary = builder::build_from_samples(
            grove, samples,
            /*threads=*/1, /*min_expression=*/-1.0f, /*absorb=*/true,
            /*min_replicates=*/0, /*fuzzy_tolerance=*/5,
            /*prune_tombstones=*/true);
        ASSERT_EQ(summary.counters.absorbed_segments, 1u);
        pruned_edge_count = grove.edge_count();
    }

    // Build B: --prune-tombstones OFF (default) — orphan edges remain
    size_t default_edge_count = 0;
    {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        std::vector<sample_info> samples;
        samples.emplace_back("ism_sample", ism_path);
        samples.back().type = "sample";
        samples.emplace_back("parent_sample", parent_path);
        samples.back().type = "sample";

        grove_type grove(3);
        auto summary = builder::build_from_samples(
            grove, samples,
            /*threads=*/1, /*min_expression=*/-1.0f, /*absorb=*/true,
            /*min_replicates=*/0, /*fuzzy_tolerance=*/5,
            /*prune_tombstones=*/false);
        ASSERT_EQ(summary.counters.absorbed_segments, 1u);
        default_edge_count = grove.edge_count();
    }

    EXPECT_LT(pruned_edge_count, default_edge_count)
        << "--prune-tombstones should have removed at least one orphan edge. "
        << "Pruned: " << pruned_edge_count
        << ", default: " << default_edge_count;
}

// ── Tier 1.4: build_summary text file contains the processing section ─────
//
// End-to-end: run the full pipeline, call write_summary, re-open the file
// and assert the "Transcript processing" section is present with the
// counter values.
TEST_F(BuilderPipelineTest, BuildSummary_WrittenFileContainsCounters) {
    auto ism_path = write_ism_gtf();
    auto parent_path = write_parent_gtf();

    std::vector<sample_info> samples;
    samples.emplace_back("ism_sample", ism_path);
    samples.back().type = "sample";
    samples.emplace_back("parent_sample", parent_path);
    samples.back().type = "sample";

    grove_type grove(3);
    auto summary = builder::build_from_samples(
        grove, samples, 1, -1.0f, true, 0, 5);

    fs::path summary_path = tmp_dir / "test.ggx.summary";
    summary.write_summary(summary_path.string());

    ASSERT_TRUE(fs::exists(summary_path));

    std::ifstream in(summary_path);
    std::stringstream buf;
    buf << in.rdbuf();
    std::string contents = buf.str();

    EXPECT_NE(contents.find("Transcript processing:"), std::string::npos);
    EXPECT_NE(contents.find("Input:"), std::string::npos);
    EXPECT_NE(contents.find("Merged:"), std::string::npos);
    EXPECT_NE(contents.find("Absorbed:"), std::string::npos);
    EXPECT_NE(contents.find("Discarded:"), std::string::npos);

    // Values should appear verbatim next to their labels
    EXPECT_NE(contents.find("Input:            2"), std::string::npos);
    EXPECT_NE(contents.find("Absorbed:         1"), std::string::npos);
}
