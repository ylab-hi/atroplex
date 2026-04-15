/*
 * Tests for the full builder pipeline: builder::build_from_samples,
 * remove_tombstones(), build_counters aggregation, and build_summary::collect /
 * write_summary.
 *
 * These exercise the orchestration layer that sits on top of build_gff / build_bam,
 * which the other test suites bypass by calling build_gff::build directly.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "build_gff.hpp"
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

        tmp_dir = fs::temp_directory_path() / ("atroplex_builder_test_" +
            std::to_string(::testing::UnitTest::GetInstance()->random_seed()));
        fs::create_directories(tmp_dir);
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
    }

    fs::path tmp_dir;

    // Write a GTF containing a 3-exon ISM (G1: ISM_FIRST)
    fs::path write_ism_gtf() {
        fs::path p = tmp_dir / "ism.gtf";
        std::ofstream out(p);
        out << "chr1\tTEST\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n";
        out << "chr1\tTEST\ttranscript\t2000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n";
        out << "chr1\tTEST\texon\t2000\t2300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"1\";\n";
        out << "chr1\tTEST\texon\t3500\t3800\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"2\";\n";
        out << "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"3\";\n";
        return p;
    }

    // Write a GTF containing the 4-exon parent that absorbs the ISM above
    fs::path write_parent_gtf() {
        fs::path p = tmp_dir / "parent.gtf";
        std::ofstream out(p);
        out << "chr1\tTEST\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n";
        out << "chr1\tTEST\ttranscript\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n";
        out << "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"1\";\n";
        out << "chr1\tTEST\texon\t2000\t2300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"2\";\n";
        out << "chr1\tTEST\texon\t3500\t3800\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"3\";\n";
        out << "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"4\";\n";
        return p;
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

// ── Tier 1.2: Physical removal of tombstones from the B+ tree ──────────────
//
// After the sweep runs, no segment in live tree traversal should have
// absorbed == true. This confirms remove_tombstones() actually unlinks
// keys from the tree (rather than leaving them as stale entries).
TEST_F(BuilderPipelineTest, RemoveTombstones_PhysicallyRemoved) {
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

    ASSERT_EQ(summary.counters.absorbed_segments, 1u)
        << "Fixture must produce exactly one tombstone";

    auto live = walk_live_segments(grove);
    EXPECT_EQ(live.size(), 1u)
        << "Expected exactly one live segment (the parent) after sweep";

    for (const auto* seg : live) {
        EXPECT_FALSE(seg->absorbed)
            << "Tree traversal should not surface any absorbed segment";
    }
}

// ── Tier 1.3: Orphan EXON_TO_EXON edges are pruned by the sweep ────────────
//
// When reverse absorption tombstones a segment, the old EXON_TO_EXON edges
// carrying that segment's segment_index become orphans (they link two exons
// that are legitimately in the live parent chain, but with a dead id).
// remove_tombstones calls remove_edges_if to strip them. Verify edge_count
// drops compared to a no-sweep baseline built via the same files.
TEST_F(BuilderPipelineTest, RemoveTombstones_OrphanEdgesPruned) {
    auto ism_path = write_ism_gtf();
    auto parent_path = write_parent_gtf();

    // Build A: full builder pipeline (sweep runs)
    size_t swept_edge_count = 0;
    {
        std::vector<sample_info> samples;
        samples.emplace_back("ism_sample", ism_path);
        samples.back().type = "sample";
        samples.emplace_back("parent_sample", parent_path);
        samples.back().type = "sample";

        grove_type grove(3);
        auto summary = builder::build_from_samples(
            grove, samples, 1, -1.0f, true, 0, 5);
        ASSERT_EQ(summary.counters.absorbed_segments, 1u);
        swept_edge_count = grove.edge_count();
    }

    // Build B: bypass the sweep by calling build_gff::build directly
    // twice against the same caches. Tombstone remains in place with its
    // orphan edges.
    size_t unswept_edge_count = 0;
    {
        // Clear registries so sample IDs start at 0 again
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        grove_type grove(3);
        chromosome_exon_caches exon_caches;
        chromosome_segment_caches segment_caches;
        chromosome_gene_segment_indices gene_indices;
        size_t segment_count = 0;
        build_counters counters;

        sample_info ism_info("ism_sample", ism_path);
        ism_info.type = "sample";
        uint32_t ism_id = sample_registry::instance().register_data(ism_info);
        build_gff::build(grove, ism_path, ism_id, exon_caches, segment_caches,
                         gene_indices, segment_count, 0, -1.0f, true, 5, counters);

        sample_info parent_info("parent_sample", parent_path);
        parent_info.type = "sample";
        uint32_t parent_id = sample_registry::instance().register_data(parent_info);
        build_gff::build(grove, parent_path, parent_id, exon_caches, segment_caches,
                         gene_indices, segment_count, 0, -1.0f, true, 5, counters);

        unswept_edge_count = grove.edge_count();
    }

    EXPECT_LT(swept_edge_count, unswept_edge_count)
        << "Sweep should have removed at least one orphan edge. "
        << "Swept: " << swept_edge_count
        << ", unswept: " << unswept_edge_count;
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
