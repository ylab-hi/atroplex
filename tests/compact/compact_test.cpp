/*
 * Tests for the `compact` subcommand.
 *
 * Existing builder_pipeline_test.cpp already covers `remove_tombstones`
 * with `physical=true` against an in-memory grove (RemoveTombstones_*
 * cases). These tests exercise the full subcall flow on disk: build a
 * grove with absorbed segments → save to an input directory → invoke
 * subcall::compact end-to-end → verify the output directory contains a
 * compacted .ggx alongside the copied .qtx and .ggx.summary, and that
 * the compacted .ggx loads with only live segments.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <cxxopts.hpp>

#include "build_summary.hpp"
#include "builder.hpp"
#include "genomic_feature.hpp"
#include "grove_walk.hpp"
#include "sample_info.hpp"
#include "subcall/compact.hpp"

namespace fs = std::filesystem;

class CompactSubcallTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_compact_test_";
        dir_name += info ? info->name() : "unknown";
        tmp_dir = fs::temp_directory_path() / dir_name;

        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
        fs::create_directories(tmp_dir, ec);
        ASSERT_FALSE(ec) << "Failed to create tmp_dir " << tmp_dir
                         << ": " << ec.message();
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
    }

    fs::path tmp_dir;

    fs::path write_gtf(const std::string& name, const std::string& contents) {
        fs::path p = tmp_dir / name;
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        out << contents;
        out.close();
        EXPECT_TRUE(fs::exists(p)) << "Failed to write fixture " << p;
        return p;
    }

    // 3-exon ISM that gets reverse-absorbed by the 4-exon parent below.
    fs::path write_ism_gtf() {
        return write_gtf("ism.gtf",
            "chr1\tTEST\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\ttranscript\t2000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\texon\t2000\t2300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"1\";\n"
            "chr1\tTEST\texon\t3500\t3800\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"2\";\n"
            "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ISM_FIRST\"; exon_number \"3\";\n");
    }

    fs::path write_parent_gtf() {
        return write_gtf("parent.gtf",
            "chr1\tTEST\tgene\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\ttranscript\t1000\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; gene_name \"TestGene\"; gene_biotype \"protein_coding\";\n"
            "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"1\";\n"
            "chr1\tTEST\texon\t2000\t2300\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"2\";\n"
            "chr1\tTEST\texon\t3500\t3800\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"3\";\n"
            "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"PARENT_LATER\"; exon_number \"4\";\n");
    }

    /// Build a grove with one tombstoned segment, serialize to `index_dir/<stem>.ggx`,
    /// and emit a matching .qtx + .ggx.summary so it looks like a real index dir.
    /// Returns the stem used for the files.
    std::string build_index_with_tombstone(const fs::path& index_dir,
                                            const std::string& stem) {
        fs::create_directories(index_dir);

        auto ism_path = write_ism_gtf();
        auto parent_path = write_parent_gtf();

        std::vector<sample_info> samples;
        samples.emplace_back("ism_sample", ism_path);
        samples.back().type = "sample";
        samples.emplace_back("parent_sample", parent_path);
        samples.back().type = "sample";

        grove_type grove(3);
        build_options opts;
        // Default: prune_tombstones = false, so the .ggx contains the
        // tombstoned segment for compact to reclaim later.
        opts.qtx_path = (index_dir / (stem + ".qtx")).string();
        auto summary = builder::build_from_samples(grove, samples, opts);

        // ASSERT (not EXPECT) — a fixture without a tombstone would let
        // every test that uses this helper proceed against a malformed
        // grove and surface a confusing downstream failure.
        [&] {
            ASSERT_EQ(summary.counters.absorbed_segments, 1u)
                << "Fixture must produce exactly one tombstone";
        }();

        // Mirror subcall::save_grove
        fs::path ggx_path = index_dir / (stem + ".ggx");
        std::ofstream ofs(ggx_path, std::ios::binary);
        EXPECT_TRUE(ofs.is_open());
        ofs.write("AGRX", 4);
        uint16_t version = 1;
        ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));
        gene_registry::instance().serialize(ofs);
        source_registry::instance().serialize(ofs);
        transcript_registry::instance().serialize(ofs);
        sample_registry::instance().serialize(ofs);
        grove.serialize(ofs);
        ofs.close();

        // Drop a build summary file alongside so we can verify it's
        // copied through to the compacted output.
        fs::path summary_path = index_dir / (stem + ".ggx.summary");
        std::ofstream sum_out(summary_path);
        sum_out << "stub-summary-for-compact-test\n";
        sum_out.close();

        return stem;
    }

    /// Collect every segment key in the grove (live + absorbed). Callers
    /// inspect `seg.absorbed` themselves.
    std::vector<const segment_feature*> walk_segments(grove_type& grove) const {
        std::vector<const segment_feature*> out;
        atroplex::test::for_each_segment(grove,
            [&](const segment_feature& seg, key_ptr) { out.push_back(&seg); });
        return out;
    }

    /// Drive subcall::compact through its public CLI surface so the
    /// validate/setup_grove/execute orchestration is actually exercised.
    void run_compact(const std::vector<std::string>& cli_args) {
        subcall::compact cmd;
        cxxopts::Options options = cmd.parse_args(0, nullptr);

        std::vector<std::string> storage = cli_args;
        std::vector<char*> argv;
        argv.reserve(storage.size());
        for (auto& s : storage) argv.push_back(s.data());

        auto parsed = options.parse(static_cast<int>(argv.size()), argv.data());
        cmd.run(parsed);
    }
};

// ── End-to-end: input has 1 absorbed + 1 live segment; after compact ──
// the output index has only the live segment, and .qtx + .ggx.summary
// are copied through unchanged.
TEST_F(CompactSubcallTest, EndToEnd_RemovesTombstoneAndCopiesSidecars) {
    fs::path input_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    std::string stem = build_index_with_tombstone(input_dir, "test");

    // Sanity: the input index file exists and is non-trivial in size.
    fs::path in_ggx = input_dir / (stem + ".ggx");
    fs::path in_qtx = input_dir / (stem + ".qtx");
    fs::path in_sum = input_dir / (stem + ".ggx.summary");
    ASSERT_TRUE(fs::exists(in_ggx));
    ASSERT_TRUE(fs::exists(in_qtx))
        << "build_from_samples must have written a .qtx alongside the .ggx";
    ASSERT_TRUE(fs::exists(in_sum));
    auto in_ggx_size = fs::file_size(in_ggx);
    auto in_qtx_bytes = fs::file_size(in_qtx);

    // Reset registries so compact's load goes through deserialize cleanly,
    // not through the in-memory state left by build_from_samples.
    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    run_compact({
        "atroplex compact",
        "-g", input_dir.string(),
        "-o", output_dir.string(),
    });

    fs::path out_ggx = output_dir / (stem + ".ggx");
    fs::path out_qtx = output_dir / (stem + ".qtx");
    fs::path out_sum = output_dir / (stem + ".ggx.summary");

    EXPECT_TRUE(fs::exists(out_ggx)) << "Compacted .ggx must exist at " << out_ggx;
    EXPECT_TRUE(fs::exists(out_qtx)) << ".qtx must be copied through";
    EXPECT_TRUE(fs::exists(out_sum)) << ".ggx.summary must be copied through";

    EXPECT_LT(fs::file_size(out_ggx), in_ggx_size)
        << "Compacted .ggx should be smaller than the input "
           "(at least the tombstone's serialized payload should be gone)";
    EXPECT_EQ(fs::file_size(out_qtx), in_qtx_bytes)
        << ".qtx is copied verbatim — sizes must match exactly";

    // Reload the compacted grove and assert only the parent survives.
    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    std::ifstream ifs(out_ggx, std::ios::binary);
    ASSERT_TRUE(ifs.is_open());
    char magic[4]; ifs.read(magic, 4);
    ASSERT_EQ(std::string(magic, 4), "AGRX");
    uint16_t version; ifs.read(reinterpret_cast<char*>(&version), sizeof(version));
    ASSERT_EQ(version, 1);
    gene_registry::instance().deserialize_into(ifs);
    source_registry::instance().deserialize_into(ifs);
    transcript_registry::instance().deserialize_into(ifs);
    (void)sample_registry::deserialize(ifs);
    auto loaded = std::make_unique<grove_type>(grove_type::deserialize(ifs));

    auto live = walk_segments(*loaded);
    ASSERT_EQ(live.size(), 1u)
        << "Compacted grove must contain exactly the parent segment";
    EXPECT_FALSE(live[0]->absorbed)
        << "No segment in a compacted grove should still carry the tombstone flag";
}

// ── Validation: missing .qtx without --no-qtx must throw. ─────────────
TEST_F(CompactSubcallTest, Validate_FailsWhenQtxMissingWithoutOptOut) {
    fs::path input_dir = tmp_dir / "in";
    std::string stem = build_index_with_tombstone(input_dir, "test");

    // Delete the .qtx so the strict check trips.
    fs::remove(input_dir / (stem + ".qtx"));

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    EXPECT_THROW(run_compact({
        "atroplex compact",
        "-g", input_dir.string(),
        "-o", (tmp_dir / "out").string(),
    }), std::runtime_error)
        << "compact must reject inputs without a .qtx unless --no-qtx is set";
}

// ── Validation: --no-qtx allows compact to proceed and skips qtx copy. ─
TEST_F(CompactSubcallTest, Validate_NoQtxOptOutSucceedsAndSkipsCopy) {
    fs::path input_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    std::string stem = build_index_with_tombstone(input_dir, "test");

    fs::remove(input_dir / (stem + ".qtx"));

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    run_compact({
        "atroplex compact",
        "-g", input_dir.string(),
        "-o", output_dir.string(),
        "--no-qtx",
    });

    EXPECT_TRUE(fs::exists(output_dir / (stem + ".ggx")));
    EXPECT_FALSE(fs::exists(output_dir / (stem + ".qtx")))
        << "--no-qtx must not produce or copy a .qtx";
}

// ── Validation: --no-qtx wins when a .qtx is present — the warning ───
// fires (not asserted here) and the sidecar is intentionally NOT carried
// into the output, so the compacted index is structure-only.
TEST_F(CompactSubcallTest, Execute_NoQtxWithExistingSidecarSkipsCopy) {
    fs::path input_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    std::string stem = build_index_with_tombstone(input_dir, "test");
    ASSERT_TRUE(fs::exists(input_dir / (stem + ".qtx")))
        << "Pre-condition: input directory must have a .qtx";

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    run_compact({
        "atroplex compact",
        "-g", input_dir.string(),
        "-o", output_dir.string(),
        "--no-qtx",
    });

    EXPECT_TRUE(fs::exists(output_dir / (stem + ".ggx")));
    EXPECT_FALSE(fs::exists(output_dir / (stem + ".qtx")))
        << "--no-qtx must skip the .qtx copy even when the input has one";
    EXPECT_TRUE(fs::exists(input_dir / (stem + ".qtx")))
        << "--no-qtx is not destructive — the input .qtx must remain untouched";
}

// ── Validation: refuses to write the compacted index back into the ───
// input directory. Without this guard, a partially-written .ggx could
// stomp on the original.
TEST_F(CompactSubcallTest, Validate_RefusesSameInputAndOutputDir) {
    fs::path input_dir = tmp_dir / "in";
    [[maybe_unused]] std::string stem = build_index_with_tombstone(input_dir, "test");

    transcript_registry::reset();
    gene_registry::reset();
    source_registry::reset();
    sample_registry::reset();

    EXPECT_THROW(run_compact({
        "atroplex compact",
        "-g", input_dir.string(),
        "-o", input_dir.string(),
    }), std::runtime_error)
        << "compact must refuse to overwrite its input directory";
}
