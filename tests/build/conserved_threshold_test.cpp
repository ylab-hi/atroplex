/*
 * Tests for the conservation-threshold mechanics:
 *   - is_conserved() semantics (>= min_required, not == total)
 *   - analysis_report::set_conserved_fraction validation + clamping
 *   - End-to-end: per-sample conserved counts and conserved_segments.tsv
 *     reflect the configured fraction.
 *
 * Strategy: build a small grove from multiple sample-typed entries that all
 * declare the same exon chain, so the segment's sample_idx accumulates one
 * bit per sample. Then run analysis_report::collect() under various
 * fraction settings and assert on the resulting per-sample counters and
 * the conserved_segments.tsv file contents.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "analysis_report.hpp"
#include "builder.hpp"
#include "genomic_feature.hpp"
#include "sample_info.hpp"

namespace fs = std::filesystem;

class ConservedThresholdTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_conserved_test_";
        dir_name += info ? info->name() : "unknown";
        tmp_dir = fs::temp_directory_path() / dir_name;

        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
        fs::create_directories(tmp_dir, ec);
        ASSERT_FALSE(ec) << "Failed to create tmp_dir " << tmp_dir << ": " << ec.message();
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
    }

    fs::path tmp_dir;

    // Write a GTF with a single 3-exon transcript. The sample-specific
    // transcript_id ensures every file contributes a registry-distinct
    // tx_id; the shared exon chain causes them all to dedup into one
    // segment whose sample_idx accumulates one bit per file.
    fs::path write_shared_gtf(const std::string& tag) {
        std::ostringstream g;
        g << "chr1\tTEST\ttranscript\t1000\t5000\t.\t+\t.\t"
          << "gene_id \"G1\"; transcript_id \"TX_" << tag << "\";\n"
          << "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\t"
          << "gene_id \"G1\"; transcript_id \"TX_" << tag << "\"; exon_number \"1\";\n"
          << "chr1\tTEST\texon\t2000\t2200\t.\t+\t.\t"
          << "gene_id \"G1\"; transcript_id \"TX_" << tag << "\"; exon_number \"2\";\n"
          << "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\t"
          << "gene_id \"G1\"; transcript_id \"TX_" << tag << "\"; exon_number \"3\";\n";

        fs::path p = tmp_dir / ("shared_" + tag + ".gtf");
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        out << g.str();
        out.close();
        return p;
    }

    // Same skeleton but with a non-overlapping exon set so this sample
    // does NOT contribute to the shared segment. Used for "dropout"
    // scenarios where N-1 of N samples share a structure.
    fs::path write_distinct_gtf(const std::string& tag) {
        std::ostringstream g;
        g << "chr1\tTEST\ttranscript\t20000\t25000\t.\t+\t.\t"
          << "gene_id \"G2_" << tag << "\"; transcript_id \"DTX_" << tag << "\";\n"
          << "chr1\tTEST\texon\t20000\t20200\t.\t+\t.\t"
          << "gene_id \"G2_" << tag << "\"; transcript_id \"DTX_" << tag << "\"; exon_number \"1\";\n"
          << "chr1\tTEST\texon\t22000\t22200\t.\t+\t.\t"
          << "gene_id \"G2_" << tag << "\"; transcript_id \"DTX_" << tag << "\"; exon_number \"2\";\n"
          << "chr1\tTEST\texon\t24500\t25000\t.\t+\t.\t"
          << "gene_id \"G2_" << tag << "\"; transcript_id \"DTX_" << tag << "\"; exon_number \"3\";\n";

        fs::path p = tmp_dir / ("distinct_" + tag + ".gtf");
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        out << g.str();
        out.close();
        return p;
    }

    // Build a grove from N sample-typed entries; each gets `share_count`
    // contributions to the shared structure and `1 - shared_share` chance
    // of getting a distinct structure instead. Returns the grove.
    grove_type build_with(size_t shared_count, size_t distinct_count) {
        std::vector<sample_info> samples;
        samples.reserve(shared_count + distinct_count);
        for (size_t i = 0; i < shared_count; ++i) {
            std::string tag = "S" + std::to_string(i);
            samples.emplace_back("sample_" + tag, write_shared_gtf(tag));
            samples.back().type = "sample";
        }
        for (size_t i = 0; i < distinct_count; ++i) {
            std::string tag = "D" + std::to_string(i);
            samples.emplace_back("dropout_" + tag, write_distinct_gtf(tag));
            samples.back().type = "sample";
        }
        grove_type grove(3);
        builder::build_from_samples(grove, samples);
        return grove;
    }
};

// ── is_conserved semantics ──────────────────────────────────────────

TEST_F(ConservedThresholdTest, IsConserved_GreaterEqualThreshold) {
    // Build a 3-sample grove where all three share the structure.
    auto grove = build_with(/*shared=*/3, /*distinct=*/0);

    // Find the live shared segment.
    const segment_feature* shared_seg = nullptr;
    auto roots = grove.get_root_nodes();
    for (auto& [seqid, root] : roots) {
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
                if (!is_segment(feature)) continue;
                auto& seg = get_segment(feature);
                if (seg.absorbed) continue;
                if (seg.sample_count() == 3 && shared_seg == nullptr) {
                    shared_seg = &seg;
                }
            }
            node = node->get_next();
        }
    }
    ASSERT_NE(shared_seg, nullptr) << "Expected a shared segment with 3 samples";

    // >= threshold: 1, 2, 3 all true; 4 false.
    EXPECT_TRUE (shared_seg->is_conserved(1));
    EXPECT_TRUE (shared_seg->is_conserved(2));
    EXPECT_TRUE (shared_seg->is_conserved(3));
    EXPECT_FALSE(shared_seg->is_conserved(4));
}

// ── set_conserved_fraction validation ───────────────────────────────

TEST_F(ConservedThresholdTest, SetConservedFraction_InvalidFallsBackToOne) {
    // Out-of-range input must not throw and must not classify the
    // 3-of-4 shared structure as conserved (silent fallback to 1.0).
    auto grove = build_with(/*shared=*/3, /*distinct=*/1);

    fs::path tsv = tmp_dir / "invalid_fraction.conserved_segments.tsv";
    {
        analysis_report report;
        report.set_conserved_fraction(-0.5);  // invalid → silently 1.0
        report.begin_conserved_segment_stream(tsv.string(), false);
        ASSERT_NO_THROW(report.collect(grove));
    }

    ASSERT_TRUE(fs::exists(tsv));
    std::ifstream in(tsv);
    std::string line;
    size_t data_rows = 0;
    bool header_seen = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (!header_seen) { header_seen = true; continue; }
        ++data_rows;
    }
    EXPECT_EQ(data_rows, 0u)
        << "Invalid fraction should fall back to strict 1.0, "
        << "leaving the 3-of-4 structure NOT conserved";
}

// ── End-to-end: conserved_segments.tsv reflects the threshold ───────

TEST_F(ConservedThresholdTest, ConservedSegmentsTsv_StrictExcludesDropout) {
    // 3 of 4 samples share a structure → strict 1.0 should NOT classify
    // it as conserved → conserved_segments.tsv should have 0 data rows.
    auto grove = build_with(/*shared=*/3, /*distinct=*/1);

    fs::path tsv = tmp_dir / "strict.conserved_segments.tsv";
    {
        analysis_report report;
        // default fraction is 1.0
        report.begin_conserved_segment_stream(tsv.string(), /*emit_expression=*/false);
        report.collect(grove);
    }

    ASSERT_TRUE(fs::exists(tsv));
    std::ifstream in(tsv);
    std::string line;
    size_t data_rows = 0;
    bool header_seen = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (!header_seen) { header_seen = true; continue; }
        ++data_rows;
    }
    EXPECT_TRUE(header_seen) << "Header missing from conserved_segments.tsv";
    EXPECT_EQ(data_rows, 0u)
        << "At fraction=1.0 with 3-of-4 sharing, no segment should be conserved";
}

TEST_F(ConservedThresholdTest, ConservedSegmentsTsv_RelaxedIncludesDropout) {
    // Same grove; relaxed fraction 0.75 → ceil(4 * 0.75) = 3, the
    // shared segment (count=3) qualifies → 1 data row expected.
    auto grove = build_with(/*shared=*/3, /*distinct=*/1);

    fs::path tsv = tmp_dir / "relaxed.conserved_segments.tsv";
    {
        analysis_report report;
        report.set_conserved_fraction(0.75);
        report.begin_conserved_segment_stream(tsv.string(), /*emit_expression=*/false);
        report.collect(grove);
    }

    ASSERT_TRUE(fs::exists(tsv));
    std::ifstream in(tsv);
    std::string line;
    size_t data_rows = 0;
    bool header_seen = false;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (!header_seen) { header_seen = true; continue; }
        ++data_rows;
    }
    EXPECT_TRUE(header_seen);
    EXPECT_GE(data_rows, 1u)
        << "At fraction=0.75 the 3/4-shared structure should appear as conserved";
}