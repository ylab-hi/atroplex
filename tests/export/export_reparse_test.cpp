/*
 * Round-trip and schema validation tests for `subcall::export_gtf`
 * (partial #61). Closes the audit's "export_gtf — NOT TESTED AT ALL,
 * no validation that exported GTF can be re-parsed" gap.
 *
 * Builds a small grove from a known multi-exon GTF, writes the .ggx,
 * runs `atroplex export` end-to-end through its cxxopts CLI surface,
 * then verifies the per-sample GTF output two ways:
 *
 *   1. **Schema**: the exported file has gene + transcript + exon
 *      lines, each with the required attributes (gene_id,
 *      transcript_id), and the exon coordinates round-trip the
 *      original input.
 *
 *   2. **Re-parse**: feeding the exported GTF back through
 *      `build_gff::build` produces a grove with the same segment
 *      count as the original — the export is structurally faithful,
 *      not just textually plausible.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cxxopts.hpp>

#include "build_gff.hpp"
#include "build_summary.hpp"
#include "builder.hpp"
#include "genomic_feature.hpp"
#include "grove_walk.hpp"
#include "sample_info.hpp"
#include "subcall/export_gtf.hpp"

namespace fs = std::filesystem;

class ExportReparseTest : public ::testing::Test {
protected:
    void SetUp() override {
        reset_registries();

        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_export_reparse_";
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

    void reset_registries() {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();
    }

    /// Write a single-sample input GTF carrying one 3-exon transcript.
    /// Coords match what the export round-trip should preserve.
    fs::path write_input_gtf(const std::string& name) {
        std::ostringstream g;
        g << "chr1\tTEST\tgene\t1000\t5000\t.\t+\t.\t"
          << "gene_id \"G_RT\"; gene_name \"GeneRT\"; gene_biotype \"protein_coding\";\n"
          << "chr1\tTEST\ttranscript\t1000\t5000\t.\t+\t.\t"
          << "gene_id \"G_RT\"; transcript_id \"TX_RT\";\n"
          << "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\t"
          << "gene_id \"G_RT\"; transcript_id \"TX_RT\"; exon_number \"1\";\n"
          << "chr1\tTEST\texon\t2500\t2700\t.\t+\t.\t"
          << "gene_id \"G_RT\"; transcript_id \"TX_RT\"; exon_number \"2\";\n"
          << "chr1\tTEST\texon\t4500\t5000\t.\t+\t.\t"
          << "gene_id \"G_RT\"; transcript_id \"TX_RT\"; exon_number \"3\";\n";

        fs::path p = tmp_dir / name;
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        out << g.str();
        out.close();
        return p;
    }

    /// Build a one-sample grove from the input GTF and serialize it
    /// to `index_dir/test.ggx`. Mirrors the save_grove pattern used
    /// elsewhere in the test suite.
    void build_index(const fs::path& index_dir,
                     const fs::path& input_gtf) {
        fs::create_directories(index_dir);

        std::vector<sample_info> samples;
        samples.emplace_back("rt_sample", input_gtf);
        samples.back().type = "sample";

        grove_type grove(3);
        [[maybe_unused]] auto summary = builder::build_from_samples(grove, samples);

        fs::path ggx_path = index_dir / "test.ggx";
        std::ofstream ofs(ggx_path, std::ios::binary);
        ASSERT_TRUE(ofs.is_open());
        ofs.write("AGRX", 4);
        uint16_t version = 1;
        ofs.write(reinterpret_cast<const char*>(&version), sizeof(version));
        gene_registry::instance().serialize(ofs);
        source_registry::instance().serialize(ofs);
        transcript_registry::instance().serialize(ofs);
        sample_registry::instance().serialize(ofs);
        grove.serialize(ofs);
        ofs.close();
    }

    /// Drive subcall::export_gtf through its cxxopts CLI surface
    /// end-to-end (matches the existing pattern in
    /// export_conserved_fraction_test).
    void run_export(const std::vector<std::string>& cli_args) {
        subcall::export_gtf cmd;
        cxxopts::Options options = cmd.parse_args(0, nullptr);

        std::vector<std::string> storage = cli_args;
        std::vector<char*> argv;
        argv.reserve(storage.size());
        for (auto& s : storage) argv.push_back(s.data());

        auto parsed = options.parse(static_cast<int>(argv.size()), argv.data());
        cmd.run(parsed);
    }

    /// Return the per-sample GTF path the export *should* have written
    /// at `<output_dir>/<sample_id>.gtf`. Strict on the filename — no
    /// fallback to "any .gtf in the dir" so a future regression in the
    /// writer's filename convention surfaces as a clear test failure
    /// rather than silently passing against the wrong file.
    fs::path expected_exported_gtf(const fs::path& output_dir,
                                    const std::string& sample_id) const {
        return output_dir / (sample_id + ".gtf");
    }

    /// Slurp a file's lines, dropping blanks and `#` comment lines.
    static std::vector<std::string> non_comment_lines(const fs::path& path) {
        std::vector<std::string> out;
        std::ifstream in(path);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            out.push_back(line);
        }
        return out;
    }

    /// Count GTF lines whose feature-type column matches `kind`.
    static size_t count_records(const std::vector<std::string>& lines,
                                const std::string& kind) {
        size_t n = 0;
        for (const auto& line : lines) {
            std::istringstream is(line);
            std::string seqid, source, feature;
            if ((is >> seqid >> source >> feature) && feature == kind) ++n;
        }
        return n;
    }

    /// Run the original sample GTF + the exported GTF through
    /// `build_gff::build` independently and compare segment counts.
    /// Returns {original_segments, exported_segments}.
    std::pair<size_t, size_t> compare_segment_counts(
        const fs::path& original_gtf,
        const fs::path& exported_gtf
    ) {
        size_t orig = 0;
        {
            reset_registries();
            sample_info info("orig");
            info.type = "sample";
            uint32_t sid = sample_registry::instance().register_data(info);
            grove_type grove(3);
            chromosome_exon_caches ec;
            chromosome_segment_caches sc;
            size_t segment_count = 0;
            build_counters counters;
            build_options opts;
            opts.include_scaffolds = true;
            build_gff::build(grove, original_gtf, sid, ec, sc,
                             segment_count, opts, counters);
            orig = atroplex::test::collect_live_segments(grove).size();
        }

        size_t exp = 0;
        {
            reset_registries();
            sample_info info("exported");
            info.type = "sample";
            uint32_t sid = sample_registry::instance().register_data(info);
            grove_type grove(3);
            chromosome_exon_caches ec;
            chromosome_segment_caches sc;
            size_t segment_count = 0;
            build_counters counters;
            build_options opts;
            opts.include_scaffolds = true;
            build_gff::build(grove, exported_gtf, sid, ec, sc,
                             segment_count, opts, counters);
            exp = atroplex::test::collect_live_segments(grove).size();
        }
        return {orig, exp};
    }
};

// ── Schema: required feature lines and attributes are present ───────
//
// The export must produce gene + transcript + exon lines. Each must
// carry `gene_id`, transcript and exon must additionally carry
// `transcript_id`. This is the minimum any downstream GTF-aware tool
// (StringTie, Salmon, IsoQuant) expects.
TEST_F(ExportReparseTest, Schema_RequiredAttributesPresent) {
    fs::path input = write_input_gtf("input.gtf");
    fs::path index_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    build_index(index_dir, input);
    reset_registries();

    run_export({
        "atroplex export",
        "-g", index_dir.string(),
        "-o", output_dir.string(),
    });

    fs::path exported = expected_exported_gtf(output_dir, "rt_sample");
    ASSERT_TRUE(fs::exists(exported))
        << "Export must produce " << exported
        << " (the expected per-sample filename).";
    ASSERT_GT(fs::file_size(exported), 0u) << "Exported GTF must not be empty";

    auto lines = non_comment_lines(exported);
    EXPECT_EQ(count_records(lines, "gene"), 1u)
        << "Expected exactly one `gene` line.";
    EXPECT_EQ(count_records(lines, "transcript"), 1u)
        << "Expected exactly one `transcript` line.";
    EXPECT_EQ(count_records(lines, "exon"), 3u)
        << "Expected three `exon` lines (input fixture has a 3-exon transcript).";

    // Every gene/transcript/exon line carries gene_id; transcript and
    // exon additionally carry transcript_id. The transcript line's
    // start/end columns must match the input range (chr1:1000-5000)
    // — catches a writer regression that derives transcript bounds
    // from individual exons rather than the segment span.
    bool saw_transcript = false;
    for (const auto& line : lines) {
        std::istringstream is(line);
        std::string seqid, source, feature;
        size_t start = 0, end = 0;
        is >> seqid >> source >> feature >> start >> end;
        EXPECT_NE(line.find("gene_id \"G_RT\""), std::string::npos)
            << "Every record must carry gene_id \"G_RT\". Offending line: " << line;
        if (feature == "transcript" || feature == "exon") {
            EXPECT_NE(line.find("transcript_id \"TX_RT\""), std::string::npos)
                << feature << " line must carry transcript_id \"TX_RT\". Line: " << line;
        }
        if (feature == "transcript") {
            saw_transcript = true;
            EXPECT_EQ(start, 1000u)
                << "Transcript start must equal the input transcript range start.";
            EXPECT_EQ(end, 5000u)
                << "Transcript end must equal the input transcript range end.";
        }
    }
    EXPECT_TRUE(saw_transcript)
        << "Schema check expected to encounter at least one transcript line.";
}

// ── Schema: exon coordinates round-trip the original input ──────────
//
// The three input exons are at chr1:1000-1200, 2500-2700, 4500-5000.
// These coordinates must appear verbatim in the exported GTF. This
// catches off-by-one regressions in the writer (e.g., switching to
// 0-based half-open) and confirms that segment-stored coords resolve
// back to GTF coords correctly.
TEST_F(ExportReparseTest, Schema_ExonCoordinatesRoundTrip) {
    fs::path input = write_input_gtf("input.gtf");
    fs::path index_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    build_index(index_dir, input);
    reset_registries();

    run_export({
        "atroplex export",
        "-g", index_dir.string(),
        "-o", output_dir.string(),
    });

    fs::path exported = expected_exported_gtf(output_dir, "rt_sample");
    ASSERT_TRUE(fs::exists(exported))
        << "Export must produce " << exported;

    auto lines = non_comment_lines(exported);
    std::vector<std::pair<size_t, size_t>> exon_coords;
    for (const auto& line : lines) {
        std::istringstream is(line);
        std::string seqid, source, feature;
        size_t start, end;
        if (!(is >> seqid >> source >> feature >> start >> end)) continue;
        if (feature != "exon") continue;
        exon_coords.push_back({start, end});
    }
    ASSERT_EQ(exon_coords.size(), 3u)
        << "Expected three exon lines to extract coordinates from.";

    // Sort by start so we don't depend on emission order.
    std::sort(exon_coords.begin(), exon_coords.end());
    EXPECT_EQ(exon_coords[0], (std::pair<size_t, size_t>{1000, 1200}));
    EXPECT_EQ(exon_coords[1], (std::pair<size_t, size_t>{2500, 2700}));
    EXPECT_EQ(exon_coords[2], (std::pair<size_t, size_t>{4500, 5000}));
}

// ── Round-trip: re-parsed GTF produces the same segment count ───────
//
// The export must be structurally faithful. Feeding the exported GTF
// back through `build_gff::build` should produce a grove with the same
// number of live segments as the original input — one in this fixture.
// If a future writer change loses an exon line or splits a transcript
// across multiple `transcript_id` values, this assertion catches it.
TEST_F(ExportReparseTest, RoundTrip_BuildFromExportedMatchesOriginal) {
    fs::path input = write_input_gtf("input.gtf");
    fs::path index_dir = tmp_dir / "in";
    fs::path output_dir = tmp_dir / "out";
    build_index(index_dir, input);
    reset_registries();

    run_export({
        "atroplex export",
        "-g", index_dir.string(),
        "-o", output_dir.string(),
    });

    fs::path exported = expected_exported_gtf(output_dir, "rt_sample");
    ASSERT_TRUE(fs::exists(exported))
        << "Export must produce " << exported;

    auto [orig_segments, exp_segments] = compare_segment_counts(input, exported);
    EXPECT_EQ(orig_segments, 1u)
        << "Sanity: the input fixture must produce exactly one segment.";
    EXPECT_EQ(exp_segments, orig_segments)
        << "Re-parsing the exported GTF must produce the same segment count "
           "as the original input. Got " << exp_segments
           << " from the export, " << orig_segments << " from the input.";
}
