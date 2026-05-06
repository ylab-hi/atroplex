/*
 * Malformed-GTF input behavior tests (partial #61).
 *
 * Pins the current behavior of `build_gff::build` on degenerate or
 * malformed input — empty files, missing required attributes, invalid
 * coordinates, truncated lines, non-existent paths. Most of these
 * silently `continue` past the bad record and produce an empty (or
 * partial) grove without surfacing a warning or error. That's a
 * documented gap, not a desired behavior; if a future PR adds warning
 * counters or strict-mode rejection, the assertions in this file are
 * the regression markers to flip.
 *
 * Each test writes a small inline GTF, runs `build_gff::build`, and
 * asserts on the resulting segment count, build counters, and (for
 * file-open failures) the thrown exception.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "build_gff.hpp"
#include "build_summary.hpp"
#include "genomic_feature.hpp"
#include "grove_walk.hpp"
#include "sample_info.hpp"

namespace fs = std::filesystem;

class MalformedGtfTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_malformed_test_";
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

    /// Write a GTF with the given inline contents and return its path.
    fs::path write_gtf(const std::string& name, const std::string& contents) {
        fs::path p = tmp_dir / name;
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        out << contents;
        out.close();
        EXPECT_TRUE(fs::exists(p)) << "Failed to write " << p;
        return p;
    }

    /// Drive `build_gff::build` and report (live segment count, counters).
    /// Caller assertions check both directly. Mirrors `AbsorptionTest`'s
    /// build_fixture pattern but without the BuildResult struct because
    /// these tests rarely care about transcript counts.
    struct build_outcome {
        size_t segment_count;     ///< Running counter after build (== total inserts)
        size_t live_segments;     ///< Live segments in the live grove
        build_counters counters;  ///< input/merged/absorbed/discarded counters
    };

    build_outcome build(const fs::path& fixture,
                       bool as_annotation = false) {
        sample_info info("malformed_test_sample");
        if (as_annotation) {
            info.type = "annotation";
            info.annotation_source = "TEST";
        }
        uint32_t sample_id = sample_registry::instance().register_data(info);

        grove_type grove(3);
        chromosome_exon_caches exon_caches;
        chromosome_segment_caches segment_caches;
        size_t segment_count = 0;
        build_counters counters;
        build_options opts;
        opts.include_scaffolds = true;

        build_gff::build(grove, fixture, sample_id, exon_caches, segment_caches,
                         segment_count, opts, counters);

        size_t live = atroplex::test::collect_live_segments(grove).size();
        return {segment_count, live, counters};
    }
};

// ── Empty file ──────────────────────────────────────────────────────
//
// A truly empty file (zero bytes) must not crash; it should produce a
// build with zero segments and zero input transcripts. The genogrove
// gff reader returns from its iterator immediately when there's
// nothing to read.
TEST_F(MalformedGtfTest, EmptyFile_ProducesEmptyGrove) {
    auto path = write_gtf("empty.gtf", "");

    auto out = build(path);
    EXPECT_EQ(out.segment_count, 0u)
        << "Empty file must not produce any segments.";
    EXPECT_EQ(out.live_segments, 0u);
    EXPECT_EQ(out.counters.input_transcripts, 0u);
    EXPECT_EQ(out.counters.merged_transcripts, 0u);
    EXPECT_EQ(out.counters.discarded_transcripts, 0u);
}

// ── Comments / blank lines only ─────────────────────────────────────
//
// A file containing only header comments (`#`) and blank lines is
// equivalent to empty for build purposes — the gff reader skips
// comment and blank lines, so no records reach `build_gff`.
TEST_F(MalformedGtfTest, CommentsAndBlankLinesOnly_ProducesEmptyGrove) {
    auto path = write_gtf("comments.gtf",
        "##gff-version 3\n"
        "# this is a comment\n"
        "\n"
        "## another comment\n"
        "\n");

    auto out = build(path);
    EXPECT_EQ(out.segment_count, 0u);
    EXPECT_EQ(out.live_segments, 0u);
    EXPECT_EQ(out.counters.input_transcripts, 0u);
}

// ── Non-existent file ───────────────────────────────────────────────
//
// Opening a non-existent path is a runtime error, not a silent skip.
// The gff reader / build path surfaces it as a thrown exception.
TEST_F(MalformedGtfTest, NonExistentFile_Throws) {
    fs::path missing = tmp_dir / "this_file_does_not_exist.gtf";
    ASSERT_FALSE(fs::exists(missing));

    EXPECT_ANY_THROW(build(missing))
        << "Build against a non-existent file must throw, not silently "
           "produce an empty grove.";
}

// ── Exon line missing `gene_id` ─────────────────────────────────────
//
// `build_gff::build` (src/build_gff.cpp:76-79) calls `entry.get_gene_id()`
// for every record and `continue`s past entries without a `gene_id`. This
// silently drops the offending records. The current behavior is
// preserved here so a future PR adding a warning counter or strict-mode
// rejection has a regression marker.
TEST_F(MalformedGtfTest, ExonMissingGeneId_SilentlySkipped) {
    auto path = write_gtf("missing_gene_id.gtf",
        "chr1\tTEST\ttranscript\t1000\t3000\t.\t+\t.\ttranscript_id \"TX_NO_GENE\";\n"
        "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\ttranscript_id \"TX_NO_GENE\"; exon_number \"1\";\n"
        "chr1\tTEST\texon\t2000\t2200\t.\t+\t.\ttranscript_id \"TX_NO_GENE\"; exon_number \"2\";\n"
        "chr1\tTEST\texon\t2800\t3000\t.\t+\t.\ttranscript_id \"TX_NO_GENE\"; exon_number \"3\";\n");

    auto out = build(path);
    EXPECT_EQ(out.segment_count, 0u)
        << "All entries lack gene_id and are silently skipped — no segment "
           "should be produced. (Documented gap: ingest currently surfaces "
           "no warning. See #61.)";
    EXPECT_EQ(out.live_segments, 0u);
    EXPECT_EQ(out.counters.input_transcripts, 0u);
}

// ── Exon line missing `transcript_id` ───────────────────────────────
//
// `process_gene` (src/build_gff.cpp:120-124) only adds entries with a
// `transcript_id` to the transcripts map. Entries without one are
// silently dropped at the gene level; if every entry lacks a
// transcript_id, the gene yields zero transcripts and `process_gene`
// emits no segments. Documents the same gap as the missing-gene-id
// case above.
TEST_F(MalformedGtfTest, ExonMissingTranscriptId_SilentlySkipped) {
    auto path = write_gtf("missing_transcript_id.gtf",
        "chr1\tTEST\tgene\t1000\t3000\t.\t+\t.\tgene_id \"G1\";\n"
        "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; exon_number \"1\";\n"
        "chr1\tTEST\texon\t2000\t2200\t.\t+\t.\tgene_id \"G1\"; exon_number \"2\";\n"
        "chr1\tTEST\texon\t2800\t3000\t.\t+\t.\tgene_id \"G1\"; exon_number \"3\";\n");

    auto out = build(path);
    EXPECT_EQ(out.segment_count, 0u)
        << "Every exon lacks transcript_id, so process_gene's transcripts "
           "map stays empty and no segment is created.";
    EXPECT_EQ(out.live_segments, 0u);
    EXPECT_EQ(out.counters.input_transcripts, 0u);
}

// ── Mixed valid + missing-attribute records ─────────────────────────
//
// A mixed file with one well-formed transcript plus a parallel
// transcript missing `gene_id` should produce exactly the well-formed
// segment. Confirms that the silent-skip path doesn't poison a build
// where some records are valid.
TEST_F(MalformedGtfTest, MixedValidAndMissingGeneId_OnlyValidIngested) {
    auto path = write_gtf("mixed.gtf",
        "chr1\tTEST\tgene\t1000\t3000\t.\t+\t.\tgene_id \"G1\"; gene_name \"GeneOne\"; gene_biotype \"protein_coding\";\n"
        "chr1\tTEST\ttranscript\t1000\t3000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\";\n"
        "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\"; exon_number \"1\";\n"
        "chr1\tTEST\texon\t2000\t2200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\"; exon_number \"2\";\n"
        "chr1\tTEST\texon\t2800\t3000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\"; exon_number \"3\";\n"
        // Parallel transcript missing gene_id — silently dropped per build_gff:76-79.
        "chr1\tTEST\ttranscript\t5000\t7000\t.\t+\t.\ttranscript_id \"TX_NO_GENE\";\n"
        "chr1\tTEST\texon\t5000\t5200\t.\t+\t.\ttranscript_id \"TX_NO_GENE\"; exon_number \"1\";\n"
        "chr1\tTEST\texon\t6000\t6200\t.\t+\t.\ttranscript_id \"TX_NO_GENE\"; exon_number \"2\";\n");

    auto out = build(path);
    EXPECT_EQ(out.segment_count, 1u)
        << "Only the well-formed TX_A should land in the grove; TX_NO_GENE "
           "is silently dropped because its records lack gene_id.";
    EXPECT_EQ(out.live_segments, 1u);
    EXPECT_EQ(out.counters.input_transcripts, 1u)
        << "input_transcripts must reflect the count of records that "
           "actually entered the build pipeline (1, not 2).";
}

// ── Truncated final line ────────────────────────────────────────────
//
// Some GTFs end without a trailing newline or with a partially-written
// last line (truncated copy, interrupted download). Records that did
// finish should still ingest cleanly; the truncated tail is the gff
// reader's problem to handle. Confirms ingest doesn't abort the whole
// build on a partial last line.
TEST_F(MalformedGtfTest, TruncatedFinalLine_PartialRecordIgnored) {
    // Two complete exon records on TX_A, then a third exon line that
    // ends mid-attribute (no closing `;` and no newline).
    auto path = write_gtf("truncated.gtf",
        "chr1\tTEST\tgene\t1000\t3000\t.\t+\t.\tgene_id \"G1\"; gene_name \"GeneOne\"; gene_biotype \"protein_coding\";\n"
        "chr1\tTEST\ttranscript\t1000\t3000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\";\n"
        "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\"; exon_number \"1\";\n"
        "chr1\tTEST\texon\t2000\t2200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\"; exon_number \"2\";\n"
        "chr1\tTEST\texon\t2800\t3000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX_A\"; exon_numb");

    // We deliberately don't pin a specific segment count here — the
    // genogrove reader's exact behavior on a truncated last line is
    // version-dependent (drop, recover, throw). What we DO require is
    // that the partial record never produces a corrupt segment that
    // flows into the live grove. With only 2 well-formed exons, no
    // multi-exon segment can be reverse-absorbed into a 3-exon parent
    // either way, so the live grove must not exceed 1 segment.
    auto out = build(path);
    EXPECT_LE(out.live_segments, 1u)
        << "A truncated final exon must not produce extra live segments "
           "beyond the well-formed prefix.";
}
