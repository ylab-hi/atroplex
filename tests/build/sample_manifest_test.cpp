/*
 * Edge-case tests for `sample_manifest` (TSV parser, partial #61).
 *
 * Covers the parser-level invariants documented in
 * include/sample_manifest.hpp: required `file` column, case-insensitive
 * headers, whitespace tolerance, `.` as the empty-value marker, auto-
 * inferred ID and replicate group, expression-attribute splitting.
 *
 * Each test writes a small TSV to a temp dir and parses it. Real GTF
 * files referenced from the manifest must exist on disk (the parser
 * resolves them via `std::filesystem::canonical`); the helpers below
 * touch empty placeholder files alongside the manifest.
 */

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

#include "sample_manifest.hpp"

namespace fs = std::filesystem;

class SampleManifestTest : public ::testing::Test {
protected:
    void SetUp() override {
        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_manifest_test_";
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

    /// Create an empty file at `tmp_dir/name` so the manifest's
    /// `canonical()` resolution succeeds.
    fs::path touch(const std::string& name) {
        fs::path p = tmp_dir / name;
        std::ofstream out(p);
        out.close();
        EXPECT_TRUE(fs::exists(p)) << "Failed to create placeholder " << p;
        return p;
    }

    /// Write `contents` to `tmp_dir/name` and return the path.
    fs::path write_file(const std::string& name, const std::string& contents) {
        fs::path p = tmp_dir / name;
        std::ofstream out(p);
        out << contents;
        out.close();
        EXPECT_TRUE(fs::exists(p)) << "Failed to write " << p;
        return p;
    }
};

// ── Required column ─────────────────────────────────────────────────

TEST_F(SampleManifestTest, MissingFileColumn_Throws) {
    auto path = write_file("manifest.tsv",
        "id\ttype\nfoo\tsample\n");

    EXPECT_THROW(sample_manifest manifest(path), std::runtime_error)
        << "Manifest without a 'file' column must throw at parse time.";
}

TEST_F(SampleManifestTest, EmptyFile_Throws) {
    auto path = write_file("manifest.tsv", "");

    EXPECT_THROW(sample_manifest manifest(path), std::runtime_error)
        << "Empty manifest file (no header) must throw, not silently parse.";
}

TEST_F(SampleManifestTest, HeaderOnly_ParsesEmpty) {
    auto path = write_file("manifest.tsv", "file\n");

    sample_manifest manifest(path);
    EXPECT_TRUE(manifest.empty())
        << "Header-only manifest must parse cleanly with zero samples.";
    EXPECT_EQ(manifest.size(), 0u);
}

// ── Header tolerance ────────────────────────────────────────────────

TEST_F(SampleManifestTest, HeaderCaseAndWhitespaceTolerance) {
    touch("a.gtf");
    auto path = write_file("manifest.tsv",
        "  FILE  \t  ID  \t  Type  \n"
        "a.gtf\tSAMPLE_A\tsample\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    EXPECT_EQ(manifest.get_samples()[0].id, "SAMPLE_A")
        << "Case- and whitespace-mangled header columns must still resolve.";
}

// ── Field tolerance ─────────────────────────────────────────────────

TEST_F(SampleManifestTest, DotIsEmptyMarker) {
    touch("a.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\tassay\tspecies\n"
        "a.gtf\tSAMPLE_A\t.\t.\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    const auto& info = manifest.get_samples()[0];
    EXPECT_EQ(info.assay, "")  << "'.' in assay column must read as empty.";
    EXPECT_EQ(info.species, "") << "'.' in species column must read as empty.";
}

TEST_F(SampleManifestTest, FieldWhitespaceTrimmed) {
    touch("a.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\tassay\n"
        "a.gtf\t  SAMPLE_A  \t  RNA-seq  \n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    EXPECT_EQ(manifest.get_samples()[0].id, "SAMPLE_A")
        << "Surrounding whitespace must be stripped from id field.";
    EXPECT_EQ(manifest.get_samples()[0].assay, "RNA-seq")
        << "Surrounding whitespace must be stripped from assay field.";
}

// ── Required field absence ─────────────────────────────────────────

TEST_F(SampleManifestTest, EmptyFilePath_Throws) {
    auto path = write_file("manifest.tsv",
        "file\tid\n"
        "\tSAMPLE_A\n");

    EXPECT_THROW(sample_manifest manifest(path), std::runtime_error)
        << "A row with an empty 'file' field must throw at parse time.";
}

// ── Auto-inferred fields ───────────────────────────────────────────

TEST_F(SampleManifestTest, AutoInferIdFromFilename) {
    touch("brain_tumor_rep01.gtf");
    auto path = write_file("manifest.tsv",
        "file\n"
        "brain_tumor_rep01.gtf\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    EXPECT_EQ(manifest.get_samples()[0].id, "brain_tumor_rep01")
        << "Missing id must default to the file's stem.";
}

TEST_F(SampleManifestTest, AutoInferGroupFromRepSuffix) {
    touch("brain_rep01.gtf");
    touch("brain_rep02.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\ttype\n"
        "brain_rep01.gtf\tbrain_rep01\tsample\n"
        "brain_rep02.gtf\tbrain_rep02\tsample\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 2u);
    EXPECT_EQ(manifest.get_samples()[0].group, "brain")
        << "Sample 'brain_rep01' must auto-infer group 'brain'.";
    EXPECT_EQ(manifest.get_samples()[1].group, "brain")
        << "Sample 'brain_rep02' must auto-infer the same group 'brain'.";
}

TEST_F(SampleManifestTest, NoRepSuffix_GroupEqualsId) {
    touch("solo.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\ttype\n"
        "solo.gtf\tsolo_sample\tsample\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    EXPECT_EQ(manifest.get_samples()[0].group, "solo_sample")
        << "Sample with no _repNN suffix must default to group == id.";
}

// ── Type defaulting ────────────────────────────────────────────────

TEST_F(SampleManifestTest, DefaultTypeIsSample) {
    touch("a.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\n"
        "a.gtf\tSAMPLE_A\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    EXPECT_EQ(manifest.get_samples()[0].type, "sample")
        << "Missing 'type' column must default to 'sample'.";
}

TEST_F(SampleManifestTest, TypeCaseInsensitive) {
    touch("a.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\ttype\n"
        "a.gtf\tANN\tANNOTATION\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    EXPECT_EQ(manifest.get_samples()[0].type, "annotation")
        << "Type values must be lower-cased on read.";
}

// ── Comments and blank lines ───────────────────────────────────────

TEST_F(SampleManifestTest, CommentsAndBlankLinesSkipped) {
    touch("a.gtf");
    touch("b.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\n"
        "# this is a comment\n"
        "\n"
        "a.gtf\tSAMPLE_A\n"
        "# another comment\n"
        "\n"
        "b.gtf\tSAMPLE_B\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 2u)
        << "Comments and blank lines must be silently skipped.";
    EXPECT_EQ(manifest.get_samples()[0].id, "SAMPLE_A");
    EXPECT_EQ(manifest.get_samples()[1].id, "SAMPLE_B");
}

// ── Expression attributes ──────────────────────────────────────────

TEST_F(SampleManifestTest, ExpressionAttributesCommaSplitWithWhitespace) {
    touch("a.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\texpression_attribute\n"
        "a.gtf\tSAMPLE_A\t  cov , TPM ,counts  \n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    const auto& attrs = manifest.get_samples()[0].expression_attributes;
    ASSERT_EQ(attrs.size(), 3u);
    EXPECT_EQ(attrs[0], "cov");
    EXPECT_EQ(attrs[1], "TPM");
    EXPECT_EQ(attrs[2], "counts");
}

TEST_F(SampleManifestTest, EmptyExpressionAttributeYieldsEmptyList) {
    touch("a.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\texpression_attribute\n"
        "a.gtf\tSAMPLE_A\t.\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 1u);
    EXPECT_TRUE(manifest.get_samples()[0].expression_attributes.empty())
        << "'.' in expression_attribute column must yield an empty list.";
}

// ── Duplicate IDs (currently accepted) ─────────────────────────────
//
// The parser does not currently dedupe on `id` — duplicate IDs flow
// through to the registry and may produce surprising downstream
// behavior. This test pins the current state so a future PR adding
// validation has a regression test ready to flip; if duplicate-ID
// detection is added, this assertion should change to EXPECT_THROW.
TEST_F(SampleManifestTest, DuplicateIdsCurrentlyAccepted) {
    touch("a.gtf");
    touch("b.gtf");
    auto path = write_file("manifest.tsv",
        "file\tid\ttype\n"
        "a.gtf\tDUPLICATE\tsample\n"
        "b.gtf\tDUPLICATE\tsample\n");

    sample_manifest manifest(path);
    ASSERT_EQ(manifest.size(), 2u)
        << "Both rows must currently parse — duplicate-id detection is a "
           "documented gap (#61 follow-up).";
    EXPECT_EQ(manifest.get_samples()[0].id, "DUPLICATE");
    EXPECT_EQ(manifest.get_samples()[1].id, "DUPLICATE");
}
