/*
 * Integration test for per-hub PSI + entropy columns in splicing_hubs.tsv.
 *
 * Builds a grove with a known hub topology (one exon with 12 distinct
 * downstream targets across 13 segments in one gene), runs
 * analysis_report::collect with hub streaming armed, and verifies the
 * emitted .psi and .entropy values against hand-computed expectations.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "analysis_report.hpp"
#include "build_gff.hpp"
#include "genomic_feature.hpp"
#include "sample_info.hpp"

namespace fs = std::filesystem;

class HubPsiEntropyTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        const auto* info = ::testing::UnitTest::GetInstance()->current_test_info();
        std::string dir_name = "atroplex_hub_psi_test_";
        dir_name += info ? info->name() : "unknown";
        tmp_dir = fs::temp_directory_path() / dir_name;

        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
        fs::create_directories(tmp_dir, ec);
        ASSERT_FALSE(ec) << "Failed to create tmp_dir: " << ec.message();
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmp_dir, ec);
    }

    fs::path tmp_dir;

    fs::path write_gtf(const std::string& name, const std::string& contents) {
        fs::path p = tmp_dir / name;
        std::ofstream out(p, std::ios::binary | std::ios::trunc);
        EXPECT_TRUE(out.is_open()) << "Failed to open " << p;
        out << contents;
        out.flush();
        out.close();
        EXPECT_TRUE(out) << "Write failed for " << p;
        return p;
    }

    // Build a GTF with 12 transcripts through a hub exon E1[1000-1200],
    // each fanning to a unique downstream exon (triggering hub detection
    // at >10 unique targets), plus 1 transcript that bypasses the hub
    // entirely so PSI < 1.
    fs::path write_hub_fixture() {
        std::ostringstream gtf;

        // Gene spanning everything
        gtf << "chr1\tTEST\tgene\t1000\t25000\t.\t+\t.\t"
            << "gene_id \"G1\"; gene_name \"HubGene\"; gene_biotype \"protein_coding\";\n";

        // 12 transcripts through hub exon E1[1000-1200], each with a
        // unique second exon at positions 2000+i*1000 .. 2200+i*1000
        for (int i = 0; i < 12; ++i) {
            std::string tid = "T" + std::to_string(i + 1);
            size_t e2_start = 2000 + static_cast<size_t>(i) * 1000;
            size_t e2_end = e2_start + 200;
            size_t tx_end = e2_end;

            gtf << "chr1\tTEST\ttranscript\t1000\t" << tx_end << "\t.\t+\t.\t"
                << "gene_id \"G1\"; transcript_id \"" << tid << "\";\n";
            gtf << "chr1\tTEST\texon\t1000\t1200\t.\t+\t.\t"
                << "gene_id \"G1\"; transcript_id \"" << tid << "\"; exon_number \"1\";\n";
            gtf << "chr1\tTEST\texon\t" << e2_start << "\t" << e2_end << "\t.\t+\t.\t"
                << "gene_id \"G1\"; transcript_id \"" << tid << "\"; exon_number \"2\";\n";
        }

        // 1 transcript NOT through the hub — different start exon entirely
        gtf << "chr1\tTEST\ttranscript\t20000\t22200\t.\t+\t.\t"
            << "gene_id \"G1\"; transcript_id \"T_ALT\";\n";
        gtf << "chr1\tTEST\texon\t20000\t20200\t.\t+\t.\t"
            << "gene_id \"G1\"; transcript_id \"T_ALT\"; exon_number \"1\";\n";
        gtf << "chr1\tTEST\texon\t22000\t22200\t.\t+\t.\t"
            << "gene_id \"G1\"; transcript_id \"T_ALT\"; exon_number \"2\";\n";

        return write_gtf("hub_fixture.gtf", gtf.str());
    }

    // Parse a TSV file into header + rows. Returns {header_tokens, rows_of_tokens}.
    static std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    parse_tsv(const fs::path& path) {
        std::ifstream in(path);
        std::string line;
        std::vector<std::string> header;
        std::vector<std::vector<std::string>> rows;

        if (std::getline(in, line)) {
            std::istringstream hs(line);
            std::string tok;
            while (std::getline(hs, tok, '\t')) header.push_back(tok);
        }
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream rs(line);
            std::string tok;
            std::vector<std::string> row;
            while (std::getline(rs, tok, '\t')) row.push_back(tok);
            rows.push_back(std::move(row));
        }
        return {header, rows};
    }

    // Find column index by name, or -1 if not found.
    static int find_col(const std::vector<std::string>& header, const std::string& name) {
        for (size_t i = 0; i < header.size(); ++i) {
            if (header[i] == name) return static_cast<int>(i);
        }
        return -1;
    }
};

TEST_F(HubPsiEntropyTest, PsiAndEntropyCorrectForKnownHub) {
    auto gtf_path = write_hub_fixture();

    // Register one sample and build the grove
    sample_info info("hub_sample");
    info.type = "sample";
    uint32_t sample_id = sample_registry::instance().register_data(info);

    grove_type grove(3);
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    size_t segment_count = 0;
    build_counters counters;

    build_gff::build(grove, gtf_path, sample_id, exon_caches, segment_caches,
                     segment_count, 0, expression_filters{}, false, 5,
                     /*include_scaffolds=*/true, counters);

    // Verify we got the expected number of segments: 12 hub-using + 1 non-hub = 13
    ASSERT_EQ(segment_count, 13u) << "Expected 13 segments (12 through hub + 1 bypass)";

    // Set up analysis_report with hub streaming, scoped so the report
    // destructs (flushing the ofstream buffer to disk) before we parse.
    fs::path hubs_path = tmp_dir / "hubs.tsv";
    fs::path branches_path = tmp_dir / "branches.tsv";

    {
        analysis_report report;
        report.begin_splicing_hub_streams(hubs_path.string(), branches_path.string());
        report.collect(grove);
    } // report destructs here → hub_stream flushes + closes

    // Verify hub file was written
    ASSERT_TRUE(fs::exists(hubs_path)) << "splicing_hubs.tsv not created";
    ASSERT_GT(fs::file_size(hubs_path), 0u) << "Hub file is empty";

    // Parse the TSV
    auto [header, rows] = parse_tsv(hubs_path);

    // Should have exactly 1 hub row (only E1 qualifies with >10 targets)
    ASSERT_EQ(rows.size(), 1u) << "Expected exactly 1 hub row";

    // Find the sample's .psi and .entropy columns
    int psi_col = find_col(header, "hub_sample.psi");
    int entropy_col = find_col(header, "hub_sample.entropy");
    ASSERT_GE(psi_col, 0) << "hub_sample.psi column not found in header";
    ASSERT_GE(entropy_col, 0) << "hub_sample.entropy column not found in header";

    // Also verify branches column exists
    int branches_col = find_col(header, "hub_sample.branches");
    ASSERT_GE(branches_col, 0) << "hub_sample.branches column not found in header";

    // Extract values from the single hub row
    ASSERT_LT(static_cast<size_t>(psi_col), rows[0].size());
    ASSERT_LT(static_cast<size_t>(entropy_col), rows[0].size());
    ASSERT_LT(static_cast<size_t>(branches_col), rows[0].size());

    // Values should not be missing
    EXPECT_NE(rows[0][psi_col], ".") << "PSI should not be missing for the sample";
    EXPECT_NE(rows[0][entropy_col], ".") << "Entropy should not be missing";

    double psi = std::stod(rows[0][psi_col]);
    double entropy = std::stod(rows[0][entropy_col]);
    int branches = std::stoi(rows[0][branches_col]);

    // PSI = 12/13 ≈ 0.923 — 12 segments use the hub, 13 total in gene
    double expected_psi = 12.0 / 13.0;
    EXPECT_NEAR(psi, expected_psi, 0.01)
        << "PSI should be 12/13 (12 hub-using segments out of 13 total)";

    // Entropy = log2(12) ≈ 3.585 — all 12 branches used once (uniform)
    double expected_entropy = std::log2(12.0);
    EXPECT_NEAR(entropy, expected_entropy, 0.01)
        << "Entropy should be log2(12) (uniform distribution over 12 branches)";

    // Branches = 12 — the sample uses all 12 downstream targets
    EXPECT_EQ(branches, 12) << "Sample should see all 12 downstream targets";

    // Verify total_branches in the global columns
    int total_branches_col = find_col(header, "total_branches");
    ASSERT_GE(total_branches_col, 0);
    EXPECT_EQ(std::stoi(rows[0][total_branches_col]), 12)
        << "Hub should have 12 total downstream targets";
}