/*
 * Integration test for atroplex#53's two open QC items:
 *
 *   - `atroplex query ... --nearest` adds the four closest_* columns to
 *     query.tsv for INTERGENIC rows.
 *   - `atroplex query ...` (no --nearest) emits every classification
 *     including INTERGENIC, with no closest_* columns.
 *
 * Drives the matcher → query_result → writer chain end-to-end without
 * needing a GTF fixture: builds a synthetic grove via insert_data(),
 * runs a manual classification loop that mirrors the 4-line copy in
 * `query::classify_transcripts` (`match.closest_* → qr.closest_*`),
 * then calls `query::write_classification` and asserts the TSV schema.
 *
 * Complements:
 *   - closest_match_test.cpp     — unit tests for populate_nearest
 *   - write_classification_test.cpp — unit tests for the writer with
 *                                      synthetic query_result inputs
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "genomic_feature.hpp"
#include "read_cluster.hpp"
#include "sample_info.hpp"
#include "subcall/query.hpp"
#include "transcript_matcher.hpp"

namespace fs = std::filesystem;
namespace gdt = genogrove::data_type;

namespace {

std::string slurp(const fs::path& p) {
    std::ifstream in(p);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::vector<std::string> split_lines(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string line;
    while (std::getline(in, line)) out.push_back(line);
    return out;
}

std::vector<std::string> split_tabs(const std::string& s) {
    std::vector<std::string> out;
    std::string field;
    for (char c : s) {
        if (c == '\t') { out.push_back(field); field.clear(); }
        else            { field.push_back(c); }
    }
    out.push_back(field);
    return out;
}

}  // namespace

class NearestIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();
        grove = std::make_unique<grove_type>(3);

        tmp = fs::temp_directory_path() /
              ("nearest_integration_" +
               std::string(::testing::UnitTest::GetInstance()
                              ->current_test_info()->name()) +
               ".tsv");
    }

    void TearDown() override {
        if (fs::exists(tmp)) fs::remove(tmp);
    }

    void insert_segment(const std::string& seqid, char strand,
                        size_t start, size_t end,
                        const std::string& gene_id,
                        const std::string& gene_name) {
        segment_feature seg;
        seg.gene_idx = gene_registry::instance().intern(gene_id, gene_name, "");
        gdt::genomic_coordinate coord(strand, start, end);
        genomic_feature feature = seg;
        grove->insert_data(seqid, coord, feature);
    }

    read_cluster make_query(const std::string& tx_id, const std::string& seqid,
                             char strand, size_t start, size_t end) {
        read_cluster c;
        c.cluster_id = tx_id;
        c.seqid = seqid;
        c.strand = strand;
        c.start = start;
        c.end = end;
        size_t mid = start + (end - start) / 2;
        c.consensus_junctions.emplace_back(mid - 5, mid + 5);
        return c;
    }

    /// Mirrors the relevant portion of `query::classify_transcripts` —
    /// the 4-line copy from match_result.closest_* → query_result.closest_*
    /// and the matched-row gene_id population. Kept narrow so this test
    /// exercises the integration without depending on the GTF reader.
    subcall::query_result classify_one(const read_cluster& cluster, bool find_nearest) {
        transcript_matcher::config cfg;
        cfg.junction_tolerance = 5;
        cfg.min_junction_score = 0.8;
        cfg.min_overlap_bp = 50;
        cfg.find_nearest = find_nearest;
        transcript_matcher matcher(*grove, cfg);

        match_result match = matcher.match(cluster);

        subcall::query_result qr;
        qr.transcript_id = cluster.cluster_id;
        qr.category = match.category;
        qr.subcat = match.subcat;

        qr.closest_gene_id     = match.closest_gene_id;
        qr.closest_gene_name   = match.closest_gene_name;
        qr.closest_distance_bp = match.closest_distance_bp;
        qr.closest_direction   = match.closest_direction;

        if (match.has_match()) {
            qr.gene_id = match.reference_gene.value_or(".");
            auto* seg_key = match.matched_segments.front();
            qr.gene_name = get_segment(seg_key->get_data()).gene_name();
        }
        return qr;
    }

    std::unique_ptr<grove_type> grove;
    fs::path tmp;
};

// ── QC item: `atroplex query --nearest` adds closest_* columns for INTERGENIC

TEST_F(NearestIntegrationTest, WithNearest_TsvHasClosestColumnsForIntergenicRow) {
    // Annotation-style segment on chr1; query lands far enough away to
    // be INTERGENIC and picks up the segment as upstream neighbor.
    insert_segment("chr1", '+', 1000, 2000, "GENE_A", "GeneA");

    std::vector<subcall::query_result> results;
    results.push_back(classify_one(
        make_query("TX_INTERGENIC", "chr1", '+', 5000, 6000), /*find_nearest=*/true));

    ASSERT_EQ(results.front().category, structural_category::INTERGENIC)
        << "Query at chr1:5000-6000 with the only segment at 1000-2000 must be INTERGENIC";
    ASSERT_TRUE(results.front().closest_gene_id.has_value())
        << "find_nearest=true must populate closest_gene_id for an INTERGENIC row";

    subcall::query::write_classification(tmp.string(), results, /*qtx_reader=*/nullptr);

    auto lines = split_lines(slurp(tmp));
    ASSERT_EQ(lines.size(), 2u) << "Header + 1 row";

    auto header = split_tabs(lines.front());
    auto col = [&](const std::string& name) {
        return std::distance(header.begin(),
                             std::find(header.begin(), header.end(), name));
    };

    ASSERT_NE(std::find(header.begin(), header.end(),
                        std::string("closest_gene_id")), header.end())
        << "Header MUST include closest_* columns when --nearest yields a hit";
    ASSERT_NE(std::find(header.begin(), header.end(),
                        std::string("closest_distance_bp")), header.end());

    auto row = split_tabs(lines[1]);
    EXPECT_EQ(row[col("transcript_id")],         "TX_INTERGENIC");
    EXPECT_EQ(row[col("structural_category")],   "intergenic");
    EXPECT_EQ(row[col("closest_gene_id")],       "GENE_A");
    EXPECT_EQ(row[col("closest_gene_name")],     "GeneA");
    EXPECT_EQ(row[col("closest_direction")],     "upstream");
    // Closed-coord gap: 5000 - 2000 - 1 = 2999
    EXPECT_EQ(row[col("closest_distance_bp")],   "2999");
}

// ── QC item: `atroplex query` (no --nearest) — every classification emitted,
// ── no closest_* columns

TEST_F(NearestIntegrationTest, WithoutNearest_EveryClassificationEmittedNoClosestCols) {
    insert_segment("chr1", '+', 1000, 2000, "GENE_A", "GeneA");

    // Two queries: one that overlaps (will be scored against the
    // segment) and one that's INTERGENIC. Without --nearest, both must
    // land in the TSV and the closest_* columns must NOT appear.
    std::vector<subcall::query_result> results;
    // Note: the overlap-case won't FSM here (no exon edges via insert_data),
    // so it falls into the second INTERGENIC fall-through. The point of
    // this test is the writer's behavior: every category emitted, no
    // closest_* columns when find_nearest=false. We assert that directly.
    results.push_back(classify_one(
        make_query("TX_FAR", "chr1", '+', 5000, 6000), /*find_nearest=*/false));

    // Append an ANTISENSE row directly so we cover that category too.
    {
        subcall::query_result anti;
        anti.transcript_id = "TX_ANTISENSE";
        anti.category = structural_category::ANTISENSE;
        results.push_back(anti);
    }

    subcall::query::write_classification(tmp.string(), results, /*qtx_reader=*/nullptr);

    auto lines = split_lines(slurp(tmp));
    ASSERT_EQ(lines.size(), 3u) << "Header + 2 rows (INTERGENIC and ANTISENSE both emitted)";

    auto header = split_tabs(lines.front());
    EXPECT_EQ(std::find(header.begin(), header.end(),
                        std::string("closest_gene_id")), header.end())
        << "Header must NOT carry closest_* columns when --nearest is off";
    EXPECT_EQ(std::find(header.begin(), header.end(),
                        std::string("closest_distance_bp")), header.end());

    EXPECT_EQ(split_tabs(lines[1]).front(), "TX_FAR");
    EXPECT_EQ(split_tabs(lines[2]).front(), "TX_ANTISENSE");

    // Find the structural_category column index from the header so this
    // doesn't break if columns are reordered.
    auto cat_col = std::distance(header.begin(),
        std::find(header.begin(), header.end(), std::string("structural_category")));
    EXPECT_EQ(split_tabs(lines[1])[cat_col], "intergenic");
    EXPECT_EQ(split_tabs(lines[2])[cat_col], "antisense");
}
