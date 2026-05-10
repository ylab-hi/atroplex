/*
 * Tests for `subcall::query::write_classification` — the TSV writer for
 * `query.tsv`. Drives the writer with synthetic `vector<query_result>`
 * inputs and reads the output file back to assert schema invariants
 * around `--nearest` (atroplex#53):
 *
 *   - Without --nearest (no row has closest_*): header omits the four
 *     closest_* columns, INTERGENIC rows are skipped.
 *   - With --nearest (some row has closest_*): header includes the four
 *     closest_* columns, INTERGENIC rows that picked up a flanking
 *     neighbor are emitted, INTERGENIC rows without a neighbor are still
 *     skipped, and present-but-empty `closest_gene_name` coalesces to "."
 *     (writer's opt_str helper).
 *
 * Reader is passed as nullptr so per-sample expression columns are
 * disabled — keeps the schema we're asserting compact.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "subcall/query.hpp"
#include "transcript_matcher.hpp"

namespace fs = std::filesystem;

namespace {

/// Read the whole file into a string.
std::string slurp(const fs::path& p) {
    std::ifstream in(p);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

/// Split text into lines.
std::vector<std::string> split_lines(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string line;
    while (std::getline(in, line)) out.push_back(line);
    return out;
}

/// Split a line on `\t`.
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

class WriteClassificationTest : public ::testing::Test {
protected:
    void SetUp() override {
        sample_registry::reset();
        gene_registry::reset();
        transcript_registry::reset();
        source_registry::reset();
        tmp = fs::temp_directory_path() /
              ("write_classification_" +
               std::string(::testing::UnitTest::GetInstance()
                              ->current_test_info()->name()) +
               ".tsv");
    }

    void TearDown() override {
        if (fs::exists(tmp)) fs::remove(tmp);
    }

    /// Build a query_result for an FSM/NIC-style row that survives the
    /// writer's "skip unmatched" filter (anything other than INTERGENIC
    /// or ANTISENSE without closest is emitted).
    subcall::query_result make_classified_row(
        const std::string& tx,
        structural_category cat,
        const std::string& gene_id,
        const std::string& gene_name) {
        subcall::query_result r;
        r.transcript_id = tx;
        r.category = cat;
        r.gene_id = gene_id;
        r.gene_name = gene_name;
        return r;
    }

    /// Build an INTERGENIC row, optionally with closest_* populated.
    subcall::query_result make_intergenic_row(
        const std::string& tx,
        std::optional<std::string> closest_gene_id = std::nullopt,
        std::optional<std::string> closest_gene_name = std::nullopt,
        std::optional<size_t>      closest_distance = std::nullopt,
        std::optional<std::string> closest_direction = std::nullopt) {
        subcall::query_result r;
        r.transcript_id = tx;
        r.category = structural_category::INTERGENIC;
        r.closest_gene_id     = std::move(closest_gene_id);
        r.closest_gene_name   = std::move(closest_gene_name);
        r.closest_distance_bp = closest_distance;
        r.closest_direction   = std::move(closest_direction);
        return r;
    }

    fs::path tmp;
};

// ── --nearest off (no row has closest_*) ─────────────────────────────

TEST_F(WriteClassificationTest, WithoutNearest_NoClosestColumnsAndIntergenicSkipped) {
    std::vector<subcall::query_result> results;
    results.push_back(make_classified_row("TX_FSM", structural_category::FSM,
                                          "GENE_A", "GeneA"));
    results.push_back(make_intergenic_row("TX_INTERGENIC"));

    subcall::query::write_classification(tmp.string(), results, /*qtx_reader=*/nullptr);

    auto lines = split_lines(slurp(tmp));
    ASSERT_FALSE(lines.empty());

    auto header = split_tabs(lines.front());
    EXPECT_EQ(std::find(header.begin(), header.end(),
                        std::string("closest_gene_id")), header.end())
        << "Header must NOT include closest_* columns when no row has them";

    // Only the FSM row makes it through; INTERGENIC without closest is
    // still skipped.
    ASSERT_EQ(lines.size(), 2u) << "Expect header + 1 emitted row";
    auto fsm_fields = split_tabs(lines[1]);
    EXPECT_EQ(fsm_fields.front(), "TX_FSM");
}

// ── --nearest on, INTERGENIC with neighbor emitted ───────────────────

TEST_F(WriteClassificationTest, WithNearest_HeaderIncludesClosestColumns) {
    std::vector<subcall::query_result> results;
    results.push_back(make_intergenic_row("TX_NEAR",
                                          /*gene_id=*/std::string("GENE_NEAR"),
                                          /*gene_name=*/std::string("GeneNear"),
                                          /*distance=*/size_t(123),
                                          /*direction=*/std::string("upstream")));

    subcall::query::write_classification(tmp.string(), results, nullptr);

    auto lines = split_lines(slurp(tmp));
    ASSERT_GE(lines.size(), 2u);

    auto header = split_tabs(lines.front());
    auto find_col = [&](const std::string& name) {
        return std::find(header.begin(), header.end(), name);
    };
    ASSERT_NE(find_col("closest_gene_id"),     header.end());
    ASSERT_NE(find_col("closest_gene_name"),   header.end());
    ASSERT_NE(find_col("closest_distance_bp"), header.end());
    ASSERT_NE(find_col("closest_direction"),   header.end());

    // Row content
    auto row = split_tabs(lines[1]);
    EXPECT_EQ(row.front(), "TX_NEAR");

    auto col_index = [&](const std::string& name) {
        return std::distance(header.begin(), find_col(name));
    };
    EXPECT_EQ(row[col_index("closest_gene_id")],     "GENE_NEAR");
    EXPECT_EQ(row[col_index("closest_gene_name")],   "GeneNear");
    EXPECT_EQ(row[col_index("closest_distance_bp")], "123");
    EXPECT_EQ(row[col_index("closest_direction")],   "upstream");
}

TEST_F(WriteClassificationTest, WithNearest_IntergenicWithoutNeighborStillSkipped) {
    std::vector<subcall::query_result> results;
    // One row with closest set so emit_closest_cols flips on; one
    // INTERGENIC row WITHOUT closest must still be skipped, otherwise
    // we'd write rows where every closest_* column is just "." with no
    // useful information.
    results.push_back(make_intergenic_row("TX_HAS_NEIGHBOR",
                                          std::string("GENE_X"), std::string("GeneX"),
                                          size_t(50), std::string("downstream")));
    results.push_back(make_intergenic_row("TX_NO_NEIGHBOR"));

    subcall::query::write_classification(tmp.string(), results, nullptr);
    auto lines = split_lines(slurp(tmp));

    // Header + exactly one data row
    ASSERT_EQ(lines.size(), 2u);
    EXPECT_NE(lines[1].find("TX_HAS_NEIGHBOR"), std::string::npos);
    EXPECT_EQ(lines[1].find("TX_NO_NEIGHBOR"),  std::string::npos);
}

// ── opt_str coalescing for present-but-empty strings ─────────────────

TEST_F(WriteClassificationTest, EmptyClosestGeneNameCoalescesToDot) {
    // Real-world case: a registry entry with no `gene_name` attribute in
    // the source GFF (TALON, StringTie outputs without the attribute).
    // segment_feature::gene_name() returns "" — match_result and then
    // query_result carry std::optional<std::string>{""} (present, empty).
    // The writer must emit "." not a literal blank field, otherwise the
    // column count goes off-by-one downstream.
    std::vector<subcall::query_result> results;
    results.push_back(make_intergenic_row("TX_NO_NAME",
                                          /*gene_id=*/std::string("MSTRG.42"),
                                          /*gene_name=*/std::string(""),  // present-but-empty
                                          /*distance=*/size_t(1000),
                                          /*direction=*/std::string("upstream")));

    subcall::query::write_classification(tmp.string(), results, nullptr);
    auto lines = split_lines(slurp(tmp));
    ASSERT_EQ(lines.size(), 2u);

    auto header = split_tabs(lines.front());
    auto row    = split_tabs(lines[1]);
    auto idx    = [&](const std::string& n) {
        return std::distance(header.begin(),
                             std::find(header.begin(), header.end(), n));
    };

    EXPECT_EQ(row[idx("closest_gene_id")],   "MSTRG.42");
    EXPECT_EQ(row[idx("closest_gene_name")], ".")
        << "present-but-empty optional must coalesce to '.'";
    EXPECT_EQ(row[idx("closest_distance_bp")], "1000");
    EXPECT_EQ(row[idx("closest_direction")],   "upstream");
}

// ── column-count consistency between header and rows ────────────────

TEST_F(WriteClassificationTest, RowFieldCountMatchesHeader) {
    // Mix matched + intergenic-with-neighbor — both go through different
    // code paths (the INTERGENIC branch synthesizes "."/0 for non-closest
    // fields). All rows must have the same field count as the header.
    std::vector<subcall::query_result> results;
    results.push_back(make_classified_row("TX_FSM", structural_category::FSM,
                                          "GENE_A", "GeneA"));
    results.push_back(make_intergenic_row("TX_NEAR",
                                          std::string("GENE_NEAR"), std::string("GeneNear"),
                                          size_t(100), std::string("upstream")));

    subcall::query::write_classification(tmp.string(), results, nullptr);
    auto lines = split_lines(slurp(tmp));
    ASSERT_GE(lines.size(), 3u);

    size_t header_cols = split_tabs(lines.front()).size();
    for (size_t i = 1; i < lines.size(); ++i) {
        EXPECT_EQ(split_tabs(lines[i]).size(), header_cols)
            << "Row " << i << " field count differs from header";
    }
}
