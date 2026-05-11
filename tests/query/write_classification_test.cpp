/*
 * Tests for `subcall::query::write_classification` — the TSV writer for
 * `query.tsv`. Drives the writer with synthetic `vector<query_result>`
 * inputs and reads the output file back to assert schema invariants
 * around `--nearest` (atroplex#53):
 *
 *   - Every classified row is emitted, INTERGENIC and ANTISENSE
 *     included (the previous skip behavior contradicted #53's "rather
 *     than silently dropping it from the output").
 *   - Without --nearest (no row has closest_*): header omits the four
 *     closest_* columns.
 *   - With --nearest (some row has closest_*): header includes the four
 *     closest_* columns; rows without a flanking neighbor leave them
 *     as "." and present-but-empty `closest_gene_name` coalesces to "."
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

    /// Build a query_result for an FSM/NIC-style row (the writer emits
    /// every classification now; this helper is just for ergonomics).
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

// ── every classification is emitted ──────────────────────────────────

TEST_F(WriteClassificationTest, AllClassificationsEmitted_NoClosestCols) {
    // Without --nearest, no row carries closest_*, so the header omits
    // the four closest_* columns. Every row is still emitted — FSM,
    // INTERGENIC, and ANTISENSE all show up so users can see which
    // input transcripts landed in each category (#53's stated motivation).
    std::vector<subcall::query_result> results;
    results.push_back(make_classified_row("TX_FSM", structural_category::FSM,
                                          "GENE_A", "GeneA"));
    results.push_back(make_intergenic_row("TX_INTERGENIC"));
    {
        subcall::query_result anti;
        anti.transcript_id = "TX_ANTISENSE";
        anti.category = structural_category::ANTISENSE;
        results.push_back(anti);
    }

    subcall::query::write_classification(tmp.string(), results, /*qtx_reader=*/nullptr);

    auto lines = split_lines(slurp(tmp));
    ASSERT_EQ(lines.size(), 4u) << "Expect header + 3 emitted rows";

    auto header = split_tabs(lines.front());
    EXPECT_EQ(std::find(header.begin(), header.end(),
                        std::string("closest_gene_id")), header.end())
        << "Header must NOT include closest_* columns when no row has them";

    EXPECT_EQ(split_tabs(lines[1]).front(), "TX_FSM");
    EXPECT_EQ(split_tabs(lines[2]).front(), "TX_INTERGENIC");
    EXPECT_EQ(split_tabs(lines[3]).front(), "TX_ANTISENSE");
}

// ── --nearest on, INTERGENIC with neighbor emitted ───────────────────

TEST_F(WriteClassificationTest, WithNearest_HeaderIncludesClosestColumns) {
    std::vector<subcall::query_result> results;
    results.push_back(make_intergenic_row("TX_NEAR",
                                          /*closest_gene_id=*/std::string("GENE_NEAR"),
                                          /*closest_gene_name=*/std::string("GeneNear"),
                                          /*closest_distance=*/static_cast<size_t>(123),
                                          /*closest_direction=*/std::string("upstream")));

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

TEST_F(WriteClassificationTest, WithNearest_IntergenicWithoutNeighborStillEmitted) {
    // With --nearest on (at least one row has closest_*), the header
    // includes closest_* columns. An INTERGENIC row that did NOT pick
    // up a flanking neighbor (e.g., on a chromosome with no segments)
    // is still emitted; its closest_* cells are "." rather than the
    // row being silently dropped.
    std::vector<subcall::query_result> results;
    results.push_back(make_intergenic_row("TX_HAS_NEIGHBOR",
                                          std::string("GENE_X"), std::string("GeneX"),
                                          static_cast<size_t>(50),
                                          std::string("downstream")));
    results.push_back(make_intergenic_row("TX_NO_NEIGHBOR"));

    subcall::query::write_classification(tmp.string(), results, nullptr);
    auto lines = split_lines(slurp(tmp));

    ASSERT_EQ(lines.size(), 3u) << "Header + both INTERGENIC rows";

    auto header = split_tabs(lines.front());
    auto col_index = [&](const std::string& name) {
        return std::distance(header.begin(),
                             std::find(header.begin(), header.end(), name));
    };

    auto with_neighbor    = split_tabs(lines[1]);
    auto without_neighbor = split_tabs(lines[2]);

    EXPECT_EQ(with_neighbor.front(),    "TX_HAS_NEIGHBOR");
    EXPECT_EQ(without_neighbor.front(), "TX_NO_NEIGHBOR");

    EXPECT_EQ(with_neighbor[col_index("closest_gene_id")],     "GENE_X");
    EXPECT_EQ(without_neighbor[col_index("closest_gene_id")],  ".")
        << "Row without a flanking hit must show '.' in closest_gene_id";
    EXPECT_EQ(without_neighbor[col_index("closest_distance_bp")], ".");
    EXPECT_EQ(without_neighbor[col_index("closest_direction")],   ".");
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
                                          /*closest_gene_id=*/std::string("MSTRG.42"),
                                          /*closest_gene_name=*/std::string(""),  // present-but-empty
                                          /*closest_distance=*/static_cast<size_t>(1000),
                                          /*closest_direction=*/std::string("upstream")));

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
                                          static_cast<size_t>(100),
                                          std::string("upstream")));

    subcall::query::write_classification(tmp.string(), results, nullptr);
    auto lines = split_lines(slurp(tmp));
    ASSERT_GE(lines.size(), 3u);

    size_t header_cols = split_tabs(lines.front()).size();
    for (size_t i = 1; i < lines.size(); ++i) {
        EXPECT_EQ(split_tabs(lines[i]).size(), header_cols)
            << "Row " << i << " field count differs from header";
    }
}
