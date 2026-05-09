/*
 * Tests for nearest-neighbor reporting on INTERGENIC query results
 * (transcript_matcher::config::find_nearest, atroplex#53).
 *
 * Builds a minimal grove directly via insert_data() so each test controls
 * the exact spatial layout of segments around the query — much easier than
 * constructing GTF fixtures for every edge case.
 *
 * The matcher reads grove.flanking() (genogrove#309 / v0.22.0), picks the
 * closer of predecessor/successor under a same-strand predicate, and
 * populates match_result::closest_{gene_id,gene_name,distance_bp,direction}.
 */

#include <gtest/gtest.h>

#include "genomic_feature.hpp"
#include "read_cluster.hpp"
#include "sample_info.hpp"
#include "transcript_matcher.hpp"

namespace gdt = genogrove::data_type;

class ClosestMatchTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();
        grove = std::make_unique<grove_type>(3);
    }

    // Insert a multi-exon segment as a single genogrove key. The matcher
    // only reads gene_idx / strand / start / end / absorbed off the
    // segment when computing closest_*; we don't need exon edges for that.
    key_ptr insert_segment(const std::string& seqid, char strand,
                            size_t start, size_t end,
                            const std::string& gene_id,
                            const std::string& gene_name = "",
                            bool absorbed = false) {
        segment_feature seg;
        seg.gene_idx = gene_registry::instance().intern(gene_id, gene_name, "");
        seg.absorbed = absorbed;
        gdt::genomic_coordinate coord(strand, start, end);
        genomic_feature feature = seg;
        return grove->insert_data(seqid, coord, feature);
    }

    // Build a 2-exon read_cluster that sits at [start, end] with one
    // junction in the middle. Multi-exon ⇒ goes through the junction
    // path; single-exon would short-circuit via the empty-junctions
    // branch in classify_match().
    read_cluster make_query(const std::string& seqid, char strand,
                             size_t start, size_t end) {
        read_cluster c;
        c.seqid = seqid;
        c.strand = strand;
        c.start = start;
        c.end = end;
        c.cluster_id = "Q";
        size_t mid = start + (end - start) / 2;
        c.consensus_junctions.emplace_back(mid - 5, mid + 5);
        return c;
    }

    transcript_matcher::config config_with_nearest(bool on) {
        transcript_matcher::config cfg;
        cfg.junction_tolerance = 5;
        cfg.min_junction_score = 0.8;
        cfg.min_overlap_bp = 50;
        cfg.find_nearest = on;
        return cfg;
    }

    std::unique_ptr<grove_type> grove;
};

// ── default off: no closest_* even with neighbors present ────────────

TEST_F(ClosestMatchTest, FlagOff_NoClosestPopulated) {
    insert_segment("chr1", '+', 1000, 2000, "GENE_LEFT");
    insert_segment("chr1", '+', 5000, 6000, "GENE_RIGHT");

    transcript_matcher matcher(*grove, config_with_nearest(false));
    auto result = matcher.match(make_query("chr1", '+', 3000, 4000));

    EXPECT_EQ(result.category, structural_category::INTERGENIC);
    EXPECT_FALSE(result.closest_gene_id.has_value());
    EXPECT_FALSE(result.closest_distance_bp.has_value());
    EXPECT_FALSE(result.closest_direction.has_value());
}

// ── basic predecessor / successor selection ──────────────────────────

TEST_F(ClosestMatchTest, PicksClosestUpstream) {
    insert_segment("chr1", '+',  500, 1000, "GENE_FAR");
    insert_segment("chr1", '+', 1500, 2000, "GENE_NEAR");

    transcript_matcher matcher(*grove, config_with_nearest(true));
    auto result = matcher.match(make_query("chr1", '+', 5000, 6000));

    EXPECT_EQ(result.category, structural_category::INTERGENIC);
    ASSERT_TRUE(result.closest_gene_id.has_value());
    EXPECT_EQ(*result.closest_gene_id, "GENE_NEAR");
    EXPECT_EQ(*result.closest_direction, "upstream");
    // Closed-coord gap: query.start - pred.end - 1 = 5000 - 2000 - 1 = 2999
    EXPECT_EQ(*result.closest_distance_bp, 2999u);
}

TEST_F(ClosestMatchTest, PicksClosestDownstream) {
    insert_segment("chr1", '+', 5000, 6000, "GENE_NEAR");
    insert_segment("chr1", '+', 9000, 9500, "GENE_FAR");

    transcript_matcher matcher(*grove, config_with_nearest(true));
    auto result = matcher.match(make_query("chr1", '+', 1000, 2000));

    EXPECT_EQ(result.category, structural_category::INTERGENIC);
    ASSERT_TRUE(result.closest_gene_id.has_value());
    EXPECT_EQ(*result.closest_gene_id, "GENE_NEAR");
    EXPECT_EQ(*result.closest_direction, "downstream");
    // Closed-coord gap: succ.start - query.end - 1 = 5000 - 2000 - 1 = 2999
    EXPECT_EQ(*result.closest_distance_bp, 2999u);
}

TEST_F(ClosestMatchTest, PicksClosestOfBracketingPair) {
    auto* upstream = insert_segment("chr1", '+', 4000, 4500, "GENE_UP");
    auto* downstream = insert_segment("chr1", '+', 5500, 6000, "GENE_DOWN");
    (void)upstream; (void)downstream;

    transcript_matcher matcher(*grove, config_with_nearest(true));
    // Query [4800,5200] → upstream gap = 4800-4500-1 = 299
    //                  →   downstream gap = 5500-5200-1 = 299  (tie!)
    // Tie-break favors upstream (predecessor) for determinism.
    auto tied = matcher.match(make_query("chr1", '+', 4800, 5200));
    EXPECT_EQ(tied.closest_direction.value_or(""), "upstream");
    EXPECT_EQ(tied.closest_gene_id.value_or(""), "GENE_UP");

    // Bias the query toward the right neighbor (gap = 100 vs 499).
    auto downward = matcher.match(make_query("chr1", '+', 5001, 5399));
    EXPECT_EQ(downward.closest_direction.value_or(""), "downstream");
    EXPECT_EQ(downward.closest_gene_id.value_or(""), "GENE_DOWN");
    EXPECT_EQ(downward.closest_distance_bp.value_or(0u), 100u);
}

// ── strand filtering via flanking predicate ──────────────────────────

TEST_F(ClosestMatchTest, OppositeStrandNeighborSkipped) {
    insert_segment("chr1", '-', 4000, 4500, "GENE_MINUS_NEAR");
    insert_segment("chr1", '+', 1000, 2000, "GENE_PLUS_FAR");

    transcript_matcher matcher(*grove, config_with_nearest(true));
    auto result = matcher.match(make_query("chr1", '+', 5000, 6000));

    ASSERT_TRUE(result.closest_gene_id.has_value());
    EXPECT_EQ(*result.closest_gene_id, "GENE_PLUS_FAR")
        << "Same-strand predicate should skip the closer minus-strand segment";
}

TEST_F(ClosestMatchTest, WildcardQueryMatchesAnyStrand) {
    insert_segment("chr1", '-', 4000, 4500, "GENE_MINUS");

    transcript_matcher matcher(*grove, config_with_nearest(true));
    // Cluster strand '*' means any strand is acceptable — see populate_nearest's
    // same_strand predicate.
    auto result = matcher.match(make_query("chr1", '*', 5000, 6000));

    ASSERT_TRUE(result.closest_gene_id.has_value());
    EXPECT_EQ(*result.closest_gene_id, "GENE_MINUS");
}

// ── empty / out-of-grove edges ───────────────────────────────────────

TEST_F(ClosestMatchTest, EmptyGroveLeavesClosestUnset) {
    transcript_matcher matcher(*grove, config_with_nearest(true));
    auto result = matcher.match(make_query("chr1", '+', 5000, 6000));

    EXPECT_EQ(result.category, structural_category::INTERGENIC);
    EXPECT_FALSE(result.closest_gene_id.has_value());
}

TEST_F(ClosestMatchTest, MissingChromosomeLeavesClosestUnset) {
    insert_segment("chr1", '+', 1000, 2000, "GENE_X");

    transcript_matcher matcher(*grove, config_with_nearest(true));
    // Query lives on chr2 — no flanking neighbor on that index.
    auto result = matcher.match(make_query("chr2", '+', 5000, 6000));

    EXPECT_EQ(result.category, structural_category::INTERGENIC);
    EXPECT_FALSE(result.closest_gene_id.has_value());
}


// ── overlapping segments are skipped by flanking ─────────────────────

TEST_F(ClosestMatchTest, OverlappingSegmentNotReportedAsClosest) {
    // grove.flanking() restricts results to keys that do NOT overlap the
    // query. We rely on this: even when the matcher reaches the
    // INTERGENIC branch on a query that spatially overlaps a segment
    // (e.g., the segment has no exon edges so the candidate is dropped),
    // populate_nearest must not "rescue" the overlapping segment as a
    // closest match. Here the only segment overlaps the query, so
    // closest_* must stay empty.
    insert_segment("chr1", '+', 1000, 6000, "GENE_OVERLAPS");

    transcript_matcher matcher(*grove, config_with_nearest(true));
    auto result = matcher.match(make_query("chr1", '+', 2000, 3000));

    EXPECT_FALSE(result.closest_gene_id.has_value())
        << "flanking() must not return overlapping keys";
}
