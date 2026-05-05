/*
 * Tests for absorption rules 0-8.
 * See absorption_rules.txt for full rule documentation.
 *
 * Note: process_gene() sorts transcripts by exon count (descending), so
 * longer transcripts are always processed first within a gene. Tests check
 * final state (live/absorbed counts) not processing order.
 */

#include <gtest/gtest.h>
#include <filesystem>
#include <memory>
#include <set>

#include "genomic_feature.hpp"
#include "build_gff.hpp"
#include "sample_info.hpp"

namespace fs = std::filesystem;

class AbsorptionTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();
    }

    struct BuildResult {
        size_t segment_count;
        size_t absorbed_count;
        size_t live_segments;
        size_t total_absorbed_into;
        size_t live_transcript_count;
        std::unique_ptr<grove_type> grove;
    };

    BuildResult build_fixture(const std::string& filename,
                              bool absorb = true,
                              bool as_annotation = false) {
        fs::path fixture = fs::path(TEST_FIXTURE_DIR) / filename;
        EXPECT_TRUE(fs::exists(fixture)) << "Fixture not found: " << fixture;

        auto grove = std::make_unique<grove_type>(3);
        chromosome_exon_caches exon_caches;
        chromosome_segment_caches segment_caches;
        size_t segment_count = 0;

        sample_info info("test_sample");
        if (as_annotation) {
            info.type = "annotation";
            info.annotation_source = "GENCODE";
        }
        uint32_t sample_id = sample_registry::instance().register_data(info);

        build_counters counters;
        build_options test_opts;
        test_opts.absorb = absorb;
        test_opts.include_scaffolds = true;
        build_gff::build(*grove, fixture, sample_id, exon_caches, segment_caches,
                         segment_count, test_opts, counters);

        size_t absorbed = 0, live = 0, total_absorbed_into = 0, live_tx_count = 0;
        auto roots = grove->get_root_nodes();
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
                    if (seg.absorbed) {
                        absorbed++;
                    } else {
                        live++;
                        total_absorbed_into += seg.absorbed_count;
                        live_tx_count += seg.transcript_ids.size();
                    }
                }
                node = node->get_next();
            }
        }

        return {segment_count, absorbed, live, total_absorbed_into, live_tx_count,
                std::move(grove)};
    }

    // Build two files sequentially to test cross-file absorption with ref/sample distinction
    BuildResult build_two_fixtures(const std::string& file1, bool file1_annotation,
                                   const std::string& file2, bool file2_annotation,
                                   bool absorb = true) {
        fs::path fix1 = fs::path(TEST_FIXTURE_DIR) / file1;
        fs::path fix2 = fs::path(TEST_FIXTURE_DIR) / file2;
        EXPECT_TRUE(fs::exists(fix1)) << "Fixture not found: " << fix1;
        EXPECT_TRUE(fs::exists(fix2)) << "Fixture not found: " << fix2;

        auto grove = std::make_unique<grove_type>(3);
        chromosome_exon_caches exon_caches;
        chromosome_segment_caches segment_caches;
        size_t segment_count = 0;

        build_counters counters;

        // File 1
        sample_info info1("sample1");
        if (file1_annotation) { info1.type = "annotation"; info1.annotation_source = "GENCODE"; }
        uint32_t sid1 = sample_registry::instance().register_data(info1);
        build_options test_opts;
        test_opts.absorb = absorb;
        test_opts.include_scaffolds = true;
        build_gff::build(*grove, fix1, sid1, exon_caches, segment_caches,
                         segment_count, test_opts, counters);

        // File 2
        sample_info info2("sample2");
        if (file2_annotation) { info2.type = "annotation"; info2.annotation_source = "GENCODE"; }
        uint32_t sid2 = sample_registry::instance().register_data(info2);
        build_gff::build(*grove, fix2, sid2, exon_caches, segment_caches,
                         segment_count, test_opts, counters);

        size_t absorbed = 0, live = 0, total_absorbed_into = 0, live_tx_count = 0;
        auto roots = grove->get_root_nodes();
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
                    if (seg.absorbed) {
                        absorbed++;
                    } else {
                        live++;
                        total_absorbed_into += seg.absorbed_count;
                        live_tx_count += seg.transcript_ids.size();
                    }
                }
                node = node->get_next();
            }
        }

        return {segment_count, absorbed, live, total_absorbed_into, live_tx_count,
                std::move(grove)};
    }
};

// ── Rule 0: FSM — identical structure → merge metadata ──────────────

TEST_F(AbsorptionTest, Rule0_FSM_Merge) {
    auto result = build_fixture("rule0_fsm.gtf");

    // TX_A and TX_B have identical exon structure → one segment, both transcript IDs
    EXPECT_EQ(result.live_segments, 1) << "FSM should merge into one segment";
    EXPECT_EQ(result.live_transcript_count, 2) << "Both transcript IDs should be on the segment";
    EXPECT_EQ(result.total_absorbed_into, 0) << "FSM merges metadata, not absorption";
}

// ── Rule 1: 5' ISM — contiguous subset at 5' end → keep ────────────

TEST_F(AbsorptionTest, Rule1_5PrimeISM_Kept) {
    auto result = build_fixture("rule6_5prime_ism.gtf");

    // PARENT: 4 exons [E1, E2, E3, E4]
    // ISM_5P: first 3 exons [E1, E2, E3] — 5' ISM, missing E4 from 3' end
    EXPECT_EQ(result.live_segments, 2) << "5' ISM should be kept as separate segment";
    EXPECT_EQ(result.total_absorbed_into, 0);
}

// ── Rule 2: 3' ISM — missing 1-2 exons from 5' end → absorb ───────

TEST_F(AbsorptionTest, Rule2_3PrimeISM_Absorbed) {
    auto result = build_fixture("rule1_3prime_ism.gtf");

    // PARENT: 4 exons [E1, E2, E3, E4]
    // ISM_3P: last 3 exons [E2, E3, E4] — 3' ISM, missing 1 exon from 5' end
    EXPECT_EQ(result.live_segments, 1) << "3' ISM should be absorbed";
    EXPECT_EQ(result.total_absorbed_into, 1);
    EXPECT_EQ(result.live_transcript_count, 2);
}

// ── Rule 3: 3' degradation — missing 3+ exons → drop vs ref ────────

TEST_F(AbsorptionTest, Rule3_3PrimeDegradation_DroppedVsRef) {
    // PARENT: 5 exons, DEGRAD: last 2 exons (missing 3 from 5' end)
    // As annotation → Rule 3 drops the fragment
    auto result = build_fixture("rule3_degradation.gtf", true, true);

    EXPECT_EQ(result.live_segments, 1) << "3' degradation fragment should be dropped vs annotation";
    EXPECT_EQ(result.live_transcript_count, 1) << "Only parent transcript survives";
}

TEST_F(AbsorptionTest, Rule3_3PrimeDegradation_KeptVsSample) {
    // Same fixture but as sample → Rule 3 keeps the fragment
    auto result = build_fixture("rule3_degradation.gtf", true, false);

    EXPECT_EQ(result.live_segments, 2) << "3' degradation fragment should be kept vs sample";
}

// ── Rule 4: Internal fragment — both ends missing → drop vs ref ─────

TEST_F(AbsorptionTest, Rule4_InternalFragment_DroppedVsRef) {
    // PARENT: 5 exons, INTERNAL: middle 3 (missing from both ends)
    // As annotation → Rule 4 drops the fragment
    auto result = build_fixture("rule2_internal_fragment.gtf", true, true);

    EXPECT_EQ(result.live_segments, 1) << "Internal fragment should be dropped vs annotation";
    EXPECT_EQ(result.live_transcript_count, 1) << "Only parent transcript survives";
}

TEST_F(AbsorptionTest, Rule4_InternalFragment_KeptVsSample) {
    // Same fixture as sample → Rule 4 keeps the fragment
    auto result = build_fixture("rule2_internal_fragment.gtf", true, false);

    EXPECT_EQ(result.live_segments, 2) << "Internal fragment should be kept vs sample";
}

// ── Rule 5: Terminal variant — same intron chain, <50bp → absorb ────

TEST_F(AbsorptionTest, Rule5_TerminalVariant_Absorbed) {
    auto result = build_fixture("rule3_terminal_variant.gtf");

    // CANONICAL: exons [1000-1200, 2000-2300, 3500-3800]
    // SHIFTED_SMALL: exons [1020-1200, 2000-2300, 3500-3790]
    // Same intron chain, TSS diff=20bp, TES diff=10bp → absorb
    EXPECT_EQ(result.live_segments, 1);
    EXPECT_EQ(result.total_absorbed_into, 1);
    EXPECT_EQ(result.live_transcript_count, 2);
}

TEST_F(AbsorptionTest, Rule5_TerminalVariant_NotAbsorbed_LargeDiff) {
    auto result = build_fixture("rule3_not_absorbed.gtf");

    // TES differs by 100bp → NOT absorbed
    EXPECT_EQ(result.live_segments, 2);
    EXPECT_EQ(result.total_absorbed_into, 0);
}

// ── Rule 6: Mono-exon gene overlap → drop ───────────────────────────

TEST_F(AbsorptionTest, Rule6_MonoExonGeneOverlap_Dropped) {
    auto result = build_fixture("rule4_mono_exon.gtf");

    // PARENT: 3 exons, MONO_FRAG: [2050-2250] within PARENT's exon2 [2000-2300]
    // Doesn't cross intron → Rule 6 → drop
    EXPECT_EQ(result.live_segments, 1) << "Mono-exon should be dropped";
}

// ── Rule 7: Mono-exon intron retention → keep ───────────────────────

TEST_F(AbsorptionTest, Rule7_MonoExonIntronRetention_Kept) {
    auto result = build_fixture("rule7_intron_retention.gtf");

    // PARENT: 3 exons, MONO_IR: spans from exon1 across intron into exon2
    EXPECT_EQ(result.live_segments, 2)
        << "Intron retention mono-exon should be kept as separate segment";
}

// ── Rule 8: Mono-exon intergenic → drop ─────────────────────────────

TEST_F(AbsorptionTest, Rule8_MonoExonIntergenic_Dropped) {
    auto result = build_fixture("rule8_mono_intergenic.gtf");

    // MULTI: 3 exons → creates a segment
    // MONO_INTER: single exon at chr1:10000-10500, no overlap with G1 → Rule 8 → drop
    EXPECT_EQ(result.live_segments, 1) << "Intergenic mono-exon should be dropped";
}

// ── Fuzzy matching: ≤5bp junction shift → absorb ────────────────────

TEST_F(AbsorptionTest, FuzzyMatch_JunctionShift_Absorbed) {
    auto result = build_fixture("rule5_junction_shift.gtf");

    // CANONICAL: exons [1000-1200, 2000-2300, 3500-3800]
    // SHIFTED_3BP: all junctions shifted 3bp → fuzzy match → absorb
    EXPECT_EQ(result.live_segments, 1);
    EXPECT_EQ(result.total_absorbed_into, 1);
    EXPECT_EQ(result.live_transcript_count, 2);
}

TEST_F(AbsorptionTest, FuzzyMatch_JunctionShift_NotAbsorbed_LargeShift) {
    auto result = build_fixture("rule5_not_absorbed.gtf");

    // 10bp shift > 5bp threshold → NOT absorbed
    EXPECT_EQ(result.live_segments, 2);
    EXPECT_EQ(result.total_absorbed_into, 0);
}

// ── Reverse absorption ──────────────────────────────────────────────

TEST_F(AbsorptionTest, ReverseAbsorption_ParentArrivesLater) {
    auto result = build_fixture("reverse_absorption.gtf");

    // ISM_FIRST (3 exons) listed first, PARENT_LATER (4 exons) listed second
    // Sorted by exon count: PARENT processed first, ISM forward-absorbed
    EXPECT_EQ(result.live_segments, 1);
    EXPECT_EQ(result.total_absorbed_into, 1);
    EXPECT_EQ(result.live_transcript_count, 2);
}

// ── Absorption disabled ─────────────────────────────────────────────

TEST_F(AbsorptionTest, NoAbsorb_KeepsAll) {
    auto result = build_fixture("rule1_3prime_ism.gtf", false);

    EXPECT_EQ(result.live_segments, 2)
        << "Without absorption, all multi-exon segments should be kept";
    EXPECT_EQ(result.total_absorbed_into, 0);
}

// ── Cross-sample gene_idx inheritance via spatial overlap ───────────

TEST_F(AbsorptionTest, CrossSampleInheritance_NovelLocus) {
    // Two non-annotation samples each contribute a novel transcript with a
    // distinct gene_id ("NOVEL_A" / "NOVEL_B") at spatially overlapping loci
    // (chr1:10000-13000 and chr1:10500-12800). Exon coordinates are disjoint
    // so neither dedup nor absorption fires, but spatial intersect returns
    // A's segment as a candidate when B is inserted. With the inheritance
    // fix, B's segment adopts A's gene_idx; without the fix, the two
    // segments accumulate two separate gene_idx entries.
    auto result = build_two_fixtures("cross_sample_novel_a.gtf", false,
                                     "cross_sample_novel_b.gtf", false);

    ASSERT_EQ(result.live_segments, 2)
        << "Expected 2 distinct live segments (different exon chains, no absorption).";

    std::set<uint32_t> distinct_gene_idx;
    auto roots = result.grove->get_root_nodes();
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
                if (!seg.absorbed) distinct_gene_idx.insert(seg.gene_idx);
            }
            node = node->get_next();
        }
    }

    EXPECT_EQ(distinct_gene_idx.size(), 1u)
        << "Sample B's novel locus should inherit Sample A's gene_idx via "
        << "spatial overlap, but found " << distinct_gene_idx.size()
        << " distinct gene_idx values.";
}

// ── Mono-exon gene_idx inheritance via spatial overlap (#75) ────────
//
// Before the fix, candidate gathering was gated on
// `absorb && exon_chain.size() >= 2`, so mono-exon sample transcripts
// hit `apply_gene_idx_inheritance` with an empty candidates list and
// always created a fresh gene_idx. On a 21K-sample build that left
// ~845K mono-exon segments as singleton genes that should have inherited
// from spatially-overlapping multi-exon parents.
//
// Sample A: 2-exon segment at chr1:10000-13000 (NOVEL_DONOR), intron
//           at 11001-11999.
// Sample B: mono-exon at chr1:10500-12500 (NOVEL_RECIPIENT), spans the
//           donor's intron → classify_mono_exon returns INTRON_RETENTION
//           → the segment is created (not dropped).
//
// With the fix, B's mono-exon segment inherits NOVEL_DONOR's gene_idx;
// without the fix, the two segments live with two separate gene_idx
// entries.
TEST_F(AbsorptionTest, CrossSampleInheritance_MonoExonInheritsFromMultiExon) {
    auto result = build_two_fixtures("cross_sample_mono_donor.gtf",     false,
                                     "cross_sample_mono_recipient.gtf", false);

    ASSERT_EQ(result.live_segments, 2)
        << "Expected 2 live segments (multi-exon donor + mono-exon recipient).";

    std::set<uint32_t> distinct_gene_idx;
    auto roots = result.grove->get_root_nodes();
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
                if (!seg.absorbed) distinct_gene_idx.insert(seg.gene_idx);
            }
            node = node->get_next();
        }
    }

    EXPECT_EQ(distinct_gene_idx.size(), 1u)
        << "Mono-exon sample transcript should inherit the multi-exon "
        << "donor's gene_idx via spatial overlap, but found "
        << distinct_gene_idx.size() << " distinct gene_idx values.";
}
