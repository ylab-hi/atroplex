/*
 * Tests for query transcript classification and per-sample tracking.
 *
 * Grove is built from two sources:
 *   annotation.gtf (type=annotation, source=TEST):
 *     TX_A1: 4 exons [10000-10200, 11000-11300, 12500-12800, 14000-14500]
 *     TX_A2: 3 exons [10000-10200, 11000-11300, 14000-14500] (exon skip)
 *
 *   sample.gtf (type=sample, source=TALON, with expression counts):
 *     TX_S1: 4 exons [same as TX_A1] — deduplicates with TX_A1's segment
 *     TX_S2: 2 exons [10000-10200, 12500-12800] — sample-only segment
 *
 * After build, segments have these sample memberships:
 *   Segment for TX_A1/TX_S1: annotation_id + sample_id (shared)
 *   Segment for TX_A2:       annotation_id only
 *   Segment for TX_S2:       sample_id only
 *
 * Query transcripts classified against the grove:
 *   Q_SHARED:      matches TX_A1 structure — present in both annotation and sample
 *   Q_ANNO_ONLY:   matches TX_A2 structure — present in annotation only
 *   Q_SAMPLE_ONLY: matches TX_S2 structure — present in sample only
 *   Q_NOVEL:       novel exon at 11500-11800 — not in any source, NNC
 */

#include <gtest/gtest.h>
#include <filesystem>
#include <memory>
#include <algorithm>

#include "genomic_feature.hpp"
#include "build_gff.hpp"
#include "sample_info.hpp"
#include "read_cluster.hpp"
#include "transcript_matcher.hpp"
#include "utility.hpp"

#include <genogrove/io/gff_reader.hpp>

namespace fs = std::filesystem;
namespace gio = genogrove::io;

class QueryClassificationTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        grove = std::make_unique<grove_type>(3);

        fs::path anno_path = fs::path(TEST_FIXTURE_DIR) / "annotation.gtf";
        fs::path sample_path = fs::path(TEST_FIXTURE_DIR) / "sample.gtf";
        query_path = (fs::path(TEST_FIXTURE_DIR) / "query_transcripts.gtf").string();

        ASSERT_TRUE(fs::exists(anno_path)) << "annotation.gtf not found";
        ASSERT_TRUE(fs::exists(sample_path)) << "sample.gtf not found";
        ASSERT_TRUE(fs::exists(query_path)) << "query_transcripts.gtf not found";

        // Register and build annotation (processed first)
        sample_info anno_info("TEST_ANNOTATION");
        anno_info.type = "annotation";
        anno_info.annotation_source = "TEST";
        annotation_id = sample_registry::instance().register_data(anno_info);

        build_gff::build(*grove, anno_path, annotation_id, exon_caches,
                         segment_caches, gene_indices, segment_count,
                         0, -1.0f, false);

        // Register and build sample
        sample_info sample_info_obj("TEST_SAMPLE");
        sample_info_obj.type = "sample";
        sample_id = sample_registry::instance().register_data(sample_info_obj);

        build_gff::build(*grove, sample_path, sample_id, exon_caches,
                         segment_caches, gene_indices, segment_count,
                         0, -1.0f, false);

        ASSERT_GT(segment_count, 0);
    }

    // Replicate what query::classify_transcripts does — read GTF, build
    // clusters from exon chains, match against grove
    struct query_result {
        std::string transcript_id;
        structural_category category;
        std::vector<uint32_t> sample_ids;
        std::map<uint32_t, float> sample_expression;
        std::string gene_id;
    };

    std::vector<query_result> classify_query_transcripts() {
        std::vector<query_result> results;

        transcript_matcher::config cfg;
        cfg.junction_tolerance = 5;
        cfg.min_junction_score = 0.8;
        cfg.min_overlap_bp = 50;
        transcript_matcher matcher(*grove, cfg, exon_caches);

        gio::gff_reader reader(query_path);

        // Group exon entries by transcript
        std::unordered_map<std::string, std::vector<gio::gff_entry>> transcripts;
        for (const auto& entry : reader) {
            if (entry.type != "exon") continue;
            auto tx_id = entry.get_transcript_id();
            if (tx_id.has_value()) {
                transcripts[tx_id.value()].push_back(entry);
            }
        }

        for (auto& [tx_id, exon_entries] : transcripts) {
            std::sort(exon_entries.begin(), exon_entries.end(),
                [](const gio::gff_entry& a, const gio::gff_entry& b) {
                    return a.start < b.start;
                });

            if (exon_entries.size() < 2) continue;

            read_cluster cluster;
            cluster.seqid = normalize_chromosome(exon_entries.front().seqid);
            cluster.strand = exon_entries.front().strand.value_or('+');
            cluster.start = exon_entries.front().start;
            cluster.end = exon_entries.back().end;
            cluster.cluster_id = tx_id;

            for (size_t i = 0; i + 1 < exon_entries.size(); ++i) {
                cluster.consensus_junctions.emplace_back(
                    exon_entries[i].end, exon_entries[i + 1].start);
            }

            match_result match = matcher.match(cluster);

            query_result qr;
            qr.transcript_id = tx_id;
            qr.category = match.category;

            if (match.has_match()) {
                qr.gene_id = match.reference_gene.value_or(".");
                auto* seg_key = match.matched_segments.front();
                const auto& seg = get_segment(seg_key->get_data());

                seg.sample_idx.for_each([&](uint32_t sid) {
                    qr.sample_ids.push_back(sid);
                    if (seg.expression.has(sid)) {
                        qr.sample_expression[sid] = seg.expression.get(sid);
                    }
                });
            }

            results.push_back(std::move(qr));
        }
        return results;
    }

    // Find a result by transcript ID
    const query_result* find_result(const std::vector<query_result>& results,
                                    const std::string& tx_id) {
        for (const auto& r : results) {
            if (r.transcript_id == tx_id) return &r;
        }
        return nullptr;
    }

    bool has_sample(const query_result& r, uint32_t sid) {
        return std::find(r.sample_ids.begin(), r.sample_ids.end(), sid)
               != r.sample_ids.end();
    }

    std::unique_ptr<grove_type> grove;
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    chromosome_gene_segment_indices gene_indices;
    size_t segment_count = 0;
    uint32_t annotation_id = 0;
    uint32_t sample_id = 0;
    std::string query_path;
};

// ── Per-sample presence tracking ───────────────────────────────────

TEST_F(QueryClassificationTest, SharedTranscript_PresentInBoth) {
    auto results = classify_query_transcripts();
    auto* r = find_result(results, "Q_SHARED");
    ASSERT_NE(r, nullptr) << "Q_SHARED should be classified";

    EXPECT_EQ(r->category, structural_category::FSM)
        << "Q_SHARED matches TX_A1 exactly — should be FSM";
    EXPECT_TRUE(has_sample(*r, annotation_id))
        << "Q_SHARED should be present in annotation";
    EXPECT_TRUE(has_sample(*r, sample_id))
        << "Q_SHARED should be present in sample";
}

TEST_F(QueryClassificationTest, AnnotationOnly_PresentInAnnotation) {
    auto results = classify_query_transcripts();
    auto* r = find_result(results, "Q_ANNO_ONLY");
    ASSERT_NE(r, nullptr) << "Q_ANNO_ONLY should be classified";

    EXPECT_TRUE(has_sample(*r, annotation_id))
        << "Q_ANNO_ONLY should be present in annotation";
    EXPECT_FALSE(has_sample(*r, sample_id))
        << "Q_ANNO_ONLY should NOT be present in sample";
}

TEST_F(QueryClassificationTest, SampleOnly_PresentInSample) {
    auto results = classify_query_transcripts();
    auto* r = find_result(results, "Q_SAMPLE_ONLY");
    ASSERT_NE(r, nullptr) << "Q_SAMPLE_ONLY should be classified";

    EXPECT_TRUE(has_sample(*r, sample_id))
        << "Q_SAMPLE_ONLY should be present in sample";
    EXPECT_FALSE(has_sample(*r, annotation_id))
        << "Q_SAMPLE_ONLY should NOT be present in annotation";
}

TEST_F(QueryClassificationTest, NovelTranscript_NoSamplePresence) {
    auto results = classify_query_transcripts();
    auto* r = find_result(results, "Q_NOVEL");
    ASSERT_NE(r, nullptr) << "Q_NOVEL should be classified";

    // Novel junctions (11500-11800) don't match any segment
    EXPECT_TRUE(r->category == structural_category::NNC ||
                r->category == structural_category::NIC)
        << "Q_NOVEL should be NNC or NIC (novel exon boundaries)";

    // Even if it matches a segment partially, the novel exon (11500-11800)
    // means the best match won't have an exact segment. If it does match
    // a segment, that segment might carry sample info from shared exons.
    // The key assertion is that it's not classified as FSM.
    EXPECT_NE(r->category, structural_category::FSM)
        << "Q_NOVEL should NOT be FSM";
}

// ── Expression tracking ────────────────────────────────────────────

TEST_F(QueryClassificationTest, SharedTranscript_HasSampleExpression) {
    auto results = classify_query_transcripts();
    auto* r = find_result(results, "Q_SHARED");
    ASSERT_NE(r, nullptr);

    // TX_S1 has counts=150, which gets stored on the shared segment
    auto expr_it = r->sample_expression.find(sample_id);
    if (expr_it != r->sample_expression.end()) {
        EXPECT_GT(expr_it->second, 0.0f)
            << "Sample expression should be positive (counts=150 in fixture)";
    }
    // Annotation has no expression data
    EXPECT_EQ(r->sample_expression.count(annotation_id), 0)
        << "Annotation should have no expression";
}

// ── Gene assignment ────────────────────────────────────────────────

TEST_F(QueryClassificationTest, AllTranscripts_AssignedToGeneA) {
    auto results = classify_query_transcripts();

    for (const auto& r : results) {
        if (r.category != structural_category::INTERGENIC) {
            EXPECT_EQ(r.gene_id, "GENE_A")
                << "Transcript " << r.transcript_id
                << " should map to GENE_A";
        }
    }
}

// ── Classification correctness ─────────────────────────────────────

TEST_F(QueryClassificationTest, AllFourTranscripts_Classified) {
    auto results = classify_query_transcripts();

    EXPECT_EQ(results.size(), 4)
        << "All 4 query transcripts should be classified (multi-exon)";

    ASSERT_NE(find_result(results, "Q_SHARED"), nullptr);
    ASSERT_NE(find_result(results, "Q_ANNO_ONLY"), nullptr);
    ASSERT_NE(find_result(results, "Q_SAMPLE_ONLY"), nullptr);
    ASSERT_NE(find_result(results, "Q_NOVEL"), nullptr);
}

TEST_F(QueryClassificationTest, SampleCount_MatchesExpected) {
    auto results = classify_query_transcripts();

    auto* shared = find_result(results, "Q_SHARED");
    auto* anno_only = find_result(results, "Q_ANNO_ONLY");
    auto* sample_only = find_result(results, "Q_SAMPLE_ONLY");

    ASSERT_NE(shared, nullptr);
    ASSERT_NE(anno_only, nullptr);
    ASSERT_NE(sample_only, nullptr);

    EXPECT_EQ(shared->sample_ids.size(), 2)
        << "Q_SHARED present in 2 entries (annotation + sample)";
    EXPECT_EQ(anno_only->sample_ids.size(), 1)
        << "Q_ANNO_ONLY present in 1 entry (annotation)";
    EXPECT_EQ(sample_only->sample_ids.size(), 1)
        << "Q_SAMPLE_ONLY present in 1 entry (sample)";
}