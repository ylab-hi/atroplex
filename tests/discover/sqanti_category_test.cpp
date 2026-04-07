/*
 * Tests for transcript classification (SQANTI-like categories).
 *
 * Reference GTF (reference_transcripts.gtf) defines:
 *   GENE_A on chr22:10000-14500 with two transcripts:
 *     TX_A1: 4 exons [10000-10200, 11000-11300, 12500-12800, 14000-14500]
 *            Junctions: (10200,11000), (11300,12500), (12800,14000)
 *     TX_A2: 3 exons [10000-10200, 11000-11300, 14000-14500]  (skips exon 3)
 *            Junctions: (10200,11000), (11300,14000)
 *   GENE_B on chr22:50000-55000 with one transcript:
 *     TX_B1: 3 exons [50000-50300, 52000-52400, 54000-55000]
 *            Junctions: (50300,52000), (52400,54000)
 *
 * Known donors:    {10200, 11300, 12800, 50300, 52400}
 * Known acceptors: {11000, 12500, 14000, 52000, 54000}
 *
 * Tests construct read_cluster objects directly (no BAM I/O) and verify
 * that transcript_matcher assigns the correct structural category.
 */

#include <gtest/gtest.h>
#include <filesystem>
#include <memory>

#include "genomic_feature.hpp"
#include "build_gff.hpp"
#include "sample_info.hpp"
#include "read_cluster.hpp"
#include "transcript_matcher.hpp"

#include <genogrove/io/bam_reader.hpp>

namespace fs = std::filesystem;

class DiscoverCategoryTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        grove = std::make_unique<grove_type>(3);

        fs::path fixture = fs::path(TEST_FIXTURE_DIR) / "reference_transcripts.gtf";
        ASSERT_TRUE(fs::exists(fixture)) << "Fixture not found: " << fixture;

        sample_info info("test_ref");
        info.type = "annotation";
        info.annotation_source = "TEST";
        uint32_t sample_id = sample_registry::instance().register_data(info);

        // absorb=false to keep all segments as separate candidates
        build_gff::build(*grove, fixture, sample_id, exon_caches, segment_caches,
                         gene_indices, segment_count, 0, -1.0f, false);

        ASSERT_GT(segment_count, 0) << "No segments built from fixture";
    }

    // Helper: construct a read_cluster with specific junctions
    static read_cluster make_cluster(const std::string& id,
                                     const std::string& seqid,
                                     char strand,
                                     size_t start, size_t end,
                                     std::vector<splice_junction> junctions) {
        read_cluster cluster;
        cluster.cluster_id = id;
        cluster.seqid = seqid;
        cluster.strand = strand;
        cluster.start = start;
        cluster.end = end;
        cluster.consensus_junctions = std::move(junctions);
        return cluster;
    }

    // Helper: match a cluster and return the result
    match_result classify(const read_cluster& cluster,
                          int min_overlap = 50) {
        transcript_matcher::config cfg;
        cfg.min_overlap_bp = min_overlap;
        transcript_matcher matcher(*grove, cfg, exon_caches);
        return matcher.match(cluster);
    }

    std::unique_ptr<grove_type> grove;
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    chromosome_gene_segment_indices gene_indices;
    size_t segment_count = 0;
};

// ── FSM: all junctions match reference, same count ─────────────────

TEST_F(DiscoverCategoryTest, FSM_AllJunctionsMatch) {
    // 3 junctions matching TX_A1 exactly
    auto cluster = make_cluster("fsm", "chr22", '+', 10000, 14500,
        {{10200, 11000}, {11300, 12500}, {12800, 14000}});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::FSM)
        << "All 3 junctions match TX_A1 — should be FSM";
    EXPECT_TRUE(result.reference_gene.has_value());
    EXPECT_EQ(result.total_query_junctions, 3);
    EXPECT_EQ(result.total_ref_junctions, 3);
    EXPECT_GE(result.junction_match_score, 0.99);
}

// ── ISM: subset of reference junctions ─────────────────────────────

TEST_F(DiscoverCategoryTest, ISM_Prefix_MissingLastJunction) {
    // First 2 of 3 junctions from TX_A1 — matches first, misses last
    auto cluster = make_cluster("ism_prefix", "chr22", '+', 10000, 12800,
        {{10200, 11000}, {11300, 12500}});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::ISM)
        << "2 of 3 junctions, missing 3' end — should be ISM";
    EXPECT_EQ(result.subcat.ism, ism_subcategory::PREFIX)
        << "Matches first junction, misses last — PREFIX (5' fragment)";
}

TEST_F(DiscoverCategoryTest, ISM_Suffix_MissingFirstJunction) {
    // Last 2 of 3 junctions from TX_A1 — misses first, matches last
    auto cluster = make_cluster("ism_suffix", "chr22", '+', 11000, 14500,
        {{11300, 12500}, {12800, 14000}});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::ISM)
        << "2 of 3 junctions, missing 5' end — should be ISM";
    EXPECT_EQ(result.subcat.ism, ism_subcategory::SUFFIX)
        << "Misses first junction, matches last — SUFFIX (3' fragment)";
}

// ── NIC: novel combination of known splice sites ───────────────────

TEST_F(DiscoverCategoryTest, NIC_NovelCombinationKnownSites) {
    // Junction (10200, 12500): donor 10200 is known (exon 1 end),
    // acceptor 12500 is known (exon 3 start), but this specific pair
    // doesn't exist in any reference transcript (exon skipping)
    auto cluster = make_cluster("nic", "chr22", '+', 10000, 14500,
        {{10200, 12500}, {12800, 14000}});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::NIC)
        << "Known splice sites in novel combination — should be NIC";
    EXPECT_EQ(result.novel_donors, 0) << "All donors are known";
    EXPECT_EQ(result.novel_acceptors, 0) << "All acceptors are known";
}

// ── NNC: at least one novel splice site ────────────────────────────

TEST_F(DiscoverCategoryTest, NNC_NovelDonorSite) {
    // Junction (11500, 12500): donor 11500 is NOT near any known donor
    // (known donors: 10200, 11300, 12800 — all >10bp away)
    auto cluster = make_cluster("nnc", "chr22", '+', 10000, 14500,
        {{10200, 11000}, {11500, 12500}, {12800, 14000}});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::NNC)
        << "Novel donor at 11500 — should be NNC";
    EXPECT_GT(result.novel_donors, 0) << "Should detect novel donor";
}

TEST_F(DiscoverCategoryTest, NNC_NovelAcceptorSite) {
    // Junction (10200, 10800): acceptor 10800 is NOT near any known acceptor
    // (known acceptors: 11000, 12500, 14000 — all >10bp away)
    auto cluster = make_cluster("nnc_acc", "chr22", '+', 10000, 14500,
        {{10200, 10800}, {11300, 12500}, {12800, 14000}});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::NNC)
        << "Novel acceptor at 10800 — should be NNC";
    EXPECT_GT(result.novel_acceptors, 0) << "Should detect novel acceptor";
}

// ── GENIC_INTRON: mono-exon read within intron ─────────────────────

TEST_F(DiscoverCategoryTest, GenicIntron_MonoExonInIntron) {
    // Single-exon read at 10400-10700, within the intron between
    // exon 1 (10000-10200) and exon 2 (11000-11300) of TX_A1.
    // The segment spans 10000-14500, so spatial query finds it.
    auto cluster = make_cluster("genic_intron", "chr22", '+', 10400, 10700, {});

    auto result = classify(cluster, /*min_overlap=*/1);

    EXPECT_EQ(result.category, structural_category::GENIC_INTRON)
        << "Mono-exon read in intron, no exon overlap — should be GENIC_INTRON";
}

// ── GENIC_GENOMIC: mono-exon read overlapping exon ─────────────────

TEST_F(DiscoverCategoryTest, GenicGenomic_MonoExonOverlapsExon) {
    // Single-exon read at 10100-10400, overlapping exon 1 end (10000-10200)
    // and extending into the intron. No junctions.
    auto cluster = make_cluster("genic_genomic", "chr22", '+', 10100, 10400, {});

    auto result = classify(cluster, /*min_overlap=*/1);

    EXPECT_EQ(result.category, structural_category::GENIC_GENOMIC)
        << "Mono-exon read overlapping exon — should be GENIC_GENOMIC";
}

// ── INTERGENIC: no gene overlap ────────────────────────────────────

TEST_F(DiscoverCategoryTest, Intergenic_NoGeneOverlap) {
    // Read at 30000-30500, between GENE_A (ends 14500) and GENE_B (starts 50000)
    auto cluster = make_cluster("intergenic", "chr22", '+', 30000, 30500, {});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::INTERGENIC)
        << "Read in intergenic region — should be INTERGENIC";
}

// ── Cross-gene: verify gene assignment ─────────────────────────────

TEST_F(DiscoverCategoryTest, FSM_GeneB) {
    // Match against GENE_B's TX_B1
    auto cluster = make_cluster("fsm_b", "chr22", '+', 50000, 55000,
        {{50300, 52000}, {52400, 54000}});

    auto result = classify(cluster);

    EXPECT_EQ(result.category, structural_category::FSM)
        << "All junctions match TX_B1 — should be FSM";
    ASSERT_TRUE(result.reference_gene.has_value());
    EXPECT_EQ(result.reference_gene.value(), "GENE_B");
}

// ── Stats accumulation ─────────────────────────────────────────────

TEST_F(DiscoverCategoryTest, MatcherStats_TrackAllCategories) {
    transcript_matcher::config cfg;
    cfg.min_overlap_bp = 1;
    transcript_matcher matcher(*grove, cfg, exon_caches);

    // FSM
    matcher.match(make_cluster("fsm", "chr22", '+', 10000, 14500,
        {{10200, 11000}, {11300, 12500}, {12800, 14000}}));
    // ISM
    matcher.match(make_cluster("ism", "chr22", '+', 10000, 12800,
        {{10200, 11000}, {11300, 12500}}));
    // Intergenic
    matcher.match(make_cluster("inter", "chr22", '+', 30000, 30500, {}));

    const auto& stats = matcher.get_stats();
    EXPECT_EQ(stats.total_matches, 3);
    EXPECT_EQ(stats.fsm_matches, 1);
    EXPECT_EQ(stats.ism_matches, 1);
    EXPECT_EQ(stats.intergenic_matches, 1);
}

// ========================================================================
// Integration tests: full pipeline from SAM file
//
// query_reads.sam contains 23 reads across 8 categories (3 per multi-exon
// category, 2-3 for single-exon). CIGAR strings produce specific junctions
// when parsed by htslib (SAM POS 1-based, htslib converts to 0-based).
//
// Multi-exon (3 reads each):
//   fsm:       200M800N300M1200N300M1200N500M             -> (10200,11000),(11300,12500),(12800,14000)
//   ism_pre:   200M800N300M1200N300M                      -> (10200,11000),(11300,12500)
//   ism_suf:   300M1200N300M1200N500M (pos 11001)         -> (11300,12500),(12800,14000)
//   nic:       200M2300N300M1200N500M                     -> (10200,12500),(12800,14000)
//   nnc:       200M800N500M1000N300M1200N500M             -> (10200,11000),(11500,12500),(12800,14000)
//
// Fuzzy (2 reads):
//   fuzzy_fsm: 202M798N300M1200N300M1200N500M             -> (10202,11000),... +2bp donor wobble
//
// Single-exon (3 reads each):
//   genic_ir:  300M/200M (pos 10401/10451)                -> intron of TX_A1
//   intergen:  500M (pos 30001)                           -> intergenic gap
// ========================================================================

namespace gio = genogrove::io;

class DiscoverIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        transcript_registry::reset();
        gene_registry::reset();
        source_registry::reset();
        sample_registry::reset();

        grove = std::make_unique<grove_type>(3);

        fs::path gtf = fs::path(TEST_FIXTURE_DIR) / "reference_transcripts.gtf";
        ASSERT_TRUE(fs::exists(gtf)) << "GTF fixture not found: " << gtf;

        sam_path = fs::path(TEST_FIXTURE_DIR) / "query_reads.sam";
        ASSERT_TRUE(fs::exists(sam_path)) << "SAM fixture not found: " << sam_path;

        sample_info info("test_ref");
        info.type = "annotation";
        info.annotation_source = "TEST";
        uint32_t sample_id = sample_registry::instance().register_data(info);

        build_gff::build(*grove, gtf, sample_id, exon_caches, segment_caches,
                         gene_indices, segment_count, 0, -1.0f, false);
    }

    std::unique_ptr<grove_type> grove;
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    chromosome_gene_segment_indices gene_indices;
    size_t segment_count = 0;
    fs::path sam_path;
};

// ── Clustering ─────────────────────────────────────────────────────

TEST_F(DiscoverIntegrationTest, Clustering_ReadCounts) {
    gio::bam_reader reader(sam_path.string());

    read_clusterer::config cfg;
    cfg.min_mapq = 20;
    read_clusterer clusterer(cfg);
    auto clusters = clusterer.cluster_reads(reader);

    const auto& stats = clusterer.get_stats();
    EXPECT_EQ(stats.total_reads, 23) << "SAM has 23 reads";
    EXPECT_EQ(stats.filtered_reads, 0) << "All reads have mapq=60, none filtered";
    EXPECT_GT(stats.total_clusters, 0);
}

TEST_F(DiscoverIntegrationTest, Clustering_IdenticalReadsClusterTogether) {
    gio::bam_reader reader(sam_path.string());

    read_clusterer::config cfg;
    cfg.min_mapq = 20;
    read_clusterer clusterer(cfg);
    auto clusters = clusterer.cluster_reads(reader);

    // Count multi-exon clusters with exactly 3 members
    // (fsm, ism_pre, ism_suf, nic, nnc each have 3 identical reads)
    int clusters_with_3 = 0;
    for (const auto& cluster : clusters) {
        if (!cluster.consensus_junctions.empty() && cluster.read_count() == 3) {
            clusters_with_3++;
        }
    }
    EXPECT_GE(clusters_with_3, 4)
        << "At least 4 multi-exon categories should form 3-read clusters";
}

TEST_F(DiscoverIntegrationTest, Clustering_MultiExonClusterCount) {
    gio::bam_reader reader(sam_path.string());

    read_clusterer::config cfg;
    cfg.min_mapq = 20;
    read_clusterer clusterer(cfg);
    auto clusters = clusterer.cluster_reads(reader);

    // fuzzy_fsm reads have a 2bp donor shift -- within junction_tolerance (5bp)
    // they may cluster with FSM reads or form their own cluster.
    // Either way, at least 5 multi-exon clusters total.
    int multi_exon_clusters = 0;
    for (const auto& cluster : clusters) {
        if (!cluster.consensus_junctions.empty()) {
            multi_exon_clusters++;
        }
    }
    EXPECT_GE(multi_exon_clusters, 5)
        << "At least 5 multi-exon clusters (fsm, ism_pre, ism_suf, nic, nnc)";
}

TEST_F(DiscoverIntegrationTest, Clustering_ConsensusJunctions) {
    gio::bam_reader reader(sam_path.string());

    read_clusterer::config cfg;
    read_clusterer clusterer(cfg);
    auto clusters = clusterer.cluster_reads(reader);

    // Find the FSM cluster: 3 junctions, second donor near 11300.
    // NNC cluster also has 3 junctions but second donor is 11500.
    // fuzzy_fsm reads (2bp wobble) merge with exact FSM reads, so
    // the FSM cluster may have 3-5 reads — don't filter by count.
    for (const auto& cluster : clusters) {
        if (cluster.consensus_junctions.size() != 3) continue;
        auto& j = cluster.consensus_junctions;
        if (std::abs(static_cast<double>(j[1].donor) - 11300.0) > 10.0)
            continue;
        EXPECT_GE(cluster.read_count(), 3);
        EXPECT_NEAR(static_cast<double>(j[0].donor),    10200.0, 5.0);
        EXPECT_NEAR(static_cast<double>(j[0].acceptor), 11000.0, 5.0);
        EXPECT_NEAR(static_cast<double>(j[1].donor),    11300.0, 5.0);
        EXPECT_NEAR(static_cast<double>(j[1].acceptor), 12500.0, 5.0);
        EXPECT_NEAR(static_cast<double>(j[2].donor),    12800.0, 5.0);
        EXPECT_NEAR(static_cast<double>(j[2].acceptor), 14000.0, 5.0);
        return;
    }
    ADD_FAILURE() << "Expected a 3-junction cluster with second donor near 11300 (FSM)";
}

// ── Full pipeline: cluster -> match -> classify ────────────────────

TEST_F(DiscoverIntegrationTest, FullPipeline_AllCategories) {
    gio::bam_reader reader(sam_path.string());

    read_clusterer::config cluster_cfg;
    cluster_cfg.min_mapq = 20;
    read_clusterer clusterer(cluster_cfg);
    auto clusters = clusterer.cluster_reads(reader);

    transcript_matcher::config match_cfg;
    match_cfg.min_overlap_bp = 1;
    transcript_matcher matcher(*grove, match_cfg, exon_caches);
    auto results = matcher.match_batch(clusters);

    ASSERT_EQ(results.size(), clusters.size());

    std::map<structural_category, int> counts;
    for (const auto& r : results) {
        counts[r.category]++;
    }

    EXPECT_GE(counts[structural_category::FSM], 1)
        << "fsm reads should produce FSM cluster(s)";
    EXPECT_GE(counts[structural_category::ISM], 1)
        << "ism_pre/ism_suf reads should produce ISM cluster(s)";
    EXPECT_GE(counts[structural_category::NIC], 1)
        << "nic reads should produce NIC cluster(s)";
    EXPECT_GE(counts[structural_category::NNC], 1)
        << "nnc reads should produce NNC cluster(s)";
    EXPECT_GE(counts[structural_category::INTERGENIC], 1)
        << "intergen reads should produce INTERGENIC cluster(s)";

    const auto& mstats = matcher.get_stats();
    EXPECT_EQ(mstats.total_matches, results.size());
    EXPECT_GE(mstats.fsm_matches, 1);
    EXPECT_GE(mstats.ism_matches, 1);
    EXPECT_GE(mstats.nic_matches, 1);
    EXPECT_GE(mstats.nnc_matches, 1);
    EXPECT_GE(mstats.intergenic_matches, 1);
}

TEST_F(DiscoverIntegrationTest, FullPipeline_NovelJunctionsDetected) {
    gio::bam_reader reader(sam_path.string());

    read_clusterer::config cluster_cfg;
    read_clusterer clusterer(cluster_cfg);
    auto clusters = clusterer.cluster_reads(reader);

    transcript_matcher::config match_cfg;
    match_cfg.min_overlap_bp = 1;
    transcript_matcher matcher(*grove, match_cfg, exon_caches);
    auto results = matcher.match_batch(clusters);

    // NNC cluster should report novel donors (11500 is not near any known donor)
    bool found_novel_donor = false;
    for (const auto& r : results) {
        if (r.category == structural_category::NNC && r.novel_donors > 0) {
            found_novel_donor = true;
            break;
        }
    }
    EXPECT_TRUE(found_novel_donor)
        << "NNC cluster from nnc reads should detect novel donor at 11500";
}
