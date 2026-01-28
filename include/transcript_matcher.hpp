/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_TRANSCRIPT_MATCHER_HPP
#define ATROPLEX_TRANSCRIPT_MATCHER_HPP

#include <vector>
#include <string>
#include <unordered_set>
#include <optional>

#include "genomic_feature.hpp"
#include "read_cluster.hpp"

/**
 * SQANTI-like structural categories for transcript classification
 */
enum class structural_category {
    FSM,            // Full Splice Match - all junctions match reference transcript
    ISM,            // Incomplete Splice Match - subset of junctions from reference
    NIC,            // Novel In Catalog - novel combination of known splice sites
    NNC,            // Novel Not in Catalog - at least one novel splice site
    GENIC_INTRON,   // Mono-exon, entirely within intron
    GENIC_GENOMIC,  // Overlaps both intron and exon regions
    ANTISENSE,      // Overlaps gene on opposite strand
    INTERGENIC,     // No gene overlap
    AMBIGUOUS       // Multiple equally good matches (internal use)
};

/**
 * Subcategories for ISM (Incomplete Splice Match)
 */
enum class ism_subcategory {
    NONE,
    PREFIX,         // 5' fragment - missing 3' end junctions
    SUFFIX,         // 3' fragment - missing 5' end junctions
    INTERNAL,       // Internal fragment - missing both ends
    MONO_EXON       // Single exon matching part of multi-exon transcript
};

/**
 * Subcategories for NIC (Novel In Catalog)
 */
enum class nic_subcategory {
    NONE,
    COMBINATION,        // Novel combination of known junctions
    INTRON_RETENTION,   // Intron retention event
    ALT_3PRIME,         // Alternative 3' splice site (known acceptor)
    ALT_5PRIME,         // Alternative 5' splice site (known donor)
    EXON_SKIPPING       // Known exon(s) skipped
};

/**
 * Subcategories for NNC (Novel Not in Catalog)
 */
enum class nnc_subcategory {
    NONE,
    NOVEL_DONOR,        // At least one novel donor site
    NOVEL_ACCEPTOR,     // At least one novel acceptor site
    NOVEL_BOTH,         // Both novel donor and acceptor sites
    NOVEL_EXON          // Novel internal exon
};

/**
 * Subcategory wrapper - holds the appropriate subcategory based on structural category
 */
struct subcategory {
    ism_subcategory ism = ism_subcategory::NONE;
    nic_subcategory nic = nic_subcategory::NONE;
    nnc_subcategory nnc = nnc_subcategory::NONE;

    std::string to_string() const;
};

/**
 * Convert structural category to string
 */
std::string to_string(structural_category cat);

/**
 * Result of matching a read cluster against reference segments
 */
struct match_result {
    // Primary SQANTI-like classification
    structural_category category = structural_category::INTERGENIC;
    subcategory subcat;

    // Best matching segment(s)
    std::vector<key_ptr> matched_segments;
    std::vector<std::string> matched_transcript_ids;
    std::vector<std::string> matched_gene_ids;

    // Reference transcript info (for FSM/ISM)
    std::optional<std::string> reference_transcript;
    std::optional<std::string> reference_gene;

    // Match quality metrics
    double junction_match_score = 0.0;   // Fraction of junctions that match
    int matching_junctions = 0;
    int total_query_junctions = 0;       // Junctions in query (read cluster)
    int total_ref_junctions = 0;         // Junctions in reference transcript
    int exon_boundary_mismatches = 0;    // Total bp off from exon boundaries

    // Splice site analysis (for NIC/NNC classification)
    int known_donors = 0;        // Query donors matching known donor sites
    int known_acceptors = 0;     // Query acceptors matching known acceptor sites
    int novel_donors = 0;        // Query donors not in reference
    int novel_acceptors = 0;     // Query acceptors not in reference

    // Novel junctions found
    std::vector<splice_junction> novel_junctions;

    // Helpers
    bool is_fsm() const { return category == structural_category::FSM; }
    bool is_ism() const { return category == structural_category::ISM; }
    bool is_nic() const { return category == structural_category::NIC; }
    bool is_nnc() const { return category == structural_category::NNC; }
    bool is_novel() const { return is_nic() || is_nnc(); }
    bool is_genic() const {
        return category == structural_category::GENIC_INTRON ||
               category == structural_category::GENIC_GENOMIC;
    }
    bool has_match() const { return !matched_segments.empty(); }

    // Get full category string (e.g., "ISM_prefix", "NNC_novel_donor")
    std::string category_string() const;
};

/**
 * Transcript matcher - matches read clusters to reference segments
 *
 * Uses coarse-to-fine query strategy:
 * 1. Spatial query to find candidate segments
 * 2. Drill down to exon chains via graph edges
 * 3. Score junction match quality
 * 4. Classify using SQANTI-like categories
 */
class transcript_matcher {
public:
    /**
     * Configuration for matching behavior
     */
    struct config {
        int junction_tolerance = 5;      // Max bp difference for junction match
        double min_junction_score = 0.8; // Min fraction of junctions to match for FSM/ISM
        int min_overlap_bp = 50;         // Minimum overlap with segment
        bool allow_partial = true;       // Allow partial transcript matches (ISM)
        bool track_novel = true;         // Track novel junctions
        int splice_site_window = 10;     // Window for matching splice sites as "known"
    };

    /**
     * Matching statistics (SQANTI categories)
     */
    struct stats {
        size_t total_matches = 0;
        size_t fsm_matches = 0;
        size_t ism_matches = 0;
        size_t nic_matches = 0;
        size_t nnc_matches = 0;
        size_t genic_intron_matches = 0;
        size_t genic_genomic_matches = 0;
        size_t antisense_matches = 0;
        size_t intergenic_matches = 0;
        size_t ambiguous_matches = 0;

        // ISM subcategory counts
        size_t ism_prefix = 0;
        size_t ism_suffix = 0;
        size_t ism_internal = 0;
        size_t ism_mono_exon = 0;

        // NIC subcategory counts
        size_t nic_combination = 0;
        size_t nic_intron_retention = 0;
        size_t nic_alt_3prime = 0;
        size_t nic_alt_5prime = 0;
        size_t nic_exon_skipping = 0;

        // NNC subcategory counts
        size_t nnc_novel_donor = 0;
        size_t nnc_novel_acceptor = 0;
        size_t nnc_novel_both = 0;
        size_t nnc_novel_exon = 0;

        // Grove updates
        size_t segments_updated = 0;
        size_t segments_created = 0;
    };

    explicit transcript_matcher(grove_type& grove) : transcript_matcher(grove, config{}) {}
    transcript_matcher(grove_type& grove, const config& cfg);

    /**
     * Match a read cluster to reference transcripts
     */
    match_result match(const read_cluster& cluster);

    /**
     * Batch match multiple clusters
     */
    std::vector<match_result> match_batch(const std::vector<read_cluster>& clusters);

    /**
     * Update grove with match results
     * @param cluster The read cluster
     * @param result The match result
     * @param add_novel Whether to add novel transcripts to grove
     */
    void update_grove(const read_cluster& cluster, const match_result& result,
                      bool add_novel = false);

    const stats& get_stats() const { return stats_; }

    /**
     * Write match results to SQANTI-like TSV file
     */
    static void write_results(const std::string& filepath,
                              const std::vector<read_cluster>& clusters,
                              const std::vector<match_result>& results);

    /**
     * Write summary statistics to file
     */
    void write_summary(const std::string& filepath) const;

private:
    grove_type& grove_;
    config cfg_;
    stats stats_;

    // Known splice sites from reference (populated during construction)
    std::unordered_set<size_t> known_donor_sites_;    // 5' splice sites
    std::unordered_set<size_t> known_acceptor_sites_; // 3' splice sites

    /**
     * Build index of known splice sites from grove
     */
    void index_splice_sites();

    /**
     * Check if a splice site is known (within tolerance)
     */
    bool is_known_donor(const std::string& seqid, size_t position) const;
    bool is_known_acceptor(const std::string& seqid, size_t position) const;

    /**
     * Find candidate segments via spatial query
     */
    std::vector<key_ptr> find_candidate_segments(const read_cluster& cluster);

    /**
     * Get exon chain for a segment/transcript
     */
    std::vector<key_ptr> get_exon_chain(key_ptr segment, const std::string& transcript_id);

    /**
     * Score how well cluster junctions match exon chain
     */
    double score_junction_match(const read_cluster& cluster,
                                const std::vector<key_ptr>& exon_chain);

    /**
     * Classify match using SQANTI categories
     */
    void classify_match(match_result& result,
                        const read_cluster& cluster,
                        const std::vector<key_ptr>& best_chain);

    /**
     * Determine ISM subcategory based on which junctions are missing
     */
    ism_subcategory classify_ism(const read_cluster& cluster,
                                  const std::vector<key_ptr>& ref_chain);

    /**
     * Determine NIC subcategory
     */
    nic_subcategory classify_nic(const read_cluster& cluster,
                                  const std::vector<key_ptr>& ref_chain);

    /**
     * Determine NNC subcategory based on novel splice sites
     */
    nnc_subcategory classify_nnc(const match_result& result);

    /**
     * Analyze splice sites in query for known/novel classification
     */
    void analyze_splice_sites(match_result& result, const read_cluster& cluster);

    /**
     * Check for intron retention events
     */
    bool has_intron_retention(const read_cluster& cluster,
                              const std::vector<key_ptr>& ref_chain);

    /**
     * Find novel junctions not present in reference
     */
    std::vector<splice_junction> find_novel_junctions(
        const read_cluster& cluster,
        const std::vector<key_ptr>& exon_chain);

    /**
     * Create a new discovered segment for novel patterns
     */
    key_ptr create_discovered_segment(const read_cluster& cluster);

    /**
     * Extract reference junctions from exon chain
     */
    std::vector<splice_junction> extract_reference_junctions(
        const std::vector<key_ptr>& exon_chain);

    /**
     * Update stats based on match result
     */
    void update_stats(const match_result& result);
};

#endif // ATROPLEX_TRANSCRIPT_MATCHER_HPP