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

// standard
#include <vector>
#include <string>
#include <unordered_set>

// atroplex
#include "genomic_feature.hpp"
#include "read_cluster.hpp"

/**
 * Result of matching a read cluster against reference segments
 */
struct match_result {
    /**
     * Classification of the match
     */
    enum class match_type {
        EXACT,           // All junctions match a known transcript exactly
        COMPATIBLE,      // Junctions are subset of known transcript (partial coverage)
        NOVEL_JUNCTION,  // Has splice junction not in any known transcript
        NOVEL_EXON,      // Extends beyond known exon boundaries
        INTERGENIC,      // No overlap with any known genes
        AMBIGUOUS        // Multiple equally good matches
    };

    match_type type = match_type::INTERGENIC;

    // Best matching segment(s)
    std::vector<key_ptr> matched_segments;
    std::vector<std::string> matched_transcript_ids;

    // Match quality metrics
    double junction_match_score = 0.0;   // Fraction of junctions that match
    int matching_junctions = 0;
    int total_junctions = 0;
    int exon_boundary_mismatches = 0;    // Total bp off from exon boundaries

    // For NOVEL_JUNCTION/NOVEL_EXON types
    std::vector<splice_junction> novel_junctions;

    // Helpers
    bool is_exact() const { return type == match_type::EXACT; }
    bool is_novel() const {
        return type == match_type::NOVEL_JUNCTION || type == match_type::NOVEL_EXON;
    }
    bool has_match() const { return !matched_segments.empty(); }
};

/**
 * Transcript matcher - matches read clusters to reference segments
 *
 * Uses coarse-to-fine query strategy:
 * 1. Spatial query to find candidate segments
 * 2. Drill down to exon chains via graph edges
 * 3. Score junction match quality
 * 4. Classify and optionally update grove
 */
class transcript_matcher {
public:
    /**
     * Configuration for matching behavior
     */
    struct config {
        int junction_tolerance = 5;      // Max bp difference for junction match
        double min_junction_score = 0.8; // Min fraction of junctions to match
        int min_overlap_bp = 50;         // Minimum overlap with segment
        bool allow_partial = true;       // Allow partial transcript matches
        bool track_novel = true;         // Track novel junctions
    };

    /**
     * Construct matcher with reference to grove
     * @param grove Reference to grove containing segments and exons
     * @param cfg Configuration options
     */
    explicit transcript_matcher(grove_type& grove) : transcript_matcher(grove, config{}) {}
    transcript_matcher(grove_type& grove, const config& cfg);

    /**
     * Match a read cluster to reference transcripts
     * @param cluster The read cluster to match
     * @return Match result with classification and matched transcripts
     */
    match_result match(const read_cluster& cluster);

    /**
     * Batch match multiple clusters
     * @param clusters Vector of clusters to match
     * @return Vector of match results (same order as input)
     */
    std::vector<match_result> match_batch(const std::vector<read_cluster>& clusters);

    /**
     * Update grove with match results
     * Adds read support to matched segments, optionally creates discovered segments
     * @param cluster The read cluster
     * @param result The match result
     */
    void update_grove(const read_cluster& cluster, const match_result& result);

    /**
     * Matching statistics
     */
    struct stats {
        size_t total_matches = 0;
        size_t exact_matches = 0;
        size_t compatible_matches = 0;
        size_t novel_junction_matches = 0;
        size_t novel_exon_matches = 0;
        size_t intergenic_matches = 0;
        size_t ambiguous_matches = 0;
        size_t segments_updated = 0;
        size_t segments_created = 0;
    };

    const stats& get_stats() const { return stats_; }

    /**
     * Write match results to TSV file
     * @param filepath Output file path
     * @param clusters The read clusters that were matched
     * @param results The match results (same order as clusters)
     */
    static void write_results(const std::string& filepath,
                              const std::vector<read_cluster>& clusters,
                              const std::vector<match_result>& results);

    /**
     * Write summary statistics to file
     * @param filepath Output file path
     */
    void write_summary(const std::string& filepath) const;

private:
    grove_type& grove_;
    config cfg_;
    stats stats_;

    /**
     * Step 1: Find candidate segments via spatial query
     * @param cluster Read cluster to query
     * @return Vector of segment key pointers that overlap
     */
    std::vector<key_ptr> find_candidate_segments(const read_cluster& cluster);

    /**
     * Step 2: Get exon chain for a segment/transcript
     * Traverses SEGMENT_TO_EXON and EXON_TO_EXON edges
     * @param segment Segment to drill down from
     * @param transcript_id Transcript to follow
     * @return Ordered vector of exon key pointers
     */
    std::vector<key_ptr> get_exon_chain(key_ptr segment, const std::string& transcript_id);

    /**
     * Step 3: Score how well cluster junctions match exon chain
     * @param cluster Read cluster with junctions
     * @param exon_chain Ordered exon keys
     * @return Fraction of matching junctions [0.0, 1.0]
     */
    double score_junction_match(const read_cluster& cluster,
                                const std::vector<key_ptr>& exon_chain);

    /**
     * Step 4: Classify match based on score and junction analysis
     * @param cluster The read cluster
     * @param best_chain Best matching exon chain
     * @param best_score Score from score_junction_match
     * @return Classified match type
     */
    match_result::match_type classify_match(
        const read_cluster& cluster,
        const std::vector<key_ptr>& best_chain,
        double best_score);

    /**
     * Find novel junctions not present in reference
     * @param cluster Read cluster with junctions
     * @param exon_chain Best matching exon chain
     * @return Vector of novel junctions
     */
    std::vector<splice_junction> find_novel_junctions(
        const read_cluster& cluster,
        const std::vector<key_ptr>& exon_chain);

    /**
     * Create a new discovered segment for novel junction patterns
     * @param cluster Read cluster defining the segment
     * @return Key pointer to new segment (or nullptr if creation failed)
     */
    key_ptr create_discovered_segment(const read_cluster& cluster);

    /**
     * Extract reference junctions from exon chain
     * Junctions are boundaries between consecutive exons
     */
    std::vector<splice_junction> extract_reference_junctions(
        const std::vector<key_ptr>& exon_chain);
};

#endif // ATROPLEX_TRANSCRIPT_MATCHER_HPP