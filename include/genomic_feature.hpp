/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_GENOMIC_FEATURE_HPP
#define ATROPLEX_GENOMIC_FEATURE_HPP

// standard
#include <cstdint>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <variant>
#include <map>

// genogrove
#include <genogrove/data_type/interval.hpp>
#include <genogrove/data_type/genomic_coordinate.hpp>
#include <genogrove/structure/grove/grove.hpp>

namespace gdt = genogrove::data_type;
namespace gst = genogrove::structure;

/**
 * Edge metadata for grove graph
 * Tracks which transcript/segment an edge belongs to
 */
struct edge_metadata {
    // Primary identifier for this edge (transcript or segment ID)
    std::string id;

    // Type of edge
    enum class edge_type {
        EXON_TO_EXON,      // Exon chain within a segment
        SEGMENT_TO_EXON,   // Segment links to its first exon
        SEGMENT_TO_SEGMENT // Transcript path through segments
    } type;

    edge_metadata() : type(edge_type::EXON_TO_EXON) {}

    edge_metadata(std::string id, edge_type t)
        : id(std::move(id)), type(t) {}
};

/**
 * Exon feature: spatial interval with biological annotations
 * Used for read overlap queries and CDS/UTR disambiguation
 */
struct exon_feature {
    // Core identifiers
    std::string id;                  // Exon ID
    std::string gene_id;             // Gene identifier
    std::string gene_name;           // Gene name/symbol
    std::string gene_biotype;        // Gene biotype (protein_coding, lncRNA, etc.)
    std::string coordinate;          // Genomic coordinate e.g., chr1:+:100-200

    // Transcript associations
    std::unordered_set<std::string> transcript_ids;  // Reference transcripts using this exon

    // Source tracking (GFF column 2)
    // Multiple sources can call the same exon (e.g., HAVANA, ENSEMBL, StringTie)
    std::unordered_set<std::string> sources;

    // Sample tracking (registry IDs from sample_registry)
    // Links to sample_info entries that contain this exon
    // Use sample_registry::instance().get(id) to retrieve full sample_info
    std::unordered_set<uint32_t> sample_idx;

    // Expression quantification (optional, for samples with expression data)
    // Maps sample registry ID -> expression value (TPM, FPKM, count, etc.)
    std::unordered_map<uint32_t, float> expression;

    // Additional annotations
    int exon_number;                 // Exon number within transcript (-1 for unknown)

    // Overlapping feature annotations (CDS, UTR, codons)
    // Maps feature type (e.g., "CDS", "five_prime_UTR", "start_codon") to intervals
    // Example: exon [100-300] might have:
    //   "CDS" -> {[120-300]}
    //   "five_prime_UTR" -> {[100-119]}
    //   "start_codon" -> {[120-122]}
    std::unordered_map<std::string, std::vector<gdt::interval>> overlapping_features;

    // Constructors
    exon_feature() : exon_number(-1) {}

    // Create exon from GFF entry
    static exon_feature from_gff_entry(
        const std::map<std::string, std::string>& attributes,
        const std::string& seqid,
        const gdt::interval& interval,
        char strand
    );

    // Check if exon overlaps with a specific feature type
    bool has_overlapping_feature(const std::string& feature_type) const {
        return overlapping_features.find(feature_type) != overlapping_features.end();
    }

    // Get intervals where a specific feature type overlaps
    const std::vector<gdt::interval>* get_overlapping_intervals(const std::string& feature_type) const {
        auto it = overlapping_features.find(feature_type);
        return it != overlapping_features.end() ? &it->second : nullptr;
    }

    // --- Sample tracking methods ---

    // Add sample reference (registry ID from sample_registry)
    void add_sample(uint32_t sample_id) {
        sample_idx.insert(sample_id);
    }

    // Add source (GFF column 2 value like HAVANA, ENSEMBL, StringTie)
    void add_source(const std::string& source) {
        sources.insert(source);
    }

    // Check if feature exists in a specific sample
    bool in_sample(uint32_t sample_id) const {
        return sample_idx.contains(sample_id);
    }

    // Check if feature has a specific source
    bool has_source(const std::string& source) const {
        return sources.contains(source);
    }

    // Get number of samples containing this feature
    size_t sample_count() const {
        return sample_idx.size();
    }

    // --- Expression methods ---

    // Set expression value for a sample
    void set_expression(uint32_t sample_id, float value) {
        expression[sample_id] = value;
    }

    // Get expression value for a sample (returns 0 if not present)
    float get_expression(uint32_t sample_id) const {
        auto it = expression.find(sample_id);
        return it != expression.end() ? it->second : 0.0f;
    }

    // Check if expression data exists for a sample
    bool has_expression(uint32_t sample_id) const {
        return expression.contains(sample_id);
    }

    // --- Pan-transcriptome statistics ---

    // Calculate frequency across samples
    float frequency(size_t total_samples) const {
        if (total_samples == 0) return 0.0f;
        return static_cast<float>(sample_idx.size()) / static_cast<float>(total_samples);
    }

    // Check if feature is conserved (present in all samples)
    bool is_conserved(size_t total_samples) const {
        return sample_idx.size() == total_samples;
    }

    // Check if feature is sample-specific (only in one sample)
    bool is_sample_specific() const {
        return sample_idx.size() == 1;
    }
};

/**
 * Segment feature: spans consecutive exons, represents transcript paths
 * Used for graph traversal and transcript assignment
 */
struct segment_feature {
    // Core identifiers
    std::string id;                  // Segment ID
    std::string gene_id;             // Gene identifier
    std::string gene_name;           // Gene name/symbol
    std::string gene_biotype;        // Gene biotype
    std::string coordinate;          // Genomic coordinate e.g., chr1:+:100-500

    // Transcript associations
    std::unordered_set<std::string> transcript_ids;  // Transcripts using this segment

    // Source tracking (GFF column 2)
    // Multiple sources can call the same segment (e.g., HAVANA, ENSEMBL, StringTie)
    std::unordered_set<std::string> sources;

    // Sample tracking (registry IDs from sample_registry)
    // Links to sample_info entries that contain this segment
    // Use sample_registry::instance().get(id) to retrieve full sample_info
    std::unordered_set<uint32_t> sample_idx;

    // Expression quantification (optional, for samples with expression data)
    // Maps sample registry ID -> expression value (TPM, FPKM, count, etc.)
    std::unordered_map<uint32_t, float> expression;

    // Transcript biotype mapping (transcript_id -> biotype)
    // e.g., protein_coding, retained_intron, nonsense_mediated_decay
    std::unordered_map<std::string, std::string> transcript_biotypes;

    // Segment structure
    size_t segment_index;            // Unique numeric index (used as edge ID for traversal)
    int exon_count;                  // Number of exons in this segment

    // Read support metadata (populated during discovery phase)
    size_t read_coverage;            // Number of reads supporting this segment
    std::vector<std::string> supporting_reads;  // Read IDs that support this segment

    // Constructors
    segment_feature()
        : segment_index(0), exon_count(0), read_coverage(0) {}

    // --- Sample tracking methods ---

    // Add sample reference (registry ID from sample_registry)
    void add_sample(uint32_t sample_id) {
        sample_idx.insert(sample_id);
    }

    // Add source (GFF column 2 value like HAVANA, ENSEMBL, StringTie)
    void add_source(const std::string& source) {
        sources.insert(source);
    }

    // Check if feature exists in a specific sample
    bool in_sample(uint32_t sample_id) const {
        return sample_idx.contains(sample_id);
    }

    // Check if feature has a specific source
    bool has_source(const std::string& source) const {
        return sources.contains(source);
    }

    // Get number of samples containing this feature
    size_t sample_count() const {
        return sample_idx.size();
    }

    // --- Expression methods ---

    // Set expression value for a sample
    void set_expression(uint32_t sample_id, float value) {
        expression[sample_id] = value;
    }

    // Get expression value for a sample (returns 0 if not present)
    float get_expression(uint32_t sample_id) const {
        auto it = expression.find(sample_id);
        return it != expression.end() ? it->second : 0.0f;
    }

    // Check if expression data exists for a sample
    bool has_expression(uint32_t sample_id) const {
        return expression.contains(sample_id);
    }

    // --- Read support methods ---

    // Add read support
    void add_read_support(const std::string& read_id) {
        supporting_reads.push_back(read_id);
        read_coverage++;
    }

    // --- Pan-transcriptome statistics ---

    // Calculate frequency across samples
    float frequency(size_t total_samples) const {
        if (total_samples == 0) return 0.0f;
        return static_cast<float>(sample_idx.size()) / static_cast<float>(total_samples);
    }

    // Check if feature is conserved (present in all samples)
    bool is_conserved(size_t total_samples) const {
        return sample_idx.size() == total_samples;
    }

    // Check if feature is sample-specific (only in one sample)
    bool is_sample_specific() const {
        return sample_idx.size() == 1;
    }
};

/**
 * Genomic feature variant: can be either exon or segment
 * Use std::holds_alternative<T> to check type
 * Use std::get<T> to access specific type
 */
using genomic_feature = std::variant<exon_feature, segment_feature>;

// Helper functions to check variant type
inline bool is_exon(const genomic_feature& feature) {
    return std::holds_alternative<exon_feature>(feature);
}

inline bool is_segment(const genomic_feature& feature) {
    return std::holds_alternative<segment_feature>(feature);
}

// Helper functions to get variant value
inline exon_feature& get_exon(genomic_feature& feature) {
    return std::get<exon_feature>(feature);
}

inline const exon_feature& get_exon(const genomic_feature& feature) {
    return std::get<exon_feature>(feature);
}

inline segment_feature& get_segment(genomic_feature& feature) {
    return std::get<segment_feature>(feature);
}

inline const segment_feature& get_segment(const genomic_feature& feature) {
    return std::get<segment_feature>(feature);
}

// ========== Grove type aliases ==========
// Canonical definitions used throughout the codebase

using grove_type = gst::grove<gdt::genomic_coordinate, genomic_feature, edge_metadata>;
using key_ptr = gdt::key<gdt::genomic_coordinate, genomic_feature>*;

// Cache types for deduplication across files (pan-transcriptome)
using exon_cache_type = std::map<gdt::genomic_coordinate, key_ptr>;
using segment_cache_type = std::unordered_map<std::string, key_ptr>;
using chromosome_exon_caches = std::map<std::string, exon_cache_type>;
using chromosome_segment_caches = std::map<std::string, segment_cache_type>;

#endif //ATROPLEX_GENOMIC_FEATURE_HPP