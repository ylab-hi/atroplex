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
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <variant>
#include <map>

// genogrove
#include <genogrove/data_type/interval.hpp>

namespace gdt = genogrove::data_type;

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

    // Additional annotations
    std::string source;              // Source of annotation (e.g., GENCODE, Ensembl)
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

    // Segment structure
    int exon_count;                  // Number of exons in this segment

    // Read support metadata (populated during discovery phase)
    size_t read_coverage;            // Number of reads supporting this segment
    std::vector<std::string> supporting_reads;  // Read IDs that support this segment

    // Provenance tracking
    enum class source_type {
        REFERENCE,   // From GFF/GTF annotation
        DISCOVERED,  // Discovered from read analysis
        FUSION       // User-defined or detected fusion
    } source;

    // Constructors
    segment_feature()
        : exon_count(0), read_coverage(0), source(source_type::REFERENCE) {}

    segment_feature(source_type src)
        : exon_count(0), read_coverage(0), source(src) {}

    // Check if this is a reference segment
    bool is_reference() const { return source == source_type::REFERENCE; }

    // Check if this is a discovered segment
    bool is_discovered() const { return source == source_type::DISCOVERED; }

    // Check if this is a fusion segment
    bool is_fusion() const { return source == source_type::FUSION; }

    // Add read support
    void add_read_support(const std::string& read_id) {
        supporting_reads.push_back(read_id);
        read_coverage++;
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

#endif //ATROPLEX_GENOMIC_FEATURE_HPP