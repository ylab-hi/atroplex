/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_GENOMIC_FEATURE_HPP
#define ATROPLEX_GENOMIC_FEATURE_HPP

#include <string>
#include <unordered_set>

/**
 * Type of genomic feature stored in the grove
 */
enum class feature_type {
    EXON,
    SEGMENT
};

/**
 * Genomic feature data stored in the grove
 * Represents both exons and segments with biological metadata
 * Links to the unified transcript graph via graph_node_id
 */
struct genomic_feature {
    // Reference to unified transcript graph
    size_t graph_node_id;

    // Feature classification
    feature_type type;

    // Core identifiers
    std::string id;              // Feature ID (exon_id or segment_id)
    std::string gene_id;         // Gene identifier
    std::string gene_name;       // Gene name/symbol
    std::string gene_type;       // Gene biotype (protein_coding, lncRNA, etc.)

    // Transcript information
    std::unordered_set<std::string> transcript_ids;  // All transcripts using this feature
    std::string transcript_type;  // Transcript biotype

    // Additional annotations
    std::string source;          // Source of annotation (e.g., GENCODE, Ensembl)
    int exon_number;            // Exon number within transcript (-1 for segments/unknown)

    // Constructors
    genomic_feature()
        : graph_node_id(0), type(feature_type::EXON), exon_number(-1) {}

    genomic_feature(size_t node_id, feature_type t)
        : graph_node_id(node_id), type(t), exon_number(-1) {}

    // Create feature from GFF entry
    static genomic_feature from_gff_entry(
        const std::string& attributes,
        size_t node_id,
        feature_type type
    );

    // Helper to check if this is an exon
    bool is_exon() const { return type == feature_type::EXON; }

    // Helper to check if this is a segment
    bool is_segment() const { return type == feature_type::SEGMENT; }
};

#endif //ATROPLEX_GENOMIC_FEATURE_HPP