/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_TRANSCRIPT_GRAPH_HPP
#define ATROPLEX_TRANSCRIPT_GRAPH_HPP

// standard
#include <string>
#include <vector>
#include <unordered_set>
#include <memory>
#include <variant>

// genogrove
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/gff_reader.hpp>

namespace gdt = genogrove::data_type;
namespace ggs = genogrove::structure;
namespace gio = genogrove::io;

/**
 * Represents an exon in the unified graph
 * Exons are deduplicated: same coordinates + same strand = same node
 */
struct exon_data {
    std::string seqid;
    gdt::interval interval;
    char strand;  // '+', '-', or '.'
    std::unordered_set<std::string> transcript_ids;  // All transcripts containing this exon
    std::string attributes;  // GFF attributes string for metadata extraction

    // Unique identifier for deduplication
    std::string get_key() const {
        return seqid + ":" + std::to_string(interval.start) + "-" +
               std::to_string(interval.end) + ":" + strand;
    }

    bool operator==(const exon_data& other) const {
        return seqid == other.seqid &&
               interval.start == other.interval.start &&
               interval.end == other.interval.end &&
               strand == other.strand;
    }
};

/**
 * Represents a segment in the unified graph
 * A segment is a linear path through exons (consecutive in transcript space)
 */
struct segment_data {
    std::string segment_id;  // Unique identifier
    std::unordered_set<std::string> transcript_ids;  // Transcripts using this segment

    // Genomic extent (for quick queries)
    std::string seqid;
    gdt::coordinate start;
    gdt::coordinate end;
    char strand;

    // For merging identical segments (based on exon composition)
    std::vector<size_t> exon_node_ids;  // Node IDs of exons in this segment

    std::string get_key() const {
        // Key based on exon sequence for merging identical segments
        std::string key;
        for (size_t id : exon_node_ids) {
            key += std::to_string(id) + ",";
        }
        return key;
    }
};

/**
 * Node type discriminator
 */
enum class node_type {
    EXON,
    SEGMENT
};

/**
 * Unified node that can be either an exon or a segment
 */
struct transcript_node {
    node_type type;
    std::variant<exon_data, segment_data> data;

    // Helper constructors
    static transcript_node make_exon(const exon_data& exon) {
        return transcript_node{node_type::EXON, exon};
    }

    static transcript_node make_segment(const segment_data& segment) {
        return transcript_node{node_type::SEGMENT, segment};
    }

    // Accessors
    exon_data* as_exon() {
        return type == node_type::EXON ? &std::get<exon_data>(data) : nullptr;
    }

    segment_data* as_segment() {
        return type == node_type::SEGMENT ? &std::get<segment_data>(data) : nullptr;
    }

    const exon_data* as_exon() const {
        return type == node_type::EXON ? &std::get<exon_data>(data) : nullptr;
    }

    const segment_data* as_segment() const {
        return type == node_type::SEGMENT ? &std::get<segment_data>(data) : nullptr;
    }
};

/**
 * Edge type discriminator
 */
enum class edge_type {
    EXON_TO_EXON,           // Connection between exons within a transcript
    SEGMENT_TO_EXON,        // Segment pointing to its first exon
    SEGMENT_TO_SEGMENT_NORMAL,  // Linear continuation between segments
    SEGMENT_TO_SEGMENT_FUSION   // Fusion/breakpoint between segments
};

/**
 * Unified edge data
 */
struct transcript_edge {
    edge_type type;
    std::unordered_set<std::string> transcript_ids;  // Transcripts using this connection

    static transcript_edge make_exon_edge(const std::unordered_set<std::string>& transcripts) {
        return transcript_edge{edge_type::EXON_TO_EXON, transcripts};
    }

    static transcript_edge make_segment_to_exon() {
        return transcript_edge{edge_type::SEGMENT_TO_EXON, {}};
    }

    static transcript_edge make_segment_normal(const std::unordered_set<std::string>& transcripts) {
        return transcript_edge{edge_type::SEGMENT_TO_SEGMENT_NORMAL, transcripts};
    }

    static transcript_edge make_segment_fusion(const std::unordered_set<std::string>& transcripts) {
        return transcript_edge{edge_type::SEGMENT_TO_SEGMENT_FUSION, transcripts};
    }
};

/**
 * Unified graph structure for transcript analysis
 *
 * Contains both exon nodes and segment nodes in a single graph:
 * - Exon nodes: Deduplicated exons with transcript connections (exon → exon edges)
 * - Segment nodes: Linear exon paths (segment → exon, segment → segment edges)
 * - Segments link to their first exon via SEGMENT_TO_EXON edges
 * - Fusion events represented as SEGMENT_TO_SEGMENT_FUSION edges
 */
class transcript_graph {
public:
    transcript_graph();
    ~transcript_graph();

    /**
     * Build unified graph from GFF entries
     * Creates exon nodes, segment nodes, and all edges
     * @param entries Vector of GFF entries (should be sorted by transcript)
     */
    void build_from_gff_entries(const std::vector<gio::gff_entry>& entries);

    /**
     * Get the unified graph
     */
    ggs::graph<transcript_node, transcript_edge>* get_graph() const { return graph; }

    /**
     * Get node by index
     */
    const transcript_node* get_node(size_t idx) const;

    /**
     * Get all exon node indices
     */
    std::vector<size_t> get_exon_nodes() const;

    /**
     * Get all segment node indices
     */
    std::vector<size_t> get_segment_nodes() const;

private:
    ggs::graph<transcript_node, transcript_edge>* graph;

    // Track node indices for different types
    std::vector<size_t> exon_node_indices;
    std::vector<size_t> segment_node_indices;

    /**
     * Build exon nodes and edges from GFF entries
     * Deduplicates exons and creates EXON_TO_EXON edges
     */
    void build_exon_nodes(const std::vector<gio::gff_entry>& entries);

    /**
     * Build segment nodes and edges
     * Creates segments as linear exon paths, adds SEGMENT_TO_EXON and SEGMENT_TO_SEGMENT edges
     */
    void build_segment_nodes();

    /**
     * Determine if connection between two exons represents a fusion/breakpoint
     */
    bool is_fusion_point(const exon_data& from, const exon_data& to) const;

    /**
     * Helper: Group GFF entries by transcript ID
     */
    std::unordered_map<std::string, std::vector<gio::gff_entry>>
        group_by_transcript(const std::vector<gio::gff_entry>& entries);
};

#endif //ATROPLEX_TRANSCRIPT_GRAPH_HPP
