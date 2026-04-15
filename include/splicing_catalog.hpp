/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SPLICING_CATALOG_HPP
#define ATROPLEX_SPLICING_CATALOG_HPP

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <cstdint>

#include "genomic_feature.hpp"
#include "sample_info.hpp"

/// Alternative splicing event types
enum class splicing_event_type {
    CASSETTE_EXON,       ///< Exon included in some paths, skipped in others
    MUTUALLY_EXCLUSIVE,  ///< Two exons never co-occur in the same path
    ALT_5PRIME,          ///< Same acceptor (exon start), different donors (exon end)
    ALT_3PRIME,          ///< Same donor (exon end), different acceptors (exon start)
    INTRON_RETENTION,    ///< Intron retained — exon spans intron boundary
    ALT_FIRST_EXON,      ///< Different first exons across transcript paths
    ALT_LAST_EXON        ///< Different last exons across transcript paths
};

/// A single alternative splicing event detected in the exon graph
struct splicing_event {
    splicing_event_type type;
    std::string gene_id;
    std::string gene_name;
    std::string chromosome;

    // Coordinates of the event
    key_ptr upstream_exon = nullptr;    ///< Exon upstream of the event (branch point)
    key_ptr downstream_exon = nullptr;  ///< Exon downstream (rejoin point, if applicable)
    std::vector<key_ptr> affected_exons; ///< Exon(s) involved in the event (cassette, MXE targets, etc.)

    // Per-sample inclusion: sample_id → PSI (percent spliced in)
    std::unordered_map<uint32_t, float> sample_psi;

    // Per-sample transcript counts at this event
    std::unordered_map<uint32_t, size_t> sample_included;  ///< Transcripts including the event
    std::unordered_map<uint32_t, size_t> sample_total;     ///< Total transcripts at the locus
};

/**
 * Splicing event catalog — detects and classifies alternative splicing events
 * from the exon graph structure.
 *
 * Walks the per-gene exon graph, identifies branch points, and classifies
 * each event using topology analysis:
 * - Compare exon chains across all segments (transcript paths) within a gene
 * - At each branch point, determine the event type from the graph topology
 * - Track per-sample PSI for each event
 */
class splicing_catalog {
public:
    /**
     * Build the catalog from builder caches (available during build).
     * Uses pre-populated gene_segment_index for efficient per-gene traversal.
     */
    static std::vector<splicing_event> collect(
        const chromosome_gene_segment_indices& gene_indices,
        grove_type& grove
    );

    /**
     * Build the catalog by walking the grove directly (no builder caches needed).
     * Reconstructs per-gene segment chains from the grove's B+ tree and graph edges.
     * Use this when loading from a .ggx file where builder caches are unavailable.
     */
    static std::vector<splicing_event> collect_from_grove(grove_type& grove);

    /**
     * Write the splicing event catalog to TSV.
     * Columns: gene_id, gene_name, event_type, chromosome, upstream_exon,
     *          downstream_exon, affected_exon(s), per-sample PSI columns.
     */
    static void write_events_tsv(
        const std::string& filepath,
        const std::vector<splicing_event>& events
    );

    /**
     * Write a summary of event counts by type.
     */
    static void write_summary(
        const std::string& filepath,
        const std::vector<splicing_event>& events
    );

    /// Detect splicing events for a single gene by comparing all segment exon chains.
    /// Exposed publicly so streaming callers (analysis_report::collect) can run
    /// event detection per-gene at chromosome-bounded gene finalization and write
    /// rows inline without holding the full catalog in memory.
    static std::vector<splicing_event> detect_gene_events(
        const std::string& gene_id,
        const std::string& chromosome,
        const std::vector<segment_chain_entry>& segments,
        grove_type& grove
    );

    /// Convert a splicing_event_type to its TSV-friendly string form.
    /// Exposed publicly for streaming row writers.
    static std::string event_type_str(splicing_event_type type);

    /// Format an exon key_ptr as its coordinate string for TSV output.
    /// Exposed publicly for streaming row writers.
    static std::string format_exon(key_ptr exon);

private:

    /// Detect cassette exon events: exons present in some paths, absent in others
    static void detect_cassette_exons(
        const std::string& gene_id,
        const std::string& chromosome,
        const std::vector<segment_chain_entry>& segments,
        std::vector<splicing_event>& events
    );

    /// Detect alternative 5'/3' splice site events
    static void detect_alt_splice_sites(
        const std::string& gene_id,
        const std::string& chromosome,
        const std::vector<segment_chain_entry>& segments,
        std::vector<splicing_event>& events
    );

    /// Detect intron retention: one segment has two adjacent exons defining an intron,
    /// another segment has a single exon spanning that intron region
    static void detect_intron_retention(
        const std::string& gene_id,
        const std::string& chromosome,
        const std::vector<segment_chain_entry>& segments,
        std::vector<splicing_event>& events
    );

    /// Detect alternative first/last exon events
    static void detect_alt_terminal_exons(
        const std::string& gene_id,
        const std::string& chromosome,
        const std::vector<segment_chain_entry>& segments,
        std::vector<splicing_event>& events
    );

    /// Detect mutually exclusive exon events
    static void detect_mutually_exclusive(
        const std::string& gene_id,
        const std::string& chromosome,
        const std::vector<segment_chain_entry>& segments,
        std::vector<splicing_event>& events
    );

    /// Compute PSI for an event based on segment sample membership and expression
    static void compute_psi(
        splicing_event& event,
        const std::vector<segment_chain_entry>& segments
    );
};

#endif //ATROPLEX_SPLICING_CATALOG_HPP
