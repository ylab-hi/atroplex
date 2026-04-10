/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_ANALYSIS_REPORT_HPP
#define ATROPLEX_ANALYSIS_REPORT_HPP

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "genomic_feature.hpp"

/**
 * Streaming analysis report — single-pass grove traversal.
 *
 * Data is stored at three levels:
 *   - per_chromosome: genes, segments, exons per chromosome
 *   - per_sample:     genes, segments, exons, transcripts per sample
 *   - distributions:  transcripts_per_gene, exons_per_segment (for median/max)
 *
 * Global totals are derived from these at output time — not tracked separately.
 *
 * Usage:
 *   analysis_report report;
 *   report.collect(grove);
 *   report.write_overview("output.tsv");
 */
struct analysis_report {

    // ── Per-sample (flat vector indexed by sample_id) ───────────────
    struct sample_counters {
        // Gene / transcript
        size_t genes = 0;
        size_t transcripts = 0;

        // Segments
        size_t segments = 0;
        size_t exclusive_segments = 0;
        size_t shared_segments = 0;
        size_t conserved_segments = 0;
        size_t single_exon_segments = 0;

        // Exons
        size_t exons = 0;
        size_t exclusive_exons = 0;
        size_t shared_exons = 0;
        size_t conserved_exons = 0;

        // Expression
        double expression_sum = 0;
        size_t expressed_segments = 0;

        // Biotypes
        std::map<std::string, size_t> genes_by_biotype;
        std::map<std::string, size_t> transcripts_by_biotype;
    };
    std::vector<sample_counters> per_sample;  // indexed by sample_id

    // ── Distributions (for median/max, bounded by feature count) ────
    std::vector<size_t> transcripts_per_gene;
    std::vector<size_t> exons_per_segment;

    // ── Counters not derivable from per-sample ────────────────────────
    size_t total_transcripts = 0;  // sum of transcript_ids.size() across all segments
    size_t total_exons = 0;        // unique exons (pointer-deduplicated)
    size_t absorbed_segments = 0;
    size_t total_edges = 0;

    // ── Collection ──────────────────────────────────────────────────

    void collect(grove_type& grove);

    // ── Output ──────────────────────────────────────────────────────

    void write_overview(const std::string& path) const;
    void write_per_sample(const std::string& path) const;
};

#endif // ATROPLEX_ANALYSIS_REPORT_HPP