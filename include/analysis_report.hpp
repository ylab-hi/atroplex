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
 * Data is stored at two levels:
 *   - per_sample:    flat vector<sample_counters> indexed by sample_id
 *                    (genes, segments, exons, transcripts, diversity)
 *   - distributions: transcripts_per_gene, exons_per_segment
 *                    (for median/max, one entry per feature — bounded)
 *
 * Global totals (total_genes, total_segments, total_exons) are derived from
 * these at output time, not tracked separately.
 *
 * Per-gene accumulators live in a local `active_genes` map inside collect()
 * and are cleared at every chromosome boundary, so memory scales with
 * feature count, never with sample count.
 *
 * Usage:
 *   analysis_report report;
 *   report.collect(grove);
 *   report.write_overview("overview.tsv");
 *   report.write_per_sample("per_sample.tsv");
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
    size_t total_edges = 0;
    // Note: the absorbed-segment count is part of build_summary::counters —
    // tombstones are physically removed by builder::remove_tombstones before
    // analyze runs, so analysis_report can never observe one.

    // ── Collection ──────────────────────────────────────────────────

    /**
     * Traverse the grove and populate per-sample counters + distributions.
     *
     * **Precondition:** the instance must be default-constructed (or empty).
     * collect() appends to its internal per-sample counters and distribution
     * vectors without clearing them first, so calling it twice on the same
     * instance would double-count every metric. Use one analysis_report per
     * analyze invocation.
     */
    void collect(grove_type& grove);

    // ── Output ──────────────────────────────────────────────────────

    void write_overview(const std::string& path) const;
    void write_per_sample(const std::string& path) const;
};

#endif // ATROPLEX_ANALYSIS_REPORT_HPP
