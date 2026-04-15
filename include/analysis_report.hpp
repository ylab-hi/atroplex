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
#include <fstream>
#include <map>
#include <memory>
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
 * these at output time, not tracked separately. Per-feature detail outputs
 * (conserved exons, splicing hubs, branch details) stream directly to disk
 * during the single traversal — never accumulated in memory.
 *
 * Per-gene accumulators live in a local `active_genes` map inside collect()
 * and are cleared at every chromosome boundary, so memory scales with
 * feature count, never with sample count.
 *
 * Usage:
 *   analysis_report report;
 *   report.begin_conserved_exon_stream("conserved.tsv");  // optional
 *   report.begin_splicing_hub_streams("hubs.tsv", "branches.tsv");  // optional
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

        // Phase 8.4 diversity — constitutive vs alternative exon classification
        size_t constitutive_exons = 0;   // exons used in every segment of their gene
        size_t alternative_exons = 0;    // exons used in some segments of their gene

        // Phase 8.4 diversity — per-gene entropy / effective isoforms.
        // Sums accumulated over multi-segment genes this sample participates in;
        // divide by multi_segment_genes at output time for the mean.
        double gene_exon_entropy_sum = 0;  // Σ H(exon usage) across multi-seg genes
        double effective_isoforms_sum = 0; // Σ 2^H(segment→tx dist)
        size_t multi_segment_genes = 0;

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
    // Note: the absorbed-segment count lives on build_summary::counters,
    // not here. By default tombstones stay physically in the grove and
    // analysis_report skips them via the `if (seg.absorbed) continue;`
    // filter in collect(). Under `--prune-tombstones` the tree no longer
    // contains them at all, so the filter is a no-op.

    // ── Streaming outputs (set before collect()) ────────────────────
    //
    // Per-feature detail files are opened up-front and rows are written
    // inline during the single grove traversal in collect(). Nothing is
    // accumulated in memory beyond the flat counters above.
    std::unique_ptr<std::ofstream> conserved_exon_stream;
    std::vector<uint32_t> conserved_stream_sample_ids;  // sample IDs to emit columns for
    std::vector<bool> conserved_stream_is_sample;       // parallel: true if type=="sample"

    /**
     * Open the conserved-exons TSV, write its header, and arm inline
     * streaming inside the next collect() call. Exons whose sample_count
     * equals the number of non-replicate sample/annotation entries are
     * written as soon as they are first visited during traversal.
     *
     * Must be called after sample_registry is populated and before
     * collect(). Each row carries per-sample expression columns for
     * entries with type=="sample".
     */
    void begin_conserved_exon_stream(const std::string& path);

    // ── Phase 8.3: splicing hub streams ─────────────────────────────
    //
    // Hub exons (>MIN_HUB_BRANCHES unique downstream exons) are detected
    // on the exon chain walk and buffered per-gene in gene_acc; rows are
    // written at gene finalization. Temporary per-sample vectors are
    // reused across hubs so peak memory stays at O(num_samples), not
    // O(num_hubs × num_samples).
    std::unique_ptr<std::ofstream> hub_stream;
    std::unique_ptr<std::ofstream> branch_stream;
    std::vector<uint32_t> hub_stream_sample_ids;   // parallel labels
    std::vector<bool> hub_stream_is_sample;
    static constexpr size_t MIN_HUB_BRANCHES = 10;

    /**
     * Open the splicing hubs TSV and branch details TSV, write their
     * headers, and arm inline streaming. Hub rows are emitted at gene
     * finalization; branch rows are emitted per (hub × target) pair.
     * Must be called after sample_registry is populated and before
     * collect().
     */
    void begin_splicing_hub_streams(const std::string& hubs_path,
                                    const std::string& branches_path);

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

    /**
     * Write exon sharing TSV: metrics as rows, samples as columns.
     * Columns: metric, total, <sample1>, <sample2>, ...
     * Rows: total, exclusive, shared, conserved.
     * Formats existing per_sample counters — no extra traversal.
     */
    void write_exon_sharing(const std::string& path) const;

    /**
     * Write segment sharing TSV: same layout as exon_sharing.
     * Rows: total, exclusive, shared, conserved.
     */
    void write_segment_sharing(const std::string& path) const;
};

#endif // ATROPLEX_ANALYSIS_REPORT_HPP
