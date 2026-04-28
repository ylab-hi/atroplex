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
#include "quant_sidecar.hpp"

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

        // Quantification (populated only when collect() receives a non-null
        // qtx_reader). expression_sum is the sum of per-segment values for
        // every (segment, sample) record in the sidecar; expressed_segments
        // is the count of distinct segments with a value for this sample.
        // mean_expression at output time is expression_sum / expressed_segments.
        double expression_sum = 0;
        size_t expressed_segments = 0;

        // Biotypes
        std::map<std::string, size_t> genes_by_biotype;
        std::map<std::string, size_t> transcripts_by_biotype;
    };
    std::vector<sample_counters> per_sample;  // indexed by sample_id

    // ── Per-source (one row per GFF column-2 source string) ────────
    //
    // Bounded by source_registry capacity (16 sources max). Each
    // entry counts segments/exons that carry the source's bit, plus
    // the "exclusive" subset where ONLY this source contributed
    // (source_count() == 1). genes is the count of distinct genes
    // that have any feature carrying this source — accumulated at
    // gene finalization from gene_acc::sources_seen_in_gene.
    struct source_stats {
        size_t segments = 0;
        size_t exclusive_segments = 0;
        size_t exons = 0;
        size_t exclusive_exons = 0;
        size_t genes = 0;
    };
    std::map<std::string, source_stats> per_source;

    // ── Distributions (for median/max, bounded by feature count) ────
    std::vector<size_t> transcripts_per_gene;
    std::vector<size_t> exons_per_segment;

    // ── Global biotype breakdown (matches sample_counters versions) ──
    // Top-level distinct gene / transcript counts per biotype, used for
    // the `total` column of biotype.tsv. Each transcript_id appears in
    // exactly one segment's transcript_biotypes after dedup, so simple
    // bumps during the segment loop give the global count without any
    // per-feature deduplication state.
    std::map<std::string, size_t> genes_by_biotype;
    std::map<std::string, size_t> transcripts_by_biotype;

    // ── Counters not derivable from per-sample ────────────────────────
    size_t total_transcripts = 0;  // sum of transcript_ids.size() across all segments
    size_t total_exons = 0;        // unique exons (pointer-deduplicated)
    size_t total_edges = 0;            // grove-level total (unfiltered)
    size_t segment_to_exon_edges = 0;   // one per segment with an exon chain
    size_t exon_to_exon_edges = 0;      // exon chain links
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
    std::unique_ptr<std::ofstream> conserved_segment_stream;
    std::vector<uint32_t> conserved_stream_sample_ids;  // sample IDs to emit columns for
    std::vector<bool> conserved_stream_is_sample;       // parallel: true if type=="sample"
    bool conserved_emit_expression = false;             // emit per-sample expression columns when qtx reader available

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
    void begin_conserved_exon_stream(const std::string& path,
                                     bool emit_expression_columns = false);

    /**
     * Open the conserved-segments TSV, write its header, and arm inline
     * streaming inside the next collect() call. Segments whose sample_count
     * equals the number of non-replicate sample/annotation entries are
     * written when first visited during the grove traversal. Reuses the
     * sample-id vectors prepared by begin_conserved_exon_stream() if it
     * was called first; otherwise builds them now.
     *
     * Each row carries gene info, span coordinate, exon count, transcript
     * count, sample count, sources, biotype, and (if a .qtx reader is
     * attached) per-sample expression columns for entries with
     * type=="sample".
     */
    void begin_conserved_segment_stream(const std::string& path,
                                        bool emit_expression_columns = false);

    /**
     * Set the conservation threshold as a fraction in (0, 1] of sample-typed
     * registry entries. The minimum sample count required for a segment or
     * exon to be classified as conserved is ceil(total_sample_entries × frac).
     * Default is 1.0 (strict: present in every sample). Must be called
     * before collect(). Annotations don't count toward the denominator.
     */
    void set_conserved_fraction(double fraction);

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
    bool hub_emit_expression = false;              // emit per-sample expression columns when qtx reader available
    static constexpr size_t MIN_HUB_BRANCHES = 10;

    /**
     * Open the splicing hubs TSV and branch details TSV, write their
     * headers, and arm inline streaming. Hub rows are emitted at gene
     * finalization; branch rows are emitted per (hub × target) pair.
     * Must be called after sample_registry is populated and before
     * collect().
     */
    void begin_splicing_hub_streams(const std::string& hubs_path,
                                    const std::string& branches_path,
                                    bool emit_expression_columns = false);

    // ── Phase 8.6: splicing event stream ────────────────────────────
    //
    // At gene finalization we have the ordered exon chains captured
    // during the segment walk; the existing `splicing_catalog`
    // detectors classify cassette / alt_5' / alt_3' / intron retention
    // / alt-terminal / mutually-exclusive events from that. Rows are
    // written inline from inside `finalize_gene`, and the per-gene
    // `segments_in_gene` vector dies with the gene at chromosome
    // boundary via `active_genes.clear()`.
    std::unique_ptr<std::ofstream> splicing_events_stream;
    std::vector<uint32_t> splicing_events_sample_ids;

    /**
     * Open the splicing events TSV, write its header, and arm inline
     * streaming. Must be called after sample_registry is populated and
     * before collect().
     */
    void begin_splicing_events_stream(const std::string& path);

    // ── Collection ──────────────────────────────────────────────────

    /**
     * Traverse the grove and populate per-sample counters + distributions.
     *
     * **Precondition:** the instance must be default-constructed (or empty).
     * collect() appends to its internal per-sample counters and distribution
     * vectors without clearing them first, so calling it twice on the same
     * instance would double-count every metric. Use one analysis_report per
     * inspect invocation.
     */
    void collect(grove_type& grove,
                 quant_sidecar::Reader* qtx_reader = nullptr,
                 size_t min_samples = 0);

    // ── Output ──────────────────────────────────────────────────────

    void write_overview(const std::string& path) const;
    void write_per_sample(const std::string& path) const;

    /**
     * Write per-source TSV: one row per GFF column-2 source.
     * Columns: source, genes, segments, exclusive_segments, exons, exclusive_exons.
     * Bounded by source_registry capacity (16 entries max).
     */
    void write_per_source(const std::string& path) const;

    /**
     * Write biotype TSV: long-form gene + transcript biotype counts.
     * Columns: level, biotype, total, <sample1>, <sample2>, ...
     * Two sections: `level == "gene"` rows then `level == "transcript"`
     * rows, each sorted by global count descending. Empty file if there
     * are no biotypes (e.g., TALON-only builds with no biotype attrs).
     */
    void write_biotype(const std::string& path) const;

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

private:
    // ── Internal types (moved from collect() for method extraction) ──

    struct pending_hub {
        key_ptr exon = nullptr;
        std::vector<key_ptr> targets;
        size_t chain_pos = 0;
        size_t chain_total = 0;
    };

    struct gene_acc {
        uint32_t gene_idx = 0;
        std::string biotype;
        size_t segment_count = 0;
        sample_bitset sample_bits;
        std::vector<size_t> sample_tx;
        std::vector<size_t> segment_tx_counts;
        std::unordered_map<key_ptr, size_t> exon_seg_counts;
        std::vector<key_ptr> first_visit_exons;
        std::vector<pending_hub> pending_hubs;
        size_t max_hub_targets_in_gene = 0;
        std::vector<segment_chain_entry> segments_in_gene;
        uint16_t sources_seen_in_gene = 0;
    };

    using exon_expr_map_t = std::unordered_map<key_ptr,
                                               std::unordered_map<uint32_t, float>>;

    // ── Per-chromosome state (cleared at chromosome boundary) ───────
    std::unordered_map<uint32_t, gene_acc> active_genes_;
    exon_expr_map_t exon_expr_sum_;

    // ── Hub tally buffers (allocated once, reused across hubs) ───────
    std::vector<size_t> hub_branches_;
    std::vector<size_t> hub_shared_;
    std::vector<size_t> hub_unique_;
    std::vector<size_t> hub_psi_num_;
    std::vector<size_t> hub_psi_den_;
    std::vector<size_t> branch_counts_;

    // ── Traversal state ─────────────────────────────────────────────
    std::unordered_set<const void*> visited_exons_;
    size_t num_samples_ = 0;
    size_t total_samples_for_conserved_ = 0;
    size_t min_required_for_conserved_ = 0;  // ceil(total * conserved_fraction_); set in collect()
    double conserved_fraction_ = 1.0;        // (0,1]; default 1.0 = strict "in every sample"
    quant_sidecar::Reader* qtx_reader_ = nullptr;

    // ── Extracted helpers from collect() ─────────────────────────────

    void finalize_gene_stats(gene_acc& acc, const std::string& seqid);
    void finalize_gene_diversity(gene_acc& acc);
    void emit_hub_rows(gene_acc& acc, const std::string& seqid);
    void emit_event_rows(gene_acc& acc, const std::string& seqid,
                         grove_type& grove);

    void accumulate_segment_stats(const segment_feature& seg);
    std::vector<quant_sidecar::Reader::ValueRecord> lookup_segment_expression(
        const segment_feature& seg);
    void accumulate_gene(const segment_feature& seg, gene_acc& acc);
    void process_exon_visit(key_ptr exon_key, const segment_feature& seg,
                            const std::string& seqid, gene_acc& acc,
                            const std::vector<quant_sidecar::Reader::ValueRecord>& expr_records,
                            size_t chain_pos, size_t chain_total,
                            grove_type& grove);

    void stream_conserved_segment_row(
        const segment_feature& seg,
        const std::string& seqid,
        key_ptr seg_key,
        const std::vector<quant_sidecar::Reader::ValueRecord>& seg_expr_records);
};

#endif // ATROPLEX_ANALYSIS_REPORT_HPP
