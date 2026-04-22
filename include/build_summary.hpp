/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_BUILD_SUMMARY_HPP
#define ATROPLEX_BUILD_SUMMARY_HPP

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "genomic_feature.hpp"

/**
 * Lightweight counters tracking transcript disposition during grove build.
 * A single instance is threaded through the build path (build_gff / build_bam
 * / segment_builder) via a pointer so that filter and absorption sites can
 * increment categories without restructuring return types. All counts are
 * later folded into the final build_summary.
 *
 * Categories:
 *   input_transcripts  — raw transcripts observed across all input files
 *                        (counted AFTER the scaffold filter, so this reflects
 *                        only transcripts on retained seqids)
 *   merged_transcripts — transcripts folded into an existing segment without
 *                        creating a new one (Rule 0 FSM, Rule 5 terminal
 *                        variant, forward Rules 1-4 ISM subsequence,
 *                        fuzzy-FSM subsequence match)
 *   absorbed_segments  — existing segments tombstoned by reverse absorption
 *                        and physically removed by the tombstone sweep
 *   discarded_transcripts — transcripts dropped entirely (expression filter,
 *                        mono-exon Rules 6/8, Rules 3/4 fragment drops against
 *                        reference parents)
 *   scaffold_filtered_transcripts — transcripts skipped at ingest because
 *                        their seqid is not a canonical main chromosome
 *                        (chr1..chr22, chrX, chrY, chrM). Disabled with
 *                        --include-scaffolds. Counted before any other
 *                        processing so scaffold features never enter the grove.
 */
struct build_counters {
    size_t input_transcripts = 0;
    size_t merged_transcripts = 0;
    size_t absorbed_segments = 0;
    size_t discarded_transcripts = 0;
    size_t scaffold_filtered_transcripts = 0;
};

/**
 * Build-phase summary — the data backing `.ggx.summary`.
 *
 * Contains only information derivable from the build caches (segment_caches,
 * exon_caches, gene_indices) plus the counters accumulated during the build
 * itself. No tree traversal is required: the builder populates this directly
 * from the caches it already holds in memory.
 *
 * For per-sample / sharing / splicing-hub analysis run `atroplex inspect`,
 * which uses `analysis_report` instead.
 */
struct build_summary {
    // ── Inputs ─────────────────────────────────────────────────────
    std::vector<std::string> annotation_sources;   // registry IDs of all entries
    size_t total_entries = 0;                      // excluding replicates

    // ── Basic counts ───────────────────────────────────────────────
    size_t total_chromosomes = 0;
    size_t total_genes = 0;
    size_t total_transcripts = 0;
    size_t total_segments = 0;
    size_t total_exons = 0;
    size_t total_edges = 0;

    // ── Biotypes ───────────────────────────────────────────────────
    std::map<std::string, size_t> genes_by_biotype;
    std::map<std::string, size_t> transcripts_by_biotype;

    // ── Transcripts per gene ──────────────────────────────────────
    double mean_transcripts_per_gene = 0;
    double median_transcripts_per_gene = 0;
    size_t max_transcripts_per_gene = 0;
    std::string max_transcripts_gene_id;
    size_t single_isoform_genes = 0;
    size_t multi_isoform_genes = 0;

    // ── Exons per segment ─────────────────────────────────────────
    double mean_exons_per_segment = 0;
    double median_exons_per_segment = 0;
    size_t max_exons_per_segment = 0;
    size_t single_exon_segments = 0;

    // ── Build processing counters (from build_counters) ───────────
    build_counters counters;

    // ── B+ tree structure ─────────────────────────────────────────
    int tree_order = 0;
    std::map<std::string, int> tree_depth_per_chromosome;

    // ── Per-chromosome summary ────────────────────────────────────
    struct chromosome_stats {
        size_t genes = 0;
        size_t segments = 0;
        size_t exons = 0;
    };
    std::map<std::string, chromosome_stats> per_chromosome;

    // ── Build timing ──────────────────────────────────────────────
    double build_time_seconds = 0;

    /**
     * Populate this build_summary from builder caches post-sweep.
     * Caller must have already run tombstone cleanup and replicate merging.
     * Mirrors analysis_report::collect — both are member mutators so callers
     * can configure instance-level state before collection.
     *
     * **Precondition:** the instance must be default-constructed (or
     * otherwise empty). collect() appends to its internal containers without
     * clearing them first, so calling it twice on the same instance would
     * double-count annotation_sources, biotype maps, per_chromosome entries,
     * and the distribution aggregates. Use one build_summary per build.
     *
     * @param grove          Grove reference (used for edge count + tree depth)
     * @param segment_caches Chromosome -> segment cache (post-sweep)
     * @param exon_caches    Chromosome -> exon cache
     * @param segment_count  Running segment counter (pre-sweep total)
     * @param counters_in    Build counters accumulated during build
     */
    void collect(
        grove_type& grove,
        const chromosome_segment_caches& segment_caches,
        const chromosome_exon_caches& exon_caches,
        size_t segment_count,
        const build_counters& counters_in
    );

    /**
     * Write `.ggx.summary` text report.
     */
    void write_summary(const std::string& path) const;
};

#endif // ATROPLEX_BUILD_SUMMARY_HPP
