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
 * Walks the B+ tree once, processing segments and their exon chains inline.
 * No intermediate maps stored — only flat counters accumulated.
 *
 * Usage:
 *   analysis_report report;
 *   report.collect(grove);
 *   report.write_overview("output.tsv");
 */
struct analysis_report {

    // ── Global counters ─────────────────────────────────────────────
    size_t total_chromosomes = 0;
    size_t total_genes = 0;
    size_t total_transcripts = 0;
    size_t total_segments = 0;
    size_t total_exons = 0;
    size_t total_edges = 0;
    size_t absorbed_segments = 0;
    size_t single_exon_segments = 0;

    // Transcripts per gene distribution
    double mean_transcripts_per_gene = 0;
    double median_transcripts_per_gene = 0;
    size_t max_transcripts_per_gene = 0;
    std::string max_transcripts_gene_id;
    size_t single_isoform_genes = 0;
    size_t multi_isoform_genes = 0;

    // Exons per segment distribution
    double mean_exons_per_segment = 0;
    double median_exons_per_segment = 0;
    size_t max_exons_per_segment = 0;

    // Biotype distributions
    std::map<std::string, size_t> genes_by_biotype;
    std::map<std::string, size_t> transcripts_by_biotype;

    // Per-chromosome
    struct chromosome_stats {
        size_t genes = 0;
        size_t segments = 0;
        size_t exons = 0;
    };
    std::map<std::string, chromosome_stats> per_chromosome;

    // Dedup ratio
    double deduplication_ratio = 0;

    // ── Per-sample counters (flat vector indexed by sample_id) ───────

    struct sample_counters {
        size_t genes = 0;
        size_t transcripts = 0;
        size_t segments = 0;
        size_t exons = 0;
        size_t single_exon_segments = 0;
        double deduplication_ratio = 0;
        std::map<std::string, size_t> genes_by_biotype;
        std::map<std::string, size_t> transcripts_by_biotype;
    };
    std::vector<sample_counters> per_sample;  // indexed by sample_id

    // ── Collection ──────────────────────────────────────────────────

    /**
     * Collect basic index statistics by walking the grove once.
     * Phase 1: segment + exon counting, per-chromosome, biotypes,
     * per-sample counters.
     */
    void collect(grove_type& grove);

    // ── Output ──────────────────────────────────────────────────────

    void write_overview(const std::string& path) const;
    void write_per_chromosome(const std::string& path) const;
    void write_biotypes(const std::string& path) const;
};

#endif // ATROPLEX_ANALYSIS_REPORT_HPP
