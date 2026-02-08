/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_INDEX_STATS_HPP
#define ATROPLEX_INDEX_STATS_HPP

#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

#include "genomic_feature.hpp"

/**
 * Statistics collected from a built grove index.
 * Provides a summary of the pan-transcriptome structure.
 */
struct index_stats {
    // Annotation sources
    std::vector<std::string> annotation_sources;

    // Basic counts
    size_t total_chromosomes = 0;
    size_t total_genes = 0;
    size_t total_transcripts = 0;
    size_t total_segments = 0;
    size_t total_exons = 0;
    size_t total_edges = 0;
    size_t total_samples = 0;       // entries with type == "sample" (excludes annotations)

    // Gene biotype breakdown
    std::map<std::string, size_t> genes_by_biotype;

    // Transcripts per gene
    double mean_transcripts_per_gene = 0;
    double median_transcripts_per_gene = 0;
    size_t max_transcripts_per_gene = 0;
    std::string max_transcripts_gene_id;
    size_t single_isoform_genes = 0;
    size_t multi_isoform_genes = 0;

    // Exons per segment
    double mean_exons_per_segment = 0;
    double median_exons_per_segment = 0;
    size_t max_exons_per_segment = 0;
    size_t single_exon_segments = 0;

    // Exon sharing
    size_t shared_exons = 0;          // exons used by >1 transcript
    double mean_transcripts_per_exon = 0;
    size_t max_transcripts_per_exon = 0;
    size_t constitutive_exons = 0;    // present in all transcripts of their gene
    size_t alternative_exons = 0;     // present in only some transcripts

    // Segment sharing distribution across samples
    size_t conserved_segments = 0;       // segments present in ALL samples
    size_t shared_segments = 0;          // segments in 2+ but not all samples
    size_t sample_specific_segments = 0; // segments in exactly 1 sample

    // Isoform diversity
    double isoform_diversity = 0;     // mean pairwise Jaccard distance of exon sets per gene (0=identical, 1=disjoint)
    double deduplication_ratio = 0;   // segments / transcripts (1.0 = no dedup, lower = more sharing)
    size_t multi_segment_genes = 0;   // genes with >1 distinct segment (included in diversity calculation)

    // Graph structure
    size_t branching_exons = 0;       // exons with >1 unique downstream exon

    // Splicing hubs (branching exons with >MIN_HUB_BRANCHES downstream targets)
    struct branching_exon_info {
        std::string gene_name;
        std::string gene_id;
        std::string exon_id;
        std::string coordinate;
        size_t branches = 0;          // total unique downstream exons
        size_t transcripts = 0;       // transcripts using this exon
        std::unordered_set<uint32_t> sample_idx;              // samples containing this exon
        std::unordered_map<uint32_t, size_t> sample_branches; // per-sample downstream branch count
    };
    std::vector<branching_exon_info> splicing_hubs;
    static constexpr size_t MIN_HUB_BRANCHES = 10;
    static constexpr size_t MAX_DISPLAY_HUBS = 20;

    // Per-chromosome summary
    struct chromosome_stats {
        size_t genes = 0;
        size_t segments = 0;
        size_t exons = 0;
    };
    std::map<std::string, chromosome_stats> per_chromosome;

    // Per-sample summary (keyed by sample registry ID)
    struct sample_stats {
        size_t segments = 0;            // segments this sample contributes to
        size_t exclusive_segments = 0;  // segments only in this sample
        size_t exons = 0;              // exons this sample contributes to
        size_t exclusive_exons = 0;    // exons only in this sample
        size_t genes = 0;             // genes this sample has features in
        size_t transcripts = 0;       // transcript count
        double isoform_diversity = 0; // mean pairwise Jaccard distance of exon sets per gene
        double deduplication_ratio = 0; // segments / transcripts
        double mean_expression = 0;   // mean expression across segments (if available)
        size_t expressed_segments = 0; // segments with expression data
    };
    std::unordered_map<uint32_t, sample_stats> per_sample;

    // Per-source summary (keyed by GFF column 2 value: HAVANA, ENSEMBL, TALON, etc.)
    struct source_stats {
        size_t segments = 0;            // segments with this source
        size_t exclusive_segments = 0;  // segments only from this source
        size_t exons = 0;              // exons with this source
        size_t exclusive_exons = 0;    // exons only from this source
        size_t genes = 0;             // genes that have features from this source
    };
    std::map<std::string, source_stats> per_source;

    /**
     * Collect statistics by traversing a built grove.
     * @param detailed  If true, run expensive Phase 4b (Jaccard diversity).
     */
    static index_stats collect(grove_type& grove, bool detailed = true);

    /**
     * Write full analysis report to a text file.
     * Includes all sections: per-sample, per-source, exon sharing, etc.
     */
    void write(const std::string& path) const;

    /**
     * Write compact summary (for --stats flag).
     * Includes: overview, biotype, per-chromosome. No per-sample analysis.
     */
    void write_summary(const std::string& path) const;

    /**
     * Write per-sample statistics as CSV.
     * Samples are columns, metrics are rows.
     * First column is the metric name, followed by one column per sample.
     */
    void write_sample_csv(const std::string& path) const;

    /**
     * Write per-source statistics as CSV.
     * Sources (GFF column 2) are columns, metrics are rows.
     */
    void write_source_csv(const std::string& path) const;

    /**
     * Write splicing hubs as TSV (exons with >MIN_HUB_BRANCHES downstream targets).
     */
    void write_splicing_hubs_tsv(const std::string& path) const;
};

#endif // ATROPLEX_INDEX_STATS_HPP