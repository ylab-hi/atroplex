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

    // Transcript biotype breakdown
    std::map<std::string, size_t> transcripts_by_biotype;

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

    // ISM absorption
    size_t absorbed_segments = 0;    // ISM segments absorbed into longer parents

    // Exon sharing
    size_t shared_exons = 0;          // exons used by >1 transcript
    double mean_transcripts_per_exon = 0;
    size_t max_transcripts_per_exon = 0;
    size_t constitutive_exons = 0;    // present in all transcripts of their gene
    size_t alternative_exons = 0;     // present in only some transcripts

    // Conserved exon detail (exons present in ALL samples)
    struct conserved_exon_entry {
        std::string exon_id;
        std::string gene_name;
        std::string gene_id;
        std::string chromosome;
        std::string coordinate;
        size_t n_transcripts = 0;     // total transcripts using this exon
        bool constitutive = false;    // present in all transcripts of its gene
        std::unordered_map<uint32_t, size_t> sample_transcripts;  // per-sample transcript count
        std::unordered_map<uint32_t, float> sample_expression;    // per-sample expression (if available)
    };
    std::vector<conserved_exon_entry> conserved_exon_details;

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
        size_t exon_number = 0;       // 1-based position in representative transcript chain
        size_t total_exons = 0;       // total exons in that representative segment
        std::unordered_set<uint32_t> sample_idx;              // samples containing this exon
        std::unordered_map<uint32_t, size_t> sample_branches;    // per-sample downstream branch count
        std::unordered_map<uint32_t, size_t> sample_shared;      // branches also in ≥1 other sample
        std::unordered_map<uint32_t, size_t> sample_unique;      // branches only in this sample
        std::unordered_map<uint32_t, size_t> sample_transcripts; // per-sample transcript count

        std::unordered_map<uint32_t, double> sample_psi;          // traditional PSI: hub tx / gene tx per sample
        std::unordered_map<uint32_t, double> sample_entropy;     // Shannon entropy of branch usage
        std::unordered_map<uint32_t, float> sample_expression;   // per-sample expression at hub exon

        // Individual branch targets (for detail file)
        struct branch_target {
            std::string exon_id;
            std::string coordinate;
            std::unordered_set<uint32_t> sample_idx;
            std::unordered_map<uint32_t, size_t> sample_transcripts; // per-sample tx count through this target
            std::unordered_map<uint32_t, float> sample_expression;   // per-sample expression at target exon
        };
        std::vector<branch_target> targets;
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
        size_t shared_segments = 0;     // segments in 2+ but not all samples
        size_t conserved_segments = 0;  // segments in ALL samples
        size_t exons = 0;              // exons this sample contributes to
        size_t exclusive_exons = 0;    // exons only in this sample
        size_t shared_exons = 0;       // exons in 2+ but not all samples
        size_t conserved_exons = 0;    // exons in ALL samples
        size_t constitutive_exons = 0; // exons in all transcripts of their gene
        size_t alternative_exons = 0;  // exons in some transcripts of their gene
        size_t genes = 0;             // genes this sample has features in
        size_t transcripts = 0;       // transcript count
        double isoform_diversity = 0; // mean pairwise Jaccard distance of exon sets per gene
        double deduplication_ratio = 0; // segments / transcripts
        double mean_expression = 0;   // mean expression across segments (if available)
        size_t expressed_segments = 0; // segments with expression data
        std::map<std::string, size_t> genes_by_biotype;        // per-sample gene biotype
        std::map<std::string, size_t> transcripts_by_biotype;  // per-sample transcript biotype
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
     * Options for collect() — controls detail level and streaming output.
     */
    struct collect_options {
        bool detailed = true;       // run expensive Phase 4b (Jaccard diversity)
        std::string output_dir;     // if set, stream large outputs to files during traversal
        std::string basename;       // file prefix for streamed outputs
    };

    /**
     * Collect statistics by traversing a built grove.
     * When output_dir is set, conserved exons and branch details are streamed
     * directly to files to reduce peak memory usage.
     */
    static index_stats collect(grove_type& grove, const collect_options& opts);
    static index_stats collect(grove_type& grove) {
        return collect(grove, collect_options{});
    }

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

    /**
     * Write per-branch detail file (one row per hub × downstream target).
     * Shows which samples each individual branch target occurs in.
     */
    void write_branch_details_tsv(const std::string& path) const;

    /**
     * Write exon sharing statistics as TSV.
     * Metrics as rows, samples as columns. Covers cross-sample and transcript-level sharing.
     */
    void write_exon_sharing_tsv(const std::string& path) const;

    /**
     * Write segment sharing statistics as TSV.
     * Metrics as rows, samples as columns. Covers cross-sample sharing distribution.
     */
    void write_segment_sharing_tsv(const std::string& path) const;

    /**
     * Write conserved exons (present in ALL samples) as TSV.
     * One row per exon with per-sample transcript counts.
     */
    void write_conserved_exons_tsv(const std::string& path) const;

    /**
     * Write global overview statistics as TSV (metric/value pairs).
     */
    void write_overview_tsv(const std::string& path) const;

    /**
     * Write per-chromosome summary as TSV.
     * One row per chromosome with genes, segments, exons.
     */
    void write_per_chromosome_tsv(const std::string& path) const;

    /**
     * Write per-source summary as TSV.
     * One row per annotation source (GFF column 2) with genes, segments, exons.
     */
    void write_per_source_tsv(const std::string& path) const;

    /**
     * Write per-sample summary as TSV.
     * One row per sample with genes, segments, exons, diversity, etc.
     */
    void write_per_sample_tsv(const std::string& path) const;

    /**
     * Write gene and transcript biotype breakdown as TSV.
     * One row per (level, biotype) pair with per-sample counts.
     * Level is "gene" or "transcript".
     */
    void write_biotype_tsv(const std::string& path) const;
};

#endif // ATROPLEX_INDEX_STATS_HPP