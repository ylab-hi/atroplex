/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_BUILD_GFF_HPP
#define ATROPLEX_BUILD_GFF_HPP

// standard
#include <string>
#include <optional>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <mutex>

// genogrove
#include <genogrove/io/gff_reader.hpp>

// class
#include "build_summary.hpp"
#include "genomic_feature.hpp"
#include "quant_sidecar.hpp"
#include "sample_info.hpp"
#include "segment_builder.hpp"

namespace gio = genogrove::io;

/**
 * GFF/GTF-specific grove builder
 * Handles construction of grove structures from GFF/GTF annotation files
 *
 * Processing strategy:
 * - Gene-by-gene processing for memory efficiency
 * - Creates exon features with CDS/UTR overlaps
 * - Creates segment features representing transcript paths
 * - Links exons within segments via graph edges
 * - Tags features with sample_id for pan-transcriptome tracking
 */
class build_gff {
public:
    /**
     * Build grove from single GFF/GTF file
     * @param grove Grove to populate
     * @param filepath Path to GFF/GTF file
     * @param sample_id Sample registry ID for pan-transcriptome tracking
     * @param exon_caches Chromosome-level exon caches (for cross-file deduplication)
     * @param segment_caches Chromosome-level segment caches (for cross-file deduplication)
     * @param segment_count Reference to segment counter (incremented for new segments)
     * @param num_threads Number of worker threads (0 = auto-detect)
     */
    static void build(grove_type& grove,
                      const std::filesystem::path& filepath,
                      std::optional<uint32_t> sample_id,
                      chromosome_exon_caches& exon_caches,
                      chromosome_segment_caches& segment_caches,
                      size_t& segment_count,
                      const build_options& opts,
                      build_counters& counters,
                      quant_sidecar::SampleStreamWriter* sidecar_writer = nullptr);

    /**
     * Parse GTF/GFF header to extract sample metadata
     * Looks for ##property: value lines (GFF comment convention)
     *
     * Supported header properties (case-insensitive):
     *   ##id: sample_001
     *   ##description: Evidence-based annotation of the human genome
     *   ##assay: RNA-seq
     *   ##biosample: brain
     *   ##biosample_type: tissue
     *   ##condition: tumor
     *   ##treatment: dexamethasone 100nM
     *   ##species: Homo sapiens
     *   ##replication_type: biological
     *   ##platform: PacBio Sequel II
     *   ##pipeline: StringTie
     *   ##pipeline_version: v2.2.1
     *   ##annotation_source: GENCODE
     *   ##annotation_version: v44
     *   ##source_url: https://...
     *   ##publication: DOI:...
     *
     * Unknown properties are stored in sample_info.attributes.
     * If no id is provided, defaults to filename stem.
     *
     * @param filepath Path to GFF/GTF file
     * @return sample_info populated with extracted metadata
     */
    static sample_info parse_header(const std::filesystem::path& filepath);

    /**
     * Process all transcripts from a single gene
     * Called by builder for parallel chromosome processing
     * @param grove Grove to add entries to
     * @param grove_mutex Mutex protecting grove operations (for thread-safety)
     * @param gene_entries All GFF entries for this gene (exons, CDS, UTR, codons)
     * @param exon_cache Chromosome-level exon cache (for deduplication)
     * @param segment_cache Chromosome-level segment cache (for deduplication)
     * @param sample_id Sample registry ID for pan-transcriptome tracking
     * @param segment_count Reference to segment counter (incremented for each segment created)
     */
    static void process_gene(
        grove_type& grove,
        std::mutex& grove_mutex,
        const std::vector<gio::gff_entry>& gene_entries,
        exon_cache_type& exon_cache,
        segment_cache_type& segment_cache,
        std::optional<uint32_t> sample_id,
        size_t& segment_count,
        const build_options& opts,
        build_counters& counters,
        quant_sidecar::SampleStreamWriter* sidecar_writer = nullptr
    );

private:
    static void process_transcript(
        grove_type& grove,
        std::mutex& grove_mutex,
        const std::string& transcript_id,
        const std::vector<gio::gff_entry>& all_entries,
        std::map<gdt::genomic_coordinate, key_ptr>& exon_keys,
        std::unordered_map<std::string, key_ptr>& segment_keys,
        std::optional<uint32_t> sample_id,
        size_t& segment_count,
        const build_options& opts,
        build_counters& counters,
        quant_sidecar::SampleStreamWriter* sidecar_writer = nullptr
    );

    /**
     * Annotate exons with overlapping features (CDS, UTR, codons)
     * @param grove Grove containing the exons
     * @param exon_keys Map of exon intervals to their keys
     * @param annotations Vector of annotation entries (CDS, UTR, etc.)
     */
    static void annotate_exons(
        grove_type& grove,
        const std::map<gdt::genomic_coordinate, key_ptr>& exon_keys,
        const std::vector<gio::gff_entry>& annotations
    );

    // ========== Transcript processing helpers ==========

    /**
     * Extract exon entries and sort in 5'→3' biological order
     * @return Sorted exon entries (empty if no exons found)
     */
    static std::vector<gio::gff_entry> extract_sorted_exons(
        const std::vector<gio::gff_entry>& transcript_entries
    );

    /**
     * Insert new exon or reuse existing one with same coordinates
     * @param gff_source GFF column 2 value (e.g., HAVANA, ENSEMBL, StringTie)
     * @return Key pointer to the exon (new or existing)
     */
    static key_ptr insert_exon(
        grove_type& grove,
        std::mutex& grove_mutex,
        const gio::gff_entry& exon_entry,
        const std::string& transcript_id,
        std::map<gdt::genomic_coordinate, key_ptr>& exon_cache,
        std::optional<uint32_t> sample_id,
        const std::string& gff_source = ""
    );

    // ========== Attribute extraction helpers ==========

    static std::optional<std::string> extract_gene_id(
        const std::map<std::string, std::string, std::less<>>& attributes);

    static std::optional<std::string> extract_transcript_id(
        const std::map<std::string, std::string, std::less<>>& attributes);

    static std::optional<std::string> extract_attribute(
        const std::map<std::string, std::string, std::less<>>& attributes,
        const std::string& key
    );
};

#endif //ATROPLEX_BUILD_GFF_HPP