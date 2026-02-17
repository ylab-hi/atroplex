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
#include <filesystem>
#include <mutex>

// genogrove
#include <genogrove/io/gff_reader.hpp>

// class
#include "genomic_feature.hpp"
#include "sample_info.hpp"

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
                      chromosome_gene_segment_indices& gene_indices,
                      size_t& segment_count,
                      uint32_t num_threads = 0,
                      float min_expression = -1.0f);

    /**
     * Generate a structure key from ordered exon coordinates
     * Format: seqid:strand:start-end,start-end,...
     * Example: chr1:+:100-200,300-400,500-600
     * @param seqid Chromosome/sequence ID
     * @param exon_coords Ordered vector of exon coordinates (must be in 5'→3' order)
     * @return Structure key string
     */
    static std::string make_exon_structure_key(
        const std::string& seqid,
        const std::vector<gdt::genomic_coordinate>& exon_coords
    );

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
        gene_segment_index_type& gene_index,
        std::optional<uint32_t> sample_id,
        size_t& segment_count,
        float min_expression = -1.0f
    );

private:
    /**
     * Process a single transcript: create exon chain and segment
     * @param grove Grove to add entries to
     * @param grove_mutex Mutex protecting grove operations
     * @param transcript_id Transcript ID
     * @param all_entries All entries for this transcript (exons + annotations)
     * @param exon_keys Map of exon coordinates to their keys (for reuse across transcripts)
     * @param segment_keys Map of exon structure keys to segment keys (for deduplication)
     * @param sample_id Sample registry ID for pan-transcriptome tracking
     * @param segment_count Reference to segment counter (incremented when segment is created)
     */
    static void process_transcript(
        grove_type& grove,
        std::mutex& grove_mutex,
        const std::string& transcript_id,
        const std::vector<gio::gff_entry>& all_entries,
        std::map<gdt::genomic_coordinate, key_ptr>& exon_keys,
        std::unordered_map<std::string, key_ptr>& segment_keys,
        gene_segment_index_type& gene_index,
        std::optional<uint32_t> sample_id,
        size_t& segment_count,
        float min_expression = -1.0f
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

    /**
     * Create new segment or reuse existing one with same exon structure
     * @param gff_source GFF column 2 value (e.g., HAVANA, ENSEMBL, StringTie)
     * @param expression_value Expression value from transcript entry (-1 if not available)
     * @param transcript_biotype Transcript biotype (e.g., protein_coding, retained_intron)
     */
    static void create_segment(
        grove_type& grove,
        std::mutex& grove_mutex,
        const std::string& transcript_id,
        const std::vector<gio::gff_entry>& sorted_exons,
        const std::vector<gdt::genomic_coordinate>& exon_coords,
        const std::vector<key_ptr>& exon_chain,
        std::unordered_map<std::string, key_ptr>& segment_cache,
        gene_segment_index_type& gene_index,
        std::optional<uint32_t> sample_id,
        const std::string& gff_source,
        size_t& segment_count,
        float expression_value = -1.0f,
        const std::string& transcript_biotype = ""
    );

    // ========== ISM absorption helpers ==========

    /**
     * Check if sub is a contiguous subsequence of parent (pointer comparison)
     * Used to detect ISM transcripts that are truncated versions of full-length segments
     */
    static bool is_contiguous_subsequence(
        const std::vector<key_ptr>& sub,
        const std::vector<key_ptr>& parent
    );

    /**
     * Merge transcript metadata into an existing segment (shared by dedup + absorption)
     */
    static void merge_into_segment(
        key_ptr target_seg,
        const std::string& transcript_id,
        std::optional<uint32_t> sample_id,
        const std::string& gff_source,
        float expression_value,
        const std::string& transcript_biotype
    );

    /**
     * Reverse absorption: absorb existing shorter segments into a new longer one
     */
    static void try_reverse_absorption(
        gene_segment_index_type& gene_index,
        const std::string& gene_id,
        key_ptr new_seg,
        const std::vector<key_ptr>& new_exon_chain,
        segment_cache_type& segment_cache
    );

    // ========== Attribute extraction helpers ==========

    static std::optional<std::string> extract_gene_id(
        const std::map<std::string, std::string>& attributes);

    static std::optional<std::string> extract_transcript_id(
        const std::map<std::string, std::string>& attributes);

    static std::optional<std::string> extract_attribute(
        const std::map<std::string, std::string>& attributes,
        const std::string& key
    );
};

#endif //ATROPLEX_BUILD_GFF_HPP