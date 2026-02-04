/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SAMPLE_MANIFEST_HPP
#define ATROPLEX_SAMPLE_MANIFEST_HPP

#include <filesystem>
#include <string>
#include <vector>
#include <optional>
#include <unordered_map>

#include "sample_info.hpp"

/**
 * Parser for sample manifest files (TSV format)
 *
 * Manifest format:
 *   - Tab-separated values
 *   - First row is header with column names
 *   - Required column: "file" (path to annotation/sample file)
 *   - All other columns are optional
 *   - Use "." for empty/missing values (VCF convention)
 *
 * Supported columns (case-insensitive, ENCODE-aligned):
 *
 *   REQUIRED:
 *     file              Path to GFF/GTF file (relative to manifest or absolute)
 *
 *   IDENTIFIERS:
 *     id                Sample identifier (auto-generated from filename if empty)
 *
 *   BIOLOGICAL (ENCODE-aligned):
 *     assay             Experimental assay (RNA-seq, CAGE, polyA RNA-seq, etc.)
 *     biosample_type    Sample type (tissue, cell line, primary cell, organoid)
 *     biosample         Specific biosample (brain, liver, HeLa, GM12878, etc.)
 *     condition         Experimental condition (tumor, healthy, treated, control)
 *     treatment         Treatment details (dexamethasone 100nM 1h, etc.)
 *     species           Species (Homo sapiens, Mus musculus, etc.)
 *     replication_type  Replicate type (biological, technical, isogenic)
 *
 *   TECHNICAL:
 *     platform          Sequencing platform (PacBio Sequel II, ONT PromethION, etc.)
 *     pipeline          Analysis pipeline (StringTie, FLAIR, IsoQuant, etc.)
 *     pipeline_version  Pipeline version (v2.2.1, etc.)
 *
 *   ANNOTATION (for reference files):
 *     annotation_source Source database (GENCODE, Ensembl, RefSeq)
 *     annotation_version Version string (v44, release 110, etc.)
 *
 *   LINKS:
 *     source_url        URL to original data (ENCODE portal, GEO, etc.)
 *     publication       Publication DOI or reference
 *
 *
 * Example manifest (samples.tsv) - empty values are allowed:
 *
 *   file	id	assay	biosample_type	biosample	condition	species	platform	pipeline
 *   gencode.v44.gtf	GENCODE_v44					Homo sapiens
 *   brain_tumor.gtf	BT001	RNA-seq	tissue	brain	tumor	Homo sapiens	PacBio Sequel II	StringTie
 *   brain_normal.gtf	BN001	RNA-seq	tissue	brain	healthy	Homo sapiens	PacBio Sequel II	StringTie
 *   hela_treated.gtf	HT001	RNA-seq	cell line	HeLa	treated	Homo sapiens	ONT PromethION	FLAIR
 *
 * Minimal manifest (just file paths):
 *
 *   file
 *   sample1.gtf
 *   sample2.gtf
 *
 * Usage:
 *   sample_manifest manifest("samples.tsv");
 *   for (const auto& info : manifest.get_samples()) {
 *       // info.source_file contains resolved path
 *       // info has all metadata from manifest
 *   }
 */
class sample_manifest {
public:
    /**
     * Parse manifest from file
     * @param manifest_path Path to TSV manifest file
     * @throws std::runtime_error if file cannot be read or has invalid format
     */
    explicit sample_manifest(const std::filesystem::path& manifest_path);

    /**
     * Get all parsed sample entries
     * @return Vector of sample_info with populated metadata
     */
    const std::vector<sample_info>& get_samples() const { return samples_; }

    /**
     * Get number of samples in manifest
     */
    size_t size() const { return samples_.size(); }

    /**
     * Check if manifest is empty
     */
    bool empty() const { return samples_.empty(); }

    /**
     * Get manifest file path
     */
    const std::filesystem::path& path() const { return manifest_path_; }

    /**
     * Iterator support
     */
    auto begin() const { return samples_.begin(); }
    auto end() const { return samples_.end(); }

    /**
     * Create a template manifest with all supported columns
     * @param output_path Path to write template file
     * @param include_examples If true, include example data rows
     */
    static void write_template(const std::filesystem::path& output_path,
                               bool include_examples = true);

private:
    std::filesystem::path manifest_path_;
    std::vector<sample_info> samples_;

    /**
     * Parse header row to get column indices
     * @return Map of column name (lowercase) to column index
     */
    static std::unordered_map<std::string, size_t> parse_header(
        const std::string& header_line);

    /**
     * Parse a data row into sample_info
     * @param line Tab-separated data line
     * @param col_indices Column name to index mapping
     * @param manifest_dir Directory of manifest file (for resolving relative paths)
     * @return Populated sample_info
     */
    static sample_info parse_row(
        const std::string& line,
        const std::unordered_map<std::string, size_t>& col_indices,
        const std::filesystem::path& manifest_dir);

    /**
     * Split string by delimiter
     */
    static std::vector<std::string> split(const std::string& str, char delim);

    /**
     * Convert string to lowercase
     */
    static std::string to_lower(const std::string& str);
};

#endif //ATROPLEX_SAMPLE_MANIFEST_HPP
