/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_GENOGROVE_BUILDER_HPP
#define ATROPLEX_GENOGROVE_BUILDER_HPP

// standard
#include <string>
#include <vector>

// genogrove
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/gff_reader.hpp>
#include <genogrove/io/filetype_detector.hpp>

// class
#include "genomic_feature.hpp"

namespace gdt = genogrove::data_type;
namespace gst = genogrove::structure;
namespace gio = genogrove::io;

// Type aliases for grove structure
using grove_type = gst::grove<gdt::interval, genomic_feature, edge_metadata>;
using key_ptr = gdt::key<gdt::interval, genomic_feature>*;

/**
 * GenogroveBuilder handles the creation of genogrove structures from various file types.
 * Processes gene-by-gene to minimize memory usage.
 * Workflow per gene:
 * 1. Insert exons into grove (spatial intervals)
 * 2. For each transcript: create segments and link exons
 * 3. Annotate exons with overlapping features (CDS, UTR, codons)
 *
 * Supports:
 * - GFF/GTF annotation files
 * - Future: BED, BAM, or other interval-based formats
 */
class genogrove_builder {
public:
    /**
     * Build genogrove from multiple input files
     * @param files Vector of file paths to process
     * @param order Order of the genogrove structure
     * @return Pointer to grove (caller owns the memory)
     */
    static grove_type* build_from_files(
        const std::vector<std::string>& files,
        int order
    );

    /**
     * Process a single GFF/GTF file and add to grove
     * Processes gene-by-gene for memory efficiency
     * @param grove Grove to add entries to
     * @param filepath Path to GFF/GTF file
     */
    static void build_from_gff(
        grove_type& grove,
        const std::filesystem::path& filepath
    );

private:
    /**
     * Process all transcripts from a single gene
     * @param grove Grove to add entries to
     * @param gene_entries All GFF entries for this gene (exons, CDS, UTR, codons)
     */
    static void process_gene(
        grove_type& grove,
        const std::vector<gio::gff_entry>& gene_entries
    );

    /**
     * Process a single transcript: create exon chain and segment
     * @param grove Grove to add entries to
     * @param transcript_id Transcript ID
     * @param all_entries All entries for this transcript (exons + annotations)
     * @param exon_keys Map of exon intervals to their keys (for reuse across transcripts)
     */
    static void process_transcript(
        grove_type& grove,
        const std::string& transcript_id,
        const std::vector<gio::gff_entry>& all_entries,
        std::unordered_map<gdt::interval, key_ptr>& exon_keys
    );

    /**
     * Annotate exons with overlapping features (CDS, UTR, codons)
     * @param grove Grove containing the exons
     * @param exon_keys Map of exon intervals to their keys
     * @param annotations Vector of annotation entries (CDS, UTR, etc.)
     */
    static void annotate_exons(
        grove_type& grove,
        const std::unordered_map<gdt::interval, key_ptr>& exon_keys,
        const std::vector<gio::gff_entry>& annotations
    );

    /**
     * Extract gene_id from GFF attributes
     */
    static std::string extract_gene_id(const std::map<std::string, std::string>& attributes);

    /**
     * Extract transcript_id from GFF attributes
     */
    static std::string extract_transcript_id(const std::string& attributes);
};

#endif //ATROPLEX_GENOGROVE_BUILDER_HPP