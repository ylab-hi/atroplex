/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the MIT
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
#include <filesystem>

// genogrove
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/gff_reader.hpp>
#include <genogrove/data_type/genomic_coordinate.hpp>

// class
#include "genomic_feature.hpp"

namespace gdt = genogrove::data_type;
namespace gst = genogrove::structure;
namespace gio = genogrove::io;

// Type aliases
using grove_type = gst::grove<gdt::interval, genomic_feature, edge_metadata>;
using key_ptr = gdt::key<gdt::interval, genomic_feature>*;

/**
 * GFF/GTF-specific grove builder
 * Handles construction of grove structures from GFF/GTF annotation files
 *
 * Processing strategy:
 * - Gene-by-gene processing for memory efficiency
 * - Creates exon features with CDS/UTR overlaps
 * - Creates segment features representing transcript paths
 * - Links exons within segments via graph edges
 */
class build_gff {
public:
    /**
     * Build grove from single GFF/GTF file
     * @param grove Grove to populate
     * @param filepath Path to GFF/GTF file
     */
    static void build(grove_type& grove, const std::filesystem::path& filepath);

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
        std::map<gdt::interval, key_ptr>& exon_keys
    );

    /**
     * Annotate exons with overlapping features (CDS, UTR, codons)
     * @param grove Grove containing the exons
     * @param exon_keys Map of exon intervals to their keys
     * @param annotations Vector of annotation entries (CDS, UTR, etc.)
     */
    static void annotate_exons(
        grove_type& grove,
        const std::map<gdt::interval, key_ptr>& exon_keys,
        const std::vector<gio::gff_entry>& annotations
    );

    /**
     * Extract gene_id from GFF attributes
     */
    static std::optional<std::string> extract_gene_id(
        const std::map<std::string, std::string>& attributes);

    /**
     * Extract transcript_id from GFF attributes
     */
    static std::optional<std::string> extract_transcript_id(
        const std::map<std::string, std::string>& attributes);

    /**
     * Extract generic attribute from GFF attributes map
     */
    static std::optional<std::string> extract_attribute(
        const std::map<std::string, std::string>& attributes,
        const std::string& key
    );
};

#endif //ATROPLEX_BUILD_GFF_HPP