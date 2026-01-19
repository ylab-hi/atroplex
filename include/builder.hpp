/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_GENOGROVE_BUILDER_HPP
#define ATROPLEX_GENOGROVE_BUILDER_HPP

// standard
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// genogrove
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/filetype_detector.hpp>

// class
#include "genomic_feature.hpp"
#include "sample_info.hpp"

namespace gdt = genogrove::data_type;
namespace gst = genogrove::structure;
namespace gio = genogrove::io;

// Type aliases for grove structure
using grove_type = gst::grove<gdt::interval, genomic_feature, edge_metadata>;
using key_ptr = gdt::key<gdt::interval, genomic_feature>*;

/**
 * GenogroveBuilder - Toolkit for grove construction and modification
 *
 * Stateless builder that operates on grove references.
 * Supports:
 * - Initial construction from GFF/GTF annotation files
 * - Automatic metadata extraction from file headers (GENCODE, Ensembl)
 * - Pan-transcriptome builds with sample tracking
 * - Extension with discovered transcripts (Phase 3)
 * - Addition of fusion segments (Phase 3)
 * - Future: BED, BAM, or other interval-based formats
 */
class builder {
public:
    /**
     * Build genogrove from multiple input files
     * Automatically:
     * - Detects file types and dispatches to appropriate builders
     * - Parses file headers to extract metadata (GENCODE version, genome build, etc.)
     * - Registers sample_info in the registry for provenance tracking
     *
     * @param grove Reference to grove to populate
     * @param files Vector of file paths to process
     */
    static void build_from_files(
        grove_type& grove,
        const std::vector<std::string>& files
    );

    // Future extension methods for Phase 2/3:
    // static void add_discovered_transcripts(grove_type& grove, const std::vector<alignment_entry>& reads);
    // static void add_fusion_segments(grove_type& grove, const fusion_data& fusions);
};

#endif //ATROPLEX_GENOGROVE_BUILDER_HPP