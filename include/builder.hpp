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
 * GenogroveBuilder - Main interface for grove construction
 *
 * Dispatches to format-specific builders based on file type detection.
 * Supports:
 * - GFF/GTF annotation files (via gff_grove_builder)
 * - Future: BED, BAM, or other interval-based formats
 */
class builder {
public:
    /**
     * Build genogrove from multiple input files
     * Automatically detects file types and dispatches to appropriate builders
     *
     * @param files Vector of file paths to process
     * @param order Order of the genogrove structure
     * @return Pointer to grove (caller owns the memory)
     */
    static grove_type* build_from_files(
        const std::vector<std::string>& files,
        int order
    );
};

#endif //ATROPLEX_GENOGROVE_BUILDER_HPP