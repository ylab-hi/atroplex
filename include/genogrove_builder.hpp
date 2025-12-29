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
#include "transcript_graph.hpp"
#include "genomic_feature.hpp"

namespace gdt = genogrove::data_type;
namespace gst = genogrove::structure;
namespace gio = genogrove::io;

/**
 * Container for genomic data structures
 * - grove: Spatial index mapping intervals to genomic features
 * - graph: Unified graph with exon and segment nodes
 */
struct genomic_structures {
    ggs::grove<gdt::interval, genomic_feature>* grove;
    transcript_graph* graph;

    genomic_structures() : grove(nullptr), graph(nullptr) {}
    ~genomic_structures() {
        delete grove;
        delete graph;
    }
};

/**
 * GenogroveBuilder handles the creation of genogrove structures from various file types.
 * Builds both interval tree (grove) and two-layer transcript graph.
 * Supports:
 * - GFF/GTF annotation files
 * - Future: BED, BAM, or other interval-based formats
 */
class genogrove_builder {
    public:
        /**
         * Build genogrove and transcript graph from multiple input files
         * @param files Vector of file paths to process
         * @param order Order of the genogrove structure
         * @return Pointer to genomic_structures (caller owns the memory)
         */
        static gst::grove<gdt::interval, genomic_feature>* build_from_files(
            const std::vector<std::string>& files,
            int order
        );
        build_from_gff()

};

#endif //ATROPLEX_GENOGROVE_BUILDER_HPP