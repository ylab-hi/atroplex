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
#include <string>
#include <vector>
#include <cstdint>

// genogrove
#include <genogrove/io/filetype_detector.hpp>

// class
#include "genomic_feature.hpp"
#include "index_stats.hpp"
#include "sample_info.hpp"

namespace gio = genogrove::io;

/**
 * GenogroveBuilder - Toolkit for grove construction and modification
 *
 * Stateless builder that operates on grove references.
 * Supports:
 * - Initial construction from GFF/GTF annotation files
 * - Pan-transcriptome builds with sample tracking via manifest
 * - Extension with discovered transcripts (Phase 3)
 * - Addition of fusion segments (Phase 3)
 * - Future: BED, BAM, or other interval-based formats
 */
class builder {
public:
    /**
     * Build genogrove from sample_info entries (primary entry point)
     * Use this with sample_manifest for full metadata support.
     *
     * @param grove Reference to grove to populate
     * @param samples Vector of sample_info with metadata (from manifest or programmatic)
     * @param threads Number of threads (reserved for future use, currently ignored)
     */
    static index_stats build_from_samples(
        grove_type& grove,
        const std::vector<sample_info>& samples,
        uint32_t threads = 1,
        float min_expression = -1.0f
    );

    /**
     * Build genogrove from file paths (convenience wrapper)
     * Creates minimal sample_info with just source_file and auto-generated ID.
     * For full metadata support, use build_from_samples() with a manifest.
     *
     * @param grove Reference to grove to populate
     * @param files Vector of file paths to process
     * @param threads Number of threads (reserved for future use, currently ignored)
     */
    static index_stats build_from_files(
        grove_type& grove,
        const std::vector<std::string>& files,
        uint32_t threads = 1
    );

    // Future extension methods for Phase 2/3:
    // static void add_discovered_transcripts(grove_type& grove, const std::vector<alignment_entry>& reads);
    // static void add_fusion_segments(grove_type& grove, const fusion_data& fusions);
};

#endif //ATROPLEX_GENOGROVE_BUILDER_HPP
