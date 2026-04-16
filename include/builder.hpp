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
#include "build_summary.hpp"
#include "genomic_feature.hpp"
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
    static build_summary build_from_samples(
        grove_type& grove,
        const std::vector<sample_info>& samples,
        uint32_t threads,
        const expression_filters& filters,
        bool absorb,
        int min_replicates,
        size_t fuzzy_tolerance,
        bool prune_tombstones,
        bool include_scaffolds = false,
        /// Final .qtx output path. Empty means no quantification sidecar
        /// is written. A temp directory is derived from this path
        /// (`<qtx_path>.tmp/`) for per-sample streams that are K-way
        /// merged into the final single file at end of build.
        const std::string& qtx_path = "",
        chromosome_exon_caches* out_exon_caches = nullptr,
        chromosome_gene_segment_indices* out_gene_indices = nullptr
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
    static build_summary build_from_files(
        grove_type& grove,
        const std::vector<std::string>& files,
        uint32_t threads = 1,
        chromosome_exon_caches* out_exon_caches = nullptr
    );

    // Future extension methods for Phase 2/3:
    // static void add_discovered_transcripts(grove_type& grove, const std::vector<alignment_entry>& reads);
    // static void add_fusion_segments(grove_type& grove, const fusion_data& fusions);

private:
    /// Post-build merge of biological replicates within groups.
    /// Iterates exon and segment caches (no grove traversal).
    /// Returns the number of replicate entries collapsed into merged groups.
    static size_t merge_replicates(
        chromosome_exon_caches& exon_caches,
        chromosome_segment_caches& segment_caches,
        int min_replicates
    );

    /// Handle absorbed (tombstoned) segments after the build.
    /// Always counts tombstones, prunes them from segment_caches and
    /// gene_indices, and returns the count.
    ///
    /// When `physical == true`, additionally unlinks the tombstoned keys
    /// from the B+ tree via grove.remove_key() and drops orphan
    /// EXON_TO_EXON chain edges. This produces a smaller .ggx at the
    /// cost of an O(N × E) sweep (genogrove's graph_overlay::remove_edges_to
    /// is O(E) per call), which can block for minutes to hours on
    /// realistic indices. Default `false` — tombstones stay in the tree
    /// and every consumer already filters them defensively.
    static size_t remove_tombstones(
        grove_type& grove,
        chromosome_segment_caches& segment_caches,
        chromosome_gene_segment_indices& gene_indices,
        bool physical
    );
};

#endif //ATROPLEX_GENOGROVE_BUILDER_HPP
