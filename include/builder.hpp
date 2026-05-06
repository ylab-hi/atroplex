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
#include <unordered_map>
#include <unordered_set>
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
        const build_options& opts = {},
        chromosome_exon_caches* out_exon_caches = nullptr
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
    ///
    /// `segment_caches` is the build-time cache used by build_from_samples
    /// to dedupe segments by ordered-exon-coordinate hash. The `compact`
    /// subcommand operates on a loaded grove without that cache; pass an
    /// empty map and the cache-cleanup phase becomes a no-op.
    ///
    /// When `out_tombstoned_segment_indices` is non-null, the set of
    /// segment_index values for tombstoned segments is written there
    /// (used by the sidecar merge to exclude their records from the
    /// final .qtx).
    ///
    /// When `out_tombstone_remap` is non-null, tombstones that carry an
    /// `absorbed_into_idx` (set by `absorb_into_parent`) produce an entry
    /// `{candidate_idx -> parent_idx}` — with transitive chains resolved
    /// via path compression. `merge_to_qtx` consumes this to rewrite each
    /// record's segment_index at emission time instead of dropping it, so
    /// expression from samples whose transcripts were processed before
    /// the parent existed is preserved on the live parent's `.qtx` block.
    /// Tombstones without `absorbed_into_idx` (Rule 3/4 `tombstone_candidate`
    /// drops) are **not** added to the remap — they remain in the drop set.
    static size_t remove_tombstones(
        grove_type& grove,
        chromosome_segment_caches& segment_caches,
        bool physical,
        std::unordered_set<uint64_t>* out_tombstoned_segment_indices = nullptr,
        std::unordered_map<uint64_t, uint64_t>* out_tombstone_remap = nullptr
    );
};

#endif //ATROPLEX_GENOGROVE_BUILDER_HPP
