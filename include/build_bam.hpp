/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_BUILD_BAM_HPP
#define ATROPLEX_BUILD_BAM_HPP

#include <string>
#include <optional>
#include <vector>
#include <map>
#include <unordered_map>
#include <filesystem>
#include <mutex>

#include "genomic_feature.hpp"
#include "sample_info.hpp"
#include "read_cluster.hpp"

/**
 * BAM-specific grove builder.
 *
 * Clusters aligned reads by splice junction pattern, derives exon coordinates
 * from consensus junctions, and indexes them via segment_builder.
 *
 * Uses the same deduplication and absorption rules as build_gff — the
 * segment_builder is format-agnostic.
 */
class build_bam {
public:
    /**
     * Build grove from a BAM file.
     * @param grove Grove to populate
     * @param filepath Path to BAM file
     * @param sample_id Sample registry ID for pan-transcriptome tracking
     * @param exon_caches Chromosome-level exon caches (for cross-file deduplication)
     * @param segment_caches Chromosome-level segment caches
     * @param gene_indices Per-gene segment index (for absorption)
     * @param segment_count Running segment counter
     * @param min_expression Minimum read count to include a cluster (default: -1 = no filter)
     * @param absorb Enable absorption rules
     * @param fuzzy_tolerance Fuzzy matching tolerance in bp
     */
    static void build(grove_type& grove,
                      const std::filesystem::path& filepath,
                      std::optional<uint32_t> sample_id,
                      chromosome_exon_caches& exon_caches,
                      chromosome_segment_caches& segment_caches,
                      chromosome_gene_segment_indices& gene_indices,
                      size_t& segment_count,
                      float min_expression = -1.0f,
                      bool absorb = true,
                      size_t fuzzy_tolerance = 5);

    /**
     * Parse BAM header to extract sample metadata.
     * Reads @RG and @PG header lines for sample ID, platform, pipeline info.
     * Falls back to filename stem if no metadata found.
     */
    static sample_info parse_header(const std::filesystem::path& filepath);

private:
    /**
     * Insert an exon by coordinate (no GFF entry needed).
     * Reuses existing exon if one with identical coordinates exists in the cache.
     */
    static key_ptr insert_exon(
        grove_type& grove,
        std::mutex& grove_mutex,
        const std::string& seqid,
        char strand,
        size_t start,
        size_t end,
        const std::string& gene_id,
        const std::string& gene_name,
        const std::string& gene_biotype,
        const std::string& transcript_id,
        exon_cache_type& exon_cache,
        std::optional<uint32_t> sample_id,
        const std::string& source = "BAM"
    );

    /**
     * Convert a read cluster's consensus junctions into exon intervals.
     * N junctions → N+1 exon intervals.
     */
    static std::vector<std::pair<size_t, size_t>> cluster_to_exon_coords(
        const read_cluster& cluster
    );

    /**
     * Resolve gene assignment for a cluster by spatial query.
     * Returns (gene_id, gene_name, gene_biotype) of the best overlapping gene,
     * or a synthetic gene ID if no overlap found.
     */
    static std::tuple<std::string, std::string, std::string> resolve_gene(
        grove_type& grove,
        const read_cluster& cluster
    );
};

#endif //ATROPLEX_BUILD_BAM_HPP