/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SEGMENT_BUILDER_HPP
#define ATROPLEX_SEGMENT_BUILDER_HPP

#include <string>
#include <optional>
#include <vector>
#include <map>
#include <unordered_map>
#include <mutex>

#include "build_summary.hpp"
#include "genomic_feature.hpp"
#include "quant_sidecar.hpp"
#include "sample_info.hpp"

/// Classification of how a subsequence aligns within a parent exon chain.
/// See absorption_rules.txt for full documentation.
enum class subsequence_type {
    NONE,               ///< Not a contiguous subsequence
    FSM,                ///< Rule 0: full-length match (same exon count, all boundaries within fuzzy tolerance) → ABSORB
    ISM_5PRIME,         ///< Rule 1: 5' intact, 3' truncated → KEEP
    ISM_3PRIME,         ///< Rule 2: 5' truncated (1-2 exons), 3' intact → ABSORB
    DEGRADATION_3PRIME, ///< Rule 3: 5' truncated (3+ exons), 3' intact → DROP/KEEP
    INTERNAL_FRAGMENT   ///< Rule 4: both ends truncated → DROP/KEEP
};

/**
 * Format-agnostic segment creation and absorption logic.
 *
 * Shared by build_gff (GTF/GFF input) and future build_bam (BAM input).
 * All methods are static — no instance state.
 *
 * Absorption rules 0-8 are implemented here. See absorption_rules.txt
 * for full documentation of each rule.
 */
class segment_builder {
public:
    /**
     * Generate a structure key from ordered exon coordinates.
     * Format: seqid:strand:start-end,start-end,...
     */
    [[nodiscard]] static std::string make_exon_structure_key(
        const std::string& seqid,
        const std::vector<gdt::genomic_coordinate>& exon_coords
    );

    /**
     * Walk a segment's exon chain via SEGMENT_TO_EXON → EXON_TO_EXON
     * graph edges. Returns the ordered exon key_ptrs.
     */
    [[nodiscard]] static std::vector<key_ptr> walk_exon_chain(
        const grove_type& grove,
        key_ptr segment_key,
        size_t segment_index
    );

    /**
     * Create new segment or match against existing ones via absorption rules.
     *
     * Absorption candidates are found via grove.intersect() — spatially,
     * not by gene_id. This enables cross-annotation absorption and
     * eliminates gene_id as a structural dependency (issue #40).
     */
    static void create_segment(
        grove_type& grove,
        std::mutex& grove_mutex,
        const std::string& transcript_id,
        const std::string& seqid,
        char strand,
        size_t span_start,
        size_t span_end,
        int exon_count,
        const std::vector<gdt::genomic_coordinate>& exon_coords,
        const std::vector<key_ptr>& exon_chain,
        segment_cache_type& segment_cache,
        uint32_t gene_idx,
        std::optional<uint32_t> sample_id,
        const std::string& gff_source,
        size_t& segment_count,
        float expression_value,
        const std::string& transcript_biotype,
        bool absorb,
        size_t fuzzy_tolerance,
        build_counters& counters,
        quant_sidecar::SampleStreamWriter* sidecar_writer = nullptr,
        bool annotated_loci_only = false
    );

    /**
     * Classify how sub aligns within parent's exon chain (pointer comparison).
     * Distinguishes ISM_3PRIME (1-2 missing) from DEGRADATION_3PRIME (3+ missing).
     */
    [[nodiscard]] static subsequence_type classify_subsequence(
        const std::vector<key_ptr>& sub,
        const std::vector<key_ptr>& parent
    );

    /**
     * Fuzzy subsequence matching (fallback when pointer identity fails).
     * Slides a window over the parent chain, checks all exon boundaries within tolerance.
     */
    [[nodiscard]] static subsequence_type fuzzy_classify_subsequence(
        const std::vector<key_ptr>& sub,
        const std::vector<key_ptr>& parent,
        size_t tolerance_bp = 5
    );

    /**
     * Merge transcript metadata into an existing segment (shared by dedup + absorption).
     *
     * When `sidecar_writer` is non-null and `expression_value >= 0`, appends
     * a record for the target segment's index to the sample's stream — so
     * expression survives the common path where a sample's transcript matches
     * a previously created segment (e.g. annotations processed first, then
     * samples merging in).
     */
    static void merge_into_segment(
        key_ptr target_seg,
        const std::string& transcript_id,
        std::optional<uint32_t> sample_id,
        const std::string& gff_source,
        float expression_value,
        const std::string& transcript_biotype,
        quant_sidecar::SampleStreamWriter* sidecar_writer = nullptr
    );

    // ── Absorption rule helpers ─────────────────────────────────────

    /** Rule 5: Check if two exon chains share the same intron chain */
    [[nodiscard]] static bool has_same_intron_chain(
        const std::vector<key_ptr>& chain_a,
        const std::vector<key_ptr>& chain_b
    );

    /** Rule 5: Check if terminal exon boundaries are within tolerance */
    [[nodiscard]] static bool terminal_boundaries_within_tolerance(
        const std::vector<key_ptr>& chain_a,
        const std::vector<key_ptr>& chain_b,
        size_t tolerance_bp = 50
    );

    /** Mono-exon classification labels (Rules 6/7/8) */
    enum class mono_exon_class { GENE_OVERLAP, INTRON_RETENTION, INTERGENIC };

    /** Check if the given sample_id corresponds to a reference annotation */
    [[nodiscard]] static bool is_annotation_sample(std::optional<uint32_t> sample_id);

    /** Check if the parent segment belongs to a reference annotation */
    [[nodiscard]] static bool is_parent_annotation(key_ptr parent_seg);

private:
    struct spatial_candidate {
        key_ptr segment;
        std::vector<key_ptr> exon_chain;
    };

    /**
     * Reverse absorption: apply Rules 0/2/3/4/5 to existing segments
     * when a new (potentially longer) segment is created. Iterates the
     * `candidates` vector already gathered by `create_segment` rather
     * than re-issuing a `grove.intersect()` post-insertion — the new
     * segment is identified by `new_seg` and skipped during the sweep.
     */
    static void try_reverse_absorption(
        grove_type& grove,
        key_ptr new_seg,
        const std::vector<key_ptr>& new_exon_chain,
        const std::vector<spatial_candidate>& candidates,
        segment_cache_type& segment_cache,
        const std::string& seqid,
        size_t fuzzy_tolerance = 5
    );

    /**
     * Rules 6/7/8: classify a mono-exon transcript using already-gathered
     * spatial candidates. For mono-exon transcripts the create_segment
     * gather coord equals the mono coord, so the candidate set is
     * exactly the overlap set this classifier needs — no separate
     * `grove.intersect()` is required.
     */
    [[nodiscard]] static mono_exon_class classify_mono_exon(
        const gdt::genomic_coordinate& mono_coord,
        const std::vector<spatial_candidate>& candidates
    );

    /** Rule 5: try merging as terminal variant (same intron chain, TSS/TES <50bp) */
    [[nodiscard]] static bool try_terminal_variant_merge(
        const std::vector<key_ptr>& exon_chain,
        const std::vector<spatial_candidate>& candidates,
        const std::string& transcript_id, std::optional<uint32_t> sample_id,
        const std::string& gff_source, float expression_value,
        const std::string& transcript_biotype,
        quant_sidecar::SampleStreamWriter* sidecar_writer,
        build_counters& counters);

    /** Rules 6/7/8: handle mono-exon sample transcripts. Returns true if handled (dropped or kept). */
    [[nodiscard]] static bool handle_mono_exon(
        const std::vector<key_ptr>& exon_chain,
        std::optional<uint32_t> sample_id, bool absorb,
        const std::vector<spatial_candidate>& candidates,
        build_counters& counters);

    /** Rules 0(fuzzy)/1/2/3/4: try subsequence absorption. Returns true if absorbed. */
    [[nodiscard]] static bool try_subsequence_absorption(
        const std::vector<key_ptr>& exon_chain,
        const std::vector<spatial_candidate>& candidates,
        size_t span_start, size_t span_end, size_t fuzzy_tolerance,
        const std::string& transcript_id, std::optional<uint32_t> sample_id,
        const std::string& gff_source, float expression_value,
        const std::string& transcript_biotype,
        quant_sidecar::SampleStreamWriter* sidecar_writer,
        build_counters& counters);

    /** Gene_idx inheritance + annotated-loci-only filter. Returns true if discarded. */
    [[nodiscard]] static bool apply_gene_idx_inheritance(
        const std::vector<spatial_candidate>& candidates,
        std::optional<uint32_t> sample_id, bool annotated_loci_only,
        uint32_t& gene_idx, build_counters& counters);
};

#endif //ATROPLEX_SEGMENT_BUILDER_HPP