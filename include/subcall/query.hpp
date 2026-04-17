/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_QUERY_HPP
#define ATROPLEX_SUBCALL_QUERY_HPP

#include "subcall/subcall.hpp"
#include "transcript_matcher.hpp"

#include <map>
#include <set>
#include <vector>

namespace subcall {

/**
 * Contrast specification for differential transcript usage
 */
struct query_contrast {
    std::string group_a;
    std::string group_b;

    static query_contrast parse(const std::string& spec);
};

/**
 * Per-transcript classification result with per-sample context
 */
struct query_result {
    std::string transcript_id;
    std::string gene_id;
    std::string gene_name;

    // Structural classification (SQANTI-like)
    structural_category category = structural_category::INTERGENIC;
    subcategory subcat;

    // Per-sample presence (sample_id -> present)
    std::vector<uint32_t> sample_ids;

    // Per-sample expression (from .qtx sidecar, populated at query time
    // when a sidecar is loaded). Empty when no sidecar is available or
    // the matched segment has no records in the sidecar.
    std::map<uint32_t, float> sample_expression;

    // segment_index of the best-matching catalog segment (only set when
    // the query has a match). Used by the DTU path to look up per-sample
    // values from the .qtx sidecar without repeating the grove query.
    std::optional<uint64_t> matched_segment_index;

    // Match quality
    double junction_match_score = 0.0;
    int matching_junctions = 0;
    int total_query_junctions = 0;
    int total_ref_junctions = 0;

    // Novel splice sites
    int known_donors = 0;
    int known_acceptors = 0;
    int novel_donors = 0;
    int novel_acceptors = 0;
};

/**
 * DTU result for a single transcript (when --contrast is used)
 */
struct dtu_result {
    std::string gene_id;
    std::string gene_name;
    std::string transcript_id;
    double prop_group_a;
    double prop_group_b;
    double delta_proportion;
    double p_value;
    double fdr;
    bool significant;
};

/**
 * Query subcommand: classify transcripts against the pan-transcriptome index.
 *
 * Pipeline:
 * 1. Build/load grove from annotations (via common options)
 * 2. Read input transcripts (GTF/GFF)
 * 3. Classify each against the grove (structural category + per-sample context)
 * 4. Optionally: run differential transcript usage between groups (--contrast)
 *
 * Output:
 * - Per-transcript classification TSV (structural category, per-sample presence/expression)
 * - Summary statistics
 * - DTU results TSV (when --contrast is provided)
 */
class query : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "query"; }
    std::string description() const override {
        return "Classify transcripts against the pan-transcriptome index";
    }

private:
    std::vector<query_contrast> contrasts_;

    /**
     * Classify input transcripts from GTF/GFF against the grove
     */
    std::vector<query_result> classify_transcripts(const std::string& input_path);

    /**
     * Run differential transcript usage for a single contrast.
     *
     * Groups samples by `sample_info::group`, computes per-transcript
     * proportions within each group, and applies a chi-squared test on
     * the 2×k contingency table (k transcripts per gene). Requires a
     * `.qtx` sidecar — per-sample expression comes from `qr.sample_expression`
     * populated in `classify_transcripts` via `qtx_reader_ptr()->lookup()`.
     *
     * Samples whose `group` field is empty or `.` are skipped.
     */
    std::vector<dtu_result> run_dtu(
        const query_contrast& contrast,
        const std::vector<query_result>& results,
        const cxxopts::ParseResult& args);

    /**
     * Benjamini–Hochberg FDR correction. Mutates `results.fdr` and
     * `results.significant` in place using `fdr_threshold`.
     */
    void apply_fdr_correction(std::vector<dtu_result>& results,
                              double fdr_threshold);

    // Output
    void write_classification(const std::string& path,
                              const std::vector<query_result>& results);
    void write_dtu_results(const std::string& path,
                           const query_contrast& contrast,
                           const std::vector<dtu_result>& results);
    void write_summary(const std::string& path,
                       const std::vector<query_result>& classification,
                       const std::vector<std::pair<query_contrast,
                           std::vector<dtu_result>>>& dtu_results);
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_QUERY_HPP