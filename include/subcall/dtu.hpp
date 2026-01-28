/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_DTU_HPP
#define ATROPLEX_SUBCALL_DTU_HPP

#include "subcall/subcall.hpp"

#include <map>
#include <set>
#include <vector>

namespace subcall {

/**
 * Sample metadata for DTU analysis
 */
struct dtu_sample {
    std::string sample_id;
    std::string file_path;
    std::string condition;
    std::string tissue;
    std::map<std::string, std::string> attributes;
};

/**
 * Per-sample transcript quantification
 */
struct transcript_count {
    std::string sample_id;
    std::string gene_id;
    std::string transcript_id;
    uint32_t read_count;
    double proportion;
};

/**
 * DTU result for a single transcript
 */
struct dtu_result {
    std::string gene_id;
    std::string gene_name;
    std::string transcript_id;
    double prop_condition1;
    double prop_condition2;
    double delta_proportion;
    double p_value;
    double fdr;
    bool significant;
};

/**
 * Contrast specification for DTU comparison
 */
struct dtu_contrast {
    std::string condition1;
    std::string condition2;

    static dtu_contrast parse(const std::string& spec);
};

/**
 * DTU subcommand: differential transcript usage analysis.
 *
 * Compares transcript proportions between conditions to identify
 * genes with differential transcript usage.
 */
class dtu : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "dtu"; }
    std::string description() const override {
        return "Differential transcript usage analysis";
    }

private:
    std::vector<dtu_sample> samples;
    std::map<std::string, std::vector<transcript_count>> counts_by_sample;
    std::vector<dtu_contrast> contrasts;

    std::vector<dtu_sample> load_manifest(const std::string& path);
    void quantify_samples(const cxxopts::ParseResult& args);
    void load_quantifications(const std::string& quant_dir);
    std::vector<dtu_result> run_dtu_analysis(const dtu_contrast& contrast,
                                              const cxxopts::ParseResult& args);
    void apply_fdr_correction(std::vector<dtu_result>& results, double fdr_threshold);

    void write_results(const std::string& path,
                       const dtu_contrast& contrast,
                       const std::vector<dtu_result>& results);
    void write_proportions(const std::string& path);
    void write_summary(const std::string& path,
                       const std::vector<std::pair<dtu_contrast, std::vector<dtu_result>>>& all_results,
                       const cxxopts::ParseResult& args);
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_DTU_HPP