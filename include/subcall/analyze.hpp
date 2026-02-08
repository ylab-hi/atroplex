/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_ANALYZE_HPP
#define ATROPLEX_SUBCALL_ANALYZE_HPP

#include "subcall/subcall.hpp"

namespace subcall {

/**
 * Analyze subcommand: full pan-transcriptome analysis with per-sample
 * statistics, isoform diversity (Jaccard), and CSV output.
 */
class analyze : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "analyze"; }
    std::string description() const override {
        return "Full pan-transcriptome analysis with per-sample statistics";
    }
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_ANALYZE_HPP