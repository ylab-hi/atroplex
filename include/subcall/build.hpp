/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_BUILD_HPP
#define ATROPLEX_SUBCALL_BUILD_HPP

#include "subcall/subcall.hpp"

namespace subcall {

/**
 * Build subcommand: create genogrove index from annotation files.
 *
 * Creates a .gg index file that can be reused across multiple
 * analysis runs, avoiding repeated parsing of annotation files.
 */
class build : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "build"; }
    std::string description() const override {
        return "Build genogrove index from annotation files";
    }
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_BUILD_HPP