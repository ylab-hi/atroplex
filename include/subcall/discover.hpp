/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_DISCOVER_HPP
#define ATROPLEX_SUBCALL_DISCOVER_HPP

#include "subcall/subcall.hpp"

#include <genogrove/io/bam_reader.hpp>
#include <genogrove/io/filetype_detector.hpp>

#include "read_cluster.hpp"
#include "transcript_matcher.hpp"

namespace gio = genogrove::io;

namespace subcall {

/**
 * Discover subcommand: detect transcripts from long-read sequencing data.
 *
 * Pipeline:
 * 1. Loads/builds grove from annotations
 * 2. Clusters reads by splice junction signature
 * 3. Matches clusters to reference transcripts
 * 4. Discovers novel transcripts
 */
class discover : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "discover"; }
    std::string description() const override {
        return "Detect transcripts from long-read sequencing data";
    }

private:
    void process_reads(const cxxopts::ParseResult& args);
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_DISCOVER_HPP