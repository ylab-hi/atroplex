/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_INSPECT_HPP
#define ATROPLEX_SUBCALL_INSPECT_HPP

#include "subcall/subcall.hpp"

namespace subcall {

/**
 * Inspect subcommand: streaming pan-transcriptome analysis — overview,
 * per-sample / per-source / biotype breakdowns, exon + segment sharing,
 * conserved-exon details, splicing hubs, splicing events. Single-pass
 * grove traversal, memory scales with feature count (not sample count).
 *
 * Outputs land directly under the common `-o/--output-dir`, grouped by
 * category subdirectories (`overview/`, `sharing/`, `splicing_hubs/`,
 * `splicing_events/`) — no extra wrapper folder.
 */
class inspect : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "inspect"; }
    std::string description() const override {
        return "Pan-transcriptome inspection: overview, sharing, splicing hubs and events";
    }
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_INSPECT_HPP