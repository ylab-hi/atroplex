/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_COMPACT_HPP
#define ATROPLEX_SUBCALL_COMPACT_HPP

#include "subcall/subcall.hpp"

namespace subcall {

/**
 * Compact subcommand: physically remove absorbed (tombstoned) segments
 * from an existing genogrove index, producing a smaller .ggx without
 * changing query semantics.
 *
 * Loads a pre-built .ggx (and its companion .qtx, when present),
 * unlinks every segment with `absorbed == true` from the B+ tree via
 * grove.remove_key(), drops orphan EXON_TO_EXON chain edges, and
 * re-serializes the grove to the output directory. The .qtx written
 * during build is already remapped against live segments
 * (`merge_to_qtx` consumes the tombstone remap at build time), so it
 * is copied through unchanged alongside the new .ggx and .ggx.summary.
 *
 * By default `compact` fails if no .qtx exists alongside the input
 * .ggx — pass `--no-qtx` to opt out when operating on a structure-only
 * index.
 */
class compact : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "compact"; }
    std::string description() const override {
        return "Compact a genogrove index by physically removing absorbed segments";
    }
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_COMPACT_HPP