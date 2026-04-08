/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SUBCALL_EXPORT_GTF_HPP
#define ATROPLEX_SUBCALL_EXPORT_GTF_HPP

#include <optional>
#include <unordered_set>

#include "subcall/subcall.hpp"

namespace subcall {

/**
 * Export subcommand: reconstruct per-sample GTF files from a .ggx index.
 *
 * Walks the grove, traverses exon chains via graph edges, resolves metadata
 * from registries, and writes one GTF per sample. Optional filters restrict
 * output to specific samples, genes, regions, biotypes, or sources.
 */
class export_gtf : public subcall {
public:
    cxxopts::Options parse_args(int argc, char** argv) override;
    void validate(const cxxopts::ParseResult& args) override;
    void execute(const cxxopts::ParseResult& args) override;

    std::string name() const override { return "export"; }
    std::string description() const override {
        return "Export per-sample GTF files from a genogrove index";
    }

private:
    struct region_spec {
        std::string seqid;
        size_t start;
        size_t end;
    };

    struct export_filters {
        std::unordered_set<std::string> sample_ids;
        std::unordered_set<std::string> gene_ids;
        std::optional<region_spec> region;
        size_t min_samples = 0;
        bool conserved_only = false;
        std::unordered_set<std::string> biotypes;
        std::unordered_set<std::string> sources;
    };

    export_filters filters_;

    void parse_filters(const cxxopts::ParseResult& args);

    [[nodiscard]] static region_spec parse_region(const std::string& region_str);
};

} // namespace subcall

#endif // ATROPLEX_SUBCALL_EXPORT_GTF_HPP