/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

// standard
#include <filesystem>
#include <algorithm>

// class
#include "utility.hpp"
#include "builder.hpp"
#include "build_gff.hpp"

// Natural chromosome sort comparator
static bool chromosome_compare(const std::string& a, const std::string& b) {
    auto get_chr_value = [](const std::string& s) -> std::pair<int, std::string> {
        if (s.size() > 3 && s.substr(0, 3) == "chr") {
            std::string suffix = s.substr(3);
            try {
                size_t pos;
                int num = std::stoi(suffix, &pos);
                if (pos == suffix.size()) {
                    return {num, ""};  // Pure numeric (chr1, chr22)
                }
            } catch (...) {}
            return {1000, suffix};  // Non-numeric (chrX, chrY, chrM)
        }
        return {2000, s};  // Non-chr prefixed
    };

    auto [num_a, str_a] = get_chr_value(a);
    auto [num_b, str_b] = get_chr_value(b);

    if (num_a != num_b) return num_a < num_b;
    return str_a < str_b;
}

void builder::build_from_samples(grove_type& grove,
                                  const std::vector<sample_info>& samples,
                                  uint32_t threads) {
    if (samples.empty()) {
        logging::warning("No samples provided to build genogrove");
        return;
    }

    if (threads > 1) {
        logging::info("Note: Multi-threading for build not yet optimized, using single thread");
    }

    logging::info("Populating grove from " + std::to_string(samples.size()) + " sample(s)");

    // Chromosome-level caches for deduplication across files
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    size_t segment_count = 0;

    // Process each sample sequentially
    for (const auto& info : samples) {
        const auto& filepath = info.source_file;

        if (!std::filesystem::exists(filepath)) {
            logging::error("File not found: " + filepath.string());
            continue;
        }

        // Detect file type
        gio::filetype_detector detector;
        auto [ftype, is_gzipped] = detector.detect_filetype(filepath);

        // Dispatch to appropriate builder based on file type
        if (ftype == gio::filetype::GFF || ftype == gio::filetype::GTF) {
            logging::info("Processing GFF/GTF file: " + filepath.string() +
                         (info.id.empty() ? "" : " (id: " + info.id + ")"));

            // Register sample_info in the registry to get uint32_t ID
            uint32_t registry_id = sample_registry::instance().register_data(info);

            // Build with persistent caches for cross-file deduplication
            build_gff::build(grove, filepath, registry_id, exon_caches, segment_caches, segment_count, threads);
        } else {
            logging::warning("Unsupported file type for: " + filepath.string());
        }
    }

    logging::info("Grove construction complete: " + std::to_string(segment_count) + " segments");

    // Collect and sort chromosome names naturally
    std::vector<std::string> seqids;
    for (const auto& [seqid, _] : segment_caches) {
        seqids.push_back(seqid);
    }
    std::sort(seqids.begin(), seqids.end(), chromosome_compare);

    // Output per-chromosome summary (exons and segments from caches)
    logging::info("Per-chromosome summary:");
    for (const auto& seqid : seqids) {
        size_t exon_count = exon_caches.count(seqid) ? exon_caches[seqid].size() : 0;
        size_t seg_count = segment_caches[seqid].size();
        logging::info("  " + seqid + ": " +
            std::to_string(exon_count) + " exons, " +
            std::to_string(seg_count) + " segments");
    }
}

void builder::build_from_files(grove_type& grove,
                                const std::vector<std::string>& files,
                                uint32_t threads) {
    // Parse headers to extract metadata from each file
    std::vector<sample_info> samples;
    samples.reserve(files.size());

    for (const auto& filepath : files) {
        // Parse metadata from GFF/GTF header (#property: value format)
        sample_info info = build_gff::parse_header(filepath);
        samples.push_back(std::move(info));
    }

    build_from_samples(grove, samples, threads);
}
