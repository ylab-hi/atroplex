/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

// standard
#include <filesystem>

// class
#include "utility.hpp"
#include "builder.hpp"
#include "build_gff.hpp"

void builder::build_from_files(grove_type& grove, const std::vector<std::string>& files) {
    if (files.empty()) {
        logging::warning("No files provided to build genogrove");
        return;
    }

    logging::info("Populating grove from " + std::to_string(files.size()) + " file(s)");
    if (!sample_id.empty()) {
        logging::info("Tagging features with sample: " + sample_id);
    }

    // Process each file
    for (const auto& filepath : files) {
        if (!std::filesystem::exists(filepath)) {
            logging::error("File not found: " + filepath);
            continue;
        }

        // Detect file type
        gio::filetype_detector detector;
        auto [ftype, is_gzipped] = detector.detect_filetype(filepath);

        // Dispatch to appropriate builder based on file type
        if (ftype == gio::filetype::GFF || ftype == gio::filetype::GTF) {
            logging::info("Processing GFF/GTF file: " + filepath);
            build_gff::build(grove, filepath, sample_id);
        } else {
            logging::warning("Unsupported file type for: " + filepath);
        }

        // Future format support:
        // else if (ftype == gio::filetype::BED) {
        //     build_bed::build(grove, filepath, sample_id);
        // }
        // else if (ftype == gio::filetype::BAM) {
        //     build_bam::build(grove, filepath, sample_id);
        // }
    }

    logging::info("Grove construction complete");
}

void builder::build_from_samples(grove_type& grove, const std::vector<sample_info>& samples) {
    if (samples.empty()) {
        logging::warning("No samples provided to build genogrove");
        return;
    }

    logging::info("Building pan-transcriptome from " + std::to_string(samples.size()) + " sample(s)");

    for (const auto& sample : samples) {
        if (sample.source_file.empty()) {
            logging::warning("Sample " + sample.id + " has no source file, skipping");
            continue;
        }

        if (!std::filesystem::exists(sample.source_file)) {
            logging::error("File not found for sample " + sample.id + ": " + sample.source_file.string());
            continue;
        }

        logging::info("Processing sample: " + sample.id + " (" + sample.source_file.string() + ")");

        // Detect file type
        gio::filetype_detector detector;
        auto [ftype, is_gzipped] = detector.detect_filetype(sample.source_file.string());

        // Dispatch to appropriate builder based on file type
        if (ftype == gio::filetype::GFF || ftype == gio::filetype::GTF) {
            build_gff::build(grove, sample.source_file, sample.id);
        } else {
            logging::warning("Unsupported file type for sample " + sample.id);
        }
    }

    logging::info("Pan-transcriptome construction complete with " + std::to_string(samples.size()) + " samples");
}