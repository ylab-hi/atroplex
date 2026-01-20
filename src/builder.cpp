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

            // Parse header to extract metadata and create sample_info
            sample_info info = build_gff::parse_header(filepath);

            // Register in the registry to get uint32_t ID
            uint32_t registry_id = sample_registry::instance().register_data(std::move(info));

            // Build with the registry ID for provenance tracking
            build_gff::build(grove, filepath, registry_id);
        } else {
            logging::warning("Unsupported file type for: " + filepath);
        }
    }

    logging::info("Grove construction complete");
}