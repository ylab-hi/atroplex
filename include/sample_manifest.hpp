/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#ifndef ATROPLEX_SAMPLE_MANIFEST_HPP
#define ATROPLEX_SAMPLE_MANIFEST_HPP

#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <optional>

#include "sample_info.hpp"

/**
 * Parser for sample manifest files (TSV format)
 *
 * Manifest format:
 * #sample_id    file    tissue    condition    [expression_file]    [key=value ...]
 *
 * Lines starting with # are treated as comments (first # line can define headers)
 * Columns after expression_file are parsed as key=value attributes
 *
 * Example:
 * #sample_id    file    tissue    condition    expression_file    species
 * brain_tumor   brain_tumor.gtf    brain    tumor    brain_tpm.tsv    human
 * liver_ctrl    liver_ctrl.gtf    liver    control    liver_tpm.tsv    human
 */
class sample_manifest {
public:
    /**
     * Parse manifest file and return vector of sample_info
     * @param filepath Path to manifest TSV file
     * @return Vector of parsed sample_info entries
     * @throws std::runtime_error if file cannot be read or is malformed
     */
    static std::vector<sample_info> parse(const std::filesystem::path& filepath) {
        std::vector<sample_info> samples;

        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open manifest file: " + filepath.string());
        }

        std::vector<std::string> headers;
        std::string line;
        size_t line_num = 0;

        while (std::getline(file, line)) {
            line_num++;

            // Skip empty lines
            if (line.empty() || std::all_of(line.begin(), line.end(), ::isspace)) {
                continue;
            }

            // Handle header line (first comment line starting with #)
            if (line[0] == '#') {
                if (headers.empty()) {
                    headers = parse_header(line.substr(1));
                }
                continue;
            }

            // Parse data line
            auto fields = split_tsv(line);
            if (fields.size() < 2) {
                throw std::runtime_error(
                    "Manifest line " + std::to_string(line_num) +
                    " must have at least sample_id and file columns"
                );
            }

            sample_info sample(fields[0]);
            sample.source_file = fields[1];

            // Parse remaining fields based on headers or position
            if (!headers.empty()) {
                parse_with_headers(sample, fields, headers);
            } else {
                parse_positional(sample, fields);
            }

            samples.push_back(std::move(sample));
        }

        return samples;
    }

    /**
     * Create a sample_info from CLI arguments
     * @param sample_id Sample identifier
     * @param file Source annotation file
     * @param tissue Optional tissue type
     * @param condition Optional condition
     * @return Constructed sample_info
     */
    static sample_info from_cli(
        const std::string& sample_id,
        const std::filesystem::path& file,
        const std::optional<std::string>& tissue = std::nullopt,
        const std::optional<std::string>& condition = std::nullopt
    ) {
        sample_info sample(sample_id, file);
        if (tissue) sample.tissue = *tissue;
        if (condition) sample.condition = *condition;
        return sample;
    }

private:
    static std::vector<std::string> parse_header(const std::string& header_line) {
        return split_tsv(header_line);
    }

    static std::vector<std::string> split_tsv(const std::string& line) {
        std::vector<std::string> fields;
        std::istringstream ss(line);
        std::string field;

        while (std::getline(ss, field, '\t')) {
            // Trim whitespace
            size_t start = field.find_first_not_of(" \t\r\n");
            size_t end = field.find_last_not_of(" \t\r\n");
            if (start != std::string::npos) {
                fields.push_back(field.substr(start, end - start + 1));
            } else {
                fields.push_back("");
            }
        }

        return fields;
    }

    static void parse_with_headers(
        sample_info& sample,
        const std::vector<std::string>& fields,
        const std::vector<std::string>& headers
    ) {
        for (size_t i = 0; i < fields.size() && i < headers.size(); i++) {
            const auto& header = headers[i];
            const auto& value = fields[i];

            if (value.empty()) continue;

            if (header == "sample_id") {
                // Already set
            } else if (header == "file") {
                sample.source_file = value;
            } else if (header == "tissue") {
                sample.tissue = value;
            } else if (header == "condition") {
                sample.condition = value;
            } else if (header == "species") {
                sample.species = value;
            } else if (header == "expression_file") {
                sample.expression_file = value;
            } else if (header == "annotation_source") {
                sample.annotation_source = value;
            } else if (header == "annotation_version") {
                sample.annotation_version = value;
            } else {
                // Store as custom attribute
                sample.attributes[header] = value;
            }
        }
    }

    static void parse_positional(
        sample_info& sample,
        const std::vector<std::string>& fields
    ) {
        // Positional: sample_id, file, tissue, condition, expression_file, ...
        // Index 0 (sample_id) and 1 (file) already set
        if (fields.size() > 2 && !fields[2].empty()) {
            sample.tissue = fields[2];
        }
        if (fields.size() > 3 && !fields[3].empty()) {
            sample.condition = fields[3];
        }
        if (fields.size() > 4 && !fields[4].empty()) {
            sample.expression_file = fields[4];
        }
        // Additional fields as key=value pairs
        for (size_t i = 5; i < fields.size(); i++) {
            auto eq_pos = fields[i].find('=');
            if (eq_pos != std::string::npos) {
                std::string key = fields[i].substr(0, eq_pos);
                std::string value = fields[i].substr(eq_pos + 1);
                sample.attributes[key] = value;
            }
        }
    }
};

#endif //ATROPLEX_SAMPLE_MANIFEST_HPP