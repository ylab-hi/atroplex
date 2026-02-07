/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "sample_manifest.hpp"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

sample_manifest::sample_manifest(const std::filesystem::path& manifest_path)
    : manifest_path_(manifest_path) {

    std::ifstream file(manifest_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open manifest file: " + manifest_path.string());
    }

    std::string line;

    // Read header
    if (!std::getline(file, line)) {
        throw std::runtime_error("Manifest file is empty: " + manifest_path.string());
    }

    auto col_indices = parse_header(line);

    // Check required column
    if (col_indices.find("file") == col_indices.end()) {
        throw std::runtime_error("Manifest missing required 'file' column: " + manifest_path.string());
    }

    // Get manifest directory for resolving relative paths
    std::filesystem::path manifest_dir = manifest_path.parent_path();
    if (manifest_dir.empty()) {
        manifest_dir = ".";
    }

    // Read data rows
    size_t line_num = 1;
    while (std::getline(file, line)) {
        line_num++;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        try {
            samples_.push_back(parse_row(line, col_indices, manifest_dir));
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing manifest line " +
                std::to_string(line_num) + ": " + e.what());
        }
    }
}

std::unordered_map<std::string, size_t> sample_manifest::parse_header(
    const std::string& header_line) {

    std::unordered_map<std::string, size_t> indices;
    auto columns = split(header_line, '\t');

    for (size_t i = 0; i < columns.size(); ++i) {
        std::string col_name = to_lower(columns[i]);
        // Trim whitespace
        col_name.erase(0, col_name.find_first_not_of(" \t\r\n"));
        col_name.erase(col_name.find_last_not_of(" \t\r\n") + 1);

        if (!col_name.empty()) {
            indices[col_name] = i;
        }
    }

    return indices;
}

sample_info sample_manifest::parse_row(
    const std::string& line,
    const std::unordered_map<std::string, size_t>& col_indices,
    const std::filesystem::path& manifest_dir) {

    auto fields = split(line, '\t');

    // Helper to get field value (empty string or "." means not set)
    auto get_field = [&](const std::string& col_name) -> std::string {
        auto it = col_indices.find(col_name);
        if (it == col_indices.end() || it->second >= fields.size()) {
            return "";
        }
        std::string value = fields[it->second];
        // Trim whitespace
        value.erase(0, value.find_first_not_of(" \t\r\n"));
        value.erase(value.find_last_not_of(" \t\r\n") + 1);
        // "." means empty/missing
        if (value == ".") {
            return "";
        }
        return value;
    };

    sample_info info;

    // Required: file path
    std::string file_str = get_field("file");
    if (file_str.empty()) {
        throw std::runtime_error("Missing required 'file' field");
    }

    // Resolve path relative to manifest directory
    std::filesystem::path file_path(file_str);
    if (file_path.is_relative()) {
        file_path = manifest_dir / file_path;
    }
    info.source_file = std::filesystem::canonical(file_path);

    // ID (auto-generate from filename if not provided)
    info.id = get_field("id");
    if (info.id.empty()) {
        info.id = info.source_file.stem().string();
    }

    // Type (default: "sample")
    std::string type_val = get_field("type");
    if (!type_val.empty()) {
        info.type = to_lower(type_val);
    }

    // Description
    info.description = get_field("description");

    // Biological metadata
    info.assay = get_field("assay");
    info.biosample_type = get_field("biosample_type");
    info.biosample = get_field("biosample");
    info.condition = get_field("condition");
    info.treatment = get_field("treatment");
    info.species = get_field("species");
    info.replication_type = get_field("replication_type");

    // Technical metadata
    info.platform = get_field("platform");
    info.pipeline = get_field("pipeline");
    info.pipeline_version = get_field("pipeline_version");

    // Annotation metadata
    info.annotation_source = get_field("annotation_source");
    info.annotation_version = get_field("annotation_version");

    // Links
    info.source_url = get_field("source_url");
    info.publication = get_field("publication");

    return info;
}

std::vector<std::string> sample_manifest::split(const std::string& str, char delim) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delim)) {
        tokens.push_back(token);
    }

    return tokens;
}

std::string sample_manifest::to_lower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
        [](unsigned char c) { return std::tolower(c); });
    return result;
}

void sample_manifest::write_template(const std::filesystem::path& output_path,
                                     bool include_examples) {
    std::ofstream file(output_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create template file: " + output_path.string());
    }

    // Write header with all supported columns
    file << "file\tid\ttype\tassay\tbiosample_type\tbiosample\tcondition\ttreatment\t"
         << "species\treplication_type\tplatform\tpipeline\tpipeline_version\t"
         << "annotation_source\tannotation_version\tsource_url\tpublication\tdescription\n";

    if (include_examples) {
        // Reference annotation example
        file << "gencode.v44.gtf\tGENCODE_v44\tannotation\t.\t.\t.\t.\t.\t"
             << "Homo sapiens\t.\t.\t.\t.\t"
             << "GENCODE\tv44\thttps://www.gencodegenes.org\t.\t"
             << "Evidence-based annotation of the human genome (GRCh38)\n";

        // Sample assembly examples
        file << "brain_tumor_rep1.gtf\tBT_rep1\tsample\tRNA-seq\ttissue\tbrain\ttumor\t.\t"
             << "Homo sapiens\tbiological\tPacBio Sequel II\tStringTie\tv2.2.1\t"
             << ".\t.\t.\t.\t.\n";

        file << "brain_normal_rep1.gtf\tBN_rep1\tsample\tRNA-seq\ttissue\tbrain\thealthy\t.\t"
             << "Homo sapiens\tbiological\tPacBio Sequel II\tStringTie\tv2.2.1\t"
             << ".\t.\t.\t.\t.\n";

        file << "hela_treated.gtf\tHeLa_dex\tsample\tRNA-seq\tcell line\tHeLa\ttreated\tdexamethasone 100nM 1h\t"
             << "Homo sapiens\t.\tONT PromethION\tFLAIR\tv2.0\t"
             << ".\t.\t.\t.\t.\n";
    }
}