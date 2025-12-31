/*
 * SPDX-License-Identifier: MIT
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the MIT
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "genomic_feature.hpp"

#include <sstream>

/**
 * Parse GFF/GTF attributes string to extract key-value pairs
 */
static std::string extract_attribute(const std::string& attributes, const std::string& key) {
    // Try GTF format: key "value"
    size_t pos = attributes.find(key);
    if (pos != std::string::npos) {
        size_t start = attributes.find('"', pos);
        size_t end = attributes.find('"', start + 1);
        if (start != std::string::npos && end != std::string::npos) {
            return attributes.substr(start + 1, end - start - 1);
        }
    }

    // Try GFF3 format: key=value
    std::string key_equals = key + "=";
    pos = attributes.find(key_equals);
    if (pos != std::string::npos) {
        size_t start = pos + key_equals.length();
        size_t end = attributes.find(';', start);
        if (end == std::string::npos) end = attributes.length();
        return attributes.substr(start, end - start);
    }

    return "";
}

genomic_feature genomic_feature::from_gff_entry(
    const std::string& attributes,
    size_t node_id,
    feature_type type
) {
    genomic_feature feature(node_id, type);

    // Extract common attributes
    feature.id = extract_attribute(attributes, "ID");
    if (feature.id.empty()) {
        feature.id = extract_attribute(attributes, "exon_id");
    }

    feature.gene_id = extract_attribute(attributes, "gene_id");
    feature.gene_name = extract_attribute(attributes, "gene_name");
    feature.gene_type = extract_attribute(attributes, "gene_type");

    if (feature.gene_type.empty()) {
        feature.gene_type = extract_attribute(attributes, "gene_biotype");
    }

    // Transcript information
    std::string transcript_id = extract_attribute(attributes, "transcript_id");
    if (transcript_id.empty()) {
        transcript_id = extract_attribute(attributes, "Parent");
    }
    if (!transcript_id.empty()) {
        feature.transcript_ids.insert(transcript_id);
    }

    feature.transcript_type = extract_attribute(attributes, "transcript_type");
    if (feature.transcript_type.empty()) {
        feature.transcript_type = extract_attribute(attributes, "transcript_biotype");
    }

    // Source
    feature.source = extract_attribute(attributes, "source");

    // Exon number
    std::string exon_num_str = extract_attribute(attributes, "exon_number");
    if (!exon_num_str.empty()) {
        try {
            feature.exon_number = std::stoi(exon_num_str);
        } catch (...) {
            feature.exon_number = -1;
        }
    }

    return feature;
}