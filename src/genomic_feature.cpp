/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "genomic_feature.hpp"

#include <sstream>

/**
 * Helper to extract attribute from map, return empty string if not found
 */
static std::string get_attribute(
    const std::map<std::string, std::string>& attributes,
    const std::string& key
) {
    auto it = attributes.find(key);
    return (it != attributes.end()) ? it->second : "";
}

exon_feature exon_feature::from_gff_entry(
    const std::map<std::string, std::string>& attributes,
    const std::string& seqid,
    const gdt::interval& interval,
    char strand
) {
    exon_feature feature;

    // Extract common attributes
    feature.id = get_attribute(attributes, "ID");
    if (feature.id.empty()) {
        feature.id = get_attribute(attributes, "exon_id");
    }

    std::string gid = get_attribute(attributes, "gene_id");
    std::string gname = get_attribute(attributes, "gene_name");
    std::string gbiotype = get_attribute(attributes, "gene_type");
    if (gbiotype.empty()) {
        gbiotype = get_attribute(attributes, "gene_biotype");
    }
    feature.gene_idx = gene_registry::instance().intern(gid, gname, gbiotype);

    // Transcript information
    std::string transcript_id = get_attribute(attributes, "transcript_id");
    if (transcript_id.empty()) {
        transcript_id = get_attribute(attributes, "Parent");
    }
    if (!transcript_id.empty()) {
        feature.transcript_ids.insert(transcript_registry::instance().intern(transcript_id));
    }

    // Note: GFF source (column 2) is added separately via add_source() in build_gff

    // Exon number
    std::string exon_num_str = get_attribute(attributes, "exon_number");
    if (!exon_num_str.empty()) {
        try {
            feature.exon_number = std::stoi(exon_num_str);
        } catch (...) {
            feature.exon_number = -1;
        }
    }

    return feature;
}