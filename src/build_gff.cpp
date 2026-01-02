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
#include <map>
#include <unordered_map>
#include <algorithm>

// class
#include "utility.hpp"
#include "build_gff.hpp"

void build_gff::build(
    grove_type& grove,
    const std::filesystem::path& filepath
) {
    gio::gff_reader reader(filepath.string());

    // Buffer for gene entries - process one gene at a time
    std::unordered_map<std::string, std::vector<gio::gff_entry>> genes;
    std::string current_gene_id;
    size_t line_count = 0;

    for(const auto& entry : reader) {
        line_count++;

        // Process relevant features: exon, CDS, UTRs, codons
        // We need all of these to properly annotate exons with overlapping features
        if (entry.type != "exon" &&
            entry.type != "CDS" &&
            entry.type != "five_prime_UTR" &&
            entry.type != "three_prime_UTR" &&
            entry.type != "start_codon" &&
            entry.type != "stop_codon") {
            continue;
        }

        std::string gene_id = extract_gene_id(entry.attributes);
        if (gene_id.empty()) {
            continue; // Skip entries without gene_id
        }

        // If we've moved to a new gene, process the previous one
        if (!current_gene_id.empty() && gene_id != current_gene_id) {
            process_gene(grove, genes[current_gene_id]);
            genes.erase(current_gene_id); // Free memory
        }

        current_gene_id = gene_id;
        genes[gene_id].push_back(entry);
    }

    // Process last gene
    if (!current_gene_id.empty()) {
        process_gene(grove, genes[current_gene_id]);
    }

    logging::info("Processed " + std::to_string(line_count) + " lines from " + filepath.filename().string());
}

void build_gff::process_gene(
    grove_type& grove,
    const std::vector<gio::gff_entry>& gene_entries
) {
    // Separate exons from annotation features
    std::vector<gio::gff_entry> annotations;

    // Group entries by transcript
    std::unordered_map<std::string, std::vector<gio::gff_entry>> transcripts;

    for (const auto& entry : gene_entries) {
        std::string transcript_id = extract_transcript_id(entry.attributes);
        if (!transcript_id.empty()) {
            transcripts[transcript_id].push_back(entry);
        }

        // Collect annotations for later
        if (entry.type != "exon") {
            annotations.push_back(entry);
        }
    }

    // Map to track exon keys (exons can be shared across transcripts)
    // Using std::map because gdt::interval doesn't have a hash function
    std::map<gdt::interval, key_ptr> exon_keys;

    // Process each transcript: create exon chains and segments
    for (auto& [transcript_id, entries] : transcripts) {
        process_transcript(grove, transcript_id, entries, exon_keys);
    }

    // Annotate all exons with overlapping features
    annotate_exons(grove, exon_keys, annotations);
}

void build_gff::process_transcript(
    grove_type& grove,
    const std::string& transcript_id,
    const std::vector<gio::gff_entry>& all_entries,
    std::map<gdt::interval, key_ptr>& exon_keys
) {
    // Extract only exons
    std::vector<gio::gff_entry> exons;
    for (const auto& entry : all_entries) {
        if (entry.type == "exon") {
            exons.push_back(entry);
        }
    }

    if (exons.empty()) {
        return;
    }

    // Sort exons by genomic position
    std::sort(exons.begin(), exons.end(),
             [](const gio::gff_entry& a, const gio::gff_entry& b) {
                 if (a.seqid != b.seqid) return a.seqid < b.seqid;
                 return a.interval.get_start() < b.interval.get_start();
             });

    // Collect exon keys for this transcript
    std::vector<key_ptr> transcript_exon_keys;
    std::string seqid = exons[0].seqid;
    char strand = exons[0].strand.value_or('+');  // Extract from optional

    // Step 1: Insert or reuse exons
    for (const auto& exon_entry : exons) {
        key_ptr exon_key;

        // Check if this exon already exists
        auto it = exon_keys.find(exon_entry.interval);
        if (it != exon_keys.end()) {
            // Reuse existing exon, add this transcript to it
            exon_key = it->second;
            // Update the feature to include this transcript using public getter
            auto& feature_data = exon_key->get_data();
            if (is_exon(feature_data)) {
                get_exon(feature_data).transcript_ids.insert(transcript_id);
            }
        } else {
            // Create new exon feature
            exon_feature exon_feat = exon_feature::from_gff_entry(
                extract_attribute(exon_entry.attributes, ""),
                exon_entry.seqid,
                exon_entry.interval,
                exon_entry.strand.value_or('+')  // Extract from optional
            );
            exon_feat.transcript_ids.insert(transcript_id);

            // Wrap in variant and insert
            genomic_feature feature = exon_feat;
            exon_key = grove.insert_data(exon_entry.seqid, exon_entry.interval, feature);
            exon_keys[exon_entry.interval] = exon_key;
        }

        transcript_exon_keys.push_back(exon_key);
    }

    // Step 2: Link exons with edges (exon chain for this transcript)
    for (size_t i = 0; i < transcript_exon_keys.size() - 1; ++i) {
        edge_metadata edge_meta(transcript_id, edge_metadata::edge_type::EXON_TO_EXON);
        grove.add_edge(transcript_exon_keys[i], transcript_exon_keys[i + 1], edge_meta);
    }

    // Step 3: Create segment spanning all exons in this transcript
    if (!transcript_exon_keys.empty()) {
        // Segment interval: first exon start → last exon end
        gdt::interval segment_interval;
        segment_interval.set_start(exons.front().interval.get_start());
        segment_interval.set_end(exons.back().interval.get_end());

        // Create segment feature
        segment_feature seg_feat(segment_feature::source_type::REFERENCE);
        seg_feat.id = transcript_id + "_segment";
        seg_feat.transcript_ids.insert(transcript_id);
        seg_feat.exon_count = static_cast<int>(transcript_exon_keys.size());

        // Copy gene metadata from first exon using public getter
        auto& first_feature = transcript_exon_keys[0]->get_data();
        if (is_exon(first_feature)) {
            const auto& first_exon = get_exon(first_feature);
            seg_feat.gene_id = first_exon.gene_id;
            seg_feat.gene_name = first_exon.gene_name;
            seg_feat.gene_biotype = first_exon.gene_biotype;
        }

        // Build coordinate string: chr:strand:start-end
        seg_feat.coordinate = seqid + ":" + strand + ":" +
                             std::to_string(segment_interval.get_start()) + "-" +
                             std::to_string(segment_interval.get_end());

        // Wrap in variant and insert
        genomic_feature feature = seg_feat;
        key_ptr segment_key = grove.insert_data(seqid, segment_interval, feature);

        // Link segment to its first exon
        edge_metadata segment_edge(transcript_id, edge_metadata::edge_type::SEGMENT_TO_EXON);
        grove.add_edge(segment_key, transcript_exon_keys[0], segment_edge);
    }
}

void build_gff::annotate_exons(
    grove_type& grove,
    const std::map<gdt::interval, key_ptr>& exon_keys,
    const std::vector<gio::gff_entry>& annotations
) {
    // For each annotation feature (CDS, UTR, codon)
    for (const auto& annot : annotations) {
        // Find overlapping exons
        for (const auto& [exon_interval, exon_key] : exon_keys) {
            // Check for overlap
            if (exon_interval.get_end() < annot.interval.get_start() ||
                exon_interval.get_start() > annot.interval.get_end()) {
                continue; // No overlap
            }

            // Compute overlapping interval
            gdt::interval overlap;
            overlap.set_start(std::max(exon_interval.get_start(), annot.interval.get_start()));
            overlap.set_end(std::min(exon_interval.get_end(), annot.interval.get_end()));

            // Add to exon's overlapping features using public getter
            auto& feature_data = exon_key->get_data();
            if (is_exon(feature_data)) {
                get_exon(feature_data).overlapping_features[annot.type].push_back(overlap);
            }
        }
    }
}

std::string build_gff::extract_gene_id(const std::map<std::string, std::string>& attributes) {
    // Try common attribute names
    auto it = attributes.find("gene_id");
    if (it != attributes.end()) {
        return it->second;
    }

    it = attributes.find("GeneID");
    if (it != attributes.end()) {
        return it->second;
    }

    return "";
}

std::string build_gff::extract_transcript_id(const std::map<std::string, std::string>& attributes) {
    // Try GTF format: transcript_id
    auto it = attributes.find("transcript_id");
    if (it != attributes.end()) {
        return it->second;
    }

    // Try GFF3 format: Parent
    it = attributes.find("Parent");
    if (it != attributes.end()) {
        return it->second;
    }

    return "";
}

std::string build_gff::extract_attribute(
    const std::map<std::string, std::string>& attributes,
    const std::string& key
) {
    auto it = attributes.find(key);
    if (it != attributes.end()) {
        return it->second;
    }
    return "";
}