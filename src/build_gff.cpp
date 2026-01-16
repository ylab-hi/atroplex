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
#include <map>
#include <unordered_map>
#include <algorithm>

// class
#include "utility.hpp"
#include "build_gff.hpp"

// genogrove

void build_gff::build(grove_type& grove,
    const std::filesystem::path& filepath,
    const std::string& sample_id) {
    gio::gff_reader reader(filepath.string());

    // Buffer for gene entries - process one gene at a time
    std::unordered_map<std::string, std::vector<gio::gff_entry>> genes;
    std::string current_gene_id;
    size_t line_count = 0;

    for(const auto& entry : reader) {
        line_count++;

        // Extract gene_id to group entries by gene
        std::optional<std::string> gene_id = entry.get_gene_id();
        if (!gene_id.has_value()) {
            continue; // Skip entries without gene_id
        }

        // If we've moved to a new gene, process the previous one
        if (!current_gene_id.empty() && gene_id.value() != current_gene_id) {
            process_gene(grove, genes[current_gene_id], sample_id);
            genes.erase(current_gene_id); // Free memory
        }

        current_gene_id = gene_id.value();
        genes[gene_id.value()].push_back(entry);
    }

    // Process last gene
    if (!current_gene_id.empty()) {
        process_gene(grove, genes[current_gene_id], sample_id);
    }

    if(line_count % 10000 == 0) {
        logging::info("Processed " + std::to_string(line_count) + " lines from " + filepath.filename().string());
    }
}

void build_gff::process_gene(
    grove_type& grove,
    const std::vector<gio::gff_entry>& gene_entries,
    const std::string& sample_id
) {
    // Group entries by transcript
    std::unordered_map<std::string, std::vector<gio::gff_entry>> transcripts;

    for (const auto& entry : gene_entries) {
        // Extract transcript_id to group by transcript
        std::optional<std::string> transcript_id = entry.get_transcript_id();
        if (transcript_id.has_value()) {
            transcripts[transcript_id.value()].push_back(entry);
        }
    }

    // Map to track exon keys (exons can be shared across transcripts)
    // Using std::map because gdt::interval doesn't have a hash function
    std::map<gdt::interval, key_ptr> exon_keys;

    // Process each transcript: create exon chains (no segments yet)
    for (auto& [transcript_id, entries] : transcripts) {
        process_transcript(grove, transcript_id, entries, exon_keys, sample_id);
    }

    // TODO: Annotate all exons with overlapping features (CDS, UTR, codons)
    // annotate_exons(grove, exon_keys, annotations);
}

void build_gff::process_transcript(
    grove_type& grove,
    const std::string& transcript_id,
    const std::vector<gio::gff_entry>& all_entries,
    std::map<gdt::interval, key_ptr>& exon_keys,
    const std::string& sample_id
) {
    // Extract only exons from all entries
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

    // Step 1: Insert or reuse exons
    for (const auto& exon_entry : exons) {
        key_ptr exon_key;

        // Check if this exon already exists (same interval)
        auto it = exon_keys.find(exon_entry.interval);
        if (it != exon_keys.end()) {
            // Reuse existing exon, add this transcript to it
            exon_key = it->second;
            auto& feature_data = exon_key->get_data();
            if (is_exon(feature_data)) {
                auto& exon = get_exon(feature_data);
                exon.transcript_ids.insert(transcript_id);
                // Pan-transcriptome: add sample if provided
                if (!sample_id.empty()) {
                    exon.add_sample(sample_id);
                }
            }
        } else {
            // Create new exon feature from GFF entry
            exon_feature exon_feat = exon_feature::from_gff_entry(
                exon_entry.attributes,
                exon_entry.seqid,
                exon_entry.interval,
                exon_entry.strand.value_or('+')
            );
            exon_feat.transcript_ids.insert(transcript_id);

            // Pan-transcriptome: add sample if provided
            if (!sample_id.empty()) {
                exon_feat.add_sample(sample_id);
            }

            // Wrap in variant and insert into grove
            genomic_feature feature = exon_feat;
            exon_key = grove.insert_data(exon_entry.seqid, exon_entry.interval, feature);
            exon_keys[exon_entry.interval] = exon_key;
        }

        transcript_exon_keys.push_back(exon_key);
    }

    // Step 2: Link exons with edges (create exon chain for this transcript)
    for (size_t i = 0; i < transcript_exon_keys.size() - 1; ++i) {
        edge_metadata edge_meta(transcript_id, edge_metadata::edge_type::EXON_TO_EXON);
        grove.add_edge(transcript_exon_keys[i], transcript_exon_keys[i + 1], edge_meta);
    }

    // TODO: Step 3 - Create segments later
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

std::optional<std::string> build_gff::extract_gene_id(
    const std::map<std::string, std::string>& attributes) {
    auto it = attributes.find("gene_id");
    if (it != attributes.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<std::string> build_gff::extract_transcript_id(
    const std::map<std::string, std::string>& attributes) {
    auto it = attributes.find("transcript_id");
    if (it != attributes.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<std::string> build_gff::extract_attribute(
    const std::map<std::string, std::string>& attributes,
    const std::string& key
) {
    auto it = attributes.find(key);
    if (it != attributes.end()) {
        return it->second;
    }
    return std::nullopt;
}