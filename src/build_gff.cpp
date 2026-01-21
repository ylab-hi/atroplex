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
#include <fstream>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <regex>

// class
#include "utility.hpp"
#include "build_gff.hpp"

// genogrove

void build_gff::build(grove_type& grove,
    const std::filesystem::path& filepath,
    std::optional<uint32_t> sample_id) {
    gio::gff_reader reader(filepath.string());

    // Buffer for gene entries - process one gene at a time
    std::unordered_map<std::string, std::vector<gio::gff_entry>> genes;
    std::string current_gene_id;
    size_t line_count = 0;
    size_t segment_count = 0;

    std::string progress_prefix = "Processing " + filepath.filename().string();
    logging::progress_start();

    for(const auto& entry : reader) {
        line_count++;

        // Update progress every 1000 entries
        if (line_count % 1000 == 0) {
            logging::progress(line_count, progress_prefix);
        }

        // Extract gene_id to group entries by gene
        std::optional<std::string> gene_id = entry.get_gene_id();
        if (!gene_id.has_value()) {
            continue; // Skip entries without gene_id
        }

        // If we've moved to a new gene, process the previous one
        if (!current_gene_id.empty() && gene_id.value() != current_gene_id) {
            process_gene(grove, genes[current_gene_id], sample_id, segment_count);
            genes.erase(current_gene_id); // Free memory
        }

        current_gene_id = gene_id.value();
        genes[gene_id.value()].push_back(entry);
    }

    // Process last gene
    if (!current_gene_id.empty()) {
        process_gene(grove, genes[current_gene_id], sample_id, segment_count);
    }

    logging::progress_done(segment_count, "Processed " + filepath.filename().string());
}

void build_gff::process_gene(
    grove_type& grove,
    const std::vector<gio::gff_entry>& gene_entries,
    std::optional<uint32_t> sample_id,
    size_t& segment_count
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
    std::map<gdt::genomic_coordinate, key_ptr> exon_keys;

    // Map to track segments (so they can be deduplicated) - according to exon structure
    // Key: string representation of the exon chain: "chr1:+:100-200,300-400,500-600"
    std::unordered_map<std::string, key_ptr> segment_keys;

    // Process each transcript: create exon chains and segments
    for (auto& [transcript_id, entries] : transcripts) {
        process_transcript(grove, transcript_id, entries, exon_keys, sample_id, segment_count);
    }

    // TODO: Annotate all exons with overlapping features (CDS, UTR, codons)
    // annotate_exons(grove, exon_keys, annotations);
}

void build_gff::process_transcript(
    grove_type& grove,
    const std::string& transcript_id,
    const std::vector<gio::gff_entry>& all_entries,
    std::map<gdt::genomic_coordinate, key_ptr>& exon_keys,
    std::optional<uint32_t> sample_id,
    size_t& segment_count
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

    // Determine strand (use first exon's strand, default to +)
    char strand = exons.front().strand.value_or('+');

    // Sort exons in 5'→3' biological order:
    // - On + strand: ascending positions (5' = low coordinate, 3' = high coordinate)
    // - On - strand: descending positions (5' = high coordinate, 3' = low coordinate)
    std::sort(exons.begin(), exons.end(),
             [strand](const gio::gff_entry& a, const gio::gff_entry& b) {
                 if (a.seqid != b.seqid) return a.seqid < b.seqid;
                 if (strand == '-') {
                     return a.interval.get_start() > b.interval.get_start();  // Descending for - strand
                 }
                 return a.interval.get_start() < b.interval.get_start();  // Ascending for + strand
             });

    // Collect exon keys for this transcript
    std::vector<gdt::genomic_coordinate> exon_coords;
    std::vector<key_ptr> transcript_exon_keys;

    // Step 1: insert or reuse exons
    for (const auto& exon_entry : exons) {
        key_ptr exon_key;

        // Create genomic_coordinate for this exon (includes strand)
        char exon_strand = exon_entry.strand.value_or('+');
        gdt::genomic_coordinate exon_coord(
            exon_strand,
            exon_entry.interval.get_start(),
            exon_entry.interval.get_end()
        );

        // Check if this exon already exists (same coordinate including strand)
        auto it = exon_keys.find(exon_coord);
        if (it != exon_keys.end()) {
            // Reuse existing exon, add this transcript to it
            exon_key = it->second;
            auto& feature_data = exon_key->get_data();
            if (is_exon(feature_data)) {
                auto& exon = get_exon(feature_data);
                exon.transcript_ids.insert(transcript_id);
                // Pan-transcriptome: add sample if provided
                if (sample_id.has_value()) {
                    exon.add_sample(*sample_id);
                }
            }
        } else {
            // Create new exon feature from GFF entry
            exon_feature exon_feat = exon_feature::from_gff_entry(
                exon_entry.attributes,
                exon_entry.seqid,
                exon_entry.interval,
                exon_strand
            );
            exon_feat.transcript_ids.insert(transcript_id);

            // Pan-transcriptome: add sample if provided
            if (sample_id.has_value()) {
                exon_feat.add_sample(*sample_id);
            }

            // Wrap in variant and add as external key (graph-only, not spatially indexed)
            genomic_feature feature = exon_feat;
            exon_key = grove.add_external_key(exon_coord, feature);
            exon_keys[exon_coord] = exon_key;
        }

        transcript_exon_keys.push_back(exon_key);
    }

    // Step 2: Link exons with edges (create exon chain for this transcript)
    for (size_t i = 0; i < transcript_exon_keys.size() - 1; ++i) {
        edge_metadata edge_meta(transcript_id, edge_metadata::edge_type::EXON_TO_EXON);
        grove.add_edge(transcript_exon_keys[i], transcript_exon_keys[i + 1], edge_meta);
    }

    // Step 3: Create segment (spatially indexed) spanning transcript extent
    // Use min/max of all exon positions (not first/last, since exons are sorted by biological order)
    auto [min_it, max_it] = std::minmax_element(exons.begin(), exons.end(),
        [](const gio::gff_entry& a, const gio::gff_entry& b) {
            return a.interval.get_start() < b.interval.get_start();
        });

    // Create genomic_coordinate for the segment (includes strand)
    gdt::genomic_coordinate segment_coord(
        strand,
        min_it->interval.get_start(),
        max_it->interval.get_end()
    );

    // Build coordinate string (strand already defined above)
    const auto& first_exon = exons.front();  // First in biological order (5' end)
    std::string coordinate = first_exon.seqid + ":" + strand + ":" +
                             std::to_string(segment_coord.get_start()) + "-" +
                             std::to_string(segment_coord.get_end());

    // Create segment feature
    segment_feature seg_feat(segment_feature::source_type::REFERENCE);
    seg_feat.id = transcript_id;
    seg_feat.transcript_ids.insert(transcript_id);
    seg_feat.exon_count = static_cast<int>(exons.size());
    seg_feat.coordinate = coordinate;

    // Copy gene info from first exon
    const auto& first_exon_key = transcript_exon_keys.front();
    if (is_exon(first_exon_key->get_data())) {
        const auto& exon_data = get_exon(first_exon_key->get_data());
        seg_feat.gene_id = exon_data.gene_id;
        seg_feat.gene_name = exon_data.gene_name;
        seg_feat.gene_biotype = exon_data.gene_biotype;
    }

    // Pan-transcriptome: add sample if provided
    if (sample_id.has_value()) {
        seg_feat.add_sample(*sample_id);
    }

    // Insert segment (spatially indexed for transcript-level queries)
    genomic_feature seg_feature = seg_feat;
    key_ptr seg_key = grove.insert_data(first_exon.seqid, segment_coord, seg_feature);
    segment_count++;

    // Step 4: Link segment to first exon
    edge_metadata seg_edge(transcript_id, edge_metadata::edge_type::SEGMENT_TO_EXON);
    grove.add_edge(seg_key, transcript_exon_keys.front(), seg_edge);
}

void build_gff::annotate_exons(
    grove_type& grove,
    const std::map<gdt::genomic_coordinate, key_ptr>& exon_keys,
    const std::vector<gio::gff_entry>& annotations
) {
    // For each annotation feature (CDS, UTR, codon)
    for (const auto& annot : annotations) {
        // Find overlapping exons
        for (const auto& [exon_coord, exon_key] : exon_keys) {
            // Check for overlap (ignoring strand for now - annotations should match exon strand)
            if (exon_coord.get_end() < annot.interval.get_start() ||
                exon_coord.get_start() > annot.interval.get_end()) {
                continue; // No overlap
            }

            // Compute overlapping interval
            gdt::interval overlap;
            overlap.set_start(std::max(exon_coord.get_start(), annot.interval.get_start()));
            overlap.set_end(std::min(exon_coord.get_end(), annot.interval.get_end()));

            // Add to exon's overlapping features using public getter
            auto& feature_data = exon_key->get_data();
            if (is_exon(feature_data)) {
                get_exon(feature_data).overlapping_features[annot.type].push_back(overlap);
            }
        }
    }
}

std::string build_gff::make_exon_structure_key(
    const std::string& seqid,
    const std::vector<gdt::genomic_coordinate>& exon_coords
) {
    if (exon_coords.empty()) {
        return "";
    }

    std::string key;
    key.reserve(seqid.size() + 3 + exon_coords.size() * 24);  // Pre-allocate estimate

    key += seqid;
    key += ':';
    key += exon_coords.front().get_strand();
    key += ':';

    for (size_t i = 0; i < exon_coords.size(); ++i) {
        if (i > 0) {
            key += ',';
        }
        key += std::to_string(exon_coords[i].get_start());
        key += '-';
        key += std::to_string(exon_coords[i].get_end());
    }

    return key;
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

sample_info build_gff::parse_header(const std::filesystem::path& filepath) {
    sample_info info;
    info.source_file = filepath;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        logging::warning("Could not open file for header parsing: " + filepath.string());
        info.id = filepath.stem().string();
        return info;
    }

    std::string line;
    bool found_provider = false;

    // Read header lines (start with #)
    while (std::getline(file, line)) {
        // Stop at first non-header line
        if (line.empty() || line[0] != '#') {
            break;
        }

        // Skip ##gff-version and ##sequence-region directives
        if (line.starts_with("##gff-version") || line.starts_with("##sequence-region")) {
            continue;
        }

        // GENCODE format: #key: value (single #)
        // Example:
        //   #description: evidence-based annotation of the human genome (GRCh38), version 49
        //   #provider: GENCODE
        //   #date: 2025-07-08
        if (line.starts_with("#") && !line.starts_with("##") && !line.starts_with("#!")) {
            // Remove leading #
            std::string content = line.substr(1);

            // Find colon separator
            size_t colon_pos = content.find(':');
            if (colon_pos != std::string::npos) {
                std::string key = content.substr(0, colon_pos);
                std::string value = content.substr(colon_pos + 1);

                // Trim whitespace from value
                size_t start = value.find_first_not_of(" \t");
                if (start != std::string::npos) {
                    value = value.substr(start);
                }

                if (key == "provider") {
                    info.annotation_source = value;
                    found_provider = true;
                } else if (key == "description") {
                    // Parse version from description
                    // "evidence-based annotation of the human genome (GRCh38), version 49 (Ensembl 115)"
                    std::regex version_regex(R"(version\s+(\d+))");
                    std::smatch match;
                    if (std::regex_search(value, match, version_regex)) {
                        info.annotation_version = "v" + match[1].str();
                    }

                    // Parse genome build from description (first parenthetical)
                    std::regex genome_regex(R"(\(([^)]+)\))");
                    if (std::regex_search(value, match, genome_regex)) {
                        info.with_attribute("genome_build", match[1].str());
                    }

                    // Store full description
                    info.with_attribute("description", value);
                } else if (key == "date") {
                    info.with_attribute("date", value);
                } else if (key == "format") {
                    info.with_attribute("format", value);
                } else if (key == "contact") {
                    info.with_attribute("contact", value);
                }
            }
        }
        // Ensembl format: #!key value
        else if (line.starts_with("#!")) {
            std::string content = line.substr(2);

            // Find first space separator
            size_t space_pos = content.find(' ');
            if (space_pos != std::string::npos) {
                std::string key = content.substr(0, space_pos);
                std::string value = content.substr(space_pos + 1);

                // Trim whitespace from value
                size_t start = value.find_first_not_of(" \t");
                if (start != std::string::npos) {
                    value = value.substr(start);
                }

                if (key == "genome-build") {
                    info.with_attribute("genome_build", value);
                } else if (key == "genome-version") {
                    info.with_attribute("genome_version", value);
                } else if (key == "genebuild-last-updated") {
                    info.with_attribute("date", value);
                } else if (key == "genome-build-accession") {
                    info.with_attribute("genome_accession", value);
                }

                // If we see Ensembl-style headers, set provider
                if (!found_provider) {
                    info.annotation_source = "Ensembl";
                    found_provider = true;
                }
            }
        }
    }

    // Build sample ID from metadata
    if (!info.annotation_source.empty()) {
        info.id = info.annotation_source;
        if (!info.annotation_version.empty()) {
            info.id += "_" + info.annotation_version;
        }
        // Add genome build if available
        std::string genome_build = info.get_attribute("genome_build");
        if (!genome_build.empty()) {
            info.id += "_" + genome_build;
        }
    } else {
        // Fallback to filename
        info.id = filepath.stem().string();
    }

    logging::info("Parsed header metadata: id=" + info.id +
                  ", source=" + info.annotation_source +
                  ", version=" + info.annotation_version);

    return info;
}