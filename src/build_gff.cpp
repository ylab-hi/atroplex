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
#include <mutex>

// class
#include "utility.hpp"
#include "build_gff.hpp"

// genogrove

void build_gff::build(grove_type& grove,
    const std::filesystem::path& filepath,
    std::optional<uint32_t> sample_id,
    chromosome_exon_caches& exon_caches,
    chromosome_segment_caches& segment_caches,
    size_t& segment_count,
    uint32_t /*num_threads*/) {

    gio::gff_reader reader(filepath.string());

    // Dummy mutex (unused in single-threaded mode, but required by process_gene signature)
    std::mutex grove_mutex;

    // Buffer for current gene's entries
    std::vector<gio::gff_entry> current_gene_entries;
    std::string current_gene_id;
    std::string current_chrom;
    size_t line_count = 0;

    std::string progress_prefix = "Processing " + filepath.filename().string();
    logging::progress_start();

    for (const auto& entry : reader) {
        line_count++;

        if (line_count % 50000 == 0) {
            logging::progress(line_count, progress_prefix);
        }

        std::optional<std::string> gene_id = entry.get_gene_id();
        if (!gene_id.has_value()) {
            continue;
        }

        std::string seqid = normalize_chromosome(entry.seqid);

        // Gene changed - process previous gene
        if (!current_gene_id.empty() && gene_id.value() != current_gene_id) {
            process_gene(grove, grove_mutex, current_gene_entries,
                exon_caches[current_chrom], segment_caches[current_chrom],
                sample_id, segment_count);
            current_gene_entries.clear();
        }

        current_gene_id = gene_id.value();
        current_chrom = seqid;
        current_gene_entries.push_back(entry);
    }

    // Process final gene
    if (!current_gene_entries.empty()) {
        process_gene(grove, grove_mutex, current_gene_entries,
            exon_caches[current_chrom], segment_caches[current_chrom],
            sample_id, segment_count);
    }

    logging::progress_done(segment_count, "Processed " + filepath.filename().string());
}

void build_gff::process_gene(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::vector<gio::gff_entry>& gene_entries,
    exon_cache_type& exon_cache,
    segment_cache_type& segment_cache,
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

    // Process each transcript: create exon chains and segments
    // Uses chromosome-level caches for cross-file deduplication
    for (auto& [transcript_id, entries] : transcripts) {
        process_transcript(grove, grove_mutex, transcript_id, entries, exon_cache,
            segment_cache, sample_id, segment_count);
    }

    // TODO: Annotate all exons with overlapping features (CDS, UTR, codons)
    // annotate_exons(grove, exon_cache, annotations);
}

void build_gff::process_transcript(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::string& transcript_id,
    const std::vector<gio::gff_entry>& transcript_entries,
    std::map<gdt::genomic_coordinate, key_ptr>& exon_cache,
    std::unordered_map<std::string, key_ptr>& segment_cache,
    std::optional<uint32_t> sample_id,
    size_t& segment_count
) {
    // Step 1: Extract and sort exons in 5'→3' biological order
    std::vector<gio::gff_entry> sorted_exons = extract_sorted_exons(transcript_entries);
    if (sorted_exons.empty()) {
        return;
    }

    // Step 2: Insert or reuse exons, collecting keys and coordinates
    std::vector<key_ptr> exon_chain;
    std::vector<gdt::genomic_coordinate> exon_coords;

    for (const auto& exon_entry : sorted_exons) {
        key_ptr exon_key = insert_exon(
            grove, grove_mutex, exon_entry, transcript_id, exon_cache, sample_id
        );
        exon_chain.push_back(exon_key);

        char strand = exon_entry.strand.value_or('+');
        exon_coords.emplace_back(strand,
            exon_entry.interval.get_start(),
            exon_entry.interval.get_end()
        );
    }

    // Step 3: Link exons into a chain
    link_exon_chain(grove, grove_mutex, exon_chain, transcript_id);

    // Step 4: Create or reuse segment
    create_segment(
        grove, grove_mutex, transcript_id, sorted_exons, exon_coords,
        exon_chain, segment_cache, sample_id, segment_count
    );
}

std::vector<gio::gff_entry> build_gff::extract_sorted_exons(
    const std::vector<gio::gff_entry>& transcript_entries
) {
    std::vector<gio::gff_entry> exons;
    for (const auto& entry : transcript_entries) {
        if (entry.type == "exon") {
            exons.push_back(entry);
        }
    }

    if (exons.empty()) {
        return exons;
    }

    // Sort in 5'→3' biological order based on strand
    char strand = exons.front().strand.value_or('+');
    std::sort(exons.begin(), exons.end(),
        [strand](const gio::gff_entry& a, const gio::gff_entry& b) {
            if (a.seqid != b.seqid) return a.seqid < b.seqid;
            if (strand == '-') {
                return a.interval.get_start() > b.interval.get_start();
            }
            return a.interval.get_start() < b.interval.get_start();
        });

    return exons;
}

key_ptr build_gff::insert_exon(
    grove_type& grove,
    std::mutex& grove_mutex,
    const gio::gff_entry& exon_entry,
    const std::string& transcript_id,
    std::map<gdt::genomic_coordinate, key_ptr>& exon_cache,
    std::optional<uint32_t> sample_id
) {
    char strand = exon_entry.strand.value_or('+');
    gdt::genomic_coordinate coord(
        strand,
        exon_entry.interval.get_start(),
        exon_entry.interval.get_end()
    );

    auto cached = exon_cache.find(coord);
    if (cached != exon_cache.end()) {
        // Reuse existing exon
        key_ptr exon_key = cached->second;
        auto& exon = get_exon(exon_key->get_data());
        exon.transcript_ids.insert(transcript_id);
        if (sample_id.has_value()) {
            exon.add_sample(*sample_id);
        }
        return exon_key;
    }

    // Create new exon (normalize seqid for coordinate string)
    std::string normalized_seqid = normalize_chromosome(exon_entry.seqid);
    exon_feature new_exon = exon_feature::from_gff_entry(
        exon_entry.attributes,
        normalized_seqid,
        exon_entry.interval,
        strand
    );
    new_exon.transcript_ids.insert(transcript_id);
    if (sample_id.has_value()) {
        new_exon.add_sample(*sample_id);
    }

    genomic_feature feature = new_exon;
    key_ptr exon_key;
    {
        std::lock_guard<std::mutex> lock(grove_mutex);
        exon_key = grove.add_external_key(coord, feature);
    }
    exon_cache[coord] = exon_key;

    return exon_key;
}

void build_gff::link_exon_chain(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::vector<key_ptr>& exon_chain,
    const std::string& transcript_id
) {
    std::lock_guard<std::mutex> lock(grove_mutex);
    for (size_t i = 0; i + 1 < exon_chain.size(); ++i) {
        edge_metadata edge(transcript_id, edge_metadata::edge_type::EXON_TO_EXON);
        grove.add_edge(exon_chain[i], exon_chain[i + 1], edge);
    }
}

void build_gff::create_segment(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::string& transcript_id,
    const std::vector<gio::gff_entry>& sorted_exons,
    const std::vector<gdt::genomic_coordinate>& exon_coords,
    const std::vector<key_ptr>& exon_chain,
    std::unordered_map<std::string, key_ptr>& segment_cache,
    std::optional<uint32_t> sample_id,
    size_t& segment_count
) {
    std::string seqid = normalize_chromosome(sorted_exons.front().seqid);
    std::string structure_key = make_exon_structure_key(seqid, exon_coords);

    auto cached = segment_cache.find(structure_key);
    if (cached != segment_cache.end()) {
        // Reuse existing segment
        key_ptr seg_key = cached->second;
        auto& seg = get_segment(seg_key->get_data());
        seg.transcript_ids.insert(transcript_id);
        if (sample_id.has_value()) {
            seg.add_sample(*sample_id);
        }

        std::lock_guard<std::mutex> lock(grove_mutex);
        edge_metadata edge(transcript_id, edge_metadata::edge_type::SEGMENT_TO_EXON);
        grove.add_edge(seg_key, exon_chain.front(), edge);
        return;
    }

    // Compute segment span (min/max positions across all exons)
    auto [min_it, max_it] = std::minmax_element(sorted_exons.begin(), sorted_exons.end(),
        [](const gio::gff_entry& a, const gio::gff_entry& b) {
            return a.interval.get_start() < b.interval.get_start();
        });

    char strand = sorted_exons.front().strand.value_or('+');
    gdt::genomic_coordinate segment_coord(
        strand,
        min_it->interval.get_start(),
        max_it->interval.get_end()
    );

    std::string coordinate_str = seqid + ":" + strand + ":" +
        std::to_string(segment_coord.get_start()) + "-" +
        std::to_string(segment_coord.get_end());

    // Build segment feature
    segment_feature new_segment(segment_feature::source_type::REFERENCE);
    new_segment.id = structure_key;
    new_segment.transcript_ids.insert(transcript_id);
    new_segment.exon_count = static_cast<int>(sorted_exons.size());
    new_segment.coordinate = coordinate_str;

    // Copy gene info from first exon
    const auto& first_exon_data = get_exon(exon_chain.front()->get_data());
    new_segment.gene_id = first_exon_data.gene_id;
    new_segment.gene_name = first_exon_data.gene_name;
    new_segment.gene_biotype = first_exon_data.gene_biotype;

    if (sample_id.has_value()) {
        new_segment.add_sample(*sample_id);
    }

    // Insert into grove (protected by mutex)
    genomic_feature feature = new_segment;
    key_ptr seg_key;
    {
        std::lock_guard<std::mutex> lock(grove_mutex);
        seg_key = grove.insert_data(seqid, segment_coord, feature);
        edge_metadata edge(transcript_id, edge_metadata::edge_type::SEGMENT_TO_EXON);
        grove.add_edge(seg_key, exon_chain.front(), edge);
    }
    segment_count++;

    segment_cache[structure_key] = seg_key;
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