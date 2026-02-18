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
    chromosome_gene_segment_indices& gene_indices,
    size_t& segment_count,
    uint32_t /*num_threads*/,
    float min_expression) {

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
                gene_indices[current_chrom],
                sample_id, segment_count, min_expression);
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
            gene_indices[current_chrom],
            sample_id, segment_count, min_expression);
    }

    logging::progress_done(segment_count, "Processed " + filepath.filename().string());
}

void build_gff::process_gene(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::vector<gio::gff_entry>& gene_entries,
    exon_cache_type& exon_cache,
    segment_cache_type& segment_cache,
    gene_segment_index_type& gene_index,
    std::optional<uint32_t> sample_id,
    size_t& segment_count,
    float min_expression
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
            segment_cache, gene_index, sample_id, segment_count, min_expression);
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
    gene_segment_index_type& gene_index,
    std::optional<uint32_t> sample_id,
    size_t& segment_count,
    float min_expression
) {
    // Step 1: Extract and sort exons in 5'→3' biological order
    std::vector<gio::gff_entry> sorted_exons = extract_sorted_exons(transcript_entries);
    if (sorted_exons.empty()) {
        return;
    }

    // Extract expression and transcript biotype from transcript-level entry
    float expression_value = -1.0f;
    std::string transcript_biotype;
    for (const auto& entry : transcript_entries) {
        if (entry.type == "transcript") {
            // Auto-detect expression attribute (priority order)
            static const std::pair<const char*, sample_info::expression_type> expr_attrs[] = {
                {"counts", sample_info::expression_type::COUNTS},
                {"TPM",    sample_info::expression_type::TPM},
                {"FPKM",   sample_info::expression_type::FPKM},
                {"RPKM",   sample_info::expression_type::RPKM},
                {"cov",    sample_info::expression_type::UNKNOWN},
            };
            for (const auto& [attr_name, expr_type] : expr_attrs) {
                auto it = entry.attributes.find(attr_name);
                if (it != entry.attributes.end()) {
                    try {
                        expression_value = std::stof(it->second);
                        // Set expr_type on sample_info if not already set
                        if (sample_id.has_value()) {
                            auto& registry = sample_registry::instance();
                            auto* info = registry.get(*sample_id);
                            if (info && !info->has_expression_type()) {
                                info->expr_type = expr_type;
                            }
                        }
                    } catch (...) {}
                    break;
                }
            }

            // Extract transcript biotype (GENCODE: transcript_type, Ensembl: transcript_biotype)
            auto bt_it = entry.attributes.find("transcript_type");
            if (bt_it == entry.attributes.end()) {
                bt_it = entry.attributes.find("transcript_biotype");
            }
            if (bt_it != entry.attributes.end()) {
                transcript_biotype = bt_it->second;
            }

            break;
        }
    }

    // Filter: skip transcripts below expression threshold
    // Only applies when both the filter is active (min_expression >= 0) and the
    // transcript has expression data (expression_value >= 0)
    if (min_expression >= 0 && expression_value >= 0 && expression_value < min_expression) {
        return;
    }

    // Step 2: Insert or reuse exons, collecting keys and coordinates
    std::vector<key_ptr> exon_chain;
    std::vector<gdt::genomic_coordinate> exon_coords;

    // Extract GFF source from first exon (column 2 - e.g., HAVANA, ENSEMBL, StringTie)
    std::string gff_source = sorted_exons.front().source;

    for (const auto& exon_entry : sorted_exons) {
        key_ptr exon_key = insert_exon(
            grove, grove_mutex, exon_entry, transcript_id, exon_cache, sample_id, gff_source
        );
        exon_chain.push_back(exon_key);

        char strand = exon_entry.strand.value_or('+');
        exon_coords.emplace_back(strand,
            exon_entry.interval.get_start(),
            exon_entry.interval.get_end()
        );
    }

    // Step 3: Create or reuse segment (edges created inside for new segments only)
    create_segment(
        grove, grove_mutex, transcript_id, sorted_exons, exon_coords,
        exon_chain, segment_cache, gene_index, sample_id, gff_source,
        segment_count, expression_value, transcript_biotype
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
    std::optional<uint32_t> sample_id,
    const std::string& gff_source
) {
    char strand = exon_entry.strand.value_or('+');
    gdt::genomic_coordinate coord(
        strand,
        exon_entry.interval.get_start(),
        exon_entry.interval.get_end()
    );

    auto cached = exon_cache.find(coord);
    if (cached != exon_cache.end()) {
        // Reuse existing exon - add transcript, sample, and source
        key_ptr exon_key = cached->second;
        auto& exon = get_exon(exon_key->get_data());
        exon.transcript_ids.insert(transcript_registry::instance().intern(transcript_id));
        if (sample_id.has_value()) {
            exon.add_sample(*sample_id);
        }
        if (!gff_source.empty()) {
            exon.add_source(gff_source);
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
    new_exon.transcript_ids.insert(transcript_registry::instance().intern(transcript_id));
    if (sample_id.has_value()) {
        new_exon.add_sample(*sample_id);
    }
    if (!gff_source.empty()) {
        new_exon.add_source(gff_source);
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

void build_gff::create_segment(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::string& transcript_id,
    const std::vector<gio::gff_entry>& sorted_exons,
    const std::vector<gdt::genomic_coordinate>& exon_coords,
    const std::vector<key_ptr>& exon_chain,
    std::unordered_map<std::string, key_ptr>& segment_cache,
    gene_segment_index_type& gene_index,
    std::optional<uint32_t> sample_id,
    const std::string& gff_source,
    size_t& segment_count,
    float expression_value,
    const std::string& transcript_biotype
) {
    std::string seqid = normalize_chromosome(sorted_exons.front().seqid);
    std::string structure_key = make_exon_structure_key(seqid, exon_coords);

    // Step 1: Exact match deduplication (unchanged)
    auto cached = segment_cache.find(structure_key);
    if (cached != segment_cache.end()) {
        merge_into_segment(cached->second, transcript_id, sample_id,
                          gff_source, expression_value, transcript_biotype);
        return;
    }

    // Step 2: Skip mono-exon transcripts (ambiguous parent, noise in long-read data)
    if (exon_chain.size() == 1) {
        return;
    }

    // Step 3: Forward absorption — check if a longer parent segment exists
    const std::string& gene_id = get_exon(exon_chain.front()->get_data()).gene_id();
    auto gene_it = gene_index.find(gene_id);

    if (gene_it != gene_index.end()) {
        key_ptr best_parent = nullptr;
        size_t best_exon_count = 0;

        for (const auto& entry : gene_it->second) {
            auto& candidate_seg = get_segment(entry.segment->get_data());
            if (candidate_seg.absorbed) continue;
            if (entry.exon_chain.size() <= exon_chain.size()) continue;

            if (is_contiguous_subsequence(exon_chain, entry.exon_chain)) {
                if (entry.exon_chain.size() > best_exon_count) {
                    best_parent = entry.segment;
                    best_exon_count = entry.exon_chain.size();
                }
            }
        }

        if (best_parent != nullptr) {
            // Absorb: merge this ISM's metadata into the parent segment
            merge_into_segment(best_parent, transcript_id, sample_id,
                              gff_source, expression_value, transcript_biotype);
            get_segment(best_parent->get_data()).absorbed_count++;
            return;
        }
    }

    // Step 4: Create new segment (no exact match, no parent found)
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

    segment_feature new_segment;
    new_segment.segment_index = segment_count;
    uint32_t tx_id = transcript_registry::instance().intern(transcript_id);
    new_segment.transcript_ids.insert(tx_id);
    if (!transcript_biotype.empty()) {
        new_segment.transcript_biotypes[tx_id] = transcript_biotype;
    }
    new_segment.exon_count = static_cast<int>(sorted_exons.size());

    const auto& first_exon_data = get_exon(exon_chain.front()->get_data());
    new_segment.gene_idx = first_exon_data.gene_idx;

    if (sample_id.has_value()) {
        new_segment.add_sample(*sample_id);
        if (expression_value >= 0.0f) {
            new_segment.set_expression(*sample_id, expression_value);
        }
    }
    if (!gff_source.empty()) {
        new_segment.add_source(gff_source);
    }

    genomic_feature feature = new_segment;
    key_ptr seg_key;
    {
        std::lock_guard<std::mutex> lock(grove_mutex);
        seg_key = grove.insert_data(seqid, segment_coord, feature);

        edge_metadata seg_edge(segment_count, edge_metadata::edge_type::SEGMENT_TO_EXON);
        grove.add_edge(seg_key, exon_chain.front(), seg_edge);

        for (size_t i = 0; i + 1 < exon_chain.size(); ++i) {
            edge_metadata exon_edge(segment_count, edge_metadata::edge_type::EXON_TO_EXON);
            grove.add_edge(exon_chain[i], exon_chain[i + 1], exon_edge);
        }
    }
    segment_count++;

    segment_cache[structure_key] = seg_key;

    // Step 5: Register in gene segment index
    gene_index[gene_id].push_back({seg_key, exon_chain, structure_key});

    // Step 6: Reverse absorption — absorb existing shorter segments into this new one
    try_reverse_absorption(gene_index, gene_id, seg_key, exon_chain, segment_cache);
}

bool build_gff::is_contiguous_subsequence(
    const std::vector<key_ptr>& sub,
    const std::vector<key_ptr>& parent
) {
    if (sub.empty() || sub.size() >= parent.size()) return false;

    for (size_t start = 0; start <= parent.size() - sub.size(); ++start) {
        if (parent[start] == sub[0]) {
            bool all_match = true;
            for (size_t j = 1; j < sub.size(); ++j) {
                if (parent[start + j] != sub[j]) {
                    all_match = false;
                    break;
                }
            }
            if (all_match) return true;
        }
    }
    return false;
}

void build_gff::merge_into_segment(
    key_ptr target_seg,
    const std::string& transcript_id,
    std::optional<uint32_t> sample_id,
    const std::string& gff_source,
    float expression_value,
    const std::string& transcript_biotype
) {
    auto& seg = get_segment(target_seg->get_data());
    uint32_t tx_id = transcript_registry::instance().intern(transcript_id);
    seg.transcript_ids.insert(tx_id);
    if (!transcript_biotype.empty()) {
        seg.transcript_biotypes[tx_id] = transcript_biotype;
    }
    if (sample_id.has_value()) {
        seg.add_sample(*sample_id);
        if (expression_value >= 0.0f) {
            seg.expression.accumulate(*sample_id, expression_value);
        }
    }
    if (!gff_source.empty()) {
        seg.add_source(gff_source);
    }
}

void build_gff::try_reverse_absorption(
    gene_segment_index_type& gene_index,
    const std::string& gene_id,
    key_ptr new_seg,
    const std::vector<key_ptr>& new_exon_chain,
    segment_cache_type& segment_cache
) {
    auto gene_it = gene_index.find(gene_id);
    if (gene_it == gene_index.end()) return;

    auto& parent_seg = get_segment(new_seg->get_data());

    for (auto& entry : gene_it->second) {
        if (entry.segment == new_seg) continue;

        auto& candidate_seg = get_segment(entry.segment->get_data());
        if (candidate_seg.absorbed) continue;
        if (entry.exon_chain.size() >= new_exon_chain.size()) continue;

        if (is_contiguous_subsequence(entry.exon_chain, new_exon_chain)) {
            // Merge candidate's metadata into new parent
            for (const auto& tx : candidate_seg.transcript_ids)
                parent_seg.transcript_ids.insert(tx);
            for (const auto& [tx, bt] : candidate_seg.transcript_biotypes)
                parent_seg.transcript_biotypes[tx] = bt;
            parent_seg.sample_idx.merge(candidate_seg.sample_idx);
            parent_seg.expression.merge(candidate_seg.expression);
            parent_seg.merge_sources(candidate_seg.sources);
            parent_seg.absorbed_count += candidate_seg.absorbed_count + 1;

            // Tombstone the absorbed segment
            candidate_seg.absorbed = true;

            // Remove from segment_cache so no future exact matches find it
            segment_cache.erase(entry.structure_key);
        }
    }
}

void build_gff::annotate_exons(
    grove_type& /*grove*/,
    const std::map<gdt::genomic_coordinate, key_ptr>& /*exon_keys*/,
    const std::vector<gio::gff_entry>& /*annotations*/
) {
    // TODO: CDS/UTR overlap annotation not yet implemented.
    // Requires adding an appropriate storage field to exon_feature.
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

    // Helper to convert string to lowercase
    auto to_lower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(),
            [](unsigned char c) { return std::tolower(c); });
        return s;
    };

    // Helper to trim whitespace
    auto trim = [](const std::string& s) -> std::string {
        size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        size_t end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    };

    std::string line;
    while (std::getline(file, line)) {
        // Stop at first non-comment line
        if (line.empty() || line[0] != '#') {
            break;
        }

        // Must be ## directive (GFF convention)
        if (line.size() < 2 || line[1] != '#') {
            continue;
        }

        // Parse ##property: value format
        std::string content = line.substr(2);  // Remove leading ##
        size_t colon_pos = content.find(':');
        if (colon_pos == std::string::npos) {
            continue;
        }

        std::string key = to_lower(trim(content.substr(0, colon_pos)));
        std::string value = trim(content.substr(colon_pos + 1));

        if (value.empty()) {
            continue;
        }

        // Map to sample_info fields
        if (key == "id") {
            info.id = value;
        } else if (key == "description") {
            info.description = value;
        } else if (key == "assay") {
            info.assay = value;
        } else if (key == "biosample_type") {
            info.biosample_type = value;
        } else if (key == "biosample") {
            info.biosample = value;
        } else if (key == "condition") {
            info.condition = value;
        } else if (key == "treatment") {
            info.treatment = value;
        } else if (key == "species") {
            info.species = value;
        } else if (key == "replication_type") {
            info.replication_type = value;
        } else if (key == "platform") {
            info.platform = value;
        } else if (key == "pipeline") {
            info.pipeline = value;
        } else if (key == "pipeline_version") {
            info.pipeline_version = value;
        } else if (key == "annotation_source") {
            info.annotation_source = value;
        } else if (key == "annotation_version") {
            info.annotation_version = value;
        } else if (key == "source_url") {
            info.source_url = value;
        } else if (key == "publication") {
            info.publication = value;
        } else {
            // Unknown properties go to attributes
            info.attributes[key] = value;
        }
    }

    // Default ID to filename if not provided
    if (info.id.empty()) {
        info.id = filepath.stem().string();
    }

    return info;
}