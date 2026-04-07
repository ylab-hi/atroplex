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
    float min_expression,
    bool absorb,
    size_t fuzzy_tolerance) {

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
                sample_id, segment_count, min_expression, absorb, fuzzy_tolerance);
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
            sample_id, segment_count, min_expression, absorb, fuzzy_tolerance);
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
    float min_expression,
    bool absorb,
    size_t fuzzy_tolerance
) {
    // Group entries by transcript
    std::unordered_map<std::string, std::vector<gio::gff_entry>> transcripts;
    for (const auto& entry : gene_entries) {
        std::optional<std::string> transcript_id = entry.get_transcript_id();
        if (transcript_id.has_value()) {
            transcripts[transcript_id.value()].push_back(entry);
        }
    }

    // Sort transcript keys by exon count (descending) so multi-exon transcripts
    // are processed first — ensures mono-exon fragments have a parent to check against.
    auto count_exons = [](const std::vector<gio::gff_entry>& entries) {
        return std::count_if(entries.begin(), entries.end(),
            [](const gio::gff_entry& e) { return e.type == "exon"; });
    };

    std::vector<std::string> tx_order;
    tx_order.reserve(transcripts.size());
    for (const auto& [tx_id, _] : transcripts) {
        tx_order.push_back(tx_id);
    }
    std::sort(tx_order.begin(), tx_order.end(),
        [&](const auto& a, const auto& b) {
            return count_exons(transcripts[a]) > count_exons(transcripts[b]);
        });

    for (const auto& transcript_id : tx_order) {
        process_transcript(grove, grove_mutex, transcript_id, transcripts[transcript_id],
            exon_cache, segment_cache, gene_index, sample_id, segment_count,
            min_expression, absorb, fuzzy_tolerance);
    }
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
    float min_expression,
    bool absorb,
    size_t fuzzy_tolerance
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
                        if (sample_id.has_value()) {
                            auto& registry = sample_registry::instance();
                            if (registry.contains(*sample_id)) {
                                auto& info = registry.get(*sample_id);
                                if (!info.has_expression_type()) {
                                    info.expr_type = expr_type;
                                }
                            }
                        }
                    } catch (...) {}
                    break;
                }
            }

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
    if (min_expression >= 0 && expression_value >= 0 && expression_value < min_expression) {
        return;
    }

    // Step 2: Insert or reuse exons, collecting keys and coordinates
    std::vector<key_ptr> exon_chain;
    std::vector<gdt::genomic_coordinate> exon_coords;
    std::string gff_source = sorted_exons.front().source;

    for (const auto& exon_entry : sorted_exons) {
        key_ptr exon_key = insert_exon(
            grove, grove_mutex, exon_entry, transcript_id, exon_cache, sample_id, gff_source
        );
        exon_chain.push_back(exon_key);

        char strand = exon_entry.strand.value_or('+');
        exon_coords.emplace_back(strand, exon_entry.start, exon_entry.end);
    }

    // Step 3: Create or reuse segment via segment_builder (format-agnostic)
    std::string seqid = normalize_chromosome(sorted_exons.front().seqid);
    char strand = sorted_exons.front().strand.value_or('+');
    auto [min_it, max_it] = std::minmax_element(sorted_exons.begin(), sorted_exons.end(),
        [](const gio::gff_entry& a, const gio::gff_entry& b) {
            return a.start < b.start;
        });

    segment_builder::create_segment(
        grove, grove_mutex, transcript_id, seqid, strand,
        min_it->start, max_it->end, static_cast<int>(sorted_exons.size()),
        exon_coords, exon_chain, segment_cache, gene_index, sample_id,
        gff_source, segment_count, expression_value, transcript_biotype, absorb, fuzzy_tolerance
    );
}

// ========== GFF-specific helpers ==========

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

    char strand = exons.front().strand.value_or('+');
    std::sort(exons.begin(), exons.end(),
        [strand](const gio::gff_entry& a, const gio::gff_entry& b) {
            if (a.seqid != b.seqid) return a.seqid < b.seqid;
            if (strand == '-') {
                return a.start > b.start;
            }
            return a.start < b.start;
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
    gdt::genomic_coordinate coord(strand, exon_entry.start, exon_entry.end);

    auto cached = exon_cache.find(coord);
    if (cached != exon_cache.end()) {
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

    std::string normalized_seqid = normalize_chromosome(exon_entry.seqid);
    exon_feature new_exon = exon_feature::from_gff_entry(
        exon_entry.attributes,
        normalized_seqid,
        gdt::interval(exon_entry.start, exon_entry.end),
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

void build_gff::annotate_exons(
    grove_type& /*grove*/,
    const std::map<gdt::genomic_coordinate, key_ptr>& /*exon_keys*/,
    const std::vector<gio::gff_entry>& /*annotations*/
) {
    // TODO: CDS/UTR overlap annotation not yet implemented.
}

// ========== Attribute extraction helpers ==========

std::optional<std::string> build_gff::extract_gene_id(
    const std::map<std::string, std::string, std::less<>>& attributes) {
    auto it = attributes.find("gene_id");
    if (it != attributes.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<std::string> build_gff::extract_transcript_id(
    const std::map<std::string, std::string, std::less<>>& attributes) {
    auto it = attributes.find("transcript_id");
    if (it != attributes.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<std::string> build_gff::extract_attribute(
    const std::map<std::string, std::string, std::less<>>& attributes,
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

    auto to_lower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(),
            [](unsigned char c) { return std::tolower(c); });
        return s;
    };

    auto trim = [](const std::string& s) -> std::string {
        size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        size_t end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    };

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] != '#') {
            break;
        }

        if (line.size() < 2 || line[1] != '#') {
            continue;
        }

        std::string content = line.substr(2);
        size_t colon_pos = content.find(':');
        if (colon_pos == std::string::npos) {
            continue;
        }

        std::string key = to_lower(trim(content.substr(0, colon_pos)));
        std::string value = trim(content.substr(colon_pos + 1));

        if (value.empty()) {
            continue;
        }

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
            info.attributes[key] = value;
        }
    }

    if (info.id.empty()) {
        info.id = filepath.stem().string();
    }

    return info;
}