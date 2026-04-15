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
    const expression_filters& filters,
    bool absorb,
    size_t fuzzy_tolerance,
    bool include_scaffolds,
    build_counters& counters) {

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

        // Scaffold filter: skip entries on unplaced contigs, alt haplotypes,
        // fix patches, decoys, etc. Default is main chromosomes only
        // (chr1..chr22, chrX, chrY, chrM); --include-scaffolds disables
        // the filter for non-human/non-mouse use cases. We only count the
        // transcript-level entries so `scaffold_filtered_transcripts`
        // mirrors the `input_transcripts` semantics.
        if (!is_main_chromosome(entry.seqid, include_scaffolds)) {
            if (entry.type == "transcript") {
                counters.scaffold_filtered_transcripts++;
            }
            continue;
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
                sample_id, segment_count, filters, absorb, fuzzy_tolerance, counters);
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
            sample_id, segment_count, filters, absorb, fuzzy_tolerance, counters);
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
    const expression_filters& filters,
    bool absorb,
    size_t fuzzy_tolerance,
    build_counters& counters
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
    // Precompute exon counts once per transcript to avoid O(T·E·log T) rescans
    // inside the sort comparator.
    std::vector<std::pair<std::string, size_t>> tx_exon_counts;
    tx_exon_counts.reserve(transcripts.size());
    for (const auto& [tx_id, entries] : transcripts) {
        size_t n = 0;
        for (const auto& e : entries) {
            if (e.type == "exon") ++n;
        }
        tx_exon_counts.emplace_back(tx_id, n);
    }
    std::sort(tx_exon_counts.begin(), tx_exon_counts.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });

    std::vector<std::string> tx_order;
    tx_order.reserve(tx_exon_counts.size());
    for (auto& [tx_id, _] : tx_exon_counts) {
        tx_order.push_back(std::move(tx_id));
    }

    for (const auto& transcript_id : tx_order) {
        process_transcript(grove, grove_mutex, transcript_id, transcripts[transcript_id],
            exon_cache, segment_cache, gene_index, sample_id, segment_count,
            filters, absorb, fuzzy_tolerance, counters);
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
    const expression_filters& filters,
    bool absorb,
    size_t fuzzy_tolerance,
    build_counters& counters
) {
    // Step 1: Extract and sort exons in 5'→3' biological order
    std::vector<gio::gff_entry> sorted_exons = extract_sorted_exons(transcript_entries);
    if (sorted_exons.empty()) {
        return;
    }
    counters.input_transcripts++;

    // Resolve the sample's declared expression_attributes (set from the
    // manifest's optional `expression_attribute` column). Empty vector
    // means no filtering and no expression storage for this sample.
    const std::vector<std::string>* declared_attrs = nullptr;
    if (sample_id.has_value()) {
        auto& registry = sample_registry::instance();
        if (registry.contains(*sample_id)) {
            declared_attrs = &registry.get(*sample_id).expression_attributes;
        }
    }

    // Extract transcript biotype and any declared-attribute values from
    // the transcript-level entry. expression_value = the FIRST declared
    // attribute's value (used for storage on the segment); the filter
    // loop below walks ALL declared attributes and evaluates each one
    // against its corresponding threshold (AND semantics).
    float expression_value = -1.0f;
    bool drop_by_filter = false;
    std::string transcript_biotype;
    for (const auto& entry : transcript_entries) {
        if (entry.type != "transcript") continue;

        if (declared_attrs && !declared_attrs->empty()) {
            bool first = true;
            for (const auto& attr : *declared_attrs) {
                auto it = entry.attributes.find(attr);
                if (it == entry.attributes.end()) {
                    continue;  // pass-through: nothing to evaluate for this attr
                }
                float value = -1.0f;
                try {
                    value = std::stof(it->second);
                } catch (const std::exception& e) {
                    logging::warning("Invalid " + attr + " value '" + it->second +
                        "' for transcript " + transcript_id + ": " + e.what());
                    continue;
                }
                if (first) {
                    // First declared attribute becomes the stored value +
                    // display label for per_sample outputs. Tag the sample
                    // with its expression_type once (idempotent).
                    expression_value = value;
                    first = false;
                    if (sample_id.has_value()) {
                        auto& registry = sample_registry::instance();
                        if (registry.contains(*sample_id)) {
                            auto& info = registry.get(*sample_id);
                            if (!info.has_expression_type()) {
                                if      (attr == "counts") info.expr_type = sample_info::expression_type::COUNTS;
                                else if (attr == "TPM")    info.expr_type = sample_info::expression_type::TPM;
                                else if (attr == "FPKM")   info.expr_type = sample_info::expression_type::FPKM;
                                else if (attr == "RPKM")   info.expr_type = sample_info::expression_type::RPKM;
                                // `cov` stays UNKNOWN to match the legacy header convention
                            }
                        }
                    }
                }
                // Evaluate the filter for this attribute (AND semantics):
                // skip if the CLI threshold is set and the transcript's
                // value is below it.
                float threshold = filters.for_attribute(attr);
                if (threshold >= 0 && value < threshold) {
                    drop_by_filter = true;
                }
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

    if (drop_by_filter) {
        counters.discarded_transcripts++;
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
        gff_source, segment_count, expression_value, transcript_biotype, absorb, fuzzy_tolerance,
        counters
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

    auto trim = [](std::string_view s) -> std::string {
        size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string_view::npos) return "";
        size_t end = s.find_last_not_of(" \t\r\n");
        return std::string(s.substr(start, end - start + 1));
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