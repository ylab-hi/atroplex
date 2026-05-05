/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/export_gtf.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "gtf_writer.hpp"
#include "utility.hpp"

namespace fs = std::filesystem;

namespace {

// Per-transcript data collected during grove walk
struct transcript_entry {
    std::string transcript_id;
    std::string transcript_biotype;
    std::string seqid;
    std::string source;
    char strand;
    size_t start;
    size_t end;
    int exon_count;
    std::vector<std::pair<size_t, size_t>> exons;  // (start, end) from chain walk
};

// Per-gene data: bounds + transcripts
struct gene_entry {
    std::string gene_id;
    std::string gene_name;
    std::string gene_biotype;
    std::string seqid;
    std::string source;
    char strand;
    size_t start = std::numeric_limits<size_t>::max();
    size_t end = 0;
    std::vector<transcript_entry> transcripts;

    void update_bounds(size_t s, size_t e) {
        start = std::min(start, s);
        end = std::max(end, e);
    }
};

// Gene sort key: (seqid, start)
struct gene_sort_key {
    std::string seqid;
    size_t start;
    bool operator<(const gene_sort_key& o) const {
        if (seqid != o.seqid) return seqid < o.seqid;
        return start < o.start;
    }
};

} // anonymous namespace

namespace subcall {

cxxopts::Options export_gtf::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex export",
        "Export per-sample GTF files from a genogrove index");

    add_common_options(options);

    options.add_options("Export")
        ("sample", "Export only specific sample(s) by ID",
            cxxopts::value<std::vector<std::string>>())
        ("gene", "Filter by gene ID or gene name",
            cxxopts::value<std::vector<std::string>>())
        ("region", "Filter by genomic region (chr:start-end)",
            cxxopts::value<std::string>())
        ("min-samples", "Only export features present in >= N samples",
            cxxopts::value<size_t>()->default_value("0"))
        ("conserved-only", "Only export features satisfying the conservation "
            "threshold. Pair with --conserved-fraction to relax the strict "
            "default. Without --conserved-fraction this means present in every "
            "sample-typed entry (annotations excluded).")
        ("conserved-fraction", "Fraction of sample-typed entries (0,1] a segment "
            "must appear in to satisfy --conserved-only. Default 1.0 = strict. "
            "Mirrors the inspect convention so the two subcommands agree on "
            "what 'conserved' means.",
            cxxopts::value<double>()->default_value("1.0"))
        ("biotype", "Filter by gene biotype (e.g., protein_coding)",
            cxxopts::value<std::vector<std::string>>())
        ("source", "Filter by GFF source (e.g., HAVANA, TALON)",
            cxxopts::value<std::vector<std::string>>())
        ;

    return options;
}

void export_gtf::validate(const cxxopts::ParseResult& args) {
    if (!args.count("genogrove") && !args.count("build-from") && !args.count("manifest")) {
        throw std::runtime_error(
            "No input specified. Use -g/--genogrove, -m/--manifest, or -b/--build-from");
    }

    if (args.count("region")) {
        [[maybe_unused]] auto _ = parse_region(args["region"].as<std::string>());
    }

    if (args.count("conserved-fraction")) {
        double f = args["conserved-fraction"].as<double>();
        if (!(f > 0.0 && f <= 1.0)) {
            throw std::runtime_error(
                "--conserved-fraction must be in (0, 1]; got " + std::to_string(f));
        }
    }
}

void export_gtf::parse_filters(const cxxopts::ParseResult& args) {
    if (args.count("sample")) {
        auto ids = args["sample"].as<std::vector<std::string>>();
        filters_.sample_ids.insert(ids.begin(), ids.end());
    }
    if (args.count("gene")) {
        auto ids = args["gene"].as<std::vector<std::string>>();
        filters_.gene_ids.insert(ids.begin(), ids.end());
    }
    if (args.count("region")) {
        filters_.region = parse_region(args["region"].as<std::string>());
    }
    filters_.min_samples = args["min-samples"].as<size_t>();
    filters_.conserved_only = args.count("conserved-only") > 0;
    filters_.conserved_fraction = args["conserved-fraction"].as<double>();
    if (args.count("biotype")) {
        auto bt = args["biotype"].as<std::vector<std::string>>();
        filters_.biotypes.insert(bt.begin(), bt.end());
    }
    if (args.count("source")) {
        auto src = args["source"].as<std::vector<std::string>>();
        filters_.sources.insert(src.begin(), src.end());
    }
}

export_gtf::region_spec export_gtf::parse_region(const std::string& region_str) {
    // Format: chr:start-end
    auto colon = region_str.find(':');
    if (colon == std::string::npos) {
        throw std::runtime_error("Invalid region format (expected chr:start-end): " + region_str);
    }
    auto dash = region_str.find('-', colon + 1);
    if (dash == std::string::npos) {
        throw std::runtime_error("Invalid region format (expected chr:start-end): " + region_str);
    }

    region_spec spec;
    spec.seqid = region_str.substr(0, colon);
    spec.start = std::stoull(region_str.substr(colon + 1, dash - colon - 1));
    spec.end = std::stoull(region_str.substr(dash + 1));

    if (spec.start >= spec.end) {
        throw std::runtime_error("Invalid region: start must be < end: " + region_str);
    }
    return spec;
}

void export_gtf::execute(const cxxopts::ParseResult& args) {
    parse_filters(args);

    // Determine output directory
    std::string fallback;
    if (args.count("genogrove")) {
        fallback = args["genogrove"].as<std::string>();
    } else if (args.count("manifest")) {
        fallback = args["manifest"].as<std::string>();
    } else if (args.count("build-from")) {
        fallback = args["build-from"].as<std::vector<std::string>>().front();
    }
    auto out_dir = resolve_output_dir(args, fallback);

    // Determine which samples to export
    size_t total_samples = 0;
    // Conservation denominator — only `type == "sample"` entries count, so
    // the threshold matches what `inspect` enforces. `total_samples` (above)
    // also includes annotations and is kept separate for the export-iteration
    // log message.
    size_t sample_typed_count = 0;
    std::vector<uint32_t> export_sample_ids;
    auto& registry = sample_registry::instance();

    for (uint32_t i = 0; i < registry.size(); ++i) {
        const auto& info = registry.get(static_cast<uint32_t>(i));
        if (info.type == "replicate") continue;
        total_samples++;
        if (info.type == "sample") sample_typed_count++;

        if (!filters_.sample_ids.empty() && !filters_.sample_ids.count(info.id))
            continue;

        export_sample_ids.push_back(i);
    }

    if (export_sample_ids.empty()) {
        logging::warning("No samples match the --sample filter");
        return;
    }

    logging::info("Exporting " + std::to_string(export_sample_ids.size()) +
                  " sample(s) from index (" + std::to_string(total_samples) + " total)");

    // Translate --conserved-fraction into a concrete sample-count threshold.
    // ceil so fractions just below 1.0 still demand nearly every sample;
    // floor at 1 so an empty sample-typed registry never silently classifies
    // every segment as conserved (the default sample_idx.count() == 0 would
    // otherwise satisfy is_conserved(0)).
    size_t conserved_min_required = 0;
    if (filters_.conserved_only) {
        conserved_min_required = static_cast<size_t>(std::ceil(
            static_cast<double>(sample_typed_count) * filters_.conserved_fraction));
        if (conserved_min_required < 1) conserved_min_required = 1;
        if (filters_.conserved_fraction < 1.0) {
            logging::info("Conservation threshold: "
                + std::to_string(conserved_min_required) + " / "
                + std::to_string(sample_typed_count) + " samples (fraction="
                + std::to_string(filters_.conserved_fraction) + ")");
        }
    }

    // Collect per-sample data in one grove walk
    // sample_id -> gene_id -> gene_entry
    std::map<uint32_t, std::unordered_map<std::string, gene_entry>> per_sample_data;

    auto roots = grove->get_root_nodes();
    for (auto& [seqid, root] : roots) {
        if (!root) continue;

        // Region filter: skip entire chromosome if not matching
        if (filters_.region && seqid != filters_.region->seqid) continue;

        // Walk to first leaf
        auto* node = root;
        while (!node->get_is_leaf()) {
            auto& children = node->get_children();
            if (children.empty()) break;
            node = children[0];
        }

        while (node) {
            for (auto* key : node->get_keys()) {
                auto& feature = key->get_data();
                if (!is_segment(feature)) continue;

                auto& seg = get_segment(feature);
                if (seg.absorbed) continue;

                auto& coord = key->get_value();

                // Apply structural filters
                // Region overlap
                if (filters_.region) {
                    if (coord.get_end() <= filters_.region->start ||
                        coord.get_start() >= filters_.region->end)
                        continue;
                }

                // Gene filter
                if (!filters_.gene_ids.empty()) {
                    if (!filters_.gene_ids.count(seg.gene_id()) &&
                        !filters_.gene_ids.count(seg.gene_name()))
                        continue;
                }

                // Min-samples / conserved
                if (filters_.min_samples > 0 && seg.sample_count() < filters_.min_samples)
                    continue;
                if (filters_.conserved_only && !seg.is_conserved(conserved_min_required))
                    continue;

                // Biotype filter
                if (!filters_.biotypes.empty() && !filters_.biotypes.count(seg.gene_biotype()))
                    continue;

                // Source filter
                if (!filters_.sources.empty()) {
                    bool has_match = false;
                    source_registry::instance().for_each(seg.sources,
                        [&](const std::string& src) {
                            if (filters_.sources.count(src)) has_match = true;
                        });
                    if (!has_match) continue;
                }

                // Traverse exon chain for this segment
                std::vector<std::pair<size_t, size_t>> exon_coords;
                size_t edge_id = seg.segment_index;

                auto first_exons = grove->get_neighbors_if(key,
                    [edge_id](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON
                            && e.id == edge_id;
                    });

                if (!first_exons.empty()) {
                    auto* current = first_exons.front();
                    while (current) {
                        auto& ec = current->get_value();
                        exon_coords.emplace_back(ec.get_start(), ec.get_end());
                        auto next = grove->get_neighbors_if(current,
                            [edge_id](const edge_metadata& e) {
                                return e.type == edge_metadata::edge_type::EXON_TO_EXON
                                    && e.id == edge_id;
                            });
                        current = next.empty() ? nullptr : next.front();
                    }
                }

                // Sort exons by start position
                std::sort(exon_coords.begin(), exon_coords.end());

                // Resolve source string (first source in bitfield)
                std::string source_str = "atroplex";
                source_registry::instance().for_each(seg.sources,
                    [&source_str](const std::string& src) {
                        source_str = src;  // takes the last one; fine for single-source
                    });

                // Add to per-sample collections
                for (uint32_t sid : export_sample_ids) {
                    if (!seg.sample_idx.test(sid)) continue;

                    auto& gene = per_sample_data[sid][seg.gene_id()];
                    if (gene.gene_id.empty()) {
                        gene.gene_id = seg.gene_id();
                        gene.gene_name = seg.gene_name();
                        gene.gene_biotype = seg.gene_biotype();
                        gene.seqid = seqid;
                        gene.source = source_str;
                        gene.strand = coord.get_strand();
                    }
                    gene.update_bounds(coord.get_start(), coord.get_end());

                    // One transcript_entry per transcript_id on this segment
                    for (uint32_t tx_id : seg.transcript_ids) {
                        std::string tx_name = transcript_registry::instance().resolve(tx_id);

                        // Look up biotype for this transcript
                        std::string tx_biotype;
                        auto bt_it = seg.transcript_biotypes.find(tx_id);
                        if (bt_it != seg.transcript_biotypes.end()) {
                            tx_biotype = bt_it->second;
                        }

                        transcript_entry entry;
                        entry.transcript_id = tx_name;
                        entry.transcript_biotype = tx_biotype;
                        entry.seqid = seqid;
                        entry.source = source_str;
                        entry.strand = coord.get_strand();
                        entry.start = coord.get_start();
                        entry.end = coord.get_end();
                        entry.exon_count = static_cast<int>(exon_coords.size());
                        entry.exons = exon_coords;

                        gene.transcripts.push_back(std::move(entry));
                    }
                }
            }
            node = node->get_next();
        }
    }

    // Write per-sample GTFs
    size_t files_written = 0;
    size_t total_transcripts_written = 0;

    for (uint32_t sid : export_sample_ids) {
        auto it = per_sample_data.find(sid);
        if (it == per_sample_data.end() || it->second.empty()) {
            logging::info("  " + registry.get(sid).id + ": no features pass filters, skipping");
            continue;
        }

        auto& gene_map = it->second;
        const auto& sinfo = registry.get(sid);

        // Sort genes by (seqid, start)
        std::map<gene_sort_key, const gene_entry*> sorted_genes;
        for (auto& [gid, gene] : gene_map) {
            sorted_genes[{gene.seqid, gene.start}] = &gene;
        }

        fs::path output_path = out_dir / (sinfo.id + ".gtf");
        std::ofstream ofs(output_path.string());
        if (!ofs.is_open()) {
            throw std::runtime_error("Cannot create output file: " + output_path.string());
        }

        gtf_writer::write_header(ofs, sinfo);

        size_t sample_tx_count = 0;

        for (auto& [sort_key, gene_ptr] : sorted_genes) {
            auto& gene = *gene_ptr;

            // Write gene line
            gtf_writer::gtf_gene g;
            g.seqid = gene.seqid;
            g.source = gene.source;
            g.strand = gene.strand;
            g.start = gene.start;
            g.end = gene.end;
            g.gene_id = gene.gene_id;
            g.gene_name = gene.gene_name;
            g.gene_biotype = gene.gene_biotype;
            gtf_writer::write_gene(ofs, g);

            // Sort transcripts by start
            auto sorted_tx = gene.transcripts;
            std::sort(sorted_tx.begin(), sorted_tx.end(),
                [](const transcript_entry& a, const transcript_entry& b) {
                    return a.start < b.start;
                });

            for (auto& tx : sorted_tx) {
                // Write transcript line
                gtf_writer::gtf_transcript t;
                t.seqid = tx.seqid;
                t.source = tx.source;
                t.strand = tx.strand;
                t.start = tx.start;
                t.end = tx.end;
                t.gene_id = gene.gene_id;
                t.gene_name = gene.gene_name;
                t.transcript_id = tx.transcript_id;
                t.transcript_biotype = tx.transcript_biotype;
                t.exon_count = tx.exon_count;
                gtf_writer::write_transcript(ofs, t);

                // Write exon lines
                auto sorted_exons = tx.exons;
                if (tx.strand == '-') {
                    std::sort(sorted_exons.begin(), sorted_exons.end(),
                        std::greater<>());
                }

                int exon_num = 1;
                for (auto& [estart, eend] : sorted_exons) {
                    gtf_writer::gtf_exon e;
                    e.seqid = tx.seqid;
                    e.source = tx.source;
                    e.strand = tx.strand;
                    e.start = estart;
                    e.end = eend;
                    e.gene_id = gene.gene_id;
                    e.gene_name = gene.gene_name;
                    e.transcript_id = tx.transcript_id;
                    e.exon_number = exon_num++;
                    gtf_writer::write_exon(ofs, e);
                }

                sample_tx_count++;
            }
        }

        files_written++;
        total_transcripts_written += sample_tx_count;
        logging::info("  " + sinfo.id + ": " + std::to_string(gene_map.size()) +
                      " genes, " + std::to_string(sample_tx_count) + " transcripts -> " +
                      output_path.string());
    }

    logging::info("Export complete: " + std::to_string(files_written) + " file(s), " +
                  std::to_string(total_transcripts_written) + " total transcripts");
}

} // namespace subcall