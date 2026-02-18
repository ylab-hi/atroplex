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
#include <algorithm>
#include <filesystem>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

// class
#include "utility.hpp"
#include "builder.hpp"
#include "build_gff.hpp"

// Natural chromosome sort comparator
static bool chromosome_compare(const std::string& a, const std::string& b) {
    auto get_chr_value = [](const std::string& s) -> std::pair<int, std::string> {
        if (s.size() > 3 && s.substr(0, 3) == "chr") {
            std::string suffix = s.substr(3);
            try {
                size_t pos;
                int num = std::stoi(suffix, &pos);
                if (pos == suffix.size()) {
                    return {num, ""};  // Pure numeric (chr1, chr22)
                }
            } catch (...) {}
            return {1000, suffix};  // Non-numeric (chrX, chrY, chrM)
        }
        return {2000, s};  // Non-chr prefixed
    };

    auto [num_a, str_a] = get_chr_value(a);
    auto [num_b, str_b] = get_chr_value(b);

    if (num_a != num_b) return num_a < num_b;
    return str_a < str_b;
}

namespace {
double compute_median(std::vector<size_t>& values) {
    if (values.empty()) return 0;
    size_t n = values.size();
    auto mid = values.begin() + static_cast<long>(n / 2);
    std::nth_element(values.begin(), mid, values.end());
    if (n % 2 == 0) {
        auto max_lower = *std::max_element(values.begin(), mid);
        return (static_cast<double>(max_lower) + static_cast<double>(*mid)) / 2.0;
    }
    return static_cast<double>(*mid);
}
} // anonymous namespace

index_stats builder::build_from_samples(grove_type& grove,
                                  const std::vector<sample_info>& samples,
                                  uint32_t threads,
                                  float min_expression) {
    if (samples.empty()) {
        logging::warning("No samples provided to build genogrove");
        return {};
    }

    if (threads > 1) {
        logging::info("Note: Multi-threading for build not yet optimized, using single thread");
    }

    logging::info("Populating grove from " + std::to_string(samples.size()) + " sample(s)");

    // Chromosome-level caches for deduplication across files
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    chromosome_gene_segment_indices gene_indices;
    size_t segment_count = 0;

    // Process each sample sequentially
    for (const auto& info : samples) {
        const auto& filepath = info.source_file;

        if (!std::filesystem::exists(filepath)) {
            logging::error("File not found: " + filepath.string());
            continue;
        }

        // Detect file type
        gio::filetype_detector detector;
        auto [ftype, is_gzipped] = detector.detect_filetype(filepath);

        // Dispatch to appropriate builder based on file type
        if (ftype == gio::filetype::GFF || ftype == gio::filetype::GTF) {
            logging::info("Processing GFF/GTF file: " + filepath.string() +
                         (info.id.empty() ? "" : " (id: " + info.id + ")"));

            // Register sample_info in the registry to get uint32_t ID
            uint32_t registry_id = sample_registry::instance().register_data(info);

            // Build with persistent caches for cross-file deduplication
            build_gff::build(grove, filepath, registry_id, exon_caches, segment_caches, gene_indices, segment_count, threads, min_expression);

            // Estimate grove memory usage after each file
            size_t mem_segments = 0, mem_exons = 0, mem_gene_idx = 0;
            for (const auto& [sid, seg_cache] : segment_caches) {
                for (const auto& [key, seg_ptr] : seg_cache) {
                    auto& seg = get_segment(seg_ptr->get_data());
                    mem_segments += sizeof(segment_feature)
                        + seg.transcript_ids.size() * 12  // ~12 bytes per set entry (uint32_t + hash node)
                        + seg.sample_idx.word_count() * 8  // dynamic bitset words
                        + seg.expression.data_bytes()      // lazy flat array
                        + seg.transcript_biotypes.size() * 40;  // uint32_t key + string value
                }
            }
            for (const auto& [sid, exon_cache] : exon_caches) {
                for (const auto& [coord, exon_ptr] : exon_cache) {
                    auto& exon = get_exon(exon_ptr->get_data());
                    mem_exons += sizeof(exon_feature)
                        + exon.id.capacity()
                        + exon.transcript_ids.size() * 12  // uint32_t + hash node
                        + exon.sample_idx.word_count() * 8  // dynamic bitset words
                        + exon.expression.data_bytes();     // lazy flat array
                }
            }
            for (const auto& [sid, gidx] : gene_indices) {
                for (const auto& [gene_id, entries] : gidx) {
                    for (const auto& e : entries) {
                        mem_gene_idx += sizeof(segment_chain_entry)
                            + e.exon_chain.capacity() * sizeof(key_ptr)
                            + e.structure_key.capacity();
                    }
                }
            }
            size_t mem_edges = grove.edge_count() * (sizeof(edge_metadata) + 32);
            size_t mem_total = mem_segments + mem_exons + mem_edges + mem_gene_idx;

            auto fmt_mb = [](size_t bytes) {
                return std::to_string(bytes / (1024 * 1024)) + " MB";
            };
            logging::info("Memory estimate: " + fmt_mb(mem_total) +
                " (segments: " + fmt_mb(mem_segments) +
                ", exons: " + fmt_mb(mem_exons) +
                ", edges: " + fmt_mb(mem_edges) +
                ", gene_idx: " + fmt_mb(mem_gene_idx) + ")");
        } else {
            logging::warning("Unsupported file type for: " + filepath.string());
        }
    }

    logging::info("Grove construction complete: " + std::to_string(segment_count) + " segments");

    // --- Compute basic index statistics from caches ---
    index_stats stats;

    // Annotation sources and sample count from registry
    auto& registry = sample_registry::instance();
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto* info = registry.get(static_cast<uint32_t>(i));
        if (info && !info->id.empty()) {
            stats.annotation_sources.push_back(info->id);
        }
        if (info && info->type == "sample") {
            stats.total_samples++;
        }
    }

    stats.total_chromosomes = segment_caches.size();
    stats.total_segments = segment_count;

    // Gene info from segment features
    struct gene_stats_info {
        std::unordered_set<uint32_t> transcript_ids;
        std::string biotype;
    };
    std::unordered_map<std::string, gene_stats_info> genes;
    std::unordered_map<std::string, std::unordered_set<std::string>> chr_gene_ids;
    std::vector<size_t> exons_per_segment;

    for (const auto& [seqid, seg_cache] : segment_caches) {
        stats.per_chromosome[seqid].segments = seg_cache.size();

        for (const auto& [key, seg_ptr] : seg_cache) {
            auto& feature = seg_ptr->get_data();
            if (!is_segment(feature)) continue;

            auto& seg = get_segment(feature);
            if (seg.absorbed) continue;  // Skip tombstoned segments

            auto& gi = genes[seg.gene_id()];
            gi.biotype = seg.gene_biotype();
            for (const auto& tx : seg.transcript_ids) {
                gi.transcript_ids.insert(tx);
            }
            chr_gene_ids[seqid].insert(seg.gene_id());

            // Collect transcript biotypes
            for (const auto& [tx_id, biotype] : seg.transcript_biotypes) {
                if (!biotype.empty()) {
                    stats.transcripts_by_biotype[biotype]++;
                }
            }

            exons_per_segment.push_back(static_cast<size_t>(seg.exon_count));
            if (seg.exon_count == 1) {
                stats.single_exon_segments++;
            }
        }
    }

    // Count absorbed segments from gene indices
    for (const auto& [seqid, gene_idx] : gene_indices) {
        for (const auto& [gene_id, entries] : gene_idx) {
            for (const auto& entry : entries) {
                if (get_segment(entry.segment->get_data()).absorbed) {
                    stats.absorbed_segments++;
                }
            }
        }
    }

    // Exon counts from exon caches
    for (const auto& [seqid, exon_cache] : exon_caches) {
        stats.per_chromosome[seqid].exons = exon_cache.size();
        stats.total_exons += exon_cache.size();
    }

    // Gene-level statistics
    stats.total_genes = genes.size();
    std::vector<size_t> transcripts_per_gene;
    transcripts_per_gene.reserve(genes.size());
    size_t total_tx = 0;

    for (const auto& [gene_id, gi] : genes) {
        size_t tx_count = gi.transcript_ids.size();
        total_tx += tx_count;
        transcripts_per_gene.push_back(tx_count);

        if (!gi.biotype.empty()) {
            stats.genes_by_biotype[gi.biotype]++;
        }
        if (tx_count > stats.max_transcripts_per_gene) {
            stats.max_transcripts_per_gene = tx_count;
            stats.max_transcripts_gene_id = gene_id;
        }
        if (tx_count == 1) {
            stats.single_isoform_genes++;
        } else {
            stats.multi_isoform_genes++;
        }
    }

    stats.total_transcripts = total_tx;
    if (!transcripts_per_gene.empty()) {
        stats.mean_transcripts_per_gene =
            static_cast<double>(total_tx) / static_cast<double>(transcripts_per_gene.size());
        stats.median_transcripts_per_gene = compute_median(transcripts_per_gene);
    }

    // Exons per segment distribution
    if (!exons_per_segment.empty()) {
        size_t total_exons_sum = std::accumulate(
            exons_per_segment.begin(), exons_per_segment.end(), size_t{0});
        stats.mean_exons_per_segment =
            static_cast<double>(total_exons_sum) / static_cast<double>(exons_per_segment.size());
        stats.median_exons_per_segment = compute_median(exons_per_segment);
        stats.max_exons_per_segment = *std::max_element(
            exons_per_segment.begin(), exons_per_segment.end());
    }

    // Deduplication ratio
    if (stats.total_transcripts > 0) {
        stats.deduplication_ratio =
            static_cast<double>(stats.total_segments) / static_cast<double>(stats.total_transcripts);
    }

    // Per-chromosome gene counts
    for (const auto& [seqid, gene_set] : chr_gene_ids) {
        stats.per_chromosome[seqid].genes = gene_set.size();
    }

    // Total edges from grove
    stats.total_edges = grove.edge_count();

    // Log per-chromosome summary
    std::vector<std::string> seqids;
    for (const auto& [seqid, _] : segment_caches) {
        seqids.push_back(seqid);
    }
    std::sort(seqids.begin(), seqids.end(), chromosome_compare);

    logging::info("Per-chromosome summary:");
    for (const auto& seqid : seqids) {
        auto& cs = stats.per_chromosome[seqid];
        logging::info("  " + seqid + ": " +
            std::to_string(cs.exons) + " exons, " +
            std::to_string(cs.segments) + " segments");
    }

    return stats;
}

index_stats builder::build_from_files(grove_type& grove,
                                const std::vector<std::string>& files,
                                uint32_t threads) {
    // Parse headers to extract metadata from each file
    std::vector<sample_info> samples;
    samples.reserve(files.size());

    for (const auto& filepath : files) {
        // Parse metadata from GFF/GTF header (#property: value format)
        sample_info info = build_gff::parse_header(filepath);
        samples.push_back(std::move(info));
    }

    return build_from_samples(grove, samples, threads);
}
