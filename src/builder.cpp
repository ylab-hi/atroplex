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
#include <charconv>
#include <filesystem>
#include <string_view>
#include <unordered_set>

// class
#include "utility.hpp"
#include "builder.hpp"
#include "build_bam.hpp"
#include "build_gff.hpp"
#include "build_summary.hpp"

// Natural chromosome sort comparator
static bool chromosome_compare(const std::string& a, const std::string& b) {
    auto get_chr_value = [](const std::string& s) -> std::pair<int, std::string_view> {
        std::string_view sv(s);
        if (sv.size() > 3 && sv.substr(0, 3) == "chr") {
            auto suffix = sv.substr(3);
            int num = 0;
            auto [ptr, ec] = std::from_chars(suffix.data(), suffix.data() + suffix.size(), num);
            if (ec == std::errc{} && ptr == suffix.data() + suffix.size()) {
                return {num, {}};  // Pure numeric (chr1, chr22)
            }
            return {1000, suffix};  // Non-numeric (chrX, chrY, chrM)
        }
        return {2000, sv};  // Non-chr prefixed
    };

    auto [num_a, str_a] = get_chr_value(a);
    auto [num_b, str_b] = get_chr_value(b);

    if (num_a != num_b) return num_a < num_b;
    return str_a < str_b;
}

build_summary builder::build_from_samples(grove_type& grove,
                                  const std::vector<sample_info>& samples,
                                  uint32_t threads,
                                  float min_expression,
                                  bool absorb,
                                  int min_replicates,
                                  size_t fuzzy_tolerance,
                                  chromosome_exon_caches* out_exon_caches,
                                  chromosome_gene_segment_indices* out_gene_indices) {
    if (samples.empty()) {
        logging::warning("No samples provided to build genogrove");
        return {};
    }

    // Per-build counters threaded through build_gff / build_bam / segment_builder
    build_counters counters;

    if (threads > 1) {
        logging::info("Note: Multi-threading for build not yet optimized, using single thread");
    }

    logging::info("Populating grove from " + std::to_string(samples.size()) + " sample(s)");

    // Sort samples: annotations first, then samples (preserving order within each group).
    // Annotations establish the reference exon/segment structure that sample transcripts
    // are absorbed against — processing them first ensures absorption rules have a
    // curated parent to match against.
    std::vector<const sample_info*> ordered;
    ordered.reserve(samples.size());
    for (const auto& s : samples) ordered.push_back(&s);
    std::stable_sort(ordered.begin(), ordered.end(),
        [](const sample_info* a, const sample_info* b) {
            bool a_anno = (a->type == "annotation" || a->is_annotation());
            bool b_anno = (b->type == "annotation" || b->is_annotation());
            return a_anno > b_anno;  // annotations first
        });

    // Chromosome-level caches for deduplication across files
    chromosome_exon_caches exon_caches;
    chromosome_segment_caches segment_caches;
    chromosome_gene_segment_indices gene_indices;
    size_t segment_count = 0;

    // Process each sample sequentially (annotations first)
    size_t total = ordered.size();
    size_t current = 0;
    for (const auto* info_ptr : ordered) {
        ++current;
        const auto& info = *info_ptr;
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
            logging::info("[" + std::to_string(current) + "/" + std::to_string(total) + "] Processing: " + filepath.filename().string() +
                         (info.id.empty() ? "" : " (id: " + info.id + ")"));

            // Register sample_info in the registry to get uint32_t ID
            uint32_t registry_id = sample_registry::instance().register_data(info);

            // Build with persistent caches for cross-file deduplication
            build_gff::build(grove, filepath, registry_id, exon_caches, segment_caches, gene_indices, segment_count, threads, min_expression, absorb, fuzzy_tolerance, counters);

            // Estimate grove memory usage after each file
            size_t mem_segments = 0, mem_exons = 0, mem_gene_idx = 0;
            for (const auto& [sid, seg_cache] : segment_caches) {
                for (const auto& [key, seg_ptr] : seg_cache) {
                    auto& seg = get_segment(seg_ptr->get_data());
                    mem_segments += sizeof(segment_feature)
                        + seg.transcript_ids.data_bytes()  // sorted flat vector
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
                        + exon.transcript_ids.data_bytes()  // sorted flat vector
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
        } else if (ftype == gio::filetype::BAM || ftype == gio::filetype::SAM) {
            logging::info("[" + std::to_string(current) + "/" + std::to_string(total) + "] Processing BAM: " + filepath.filename().string() +
                         (info.id.empty() ? "" : " (id: " + info.id + ")"));

            uint32_t registry_id = sample_registry::instance().register_data(info);

            build_bam::build(grove, filepath, registry_id, exon_caches, segment_caches,
                            gene_indices, segment_count, min_expression, absorb, fuzzy_tolerance, counters);
        } else {
            logging::warning("Unsupported file type for: " + filepath.string());
        }
    }

    logging::info("Grove construction complete: " + std::to_string(segment_count) + " segments");

    // --- Post-build replicate merging ---
    if (min_replicates > 0) {
        counters.replicates_merged = merge_replicates(exon_caches, segment_caches, min_replicates);
    }

    // --- Physically remove absorbed segments from the grove ---
    counters.absorbed_segments = remove_tombstones(grove, segment_caches, gene_indices);

    // --- Collect summary statistics ---
    build_summary stats;
    stats.collect(grove, segment_caches, exon_caches, segment_count, counters);

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

    // Move caches to caller if requested
    if (out_exon_caches) {
        *out_exon_caches = std::move(exon_caches);
    }
    if (out_gene_indices) {
        *out_gene_indices = std::move(gene_indices);
    }

    return stats;
}

size_t builder::merge_replicates(
    chromosome_exon_caches& exon_caches,
    chromosome_segment_caches& segment_caches,
    int min_replicates
) {
    auto& registry = sample_registry::instance();

    // Step 1: Build group -> replicate ID mapping
    std::map<std::string, std::vector<uint32_t>> groups;
    for (size_t i = 0; i < registry.size(); ++i) {
        uint32_t rid = static_cast<uint32_t>(i);
        const auto& info = registry.get(rid);
        if (info.type != "sample" || info.group.empty()) continue;
        groups[info.group].push_back(rid);
    }

    if (groups.empty()) {
        logging::warning("No replicate groups found; skipping merge");
        return 0;
    }

    // Step 2: Register merged sample_info for each group
    std::map<std::string, uint32_t> group_merged_ids;
    for (auto& [group_name, replicate_ids] : groups) {
        const auto& first = registry.get(replicate_ids[0]);
        sample_info merged(first);
        merged.id = group_name;
        merged.group = group_name;
        merged.description = "Merged from " + std::to_string(replicate_ids.size()) + " replicates";
        merged.source_file.clear();

        uint32_t merged_id = registry.register_data(std::move(merged));
        group_merged_ids[group_name] = merged_id;

        if (replicate_ids.size() == 1 && min_replicates > 1) {
            logging::warning("Group '" + group_name + "' has only 1 replicate "
                "but min_replicates=" + std::to_string(min_replicates) +
                "; threshold capped to 1");
        }

        logging::info("Replicate group '" + group_name + "': " +
            std::to_string(replicate_ids.size()) + " replicates");
    }

    // Step 3: Walk caches and merge bits/expression
    // Lambda that works on both exon_feature and segment_feature via std::visit
    auto merge_feature = [&](genomic_feature& feature) {
        std::visit([&](auto& f) {
            for (const auto& [group_name, replicate_ids] : groups) {
                size_t rep_count = 0;
                float expr_sum = 0.0f;
                size_t expr_count = 0;

                for (uint32_t rid : replicate_ids) {
                    if (f.sample_idx.test(rid)) {
                        rep_count++;
                        if (f.expression.has(rid)) {
                            expr_sum += f.expression.get(rid);
                            expr_count++;
                        }
                    }
                }

                // Cap threshold at group size (singletons always pass)
                size_t threshold = std::min(
                    static_cast<size_t>(min_replicates),
                    replicate_ids.size());

                uint32_t merged_id = group_merged_ids.at(group_name);
                if (rep_count >= threshold) {
                    f.add_sample(merged_id);
                    if (expr_count > 0) {
                        f.set_expression(merged_id,
                            expr_sum / static_cast<float>(expr_count));
                    }
                }

                // Clear replicate bits and expression
                for (uint32_t rid : replicate_ids) {
                    f.sample_idx.clear(rid);
                    f.expression.remove(rid);
                }
            }
        }, feature);
    };

    // Walk exon caches
    for (auto& [seqid, exon_cache] : exon_caches) {
        for (auto& [coord, exon_ptr] : exon_cache) {
            merge_feature(exon_ptr->get_data());
        }
    }

    // Walk segment caches
    for (auto& [seqid, seg_cache] : segment_caches) {
        for (auto& [key, seg_ptr] : seg_cache) {
            auto& feature = seg_ptr->get_data();
            if (!is_segment(feature)) continue;
            if (get_segment(feature).absorbed) continue;
            merge_feature(feature);
        }
    }

    // Step 4: Mark replicate entries as type="replicate" so they're excluded from stats
    size_t collapsed = 0;
    for (const auto& [group_name, replicate_ids] : groups) {
        for (uint32_t rid : replicate_ids) {
            registry.get(rid).type = "replicate";
            ++collapsed;
        }
    }

    logging::info("Replicate merging complete: " +
        std::to_string(group_merged_ids.size()) + " groups, min_replicates=" +
        std::to_string(min_replicates));
    return collapsed;
}

size_t builder::remove_tombstones(
    grove_type& /*grove*/,
    chromosome_segment_caches& segment_caches,
    chromosome_gene_segment_indices& gene_indices
) {
    // Physical removal via grove.remove_key() is O(N × E) because
    // genogrove's graph_overlay::remove_edges_to scans the full adjacency
    // map on every call (even though segments are never edge targets, so
    // the scan always finds zero matches). On realistic indices this
    // blocks the build for hours after "Grove construction complete".
    //
    // For now we only COUNT tombstones and prune them from segment_caches
    // and gene_indices so downstream summary collection doesn't see them.
    // The tombstoned keys themselves stay in the B+ tree; every consumer
    // already skips them defensively via `if (seg.absorbed) continue;`.
    // This wastes a bit of memory (and bloats the serialized .ggx) but
    // keeps the build moving. See TODO: replace with a batch-removal
    // genogrove API that skips the incoming-edge scan when segments are
    // known to be edge-source-only.

    size_t removed = 0;

    // Count absorbed entries via the authoritative source (gene_indices).
    for (auto& [seqid, gene_idx] : gene_indices) {
        for (auto& [gene_id, entries] : gene_idx) {
            for (auto& entry : entries) {
                auto& feature = entry.segment->get_data();
                if (!is_segment(feature)) continue;
                if (get_segment(feature).absorbed) ++removed;
            }
        }
    }

    // Prune absorbed entries from gene_indices so that downstream stats
    // collection never sees a tombstone.
    for (auto& [seqid, gene_idx] : gene_indices) {
        for (auto& [gene_id, entries] : gene_idx) {
            std::erase_if(entries, [](const segment_chain_entry& e) {
                return get_segment(e.segment->get_data()).absorbed;
            });
        }
    }

    // Safety net: try_reverse_absorption already erases from segment_cache
    // at the moment it tombstones, but drop any stragglers defensively.
    for (auto& [seqid, seg_cache] : segment_caches) {
        for (auto it = seg_cache.begin(); it != seg_cache.end(); ) {
            auto& feature = it->second->get_data();
            if (is_segment(feature) && get_segment(feature).absorbed) {
                it = seg_cache.erase(it);
            } else {
                ++it;
            }
        }
    }

    if (removed > 0) {
        logging::info("Counted " + std::to_string(removed) +
            " absorbed (tombstoned) segments; physical removal deferred "
            "pending a batch genogrove remove_key variant");
    }

    return removed;
}

build_summary builder::build_from_files(grove_type& grove,
                                const std::vector<std::string>& files,
                                uint32_t threads,
                                chromosome_exon_caches* out_exon_caches) {
    // Parse headers to extract metadata from each file
    std::vector<sample_info> samples;
    samples.reserve(files.size());

    for (const auto& filepath : files) {
        gio::filetype_detector detector;
        auto [ftype, is_gzipped] = detector.detect_filetype(filepath);

        sample_info info;
        if (ftype == gio::filetype::BAM || ftype == gio::filetype::SAM) {
            info = build_bam::parse_header(filepath);
        } else {
            info = build_gff::parse_header(filepath);
        }
        samples.push_back(std::move(info));
    }

    return build_from_samples(grove, samples, threads, -1.0f, true, 0, 5, out_exon_caches);
}
