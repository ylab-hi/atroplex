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
#include <memory>
#include <string_view>
#include <system_error>
#include <unordered_set>

// class
#include "utility.hpp"
#include "builder.hpp"
#include "build_bam.hpp"
#include "build_gff.hpp"
#include "build_summary.hpp"
#include "quant_sidecar.hpp"

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
                                  const expression_filters& filters,
                                  bool absorb,
                                  size_t fuzzy_tolerance,
                                  bool prune_tombstones,
                                  bool include_scaffolds,
                                  const std::string& qtx_path,
                                  chromosome_exon_caches* out_exon_caches,
                                  bool annotated_loci_only,
                                  const std::unordered_set<std::string>& chromosomes_filter) {
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

    // Quantification sidecar (final segment-major .qtx + per-sample temp
    // streams). Empty `qtx_path` disables sidecar writing entirely (tests,
    // build_from_files, callers without an output directory).
    //
    // The final single-file output is written at `qtx_path`; per-sample
    // temp streams go in `{qtx_path}.streams/` and are K-way merged into
    // the final file after the main build + tombstone phases
    // complete. The streams directory is deleted on successful merge;
    // kept for debugging on failure. (`.tmp` suffix is reserved for
    // merge_to_qtx's own atomic-rename scratch file.)
    std::filesystem::path qtx_output_path;
    std::filesystem::path qtx_temp_dir;
    bool sidecar_enabled = !qtx_path.empty();
    if (sidecar_enabled) {
        qtx_output_path = qtx_path;
        // Use `.streams` suffix rather than `.tmp` to avoid collision
        // with merge_to_qtx's own atomic-rename scratch file, which
        // sits at `{qtx_path}.tmp`. The streams directory holds per-
        // sample `.stream` files during build and is deleted after a
        // successful merge.
        qtx_temp_dir    = qtx_output_path.string() + ".streams";
        std::error_code ec;
        std::filesystem::create_directories(qtx_temp_dir, ec);
        if (ec) {
            logging::warning("Failed to create sidecar temp directory '" +
                qtx_temp_dir.string() + "': " + ec.message() +
                " — sidecar writing disabled");
            sidecar_enabled = false;
        } else {
            logging::info("Quantification sidecar will be written to: " +
                qtx_output_path.string());
        }
    }

    // Accumulate per-sample stream paths + SampleMetadata in the same
    // order they are produced so the final merge can pair them up.
    std::vector<std::filesystem::path>         sidecar_stream_paths;
    std::vector<quant_sidecar::SampleMetadata> sidecar_sample_meta;

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

        // Open a per-sample temp-stream writer when a sidecar output was
        // configured. The temp file is named by numeric registry id so
        // the filename is always filesystem-safe and unambiguous. The
        // accompanying SampleMetadata carries the human-readable name
        // through to the final .qtx.
        std::unique_ptr<quant_sidecar::SampleStreamWriter> sidecar;
        std::filesystem::path              current_stream_path;
        std::optional<quant_sidecar::SampleMetadata> current_sample_meta;
        auto open_sidecar = [&](uint32_t registry_id) {
            if (!sidecar_enabled) return;
            std::string sample_label = info.id.empty()
                ? ("sample_" + std::to_string(registry_id))
                : info.id;
            current_stream_path =
                qtx_temp_dir / (std::to_string(registry_id) + ".stream");
            try {
                sidecar = std::make_unique<quant_sidecar::SampleStreamWriter>(
                    current_stream_path, registry_id);
                current_sample_meta = quant_sidecar::SampleMetadata{
                    registry_id,
                    static_cast<uint8_t>(info.expr_type),
                    sample_label
                };
            } catch (const std::exception& e) {
                logging::warning("Failed to open sidecar stream for sample '" +
                    sample_label + "': " + e.what());
                sidecar.reset();
                current_sample_meta.reset();
            }
        };

        // Dispatch to appropriate builder based on file type
        if (ftype == gio::filetype::GFF || ftype == gio::filetype::GTF) {
            logging::info("[" + std::to_string(current) + "/" + std::to_string(total) + "] Processing: " + filepath.filename().string() +
                         (info.id.empty() ? "" : " (id: " + info.id + ")"));

            // Register sample_info in the registry to get uint32_t ID
            uint32_t registry_id = sample_registry::instance().register_data(info);
            open_sidecar(registry_id);

            // Build with persistent caches for cross-file deduplication
            build_gff::build(grove, filepath, registry_id, exon_caches, segment_caches, segment_count, threads, filters, absorb, fuzzy_tolerance, include_scaffolds, counters, sidecar.get(), annotated_loci_only, chromosomes_filter);
        } else if (ftype == gio::filetype::BAM || ftype == gio::filetype::SAM) {
            logging::info("[" + std::to_string(current) + "/" + std::to_string(total) + "] Processing BAM: " + filepath.filename().string() +
                         (info.id.empty() ? "" : " (id: " + info.id + ")"));

            uint32_t registry_id = sample_registry::instance().register_data(info);
            open_sidecar(registry_id);

            build_bam::build(grove, filepath, registry_id, exon_caches, segment_caches,
                            segment_count, filters, absorb, fuzzy_tolerance, include_scaffolds, counters, sidecar.get(),
                            annotated_loci_only, chromosomes_filter);
        } else {
            logging::warning("Unsupported file type for: " + filepath.string());
        }

        // Explicitly finalize so I/O errors surface as warnings here
        // (destructor-driven flushes swallow exceptions). Empty streams
        // still produce a valid temp-stream file with zero records, which
        // the merge handles naturally.
        if (sidecar) {
            bool finalized_ok = false;
            try {
                sidecar->finalize();
                finalized_ok = true;
            } catch (const std::exception& e) {
                logging::warning("Failed to finalize sidecar stream '" +
                    sidecar->path().string() + "': " + e.what());
            }
            if (finalized_ok && current_sample_meta) {
                sidecar_stream_paths.push_back(current_stream_path);
                sidecar_sample_meta.push_back(std::move(*current_sample_meta));
            }
        }
        current_sample_meta.reset();
    }

    // --- Count (and optionally physically remove) absorbed segments ---
    // Collect the tombstoned segment_index values so the sidecar merge
    // below can skip emitting dead records for them in the final .qtx.
    // Tombstones that were absorbed into a parent (via
    // `absorb_into_parent`) populate `tombstone_remap` so merge_to_qtx
    // rewrites those records to the parent's segment_index instead of
    // dropping them. Issue #34.
    std::unordered_set<uint64_t> tombstoned_seg_indices;
    std::unordered_map<uint64_t, uint64_t> tombstone_remap;
    counters.absorbed_segments = remove_tombstones(
        grove, segment_caches, prune_tombstones,
        &tombstoned_seg_indices, &tombstone_remap);

    size_t live_segments = (segment_count >= counters.absorbed_segments)
        ? segment_count - counters.absorbed_segments
        : segment_count;
    std::string tombstone_note;
    if (counters.absorbed_segments > 0) {
        tombstone_note = " (" + std::to_string(counters.absorbed_segments)
            + (prune_tombstones ? " tombstones pruned)" : " tombstones)");
    }
    logging::info("Grove construction complete: " + std::to_string(live_segments)
        + " segments" + tombstone_note);

    // End-of-build memory estimate — walks segment_caches / exon_caches /
    // gene_indices once and sums struct-owned heap allocations. Approximate
    // (misses hash-table bucket overhead, grove tree nodes, and allocator
    // slack) but useful for sizing expectations. Moved here from the old
    // per-file logging because that made build runtime O(N_files ×
    // total_features) for a non-load-bearing log line.
    {
        size_t mem_segments = 0, mem_exons = 0, mem_gene_idx = 0;
        for (const auto& [sid, seg_cache] : segment_caches) {
            for (const auto& [key, seg_ptr] : seg_cache) {
                auto& seg = get_segment(seg_ptr->get_data());
                mem_segments += sizeof(segment_feature)
                    + seg.transcript_ids.data_bytes()
                    + seg.sample_idx.word_count() * 8
                    + seg.transcript_biotypes.size() * 40;
            }
        }
        for (const auto& [sid, exon_cache] : exon_caches) {
            for (const auto& [coord, exon_ptr] : exon_cache) {
                auto& exon = get_exon(exon_ptr->get_data());
                mem_exons += sizeof(exon_feature)
                    + exon.id.capacity()
                    + exon.transcript_ids.data_bytes()
                    + exon.sample_idx.word_count() * 8;
            }
        }
        size_t mem_edges = grove.edge_count() * (sizeof(edge_metadata) + 32);
        size_t mem_total = mem_segments + mem_exons + mem_edges;

        auto fmt_mb = [](size_t bytes) {
            return std::to_string(bytes / (1024 * 1024)) + " MB";
        };
        logging::info("Memory estimate: " + fmt_mb(mem_total) +
            " (segments: " + fmt_mb(mem_segments) +
            ", exons: " + fmt_mb(mem_exons) +
            ", edges: " + fmt_mb(mem_edges) + ")");
    }

    // --- Phase: finalize quantification sidecar (K-way merge) ---
    // Happens after remove_tombstones so that any future logic that
    // rewrites the sidecar in response to that phase can hook in
    // here. The merge itself simply consumes the
    // per-sample temp streams produced during the build loop above.
    if (sidecar_enabled && !sidecar_stream_paths.empty()) {
        try {
            quant_sidecar::merge_to_qtx(
                qtx_output_path,
                sidecar_stream_paths,
                sidecar_sample_meta,
                /*max_fds_per_pass=*/256,
                tombstoned_seg_indices,
                tombstone_remap);
            logging::info("Merged " + std::to_string(sidecar_stream_paths.size()) +
                " per-sample quantification streams into " +
                qtx_output_path.string());

            // Merge succeeded atomically — safe to delete the temp dir.
            std::error_code ec;
            std::filesystem::remove_all(qtx_temp_dir, ec);
            if (ec) {
                logging::warning("Failed to remove sidecar temp dir '" +
                    qtx_temp_dir.string() + "': " + ec.message());
            }
        } catch (const std::exception& e) {
            logging::error("Failed to merge quantification streams into '" +
                qtx_output_path.string() + "': " + e.what() +
                " (temp streams preserved at " + qtx_temp_dir.string() + ")");
            throw;
        }
    }

    // --- Collect summary statistics ---
    build_summary stats;
    stats.collect(grove, segment_caches, exon_caches, segment_count, counters);

    // Concise one-line registry summary (detailed per-chromosome breakdown
    // still lives on stats.per_chromosome and surfaces in the .ggx.summary).
    logging::info("Chromosomes: " + std::to_string(stats.per_chromosome.size()) +
                  ", Genes: " + std::to_string(stats.total_genes) +
                  ", Sources: " + std::to_string(source_registry::instance().size()) +
                  ", Transcripts: " + std::to_string(stats.total_transcripts) +
                  ", Samples: " + std::to_string(sample_registry::instance().size()));

    if (out_exon_caches) {
        *out_exon_caches = std::move(exon_caches);
    }

    return stats;
}

size_t builder::remove_tombstones(
    grove_type& grove,
    chromosome_segment_caches& segment_caches,
    bool physical,
    std::unordered_set<uint64_t>* out_tombstoned_segment_indices,
    std::unordered_map<uint64_t, uint64_t>* out_tombstone_remap
) {
    // Phase 1: walk the grove's B+ tree to find tombstoned segments.
    // This replaces the former gene_indices walk — the grove IS the
    // authoritative source of all segments (issue #40).
    std::unordered_set<size_t> tombstone_ids;
    size_t removed = 0;

    auto roots = grove.get_root_nodes();
    for (auto& [seqid, root] : roots) {
        if (!root) continue;
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
                if (!seg.absorbed) continue;

                ++removed;
                if (out_tombstoned_segment_indices) {
                    out_tombstoned_segment_indices->insert(
                        static_cast<uint64_t>(seg.segment_index));
                }
                if (out_tombstone_remap && seg.absorbed_into_idx.has_value()) {
                    (*out_tombstone_remap)[static_cast<uint64_t>(seg.segment_index)] =
                        static_cast<uint64_t>(*seg.absorbed_into_idx);
                }
                if (physical) {
                    tombstone_ids.insert(seg.segment_index);
                    grove.remove_key(seqid, key);
                }
            }
            node = node->get_next();
        }
    }

    // Phase 2 (physical only): drop orphan EXON_TO_EXON chain edges.
    if (physical && !tombstone_ids.empty()) {
        grove.remove_edges_if([&tombstone_ids](const auto& e) {
            return tombstone_ids.count(e.metadata.id) > 0;
        });
    }

    // Phase 3: safety net on segment_caches.
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

    // Phase 4: collapse absorption chains (issue #34).
    if (out_tombstone_remap) {
        for (auto& [old_idx, new_idx] : *out_tombstone_remap) {
            std::unordered_set<uint64_t> visited{old_idx};
            auto it = out_tombstone_remap->find(new_idx);
            while (it != out_tombstone_remap->end() &&
                   visited.insert(new_idx).second) {
                new_idx = it->second;
                it = out_tombstone_remap->find(new_idx);
            }
        }
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

    return build_from_samples(grove, samples, threads, expression_filters{}, true, 5, false, false, /*qtx_path=*/"", out_exon_caches);
}
