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
                                  int min_replicates,
                                  size_t fuzzy_tolerance,
                                  bool prune_tombstones,
                                  bool include_scaffolds,
                                  const std::string& qtx_path,
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

    // Quantification sidecar (final segment-major .qtx + per-sample temp
    // streams). Empty `qtx_path` disables sidecar writing entirely (tests,
    // build_from_files, callers without an output directory).
    //
    // The final single-file output is written at `qtx_path`; per-sample
    // temp streams go in `{qtx_path}.streams/` and are K-way merged into
    // the final file after the main build + tombstone + replicate phases
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
            build_gff::build(grove, filepath, registry_id, exon_caches, segment_caches, gene_indices, segment_count, threads, filters, absorb, fuzzy_tolerance, include_scaffolds, counters, sidecar.get());
        } else if (ftype == gio::filetype::BAM || ftype == gio::filetype::SAM) {
            logging::info("[" + std::to_string(current) + "/" + std::to_string(total) + "] Processing BAM: " + filepath.filename().string() +
                         (info.id.empty() ? "" : " (id: " + info.id + ")"));

            uint32_t registry_id = sample_registry::instance().register_data(info);
            open_sidecar(registry_id);

            build_bam::build(grove, filepath, registry_id, exon_caches, segment_caches,
                            gene_indices, segment_count, filters, absorb, fuzzy_tolerance, include_scaffolds, counters, sidecar.get());
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

    // --- Post-build replicate merging ---
    if (min_replicates > 0) {
        counters.replicates_merged = merge_replicates(exon_caches, segment_caches, min_replicates);
    }

    // --- Count (and optionally physically remove) absorbed segments ---
    // Collect the tombstoned segment_index values so the sidecar merge
    // below can skip emitting dead records for them in the final .qtx.
    std::unordered_set<uint64_t> tombstoned_seg_indices;
    counters.absorbed_segments = remove_tombstones(
        grove, segment_caches, gene_indices, prune_tombstones,
        &tombstoned_seg_indices);

    // Live segments = total minus tombstones that were reverse-absorbed.
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
    }

    // --- Phase: finalize quantification sidecar (K-way merge) ---
    // Happens after remove_tombstones + merge_replicates so that any
    // future logic that rewrites the sidecar in response to those
    // phases can hook in here. The merge itself simply consumes the
    // per-sample temp streams produced during the build loop above.
    if (sidecar_enabled && !sidecar_stream_paths.empty()) {
        try {
            quant_sidecar::merge_to_qtx(
                qtx_output_path,
                sidecar_stream_paths,
                sidecar_sample_meta,
                /*max_fds_per_pass=*/256,
                tombstoned_seg_indices);
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
    // NOTE on interaction with the quantification sidecar (.qtx):
    // Replicate merging ORs the per-replicate sample_idx bits into the
    // group's merged sample_id on each feature. The sidecar, however, was
    // written out during the per-sample build phase and only contains
    // records keyed by the ORIGINAL per-replicate sample_ids — the merged
    // group_id is a "phantom" in the sidecar. At query time, code that
    // wants quantitative aggregation across a replicate group has to
    // either (a) expand the group's sample_id back to the constituent
    // replicate IDs via sample_registry::group lookup and union their
    // sidecar records, or (b) report per-replicate without collapsing.
    // The current 3a PR defers this decision to the follow-up PR that
    // wires quant_sidecar::Reader into query/analyze.
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

                for (uint32_t rid : replicate_ids) {
                    if (f.sample_idx.test(rid)) {
                        rep_count++;
                    }
                }

                // Cap threshold at group size (singletons always pass)
                size_t threshold = std::min(
                    static_cast<size_t>(min_replicates),
                    replicate_ids.size());

                uint32_t merged_id = group_merged_ids.at(group_name);
                if (rep_count >= threshold) {
                    f.add_sample(merged_id);
                }

                // Clear replicate bits
                for (uint32_t rid : replicate_ids) {
                    f.sample_idx.clear(rid);
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
    grove_type& grove,
    chromosome_segment_caches& segment_caches,
    chromosome_gene_segment_indices& gene_indices,
    bool physical,
    std::unordered_set<uint64_t>* out_tombstoned_segment_indices
) {
    // Phase 1 (always): walk gene_indices, count tombstones, and (when
    // physical == true) issue grove.remove_key calls while we have the
    // seqid + key pointer handy. We also collect the tombstone
    // segment_indices so the orphan-edge sweep below can drop EXON_TO_EXON
    // chain edges that carry a dead id.
    //
    // Note on the physical path cost: each grove.remove_key() call fans
    // out into genogrove's graph_overlay::remove_edges_to, which is O(E)
    // (full adjacency scan). Segments are never edge targets, so the
    // scan always finds zero matches but still pays the full cost — that
    // makes the total sweep O(N × E) and can block builds for minutes
    // to hours on realistic indices. Gated behind --prune-tombstones so
    // callers opt in when they want a clean distributable .ggx.
    std::unordered_set<size_t> tombstone_ids;
    size_t removed = 0;

    for (auto& [seqid, gene_idx] : gene_indices) {
        for (auto& [gene_id, entries] : gene_idx) {
            for (auto& entry : entries) {
                auto& feature = entry.segment->get_data();
                if (!is_segment(feature)) continue;
                auto& seg = get_segment(feature);
                if (!seg.absorbed) continue;

                ++removed;
                if (out_tombstoned_segment_indices) {
                    out_tombstoned_segment_indices->insert(
                        static_cast<uint64_t>(seg.segment_index));
                }
                if (physical) {
                    tombstone_ids.insert(seg.segment_index);
                    grove.remove_key(seqid, entry.segment);
                }
            }
        }
    }

    // Phase 2 (physical only): one pass to drop orphan EXON_TO_EXON
    // edges carrying a tombstone's segment_index. remove_key only
    // handles edges where the segment key itself is source or target;
    // chain edges between exon keys are not touched automatically.
    if (physical && !tombstone_ids.empty()) {
        grove.remove_edges_if([&tombstone_ids](const auto& e) {
            return tombstone_ids.count(e.metadata.id) > 0;
        });
    }

    // Phase 3 (always): prune absorbed entries from gene_indices so that
    // downstream stats collection never sees a tombstone.
    for (auto& [seqid, gene_idx] : gene_indices) {
        for (auto& [gene_id, entries] : gene_idx) {
            std::erase_if(entries, [](const segment_chain_entry& e) {
                return get_segment(e.segment->get_data()).absorbed;
            });
        }
    }

    // Phase 4 (always): safety net on segment_caches. try_reverse_absorption
    // already erases from segment_cache at the moment it tombstones, but
    // drop any stragglers defensively.
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

    return build_from_samples(grove, samples, threads, expression_filters{}, true, 0, 5, false, false, /*qtx_path=*/"", out_exon_caches);
}
