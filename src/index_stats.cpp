/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "index_stats.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "sample_info.hpp"
#include "utility.hpp"

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

double jaccard_distance(const std::set<const void*>& a, const std::set<const void*>& b) {
    if (a.empty() && b.empty()) return 0.0;
    size_t intersection = 0;
    auto it_a = a.begin(), it_b = b.begin();
    while (it_a != a.end() && it_b != b.end()) {
        if (*it_a == *it_b) {
            intersection++;
            ++it_a;
            ++it_b;
        } else if (*it_a < *it_b) {
            ++it_a;
        } else {
            ++it_b;
        }
    }
    size_t union_size = a.size() + b.size() - intersection;
    if (union_size == 0) return 0.0;
    return 1.0 - static_cast<double>(intersection) / static_cast<double>(union_size);
}

} // anonymous namespace

index_stats index_stats::collect(grove_type& grove) {
    index_stats stats;

    // Collect annotation sources from sample_registry
    auto& registry = sample_registry::instance();
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto* info = registry.get(static_cast<uint32_t>(i));
        if (info && !info->id.empty()) {
            stats.annotation_sources.push_back(info->id);
        }
    }

    // --- Phase 1: Traverse B+ tree to collect all segments ---
    // gene_id -> { set of transcript_ids, biotype }
    struct gene_info {
        std::unordered_set<std::string> transcript_ids;
        std::string biotype;
    };
    std::unordered_map<std::string, gene_info> genes;
    std::vector<size_t> exons_per_segment;

    // Track per-chromosome gene sets
    std::unordered_map<std::string, std::unordered_set<std::string>> chr_gene_ids;

    // All segment key pointers (needed for exon traversal)
    std::vector<std::pair<std::string, key_ptr>> segment_keys; // (seqid, key)

    auto roots = grove.get_root_nodes();
    stats.total_chromosomes = roots.size();

    for (auto& [seqid, root] : roots) {
        if (!root) continue;

        // Walk to first leaf
        auto* node = root;
        while (!node->get_is_leaf()) {
            auto& children = node->get_children();
            if (children.empty()) break;
            node = children[0];
        }

        // Iterate all leaves
        while (node) {
            for (auto* key : node->get_keys()) {
                auto& feature = key->get_data();
                if (!is_segment(feature)) continue;

                auto& seg = get_segment(feature);
                stats.total_segments++;
                stats.per_chromosome[seqid].segments++;

                // Collect gene info
                auto& gi = genes[seg.gene_id];
                gi.biotype = seg.gene_biotype;
                for (const auto& tx : seg.transcript_ids) {
                    gi.transcript_ids.insert(tx);
                }
                chr_gene_ids[seqid].insert(seg.gene_id);

                exons_per_segment.push_back(static_cast<size_t>(seg.exon_count));
                if (seg.exon_count == 1) {
                    stats.single_exon_segments++;
                }

                segment_keys.emplace_back(seqid, key);
            }
            node = node->get_next();
        }
    }

    // --- Phase 2: Traverse graph from segments to collect exons ---
    std::set<const void*> visited_exons; // deduplicate by pointer

    for (auto& [seqid, seg_key] : segment_keys) {
        // Follow SEGMENT_TO_EXON edges
        auto first_exons = grove.get_neighbors_if(seg_key,
            [](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON;
            });

        for (auto* exon_key : first_exons) {
            // Walk exon chain via BFS
            std::vector<key_ptr> to_visit;
            to_visit.push_back(exon_key);

            while (!to_visit.empty()) {
                auto* current = to_visit.back();
                to_visit.pop_back();

                if (!visited_exons.insert(current).second) {
                    continue; // already visited
                }

                stats.total_exons++;
                stats.per_chromosome[seqid].exons++;

                // Follow EXON_TO_EXON edges
                auto next = grove.get_neighbors_if(current,
                    [](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::EXON_TO_EXON;
                    });

                // Deduplicate targets for branching detection
                std::set<const void*> unique_targets;
                for (auto* n : next) {
                    unique_targets.insert(n);
                    if (!visited_exons.count(n)) {
                        to_visit.push_back(n);
                    }
                }
                if (unique_targets.size() > 1) {
                    stats.branching_exons++;

                    // Track top branching exons
                    auto& feature = current->get_data();
                    if (is_exon(feature)) {
                        auto& exon = get_exon(feature);
                        index_stats::branching_exon_info info;
                        info.gene_name = exon.gene_name;
                        info.gene_id = exon.gene_id;
                        info.coordinate = exon.coordinate;
                        info.branches = unique_targets.size();
                        info.transcripts = exon.transcript_ids.size();

                        if (stats.top_branching_exons.size() < index_stats::MAX_TOP_BRANCHING) {
                            stats.top_branching_exons.push_back(info);
                        } else if (info.branches > stats.top_branching_exons.back().branches) {
                            stats.top_branching_exons.back() = info;
                        }
                        std::sort(stats.top_branching_exons.begin(),
                                  stats.top_branching_exons.end(),
                                  [](const auto& a, const auto& b) {
                                      return a.branches > b.branches;
                                  });
                    }
                }
            }
        }
    }

    // Initialize per-sample and per-source tracking (needed by Phase 3 and Phase 5)
    std::vector<uint32_t> sample_ids;
    for (size_t i = 0; i < registry.size(); ++i) {
        sample_ids.push_back(static_cast<uint32_t>(i));
        stats.per_sample[static_cast<uint32_t>(i)] = {};
    }

    std::unordered_map<uint32_t, std::unordered_set<std::string>> sample_gene_ids;
    std::unordered_map<uint32_t, std::unordered_set<std::string>> sample_transcript_ids;
    std::unordered_map<std::string, std::unordered_set<std::string>> source_gene_ids;

    // --- Phase 3: Compute exon sharing and per-sample/per-source exon stats ---
    // Single traversal with global deduplication for both sharing stats and per-sample/per-source counts
    std::unordered_map<std::string, size_t> gene_transcript_counts;
    for (const auto& [gene_id, gi] : genes) {
        gene_transcript_counts[gene_id] = gi.transcript_ids.size();
    }

    // Count only experimental samples (type == "sample") for conserved/sharing stats
    size_t total_samples = 0;
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto* info = registry.get(static_cast<uint32_t>(i));
        if (info && info->type == "sample") {
            total_samples++;
        }
    }
    stats.total_samples = total_samples;

    std::set<const void*> visited_exons_stats;
    size_t total_tx_per_exon = 0;

    for (auto& [seqid, seg_key] : segment_keys) {
        auto first_exons = grove.get_neighbors_if(seg_key,
            [](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON;
            });

        for (auto* exon_key : first_exons) {
            std::vector<key_ptr> to_visit;
            to_visit.push_back(exon_key);

            while (!to_visit.empty()) {
                auto* current = to_visit.back();
                to_visit.pop_back();

                if (!visited_exons_stats.insert(current).second) continue;

                auto& feature = current->get_data();
                if (!is_exon(feature)) continue;

                auto& exon = get_exon(feature);

                // Exon sharing stats
                size_t tx_count = exon.transcript_ids.size();
                total_tx_per_exon += tx_count;

                if (tx_count > 1) {
                    stats.shared_exons++;
                }
                if (tx_count > stats.max_transcripts_per_exon) {
                    stats.max_transcripts_per_exon = tx_count;
                }

                // Constitutive vs alternative
                size_t gene_total = gene_transcript_counts[exon.gene_id];
                if (gene_total > 0 && tx_count >= gene_total) {
                    stats.constitutive_exons++;
                } else {
                    stats.alternative_exons++;
                }

                // Per-sample exon stats
                for (uint32_t sid : exon.sample_idx) {
                    stats.per_sample[sid].exons++;
                    if (exon.sample_idx.size() == 1) {
                        stats.per_sample[sid].exclusive_exons++;
                    }
                }

                // Per-source exon stats
                for (const auto& src : exon.sources) {
                    stats.per_source[src].exons++;
                    if (exon.sources.size() == 1) {
                        stats.per_source[src].exclusive_exons++;
                    }
                    source_gene_ids[src].insert(exon.gene_id);
                }

                // Continue chain
                auto next = grove.get_neighbors_if(current,
                    [](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::EXON_TO_EXON;
                    });
                for (auto* n : next) {
                    if (!visited_exons_stats.count(n)) {
                        to_visit.push_back(n);
                    }
                }
            }
        }
    }

    if (stats.total_exons > 0) {
        stats.mean_transcripts_per_exon =
            static_cast<double>(total_tx_per_exon) / static_cast<double>(stats.total_exons);
    }

    // --- Phase 4: Compute gene-level statistics ---
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

    // Deduplication ratio
    if (stats.total_transcripts > 0) {
        stats.deduplication_ratio =
            static_cast<double>(stats.total_segments) / static_cast<double>(stats.total_transcripts);
    }

    // --- Phase 4b: Structural diversity (pairwise Jaccard distance of exon sets) ---
    // For each segment, walk its exon chain to collect the exon pointer set.
    // Group by gene, then compute mean pairwise Jaccard distance.
    struct segment_exon_info {
        std::unordered_set<uint32_t> sample_idx;
        std::set<const void*> exon_set;
    };
    std::unordered_map<std::string, std::vector<segment_exon_info>> gene_segment_exons;

    for (auto& [seqid, seg_key] : segment_keys) {
        auto& seg = get_segment(seg_key->get_data());
        if (seg.transcript_ids.empty()) continue;

        // Pick one transcript to follow its specific exon chain
        const std::string& tx_id = *seg.transcript_ids.begin();

        segment_exon_info sei;
        sei.sample_idx = seg.sample_idx;

        // Get first exon via SEGMENT_TO_EXON edge for this transcript
        auto first_exons = grove.get_neighbors_if(seg_key,
            [&tx_id](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON && e.id == tx_id;
            });

        if (first_exons.empty()) continue;

        // Walk exon chain following this transcript's EXON_TO_EXON edges
        key_ptr current = first_exons.front();
        sei.exon_set.insert(current);

        bool found_next = true;
        while (found_next) {
            found_next = false;
            auto next = grove.get_neighbors_if(current,
                [&tx_id](const edge_metadata& e) {
                    return e.type == edge_metadata::edge_type::EXON_TO_EXON && e.id == tx_id;
                });
            if (!next.empty()) {
                current = next.front();
                sei.exon_set.insert(current);
                found_next = true;
            }
        }

        gene_segment_exons[seg.gene_id].push_back(std::move(sei));
    }

    // Compute global mean pairwise Jaccard distance (averaged per gene, then across genes)
    {
        double global_sum = 0;
        size_t gene_count = 0;

        for (auto& [gene_id, seg_infos] : gene_segment_exons) {
            if (seg_infos.size() < 2) continue;

            double gene_sum = 0;
            size_t gene_pairs = 0;

            for (size_t i = 0; i < seg_infos.size(); ++i) {
                for (size_t j = i + 1; j < seg_infos.size(); ++j) {
                    gene_sum += jaccard_distance(seg_infos[i].exon_set, seg_infos[j].exon_set);
                    gene_pairs++;
                }
            }

            if (gene_pairs > 0) {
                global_sum += gene_sum / static_cast<double>(gene_pairs);
                gene_count++;
            }
        }

        stats.multi_segment_genes = gene_count;
        if (gene_count > 0) {
            stats.isoform_diversity = global_sum / static_cast<double>(gene_count);
        }
    }

    // Per-sample Jaccard diversity
    for (auto& [sid, ss] : stats.per_sample) {
        double sample_sum = 0;
        size_t sample_genes = 0;

        for (auto& [gene_id, seg_infos] : gene_segment_exons) {
            // Filter to segments containing this sample
            std::vector<size_t> indices;
            for (size_t k = 0; k < seg_infos.size(); ++k) {
                if (seg_infos[k].sample_idx.contains(sid)) {
                    indices.push_back(k);
                }
            }
            if (indices.size() < 2) continue;

            double gene_sum = 0;
            size_t gene_pairs = 0;
            for (size_t i = 0; i < indices.size(); ++i) {
                for (size_t j = i + 1; j < indices.size(); ++j) {
                    gene_sum += jaccard_distance(
                        seg_infos[indices[i]].exon_set,
                        seg_infos[indices[j]].exon_set);
                    gene_pairs++;
                }
            }

            if (gene_pairs > 0) {
                sample_sum += gene_sum / static_cast<double>(gene_pairs);
                sample_genes++;
            }
        }

        if (sample_genes > 0) {
            ss.isoform_diversity = sample_sum / static_cast<double>(sample_genes);
        }
    }

    // Exons per segment distribution
    if (!exons_per_segment.empty()) {
        size_t total_exons_sum = std::accumulate(exons_per_segment.begin(), exons_per_segment.end(), size_t{0});
        stats.mean_exons_per_segment =
            static_cast<double>(total_exons_sum) / static_cast<double>(exons_per_segment.size());
        stats.median_exons_per_segment = compute_median(exons_per_segment);
        stats.max_exons_per_segment = *std::max_element(exons_per_segment.begin(), exons_per_segment.end());
    }

    // Per-chromosome gene counts
    for (const auto& [seqid, gene_set] : chr_gene_ids) {
        stats.per_chromosome[seqid].genes = gene_set.size();
    }

    // --- Phase 5: Per-sample and per-source segment statistics + sharing distribution ---
    for (auto& [seqid, seg_key] : segment_keys) {
        auto& seg = get_segment(seg_key->get_data());

        // Segment sharing distribution
        if (seg.is_conserved(total_samples)) {
            stats.conserved_segments++;
        } else if (seg.is_sample_specific()) {
            stats.sample_specific_segments++;
        } else {
            stats.shared_segments++;
        }

        for (uint32_t sid : seg.sample_idx) {
            auto& ss = stats.per_sample[sid];
            ss.segments++;

            if (seg.sample_idx.size() == 1) {
                ss.exclusive_segments++;
            }

            sample_gene_ids[sid].insert(seg.gene_id);
            for (const auto& tx : seg.transcript_ids) {
                sample_transcript_ids[sid].insert(tx);
            }

            if (seg.has_expression(sid)) {
                ss.mean_expression += seg.get_expression(sid);
                ss.expressed_segments++;
            }
        }

        // Per-source segment stats
        for (const auto& src : seg.sources) {
            stats.per_source[src].segments++;
            if (seg.sources.size() == 1) {
                stats.per_source[src].exclusive_segments++;
            }
            source_gene_ids[src].insert(seg.gene_id);
        }
    }

    // Finalize per-sample stats
    for (auto& [sid, ss] : stats.per_sample) {
        ss.genes = sample_gene_ids[sid].size();
        ss.transcripts = sample_transcript_ids[sid].size();
        if (ss.transcripts > 0) {
            ss.deduplication_ratio =
                static_cast<double>(ss.segments) / static_cast<double>(ss.transcripts);
        }
        if (ss.expressed_segments > 0) {
            ss.mean_expression /= static_cast<double>(ss.expressed_segments);
        }
    }

    // Finalize per-source stats
    for (auto& [src, ss] : stats.per_source) {
        ss.genes = source_gene_ids[src].size();
    }

    // Total edges from grove
    stats.total_edges = grove.edge_count();

    return stats;
}

void index_stats::write(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open stats file: " + path);
        return;
    }

    out << std::fixed << std::setprecision(2);

    // Section header helper
    auto section = [&](const std::string& name) {
        std::string line = "== " + name + " ";
        while (line.size() < 72) line += '=';
        out << "\n" << line << "\n";
    };

    out << "Pan-Transcriptome Index Statistics\n";

    // Input entries
    if (!annotation_sources.empty()) {
        out << "\nInput entries:\n";
        auto& reg = sample_registry::instance();
        for (size_t i = 0; i < reg.size(); ++i) {
            const auto* info = reg.get(static_cast<uint32_t>(i));
            if (info && !info->id.empty()) {
                out << "  " << info->id << "  [" << info->type << "]\n";
            }
        }
    }

    // ── Overview ──
    section("Overview");
    out << "Chromosomes:        " << total_chromosomes << "\n";
    out << "Genes:              " << total_genes << "\n";
    out << "Transcripts:        " << total_transcripts
        << "  (unique transcript IDs across all inputs)\n";
    out << "Segments:           " << total_segments
        << "  (deduplicated by exon structure)\n";
    out << "Unique exons:       " << total_exons << "\n";
    out << "Graph edges:        " << total_edges << "\n";
    out << "\n";

    out << "Transcripts per gene:\n";
    out << "  Mean:             " << mean_transcripts_per_gene << "\n";
    out << "  Median:           " << median_transcripts_per_gene << "\n";
    out << "  Max:              " << max_transcripts_per_gene;
    if (!max_transcripts_gene_id.empty()) {
        out << "  (" << max_transcripts_gene_id << ")";
    }
    out << "\n";
    out << "  Single-isoform:   " << single_isoform_genes << "\n";
    out << "  Multi-isoform:    " << multi_isoform_genes << "\n";
    out << "\n";

    out << "Exons per segment:\n";
    out << "  Mean:             " << mean_exons_per_segment << "\n";
    out << "  Median:           " << median_exons_per_segment << "\n";
    out << "  Max:              " << max_exons_per_segment << "\n";
    out << "  Single-exon:      " << single_exon_segments << "\n";

    // ── Genes by biotype ──
    if (!genes_by_biotype.empty()) {
        section("Genes by biotype");

        std::vector<std::pair<std::string, size_t>> sorted_biotypes(
            genes_by_biotype.begin(), genes_by_biotype.end());
        std::sort(sorted_biotypes.begin(), sorted_biotypes.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        for (const auto& [biotype, count] : sorted_biotypes) {
            double pct = total_genes > 0 ?
                100.0 * static_cast<double>(count) / static_cast<double>(total_genes) : 0;
            out << "  " << std::left << std::setw(30) << biotype
                << std::right << std::setw(8) << count
                << "  (" << std::setprecision(1) << pct << "%)\n";
        }
        out << std::setprecision(2);
    }

    // ── Graph structure ──
    section("Graph structure");
    out << "   Edges connect exons within segments (EXON_TO_EXON) and segments to\n";
    out << "   their first exon (SEGMENT_TO_EXON). Branching exons have >1 downstream target.\n\n";
    out << "Total edges:        " << total_edges << "\n";
    out << "Branching exons:    " << branching_exons << "\n";

    if (!top_branching_exons.empty()) {
        out << "\nTop splicing hubs:\n";
        out << "   Exons with the most downstream splice choices.\n\n";
        out << std::left << std::setw(20) << "  Gene"
            << std::setw(16) << "Gene ID"
            << std::setw(30) << "Coordinate"
            << std::right << std::setw(10) << "Branches"
            << std::setw(12) << "Transcripts" << "\n";

        for (const auto& be : top_branching_exons) {
            std::string gene_label = be.gene_name.empty() ? be.gene_id : be.gene_name;
            out << "  " << std::left << std::setw(18) << gene_label
                << std::setw(16) << be.gene_id
                << std::setw(30) << be.coordinate
                << std::right << std::setw(10) << be.branches
                << std::setw(12) << be.transcripts << "\n";
        }
    }

    // ── Exon sharing ──
    section("Exon sharing");
    out << "   Exons shared across transcripts within the same gene.\n";
    out << "   Constitutive = in all transcripts of a gene; alternative = subset only.\n\n";

    out << "Shared exons:       " << shared_exons;
    if (total_exons > 0) {
        double pct = 100.0 * static_cast<double>(shared_exons) / static_cast<double>(total_exons);
        out << "  (" << std::setprecision(1) << pct << "% of " << total_exons << ")";
    }
    out << "\n" << std::setprecision(2);
    out << "Mean tx/exon:       " << mean_transcripts_per_exon << "\n";
    out << "Max tx/exon:        " << max_transcripts_per_exon << "\n";
    out << "Constitutive:       " << constitutive_exons << "\n";
    out << "Alternative:        " << alternative_exons << "\n";

    // Per-entry exon breakdown
    if (!per_sample.empty()) {
        auto& reg_exon = sample_registry::instance();
        std::vector<uint32_t> exon_ids;
        for (const auto& [sid, _] : per_sample) {
            exon_ids.push_back(sid);
        }
        std::sort(exon_ids.begin(), exon_ids.end());

        out << "\nPer entry (sharing across inputs):\n";
        out << "   Exclusive = exons only in this entry, not shared with any other.\n\n";
        for (uint32_t sid : exon_ids) {
            const auto& ss = per_sample.at(sid);
            const auto* info = reg_exon.get(sid);
            std::string label = (info && !info->id.empty()) ? info->id : std::to_string(sid);
            if (info && info->type != "sample") {
                label += " [" + info->type + "]";
            }
            double pct = ss.exons > 0 ?
                100.0 * static_cast<double>(ss.exclusive_exons) / static_cast<double>(ss.exons) : 0;
            out << "  " << std::left << std::setw(30) << label
                << std::right << std::setw(8) << ss.exons << " exons"
                << std::setw(10) << ss.exclusive_exons << " exclusive"
                << "  (" << std::setprecision(1) << pct << "%)\n";
        }
        out << std::setprecision(2);
    }

    // ── Segment sharing ──
    if (!per_sample.empty()) {
        section("Segment sharing");
        out << "   Distribution across experimental samples (annotations excluded).\n";
        out << "   Conserved = all samples; shared = 2+ but not all; sample-specific = exactly 1.\n\n";

        out << "Samples:            " << total_samples
            << "  (" << per_sample.size() << " entries total, "
            << (per_sample.size() - total_samples) << " annotations)\n";
        out << "Conserved:          " << conserved_segments;
        if (total_segments > 0) {
            double pct = 100.0 * static_cast<double>(conserved_segments) / static_cast<double>(total_segments);
            out << "  (" << std::setprecision(1) << pct << "%)";
        }
        out << "\n" << std::setprecision(2);
        out << "Shared (2+):        " << shared_segments;
        if (total_segments > 0) {
            double pct = 100.0 * static_cast<double>(shared_segments) / static_cast<double>(total_segments);
            out << "  (" << std::setprecision(1) << pct << "%)";
        }
        out << "\n" << std::setprecision(2);
        out << "Sample-specific:    " << sample_specific_segments;
        if (total_segments > 0) {
            double pct = 100.0 * static_cast<double>(sample_specific_segments) / static_cast<double>(total_segments);
            out << "  (" << std::setprecision(1) << pct << "%)";
        }
        out << "\n" << std::setprecision(2);
    }

    // ── Per sample ──
    if (!per_sample.empty()) {
        auto& registry = sample_registry::instance();

        std::vector<uint32_t> sorted_ids;
        for (const auto& [sid, _] : per_sample) {
            sorted_ids.push_back(sid);
        }
        std::sort(sorted_ids.begin(), sorted_ids.end());

        section("Per sample");
        out << "   Exclusive = features only in this entry.\n";
        out << "   Diversity = mean pairwise Jaccard distance (0=identical, 1=disjoint).\n";
        out << "   Dedup = segments/transcripts (1.0 = all unique).\n\n";

        out << std::left << std::setw(30) << "Sample"
            << std::right << std::setw(8) << "Genes"
            << std::setw(10) << "Segments"
            << std::setw(10) << "Exclusive"
            << std::setw(10) << "Exons"
            << std::setw(10) << "Exclusive"
            << std::setw(10) << "Tx"
            << std::setw(10) << "Diversity"
            << std::setw(10) << "Dedup"
            << std::setw(12) << "Expr.Segs"
            << std::setw(12) << "Mean.Expr" << "\n";

        for (uint32_t sid : sorted_ids) {
            const auto& ss = per_sample.at(sid);
            const auto* info = registry.get(sid);
            std::string label = (info && !info->id.empty()) ? info->id : std::to_string(sid);
            if (info && info->type != "sample") {
                label += " [" + info->type + "]";
            }

            out << std::left << std::setw(30) << label
                << std::right << std::setw(8) << ss.genes
                << std::setw(10) << ss.segments
                << std::setw(10) << ss.exclusive_segments
                << std::setw(10) << ss.exons
                << std::setw(10) << ss.exclusive_exons
                << std::setw(10) << ss.transcripts
                << std::setw(10) << std::setprecision(2) << ss.isoform_diversity
                << std::setw(10) << std::setprecision(3) << ss.deduplication_ratio;

            if (ss.expressed_segments > 0) {
                out << std::setw(12) << std::setprecision(2) << ss.mean_expression;
            } else {
                out << std::setw(12) << ".";
            }
            out << "\n";
        }
    }

    // ── Per source ──
    if (!per_source.empty()) {
        section("Per source");
        out << "   Features by annotation source (GFF column 2: HAVANA, ENSEMBL, TALON, etc.).\n";
        out << "   Exclusive = only annotated by this source.\n\n";

        out << std::left << std::setw(20) << "Source"
            << std::right << std::setw(8) << "Genes"
            << std::setw(10) << "Segments"
            << std::setw(10) << "Exclusive"
            << std::setw(10) << "Exons"
            << std::setw(10) << "Exclusive" << "\n";

        for (const auto& [src, ss] : per_source) {
            out << std::left << std::setw(20) << src
                << std::right << std::setw(8) << ss.genes
                << std::setw(10) << ss.segments
                << std::setw(10) << ss.exclusive_segments
                << std::setw(10) << ss.exons
                << std::setw(10) << ss.exclusive_exons << "\n";
        }
    }

    // ── Per chromosome ──
    if (!per_chromosome.empty()) {
        section("Per chromosome");
        out << std::left << std::setw(12) << "Chromosome"
            << std::right << std::setw(8) << "Genes"
            << std::setw(10) << "Segments"
            << std::setw(10) << "Exons" << "\n";

        for (const auto& [seqid, cs] : per_chromosome) {
            out << std::left << std::setw(12) << seqid
                << std::right << std::setw(8) << cs.genes
                << std::setw(10) << cs.segments
                << std::setw(10) << cs.exons << "\n";
        }
    }

    logging::info("Index statistics written to: " + path);
}

void index_stats::write_sample_csv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open CSV stats file: " + path);
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    if (sample_ids.empty()) {
        logging::warning("No per-sample stats to write");
        return;
    }

    // Header row: metric, sample1, sample2, ...
    out << "metric";
    for (uint32_t sid : sample_ids) {
        const auto* info = registry.get(sid);
        out << "," << (info && !info->id.empty() ? info->id : std::to_string(sid));
    }
    out << "\n";

    // Metadata rows
    auto write_metadata_row = [&](const std::string& label, auto getter) {
        out << label;
        for (uint32_t sid : sample_ids) {
            const auto* info = registry.get(sid);
            out << "," << (info ? getter(info) : "");
        }
        out << "\n";
    };

    write_metadata_row("type", [](const sample_info* i) -> const std::string& { return i->type; });
    write_metadata_row("source_file", [](const sample_info* i) { return i->source_file.filename().string(); });
    write_metadata_row("assay", [](const sample_info* i) -> const std::string& { return i->assay; });
    write_metadata_row("biosample", [](const sample_info* i) -> const std::string& { return i->biosample; });
    write_metadata_row("condition", [](const sample_info* i) -> const std::string& { return i->condition; });
    write_metadata_row("species", [](const sample_info* i) -> const std::string& { return i->species; });
    write_metadata_row("platform", [](const sample_info* i) -> const std::string& { return i->platform; });
    write_metadata_row("pipeline", [](const sample_info* i) -> const std::string& { return i->pipeline; });

    // Stat rows
    auto write_stat_row = [&](const std::string& label, auto getter) {
        out << label;
        for (uint32_t sid : sample_ids) {
            auto it = per_sample.find(sid);
            out << "," << (it != per_sample.end() ? std::to_string(getter(it->second)) : "0");
        }
        out << "\n";
    };

    write_stat_row("segments", [](const sample_stats& s) { return s.segments; });
    write_stat_row("exclusive_segments", [](const sample_stats& s) { return s.exclusive_segments; });
    write_stat_row("exons", [](const sample_stats& s) { return s.exons; });
    write_stat_row("exclusive_exons", [](const sample_stats& s) { return s.exclusive_exons; });
    write_stat_row("genes", [](const sample_stats& s) { return s.genes; });
    write_stat_row("transcripts", [](const sample_stats& s) { return s.transcripts; });
    write_stat_row("expressed_segments", [](const sample_stats& s) { return s.expressed_segments; });

    // Float metrics (need special handling)
    auto write_float_row = [&](const std::string& label, auto getter, int precision) {
        out << label;
        for (uint32_t sid : sample_ids) {
            auto it = per_sample.find(sid);
            if (it != per_sample.end()) {
                double val = getter(it->second);
                if (val > 0) {
                    out << "," << std::fixed << std::setprecision(precision) << val;
                } else {
                    out << ",.";
                }
            } else {
                out << ",.";
            }
        }
        out << "\n";
    };

    write_float_row("isoform_diversity", [](const sample_stats& s) { return s.isoform_diversity; }, 2);
    write_float_row("deduplication_ratio", [](const sample_stats& s) { return s.deduplication_ratio; }, 3);
    write_float_row("mean_expression", [](const sample_stats& s) -> double {
        return s.expressed_segments > 0 ? s.mean_expression : 0;
    }, 2);

    logging::info("Per-sample CSV statistics written to: " + path);
}

void index_stats::write_source_csv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open source CSV stats file: " + path);
        return;
    }

    if (per_source.empty()) {
        logging::warning("No per-source stats to write");
        return;
    }

    // Collect source names in sorted order (std::map is already sorted)
    std::vector<std::string> source_names;
    for (const auto& [src, _] : per_source) {
        source_names.push_back(src);
    }

    // Header row: metric, source1, source2, ...
    out << "metric";
    for (const auto& src : source_names) {
        out << "," << src;
    }
    out << "\n";

    // Stat rows
    auto write_row = [&](const std::string& label, auto getter) {
        out << label;
        for (const auto& src : source_names) {
            auto it = per_source.find(src);
            out << "," << (it != per_source.end() ? std::to_string(getter(it->second)) : "0");
        }
        out << "\n";
    };

    write_row("segments", [](const source_stats& s) { return s.segments; });
    write_row("exclusive_segments", [](const source_stats& s) { return s.exclusive_segments; });
    write_row("exons", [](const source_stats& s) { return s.exons; });
    write_row("exclusive_exons", [](const source_stats& s) { return s.exclusive_exons; });
    write_row("genes", [](const source_stats& s) { return s.genes; });

    logging::info("Per-source CSV statistics written to: " + path);
}