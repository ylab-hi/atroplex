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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "sample_info.hpp"
#include "utility.hpp"

namespace {

std::string expr_type_label(sample_info::expression_type t) {
    switch (t) {
        case sample_info::expression_type::COUNTS: return "counts";
        case sample_info::expression_type::TPM:    return "TPM";
        case sample_info::expression_type::FPKM:   return "FPKM";
        case sample_info::expression_type::RPKM:   return "RPKM";
        case sample_info::expression_type::CPM:    return "CPM";
        default:                                   return "expression";
    }
}

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

index_stats index_stats::collect(grove_type& grove, const collect_options& opts) {
    index_stats stats;
    bool detailed = opts.detailed;
    bool streaming = !opts.output_dir.empty();

    // Collect annotation sources from sample_registry
    auto& registry = sample_registry::instance();
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto* info = registry.get(static_cast<uint32_t>(i));
        if (info && !info->id.empty()) {
            stats.annotation_sources.push_back(info->id);
        }
    }

    // --- Phase 1: Traverse B+ tree to collect all segments ---
    // gene_id -> { set of transcript_ids (interned), biotype }
    struct gene_stats_info {
        std::unordered_set<uint32_t> transcript_ids;
        std::string biotype;
    };
    std::unordered_map<std::string, gene_stats_info> genes;
    std::vector<size_t> exons_per_segment;

    // Track per-chromosome gene sets
    std::unordered_map<std::string, std::unordered_set<std::string>> chr_gene_ids;

    // All segment key pointers (needed for exon traversal)
    std::vector<std::pair<std::string, key_ptr>> segment_keys; // (seqid, key)

    // Transcript -> sample mapping (interned transcript ID -> sample IDs)
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> tx_to_samples;

    // Transcript -> biotype mapping (interned transcript ID -> biotype)
    std::unordered_map<uint32_t, std::string> tx_biotypes;

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
                if (seg.absorbed) continue;  // Skip tombstoned segments
                stats.total_segments++;
                stats.per_chromosome[seqid].segments++;

                // Collect gene info
                auto& gi = genes[seg.gene_id()];
                gi.biotype = seg.gene_biotype();
                for (const auto& tx : seg.transcript_ids) {
                    gi.transcript_ids.insert(tx);
                    // Build transcript -> sample mapping (for per-sample hub stats)
                    for (uint32_t sid : seg.sample_idx) {
                        tx_to_samples[tx].insert(sid);
                    }
                }
                chr_gene_ids[seqid].insert(seg.gene_id());

                // Collect transcript biotypes
                for (const auto& [tx_id, biotype] : seg.transcript_biotypes) {
                    if (!biotype.empty()) {
                        tx_biotypes[tx_id] = biotype;
                    }
                }

                exons_per_segment.push_back(static_cast<size_t>(seg.exon_count));
                if (seg.exon_count == 1) {
                    stats.single_exon_segments++;
                }

                segment_keys.emplace_back(seqid, key);
            }
            node = node->get_next();
        }
    }

    // Build per-sample per-gene transcript counts (for traditional PSI)
    // gene_id -> sample_id -> transcript count
    std::unordered_map<std::string, std::unordered_map<uint32_t, size_t>> gene_sample_tx;
    for (const auto& [gene_id, gi] : genes) {
        for (const auto& tx : gi.transcript_ids) {
            auto it = tx_to_samples.find(tx);
            if (it != tx_to_samples.end()) {
                for (uint32_t sid : it->second) {
                    gene_sample_tx[gene_id][sid]++;
                }
            }
        }
    }

    // --- Phase 2: Traverse graph from segments to collect exons ---
    std::set<const void*> visited_exons; // deduplicate by pointer

    // Streaming: open branch details file if output_dir is set
    std::unique_ptr<std::ofstream> branch_stream;
    std::vector<uint32_t> stream_sample_ids;
    std::vector<std::string> stream_labels;
    std::vector<bool> stream_is_sample;
    std::vector<std::string> stream_expr_labels;

    if (streaming) {
        // Prepare sample ID ordering for streaming output
        for (size_t i = 0; i < registry.size(); ++i) {
            stream_sample_ids.push_back(static_cast<uint32_t>(i));
        }
        for (uint32_t sid : stream_sample_ids) {
            const auto* info = registry.get(sid);
            stream_labels.push_back((info && !info->id.empty()) ? info->id : std::to_string(sid));
            bool is_samp = info && info->type == "sample";
            stream_is_sample.push_back(is_samp);
            stream_expr_labels.push_back(is_samp ? expr_type_label(info->expr_type) : "");
        }

        auto hubs_dir = std::filesystem::path(opts.output_dir) / "splicing_hubs";
        std::filesystem::create_directories(hubs_dir);
        auto bd_path = hubs_dir / (opts.basename + ".branch_details.tsv");
        branch_stream = std::make_unique<std::ofstream>(bd_path.string());

        if (branch_stream->is_open()) {
            *branch_stream << "hub_gene_name\thub_gene_id\thub_exon_id\thub_coordinate"
                           << "\ttarget_exon_id\ttarget_coordinate";
            for (size_t i = 0; i < stream_sample_ids.size(); ++i) {
                *branch_stream << "\t" << stream_labels[i] << ".fraction";
                if (stream_is_sample[i]) {
                    *branch_stream << "\t" << stream_labels[i] << "." << stream_expr_labels[i];
                }
            }
            *branch_stream << "\n";
            *branch_stream << std::fixed;
        }
    }

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

                    // Collect splicing hubs (exons with many downstream targets)
                    if (unique_targets.size() > index_stats::MIN_HUB_BRANCHES) {
                        auto& feature = current->get_data();
                        if (is_exon(feature)) {
                            auto& exon = get_exon(feature);
                            index_stats::branching_exon_info info;
                            info.gene_name = exon.gene_name();
                            info.gene_id = exon.gene_id();
                            info.exon_id = exon.id;
                            info.coordinate = format_coordinate(seqid, current->get_value());
                            info.exon_key = current;
                            info.branches = unique_targets.size();
                            info.sample_idx = exon.sample_idx;

                            // Count transcripts for this gene only (exons can accumulate
                            // transcript_ids from overlapping genes via deduplication)
                            auto gene_tx_set = genes.find(exon.gene_id());
                            if (gene_tx_set != genes.end()) {
                                size_t gene_filtered = 0;
                                for (const auto& tx : exon.transcript_ids) {
                                    if (gene_tx_set->second.transcript_ids.count(tx)) gene_filtered++;
                                }
                                info.transcripts = gene_filtered;
                            } else {
                                info.transcripts = exon.transcript_ids.size();
                            }

                            // Collect downstream targets with per-sample info
                            std::set<key_ptr> seen_targets;
                            for (auto* n : next) {
                                if (!seen_targets.insert(n).second) continue;
                                auto& tf = n->get_data();
                                if (is_exon(tf)) {
                                    auto& te = get_exon(tf);

                                    // Per-sample branch counts
                                    for (uint32_t sid : te.sample_idx) {
                                        info.sample_branches[sid]++;
                                    }

                                    // Store target detail with per-sample transcript counts
                                    index_stats::branching_exon_info::branch_target bt;
                                    bt.exon_id = te.id;
                                    bt.coordinate = format_coordinate(seqid, n->get_value());
                                    bt.sample_idx = te.sample_idx;

                                    // Count transcripts through this target per sample
                                    // Only count transcripts that also use the hub exon AND belong to same gene
                                    auto gene_filter = genes.find(exon.gene_id());
                                    for (const auto& tx : te.transcript_ids) {
                                        if (!exon.transcript_ids.count(tx)) continue;
                                        if (gene_filter != genes.end() && !gene_filter->second.transcript_ids.count(tx)) continue;
                                        auto it = tx_to_samples.find(tx);
                                        if (it != tx_to_samples.end()) {
                                            for (uint32_t sid : it->second) {
                                                bt.sample_transcripts[sid]++;
                                            }
                                        }
                                    }

                                    bt.sample_expression = te.expression.to_map();
                                    info.targets.push_back(std::move(bt));
                                }
                            }

                            // Hub exon expression
                            info.sample_expression = exon.expression.to_map();

                            // Classify branches as shared/unique per sample
                            for (uint32_t sid : exon.sample_idx) {
                                for (const auto& bt : info.targets) {
                                    if (!bt.sample_idx.test(sid)) continue;
                                    if (bt.sample_idx.count() > 1) {
                                        info.sample_shared[sid]++;
                                    } else {
                                        info.sample_unique[sid]++;
                                    }
                                }
                            }

                            // Count per-sample transcripts using this hub exon
                            // Filter to only transcripts belonging to this gene (exons can
                            // accumulate transcript_ids from overlapping genes via deduplication)
                            auto gene_it = genes.find(exon.gene_id());
                            for (const auto& tx : exon.transcript_ids) {
                                if (gene_it != genes.end() && !gene_it->second.transcript_ids.count(tx)) continue;
                                auto it = tx_to_samples.find(tx);
                                if (it != tx_to_samples.end()) {
                                    for (uint32_t sid : it->second) {
                                        info.sample_transcripts[sid]++;
                                    }
                                }
                            }

                            // Compute Shannon entropy per sample from branch PSI distribution
                            for (uint32_t sid : exon.sample_idx) {
                                auto tx_it = info.sample_transcripts.find(sid);
                                if (tx_it == info.sample_transcripts.end() || tx_it->second == 0) continue;
                                double total_tx = static_cast<double>(tx_it->second);

                                double entropy = 0.0;
                                for (const auto& bt : info.targets) {
                                    auto bt_tx_it = bt.sample_transcripts.find(sid);
                                    if (bt_tx_it == bt.sample_transcripts.end() || bt_tx_it->second == 0) continue;
                                    double psi = static_cast<double>(bt_tx_it->second) / total_tx;
                                    if (psi > 0.0) {
                                        entropy -= psi * std::log2(psi);
                                    }
                                }
                                info.sample_entropy[sid] = entropy;
                            }

                            // Traditional PSI: hub transcripts / gene transcripts per sample
                            auto gene_tx_it = gene_sample_tx.find(exon.gene_id());
                            if (gene_tx_it != gene_sample_tx.end()) {
                                for (uint32_t sid : exon.sample_idx) {
                                    auto hub_tx = info.sample_transcripts.find(sid);
                                    auto gene_tx = gene_tx_it->second.find(sid);
                                    if (hub_tx != info.sample_transcripts.end()
                                        && gene_tx != gene_tx_it->second.end()
                                        && gene_tx->second > 0) {
                                        info.sample_psi[sid] =
                                            static_cast<double>(hub_tx->second)
                                            / static_cast<double>(gene_tx->second);
                                    }
                                }
                            }

                            // Stream branch detail rows to file and drop targets to save memory
                            if (branch_stream && branch_stream->is_open()) {
                                for (const auto& target : info.targets) {
                                    *branch_stream << info.gene_name << "\t"
                                        << info.gene_id << "\t"
                                        << info.exon_id << "\t"
                                        << info.coordinate << "\t"
                                        << target.exon_id << "\t"
                                        << target.coordinate;

                                    for (size_t si = 0; si < stream_sample_ids.size(); ++si) {
                                        uint32_t sid = stream_sample_ids[si];
                                        if (!target.sample_idx.test(sid)) {
                                            *branch_stream << "\t.";
                                            if (stream_is_sample[si]) *branch_stream << "\t.";
                                        } else {
                                            auto bt_tx_it = target.sample_transcripts.find(sid);
                                            auto hub_tx_it = info.sample_transcripts.find(sid);
                                            if (bt_tx_it != target.sample_transcripts.end()
                                                && hub_tx_it != info.sample_transcripts.end()
                                                && hub_tx_it->second > 0) {
                                                double frac = static_cast<double>(bt_tx_it->second)
                                                            / static_cast<double>(hub_tx_it->second);
                                                *branch_stream << "\t" << std::setprecision(3) << frac;
                                            } else {
                                                *branch_stream << "\t0";
                                            }
                                            if (stream_is_sample[si]) {
                                                auto expr_it = target.sample_expression.find(sid);
                                                if (expr_it != target.sample_expression.end()) {
                                                    *branch_stream << "\t" << std::setprecision(2) << expr_it->second;
                                                } else {
                                                    *branch_stream << "\t.";
                                                }
                                            }
                                        }
                                    }
                                    *branch_stream << "\n";
                                }
                                info.targets.clear();
                                info.targets.shrink_to_fit();
                            }

                            stats.splicing_hubs.push_back(std::move(info));
                        }
                    }
                }
            }
        }
    }

    if (branch_stream) {
        branch_stream.reset();
        logging::info("Branch details streamed to file");
    }

    // Sort splicing hubs by branch count descending
    std::sort(stats.splicing_hubs.begin(), stats.splicing_hubs.end(),
              [](const auto& a, const auto& b) { return a.branches > b.branches; });

    // Enrich hubs with exon position (walk chain from a representative segment)
    if (!stats.splicing_hubs.empty()) {
        // Build gene -> segments index
        std::unordered_map<std::string, std::vector<std::pair<std::string, key_ptr>>> gene_segments;
        for (auto& [seqid, seg_key] : segment_keys) {
            auto& seg = get_segment(seg_key->get_data());
            gene_segments[seg.gene_id()].emplace_back(seqid, seg_key);
        }

        for (auto& hub : stats.splicing_hubs) {
            auto git = gene_segments.find(hub.gene_id);
            if (git == gene_segments.end()) continue;

            bool found = false;
            for (auto& [seqid, seg_key] : git->second) {
                if (found) break;
                auto& seg = get_segment(seg_key->get_data());
                if (seg.transcript_ids.empty()) continue;

                size_t edge_id = seg.segment_index;
                auto first = grove.get_neighbors_if(seg_key,
                    [edge_id](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON
                            && e.id == edge_id;
                    });
                if (first.empty()) continue;

                // Walk chain counting position
                size_t pos = 1;
                key_ptr cur = first.front();
                while (true) {
                    auto& ef = cur->get_data();
                    if (cur == hub.exon_key) {
                        hub.exon_number = pos;
                        hub.total_exons = seg.exon_count;
                        found = true;
                        break;
                    }
                    auto nxt = grove.get_neighbors_if(cur,
                        [edge_id](const edge_metadata& e) {
                            return e.type == edge_metadata::edge_type::EXON_TO_EXON
                                && e.id == edge_id;
                        });
                    if (nxt.empty()) break;
                    cur = nxt.front();
                    pos++;
                }
            }
        }
    }

    // Free gene_sample_tx — only needed for Phase 2 PSI computation
    gene_sample_tx.clear();

    // Initialize per-sample and per-source tracking (needed by Phase 3 and Phase 5)
    std::vector<uint32_t> sample_ids;
    for (size_t i = 0; i < registry.size(); ++i) {
        sample_ids.push_back(static_cast<uint32_t>(i));
        stats.per_sample[static_cast<uint32_t>(i)] = {};
    }

    std::unordered_map<uint32_t, std::unordered_set<std::string>> sample_gene_ids;
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> sample_transcript_ids;
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

    // Compute global + per-sample transcript biotype counts
    for (const auto& [tx_id, biotype] : tx_biotypes) {
        stats.transcripts_by_biotype[biotype]++;
        auto it = tx_to_samples.find(tx_id);
        if (it != tx_to_samples.end()) {
            for (uint32_t sid : it->second) {
                stats.per_sample[sid].transcripts_by_biotype[biotype]++;
            }
        }
    }

    // Streaming: open conserved exons file if output_dir is set
    std::unique_ptr<std::ofstream> conserved_stream;
    if (streaming) {
        auto sharing_dir = std::filesystem::path(opts.output_dir) / "sharing";
        std::filesystem::create_directories(sharing_dir);
        auto ce_path = sharing_dir / (opts.basename + ".conserved_exons.tsv");
        conserved_stream = std::make_unique<std::ofstream>(ce_path.string());

        if (conserved_stream->is_open()) {
            *conserved_stream << "exon_id\tgene_name\tgene_id\tchromosome\tcoordinate\tn_transcripts\tconstitutive";
            for (size_t i = 0; i < stream_sample_ids.size(); ++i) {
                *conserved_stream << "\t" << stream_labels[i] << ".transcripts";
                if (stream_is_sample[i]) {
                    *conserved_stream << "\t" << stream_labels[i] << "." << stream_expr_labels[i];
                }
            }
            *conserved_stream << "\n";
            *conserved_stream << std::fixed;
        }
    }

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
                size_t gene_total = gene_transcript_counts[exon.gene_id()];
                if (gene_total > 0 && tx_count >= gene_total) {
                    stats.constitutive_exons++;
                } else {
                    stats.alternative_exons++;
                }

                // Per-sample exon stats (cross-sample sharing + transcript-level)
                bool is_conserved = exon.sample_idx.count() >= total_samples && total_samples > 0;
                for (uint32_t sid : exon.sample_idx) {
                    auto& ss = stats.per_sample[sid];
                    ss.exons++;
                    if (exon.sample_idx.count() == 1) {
                        ss.exclusive_exons++;
                    } else if (is_conserved) {
                        ss.conserved_exons++;
                    } else {
                        ss.shared_exons++;
                    }
                    if (gene_total > 0 && tx_count >= gene_total) {
                        ss.constitutive_exons++;
                    } else {
                        ss.alternative_exons++;
                    }
                }

                // Collect conserved exon details
                if (is_conserved) {
                    bool constitutive = (gene_total > 0 && tx_count >= gene_total);

                    // Build per-sample transcript counts for this exon
                    std::unordered_map<uint32_t, size_t> entry_sample_tx;
                    for (const auto& tx : exon.transcript_ids) {
                        auto it = tx_to_samples.find(tx);
                        if (it != tx_to_samples.end()) {
                            for (uint32_t sid : it->second) {
                                entry_sample_tx[sid]++;
                            }
                        }
                    }

                    if (conserved_stream && conserved_stream->is_open()) {
                        // Stream directly to file
                        *conserved_stream << exon.id << "\t"
                            << exon.gene_name() << "\t"
                            << exon.gene_id() << "\t"
                            << seqid << "\t"
                            << format_coordinate(seqid, current->get_value()) << "\t"
                            << tx_count << "\t"
                            << (constitutive ? "yes" : "no");

                        for (size_t si = 0; si < stream_sample_ids.size(); ++si) {
                            uint32_t sid = stream_sample_ids[si];
                            auto tx_it = entry_sample_tx.find(sid);
                            *conserved_stream << "\t" << (tx_it != entry_sample_tx.end() ? tx_it->second : size_t{0});
                            if (stream_is_sample[si]) {
                                if (exon.has_expression(sid)) {
                                    *conserved_stream << "\t" << std::setprecision(2) << exon.get_expression(sid);
                                } else {
                                    *conserved_stream << "\t.";
                                }
                            }
                        }
                        *conserved_stream << "\n";
                    } else {
                        // Accumulate in memory (non-streaming path)
                        index_stats::conserved_exon_entry entry;
                        entry.exon_id = exon.id;
                        entry.gene_name = exon.gene_name();
                        entry.gene_id = exon.gene_id();
                        entry.chromosome = seqid;
                        entry.coordinate = format_coordinate(seqid, current->get_value());
                        entry.n_transcripts = tx_count;
                        entry.constitutive = constitutive;
                        entry.sample_transcripts = std::move(entry_sample_tx);
                        entry.sample_expression = exon.expression.to_map();
                        stats.conserved_exon_details.push_back(std::move(entry));
                    }
                }

                // Per-source exon stats
                source_registry::instance().for_each(exon.sources, [&](const std::string& src) {
                    stats.per_source[src].exons++;
                    if (exon.source_count() == 1) {
                        stats.per_source[src].exclusive_exons++;
                    }
                    source_gene_ids[src].insert(exon.gene_id());
                });

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

    // Close conserved exon stream and free tx_to_samples (no longer needed)
    if (conserved_stream) {
        conserved_stream.reset();
        logging::info("Conserved exons streamed to file");
    }
    tx_to_samples.clear();
    tx_biotypes.clear();

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
    // Only computed in detailed mode (analyze subcommand). Skipped for --stats.
    if (detailed) {
    // For each segment, walk its exon chain to collect the exon pointer set.
    // Group by gene, then compute mean pairwise Jaccard distance.
    struct segment_exon_info {
        sample_bitset sample_idx;
        std::set<const void*> exon_set;
    };
    std::unordered_map<std::string, std::vector<segment_exon_info>> gene_segment_exons;

    for (auto& [seqid, seg_key] : segment_keys) {
        auto& seg = get_segment(seg_key->get_data());

        size_t edge_id = seg.segment_index;

        segment_exon_info sei;
        sei.sample_idx = seg.sample_idx;

        // Get first exon via SEGMENT_TO_EXON edge (filter by segment index)
        auto first_exons = grove.get_neighbors_if(seg_key,
            [edge_id](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON && e.id == edge_id;
            });

        if (first_exons.empty()) continue;

        // Walk exon chain following segment's EXON_TO_EXON edges
        key_ptr current = first_exons.front();
        sei.exon_set.insert(current);

        bool found_next = true;
        while (found_next) {
            found_next = false;
            auto next = grove.get_neighbors_if(current,
                [edge_id](const edge_metadata& e) {
                    return e.type == edge_metadata::edge_type::EXON_TO_EXON && e.id == edge_id;
                });
            if (!next.empty()) {
                current = next.front();
                sei.exon_set.insert(current);
                found_next = true;
            }
        }

        gene_segment_exons[seg.gene_id()].push_back(std::move(sei));
    }

    // Compute global mean pairwise Jaccard distance (averaged per gene, then across genes)
    // Cap at MAX_JACCARD_SEGMENTS to avoid O(n²) blowup on highly fragmented genes
    static constexpr size_t MAX_JACCARD_SEGMENTS = 200;
    size_t skipped_genes = 0;

    {
        double global_sum = 0;
        size_t gene_count = 0;

        for (auto& [gene_id, seg_infos] : gene_segment_exons) {
            if (seg_infos.size() < 2) continue;
            if (seg_infos.size() > MAX_JACCARD_SEGMENTS) {
                skipped_genes++;
                continue;
            }

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

    if (skipped_genes > 0) {
        logging::info("Jaccard diversity: skipped " + std::to_string(skipped_genes)
            + " genes with >" + std::to_string(MAX_JACCARD_SEGMENTS) + " segments");
    }

    // Per-sample Jaccard diversity
    for (auto& [sid, ss] : stats.per_sample) {
        double sample_sum = 0;
        size_t sample_genes = 0;

        for (auto& [gene_id, seg_infos] : gene_segment_exons) {
            if (seg_infos.size() > MAX_JACCARD_SEGMENTS) continue;

            // Filter to segments containing this sample
            std::vector<size_t> indices;
            for (size_t k = 0; k < seg_infos.size(); ++k) {
                if (seg_infos[k].sample_idx.test(sid)) {
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

    // Explicitly free Jaccard data
    gene_segment_exons.clear();
    } // end if (detailed)

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

            if (seg.sample_idx.count() == 1) {
                ss.exclusive_segments++;
            } else if (seg.is_conserved(total_samples)) {
                ss.conserved_segments++;
            } else {
                ss.shared_segments++;
            }

            sample_gene_ids[sid].insert(seg.gene_id());
            for (const auto& tx : seg.transcript_ids) {
                sample_transcript_ids[sid].insert(tx);
            }

            if (seg.has_expression(sid)) {
                ss.mean_expression += seg.get_expression(sid);
                ss.expressed_segments++;
            }
        }

        // Per-source segment stats
        source_registry::instance().for_each(seg.sources, [&](const std::string& src) {
            stats.per_source[src].segments++;
            if (seg.source_count() == 1) {
                stats.per_source[src].exclusive_segments++;
            }
            source_gene_ids[src].insert(seg.gene_id());
        });
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
        // Per-sample gene biotype counts
        for (const auto& gene_id : sample_gene_ids[sid]) {
            auto it = genes.find(gene_id);
            if (it != genes.end() && !it->second.biotype.empty()) {
                ss.genes_by_biotype[it->second.biotype]++;
            }
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
    if (absorbed_segments > 0) {
        out << "Absorbed segments:  " << absorbed_segments
            << "  (ISM fragments merged into longer parents)\n";
    }
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

    // ── Transcripts by biotype ──
    if (!transcripts_by_biotype.empty()) {
        section("Transcripts by biotype");

        std::vector<std::pair<std::string, size_t>> sorted_tx_biotypes(
            transcripts_by_biotype.begin(), transcripts_by_biotype.end());
        std::sort(sorted_tx_biotypes.begin(), sorted_tx_biotypes.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        for (const auto& [biotype, count] : sorted_tx_biotypes) {
            double pct = total_transcripts > 0 ?
                100.0 * static_cast<double>(count) / static_cast<double>(total_transcripts) : 0;
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

    if (!splicing_hubs.empty()) {
        out << "\nTop splicing hubs (>" << MIN_HUB_BRANCHES << " branches):\n";
        out << "   Exons with the most downstream splice choices.\n\n";
        out << std::left << std::setw(20) << "  Gene"
            << std::setw(16) << "Gene ID"
            << std::setw(30) << "Coordinate"
            << std::right << std::setw(10) << "Branches"
            << std::setw(12) << "Transcripts" << "\n";

        size_t display_count = std::min(splicing_hubs.size(), MAX_DISPLAY_HUBS);
        for (size_t i = 0; i < display_count; ++i) {
            const auto& be = splicing_hubs[i];
            std::string gene_label = be.gene_name.empty() ? be.gene_id : be.gene_name;
            out << "  " << std::left << std::setw(18) << gene_label
                << std::setw(16) << be.gene_id
                << std::setw(30) << be.coordinate
                << std::right << std::setw(10) << be.branches
                << std::setw(12) << be.transcripts << "\n";
        }
        if (splicing_hubs.size() > MAX_DISPLAY_HUBS) {
            out << "  ... and " << (splicing_hubs.size() - MAX_DISPLAY_HUBS)
                << " more (see splicing_hubs.tsv)\n";
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

void index_stats::write_summary(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open summary file: " + path);
        return;
    }

    out << std::fixed << std::setprecision(2);

    out << "Index Summary\n";
    out << "=============\n\n";

    // Input entries
    if (!annotation_sources.empty()) {
        auto& reg = sample_registry::instance();
        out << "Inputs: " << reg.size() << " entries (";
        out << total_samples << " samples, "
            << (reg.size() - total_samples) << " annotations)\n";
        for (size_t i = 0; i < reg.size(); ++i) {
            const auto* info = reg.get(static_cast<uint32_t>(i));
            if (info && !info->id.empty()) {
                out << "  " << info->id << "  [" << info->type << "]\n";
            }
        }
        out << "\n";
    }

    // Overview
    out << "Chromosomes:        " << total_chromosomes << "\n";
    out << "Genes:              " << total_genes << "\n";
    out << "Transcripts:        " << total_transcripts << "\n";
    out << "Segments:           " << total_segments
        << "  (deduplicated by exon structure)\n";
    if (absorbed_segments > 0) {
        out << "Absorbed segments:  " << absorbed_segments
            << "  (ISM fragments merged into longer parents)\n";
    }
    out << "Unique exons:       " << total_exons << "\n";
    out << "Graph edges:        " << total_edges << "\n";
    if (deduplication_ratio > 0) {
        out << "Dedup ratio:        " << std::setprecision(3)
            << deduplication_ratio << std::setprecision(2)
            << "  (segments/transcripts)\n";
    }
    out << "\n";

    // Transcripts per gene
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

    // Exons per segment
    out << "Exons per segment:\n";
    out << "  Mean:             " << mean_exons_per_segment << "\n";
    out << "  Median:           " << median_exons_per_segment << "\n";
    out << "  Max:              " << max_exons_per_segment << "\n";
    out << "  Single-exon:      " << single_exon_segments << "\n";
    out << "\n";

    // Genes by biotype
    if (!genes_by_biotype.empty()) {
        out << "Genes by biotype:\n";
        std::vector<std::pair<std::string, size_t>> sorted_biotypes(
            genes_by_biotype.begin(), genes_by_biotype.end());
        std::sort(sorted_biotypes.begin(), sorted_biotypes.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });
        for (const auto& [biotype, count] : sorted_biotypes) {
            double pct = total_genes > 0 ?
                100.0 * static_cast<double>(count) / static_cast<double>(total_genes) : 0;
            out << "  " << std::left << std::setw(28) << biotype
                << std::right << std::setw(8) << count
                << "  (" << std::setprecision(1) << pct << "%)\n";
        }
        out << std::setprecision(2) << "\n";
    }

    // Transcripts by biotype
    if (!transcripts_by_biotype.empty()) {
        out << "Transcripts by biotype:\n";
        std::vector<std::pair<std::string, size_t>> sorted_tx_biotypes(
            transcripts_by_biotype.begin(), transcripts_by_biotype.end());
        std::sort(sorted_tx_biotypes.begin(), sorted_tx_biotypes.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });
        for (const auto& [biotype, count] : sorted_tx_biotypes) {
            double pct = total_transcripts > 0 ?
                100.0 * static_cast<double>(count) / static_cast<double>(total_transcripts) : 0;
            out << "  " << std::left << std::setw(28) << biotype
                << std::right << std::setw(8) << count
                << "  (" << std::setprecision(1) << pct << "%)\n";
        }
        out << std::setprecision(2) << "\n";
    }

    // Per chromosome
    if (!per_chromosome.empty()) {
        out << "Per chromosome:\n";
        out << "  " << std::left << std::setw(10) << "Chrom"
            << std::right << std::setw(8) << "Genes"
            << std::setw(10) << "Segments"
            << std::setw(10) << "Exons" << "\n";
        for (const auto& [seqid, cs] : per_chromosome) {
            out << "  " << std::left << std::setw(10) << seqid
                << std::right << std::setw(8) << cs.genes
                << std::setw(10) << cs.segments
                << std::setw(10) << cs.exons << "\n";
        }
        out << "\n";
    }

    out << "For detailed per-sample analysis, run: atroplex analyze\n";

    logging::info("Index summary written to: " + path);
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

void index_stats::write_splicing_hubs_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open splicing hubs file: " + path);
        return;
    }

    if (splicing_hubs.empty()) {
        logging::warning("No splicing hubs to write");
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    // Determine which entries are samples (for expression columns)
    std::vector<bool> is_sample_type;
    std::vector<std::string> labels;
    std::vector<std::string> expr_labels;
    for (uint32_t sid : sample_ids) {
        const auto* info = registry.get(sid);
        labels.push_back((info && !info->id.empty()) ? info->id : std::to_string(sid));
        bool is_sample = info && info->type == "sample";
        is_sample_type.push_back(is_sample);
        expr_labels.push_back(is_sample ? expr_type_label(info->expr_type) : "");
    }

    // Header — per sample: branches, shared, unique, transcripts, entropy, psi + expression for samples
    out << "gene_name\tgene_id\texon_id\tcoordinate\texon_number\ttotal_exons\ttotal_branches\ttotal_transcripts";
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        out << "\t" << labels[i] << ".branches"
            << "\t" << labels[i] << ".shared"
            << "\t" << labels[i] << ".unique"
            << "\t" << labels[i] << ".transcripts"
            << "\t" << labels[i] << ".entropy"
            << "\t" << labels[i] << ".psi";
        if (is_sample_type[i]) {
            out << "\t" << labels[i] << "." << expr_labels[i];
        }
    }
    out << "\n";

    // Rows
    out << std::fixed;
    for (const auto& hub : splicing_hubs) {
        out << hub.gene_name << "\t"
            << hub.gene_id << "\t"
            << hub.exon_id << "\t"
            << hub.coordinate << "\t"
            << hub.exon_number << "\t"
            << hub.total_exons << "\t"
            << hub.branches << "\t"
            << hub.transcripts;

        for (size_t i = 0; i < sample_ids.size(); ++i) {
            uint32_t sid = sample_ids[i];
            if (!hub.sample_idx.test(sid)) {
                out << "\t.\t.\t.\t.\t.\t.";
                if (is_sample_type[i]) out << "\t.";
            } else {
                auto br_it = hub.sample_branches.find(sid);
                auto sh_it = hub.sample_shared.find(sid);
                auto uq_it = hub.sample_unique.find(sid);
                auto tx_it = hub.sample_transcripts.find(sid);
                auto en_it = hub.sample_entropy.find(sid);
                auto psi_it = hub.sample_psi.find(sid);

                out << "\t" << (br_it != hub.sample_branches.end() ? br_it->second : 0)
                    << "\t" << (sh_it != hub.sample_shared.end() ? sh_it->second : 0)
                    << "\t" << (uq_it != hub.sample_unique.end() ? uq_it->second : 0)
                    << "\t" << (tx_it != hub.sample_transcripts.end() ? tx_it->second : 0)
                    << "\t";
                if (en_it != hub.sample_entropy.end()) {
                    out << std::setprecision(3) << en_it->second;
                } else {
                    out << ".";
                }
                out << "\t";
                if (psi_it != hub.sample_psi.end()) {
                    out << std::setprecision(3) << psi_it->second;
                } else {
                    out << ".";
                }
                if (is_sample_type[i]) {
                    auto expr_it = hub.sample_expression.find(sid);
                    if (expr_it != hub.sample_expression.end()) {
                        out << "\t" << std::setprecision(2) << expr_it->second;
                    } else {
                        out << "\t.";
                    }
                }
            }
        }
        out << "\n";
    }

    logging::info("Splicing hubs (" + std::to_string(splicing_hubs.size())
        + " exons with >" + std::to_string(MIN_HUB_BRANCHES)
        + " branches) written to: " + path);
}

void index_stats::write_branch_details_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open branch details file: " + path);
        return;
    }

    if (splicing_hubs.empty()) {
        logging::warning("No splicing hubs for branch details");
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    // Determine which entries are samples (for expression columns)
    std::vector<bool> is_sample_type;
    std::vector<std::string> labels;
    std::vector<std::string> expr_labels;
    for (uint32_t sid : sample_ids) {
        const auto* info = registry.get(sid);
        labels.push_back((info && !info->id.empty()) ? info->id : std::to_string(sid));
        bool is_sample = info && info->type == "sample";
        is_sample_type.push_back(is_sample);
        expr_labels.push_back(is_sample ? expr_type_label(info->expr_type) : "");
    }

    // Header — sample columns show branch usage fraction + expression for samples
    out << "hub_gene_name\thub_gene_id\thub_exon_id\thub_coordinate"
        << "\ttarget_exon_id\ttarget_coordinate";
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        out << "\t" << labels[i] << ".fraction";
        if (is_sample_type[i]) {
            out << "\t" << labels[i] << "." << expr_labels[i];
        }
    }
    out << "\n";

    // One row per (hub, target) pair
    out << std::fixed;
    for (const auto& hub : splicing_hubs) {
        for (const auto& target : hub.targets) {
            out << hub.gene_name << "\t"
                << hub.gene_id << "\t"
                << hub.exon_id << "\t"
                << hub.coordinate << "\t"
                << target.exon_id << "\t"
                << target.coordinate;

            for (size_t i = 0; i < sample_ids.size(); ++i) {
                uint32_t sid = sample_ids[i];
                if (!target.sample_idx.test(sid)) {
                    out << "\t.";
                    if (is_sample_type[i]) out << "\t.";
                } else {
                    auto bt_tx_it = target.sample_transcripts.find(sid);
                    auto hub_tx_it = hub.sample_transcripts.find(sid);
                    if (bt_tx_it != target.sample_transcripts.end()
                        && hub_tx_it != hub.sample_transcripts.end()
                        && hub_tx_it->second > 0) {
                        double frac = static_cast<double>(bt_tx_it->second)
                                    / static_cast<double>(hub_tx_it->second);
                        out << "\t" << std::setprecision(3) << frac;
                    } else {
                        out << "\t0";
                    }
                    if (is_sample_type[i]) {
                        auto expr_it = target.sample_expression.find(sid);
                        if (expr_it != target.sample_expression.end()) {
                            out << "\t" << std::setprecision(2) << expr_it->second;
                        } else {
                            out << "\t.";
                        }
                    }
                }
            }
            out << "\n";
        }
    }

    logging::info("Branch details written to: " + path);
}

void index_stats::write_exon_sharing_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open exon sharing file: " + path);
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    // Resolve sample labels
    std::vector<std::string> labels;
    for (uint32_t sid : sample_ids) {
        const auto* info = registry.get(sid);
        labels.push_back((info && !info->id.empty()) ? info->id : std::to_string(sid));
    }

    // Header
    out << "metric\ttotal";
    for (const auto& label : labels) {
        out << "\t" << label;
    }
    out << "\n";

    // Helper to write a row
    auto write_row = [&](const std::string& metric, size_t global_val,
                         auto per_sample_fn) {
        out << metric << "\t" << global_val;
        for (size_t i = 0; i < sample_ids.size(); ++i) {
            auto it = per_sample.find(sample_ids[i]);
            if (it != per_sample.end()) {
                out << "\t" << per_sample_fn(it->second);
            } else {
                out << "\t0";
            }
        }
        out << "\n";
    };

    write_row("total", total_exons,
        [](const sample_stats& s) { return s.exons; });
    write_row("exclusive", 0,
        [](const sample_stats& s) { return s.exclusive_exons; });
    write_row("shared", shared_exons,
        [](const sample_stats& s) { return s.shared_exons; });
    write_row("conserved", 0,
        [](const sample_stats& s) { return s.conserved_exons; });
    write_row("constitutive", constitutive_exons,
        [](const sample_stats& s) { return s.constitutive_exons; });
    write_row("alternative", alternative_exons,
        [](const sample_stats& s) { return s.alternative_exons; });

    logging::info("Exon sharing statistics written to: " + path);
}

void index_stats::write_segment_sharing_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open segment sharing file: " + path);
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    // Resolve sample labels
    std::vector<std::string> labels;
    for (uint32_t sid : sample_ids) {
        const auto* info = registry.get(sid);
        labels.push_back((info && !info->id.empty()) ? info->id : std::to_string(sid));
    }

    // Header
    out << "metric\ttotal";
    for (const auto& label : labels) {
        out << "\t" << label;
    }
    out << "\n";

    // Helper to write a row
    auto write_row = [&](const std::string& metric, size_t global_val,
                         auto per_sample_fn) {
        out << metric << "\t" << global_val;
        for (size_t i = 0; i < sample_ids.size(); ++i) {
            auto it = per_sample.find(sample_ids[i]);
            if (it != per_sample.end()) {
                out << "\t" << per_sample_fn(it->second);
            } else {
                out << "\t0";
            }
        }
        out << "\n";
    };

    write_row("total", total_segments,
        [](const sample_stats& s) { return s.segments; });
    write_row("exclusive", sample_specific_segments,
        [](const sample_stats& s) { return s.exclusive_segments; });
    write_row("shared", shared_segments,
        [](const sample_stats& s) { return s.shared_segments; });
    write_row("conserved", conserved_segments,
        [](const sample_stats& s) { return s.conserved_segments; });

    logging::info("Segment sharing statistics written to: " + path);
}

void index_stats::write_conserved_exons_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open conserved exons file: " + path);
        return;
    }

    if (conserved_exon_details.empty()) {
        logging::warning("No conserved exons to write");
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    // Determine which entries are samples (have expression data possible)
    std::vector<bool> is_sample_type;
    std::vector<std::string> labels;
    std::vector<std::string> expr_labels;
    for (uint32_t sid : sample_ids) {
        const auto* info = registry.get(sid);
        labels.push_back((info && !info->id.empty()) ? info->id : std::to_string(sid));
        bool is_sample = info && info->type == "sample";
        is_sample_type.push_back(is_sample);
        expr_labels.push_back(is_sample ? expr_type_label(info->expr_type) : "");
    }

    // Header: .transcripts for all, .expression_type for sample-type entries only
    out << "exon_id\tgene_name\tgene_id\tchromosome\tcoordinate\tn_transcripts\tconstitutive";
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        out << "\t" << labels[i] << ".transcripts";
        if (is_sample_type[i]) {
            out << "\t" << labels[i] << "." << expr_labels[i];
        }
    }
    out << "\n";

    // Sort by chromosome then coordinate for stable output
    auto sorted = conserved_exon_details;
    std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) {
        if (a.chromosome != b.chromosome) return a.chromosome < b.chromosome;
        return a.coordinate < b.coordinate;
    });

    out << std::fixed;
    for (const auto& entry : sorted) {
        out << entry.exon_id << "\t"
            << entry.gene_name << "\t"
            << entry.gene_id << "\t"
            << entry.chromosome << "\t"
            << entry.coordinate << "\t"
            << entry.n_transcripts << "\t"
            << (entry.constitutive ? "yes" : "no");

        for (size_t i = 0; i < sample_ids.size(); ++i) {
            uint32_t sid = sample_ids[i];
            auto tx_it = entry.sample_transcripts.find(sid);
            out << "\t" << (tx_it != entry.sample_transcripts.end() ? tx_it->second : size_t{0});

            if (is_sample_type[i]) {
                auto expr_it = entry.sample_expression.find(sid);
                if (expr_it != entry.sample_expression.end()) {
                    out << "\t" << std::setprecision(2) << expr_it->second;
                } else {
                    out << "\t.";
                }
            }
        }
        out << "\n";
    }

    logging::info("Conserved exons (" + std::to_string(sorted.size())
        + " exons in all " + std::to_string(total_samples) + " samples) written to: " + path);
}

void index_stats::write_overview_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open overview file: " + path);
        return;
    }

    out << "metric\tvalue\n";

    out << "chromosomes\t" << total_chromosomes << "\n";
    out << "genes\t" << total_genes << "\n";
    out << "transcripts\t" << total_transcripts << "\n";
    out << "segments\t" << total_segments << "\n";
    out << "exons\t" << total_exons << "\n";
    out << "edges\t" << total_edges << "\n";
    out << "samples\t" << total_samples << "\n";

    out << std::fixed << std::setprecision(2);
    out << "mean_transcripts_per_gene\t" << mean_transcripts_per_gene << "\n";
    out << "median_transcripts_per_gene\t" << median_transcripts_per_gene << "\n";
    out << "max_transcripts_per_gene\t" << max_transcripts_per_gene << "\n";
    if (!max_transcripts_gene_id.empty()) {
        out << "max_transcripts_gene_id\t" << max_transcripts_gene_id << "\n";
    }
    out << "single_isoform_genes\t" << single_isoform_genes << "\n";
    out << "multi_isoform_genes\t" << multi_isoform_genes << "\n";

    out << "mean_exons_per_segment\t" << mean_exons_per_segment << "\n";
    out << "median_exons_per_segment\t" << median_exons_per_segment << "\n";
    out << "max_exons_per_segment\t" << max_exons_per_segment << "\n";
    out << "single_exon_segments\t" << single_exon_segments << "\n";

    out << "branching_exons\t" << branching_exons << "\n";
    out << "deduplication_ratio\t" << deduplication_ratio << "\n";

    if (isoform_diversity > 0) {
        out << "isoform_diversity\t" << isoform_diversity << "\n";
        out << "multi_segment_genes\t" << multi_segment_genes << "\n";
    }

    logging::info("Overview statistics written to: " + path);
}

void index_stats::write_per_chromosome_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open per-chromosome file: " + path);
        return;
    }

    out << "chromosome\tgenes\tsegments\texons\n";
    for (const auto& [seqid, cs] : per_chromosome) {
        out << seqid << "\t" << cs.genes << "\t" << cs.segments << "\t" << cs.exons << "\n";
    }

    logging::info("Per-chromosome statistics written to: " + path);
}

void index_stats::write_per_source_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open per-source file: " + path);
        return;
    }

    out << "source\tgenes\tsegments\texclusive_segments\texons\texclusive_exons\n";
    for (const auto& [source, ss] : per_source) {
        out << source << "\t" << ss.genes
            << "\t" << ss.segments << "\t" << ss.exclusive_segments
            << "\t" << ss.exons << "\t" << ss.exclusive_exons << "\n";
    }

    logging::info("Per-source statistics written to: " + path);
}

void index_stats::write_per_sample_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open per-sample file: " + path);
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    out << "sample\ttype\tgenes\tsegments\texclusive_segments\texons\texclusive_exons"
        << "\ttranscripts\tisoform_diversity\tdeduplication_ratio"
        << "\texpressed_segments\tmean_expression\n";

    out << std::fixed;
    for (uint32_t sid : sample_ids) {
        const auto& ss = per_sample.at(sid);
        const auto* info = registry.get(sid);
        std::string label = (info && !info->id.empty()) ? info->id : std::to_string(sid);
        std::string type = (info ? info->type : "sample");

        out << label << "\t" << type
            << "\t" << ss.genes
            << "\t" << ss.segments << "\t" << ss.exclusive_segments
            << "\t" << ss.exons << "\t" << ss.exclusive_exons
            << "\t" << ss.transcripts
            << "\t" << std::setprecision(2) << ss.isoform_diversity
            << "\t" << std::setprecision(3) << ss.deduplication_ratio;

        if (ss.expressed_segments > 0) {
            out << "\t" << ss.expressed_segments
                << "\t" << std::setprecision(2) << ss.mean_expression;
        } else {
            out << "\t.\t.";
        }
        out << "\n";
    }

    logging::info("Per-sample statistics written to: " + path);
}

void index_stats::write_biotype_tsv(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open biotype file: " + path);
        return;
    }

    if (genes_by_biotype.empty() && transcripts_by_biotype.empty()) {
        logging::warning("No biotype data to write");
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs in sorted order
    std::vector<uint32_t> sample_ids;
    for (const auto& [sid, _] : per_sample) {
        sample_ids.push_back(sid);
    }
    std::sort(sample_ids.begin(), sample_ids.end());

    // Resolve sample labels
    std::vector<std::string> labels;
    for (uint32_t sid : sample_ids) {
        const auto* info = registry.get(sid);
        labels.push_back((info && !info->id.empty()) ? info->id : std::to_string(sid));
    }

    // Header
    out << "level\tbiotype\ttotal";
    for (const auto& label : labels) {
        out << "\t" << label;
    }
    out << "\n";

    // Gene biotype rows (sorted by count descending)
    {
        std::vector<std::pair<std::string, size_t>> sorted_bt(
            genes_by_biotype.begin(), genes_by_biotype.end());
        std::sort(sorted_bt.begin(), sorted_bt.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        for (const auto& [biotype, count] : sorted_bt) {
            out << "gene\t" << biotype << "\t" << count;
            for (uint32_t sid : sample_ids) {
                auto it = per_sample.find(sid);
                if (it != per_sample.end()) {
                    auto bt_it = it->second.genes_by_biotype.find(biotype);
                    out << "\t" << (bt_it != it->second.genes_by_biotype.end() ? bt_it->second : size_t{0});
                } else {
                    out << "\t0";
                }
            }
            out << "\n";
        }
    }

    // Transcript biotype rows (sorted by count descending)
    {
        std::vector<std::pair<std::string, size_t>> sorted_bt(
            transcripts_by_biotype.begin(), transcripts_by_biotype.end());
        std::sort(sorted_bt.begin(), sorted_bt.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        for (const auto& [biotype, count] : sorted_bt) {
            out << "transcript\t" << biotype << "\t" << count;
            for (uint32_t sid : sample_ids) {
                auto it = per_sample.find(sid);
                if (it != per_sample.end()) {
                    auto bt_it = it->second.transcripts_by_biotype.find(biotype);
                    out << "\t" << (bt_it != it->second.transcripts_by_biotype.end() ? bt_it->second : size_t{0});
                } else {
                    out << "\t0";
                }
            }
            out << "\n";
        }
    }

    logging::info("Biotype statistics written to: " + path);
}