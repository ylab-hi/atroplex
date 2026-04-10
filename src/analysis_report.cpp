/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "analysis_report.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "sample_info.hpp"
#include "utility.hpp"

namespace {

double compute_median(std::vector<size_t>& values) {
    if (values.empty()) return 0;
    auto mid = values.begin() + static_cast<long>(values.size()) / 2;
    std::nth_element(values.begin(), mid, values.end());
    if (values.size() % 2 == 0) {
        auto mid2 = std::max_element(values.begin(), mid);
        return (static_cast<double>(*mid) + static_cast<double>(*mid2)) / 2.0;
    }
    return static_cast<double>(*mid);
}

} // anonymous namespace

void analysis_report::collect(grove_type& grove) {
    logging::info("Collecting analysis report...");

    size_t num_samples = sample_registry::instance().size();
    per_sample.resize(num_samples);

    // Count sample-type entries for conserved classification
    size_t total_samples_for_conserved = 0;
    for (size_t i = 0; i < num_samples; ++i) {
        const auto& info = sample_registry::instance().get(static_cast<uint32_t>(i));
        if (info.type == "sample") total_samples_for_conserved++;
    }

    // ── Per-gene accumulator (temporary, freed at end of chromosome) ─
    struct gene_acc {
        uint32_t gene_idx = 0;
        std::string biotype;
        size_t segment_count = 0;
        sample_bitset sample_bits;
        std::vector<size_t> sample_tx;  // per-sample transcript count
    };
    std::unordered_map<uint32_t, gene_acc> active_genes;

    auto finalize_gene = [&](const gene_acc& acc) {
        transcripts_per_gene.push_back(acc.segment_count);

        for (uint32_t sid : acc.sample_bits) {
            per_sample[sid].genes++;
            if (!acc.biotype.empty()) {
                per_sample[sid].genes_by_biotype[acc.biotype]++;
            }
            if (sid < acc.sample_tx.size()) {
                per_sample[sid].transcripts += acc.sample_tx[sid];
            }
        }
    };

    // Track visited exons (pointer dedup)
    std::unordered_set<const void*> visited_exons;

    // ── Traverse grove ──────────────────────────────────────────────
    auto roots = grove.get_root_nodes();

    for (auto& [seqid, root] : roots) {
        if (!root) continue;
        active_genes.clear();

        logging::info("  " + seqid);

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
                if (seg.absorbed) {
                    absorbed_segments++;
                    continue;
                }

                // ── Segment → per-sample ────────────────────────────
                size_t seg_sample_count = seg.sample_count();
                bool seg_exclusive = (seg_sample_count == 1);
                bool seg_conserved = seg.is_conserved(total_samples_for_conserved);

                exons_per_segment.push_back(static_cast<size_t>(seg.exon_count));

                for (uint32_t sid : seg.sample_idx) {
                    auto& sc = per_sample[sid];
                    sc.segments++;
                    if (seg_exclusive) sc.exclusive_segments++;
                    else if (seg_conserved) sc.conserved_segments++;
                    else sc.shared_segments++;
                    if (seg.exon_count == 1) sc.single_exon_segments++;

                    // Expression
                    if (seg.has_expression(sid)) {
                        sc.expression_sum += seg.get_expression(sid);
                        sc.expressed_segments++;
                    }
                }

                // ── Gene accumulation ───────────────────────────────
                uint32_t gidx = seg.gene_idx;
                auto& acc = active_genes[gidx];
                if (acc.segment_count == 0) {
                    acc.gene_idx = gidx;
                    acc.biotype = seg.gene_biotype();
                    acc.sample_tx.resize(num_samples, 0);
                }
                acc.segment_count++;
                acc.sample_bits.merge(seg.sample_idx);

                size_t tx_count = seg.transcript_ids.size();
                total_transcripts += tx_count;
                for (uint32_t sid : seg.sample_idx) {
                    acc.sample_tx[sid] += tx_count;
                }

                // Transcript biotypes → per-sample
                for (const auto& [tx_id, biotype] : seg.transcript_biotypes) {
                    if (!biotype.empty()) {
                        for (uint32_t sid : seg.sample_idx) {
                            per_sample[sid].transcripts_by_biotype[biotype]++;
                        }
                    }
                }

                // ── Exon chain traversal ────────────────────────────
                size_t edge_id = seg.segment_index;
                auto first_exons = grove.get_neighbors_if(key,
                    [edge_id](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON
                            && e.id == edge_id;
                    });

                if (!first_exons.empty()) {
                    auto* current = first_exons.front();
                    while (current) {
                        if (visited_exons.insert(current).second) {
                            auto& exon = get_exon(current->get_data());
                            size_t exon_sample_count = exon.sample_count();
                            bool exon_exclusive = (exon_sample_count == 1);
                            bool exon_conserved = exon.is_conserved(total_samples_for_conserved);

                            for (uint32_t sid : exon.sample_idx) {
                                auto& sc = per_sample[sid];
                                sc.exons++;
                                if (exon_exclusive) sc.exclusive_exons++;
                                else if (exon_conserved) sc.conserved_exons++;
                                else sc.shared_exons++;
                            }
                        }

                        auto next = grove.get_neighbors_if(current,
                            [edge_id](const edge_metadata& e) {
                                return e.type == edge_metadata::edge_type::EXON_TO_EXON
                                    && e.id == edge_id;
                            });
                        current = next.empty() ? nullptr : next.front();
                    }
                }
            }
            node = node->get_next();
        }

        // Finalize all genes for this chromosome
        for (auto& [gidx, acc] : active_genes) {
            finalize_gene(acc);
        }
        active_genes.clear();
    }

    total_exons = visited_exons.size();
    total_edges = grove.edge_count();

    // Summary log
    size_t total_seg = 0, total_ex = 0, total_genes = 0;
    for (const auto& sc : per_sample) {
        total_seg = std::max(total_seg, sc.segments);
        total_ex = std::max(total_ex, sc.exons);
        total_genes = std::max(total_genes, sc.genes);
    }
    logging::info("Analysis report collected: " +
                  std::to_string(per_sample.size()) + " samples, " +
                  std::to_string(exons_per_segment.size()) + " segments, " +
                  std::to_string(transcripts_per_gene.size()) + " genes");
}

// ── Output ──────────────────────────────────────────────────────────

void analysis_report::write_overview(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) return;

    size_t total_segments = exons_per_segment.size();
    size_t total_genes = transcripts_per_gene.size();

    // Distribution stats
    auto txpg = transcripts_per_gene;
    double mean_tpg = 0, med_tpg = 0;
    size_t max_tpg = 0, single_iso = 0, multi_iso = 0;
    if (!txpg.empty()) {
        size_t sum = std::accumulate(txpg.begin(), txpg.end(), size_t{0});
        mean_tpg = static_cast<double>(sum) / static_cast<double>(txpg.size());
        med_tpg = compute_median(txpg);
        max_tpg = *std::max_element(txpg.begin(), txpg.end());
        single_iso = static_cast<size_t>(std::count(txpg.begin(), txpg.end(), size_t{1}));
        multi_iso = txpg.size() - single_iso;
    }

    auto epsg = exons_per_segment;
    double mean_eps = 0, med_eps = 0;
    size_t max_eps = 0, single_exon = 0;
    if (!epsg.empty()) {
        size_t sum = std::accumulate(epsg.begin(), epsg.end(), size_t{0});
        mean_eps = static_cast<double>(sum) / static_cast<double>(epsg.size());
        med_eps = compute_median(epsg);
        max_eps = *std::max_element(epsg.begin(), epsg.end());
        single_exon = static_cast<size_t>(std::count(epsg.begin(), epsg.end(), size_t{1}));
    }

    double dedup = (total_transcripts > 0)
        ? static_cast<double>(total_segments) / static_cast<double>(total_transcripts)
        : 0;

    out << "metric\tvalue\n" << std::fixed;
    out << "samples\t" << per_sample.size() << "\n";
    out << "genes\t" << total_genes << "\n";
    out << "transcripts\t" << total_transcripts << "\n";
    out << "segments\t" << total_segments << "\n";
    out << "exons\t" << total_exons << "\n";
    out << "graph_edges\t" << total_edges << "\n";
    out << "absorbed_segments\t" << absorbed_segments << "\n";
    out << "single_exon_segments\t" << single_exon << "\n";
    out << "single_isoform_genes\t" << single_iso << "\n";
    out << "multi_isoform_genes\t" << multi_iso << "\n";
    out << "mean_segments_per_gene\t" << std::setprecision(2) << mean_tpg << "\n";
    out << "median_segments_per_gene\t" << med_tpg << "\n";
    out << "max_segments_per_gene\t" << max_tpg << "\n";
    out << "mean_exons_per_segment\t" << mean_eps << "\n";
    out << "median_exons_per_segment\t" << med_eps << "\n";
    out << "max_exons_per_segment\t" << max_eps << "\n";
    out << "deduplication_ratio\t" << std::setprecision(3) << dedup << "\n";

    logging::info("Overview written to: " + path);
}

void analysis_report::write_per_sample(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) return;

    auto& registry = sample_registry::instance();

    out << "sample\ttype\tgenes\tsegments\texclusive_segments\tshared_segments"
        << "\tconserved_segments\texons\texclusive_exons\tshared_exons"
        << "\tconserved_exons\ttranscripts\tsingle_exon_segments"
        << "\texpressed_segments\tmean_expression\tdeduplication_ratio\n";
    out << std::fixed;

    for (size_t sid = 0; sid < per_sample.size(); ++sid) {
        const auto& sc = per_sample[sid];
        if (sc.segments == 0) continue;

        const auto& info = registry.get(static_cast<uint32_t>(sid));
        std::string label = info.id.empty() ? std::to_string(sid) : info.id;

        double dedup = (sc.transcripts > 0)
            ? static_cast<double>(sc.segments) / static_cast<double>(sc.transcripts)
            : 0;
        double mean_expr = (sc.expressed_segments > 0)
            ? sc.expression_sum / static_cast<double>(sc.expressed_segments)
            : 0;

        out << label << "\t" << info.type
            << "\t" << sc.genes
            << "\t" << sc.segments
            << "\t" << sc.exclusive_segments
            << "\t" << sc.shared_segments
            << "\t" << sc.conserved_segments
            << "\t" << sc.exons
            << "\t" << sc.exclusive_exons
            << "\t" << sc.shared_exons
            << "\t" << sc.conserved_exons
            << "\t" << sc.transcripts
            << "\t" << sc.single_exon_segments
            << "\t" << sc.expressed_segments
            << "\t" << std::setprecision(2) << mean_expr
            << "\t" << std::setprecision(3) << dedup
            << "\n";
    }

    logging::info("Per-sample stats written to: " + path);
}