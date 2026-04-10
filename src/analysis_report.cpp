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
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include "sample_info.hpp"
#include "utility.hpp"

namespace gdt = genogrove::data_type;

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

    // Size per-sample counters to the number of registered samples
    size_t num_samples = sample_registry::instance().size();
    per_sample.resize(num_samples);

    // ── Per-gene accumulator ────────────────────────────────────────
    // Lightweight: only counters and a biotype string.
    // Active for overlapping genes on the current chromosome.
    // Finalized when we pass the gene's last segment coordinate.
    struct gene_acc {
        uint32_t gene_idx;
        std::string biotype;
        size_t segment_count = 0;
        size_t transcript_count = 0;       // accumulated from segment transcript_ids
        size_t last_end = 0;               // coordinate end of last segment seen
        sample_bitset sample_bits;         // union of all segments' sample_idx
        // Per-sample transcript count: keyed by sample_id.
        // We need this to compute per-sample transcripts per gene.
        // Bounded by gene size (typically <100 segments), not by sample count.
        // For each segment in the gene, we add transcript_ids.size() for each sample
        // in that segment's sample_idx. But transcripts can be shared across segments
        // (dedup merges metadata onto same segment), so each tx_id is on exactly one
        // segment — seg.transcript_ids.size() * sample_count is the correct accumulation.
        //
        // Actually: unique transcripts per sample = sum over segments of
        // (segment.transcript_ids.size()) for segments where sample is present.
        // This works because each transcript_id belongs to exactly one segment.
        std::vector<size_t> sample_tx;     // indexed by sample_id, sized on first use
    };

    // Active genes on current chromosome, keyed by gene_idx
    std::unordered_map<uint32_t, gene_acc> active_genes;

    // Finalize a gene accumulator: update global + per-sample stats
    auto finalize_gene = [&](gene_acc& acc) {
        total_genes++;
        total_transcripts += acc.transcript_count;

        // Biotype
        if (!acc.biotype.empty()) {
            genes_by_biotype[acc.biotype]++;
        }

        // Isoform distribution
        if (acc.segment_count == 1) single_isoform_genes++;
        else if (acc.segment_count > 1) multi_isoform_genes++;

        // Track for median/max computation (deferred)
        // We'll use a vector to collect these — bounded by gene count, not sample count

        // Per-sample: increment gene count for each sample in this gene
        for (uint32_t sid : acc.sample_bits) {
            per_sample[sid].genes++;
            if (!acc.biotype.empty()) {
                per_sample[sid].genes_by_biotype[acc.biotype]++;
            }
        }

        // Per-sample transcripts
        if (!acc.sample_tx.empty()) {
            for (uint32_t sid : acc.sample_bits) {
                if (sid < acc.sample_tx.size()) {
                    per_sample[sid].transcripts += acc.sample_tx[sid];
                }
            }
        }
    };

    // Distribution vectors — bounded by gene/segment count, not samples
    std::vector<size_t> transcripts_per_gene_vec;
    std::vector<size_t> exons_per_segment_vec;

    // Track visited exons globally (pointer dedup)
    std::unordered_set<const void*> visited_exons;

    // ── Traverse grove ──────────────────────────────────────────────
    auto roots = grove.get_root_nodes();
    total_chromosomes = roots.size();

    for (auto& [seqid, root] : roots) {
        if (!root) continue;
        active_genes.clear();

        logging::info("  Processing chromosome: " + seqid);

        // Walk to first leaf
        auto* node = root;
        while (!node->get_is_leaf()) {
            auto& children = node->get_children();
            if (children.empty()) break;
            node = children[0];
        }

        // Walk all leaves
        while (node) {
            for (auto* key : node->get_keys()) {
                auto& feature = key->get_data();
                if (!is_segment(feature)) continue;

                auto& seg = get_segment(feature);
                if (seg.absorbed) {
                    absorbed_segments++;
                    continue;
                }

                auto& coord = key->get_value();

                // ── Segment counting ────────────────────────────
                total_segments++;
                per_chromosome[seqid].segments++;

                if (seg.exon_count == 1) single_exon_segments++;
                exons_per_segment_vec.push_back(static_cast<size_t>(seg.exon_count));

                // Per-sample segment counting
                for (uint32_t sid : seg.sample_idx) {
                    per_sample[sid].segments++;
                    if (seg.exon_count == 1) per_sample[sid].single_exon_segments++;
                }

                // ── Gene accumulation ───────────────────────────
                uint32_t gidx = seg.gene_idx;
                auto& acc = active_genes[gidx];
                if (acc.segment_count == 0) {
                    // First segment for this gene
                    acc.gene_idx = gidx;
                    acc.biotype = seg.gene_biotype();
                    acc.sample_tx.resize(num_samples, 0);
                }
                acc.segment_count++;
                acc.last_end = std::max(acc.last_end, coord.get_end());
                acc.sample_bits.merge(seg.sample_idx);

                // Accumulate per-sample transcript count
                size_t tx_count = seg.transcript_ids.size();
                for (uint32_t sid : seg.sample_idx) {
                    acc.sample_tx[sid] += tx_count;
                }
                acc.transcript_count += tx_count;

                // Transcript biotypes
                for (const auto& [tx_id, biotype] : seg.transcript_biotypes) {
                    if (!biotype.empty()) {
                        transcripts_by_biotype[biotype]++;
                        // Per-sample transcript biotype
                        for (uint32_t sid : seg.sample_idx) {
                            per_sample[sid].transcripts_by_biotype[biotype]++;
                        }
                    }
                }

                // ── Exon chain traversal ────────────────────────
                size_t edge_id = seg.segment_index;
                auto first_exons = grove.get_neighbors_if(key,
                    [edge_id](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON
                            && e.id == edge_id;
                    });

                if (!first_exons.empty()) {
                    auto* current = first_exons.front();
                    while (current) {
                        // Count exon only on first visit (pointer dedup)
                        if (visited_exons.insert(current).second) {
                            total_exons++;
                            per_chromosome[seqid].exons++;

                            // Per-sample exon counting
                            auto& exon = get_exon(current->get_data());
                            for (uint32_t sid : exon.sample_idx) {
                                per_sample[sid].exons++;
                            }
                        }

                        // Follow chain
                        auto next = grove.get_neighbors_if(current,
                            [edge_id](const edge_metadata& e) {
                                return e.type == edge_metadata::edge_type::EXON_TO_EXON
                                    && e.id == edge_id;
                            });
                        current = next.empty() ? nullptr : next.front();
                    }
                }

                // ── Finalize genes that we've passed ────────────
                // Check if any active gene's last_end is before this segment's start
                // (meaning no more segments for that gene can appear)
                auto it = active_genes.begin();
                while (it != active_genes.end()) {
                    if (it->first != gidx && it->second.last_end < coord.get_start()) {
                        transcripts_per_gene_vec.push_back(it->second.transcript_count);
                        per_chromosome[seqid].genes++;
                        finalize_gene(it->second);
                        it = active_genes.erase(it);
                    } else {
                        ++it;
                    }
                }
            }
            node = node->get_next();
        }

        // Finalize remaining genes for this chromosome
        for (auto& [gidx, acc] : active_genes) {
            transcripts_per_gene_vec.push_back(acc.transcript_count);
            per_chromosome[seqid].genes++;
            finalize_gene(acc);
        }
        active_genes.clear();
    }

    // ── Compute distributions ───────────────────────────────────────
    if (!transcripts_per_gene_vec.empty()) {
        size_t sum = std::accumulate(transcripts_per_gene_vec.begin(),
                                      transcripts_per_gene_vec.end(), size_t{0});
        mean_transcripts_per_gene = static_cast<double>(sum) /
                                     static_cast<double>(transcripts_per_gene_vec.size());
        median_transcripts_per_gene = compute_median(transcripts_per_gene_vec);
        max_transcripts_per_gene = *std::max_element(
            transcripts_per_gene_vec.begin(), transcripts_per_gene_vec.end());
    }

    if (!exons_per_segment_vec.empty()) {
        size_t sum = std::accumulate(exons_per_segment_vec.begin(),
                                      exons_per_segment_vec.end(), size_t{0});
        mean_exons_per_segment = static_cast<double>(sum) /
                                  static_cast<double>(exons_per_segment_vec.size());
        median_exons_per_segment = compute_median(exons_per_segment_vec);
        max_exons_per_segment = *std::max_element(
            exons_per_segment_vec.begin(), exons_per_segment_vec.end());
    }

    // Dedup ratio
    if (total_transcripts > 0) {
        deduplication_ratio = static_cast<double>(total_segments) /
                               static_cast<double>(total_transcripts);
    }

    // Per-sample dedup ratio
    for (size_t sid = 0; sid < num_samples; ++sid) {
        auto& sc = per_sample[sid];
        if (sc.transcripts > 0) {
            sc.deduplication_ratio = static_cast<double>(sc.segments) /
                                     static_cast<double>(sc.transcripts);
        }
    }

    // Edge count
    total_edges = grove.edge_count();

    logging::info("Analysis report collected: " +
                  std::to_string(total_genes) + " genes, " +
                  std::to_string(total_segments) + " segments, " +
                  std::to_string(total_exons) + " exons");
}

// ── Output methods ──────────────────────────────────────────────────

void analysis_report::write_overview(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) return;

    out << "metric\tvalue\n";
    out << "chromosomes\t" << total_chromosomes << "\n";
    out << "genes\t" << total_genes << "\n";
    out << "transcripts\t" << total_transcripts << "\n";
    out << "segments\t" << total_segments << "\n";
    out << "exons\t" << total_exons << "\n";
    out << "graph_edges\t" << total_edges << "\n";
    out << "absorbed_segments\t" << absorbed_segments << "\n";
    out << "single_exon_segments\t" << single_exon_segments << "\n";
    out << "single_isoform_genes\t" << single_isoform_genes << "\n";
    out << "multi_isoform_genes\t" << multi_isoform_genes << "\n";
    out << "mean_transcripts_per_gene\t" << std::setprecision(2) << std::fixed
        << mean_transcripts_per_gene << "\n";
    out << "median_transcripts_per_gene\t" << median_transcripts_per_gene << "\n";
    out << "max_transcripts_per_gene\t" << max_transcripts_per_gene << "\n";
    out << "mean_exons_per_segment\t" << mean_exons_per_segment << "\n";
    out << "median_exons_per_segment\t" << median_exons_per_segment << "\n";
    out << "max_exons_per_segment\t" << max_exons_per_segment << "\n";
    out << "deduplication_ratio\t" << std::setprecision(3) << deduplication_ratio << "\n";

    logging::info("Overview written to: " + path);
}

void analysis_report::write_per_chromosome(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) return;

    out << "chromosome\tgenes\tsegments\texons\n";
    for (const auto& [chr, stats] : per_chromosome) {
        out << chr << "\t" << stats.genes << "\t" << stats.segments
            << "\t" << stats.exons << "\n";
    }

    logging::info("Per-chromosome stats written to: " + path);
}

void analysis_report::write_biotypes(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) return;

    out << "level\tbiotype\tcount\n";
    for (const auto& [bt, count] : genes_by_biotype) {
        out << "gene\t" << bt << "\t" << count << "\n";
    }
    for (const auto& [bt, count] : transcripts_by_biotype) {
        out << "transcript\t" << bt << "\t" << count << "\n";
    }

    logging::info("Biotype stats written to: " + path);
}