/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "build_summary.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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

} // anonymous namespace

build_summary build_summary::collect(
    grove_type& grove,
    const chromosome_segment_caches& segment_caches,
    const chromosome_exon_caches& exon_caches,
    size_t segment_count,
    const build_counters& counters
) {
    build_summary stats;
    stats.counters = counters;

    // Annotation sources and entry count from registry
    auto& registry = sample_registry::instance();
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto& info = registry.get(static_cast<uint32_t>(i));
        if (!info.id.empty()) {
            stats.annotation_sources.push_back(info.id);
        }
        if (info.type != "replicate") {
            stats.total_entries++;
        }
    }

    stats.total_chromosomes = segment_caches.size();
    // Pre-sweep running total minus tombstones physically removed
    stats.total_segments = (segment_count >= counters.absorbed_segments)
        ? segment_count - counters.absorbed_segments
        : segment_count;
    stats.tree_order = grove.get_order();

    // B+ tree depth per chromosome
    for (const auto& [index_name, root_node] : grove.get_root_nodes()) {
        int depth = 0;
        auto* current = root_node;
        while (current && !current->get_is_leaf()) {
            current = current->get_children().empty() ? nullptr : current->get_children()[0];
            depth++;
        }
        if (current) depth++;
        stats.tree_depth_per_chromosome[index_name] = depth;
    }

    // Gene / transcript aggregation from segment caches (post-sweep: no tombstones)
    struct gene_info {
        std::unordered_set<uint32_t> transcript_ids;
        std::string biotype;
    };
    std::unordered_map<std::string, gene_info> genes;
    std::unordered_map<std::string, std::unordered_set<std::string>> chr_gene_ids;
    std::vector<size_t> exons_per_segment;

    for (const auto& [seqid, seg_cache] : segment_caches) {
        stats.per_chromosome[seqid].segments = seg_cache.size();

        for (const auto& [key, seg_ptr] : seg_cache) {
            auto& feature = seg_ptr->get_data();
            if (!is_segment(feature)) continue;

            auto& seg = get_segment(feature);
            auto& gi = genes[seg.gene_id()];
            gi.biotype = seg.gene_biotype();
            for (const auto& tx : seg.transcript_ids) {
                gi.transcript_ids.insert(tx);
            }
            chr_gene_ids[seqid].insert(seg.gene_id());

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

    // Exon counts
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

    if (!exons_per_segment.empty()) {
        size_t total_exons_sum = std::accumulate(
            exons_per_segment.begin(), exons_per_segment.end(), size_t{0});
        stats.mean_exons_per_segment =
            static_cast<double>(total_exons_sum) / static_cast<double>(exons_per_segment.size());
        stats.median_exons_per_segment = compute_median(exons_per_segment);
        stats.max_exons_per_segment = *std::max_element(
            exons_per_segment.begin(), exons_per_segment.end());
    }

    for (const auto& [seqid, gene_set] : chr_gene_ids) {
        stats.per_chromosome[seqid].genes = gene_set.size();
    }

    stats.total_edges = grove.edge_count();

    return stats;
}

void build_summary::write_summary(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open summary file: " + path);
        return;
    }

    out << std::fixed << std::setprecision(2);

    out << "Index Summary\n";
    out << "=============\n\n";

    // Build time
    if (build_time_seconds > 0) {
        if (build_time_seconds >= 60) {
            int minutes = static_cast<int>(build_time_seconds) / 60;
            double seconds = build_time_seconds - minutes * 60;
            out << "Build time:         " << minutes << "m "
                << std::setprecision(1) << seconds << "s\n";
        } else {
            out << "Build time:         " << std::setprecision(1)
                << build_time_seconds << "s\n";
        }
        out << std::setprecision(2) << "\n";
    }

    // Inputs
    if (!annotation_sources.empty()) {
        auto& reg = sample_registry::instance();

        size_t n_annotations = 0, n_replicates = 0;
        for (size_t i = 0; i < reg.size(); ++i) {
            const auto& info = reg.get(static_cast<uint32_t>(i));
            if (info.type == "annotation") n_annotations++;
            else if (info.type == "replicate") n_replicates++;
        }

        out << "Inputs: " << (reg.size() - n_replicates) << " entries ("
            << total_entries << " samples, "
            << n_annotations << " annotations";
        if (n_replicates > 0) {
            out << ", " << n_replicates << " replicates merged";
        }
        out << ")\n";
        for (size_t i = 0; i < reg.size(); ++i) {
            const auto& info = reg.get(static_cast<uint32_t>(i));
            if (!info.id.empty() && info.type != "replicate") {
                out << "  " << info.id << "  [" << info.type << "]\n";
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
    out << "Unique exons:       " << total_exons << "\n";
    out << "Graph edges:        " << total_edges << "\n";
    out << "\n";

    // Processing summary — what happened to input transcripts during build
    out << "Transcript processing:\n";
    if (counters.input_transcripts > 0) {
        out << "  Input:            " << counters.input_transcripts << "\n";
    }
    out << "  Merged:           " << counters.merged_transcripts
        << "  (folded into existing segments via FSM/terminal/ISM absorption)\n";
    out << "  Absorbed:         " << counters.absorbed_segments
        << "  (segments tombstoned by reverse absorption and removed)\n";
    out << "  Discarded:        " << counters.discarded_transcripts
        << "  (--min-expression, mono-exon, fragment drops)\n";
    if (counters.replicates_merged > 0) {
        out << "  Replicates merged:" << counters.replicates_merged << "\n";
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
            double pct = total_genes > 0
                ? 100.0 * static_cast<double>(count) / static_cast<double>(total_genes) : 0;
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
            double pct = total_transcripts > 0
                ? 100.0 * static_cast<double>(count) / static_cast<double>(total_transcripts) : 0;
            out << "  " << std::left << std::setw(28) << biotype
                << std::right << std::setw(8) << count
                << "  (" << std::setprecision(1) << pct << "%)\n";
        }
        out << std::setprecision(2) << "\n";
    }

    // B+ tree structure
    if (tree_order > 0) {
        out << "B+ tree:\n";
        out << "  Order:            " << tree_order
            << "  (max " << (tree_order - 1) << " keys/node)\n";
        if (!tree_depth_per_chromosome.empty()) {
            int max_depth = 0;
            for (const auto& [chr, d] : tree_depth_per_chromosome) {
                if (d > max_depth) max_depth = d;
            }
            out << "  Max tree depth:   " << max_depth << " levels\n";
        }
        out << "\n";
    }

    // Per chromosome
    if (!per_chromosome.empty()) {
        out << "Per chromosome:\n";
        out << "  " << std::left << std::setw(10) << "Chrom"
            << std::right << std::setw(8) << "Genes"
            << std::setw(10) << "Segments"
            << std::setw(10) << "Exons"
            << std::setw(8) << "Depth" << "\n";
        for (const auto& [seqid, cs] : per_chromosome) {
            out << "  " << std::left << std::setw(10) << seqid
                << std::right << std::setw(8) << cs.genes
                << std::setw(10) << cs.segments
                << std::setw(10) << cs.exons;
            auto depth_it = tree_depth_per_chromosome.find(seqid);
            if (depth_it != tree_depth_per_chromosome.end()) {
                out << std::setw(8) << depth_it->second;
            }
            out << "\n";
        }
        out << "\n";
    }

    out << "For detailed per-sample analysis, run: atroplex analyze\n";

    logging::info("Index summary written to: " + path);
}