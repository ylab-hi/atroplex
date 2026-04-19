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
#include "splicing_catalog.hpp"
#include "utility.hpp"

namespace {

/// Render a sample's declared expression_type as the short suffix used
/// in per-sample output column headers (e.g., `SAMPLE.counts`).
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
    auto mid = values.begin() + static_cast<long>(values.size()) / 2;
    std::nth_element(values.begin(), mid, values.end());
    if (values.size() % 2 == 0) {
        auto mid2 = std::max_element(values.begin(), mid);
        return (static_cast<double>(*mid) + static_cast<double>(*mid2)) / 2.0;
    }
    return static_cast<double>(*mid);
}

} // anonymous namespace

void analysis_report::collect(grove_type& grove,
                              quant_sidecar::Reader* qtx_reader) {
    logging::info("Collecting analysis report...");
    if (qtx_reader) {
        logging::info("Quantitative columns enabled (reading from .qtx sidecar)");
    }

    size_t num_samples = sample_registry::instance().size();
    per_sample.resize(num_samples);

    // Count sample-type entries for conserved classification
    size_t total_samples_for_conserved = 0;
    for (size_t i = 0; i < num_samples; ++i) {
        const auto& info = sample_registry::instance().get(static_cast<uint32_t>(i));
        if (info.type == "sample") total_samples_for_conserved++;
    }

    // ── Per-gene accumulator (temporary, freed at end of chromosome) ─
    //
    // Grows only with features-in-this-gene, never with sample count. At
    // chromosome-boundary finalization we compute diversity metrics from
    // `segment_tx_counts` and `exon_seg_counts`, then discard the whole
    // `active_genes` map.
    struct pending_hub {
        key_ptr exon = nullptr;
        std::vector<key_ptr> targets;  // deduplicated downstream exons
        size_t chain_pos = 0;          // 1-based position in first-visit segment
        size_t chain_total = 0;        // total exons in that segment
    };
    struct gene_acc {
        uint32_t gene_idx = 0;
        std::string biotype;
        size_t segment_count = 0;
        sample_bitset sample_bits;
        // Per-sample "transcripts touched": sum over segments the sample
        // contributes to of that segment's full transcript_ids count. This
        // over-approximates for shared segments — the merged transcript_ids
        // of a segment shared by N samples get credited to every sample,
        // even though each sample only really contributed a subset. Strict
        // per-sample attribution would require a tx→sample map, which is
        // the sample×feature structure the streaming rewrite exists to
        // avoid. The TSV column is labelled `transcripts_touched` so the
        // over-approximation is visible to consumers.
        std::vector<size_t> sample_tx;
        std::vector<size_t> segment_tx_counts;    // Phase 8.4: for effective isoforms
        std::unordered_map<key_ptr, size_t> exon_seg_counts;  // Phase 8.4: segments-per-exon within gene
        std::vector<key_ptr> first_visit_exons;   // Phase 8.4: exons whose global first visit occurred in this gene
        std::vector<pending_hub> pending_hubs;    // Phase 8.3: hubs detected in this gene
        size_t max_hub_targets_in_gene = 0;       // widest fan-out among this gene's hubs (for lazy branch_counts grow)
        std::vector<segment_chain_entry> segments_in_gene;  // Phase 8.6: per-gene chains for splicing_catalog detectors
        uint16_t sources_seen_in_gene = 0;        // 8.7a: OR of all segment.sources bitfields in the gene
    };
    std::unordered_map<uint32_t, gene_acc> active_genes;

    // Per-exon × per-sample expression accumulator, cleared per-chromosome
    // alongside `active_genes`. The sidecar stores per-segment values; we
    // attribute them to every exon in the segment's chain so hub /
    // branch / conserved-exon row emission can quote a per-exon expression
    // value (= total sample-level support flowing through this exon).
    // Empty when qtx_reader is null — no allocations, no lookups, no
    // accumulator overhead in the legacy/no-sidecar path.
    using exon_expr_map_t = std::unordered_map<key_ptr,
                                               std::unordered_map<uint32_t, float>>;
    exon_expr_map_t exon_expr_sum;

    // Phase 8.3: single set of per-sample tally buffers, allocated once and
    // reused across every hub in every gene. std::fill-reset before each hub.
    // Peak = O(num_samples), never O(num_hubs × num_samples).
    std::vector<size_t> hub_branches(num_samples, 0);
    std::vector<size_t> hub_shared(num_samples, 0);
    std::vector<size_t> hub_unique(num_samples, 0);

    // Phase 8.6b: PSI + entropy buffers. hub_psi_den is recomputed once per
    // gene (same denominator for every hub in the gene); hub_psi_num and
    // branch_counts are reset per-hub. branch_counts is a flat
    // num_samples × max_targets_in_gene matrix grown lazily via resize() —
    // over the run it reaches num_samples × global_max_targets and then
    // stops reallocating. Peak stays at O(num_samples × global_max_targets),
    // independent of num_hubs.
    std::vector<size_t> hub_psi_num(num_samples, 0);
    std::vector<size_t> hub_psi_den(num_samples, 0);
    std::vector<size_t> branch_counts;  // grown per-gene in finalize_gene

    auto finalize_gene = [&](const gene_acc& acc, const std::string& seqid_for_hub) {
        transcripts_per_gene.push_back(acc.segment_count);

        // 8.7a: global gene biotype count for biotype.tsv `total` column.
        if (!acc.biotype.empty()) {
            genes_by_biotype[acc.biotype]++;
        }

        for (uint32_t sid : acc.sample_bits) {
            per_sample[sid].genes++;
            if (!acc.biotype.empty()) {
                per_sample[sid].genes_by_biotype[acc.biotype]++;
            }
            if (sid < acc.sample_tx.size()) {
                per_sample[sid].transcripts += acc.sample_tx[sid];
            }
        }

        // 8.7a: bump per_source.genes for every source that contributed to
        // any segment in this gene. sources_seen_in_gene is the OR of all
        // segment.sources bitfields collected during the segment loop.
        source_registry::instance().for_each(acc.sources_seen_in_gene,
            [&](const std::string& src) {
                per_source[src].genes++;
            });

        // ── Phase 8.4: diversity metrics for multi-segment genes ────────
        if (acc.segment_count >= 2) {
            // Effective isoforms: 2^H over the segment→transcript-count distribution
            double total_tx = 0;
            for (size_t c : acc.segment_tx_counts) total_tx += static_cast<double>(c);
            double seg_H = 0;
            if (total_tx > 0) {
                for (size_t c : acc.segment_tx_counts) {
                    if (c == 0) continue;
                    double p = static_cast<double>(c) / total_tx;
                    seg_H -= p * std::log2(p);
                }
            }
            double effective_isoforms = std::pow(2.0, seg_H);

            // Gene exon entropy: H over exon usage fractions within the gene
            double total_usage = 0;
            for (const auto& [e, cnt] : acc.exon_seg_counts) {
                total_usage += static_cast<double>(cnt);
            }
            double exon_H = 0;
            if (total_usage > 0) {
                for (const auto& [e, cnt] : acc.exon_seg_counts) {
                    if (cnt == 0) continue;
                    double p = static_cast<double>(cnt) / total_usage;
                    exon_H -= p * std::log2(p);
                }
            }

            for (uint32_t sid : acc.sample_bits) {
                per_sample[sid].effective_isoforms_sum += effective_isoforms;
                per_sample[sid].gene_exon_entropy_sum += exon_H;
                per_sample[sid].multi_segment_genes++;
            }
        }

        // ── Phase 8.4: constitutive / alternative classification ────────
        // An exon is constitutive if it participates in every segment of its
        // gene. We classify each exon exactly once (on its first-visit gene)
        // so no per-sample double counting.
        for (key_ptr exon_key : acc.first_visit_exons) {
            auto it = acc.exon_seg_counts.find(exon_key);
            if (it == acc.exon_seg_counts.end()) continue;
            bool constitutive = (it->second >= acc.segment_count);

            auto& exon = get_exon(exon_key->get_data());
            for (uint32_t sid : exon.sample_idx) {
                if (constitutive) per_sample[sid].constitutive_exons++;
                else               per_sample[sid].alternative_exons++;
            }
        }

        // ── Phase 8.3: splicing hub rows ────────────────────────────────
        // Write one row per pending hub in this gene, then drop the gene's
        // hub data. The hub_branches/hub_shared/hub_unique buffers are
        // outer-scope vectors reused across all hubs — sized O(num_samples),
        // reset with std::fill before each hub so peak is never multiplied.
        if (hub_stream && hub_stream->is_open() && !acc.pending_hubs.empty()) {
            auto& hub_out = *hub_stream;

            // ── Pass 1 (Phase 8.6b): PSI denominator, shared across all
            // hubs in this gene. hub_psi_den[S] = number of segments in
            // this gene that sample S contributes to — sample S's "path
            // budget" in gene G. Same value is reused for every hub.
            std::fill(hub_psi_den.begin(), hub_psi_den.end(), 0);
            for (const auto& seg_entry : acc.segments_in_gene) {
                if (!seg_entry.segment) continue;
                auto& seg = get_segment(seg_entry.segment->get_data());
                for (uint32_t sid : seg.sample_idx) {
                    hub_psi_den[sid]++;
                }
            }

            // Grow branch_counts lazily to fit the widest hub in this gene.
            // Over the run this stabilizes at num_samples × global_max_targets
            // and stops reallocating.
            size_t needed = num_samples * acc.max_hub_targets_in_gene;
            if (branch_counts.size() < needed) branch_counts.resize(needed);

            for (const auto& hub : acc.pending_hubs) {
                auto& hub_exon = get_exon(hub.exon->get_data());
                const size_t K = hub.targets.size();

                // Reset per-sample tallies for this hub. branch_counts is
                // indexed with THIS hub's K stride, so we only need to zero
                // [0, num_samples × K) — higher cells will be rewritten by
                // hubs with larger K in the same gene (if any).
                std::fill(hub_branches.begin(), hub_branches.end(), 0);
                std::fill(hub_shared.begin(), hub_shared.end(), 0);
                std::fill(hub_unique.begin(), hub_unique.end(), 0);
                std::fill(hub_psi_num.begin(), hub_psi_num.end(), 0);
                std::fill(branch_counts.begin(),
                          branch_counts.begin() + static_cast<long>(num_samples * K),
                          0);

                // Populate branches/shared/unique from downstream target bitsets
                for (key_ptr tgt_key : hub.targets) {
                    auto& tgt = get_exon(tgt_key->get_data());
                    size_t tgt_sample_count = tgt.sample_count();
                    bool tgt_unique = (tgt_sample_count == 1);
                    for (uint32_t sid : tgt.sample_idx) {
                        hub_branches[sid]++;
                        if (tgt_unique) hub_unique[sid]++;
                        else            hub_shared[sid]++;
                    }
                }

                // ── Pass 2 (Phase 8.6b): PSI numerator + per-target branch
                // distribution. Walk every segment in this gene, find the
                // hub's position in its chain, and fan out the segment's
                // sample_idx into both counters.
                for (const auto& seg_entry : acc.segments_in_gene) {
                    if (!seg_entry.segment) continue;
                    const auto& chain = seg_entry.exon_chain;

                    size_t pos = SIZE_MAX;
                    for (size_t i = 0; i < chain.size(); ++i) {
                        if (chain[i] == hub.exon) { pos = i; break; }
                    }
                    if (pos == SIZE_MAX) continue;  // segment doesn't touch hub

                    auto& seg = get_segment(seg_entry.segment->get_data());

                    // PSI numerator
                    for (uint32_t sid : seg.sample_idx) hub_psi_num[sid]++;

                    // Entropy: find the downstream target this segment picks
                    if (pos + 1 < chain.size()) {
                        key_ptr next_exon = chain[pos + 1];
                        size_t t = SIZE_MAX;
                        for (size_t ti = 0; ti < K; ++ti) {
                            if (hub.targets[ti] == next_exon) { t = ti; break; }
                        }
                        if (t != SIZE_MAX) {
                            for (uint32_t sid : seg.sample_idx) {
                                branch_counts[sid * K + t]++;
                            }
                        }
                    }
                }

                // Global row prefix — gene info from gene_acc, not exon
                auto& gene_info = gene_registry::instance().resolve(acc.gene_idx);
                hub_out << gene_info.gene_name << "\t"
                        << gene_info.gene_id << "\t"
                        << hub_exon.id << "\t"
                        << format_coordinate(seqid_for_hub, hub.exon->get_value()) << "\t"
                        << hub.chain_pos << "\t"
                        << hub.chain_total << "\t"
                        << hub.targets.size() << "\t"
                        << hub_exon.transcript_ids.size();

                // Per-sample columns: .branches, .shared, .unique, .psi,
                // .entropy, optional .{expr_type}. Gated on "sample is
                // at the hub" — being at the hub implies being in the
                // gene, so hub_psi_den[sid] > 0 when
                // hub_exon.sample_idx.test(sid).
                auto hub_expr_it = (qtx_reader && hub_emit_expression)
                    ? exon_expr_sum.find(hub.exon)
                    : exon_expr_sum.end();
                for (size_t i = 0; i < hub_stream_sample_ids.size(); ++i) {
                    uint32_t sid = hub_stream_sample_ids[i];
                    if (!hub_exon.sample_idx.test(sid)) {
                        hub_out << "\t.\t.\t.\t.\t.";
                        if (hub_emit_expression && hub_stream_is_sample[i]) {
                            hub_out << "\t.";
                        }
                    } else {
                        // PSI: hub_psi_den[sid] > 0 is guaranteed here
                        double psi = (hub_psi_den[sid] > 0)
                            ? static_cast<double>(hub_psi_num[sid])
                              / static_cast<double>(hub_psi_den[sid])
                            : 0.0;

                        // Entropy over this sample's branch distribution
                        size_t total = 0;
                        for (size_t ti = 0; ti < K; ++ti) {
                            total += branch_counts[sid * K + ti];
                        }
                        double entropy = 0.0;
                        if (total > 0) {
                            for (size_t ti = 0; ti < K; ++ti) {
                                size_t c = branch_counts[sid * K + ti];
                                if (c == 0) continue;
                                double p = static_cast<double>(c)
                                         / static_cast<double>(total);
                                entropy -= p * std::log2(p);
                            }
                        }

                        hub_out << "\t" << hub_branches[sid]
                                << "\t" << hub_shared[sid]
                                << "\t" << hub_unique[sid];
                        {
                            auto prev_prec = hub_out.precision(4);
                            hub_out << "\t" << psi << "\t" << entropy;
                            hub_out.precision(prev_prec);
                        }
                        if (hub_emit_expression && hub_stream_is_sample[i]) {
                            if (hub_expr_it != exon_expr_sum.end()) {
                                auto v_it = hub_expr_it->second.find(sid);
                                if (v_it != hub_expr_it->second.end()) {
                                    hub_out << "\t" << v_it->second;
                                    continue;
                                }
                            }
                            hub_out << "\t.";
                        }
                    }
                }
                hub_out << "\n";

                // Branch detail rows: one per (hub, target)
                if (branch_stream && branch_stream->is_open()) {
                    auto& br_out = *branch_stream;
                    for (key_ptr tgt_key : hub.targets) {
                        auto& tgt = get_exon(tgt_key->get_data());
                        br_out << gene_info.gene_name << "\t"
                               << gene_info.gene_id << "\t"
                               << hub_exon.id << "\t"
                               << format_coordinate(seqid_for_hub, hub.exon->get_value()) << "\t"
                               << tgt.id << "\t"
                               << format_coordinate(seqid_for_hub, tgt_key->get_value());
                        auto tgt_expr_it = (qtx_reader && hub_emit_expression)
                            ? exon_expr_sum.find(tgt_key)
                            : exon_expr_sum.end();
                        for (size_t i = 0; i < hub_stream_sample_ids.size(); ++i) {
                            uint32_t sid = hub_stream_sample_ids[i];
                            bool present = tgt.sample_idx.test(sid);
                            br_out << "\t" << (present ? "1" : ".");
                            if (hub_emit_expression && hub_stream_is_sample[i]) {
                                if (present && tgt_expr_it != exon_expr_sum.end()) {
                                    auto v_it = tgt_expr_it->second.find(sid);
                                    if (v_it != tgt_expr_it->second.end()) {
                                        br_out << "\t" << v_it->second;
                                        continue;
                                    }
                                }
                                br_out << "\t.";
                            }
                        }
                        br_out << "\n";
                    }
                }
            }
        }

        // ── Phase 8.6: splicing event rows ─────────────────────────────
        // detect_gene_events runs cassette / alt-splice / intron
        // retention / alt-terminal / mutually-exclusive detectors over
        // the segments captured in this gene. Events are written
        // inline and their per-gene vectors are dropped with `acc` at
        // the next chromosome boundary.
        if (splicing_events_stream && splicing_events_stream->is_open()
            && !acc.segments_in_gene.empty()) {
            // Resolve gene_id / gene_name from the first live segment.
            std::string gene_id;
            std::string gene_name;
            for (const auto& entry : acc.segments_in_gene) {
                if (!entry.segment) continue;
                auto& s = get_segment(entry.segment->get_data());
                gene_id = s.gene_id();
                gene_name = s.gene_name();
                break;
            }

            auto events = splicing_catalog::detect_gene_events(
                gene_id, seqid_for_hub, acc.segments_in_gene, grove);

            auto& ev_out = *splicing_events_stream;
            for (const auto& event : events) {
                ev_out << event.gene_id
                       << "\t" << event.gene_name
                       << "\t" << splicing_catalog::event_type_str(event.type)
                       << "\t" << event.chromosome
                       << "\t" << splicing_catalog::format_exon(event.upstream_exon)
                       << "\t" << splicing_catalog::format_exon(event.downstream_exon)
                       << "\t";
                for (size_t i = 0; i < event.affected_exons.size(); ++i) {
                    if (i > 0) ev_out << ",";
                    ev_out << splicing_catalog::format_exon(event.affected_exons[i]);
                }
                for (uint32_t sid : splicing_events_sample_ids) {
                    auto psi_it = event.sample_psi.find(sid);
                    auto inc_it = event.sample_included.find(sid);
                    auto tot_it = event.sample_total.find(sid);
                    if (psi_it != event.sample_psi.end()) {
                        ev_out << "\t" << psi_it->second;
                    } else {
                        ev_out << "\t.";
                    }
                    if (inc_it != event.sample_included.end()) {
                        ev_out << "\t" << inc_it->second;
                    } else {
                        ev_out << "\t.";
                    }
                    if (tot_it != event.sample_total.end()) {
                        ev_out << "\t" << tot_it->second;
                    } else {
                        ev_out << "\t.";
                    }
                }
                ev_out << "\n";
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
        exon_expr_sum.clear();

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
                // Skip tombstones and zero-attribution segments (features
                // that lost all sample bits after replicate merging).
                if (seg.absorbed) continue;
                if (seg.sample_count() == 0) continue;

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
                }

                // ── Quantification (sidecar read) ───────────────────
                // When a .qtx reader is available, look up this segment's
                // per-sample expression records once and:
                //   1. accumulate per-sample stats (expression_sum,
                //      expressed_segments) for the sample_stats output
                //   2. attribute the value to every exon in this segment's
                //      chain via exon_expr_sum[exon_key][sample_id] so
                //      hub / branch / conserved-exon emission can quote a
                //      per-exon value derived from the segments containing
                //      the exon
                // Records are pre-sorted by sample_id within a block. We
                // walk the segment's exon chain via the same SEGMENT_TO_EXON
                // / EXON_TO_EXON edges the main exon walk uses below, but
                // do it inline here so the per-segment lookup happens once.
                std::vector<quant_sidecar::Reader::ValueRecord> seg_expr_records;
                if (qtx_reader) {
                    seg_expr_records = qtx_reader->lookup(seg.segment_index);
                    // expression_sum adds every record (sum-aggregation
                    // semantic when multiple same-sample transcripts merged
                    // into this segment). expressed_segments counts distinct
                    // samples only — records within a segment block are
                    // sorted by sample_id, so duplicates are adjacent and a
                    // prev-id sentinel is enough to dedup without an
                    // auxiliary set.
                    std::optional<uint32_t> prev_sid;
                    for (const auto& rec : seg_expr_records) {
                        auto& sc = per_sample[rec.sample_id];
                        sc.expression_sum += static_cast<double>(rec.value);
                        if (prev_sid != rec.sample_id) {
                            sc.expressed_segments++;
                            prev_sid = rec.sample_id;
                        }
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
                acc.sources_seen_in_gene |= seg.sources;  // 8.7a: track sources per gene

                // 8.7a: per-source segment counts. source_count == 1 means
                // this segment was contributed by exactly one source, so it
                // is "exclusive" to that source.
                bool seg_source_exclusive = (seg.source_count() == 1);
                source_registry::instance().for_each(seg.sources,
                    [&](const std::string& src) {
                        per_source[src].segments++;
                        if (seg_source_exclusive) per_source[src].exclusive_segments++;
                    });

                size_t tx_count = seg.transcript_ids.size();
                acc.segment_tx_counts.push_back(tx_count);  // Phase 8.4: for effective isoforms
                total_transcripts += tx_count;
                for (uint32_t sid : seg.sample_idx) {
                    acc.sample_tx[sid] += tx_count;
                }

                // Transcript biotypes → global + per-sample.
                // Each (tx_id, biotype) appears in exactly one segment after
                // dedup, so the per-segment loop already gives the global
                // count without any extra deduplication state.
                for (const auto& [tx_id, biotype] : seg.transcript_biotypes) {
                    if (!biotype.empty()) {
                        transcripts_by_biotype[biotype]++;  // 8.7a: global biotype.tsv `total` column
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
                    size_t chain_pos = 1;
                    size_t chain_total = static_cast<size_t>(seg.exon_count);
                    // Phase 8.6: capture the ordered chain for this segment so
                    // splicing_catalog::detect_gene_events can run at gene
                    // finalization. Only populated when the splicing events
                    // stream is armed — avoids the per-gene chain memory for
                    // inspect runs that don't need events.
                    std::vector<key_ptr> chain_for_event_detection;
                    // Capture if either splicing events OR splicing hubs are armed.
                    // Hubs need the chain for per-hub PSI + entropy (Phase 8.6b),
                    // which walks segments_in_gene at finalize time.
                    bool capture_chain = (splicing_events_stream && splicing_events_stream->is_open())
                                      || (hub_stream && hub_stream->is_open());
                    if (capture_chain) chain_for_event_detection.reserve(chain_total);
                    while (current) {
                        if (capture_chain) chain_for_event_detection.push_back(current);
                        // Phase 8.4: every visit counts toward the per-gene
                        // exon-usage map, regardless of global first-visit
                        // dedup (the map is keyed by exon pointer so repeats
                        // within a segment don't over-count).
                        acc.exon_seg_counts[current]++;

                        // Quantification: attribute this segment's per-sample
                        // expression to every exon in its chain. The result
                        // is "total read support flowing through this exon
                        // across all segments containing it, per sample".
                        // Cleared per-chromosome alongside `active_genes`.
                        if (qtx_reader && !seg_expr_records.empty()) {
                            auto& exon_map = exon_expr_sum[current];
                            for (const auto& rec : seg_expr_records) {
                                exon_map[rec.sample_id] += rec.value;
                            }
                        }

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

                            // 8.7a: per-source exon counts. Only emitted at
                            // first global visit so each unique exon is
                            // counted once per source (matches segment-side
                            // semantics where each unique segment is one
                            // contribution).
                            bool exon_source_exclusive = (exon.source_count() == 1);
                            source_registry::instance().for_each(exon.sources,
                                [&](const std::string& src) {
                                    per_source[src].exons++;
                                    if (exon_source_exclusive) per_source[src].exclusive_exons++;
                                });

                            // Phase 8.4: remember this exon for constitutive
                            // classification at gene finalization
                            acc.first_visit_exons.push_back(current);

                            // Phase 8.3: hub detection. On first visit only,
                            // query ALL outgoing EXON_TO_EXON edges (not
                            // filtered to the current segment) so we see the
                            // full branching fan-out. Deduplicate targets,
                            // register as pending if size > MIN_HUB_BRANCHES.
                            if (hub_stream && hub_stream->is_open()) {
                                auto all_next = grove.get_neighbors_if(current,
                                    [](const edge_metadata& e) {
                                        return e.type == edge_metadata::edge_type::EXON_TO_EXON;
                                    });
                                std::unordered_set<key_ptr> unique_targets(
                                    all_next.begin(), all_next.end());
                                if (unique_targets.size() > MIN_HUB_BRANCHES) {
                                    pending_hub ph;
                                    ph.exon = current;
                                    ph.targets.assign(unique_targets.begin(), unique_targets.end());
                                    ph.chain_pos = chain_pos;
                                    ph.chain_total = chain_total;
                                    if (ph.targets.size() > acc.max_hub_targets_in_gene) {
                                        acc.max_hub_targets_in_gene = ph.targets.size();
                                    }
                                    acc.pending_hubs.push_back(std::move(ph));
                                }
                            }

                            // Phase 8.5: stream conserved exon row inline.
                            if (exon_conserved && conserved_exon_stream
                                && conserved_exon_stream->is_open()) {
                                auto& out = *conserved_exon_stream;
                                // Gene info from the parent segment
                                out << exon.id << "\t"
                                    << seg.gene_name() << "\t"
                                    << seg.gene_id() << "\t"
                                    << seqid << "\t"
                                    << format_coordinate(seqid, current->get_value()) << "\t"
                                    << exon.transcript_ids.size();
                                if (conserved_emit_expression) {
                                    auto exon_map_it = qtx_reader
                                        ? exon_expr_sum.find(current)
                                        : exon_expr_sum.end();
                                    for (size_t i = 0; i < conserved_stream_sample_ids.size(); ++i) {
                                        if (!conserved_stream_is_sample[i]) continue;
                                        uint32_t sid = conserved_stream_sample_ids[i];
                                        if (exon_map_it != exon_expr_sum.end()) {
                                            auto v_it = exon_map_it->second.find(sid);
                                            if (v_it != exon_map_it->second.end()) {
                                                out << "\t" << v_it->second;
                                                continue;
                                            }
                                        }
                                        out << "\t.";
                                    }
                                }
                                out << "\n";
                            }
                        }

                        auto next = grove.get_neighbors_if(current,
                            [edge_id](const edge_metadata& e) {
                                return e.type == edge_metadata::edge_type::EXON_TO_EXON
                                    && e.id == edge_id;
                            });
                        current = next.empty() ? nullptr : next.front();
                        chain_pos++;
                    }

                    // Phase 8.6: push the captured chain onto the per-gene
                    // list so splicing_catalog::detect_gene_events can run
                    // at finalization. structure_key left empty — detectors
                    // only read `segment` and `exon_chain`.
                    if (capture_chain) {
                        segment_chain_entry entry;
                        entry.segment = key;
                        entry.exon_chain = std::move(chain_for_event_detection);
                        acc.segments_in_gene.push_back(std::move(entry));
                    }
                }
            }
            node = node->get_next();
        }

        // Finalize all genes for this chromosome
        for (auto& [gidx, acc] : active_genes) {
            finalize_gene(acc, seqid);
        }
        active_genes.clear();
        exon_expr_sum.clear();
    }

    total_exons = visited_exons.size();
    total_edges = grove.edge_count();

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

    out << "metric\tvalue\n" << std::fixed;
    out << "samples\t" << per_sample.size() << "\n";
    out << "genes\t" << total_genes << "\n";
    out << "source_transcript_ids\t" << total_transcripts << "\n";
    out << "segments/isoforms\t" << total_segments << "\n";
    out << "exons\t" << total_exons << "\n";
    out << "graph_edges\t" << total_edges << "\n";
    out << "single_exon_segments\t" << single_exon << "\n";
    out << "single_isoform_genes\t" << single_iso << "\n";
    out << "multi_isoform_genes\t" << multi_iso << "\n";
    out << "mean_segments_per_gene\t" << std::setprecision(2) << mean_tpg << "\n";
    out << "median_segments_per_gene\t" << med_tpg << "\n";
    out << "max_segments_per_gene\t" << max_tpg << "\n";
    out << "mean_exons_per_segment\t" << mean_eps << "\n";
    out << "median_exons_per_segment\t" << med_eps << "\n";
    out << "max_exons_per_segment\t" << max_eps << "\n";

    logging::info("Overview written to: " + path);
}

void analysis_report::write_per_sample(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) return;

    auto& registry = sample_registry::instance();

    // Note: `transcripts_touched` counts transcript_ids summed over every
    // segment this sample contributes to. Shared segments over-count —
    // each sample in a segment's sample_idx is credited with the segment's
    // full transcript set, including transcripts originally contributed by
    // other samples. This is an intentional streaming-friendly approximation;
    // strict per-sample attribution would require a tx→sample map.
    // expressed_segments / mean_expression appear only when the qtx
    // reader was plugged in at collect() time. We detect this by checking
    // whether any sample has expressed_segments > 0; if so we emit both
    // columns. Otherwise we omit them so no-sidecar runs produce a clean
    // TSV that matches the legacy (pre-sidecar) shape.
    bool emit_expression_stats = false;
    for (const auto& sc : per_sample) {
        if (sc.expressed_segments > 0) { emit_expression_stats = true; break; }
    }

    out << "sample\ttype\tgenes\tsegments\texclusive_segments\tshared_segments"
        << "\tconserved_segments\texons\texclusive_exons\tshared_exons"
        << "\tconserved_exons\tconstitutive_exons\talternative_exons"
        << "\ttranscripts_touched\tsingle_exon_segments"
        << "\tmean_gene_exon_entropy\tmean_effective_isoforms";
    if (emit_expression_stats) {
        out << "\texpressed_segments\tmean_expression";
    }
    out << "\n";
    out << std::fixed;

    for (size_t sid = 0; sid < per_sample.size(); ++sid) {
        const auto& sc = per_sample[sid];
        if (sc.segments == 0) continue;

        const auto& info = registry.get(static_cast<uint32_t>(sid));
        std::string label = info.id.empty() ? std::to_string(sid) : info.id;

        double mean_exon_H = (sc.multi_segment_genes > 0)
            ? sc.gene_exon_entropy_sum / static_cast<double>(sc.multi_segment_genes)
            : 0;
        double mean_eff_iso = (sc.multi_segment_genes > 0)
            ? sc.effective_isoforms_sum / static_cast<double>(sc.multi_segment_genes)
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
            << "\t" << sc.constitutive_exons
            << "\t" << sc.alternative_exons
            << "\t" << sc.transcripts
            << "\t" << sc.single_exon_segments
            << "\t" << std::setprecision(3) << mean_exon_H
            << "\t" << std::setprecision(2) << mean_eff_iso;
        if (emit_expression_stats) {
            double mean_expr = (sc.expressed_segments > 0)
                ? sc.expression_sum / static_cast<double>(sc.expressed_segments)
                : 0.0;
            out << "\t" << sc.expressed_segments
                << "\t" << std::setprecision(3) << mean_expr;
        }
        out << "\n";
    }

    logging::info("Per-sample stats written to: " + path);
}

void analysis_report::write_per_source(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open per-source file: " + path);
        return;
    }

    out << "source\tgenes\tsegments\texclusive_segments\texons\texclusive_exons\n";
    for (const auto& [source, ss] : per_source) {
        out << source
            << "\t" << ss.genes
            << "\t" << ss.segments
            << "\t" << ss.exclusive_segments
            << "\t" << ss.exons
            << "\t" << ss.exclusive_exons
            << "\n";
    }

    logging::info("Per-source stats written to: " + path);
}

void analysis_report::write_biotype(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open biotype file: " + path);
        return;
    }

    if (genes_by_biotype.empty() && transcripts_by_biotype.empty()) {
        logging::warning("No biotype data to write (TALON-only builds carry no biotype attrs)");
        return;
    }

    auto& registry = sample_registry::instance();

    // Enumerate non-replicate samples in registry order, matching the
    // column convention used by per_sample.tsv.
    std::vector<uint32_t> sample_ids;
    std::vector<std::string> labels;
    for (size_t sid = 0; sid < per_sample.size(); ++sid) {
        if (per_sample[sid].segments == 0 && per_sample[sid].exons == 0) continue;
        const auto& info = registry.get(static_cast<uint32_t>(sid));
        sample_ids.push_back(static_cast<uint32_t>(sid));
        labels.push_back(info.id.empty() ? std::to_string(sid) : info.id);
    }

    out << "level\tbiotype\ttotal";
    for (const auto& label : labels) out << "\t" << label;
    out << "\n";

    auto write_section = [&](const std::string& level,
                             const std::map<std::string, size_t>& global,
                             auto per_sample_getter)
    {
        // Sort by global count descending so the biggest biotypes come first.
        std::vector<std::pair<std::string, size_t>> sorted_bt(global.begin(), global.end());
        std::sort(sorted_bt.begin(), sorted_bt.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        for (const auto& [biotype, count] : sorted_bt) {
            out << level << "\t" << biotype << "\t" << count;
            for (uint32_t sid : sample_ids) {
                const auto& bt_map = per_sample_getter(per_sample[sid]);
                auto it = bt_map.find(biotype);
                out << "\t" << (it != bt_map.end() ? it->second : size_t{0});
            }
            out << "\n";
        }
    };

    write_section("gene", genes_by_biotype,
        [](const sample_counters& s) -> const std::map<std::string, size_t>& {
            return s.genes_by_biotype;
        });
    write_section("transcript", transcripts_by_biotype,
        [](const sample_counters& s) -> const std::map<std::string, size_t>& {
            return s.transcripts_by_biotype;
        });

    logging::info("Biotype stats written to: " + path);
}

// ── Splicing hub streaming setup (Phase 8.3) ───────────────────────

void analysis_report::begin_splicing_hub_streams(
    const std::string& hubs_path,
    const std::string& branches_path,
    bool emit_expression_columns)
{
    hub_emit_expression = emit_expression_columns;
    auto& registry = sample_registry::instance();

    hub_stream_sample_ids.clear();
    hub_stream_is_sample.clear();
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto& info = registry.get(static_cast<uint32_t>(i));
        if (info.type == "replicate") continue;
        hub_stream_sample_ids.push_back(static_cast<uint32_t>(i));
        hub_stream_is_sample.push_back(info.type == "sample");
    }

    hub_stream = std::make_unique<std::ofstream>(hubs_path);
    branch_stream = std::make_unique<std::ofstream>(branches_path);

    if (!hub_stream->is_open()) {
        logging::error("Cannot open splicing hubs file: " + hubs_path);
        hub_stream.reset();
    }
    if (!branch_stream->is_open()) {
        logging::error("Cannot open branch details file: " + branches_path);
        branch_stream.reset();
    }
    if (!hub_stream || !hub_stream->is_open()) {
        // Without the hub stream, inline detection is disabled. Drop branch too.
        hub_stream.reset();
        branch_stream.reset();
        return;
    }

    *hub_stream << std::fixed;
    if (branch_stream && branch_stream->is_open()) *branch_stream << std::fixed;

    // splicing_hubs.tsv header
    *hub_stream << "gene_name\tgene_id\texon_id\tcoordinate\texon_number\ttotal_exons"
                << "\ttotal_branches\ttotal_transcripts";
    for (size_t i = 0; i < hub_stream_sample_ids.size(); ++i) {
        const auto& info = registry.get(hub_stream_sample_ids[i]);
        std::string label = info.id.empty()
            ? std::to_string(hub_stream_sample_ids[i]) : info.id;
        *hub_stream << "\t" << label << ".branches"
                    << "\t" << label << ".shared"
                    << "\t" << label << ".unique"
                    << "\t" << label << ".psi"
                    << "\t" << label << ".entropy";
        if (hub_emit_expression && hub_stream_is_sample[i]) {
            *hub_stream << "\t" << label << "." << expr_type_label(info.expr_type);
        }
    }
    *hub_stream << "\n";

    // branch_details.tsv header
    if (branch_stream && branch_stream->is_open()) {
        *branch_stream << "hub_gene_name\thub_gene_id\thub_exon_id\thub_coordinate"
                       << "\ttarget_exon_id\ttarget_coordinate";
        for (size_t i = 0; i < hub_stream_sample_ids.size(); ++i) {
            const auto& info = registry.get(hub_stream_sample_ids[i]);
            std::string label = info.id.empty()
                ? std::to_string(hub_stream_sample_ids[i]) : info.id;
            *branch_stream << "\t" << label << ".present";
            if (hub_emit_expression && hub_stream_is_sample[i]) {
                *branch_stream << "\t" << label << "." << expr_type_label(info.expr_type);
            }
        }
        *branch_stream << "\n";
    }

    logging::info("Splicing hubs streaming to: " + hubs_path);
    logging::info("Branch details streaming to: " + branches_path);
}

// ── Splicing event streaming setup (Phase 8.6) ─────────────────────

void analysis_report::begin_splicing_events_stream(const std::string& path) {
    auto& registry = sample_registry::instance();

    splicing_events_sample_ids.clear();
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto& info = registry.get(static_cast<uint32_t>(i));
        if (info.type == "replicate") continue;
        splicing_events_sample_ids.push_back(static_cast<uint32_t>(i));
    }

    splicing_events_stream = std::make_unique<std::ofstream>(path);
    if (!splicing_events_stream->is_open()) {
        logging::error("Cannot open splicing events file: " + path);
        splicing_events_stream.reset();
        return;
    }

    auto& out = *splicing_events_stream;
    out << std::fixed;
    out << "gene_id\tgene_name\tevent_type\tchromosome\tupstream_exon"
        << "\tdownstream_exon\taffected_exons";
    for (uint32_t sid : splicing_events_sample_ids) {
        const auto& info = registry.get(sid);
        std::string label = info.id.empty() ? std::to_string(sid) : info.id;
        out << "\t" << label << ".psi"
            << "\t" << label << ".included"
            << "\t" << label << ".total";
    }
    out << "\n";

    logging::info("Splicing events streaming to: " + path);
}

// ── Conserved exon streaming setup (Phase 8.5) ──────────────────────

void analysis_report::begin_conserved_exon_stream(const std::string& path,
                                                  bool emit_expression_columns) {
    conserved_emit_expression = emit_expression_columns;
    auto& registry = sample_registry::instance();

    conserved_stream_sample_ids.clear();
    conserved_stream_is_sample.clear();

    // Enumerate non-replicate entries in registry order
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto& info = registry.get(static_cast<uint32_t>(i));
        if (info.type == "replicate") continue;
        conserved_stream_sample_ids.push_back(static_cast<uint32_t>(i));
        conserved_stream_is_sample.push_back(info.type == "sample");
    }

    conserved_exon_stream = std::make_unique<std::ofstream>(path);
    if (!conserved_exon_stream->is_open()) {
        logging::error("Cannot open conserved exons file: " + path);
        conserved_exon_stream.reset();
        return;
    }

    auto& out = *conserved_exon_stream;
    out << std::fixed;
    out << "exon_id\tgene_name\tgene_id\tchromosome\tcoordinate\tn_transcripts";
    if (conserved_emit_expression) {
        // Per-sample expression columns for sample-typed entries only
        for (size_t i = 0; i < conserved_stream_sample_ids.size(); ++i) {
            if (!conserved_stream_is_sample[i]) continue;
            const auto& info = registry.get(conserved_stream_sample_ids[i]);
            std::string label = info.id.empty()
                ? std::to_string(conserved_stream_sample_ids[i]) : info.id;
            out << "\t" << label << "." << expr_type_label(info.expr_type);
        }
    }
    out << "\n";

    logging::info("Conserved exons streaming to: " + path);
}

// ── Sharing TSVs (Phase 8.2) ────────────────────────────────────────

void analysis_report::write_exon_sharing(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open exon sharing file: " + path);
        return;
    }

    auto& registry = sample_registry::instance();

    // Collect sample IDs that actually have features, in index order
    std::vector<uint32_t> sample_ids;
    std::vector<std::string> labels;
    for (size_t sid = 0; sid < per_sample.size(); ++sid) {
        if (per_sample[sid].segments == 0 && per_sample[sid].exons == 0) continue;
        const auto& info = registry.get(static_cast<uint32_t>(sid));
        sample_ids.push_back(static_cast<uint32_t>(sid));
        labels.push_back(info.id.empty() ? std::to_string(sid) : info.id);
    }

    out << "metric\ttotal";
    for (const auto& label : labels) out << "\t" << label;
    out << "\n";

    auto write_row = [&](const std::string& metric, size_t global,
                         auto getter) {
        out << metric << "\t" << global;
        for (uint32_t sid : sample_ids) {
            out << "\t" << getter(per_sample[sid]);
        }
        out << "\n";
    };

    write_row("total",        total_exons, [](const sample_counters& s) { return s.exons; });
    write_row("exclusive",    size_t{0},   [](const sample_counters& s) { return s.exclusive_exons; });
    write_row("shared",       size_t{0},   [](const sample_counters& s) { return s.shared_exons; });
    write_row("conserved",    size_t{0},   [](const sample_counters& s) { return s.conserved_exons; });
    write_row("constitutive", size_t{0},   [](const sample_counters& s) { return s.constitutive_exons; });
    write_row("alternative",  size_t{0},   [](const sample_counters& s) { return s.alternative_exons; });

    logging::info("Exon sharing written to: " + path);
}

void analysis_report::write_segment_sharing(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open segment sharing file: " + path);
        return;
    }

    auto& registry = sample_registry::instance();

    std::vector<uint32_t> sample_ids;
    std::vector<std::string> labels;
    for (size_t sid = 0; sid < per_sample.size(); ++sid) {
        if (per_sample[sid].segments == 0 && per_sample[sid].exons == 0) continue;
        const auto& info = registry.get(static_cast<uint32_t>(sid));
        sample_ids.push_back(static_cast<uint32_t>(sid));
        labels.push_back(info.id.empty() ? std::to_string(sid) : info.id);
    }

    out << "metric\ttotal";
    for (const auto& label : labels) out << "\t" << label;
    out << "\n";

    size_t total_segments = exons_per_segment.size();
    auto write_row = [&](const std::string& metric, size_t global,
                         auto getter) {
        out << metric << "\t" << global;
        for (uint32_t sid : sample_ids) {
            out << "\t" << getter(per_sample[sid]);
        }
        out << "\n";
    };

    write_row("total",     total_segments, [](const sample_counters& s) { return s.segments; });
    write_row("exclusive", size_t{0},      [](const sample_counters& s) { return s.exclusive_segments; });
    write_row("shared",    size_t{0},      [](const sample_counters& s) { return s.shared_segments; });
    write_row("conserved", size_t{0},      [](const sample_counters& s) { return s.conserved_segments; });

    logging::info("Segment sharing written to: " + path);
}
