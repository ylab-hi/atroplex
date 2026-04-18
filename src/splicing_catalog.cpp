/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "splicing_catalog.hpp"
#include "utility.hpp"

#include <fstream>
#include <set>
#include <algorithm>
#include <cmath>

// ========================================================================
// Splicing event detection
//
// For each gene, we compare exon chains across all segments (transcript paths).
// Events are detected by analyzing which exons are shared, skipped, or
// alternative between different paths.
//
// Key data structure: each segment_chain_entry has a pre-built exon_chain
// (vector<key_ptr>) where pointer identity means the same exon object.
// ========================================================================

std::vector<splicing_event> splicing_catalog::collect_from_grove(grove_type& grove) {
    // Build per-gene segment chains by walking the grove
    using gene_index_type = std::unordered_map<std::string, std::vector<segment_chain_entry>>;
    std::map<std::string, gene_index_type> gene_indices;

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
                if (seg.absorbed) continue;

                // Traverse exon chain
                size_t edge_id = seg.segment_index;
                std::vector<key_ptr> exon_chain;

                auto first_exons = grove.get_neighbors_if(key,
                    [edge_id](const edge_metadata& e) {
                        return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON
                            && e.id == edge_id;
                    });

                if (!first_exons.empty()) {
                    auto* current = first_exons.front();
                    while (current) {
                        exon_chain.push_back(current);
                        auto next = grove.get_neighbors_if(current,
                            [edge_id](const edge_metadata& e) {
                                return e.type == edge_metadata::edge_type::EXON_TO_EXON
                                    && e.id == edge_id;
                            });
                        current = next.empty() ? nullptr : next.front();
                    }
                }

                // Build structure key from exon coordinates
                std::string structure_key;
                if (!exon_chain.empty()) {
                    structure_key = format_coordinate(seqid, exon_chain.front()->get_value());
                    for (size_t i = 1; i < exon_chain.size(); ++i) {
                        structure_key += "," + format_coordinate(seqid, exon_chain[i]->get_value());
                    }
                }

                gene_indices[seqid][seg.gene_id()].push_back(
                    {key, std::move(exon_chain), std::move(structure_key)});
            }
            node = node->get_next();
        }
    }

    // Run event detection per gene
    std::vector<splicing_event> all_events;
    for (const auto& [chrom, gene_idx] : gene_indices) {
        for (const auto& [gene_id, segments] : gene_idx) {
            std::vector<segment_chain_entry> valid;
            for (const auto& entry : segments) {
                auto& seg = get_segment(entry.segment->get_data());
                if (seg.absorbed) continue;
                if (entry.exon_chain.size() < 2) continue;
                valid.push_back(entry);
            }
            if (valid.size() < 2) continue;

            auto gene_events = detect_gene_events(gene_id, chrom, valid, grove);
            all_events.insert(all_events.end(),
                std::make_move_iterator(gene_events.begin()),
                std::make_move_iterator(gene_events.end()));
        }
    }
    return all_events;
}

std::vector<splicing_event> splicing_catalog::detect_gene_events(
    const std::string& gene_id,
    const std::string& chromosome,
    const std::vector<segment_chain_entry>& segments,
    grove_type& /*grove*/
) {
    std::vector<splicing_event> events;

    detect_cassette_exons(gene_id, chromosome, segments, events);
    detect_alt_splice_sites(gene_id, chromosome, segments, events);
    detect_intron_retention(gene_id, chromosome, segments, events);
    detect_alt_terminal_exons(gene_id, chromosome, segments, events);
    detect_mutually_exclusive(gene_id, chromosome, segments, events);

    return events;
}

// ========================================================================
// Cassette exon detection
//
// An exon is a cassette exon if:
// 1. It is present in some transcript paths but absent in others
// 2. The flanking exons (upstream and downstream) are shared between the
//    including and skipping paths
//
// Algorithm: for each exon, check if it appears in all segments. If not,
// find paths that skip it — those must have adjacent exons that bridge
// across the skipped exon.
// ========================================================================

void splicing_catalog::detect_cassette_exons(
    const std::string& gene_id,
    const std::string& chromosome,
    const std::vector<segment_chain_entry>& segments,
    std::vector<splicing_event>& events
) {
    // Build a set of all unique exons across all segments (by pointer)
    std::set<key_ptr> all_exons;
    for (const auto& seg : segments) {
        for (auto* exon : seg.exon_chain) {
            all_exons.insert(exon);
        }
    }

    // For each exon, check if it's skipped by any segment
    for (auto* candidate : all_exons) {
        // Find segments that include vs skip this exon
        std::vector<const segment_chain_entry*> including;
        std::vector<const segment_chain_entry*> skipping;

        for (const auto& seg : segments) {
            bool found = false;
            for (auto* exon : seg.exon_chain) {
                if (exon == candidate) {
                    found = true;
                    break;
                }
            }
            if (found) {
                including.push_back(&seg);
            } else {
                skipping.push_back(&seg);
            }
        }

        // Not a cassette if present in all or none
        if (skipping.empty() || including.empty()) continue;

        // Don't count first/last exons as cassettes (those are alt terminal exons)
        bool is_first_in_any = false;
        bool is_last_in_any = false;
        for (const auto* seg : including) {
            if (seg->exon_chain.front() == candidate) is_first_in_any = true;
            if (seg->exon_chain.back() == candidate) is_last_in_any = true;
        }
        if (is_first_in_any || is_last_in_any) continue;

        // Verify that skipping paths bridge across: the exons flanking the
        // cassette in the including path must be adjacent in a skipping path
        bool has_bridging_path = false;
        for (const auto* inc_seg : including) {
            // Find position of candidate in including segment
            for (size_t i = 0; i < inc_seg->exon_chain.size(); ++i) {
                if (inc_seg->exon_chain[i] != candidate) continue;
                if (i == 0 || i == inc_seg->exon_chain.size() - 1) continue;

                key_ptr upstream = inc_seg->exon_chain[i - 1];
                key_ptr downstream = inc_seg->exon_chain[i + 1];

                // Check if any skipping path has upstream→downstream adjacency
                for (const auto* skip_seg : skipping) {
                    for (size_t j = 0; j + 1 < skip_seg->exon_chain.size(); ++j) {
                        if (skip_seg->exon_chain[j] == upstream &&
                            skip_seg->exon_chain[j + 1] == downstream) {
                            has_bridging_path = true;
                            break;
                        }
                    }
                    if (has_bridging_path) break;
                }
                if (has_bridging_path) break;
            }
            if (has_bridging_path) break;
        }

        if (!has_bridging_path) continue;

        // Create the event
        splicing_event event;
        event.type = splicing_event_type::CASSETTE_EXON;
        event.gene_id = gene_id;
        auto& exon_data = get_exon(candidate->get_data());
        event.gene_name = get_segment(segments.front().segment->get_data()).gene_name();
        event.chromosome = chromosome;
        event.affected_exons.push_back(candidate);

        // Find upstream/downstream from the first including segment
        for (const auto* inc_seg : including) {
            for (size_t i = 1; i + 1 < inc_seg->exon_chain.size(); ++i) {
                if (inc_seg->exon_chain[i] == candidate) {
                    event.upstream_exon = inc_seg->exon_chain[i - 1];
                    event.downstream_exon = inc_seg->exon_chain[i + 1];
                    break;
                }
            }
            if (event.upstream_exon) break;
        }

        compute_psi(event, segments);
        events.push_back(std::move(event));
    }
}

// ========================================================================
// Alternative 5'/3' splice site detection
//
// Alt 5' (ALT_5PRIME): two exons share the same acceptor (start) but have
//   different donors (end). This means the upstream intron boundary varies.
//
// Alt 3' (ALT_3PRIME): two exons share the same donor (end) but have
//   different acceptors (start). The downstream intron boundary varies.
//
// Algorithm: group exons by shared boundaries, find groups with >1 exon.
// ========================================================================

void splicing_catalog::detect_alt_splice_sites(
    const std::string& gene_id,
    const std::string& chromosome,
    const std::vector<segment_chain_entry>& segments,
    std::vector<splicing_event>& events
) {
    // Collect all internal exons (not first/last in any segment)
    std::set<key_ptr> first_exons, last_exons;
    for (const auto& seg : segments) {
        first_exons.insert(seg.exon_chain.front());
        last_exons.insert(seg.exon_chain.back());
    }

    // Group exons by shared acceptor (start position) → alt 5' donor
    std::map<size_t, std::vector<key_ptr>> by_acceptor;
    // Group exons by shared donor (end position) → alt 3' acceptor
    std::map<size_t, std::vector<key_ptr>> by_donor;

    std::set<key_ptr> seen;
    for (const auto& seg : segments) {
        for (auto* exon : seg.exon_chain) {
            if (!seen.insert(exon).second) continue;
            // Skip terminal exons (handled by alt first/last exon detection)
            if (first_exons.count(exon) || last_exons.count(exon)) continue;

            auto& coord = exon->get_value();
            by_acceptor[coord.get_start()].push_back(exon);
            by_donor[coord.get_end()].push_back(exon);
        }
    }

    // Alt 5' splice site: same acceptor, different donors
    for (const auto& [acceptor, exon_group] : by_acceptor) {
        if (exon_group.size() < 2) continue;

        // Verify they actually have different end positions
        std::set<size_t> donors;
        for (auto* e : exon_group) {
            donors.insert(e->get_value().get_end());
        }
        if (donors.size() < 2) continue;

        splicing_event event;
        event.type = splicing_event_type::ALT_5PRIME;
        event.gene_id = gene_id;
        auto& exon_data = get_exon(exon_group[0]->get_data());
        event.gene_name = get_segment(segments.front().segment->get_data()).gene_name();
        event.chromosome = chromosome;
        event.affected_exons = exon_group;

        compute_psi(event, segments);
        events.push_back(std::move(event));
    }

    // Alt 3' splice site: same donor, different acceptors
    for (const auto& [donor, exon_group] : by_donor) {
        if (exon_group.size() < 2) continue;

        std::set<size_t> acceptors;
        for (auto* e : exon_group) {
            acceptors.insert(e->get_value().get_start());
        }
        if (acceptors.size() < 2) continue;

        splicing_event event;
        event.type = splicing_event_type::ALT_3PRIME;
        event.gene_id = gene_id;
        auto& exon_data = get_exon(exon_group[0]->get_data());
        event.gene_name = get_segment(segments.front().segment->get_data()).gene_name();
        event.chromosome = chromosome;
        event.affected_exons = exon_group;

        compute_psi(event, segments);
        events.push_back(std::move(event));
    }
}

// ========================================================================
// Intron retention detection
//
// An intron retention event occurs when:
// 1. Some transcript paths have two adjacent exons (intron spliced out)
// 2. Another transcript path has a single exon spanning the intron region
//    (intron retained — the exon covers both flanking exons + the intron)
//
// Algorithm: for each intron (gap between adjacent exons in any segment),
// check if any other segment has an exon that spans across it.
// ========================================================================

void splicing_catalog::detect_intron_retention(
    const std::string& gene_id,
    const std::string& chromosome,
    const std::vector<segment_chain_entry>& segments,
    std::vector<splicing_event>& events
) {
    // Collect all unique introns defined by adjacent exon pairs
    // Key: (donor_end, acceptor_start) = intron boundaries
    struct intron_def {
        key_ptr upstream_exon;
        key_ptr downstream_exon;
        size_t donor;    // end of upstream exon
        size_t acceptor; // start of downstream exon
    };

    std::vector<intron_def> introns;
    std::set<std::pair<size_t, size_t>> seen_introns;

    for (const auto& seg : segments) {
        for (size_t i = 0; i + 1 < seg.exon_chain.size(); ++i) {
            auto& up_coord = seg.exon_chain[i]->get_value();
            auto& down_coord = seg.exon_chain[i + 1]->get_value();
            size_t donor = up_coord.get_end();
            size_t acceptor = down_coord.get_start();

            if (seen_introns.insert({donor, acceptor}).second) {
                introns.push_back({seg.exon_chain[i], seg.exon_chain[i + 1],
                                   donor, acceptor});
            }
        }
    }

    // For each intron, check if any exon in any segment spans across it
    for (const auto& intron : introns) {
        std::vector<const segment_chain_entry*> retaining;  // segments with retention
        std::vector<const segment_chain_entry*> splicing;   // segments with normal splicing

        for (const auto& seg : segments) {
            bool has_intron_spliced = false;
            bool has_intron_retained = false;

            for (size_t i = 0; i < seg.exon_chain.size(); ++i) {
                auto& coord = seg.exon_chain[i]->get_value();

                // Check if this exon spans the intron
                if (coord.get_start() <= intron.donor && coord.get_end() >= intron.acceptor) {
                    // This exon covers both the upstream and downstream boundaries
                    // → intron is retained in this transcript
                    has_intron_retained = true;
                }

                // Check if this segment has the normal splice (adjacent exons)
                if (i + 1 < seg.exon_chain.size()) {
                    auto& next_coord = seg.exon_chain[i + 1]->get_value();
                    if (coord.get_end() == intron.donor &&
                        next_coord.get_start() == intron.acceptor) {
                        has_intron_spliced = true;
                    }
                }
            }

            if (has_intron_retained) retaining.push_back(&seg);
            if (has_intron_spliced) splicing.push_back(&seg);
        }

        // Only report if both retained and spliced forms exist
        if (retaining.empty() || splicing.empty()) continue;

        splicing_event event;
        event.type = splicing_event_type::INTRON_RETENTION;
        event.gene_id = gene_id;
        auto& exon_data = get_exon(intron.upstream_exon->get_data());
        event.gene_name = get_segment(segments.front().segment->get_data()).gene_name();
        event.chromosome = chromosome;
        event.upstream_exon = intron.upstream_exon;
        event.downstream_exon = intron.downstream_exon;

        // Find the retaining exon(s) as affected exons
        for (const auto* seg : retaining) {
            for (auto* exon : seg->exon_chain) {
                auto& coord = exon->get_value();
                if (coord.get_start() <= intron.donor && coord.get_end() >= intron.acceptor) {
                    // Avoid duplicates
                    bool already = false;
                    for (auto* e : event.affected_exons) {
                        if (e == exon) { already = true; break; }
                    }
                    if (!already) event.affected_exons.push_back(exon);
                }
            }
        }

        compute_psi(event, segments);
        events.push_back(std::move(event));
    }
}

// ========================================================================
// Alternative first/last exon detection
//
// Alt first exon: different segments start with different exons
// Alt last exon: different segments end with different exons
// ========================================================================

void splicing_catalog::detect_alt_terminal_exons(
    const std::string& gene_id,
    const std::string& chromosome,
    const std::vector<segment_chain_entry>& segments,
    std::vector<splicing_event>& events
) {
    // Collect unique first and last exons
    std::set<key_ptr> first_exons, last_exons;
    for (const auto& seg : segments) {
        first_exons.insert(seg.exon_chain.front());
        last_exons.insert(seg.exon_chain.back());
    }

    // Alt first exon: >1 unique first exon
    if (first_exons.size() > 1) {
        splicing_event event;
        event.type = splicing_event_type::ALT_FIRST_EXON;
        event.gene_id = gene_id;
        auto& exon_data = get_exon(segments[0].exon_chain.front()->get_data());
        event.gene_name = get_segment(segments.front().segment->get_data()).gene_name();
        event.chromosome = chromosome;
        event.affected_exons.assign(first_exons.begin(), first_exons.end());

        compute_psi(event, segments);
        events.push_back(std::move(event));
    }

    // Alt last exon: >1 unique last exon
    if (last_exons.size() > 1) {
        splicing_event event;
        event.type = splicing_event_type::ALT_LAST_EXON;
        event.gene_id = gene_id;
        auto& exon_data = get_exon(segments[0].exon_chain.back()->get_data());
        event.gene_name = get_segment(segments.front().segment->get_data()).gene_name();
        event.chromosome = chromosome;
        event.affected_exons.assign(last_exons.begin(), last_exons.end());

        compute_psi(event, segments);
        events.push_back(std::move(event));
    }
}

// ========================================================================
// Mutually exclusive exon detection
//
// Two exons are mutually exclusive if:
// 1. They never co-occur in the same transcript path
// 2. They occupy the same position in the gene structure (flanked by the
//    same upstream and downstream exons)
// ========================================================================

void splicing_catalog::detect_mutually_exclusive(
    const std::string& gene_id,
    const std::string& chromosome,
    const std::vector<segment_chain_entry>& segments,
    std::vector<splicing_event>& events
) {
    // Build exon→segment membership (which segments contain each exon)
    std::map<key_ptr, std::set<size_t>> exon_to_segments;
    for (size_t s = 0; s < segments.size(); ++s) {
        for (auto* exon : segments[s].exon_chain) {
            exon_to_segments[exon].insert(s);
        }
    }

    // Find pairs of exons that never co-occur
    std::vector<key_ptr> all_exons;
    for (const auto& [exon, _] : exon_to_segments) {
        all_exons.push_back(exon);
    }

    for (size_t i = 0; i < all_exons.size(); ++i) {
        for (size_t j = i + 1; j < all_exons.size(); ++j) {
            auto& segs_i = exon_to_segments[all_exons[i]];
            auto& segs_j = exon_to_segments[all_exons[j]];

            // Check for zero overlap in segment membership
            bool co_occur = false;
            for (size_t s : segs_i) {
                if (segs_j.count(s)) {
                    co_occur = true;
                    break;
                }
            }
            if (co_occur) continue;

            // Must each be in at least one segment
            if (segs_i.empty() || segs_j.empty()) continue;

            // Verify they're flanked by the same upstream/downstream exons
            // in their respective segments
            key_ptr shared_upstream = nullptr;
            key_ptr shared_downstream = nullptr;

            for (size_t s : segs_i) {
                const auto& chain = segments[s].exon_chain;
                for (size_t k = 0; k < chain.size(); ++k) {
                    if (chain[k] != all_exons[i]) continue;
                    key_ptr up = (k > 0) ? chain[k - 1] : nullptr;
                    key_ptr down = (k + 1 < chain.size()) ? chain[k + 1] : nullptr;

                    // Check if exon j has the same flanking in any of its segments
                    for (size_t s2 : segs_j) {
                        const auto& chain2 = segments[s2].exon_chain;
                        for (size_t m = 0; m < chain2.size(); ++m) {
                            if (chain2[m] != all_exons[j]) continue;
                            key_ptr up2 = (m > 0) ? chain2[m - 1] : nullptr;
                            key_ptr down2 = (m + 1 < chain2.size()) ? chain2[m + 1] : nullptr;

                            if (up && up2 && up == up2 && down && down2 && down == down2) {
                                shared_upstream = up;
                                shared_downstream = down;
                            }
                        }
                    }
                }
            }

            if (!shared_upstream || !shared_downstream) continue;

            splicing_event event;
            event.type = splicing_event_type::MUTUALLY_EXCLUSIVE;
            event.gene_id = gene_id;
            auto& exon_data = get_exon(all_exons[i]->get_data());
            event.gene_name = get_segment(segments.front().segment->get_data()).gene_name();
            event.chromosome = chromosome;
            event.upstream_exon = shared_upstream;
            event.downstream_exon = shared_downstream;
            event.affected_exons = {all_exons[i], all_exons[j]};

            compute_psi(event, segments);
            events.push_back(std::move(event));
        }
    }
}

// ========================================================================
// PSI computation
//
// PSI = fraction of transcripts at the locus that include the event.
// For each sample, count how many segments include any affected exon
// vs total segments for the gene in that sample.
// ========================================================================

void splicing_catalog::compute_psi(
    splicing_event& event,
    const std::vector<segment_chain_entry>& segments
) {
    // Collect all samples across all segments
    std::set<uint32_t> all_samples;
    for (const auto& seg : segments) {
        auto& s = get_segment(seg.segment->get_data());
        s.sample_idx.for_each([&](uint32_t sid) {
            all_samples.insert(sid);
        });
    }

    // Build set of affected exon pointers for fast lookup
    std::set<key_ptr> affected_set(event.affected_exons.begin(), event.affected_exons.end());

    for (uint32_t sid : all_samples) {
        size_t total = 0;
        size_t included = 0;

        for (const auto& seg : segments) {
            auto& s = get_segment(seg.segment->get_data());
            if (!s.sample_idx.test(sid)) continue;
            total++;

            // Check if this segment includes any affected exon
            for (auto* exon : seg.exon_chain) {
                if (affected_set.count(exon)) {
                    included++;
                    break;
                }
            }
        }

        if (total > 0) {
            event.sample_psi[sid] = static_cast<float>(included) / static_cast<float>(total);
            event.sample_included[sid] = included;
            event.sample_total[sid] = total;
        }
    }
}

// ========================================================================
// Output
// ========================================================================

std::string splicing_catalog::event_type_str(splicing_event_type type) {
    switch (type) {
        case splicing_event_type::CASSETTE_EXON: return "cassette_exon";
        case splicing_event_type::MUTUALLY_EXCLUSIVE: return "mutually_exclusive";
        case splicing_event_type::ALT_5PRIME: return "alt_5prime";
        case splicing_event_type::ALT_3PRIME: return "alt_3prime";
        case splicing_event_type::INTRON_RETENTION: return "intron_retention";
        case splicing_event_type::ALT_FIRST_EXON: return "alt_first_exon";
        case splicing_event_type::ALT_LAST_EXON: return "alt_last_exon";
    }
    return "unknown";
}

std::string splicing_catalog::format_exon(key_ptr exon) {
    if (!exon) return ".";
    auto& coord = exon->get_value();
    return std::to_string(coord.get_start()) + "-" + std::to_string(coord.get_end());
}

void splicing_catalog::write_events_tsv(
    const std::string& filepath,
    const std::vector<splicing_event>& events
) {
    std::ofstream out(filepath);
    if (!out.is_open()) return;

    // Collect all sample IDs across all events
    std::set<uint32_t> all_samples;
    for (const auto& event : events) {
        for (const auto& [sid, _] : event.sample_psi) {
            all_samples.insert(sid);
        }
    }

    // Get sample labels
    auto& registry = sample_registry::instance();
    std::vector<uint32_t> sample_ids(all_samples.begin(), all_samples.end());

    // Header
    out << "gene_id\tgene_name\tevent_type\tchromosome\tupstream_exon\tdownstream_exon\taffected_exons";
    for (uint32_t sid : sample_ids) {
        const auto& info = registry.get(sid);
        std::string label = info.id.empty() ? std::to_string(sid) : info.id;
        out << "\t" << label << ".psi"
            << "\t" << label << ".included"
            << "\t" << label << ".total";
    }
    out << "\n";

    // Rows
    for (const auto& event : events) {
        out << event.gene_id
            << "\t" << event.gene_name
            << "\t" << event_type_str(event.type)
            << "\t" << event.chromosome
            << "\t" << format_exon(event.upstream_exon)
            << "\t" << format_exon(event.downstream_exon);

        // Affected exons (comma-separated)
        out << "\t";
        for (size_t i = 0; i < event.affected_exons.size(); ++i) {
            if (i > 0) out << ",";
            out << format_exon(event.affected_exons[i]);
        }

        // Per-sample columns
        for (uint32_t sid : sample_ids) {
            auto psi_it = event.sample_psi.find(sid);
            auto inc_it = event.sample_included.find(sid);
            auto tot_it = event.sample_total.find(sid);

            if (psi_it != event.sample_psi.end()) {
                out << "\t" << psi_it->second;
            } else {
                out << "\t.";
            }
            if (inc_it != event.sample_included.end()) {
                out << "\t" << inc_it->second;
            } else {
                out << "\t.";
            }
            if (tot_it != event.sample_total.end()) {
                out << "\t" << tot_it->second;
            } else {
                out << "\t.";
            }
        }

        out << "\n";
    }
}

void splicing_catalog::write_summary(
    const std::string& filepath,
    const std::vector<splicing_event>& events
) {
    std::ofstream out(filepath);
    if (!out.is_open()) return;

    // Count events by type
    std::map<splicing_event_type, size_t> counts;
    for (const auto& event : events) {
        counts[event.type]++;
    }

    out << "# Atroplex Splicing Event Catalog Summary\n\n";
    out << "Total events: " << events.size() << "\n\n";
    out << "## Event Types\n";

    auto print_count = [&](splicing_event_type type, const std::string& name) {
        size_t n = counts[type];
        float pct = events.empty() ? 0.0f : 100.0f * n / events.size();
        out << name << ": " << n << " (" << pct << "%)\n";
    };

    print_count(splicing_event_type::CASSETTE_EXON, "Cassette exon");
    print_count(splicing_event_type::MUTUALLY_EXCLUSIVE, "Mutually exclusive");
    print_count(splicing_event_type::ALT_5PRIME, "Alternative 5' splice site");
    print_count(splicing_event_type::ALT_3PRIME, "Alternative 3' splice site");
    print_count(splicing_event_type::INTRON_RETENTION, "Intron retention");
    print_count(splicing_event_type::ALT_FIRST_EXON, "Alternative first exon");
    print_count(splicing_event_type::ALT_LAST_EXON, "Alternative last exon");

    // Count unique genes with events
    std::set<std::string> genes_with_events;
    for (const auto& event : events) {
        genes_with_events.insert(event.gene_id);
    }
    out << "\nGenes with splicing events: " << genes_with_events.size() << "\n";
}