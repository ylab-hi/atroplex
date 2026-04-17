/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "segment_builder.hpp"

// Absorption rule thresholds (see absorption_rules.txt)
static constexpr size_t TERMINAL_TOLERANCE_BP = 50;
// fuzzy_tolerance removed — now passed as parameter through the call chain
static constexpr size_t MAX_ISM_MISSING_EXONS = 2;

static size_t abs_diff(size_t a, size_t b) {
    return (a > b) ? a - b : b - a;
}

// ========================================================================
// Absorption rules — see absorption_rules.txt for full documentation.
//
// Rule 0: FSM            — identical exon coordinates → Merge
// Rule 1: 5' ISM         — contiguous subset at 5' end → Keep
// Rule 2: 3' ISM         — contiguous subset at 3' end, 1-2 missing → Absorb
// Rule 3: 3' degradation — contiguous subset at 3' end, 3+ missing → Drop(ref)/Keep(sample)
// Rule 4: Internal frag  — contiguous subset, both ends missing → Drop(ref)/Keep(sample)
// Rule 5: Terminal var    — same intron chain, TSS/TES <50bp → Absorb
// Rule 6: Mono-exon gene — overlaps gene, no intron crossing → Drop
// Rule 7: Mono-exon IR   — spans exon-intron-exon → Keep
// Rule 8: Mono-exon IG   — intergenic, no gene overlap → Drop
//
// Matching: pointer identity first, fuzzy fallback (≤5bp) for Rules 0-4.
// Annotations are processed before samples. Default: create segment.
//
// Execution order:
//   Step 1: Rule 0        (FSM — pointer, then fuzzy)
//   Step 2: Rule 5        (terminal variant)
//   Step 3: Rules 6/7/8   (mono-exon)
//   Step 4: Rules 1/2/3/4 (subsequence — pointer, then fuzzy)
//   Step 5: Create segment + edges
//   Step 6: Reverse absorption (Steps 1, 2, 4 on existing segments)
// ========================================================================

void segment_builder::create_segment(
    grove_type& grove,
    std::mutex& grove_mutex,
    const std::string& transcript_id,
    const std::string& seqid,
    char strand,
    size_t span_start,
    size_t span_end,
    int exon_count,
    const std::vector<gdt::genomic_coordinate>& exon_coords,
    const std::vector<key_ptr>& exon_chain,
    segment_cache_type& segment_cache,
    gene_segment_index_type& gene_index,
    std::optional<uint32_t> sample_id,
    const std::string& gff_source,
    size_t& segment_count,
    float expression_value,
    const std::string& transcript_biotype,
    bool absorb,
    size_t fuzzy_tolerance,
    build_counters& counters,
    quant_sidecar::SampleStreamWriter* sidecar_writer
) {
    std::string structure_key = make_exon_structure_key(seqid, exon_coords);

    // Step 1: Rule 0 — FSM (exact structure match → merge metadata)
    auto cached = segment_cache.find(structure_key);
    if (cached != segment_cache.end()) {
        merge_into_segment(cached->second, transcript_id, sample_id,
                          gff_source, expression_value, transcript_biotype,
                          sidecar_writer);
        counters.merged_transcripts++;
        return;
    }

    const std::string& gene_id = get_exon(exon_chain.front()->get_data()).gene_id();
    auto gene_it = absorb ? gene_index.find(gene_id) : gene_index.end();
    bool has_gene_entries = (gene_it != gene_index.end());

    // Step 2: Rule 5 — Terminal variant (same intron chain, TSS/TES <50bp → absorb)
    if (absorb && exon_chain.size() >= 2 && has_gene_entries) {
        for (const auto& entry : gene_it->second) {
            auto& candidate_seg = get_segment(entry.segment->get_data());
            if (candidate_seg.absorbed) continue;

            // Rule 5 requires equal exon counts — skip before calling into the helpers.
            if (entry.exon_chain.size() != exon_chain.size()) continue;

            if (has_same_intron_chain(exon_chain, entry.exon_chain) &&
                terminal_boundaries_within_tolerance(exon_chain, entry.exon_chain, TERMINAL_TOLERANCE_BP)) {
                    merge_into_segment(entry.segment, transcript_id, sample_id,
                                      gff_source, expression_value, transcript_biotype,
                                      sidecar_writer);
                    get_segment(entry.segment->get_data()).absorbed_count++;
                    counters.merged_transcripts++;
                    return;
                }
            }
    }

    // Step 3: Rules 6/7/8 — Mono-exon handling
    if (exon_chain.size() == 1) {
        if (!absorb) return;

        auto& mono_coord = exon_chain.front()->get_value();
        auto mono_class = classify_mono_exon(mono_coord, gene_index, gene_id);

        switch (mono_class) {
            case mono_exon_class::INTRON_RETENTION:
                break;  // Rule 7: keep — fall through to segment creation
            case mono_exon_class::GENE_OVERLAP:
                counters.discarded_transcripts++;
                return; // Rule 6: drop
            case mono_exon_class::INTERGENIC:
                counters.discarded_transcripts++;
                return; // Rule 8: drop
        }
    }

    // Step 4: Rules 0(fuzzy)/1/2/3/4 — Subsequence matching (pointer identity, then fuzzy)
    if (absorb && exon_chain.size() >= 2 && has_gene_entries) {
        {
            key_ptr best_parent = nullptr;
            size_t best_exon_count = 0;
            subsequence_type best_match = subsequence_type::NONE;

            for (const auto& entry : gene_it->second) {
                auto& candidate_seg = get_segment(entry.segment->get_data());
                if (candidate_seg.absorbed) continue;

                // Early-exit: parent must have at least as many exons as child.
                // Equal-size chains can still match via fuzzy-FSM (all boundaries
                // shifted by <=tolerance) and are handled as subsequence_type::FSM.
                if (entry.exon_chain.size() < exon_chain.size()) continue;

                // Early-exit: child's coordinate span must fit within parent's span
                // (accounting for fuzzy tolerance on the boundaries).
                const auto& parent_first = entry.exon_chain.front()->get_value();
                const auto& parent_last = entry.exon_chain.back()->get_value();
                if (span_start + fuzzy_tolerance < parent_first.get_start()) continue;
                if (span_end > parent_last.get_end() + fuzzy_tolerance) continue;

                auto match = classify_subsequence(exon_chain, entry.exon_chain);
                if (match == subsequence_type::NONE) {
                    match = fuzzy_classify_subsequence(exon_chain, entry.exon_chain, fuzzy_tolerance);
                }

                if (match == subsequence_type::NONE) continue;
                if (match == subsequence_type::ISM_5PRIME) continue; // Rule 1: keep

                // Fuzzy-FSM takes priority over any ISM/fragment match — absorb immediately.
                if (match == subsequence_type::FSM) {
                    merge_into_segment(entry.segment, transcript_id, sample_id,
                                      gff_source, expression_value, transcript_biotype,
                                      sidecar_writer);
                    get_segment(entry.segment->get_data()).absorbed_count++;
                    counters.merged_transcripts++;
                    return;
                }

                if (entry.exon_chain.size() > best_exon_count) {
                    best_parent = entry.segment;
                    best_exon_count = entry.exon_chain.size();
                    best_match = match;
                }
            }

            if (best_parent != nullptr) {
                if (best_match == subsequence_type::ISM_3PRIME) {
                    // Rule 2: always absorb
                    merge_into_segment(best_parent, transcript_id, sample_id,
                                      gff_source, expression_value, transcript_biotype,
                                      sidecar_writer);
                    get_segment(best_parent->get_data()).absorbed_count++;
                    counters.merged_transcripts++;
                    return;
                }

                // Rules 3/4: drop vs ref, keep vs sample
                if (is_parent_annotation(best_parent)) {
                    counters.discarded_transcripts++;
                    return; // Drop
                }
                // vs sample: keep — fall through to creation
            }
        }
    }

    // Step 5: Create new segment
    gdt::genomic_coordinate segment_coord(strand, span_start, span_end);

    segment_feature new_segment;
    new_segment.segment_index = segment_count;
    uint32_t tx_id = transcript_registry::instance().intern(transcript_id);
    new_segment.transcript_ids.insert(tx_id);
    if (!transcript_biotype.empty()) {
        new_segment.transcript_biotypes[tx_id] = transcript_biotype;
    }
    new_segment.exon_count = exon_count;

    const auto& first_exon_data = get_exon(exon_chain.front()->get_data());
    new_segment.gene_idx = first_exon_data.gene_idx;

    if (sample_id.has_value()) {
        new_segment.add_sample(*sample_id);
    }
    if (!gff_source.empty()) {
        new_segment.add_source(gff_source);
    }

    genomic_feature feature = new_segment;
    key_ptr seg_key;
    {
        std::lock_guard<std::mutex> lock(grove_mutex);
        seg_key = grove.insert_data(seqid, segment_coord, feature);

        edge_metadata seg_edge(segment_count, edge_metadata::edge_type::SEGMENT_TO_EXON);
        grove.add_edge(seg_key, exon_chain.front(), seg_edge);

        for (size_t i = 0; i + 1 < exon_chain.size(); ++i) {
            edge_metadata exon_edge(segment_count, edge_metadata::edge_type::EXON_TO_EXON);
            grove.add_edge(exon_chain[i], exon_chain[i + 1], exon_edge);
        }
    }
    segment_count++;

    // Write expression value to sidecar (if writer armed and value present).
    // segment_index = segment_count - 1 (just assigned above).
    if (sidecar_writer && expression_value >= 0.0f) {
        sidecar_writer->append(segment_count - 1, expression_value);
    }

    segment_cache[structure_key] = seg_key;

    // Step 5b: Register in gene segment index
    gene_index[gene_id].push_back({seg_key, exon_chain, structure_key});

    // Step 6: Reverse absorption
    if (absorb) {
        try_reverse_absorption(gene_index, gene_id, seg_key, exon_chain, segment_cache, fuzzy_tolerance);
    }
}

// ========================================================================
// Subsequence classification
// ========================================================================

subsequence_type segment_builder::classify_subsequence(
    const std::vector<key_ptr>& sub,
    const std::vector<key_ptr>& parent
) {
    if (sub.size() < 2 || sub.size() >= parent.size()) return subsequence_type::NONE;

    for (size_t start = 0; start <= parent.size() - sub.size(); ++start) {
        if (parent[start] == sub[0]) {
            bool all_match = true;
            for (size_t j = 1; j < sub.size(); ++j) {
                if (parent[start + j] != sub[j]) {
                    all_match = false;
                    break;
                }
            }
            if (all_match) {
                bool shares_first = (start == 0);
                bool shares_last  = (start + sub.size() == parent.size());
                if (shares_first) return subsequence_type::ISM_5PRIME;
                if (shares_last) {
                    size_t missing_from_5prime = start;
                    return (missing_from_5prime <= MAX_ISM_MISSING_EXONS)
                        ? subsequence_type::ISM_3PRIME
                        : subsequence_type::DEGRADATION_3PRIME;
                }
                return subsequence_type::INTERNAL_FRAGMENT;
            }
        }
    }
    return subsequence_type::NONE;
}

subsequence_type segment_builder::fuzzy_classify_subsequence(
    const std::vector<key_ptr>& sub,
    const std::vector<key_ptr>& parent,
    size_t tolerance_bp
) {
    if (sub.size() < 2 || sub.size() > parent.size()) return subsequence_type::NONE;

    for (size_t start = 0; start <= parent.size() - sub.size(); ++start) {
        bool all_match = true;
        for (size_t j = 0; j < sub.size(); ++j) {
            auto& coord_s = sub[j]->get_value();
            auto& coord_p = parent[start + j]->get_value();
            if (abs_diff(coord_s.get_start(), coord_p.get_start()) > tolerance_bp ||
                abs_diff(coord_s.get_end(), coord_p.get_end()) > tolerance_bp) {
                all_match = false;
                break;
            }
        }
        if (all_match) {
            // Fuzzy FSM: same exon count, all boundaries within tolerance
            if (sub.size() == parent.size()) return subsequence_type::FSM;
            bool shares_first = (start == 0);
            bool shares_last  = (start + sub.size() == parent.size());
            if (shares_first) return subsequence_type::ISM_5PRIME;
            if (shares_last) {
                size_t missing_from_5prime = start;
                return (missing_from_5prime <= MAX_ISM_MISSING_EXONS)
                    ? subsequence_type::ISM_3PRIME
                    : subsequence_type::DEGRADATION_3PRIME;
            }
            return subsequence_type::INTERNAL_FRAGMENT;
        }
    }
    return subsequence_type::NONE;
}

// ========================================================================
// Merge and reverse absorption
// ========================================================================

void segment_builder::merge_into_segment(
    key_ptr target_seg,
    const std::string& transcript_id,
    std::optional<uint32_t> sample_id,
    const std::string& gff_source,
    float expression_value,
    const std::string& transcript_biotype,
    quant_sidecar::SampleStreamWriter* sidecar_writer
) {
    auto& seg = get_segment(target_seg->get_data());
    uint32_t tx_id = transcript_registry::instance().intern(transcript_id);
    seg.transcript_ids.insert(tx_id);
    if (!transcript_biotype.empty()) {
        seg.transcript_biotypes[tx_id] = transcript_biotype;
    }
    if (sample_id.has_value()) {
        seg.add_sample(*sample_id);
    }
    if (!gff_source.empty()) {
        seg.add_source(gff_source);
    }
    // Sidecar write on merge: the original 3a wiring only wrote inside
    // create_segment's "Step 5: Create new segment" branch, so any
    // transcript that merged into an existing segment (Rule 0 exact,
    // Rule 5 terminal variant, fuzzy-FSM, Rule 2 ISM_3PRIME, or
    // try_reverse_absorption) dropped its counts. With annotations
    // processed first, every segment is created without expression and
    // every sample transcript then merges — leaving the .qtx empty.
    if (sidecar_writer && expression_value >= 0.0f) {
        sidecar_writer->append(seg.segment_index, expression_value);
    }
}

void segment_builder::try_reverse_absorption(
    gene_segment_index_type& gene_index,
    const std::string& gene_id,
    key_ptr new_seg,
    const std::vector<key_ptr>& new_exon_chain,
    segment_cache_type& segment_cache,
    size_t fuzzy_tolerance
) {
    auto gene_it = gene_index.find(gene_id);
    if (gene_it == gene_index.end()) return;

    auto& parent_seg = get_segment(new_seg->get_data());
    bool new_seg_is_ref = is_parent_annotation(new_seg);

    auto absorb_into_parent = [&](auto& candidate_seg, auto& entry) {
        for (const auto& tx : candidate_seg.transcript_ids)
            parent_seg.transcript_ids.insert(tx);
        for (const auto& [tx, bt] : candidate_seg.transcript_biotypes)
            parent_seg.transcript_biotypes[tx] = bt;
        parent_seg.sample_idx.merge(candidate_seg.sample_idx);
        parent_seg.merge_sources(candidate_seg.sources);
        parent_seg.absorbed_count += candidate_seg.absorbed_count + 1;
        candidate_seg.absorbed = true;
        // Record the parent link so any .qtx records already written to
        // disk under the candidate's segment_index get remapped to the
        // parent at merge_to_qtx time instead of being dropped by the
        // tombstone filter. Issue #34.
        candidate_seg.absorbed_into_idx = parent_seg.segment_index;
        segment_cache.erase(entry.structure_key);
    };

    auto tombstone_candidate = [&](auto& candidate_seg, auto& entry) {
        candidate_seg.absorbed = true;
        segment_cache.erase(entry.structure_key);
    };

    // Precompute new segment's coordinate span for span filter
    const auto& new_first = new_exon_chain.front()->get_value();
    const auto& new_last = new_exon_chain.back()->get_value();
    size_t new_start = new_first.get_start();
    size_t new_end = new_last.get_end();

    for (auto& entry : gene_it->second) {
        if (entry.segment == new_seg) continue;

        auto& candidate_seg = get_segment(entry.segment->get_data());
        if (candidate_seg.absorbed) continue;

        // Rule 5: Terminal variant (requires equal exon count)
        if (entry.exon_chain.size() == new_exon_chain.size() &&
            has_same_intron_chain(entry.exon_chain, new_exon_chain) &&
            terminal_boundaries_within_tolerance(entry.exon_chain, new_exon_chain, TERMINAL_TOLERANCE_BP)) {
            absorb_into_parent(candidate_seg, entry);
            continue;
        }

        // Early-exit: candidate (child) must have at most as many exons as new
        // segment (parent). Equal-size candidates can still match via fuzzy-FSM.
        if (entry.exon_chain.size() > new_exon_chain.size()) continue;

        // Early-exit: candidate's span must fit within the new segment's span
        // (accounting for fuzzy tolerance).
        const auto& cand_first = entry.exon_chain.front()->get_value();
        const auto& cand_last = entry.exon_chain.back()->get_value();
        if (cand_first.get_start() + fuzzy_tolerance < new_start) continue;
        if (cand_last.get_end() > new_end + fuzzy_tolerance) continue;

        // Rules 0(fuzzy)/1/2/3/4: Subsequence (pointer, then fuzzy)
        auto match = classify_subsequence(entry.exon_chain, new_exon_chain);
        if (match == subsequence_type::NONE) {
            match = fuzzy_classify_subsequence(entry.exon_chain, new_exon_chain, fuzzy_tolerance);
        }

        if (match == subsequence_type::NONE) continue;
        if (match == subsequence_type::ISM_5PRIME) continue; // Rule 1: keep

        if (match == subsequence_type::FSM || match == subsequence_type::ISM_3PRIME) {
            // Fuzzy-FSM or Rule 2: always absorb candidate into new segment
            absorb_into_parent(candidate_seg, entry);
        } else {
            // Rules 3/4: drop vs ref, keep vs sample
            if (new_seg_is_ref) {
                tombstone_candidate(candidate_seg, entry);
            }
        }
    }
}

// ========================================================================
// Rule helpers
// ========================================================================

bool segment_builder::has_same_intron_chain(
    const std::vector<key_ptr>& chain_a,
    const std::vector<key_ptr>& chain_b
) {
    if (chain_a.size() != chain_b.size()) return false;
    if (chain_a.size() < 2) return false;

    for (size_t i = 1; i + 1 < chain_a.size(); ++i) {
        if (chain_a[i] != chain_b[i]) return false;
    }

    if (chain_a.front()->get_value().get_end() != chain_b.front()->get_value().get_end())
        return false;
    if (chain_a.back()->get_value().get_start() != chain_b.back()->get_value().get_start())
        return false;

    return true;
}

bool segment_builder::terminal_boundaries_within_tolerance(
    const std::vector<key_ptr>& chain_a,
    const std::vector<key_ptr>& chain_b,
    size_t tolerance_bp
) {
    size_t start_diff = abs_diff(
        chain_a.front()->get_value().get_start(),
        chain_b.front()->get_value().get_start());
    if (start_diff > tolerance_bp) return false;

    size_t end_diff = abs_diff(
        chain_a.back()->get_value().get_end(),
        chain_b.back()->get_value().get_end());
    if (end_diff > tolerance_bp) return false;

    return true;
}

segment_builder::mono_exon_class segment_builder::classify_mono_exon(
    const gdt::genomic_coordinate& mono_coord,
    const gene_segment_index_type& gene_index,
    const std::string& gene_id
) {
    auto gene_it = gene_index.find(gene_id);
    if (gene_it == gene_index.end()) return mono_exon_class::INTERGENIC;

    bool overlaps_gene = false;

    for (const auto& entry : gene_it->second) {
        auto& candidate_seg = get_segment(entry.segment->get_data());
        if (candidate_seg.absorbed) continue;
        if (entry.exon_chain.size() < 2) continue;

        // Check if mono-exon spans from one exon across an intron into the next
        for (size_t i = 0; i + 1 < entry.exon_chain.size(); ++i) {
            auto& exon_coord = entry.exon_chain[i]->get_value();
            auto& next_coord = entry.exon_chain[i + 1]->get_value();

            if (mono_coord.get_start() < exon_coord.get_end() &&
                mono_coord.get_end() > next_coord.get_start()) {
                return mono_exon_class::INTRON_RETENTION;
            }
        }

        // Check if mono-exon overlaps any exon of this segment
        if (!overlaps_gene) {
            for (const auto& exon_key : entry.exon_chain) {
                auto& exon_coord = exon_key->get_value();
                if (mono_coord.get_start() < exon_coord.get_end() &&
                    mono_coord.get_end() > exon_coord.get_start()) {
                    overlaps_gene = true;
                    break;
                }
            }
        }
    }

    return overlaps_gene ? mono_exon_class::GENE_OVERLAP
                         : mono_exon_class::INTERGENIC;
}

bool segment_builder::is_annotation_sample(std::optional<uint32_t> sample_id) {
    if (!sample_id.has_value()) return false;
    auto& registry = sample_registry::instance();
    const auto& info = registry.get(*sample_id);
    return info.type == "annotation" || info.is_annotation();
}

bool segment_builder::is_parent_annotation(key_ptr parent_seg) {
    auto& seg = get_segment(parent_seg->get_data());
    auto& registry = sample_registry::instance();
    bool is_ref = false;
    seg.sample_idx.for_each([&](uint32_t sid) {
        if (!is_ref) {
            const auto& info = registry.get(sid);
            if (info.type == "annotation" || info.is_annotation()) {
                is_ref = true;
            }
        }
    });
    return is_ref;
}

std::string segment_builder::make_exon_structure_key(
    const std::string& seqid,
    const std::vector<gdt::genomic_coordinate>& exon_coords
) {
    if (exon_coords.empty()) {
        return "";
    }

    std::string key;
    key.reserve(seqid.size() + 3 + exon_coords.size() * 24);

    key += seqid;
    key += ':';
    key += exon_coords.front().get_strand();
    key += ':';

    for (size_t i = 0; i < exon_coords.size(); ++i) {
        if (i > 0) {
            key += ',';
        }
        key += std::to_string(exon_coords[i].get_start());
        key += '-';
        key += std::to_string(exon_coords[i].get_end());
    }

    return key;
}