/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "transcript_matcher.hpp"
#include "utility.hpp"

#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>

// ============================================================================
// transcript_matcher implementation
// ============================================================================

transcript_matcher::transcript_matcher(grove_type& grove, const config& cfg)
    : grove_(grove), cfg_(cfg) {}

match_result transcript_matcher::match(const read_cluster& cluster) {
    match_result result;
    stats_.total_matches++;

    // Step 1: Find candidate segments via spatial query
    auto candidates = find_candidate_segments(cluster);

    if (candidates.empty()) {
        result.type = match_result::match_type::INTERGENIC;
        stats_.intergenic_matches++;
        return result;
    }

    // Step 2-3: For each candidate, get exon chain and score junction match
    double best_score = 0.0;
    key_ptr best_segment = nullptr;
    std::string best_transcript_id;
    std::vector<key_ptr> best_exon_chain;

    for (key_ptr seg_key : candidates) {
        if (!is_segment(seg_key->get_data())) continue;

        const auto& seg = get_segment(seg_key->get_data());

        // Try each transcript associated with this segment
        for (const auto& transcript_id : seg.transcript_ids) {
            auto exon_chain = get_exon_chain(seg_key, transcript_id);
            if (exon_chain.empty()) continue;

            double score = score_junction_match(cluster, exon_chain);

            if (score > best_score) {
                best_score = score;
                best_segment = seg_key;
                best_transcript_id = transcript_id;
                best_exon_chain = exon_chain;
            }
        }
    }

    // Step 4: Classify match
    if (best_segment == nullptr) {
        result.type = match_result::match_type::INTERGENIC;
        stats_.intergenic_matches++;
        return result;
    }

    result.matched_segments.push_back(best_segment);
    result.matched_transcript_ids.push_back(best_transcript_id);
    result.junction_match_score = best_score;
    result.total_junctions = static_cast<int>(cluster.consensus_junctions.size());
    result.matching_junctions = static_cast<int>(best_score * result.total_junctions + 0.5);

    // Classify based on score
    result.type = classify_match(cluster, best_exon_chain, best_score);

    // Find novel junctions if applicable
    if (cfg_.track_novel && result.is_novel()) {
        result.novel_junctions = find_novel_junctions(cluster, best_exon_chain);
    }

    // Check for ambiguous matches (multiple transcripts with same best score)
    int equal_matches = 0;
    for (key_ptr seg_key : candidates) {
        if (!is_segment(seg_key->get_data())) continue;
        const auto& seg = get_segment(seg_key->get_data());

        for (const auto& transcript_id : seg.transcript_ids) {
            if (seg_key == best_segment && transcript_id == best_transcript_id) continue;

            auto exon_chain = get_exon_chain(seg_key, transcript_id);
            if (exon_chain.empty()) continue;

            double score = score_junction_match(cluster, exon_chain);
            if (std::abs(score - best_score) < 0.001) {
                equal_matches++;
                result.matched_segments.push_back(seg_key);
                result.matched_transcript_ids.push_back(transcript_id);
            }
        }
    }

    if (equal_matches > 0) {
        result.type = match_result::match_type::AMBIGUOUS;
        stats_.ambiguous_matches++;
    } else {
        // Update stats based on type
        switch (result.type) {
            case match_result::match_type::EXACT:
                stats_.exact_matches++;
                break;
            case match_result::match_type::COMPATIBLE:
                stats_.compatible_matches++;
                break;
            case match_result::match_type::NOVEL_JUNCTION:
                stats_.novel_junction_matches++;
                break;
            case match_result::match_type::NOVEL_EXON:
                stats_.novel_exon_matches++;
                break;
            default:
                break;
        }
    }

    return result;
}

std::vector<match_result> transcript_matcher::match_batch(
    const std::vector<read_cluster>& clusters) {

    std::vector<match_result> results;
    results.reserve(clusters.size());

    for (const auto& cluster : clusters) {
        results.push_back(match(cluster));
    }

    return results;
}

void transcript_matcher::update_grove(const read_cluster& cluster,
                                       const match_result& result) {
    if (result.has_match()) {
        // Add read support to matched segments
        for (key_ptr seg_key : result.matched_segments) {
            if (is_segment(seg_key->get_data())) {
                auto& seg = get_segment(seg_key->get_data());
                // Add representative read ID
                if (const auto* rep = cluster.get_representative()) {
                    seg.add_read_support(rep->read_id);
                }
                // Or add cluster ID
                // seg.add_read_support(cluster.cluster_id);
                stats_.segments_updated++;
            }
        }
    } else if (result.is_novel()) {
        // Create discovered segment for novel patterns
        key_ptr new_seg = create_discovered_segment(cluster);
        if (new_seg != nullptr) {
            stats_.segments_created++;
        }
    }
}

std::vector<key_ptr> transcript_matcher::find_candidate_segments(
    const read_cluster& cluster) {

    std::vector<key_ptr> candidates;

    gdt::genomic_coordinate query = cluster.get_coordinate();
    auto result = grove_.intersect(query, cluster.seqid);

    for (auto* key : result.get_keys()) {
        if (!is_segment(key->get_data())) continue;

        // Check minimum overlap
        const auto& seg_coord = key->get_value();
        size_t overlap_start = std::max(query.get_start(), seg_coord.get_start());
        size_t overlap_end = std::min(query.get_end(), seg_coord.get_end());

        if (overlap_end > overlap_start) {
            size_t overlap = overlap_end - overlap_start;
            if (static_cast<int>(overlap) >= cfg_.min_overlap_bp) {
                candidates.push_back(key);
            }
        }
    }

    return candidates;
}

std::vector<key_ptr> transcript_matcher::get_exon_chain(
    key_ptr segment, const std::string& transcript_id) {

    std::vector<key_ptr> chain;

    // Get SEGMENT_TO_EXON edge to find first exon
    auto first_exons = grove_.get_neighbors_if(segment,
        [&transcript_id](const edge_metadata& e) {
            return e.type == edge_metadata::edge_type::SEGMENT_TO_EXON &&
                   e.id == transcript_id;
        });

    if (first_exons.empty()) return chain;

    // Traverse EXON_TO_EXON edges
    key_ptr current = first_exons[0];
    while (current != nullptr) {
        chain.push_back(current);

        auto next_exons = grove_.get_neighbors_if(current,
            [&transcript_id](const edge_metadata& e) {
                return e.type == edge_metadata::edge_type::EXON_TO_EXON &&
                       e.id == transcript_id;
            });

        current = next_exons.empty() ? nullptr : next_exons[0];
    }

    return chain;
}

double transcript_matcher::score_junction_match(
    const read_cluster& cluster,
    const std::vector<key_ptr>& exon_chain) {

    if (cluster.consensus_junctions.empty()) {
        // Single-exon read - check if it fits within one exon
        // or spans appropriately
        return exon_chain.size() >= 1 ? 1.0 : 0.0;
    }

    // Extract reference junctions from exon chain
    auto ref_junctions = extract_reference_junctions(exon_chain);

    if (ref_junctions.empty()) {
        return 0.0;
    }

    // Count matching junctions
    int matches = 0;
    for (const auto& read_junc : cluster.consensus_junctions) {
        for (const auto& ref_junc : ref_junctions) {
            if (read_junc.matches(ref_junc, cfg_.junction_tolerance)) {
                matches++;
                break;
            }
        }
    }

    return static_cast<double>(matches) /
           static_cast<double>(cluster.consensus_junctions.size());
}

match_result::match_type transcript_matcher::classify_match(
    const read_cluster& cluster,
    const std::vector<key_ptr>& best_chain,
    double best_score) {

    // Perfect junction match
    if (best_score >= 0.999) {
        return match_result::match_type::EXACT;
    }

    // Good but not perfect - check for novel junctions
    if (best_score >= cfg_.min_junction_score) {
        auto novel = find_novel_junctions(cluster, best_chain);
        if (!novel.empty()) {
            return match_result::match_type::NOVEL_JUNCTION;
        }
        return match_result::match_type::COMPATIBLE;
    }

    // Lower score - check if boundaries extend beyond known exons
    if (best_score > 0.0) {
        // Check if cluster extends beyond exon boundaries
        if (!best_chain.empty()) {
            const auto& first_exon = best_chain.front()->get_value();
            const auto& last_exon = best_chain.back()->get_value();

            if (cluster.start < first_exon.get_start() ||
                cluster.end > last_exon.get_end()) {
                return match_result::match_type::NOVEL_EXON;
            }
        }
        return match_result::match_type::COMPATIBLE;
    }

    return match_result::match_type::INTERGENIC;
}

std::vector<splice_junction> transcript_matcher::find_novel_junctions(
    const read_cluster& cluster,
    const std::vector<key_ptr>& exon_chain) {

    std::vector<splice_junction> novel;
    auto ref_junctions = extract_reference_junctions(exon_chain);

    for (const auto& read_junc : cluster.consensus_junctions) {
        bool found = false;
        for (const auto& ref_junc : ref_junctions) {
            if (read_junc.matches(ref_junc, cfg_.junction_tolerance)) {
                found = true;
                break;
            }
        }
        if (!found) {
            novel.push_back(read_junc);
        }
    }

    return novel;
}

key_ptr transcript_matcher::create_discovered_segment(const read_cluster& cluster) {
    // Build segment feature for discovered pattern
    segment_feature new_segment(segment_feature::source_type::DISCOVERED);

    std::ostringstream id_ss;
    id_ss << "DISC_" << cluster.seqid << "_" << cluster.start << "_" << cluster.end;
    new_segment.id = id_ss.str();

    std::ostringstream coord_ss;
    coord_ss << cluster.seqid << ":" << cluster.strand << ":"
             << cluster.start << "-" << cluster.end;
    new_segment.coordinate = coord_ss.str();

    new_segment.exon_count = static_cast<int>(cluster.consensus_junctions.size() + 1);

    // Add read support
    if (const auto* rep = cluster.get_representative()) {
        new_segment.add_read_support(rep->read_id);
    }
    new_segment.read_coverage = cluster.read_count();

    // Insert into grove
    gdt::genomic_coordinate coord(cluster.strand, cluster.start, cluster.end);
    genomic_feature feature = new_segment;

    key_ptr seg_key = grove_.insert_data(cluster.seqid, coord, feature);

    return seg_key;
}

std::vector<splice_junction> transcript_matcher::extract_reference_junctions(
    const std::vector<key_ptr>& exon_chain) {

    std::vector<splice_junction> junctions;

    for (size_t i = 0; i + 1 < exon_chain.size(); ++i) {
        const auto& exon_coord = exon_chain[i]->get_value();
        const auto& next_coord = exon_chain[i + 1]->get_value();

        // Junction is at end of current exon to start of next exon
        junctions.emplace_back(exon_coord.get_end(), next_coord.get_start());
    }

    return junctions;
}

// ============================================================================
// Output methods
// ============================================================================

namespace {
    std::string match_type_to_string(match_result::match_type type) {
        switch (type) {
            case match_result::match_type::EXACT: return "EXACT";
            case match_result::match_type::COMPATIBLE: return "COMPATIBLE";
            case match_result::match_type::NOVEL_JUNCTION: return "NOVEL_JUNCTION";
            case match_result::match_type::NOVEL_EXON: return "NOVEL_EXON";
            case match_result::match_type::INTERGENIC: return "INTERGENIC";
            case match_result::match_type::AMBIGUOUS: return "AMBIGUOUS";
            default: return "UNKNOWN";
        }
    }

    std::string join_strings(const std::vector<std::string>& strs, const std::string& sep) {
        if (strs.empty()) return "";
        std::ostringstream ss;
        for (size_t i = 0; i < strs.size(); ++i) {
            if (i > 0) ss << sep;
            ss << strs[i];
        }
        return ss.str();
    }

    std::string junctions_to_string(const std::vector<splice_junction>& junctions) {
        if (junctions.empty()) return "";
        std::ostringstream ss;
        for (size_t i = 0; i < junctions.size(); ++i) {
            if (i > 0) ss << ",";
            ss << junctions[i].donor << "-" << junctions[i].acceptor;
        }
        return ss.str();
    }
}

void transcript_matcher::write_results(const std::string& filepath,
                                        const std::vector<read_cluster>& clusters,
                                        const std::vector<match_result>& results) {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        logging::error("Could not open output file: " + filepath);
        return;
    }

    // Header
    out << "cluster_id\t"
        << "seqid\t"
        << "strand\t"
        << "start\t"
        << "end\t"
        << "read_count\t"
        << "exon_count\t"
        << "junctions\t"
        << "match_type\t"
        << "junction_score\t"
        << "matching_junctions\t"
        << "total_junctions\t"
        << "transcript_ids\t"
        << "novel_junctions\t"
        << "representative_read\n";

    for (size_t i = 0; i < clusters.size(); ++i) {
        const auto& cluster = clusters[i];
        const auto& result = results[i];

        const auto* rep = cluster.get_representative();
        std::string rep_id = rep ? rep->read_id : "";

        out << cluster.cluster_id << "\t"
            << cluster.seqid << "\t"
            << cluster.strand << "\t"
            << cluster.start << "\t"
            << cluster.end << "\t"
            << cluster.read_count() << "\t"
            << (cluster.consensus_junctions.size() + 1) << "\t"
            << junctions_to_string(cluster.consensus_junctions) << "\t"
            << match_type_to_string(result.type) << "\t"
            << result.junction_match_score << "\t"
            << result.matching_junctions << "\t"
            << result.total_junctions << "\t"
            << join_strings(result.matched_transcript_ids, ",") << "\t"
            << junctions_to_string(result.novel_junctions) << "\t"
            << rep_id << "\n";
    }

    out.close();
    logging::info("Wrote results to: " + filepath);
}

void transcript_matcher::write_summary(const std::string& filepath) const {
    std::ofstream out(filepath);
    if (!out.is_open()) {
        logging::error("Could not open summary file: " + filepath);
        return;
    }

    out << "# Atroplex Transcript Matching Summary\n\n";

    out << "## Match Statistics\n";
    out << "Total clusters matched: " << stats_.total_matches << "\n";
    out << "Exact matches: " << stats_.exact_matches << "\n";
    out << "Compatible matches: " << stats_.compatible_matches << "\n";
    out << "Novel junction: " << stats_.novel_junction_matches << "\n";
    out << "Novel exon: " << stats_.novel_exon_matches << "\n";
    out << "Intergenic: " << stats_.intergenic_matches << "\n";
    out << "Ambiguous: " << stats_.ambiguous_matches << "\n";
    out << "\n";

    out << "## Grove Updates\n";
    out << "Segments updated with read support: " << stats_.segments_updated << "\n";
    out << "Novel segments created: " << stats_.segments_created << "\n";

    // Calculate percentages
    if (stats_.total_matches > 0) {
        out << "\n## Match Type Distribution\n";
        auto pct = [this](size_t count) {
            return 100.0 * static_cast<double>(count) / static_cast<double>(stats_.total_matches);
        };
        out << "Exact: " << pct(stats_.exact_matches) << "%\n";
        out << "Compatible: " << pct(stats_.compatible_matches) << "%\n";
        out << "Novel junction: " << pct(stats_.novel_junction_matches) << "%\n";
        out << "Novel exon: " << pct(stats_.novel_exon_matches) << "%\n";
        out << "Intergenic: " << pct(stats_.intergenic_matches) << "%\n";
        out << "Ambiguous: " << pct(stats_.ambiguous_matches) << "%\n";
    }

    out.close();
    logging::info("Wrote summary to: " + filepath);
}