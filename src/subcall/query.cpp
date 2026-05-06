/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/query.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>

#include <genogrove/io/gff_reader.hpp>

#include "sample_info.hpp"
#include "utility.hpp"

namespace gio = genogrove::io;

namespace {

/// Render a sample's declared expression_type as the short suffix used
/// in output column headers. Matches analysis_report's copy of the same
/// helper — any change here must be mirrored there.
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

} // anonymous namespace

namespace subcall {

// ============================================================================
// Contrast parsing
// ============================================================================

query_contrast query_contrast::parse(const std::string& spec) {
    query_contrast c;
    size_t pos = spec.find(':');
    if (pos == std::string::npos) {
        throw std::invalid_argument("Invalid contrast format: " + spec +
                                    " (expected group_a:group_b)");
    }
    c.group_a = spec.substr(0, pos);
    c.group_b = spec.substr(pos + 1);

    if (c.group_a.empty() || c.group_b.empty()) {
        throw std::invalid_argument("Empty group in contrast: " + spec);
    }
    return c;
}

// ============================================================================
// CLI options
// ============================================================================

cxxopts::Options query::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex query",
        "Classify transcripts against the pan-transcriptome index");

    options.add_options("Input/Output")
        ("i,input", "Input GTF/GFF with transcripts to classify",
            cxxopts::value<std::string>())
        ;

    options.add_options("Classification")
        ("junction-tolerance", "Max bp difference for junction match",
            cxxopts::value<int>()->default_value("5"))
        ("min-junction-score", "Minimum junction match score (0-1)",
            cxxopts::value<double>()->default_value("0.8"))
        ("min-overlap", "Minimum overlap with reference segment (bp)",
            cxxopts::value<int>()->default_value("50"))
        ;

    options.add_options("Differential usage")
        ("c,contrast",
            "Group contrast for DTU (format: group_a:group_b, repeatable). "
            "group_a and group_b must be values of sample_info::group "
            "(manifest 'group' column or auto-inferred from '_repNN' suffix). "
            "Requires a .qtx sidecar alongside the grove.",
            cxxopts::value<std::vector<std::string>>())
        ("fdr", "FDR threshold for significance",
            cxxopts::value<double>()->default_value("0.05"))
        ;

    add_common_options(options);

    return options;
}

void query::validate(const cxxopts::ParseResult& args) {
    // Must have a grove source
    if (!args.count("genogrove") && !args.count("build-from") && !args.count("manifest")) {
        throw std::runtime_error("Must provide -g/--genogrove, -m/--manifest, or -b/--build-from");
    }

    // Must have input transcripts
    if (!args.count("input")) {
        throw std::runtime_error("Must provide -i/--input with transcripts to classify");
    }

    std::string input = args["input"].as<std::string>();
    if (!std::filesystem::exists(input)) {
        throw std::runtime_error("Input file not found: " + input);
    }
}

// ============================================================================
// Main execution
// ============================================================================

void query::execute(const cxxopts::ParseResult& args) {
    // DTU requires per-sample expression from the .qtx sidecar; abort
    // early if the user asked for it but the sidecar didn't load.
    if (args.count("contrast") && qtx_reader_ptr() == nullptr) {
        throw std::runtime_error(
            "--contrast requires a .qtx sidecar alongside the grove, but none "
            "was found. Rebuild with expression data or drop --contrast to "
            "run structural classification only."
        );
    }

    std::string input_path = args["input"].as<std::string>();
    logging::info("Classifying transcripts from: " + input_path);

    // Step 1: Classify input transcripts against the grove
    auto results = classify_transcripts(input_path);

    // Log classification summary
    std::map<structural_category, size_t> category_counts;
    for (const auto& r : results) {
        category_counts[r.category]++;
    }

    logging::info("Classification complete: " + std::to_string(results.size()) + " transcripts");
    for (const auto& [cat, count] : category_counts) {
        logging::info("  " + to_string(cat) + ": " + std::to_string(count));
    }

    // Step 2: Output classification
    auto out_dir = resolve_output_dir(args, input_path);
    std::string basename = std::filesystem::path(input_path).stem().string();

    std::string class_path = (out_dir / (basename + ".query.tsv")).string();
    write_classification(class_path, results);

    // Step 3: Optional DTU analysis (one output file per contrast)
    std::vector<std::pair<query_contrast, std::vector<dtu_result>>> all_dtu;

    if (args.count("contrast")) {
        double fdr_threshold = args["fdr"].as<double>();
        auto contrast_strs = args["contrast"].as<std::vector<std::string>>();

        for (const auto& spec : contrast_strs) {
            auto contrast = query_contrast::parse(spec);
            logging::info("Running DTU: " + contrast.group_a + " vs " + contrast.group_b);

            auto dtu = run_dtu(contrast, results);
            apply_fdr_correction(dtu, fdr_threshold);

            size_t sig = std::count_if(dtu.begin(), dtu.end(),
                [](const dtu_result& r) { return r.significant; });
            logging::info("  Significant: " + std::to_string(sig) +
                          "/" + std::to_string(dtu.size()));

            std::string dtu_path = (out_dir / (basename + "." +
                contrast.group_a + "_vs_" + contrast.group_b + ".dtu.tsv")).string();
            write_dtu_results(dtu_path, contrast, dtu);

            all_dtu.emplace_back(contrast, std::move(dtu));
        }
    }

    // Step 4: Summary
    std::string summary_path = (out_dir / (basename + ".query.summary.txt")).string();
    write_summary(summary_path, results, all_dtu);

    logging::info("Results written to: " + out_dir.string());
}

// ============================================================================
// Transcript classification
// ============================================================================

std::vector<query_result> query::classify_transcripts(const std::string& input_path) {
    std::vector<query_result> results;

    // Configure matcher
    transcript_matcher::config match_cfg;
    match_cfg.junction_tolerance = 5;
    match_cfg.min_junction_score = 0.8;
    match_cfg.min_overlap_bp = 50;

    transcript_matcher matcher(*grove, match_cfg);

    // Read input GTF/GFF gene-by-gene and classify each transcript
    gio::gff_reader reader(input_path);

    // Group entries by transcript
    std::unordered_map<std::string, std::vector<gio::gff_entry>> transcripts;
    for (const auto& entry : reader) {
        if (entry.type != "exon") continue;
        auto tx_id = entry.get_transcript_id();
        if (tx_id.has_value()) {
            transcripts[tx_id.value()].push_back(entry);
        }
    }

    logging::info("Read " + std::to_string(transcripts.size()) +
                  " transcripts from input");

    for (auto& [tx_id, exon_entries] : transcripts) {
        // Sort exons in genomic order
        std::sort(exon_entries.begin(), exon_entries.end(),
            [](const gio::gff_entry& a, const gio::gff_entry& b) {
                return a.start < b.start;
            });

        // Skip mono-exon transcripts (consistent with build pipeline)
        if (exon_entries.size() < 2) continue;

        // Build a read_cluster from the exon chain (reuse matcher infrastructure)
        read_cluster cluster;
        cluster.seqid = normalize_chromosome(exon_entries.front().seqid);
        cluster.strand = exon_entries.front().strand.value_or('+');
        cluster.start = exon_entries.front().start;
        cluster.end = exon_entries.back().end;

        // Extract splice junctions from consecutive exons
        for (size_t i = 0; i + 1 < exon_entries.size(); ++i) {
            splice_junction junc(exon_entries[i].end, exon_entries[i + 1].start);
            cluster.consensus_junctions.push_back(junc);
        }

        cluster.cluster_id = tx_id;

        // Classify against the grove
        match_result match = matcher.match(cluster);

        // Build query result with per-sample context
        query_result qr;
        qr.transcript_id = tx_id;
        qr.category = match.category;
        qr.subcat = match.subcat;
        qr.junction_match_score = match.junction_match_score;
        qr.matching_junctions = match.matching_junctions;
        qr.total_query_junctions = match.total_query_junctions;
        qr.total_ref_junctions = match.total_ref_junctions;
        qr.known_donors = match.known_donors;
        qr.known_acceptors = match.known_acceptors;
        qr.novel_donors = match.novel_donors;
        qr.novel_acceptors = match.novel_acceptors;

        if (match.has_match()) {
            qr.gene_id = match.reference_gene.value_or(".");
            qr.gene_name = ".";

            // Collect per-sample presence and expression from matched segment
            auto* best_seg_key = match.matched_segments.front();
            const auto& seg = get_segment(best_seg_key->get_data());

            qr.gene_name = seg.gene_name();
            qr.matched_segment_index = seg.segment_index;

            seg.sample_idx.for_each([&](uint32_t sid) {
                qr.sample_ids.push_back(sid);
            });

            // Sidecar-backed per-sample expression (3a removed expression
            // from segment features; values now live in the .qtx sidecar).
            // One lookup per matched segment — bounded by sample_idx size.
            if (auto* reader = qtx_reader_ptr()) {
                auto records = reader->lookup(seg.segment_index);
                for (const auto& rec : records) {
                    qr.sample_expression[rec.sample_id] = rec.value;
                }
            }
        }

        results.push_back(std::move(qr));
    }

    return results;
}

// ============================================================================
// Output
// ============================================================================

void query::write_classification(const std::string& path,
                                  const std::vector<query_result>& results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open output file: " + path);
        return;
    }

    // Collect all sample entries for column headers
    auto& registry = sample_registry::instance();
    std::vector<uint32_t> sample_ids;
    for (size_t i = 0; i < registry.size(); ++i) {
        const auto& info = registry.get(static_cast<uint32_t>(i));
        if (info.type != "replicate") {
            sample_ids.push_back(static_cast<uint32_t>(i));
        }
    }

    // Emit per-sample expression columns only when a .qtx sidecar was
    // loaded alongside the grove. Without it every row would be `.` —
    // noisy columns that say "no data" for every value.
    const bool emit_expression_cols = qtx_reader_ptr() != nullptr;

    // Header
    out << "transcript_id\tgene_id\tgene_name\t"
        << "structural_category\tsubcategory\t"
        << "junction_match_score\tmatching_junctions\t"
        << "query_junctions\tref_junctions\t"
        << "known_donors\tknown_acceptors\t"
        << "novel_donors\tnovel_acceptors\t"
        << "n_samples";

    // Per-sample presence columns
    for (uint32_t sid : sample_ids) {
        const auto& info = registry.get(sid);
        out << "\t" << info.id << ".present";
    }

    // Per-sample expression columns (sample-typed entries only; annotations
    // never carry expression). Column label encodes the declared expr_type,
    // matching the hub / sample_stats convention in analysis_report.
    if (emit_expression_cols) {
        for (uint32_t sid : sample_ids) {
            const auto& info = registry.get(sid);
            if (info.type != "sample") continue;
            out << "\t" << info.id << "." << expr_type_label(info.expr_type);
        }
    }
    out << "\n";

    // Rows — skip unmatched transcripts (intergenic/antisense/genic with no
    // useful per-sample data; counts are in the summary file)
    for (const auto& r : results) {
        if (r.category == structural_category::INTERGENIC ||
            r.category == structural_category::ANTISENSE) continue;

        out << r.transcript_id << "\t"
            << r.gene_id << "\t"
            << r.gene_name << "\t"
            << to_string(r.category) << "\t"
            << r.subcat.to_string() << "\t"
            << r.junction_match_score << "\t"
            << r.matching_junctions << "\t"
            << r.total_query_junctions << "\t"
            << r.total_ref_junctions << "\t"
            << r.known_donors << "\t"
            << r.known_acceptors << "\t"
            << r.novel_donors << "\t"
            << r.novel_acceptors << "\t"
            << r.sample_ids.size();

        // Per-sample presence
        std::set<uint32_t> present_set(r.sample_ids.begin(), r.sample_ids.end());
        for (uint32_t sid : sample_ids) {
            out << "\t" << (present_set.count(sid) ? "1" : "0");
        }

        // Per-sample expression (sample-typed only; `.` when no value
        // for this (transcript, sample) pair in the sidecar)
        if (emit_expression_cols) {
            for (uint32_t sid : sample_ids) {
                const auto& info = registry.get(sid);
                if (info.type != "sample") continue;
                auto it = r.sample_expression.find(sid);
                if (it != r.sample_expression.end()) {
                    out << "\t" << it->second;
                } else {
                    out << "\t.";
                }
            }
        }
        out << "\n";
    }

    logging::info("Classification written to: " + path);
}

// ============================================================================
// DTU analysis
// ============================================================================

std::vector<dtu_result> query::run_dtu(
    const query_contrast& contrast,
    const std::vector<query_result>& results) {

    auto& registry = sample_registry::instance();

    // Partition samples by sample_info::group. Samples with empty or "."
    // group are skipped entirely — per the product decision that group
    // is the one canonical grouping axis (no --group-by field selector).
    std::set<uint32_t> group_a_samples, group_b_samples;
    for (size_t i = 0; i < registry.size(); ++i) {
        auto sid = static_cast<uint32_t>(i);
        const auto& info = registry.get(sid);
        if (info.type != "sample") continue;
        if (info.group.empty() || info.group == ".") continue;

        if      (info.group == contrast.group_a) group_a_samples.insert(sid);
        else if (info.group == contrast.group_b) group_b_samples.insert(sid);
    }

    if (group_a_samples.empty() || group_b_samples.empty()) {
        logging::warning("No samples found for contrast " +
                         contrast.group_a + " vs " + contrast.group_b +
                         " (check sample_info::group values)");
        return {};
    }

    logging::info("  Group " + contrast.group_a + ": " +
                  std::to_string(group_a_samples.size()) + " samples");
    logging::info("  Group " + contrast.group_b + ": " +
                  std::to_string(group_b_samples.size()) + " samples");

    // Group classified transcripts by gene
    std::map<std::string, std::vector<const query_result*>> by_gene;
    for (const auto& r : results) {
        if (r.gene_id.empty() || r.gene_id == ".") continue;
        by_gene[r.gene_id].push_back(&r);
    }

    std::vector<dtu_result> dtu_results;

    // For each gene, compute per-transcript proportions and test
    for (const auto& [gene_id, gene_transcripts] : by_gene) {
        if (gene_transcripts.size() < 2) continue;

        struct tx_expr {
            std::string transcript_id;
            std::string gene_name;
            double mean_a = 0.0;
            double mean_b = 0.0;
        };

        std::vector<tx_expr> tx_data;
        double gene_total_a = 0.0;
        double gene_total_b = 0.0;

        for (const auto* qr : gene_transcripts) {
            tx_expr te;
            te.transcript_id = qr->transcript_id;
            te.gene_name = qr->gene_name;

            double sum_a = 0.0, sum_b = 0.0;
            size_t n_a = 0, n_b = 0;

            for (const auto& [sid, expr] : qr->sample_expression) {
                if (group_a_samples.count(sid)) { sum_a += expr; n_a++; }
                else if (group_b_samples.count(sid)) { sum_b += expr; n_b++; }
            }

            te.mean_a = (n_a > 0) ? sum_a / n_a : 0.0;
            te.mean_b = (n_b > 0) ? sum_b / n_b : 0.0;

            gene_total_a += te.mean_a;
            gene_total_b += te.mean_b;

            tx_data.push_back(te);
        }

        if (gene_total_a <= 0.0 && gene_total_b <= 0.0) continue;

        // Chi-squared on 2×k contingency table (2 groups × k transcripts).
        // H0: transcript proportions are equal across groups.
        double chi2 = 0.0;
        int df = 0;
        size_t first_dtu_idx = dtu_results.size();

        for (auto& te : tx_data) {
            double prop_a = (gene_total_a > 0) ? te.mean_a / gene_total_a : 0.0;
            double prop_b = (gene_total_b > 0) ? te.mean_b / gene_total_b : 0.0;

            double total = te.mean_a + te.mean_b;
            double row_total = gene_total_a + gene_total_b;

            if (total <= 0.0 || row_total <= 0.0) continue;

            double expected_a = total * gene_total_a / row_total;
            double expected_b = total * gene_total_b / row_total;

            if (expected_a > 0.0) {
                chi2 += (te.mean_a - expected_a) * (te.mean_a - expected_a) / expected_a;
            }
            if (expected_b > 0.0) {
                chi2 += (te.mean_b - expected_b) * (te.mean_b - expected_b) / expected_b;
            }

            dtu_result dr;
            dr.gene_id = gene_id;
            dr.gene_name = te.gene_name;
            dr.transcript_id = te.transcript_id;
            dr.prop_group_a = prop_a;
            dr.prop_group_b = prop_b;
            dr.delta_proportion = prop_a - prop_b;
            dr.p_value = 1.0;
            dr.fdr = 1.0;
            dr.significant = false;

            dtu_results.push_back(dr);
            df++;
        }

        // Wilson-Hilferty approximation for the chi-squared CDF. Returns
        // P(X > chi2) with df = k - 1. Applied to every transcript in the
        // gene (gene-level p-value — the test is joint across transcripts).
        if (df > 1 && chi2 > 0.0) {
            int k = df - 1;
            double z = std::pow(chi2 / k, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * k));
            z /= std::sqrt(2.0 / (9.0 * k));
            double p = 0.5 * std::erfc(z / std::sqrt(2.0));
            p = std::max(p, 1e-300);

            for (size_t i = first_dtu_idx; i < dtu_results.size(); ++i) {
                dtu_results[i].p_value = p;
            }
        }
    }

    return dtu_results;
}

void query::apply_fdr_correction(std::vector<dtu_result>& results,
                                  double fdr_threshold) {
    if (results.empty()) return;

    std::vector<size_t> indices(results.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
              [&results](size_t a, size_t b) {
                  return results[a].p_value < results[b].p_value;
              });

    size_t n = results.size();
    double prev_fdr = 1.0;

    // Benjamini–Hochberg: fdr_i = p_i * n / rank_i, enforcing monotonicity
    // from largest rank down so fdr stays non-decreasing with rank.
    for (size_t i = n; i > 0; --i) {
        size_t idx = indices[i - 1];
        double fdr = results[idx].p_value * static_cast<double>(n)
                                          / static_cast<double>(i);
        fdr = std::min(fdr, prev_fdr);
        fdr = std::min(fdr, 1.0);
        results[idx].fdr = fdr;
        results[idx].significant = (fdr <= fdr_threshold);
        prev_fdr = fdr;
    }
}

// ============================================================================
// DTU output
// ============================================================================

void query::write_dtu_results(const std::string& path,
                               const query_contrast& contrast,
                               const std::vector<dtu_result>& results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open output file: " + path);
        return;
    }

    out << "gene_id\tgene_name\ttranscript_id\t"
        << "prop_" << contrast.group_a << "\t"
        << "prop_" << contrast.group_b << "\t"
        << "delta_proportion\tp_value\tfdr\tsignificant\n";

    for (const auto& r : results) {
        out << r.gene_id << "\t"
            << r.gene_name << "\t"
            << r.transcript_id << "\t"
            << r.prop_group_a << "\t"
            << r.prop_group_b << "\t"
            << r.delta_proportion << "\t"
            << r.p_value << "\t"
            << r.fdr << "\t"
            << (r.significant ? "yes" : "no") << "\n";
    }

    logging::info("DTU results written to: " + path);
}

void query::write_summary(const std::string& path,
                           const std::vector<query_result>& classification,
                           const std::vector<std::pair<query_contrast,
                               std::vector<dtu_result>>>& dtu_results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open summary file: " + path);
        return;
    }

    out << "# Atroplex Query Summary\n\n";

    // Classification summary
    out << "## Structural Classification\n";
    out << "Total transcripts: " << classification.size() << "\n\n";

    std::map<structural_category, size_t> counts;
    for (const auto& r : classification) {
        counts[r.category]++;
    }

    auto pct = [&](size_t count) {
        if (classification.empty()) return 0.0;
        return 100.0 * static_cast<double>(count) /
               static_cast<double>(classification.size());
    };

    for (const auto& [cat, count] : counts) {
        out << to_string(cat) << ": " << count
            << " (" << pct(count) << "%)\n";
    }

    // DTU summary
    if (!dtu_results.empty()) {
        out << "\n## Differential Transcript Usage\n";
        for (const auto& [contrast, results] : dtu_results) {
            size_t sig = std::count_if(results.begin(), results.end(),
                [](const dtu_result& r) { return r.significant; });
            out << contrast.group_a << " vs " << contrast.group_b << ": "
                << sig << "/" << results.size() << " significant\n";
        }
    }

    logging::info("Summary written to: " + path);
}

} // namespace subcall