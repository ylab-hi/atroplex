/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "subcall/dtu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <filesystem>

#include "utility.hpp"

namespace subcall {

dtu_contrast dtu_contrast::parse(const std::string& spec) {
    dtu_contrast c;
    size_t pos = spec.find(':');
    if (pos == std::string::npos) {
        throw std::invalid_argument("Invalid contrast format: " + spec +
                                    " (expected condition1:condition2)");
    }
    c.condition1 = spec.substr(0, pos);
    c.condition2 = spec.substr(pos + 1);

    if (c.condition1.empty() || c.condition2.empty()) {
        throw std::invalid_argument("Empty condition in contrast: " + spec);
    }
    return c;
}

cxxopts::Options dtu::parse_args(int argc, char** argv) {
    cxxopts::Options options("atroplex dtu",
        "Differential transcript usage analysis");

    options.add_options("Input/Output")
        ("m,manifest", "Sample manifest TSV (sample_id, file, condition, [tissue])",
            cxxopts::value<std::string>())
        ("q,quantifications", "Directory with pre-computed quantification files",
            cxxopts::value<std::string>())
        ("o,output", "Output prefix for results",
            cxxopts::value<std::string>()->default_value("dtu_results"))
        ;

    options.add_options("Analysis")
        ("c,contrast", "Condition contrast (format: condition1:condition2, can specify multiple)",
            cxxopts::value<std::vector<std::string>>())
        ("method", "Statistical method (ratio)",
            cxxopts::value<std::string>()->default_value("ratio"))
        ;

    options.add_options("Filtering")
        ("min-reads-transcript", "Minimum reads per transcript",
            cxxopts::value<uint32_t>()->default_value("10"))
        ("min-reads-gene", "Minimum reads per gene",
            cxxopts::value<uint32_t>()->default_value("20"))
        ("min-samples", "Minimum samples per condition",
            cxxopts::value<uint32_t>()->default_value("2"))
        ("fdr", "FDR threshold for significance",
            cxxopts::value<double>()->default_value("0.05"))
        ;

    options.add_options("Grove (for quantification)")
        ("g,genogrove", "Pre-built genogrove index (.gg)",
            cxxopts::value<std::string>())
        ("b,build-from", "Build grove from annotation (GFF/GTF)",
            cxxopts::value<std::vector<std::string>>())
        ("k,order", "Genogrove tree order",
            cxxopts::value<int>()->default_value("3"))
        ;

    add_common_options(options);

    return options;
}

void dtu::validate(const cxxopts::ParseResult& args) {
    if (!args.count("manifest")) {
        throw std::runtime_error("Sample manifest required. Use -m/--manifest");
    }

    if (!args.count("contrast")) {
        throw std::runtime_error("At least one contrast required. Use -c/--contrast");
    }

    std::string manifest = args["manifest"].as<std::string>();
    if (!std::filesystem::exists(manifest)) {
        throw std::runtime_error("Manifest file not found: " + manifest);
    }
}

void dtu::execute(const cxxopts::ParseResult& args) {
    apply_common_options(args);
    logging::info("Starting DTU analysis...");

    std::string manifest_path = args["manifest"].as<std::string>();
    std::string output_prefix = args["output"].as<std::string>();
    double fdr_threshold = args["fdr"].as<double>();

    // Parse contrasts
    auto contrast_strs = args["contrast"].as<std::vector<std::string>>();
    for (const auto& spec : contrast_strs) {
        contrasts.push_back(dtu_contrast::parse(spec));
    }

    // Load sample manifest
    logging::info("Loading sample manifest: " + manifest_path);
    samples = load_manifest(manifest_path);
    logging::info("Loaded " + std::to_string(samples.size()) + " samples");

    // Validate contrasts against available conditions
    std::set<std::string> conditions;
    for (const auto& s : samples) {
        conditions.insert(s.condition);
    }

    for (const auto& c : contrasts) {
        if (conditions.find(c.condition1) == conditions.end()) {
            throw std::runtime_error("Condition not found in manifest: " + c.condition1);
        }
        if (conditions.find(c.condition2) == conditions.end()) {
            throw std::runtime_error("Condition not found in manifest: " + c.condition2);
        }
    }

    // Load or compute quantifications
    if (args.count("quantifications")) {
        std::string quant_dir = args["quantifications"].as<std::string>();
        logging::info("Loading pre-computed quantifications from: " + quant_dir);
        load_quantifications(quant_dir);
    } else {
        logging::info("Quantifying transcripts for each sample...");
        quantify_samples(args);
    }

    // Run DTU analysis for each contrast
    std::vector<std::pair<dtu_contrast, std::vector<dtu_result>>> all_results;

    for (const auto& contrast : contrasts) {
        logging::info("Analyzing contrast: " + contrast.condition1 + " vs " + contrast.condition2);
        auto results = run_dtu_analysis(contrast, args);
        apply_fdr_correction(results, fdr_threshold);

        size_t sig_count = std::count_if(results.begin(), results.end(),
            [](const dtu_result& r) { return r.significant; });
        logging::info("  Significant DTU events: " + std::to_string(sig_count) +
                      "/" + std::to_string(results.size()));

        std::string contrast_prefix = output_prefix + "." +
                                      contrast.condition1 + "_vs_" + contrast.condition2;
        write_results(contrast_prefix + ".dtu.tsv", contrast, results);

        all_results.emplace_back(contrast, std::move(results));
    }

    write_proportions(output_prefix + ".proportions.tsv");
    write_summary(output_prefix + ".summary.txt", all_results, args);

    logging::info("DTU analysis complete");
}

std::vector<dtu_sample> dtu::load_manifest(const std::string& path) {
    std::vector<dtu_sample> result;
    std::ifstream file(path);

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open manifest file: " + path);
    }

    std::string line;
    std::vector<std::string> header;

    if (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string col;
        while (std::getline(iss, col, '\t')) {
            header.push_back(col);
        }
    }

    auto find_col = [&header](const std::string& name) -> int {
        auto it = std::find(header.begin(), header.end(), name);
        return (it != header.end()) ? std::distance(header.begin(), it) : -1;
    };

    int idx_id = find_col("sample_id");
    int idx_file = find_col("file");
    int idx_condition = find_col("condition");
    int idx_tissue = find_col("tissue");

    if (idx_id < 0 || idx_file < 0 || idx_condition < 0) {
        throw std::runtime_error("Manifest must have columns: sample_id, file, condition");
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::vector<std::string> fields;
        std::istringstream iss(line);
        std::string field;
        while (std::getline(iss, field, '\t')) {
            fields.push_back(field);
        }

        if (fields.size() < 3) continue;

        dtu_sample s;
        s.sample_id = fields[idx_id];
        s.file_path = fields[idx_file];
        s.condition = fields[idx_condition];
        if (idx_tissue >= 0 && static_cast<size_t>(idx_tissue) < fields.size()) {
            s.tissue = fields[idx_tissue];
        }

        for (size_t i = 0; i < header.size() && i < fields.size(); ++i) {
            if (static_cast<int>(i) != idx_id &&
                static_cast<int>(i) != idx_file &&
                static_cast<int>(i) != idx_condition &&
                static_cast<int>(i) != idx_tissue) {
                s.attributes[header[i]] = fields[i];
            }
        }

        result.push_back(std::move(s));
    }

    return result;
}

void dtu::quantify_samples(const cxxopts::ParseResult& args) {
    // Load grove if needed
    if (!grove) {
        if (args.count("genogrove")) {
            std::string gg_path = args["genogrove"].as<std::string>();
            load_grove(gg_path);
        } else if (args.count("build-from")) {
            auto build_files = args["build-from"].as<std::vector<std::string>>();
            int order = args["order"].as<int>();
            uint32_t threads = args["threads"].as<uint32_t>();
            build_grove(build_files, order, threads);
        } else {
            throw std::runtime_error("Grove required for quantification. "
                                     "Use --genogrove or --build-from");
        }
    }

    // TODO: Implement per-sample quantification using transcript_matcher
    logging::warning("Quantification not yet implemented - use pre-computed files");
    throw std::runtime_error("Per-sample quantification not yet implemented. "
                             "Run 'atroplex discover' on each sample first, "
                             "then provide --quantifications directory");
}

void dtu::load_quantifications(const std::string& quant_dir) {
    for (const auto& sample : samples) {
        std::string quant_path = quant_dir + "/" + sample.sample_id + ".atroplex.tsv";

        if (!std::filesystem::exists(quant_path)) {
            throw std::runtime_error("Quantification file not found: " + quant_path);
        }

        std::ifstream file(quant_path);
        std::string line;
        std::getline(file, line);  // Skip header

        std::vector<transcript_count> sample_counts;

        while (std::getline(file, line)) {
            if (line.empty()) continue;

            std::istringstream iss(line);
            std::string cluster_id, seqid, strand;
            size_t start, end, read_count, exon_count;
            std::string junctions, match_type;
            double junction_score;
            size_t matching_junctions, total_junctions;
            std::string transcript_ids;

            if (iss >> cluster_id >> seqid >> strand >> start >> end
                    >> read_count >> exon_count) {
                std::getline(iss, junctions, '\t');
                iss >> match_type >> junction_score
                    >> matching_junctions >> total_junctions;
                std::getline(iss, transcript_ids, '\t');

                if (!transcript_ids.empty() && transcript_ids != "." &&
                    match_type != "INTERGENIC") {
                    transcript_count tc;
                    tc.sample_id = sample.sample_id;
                    tc.transcript_id = transcript_ids;
                    tc.read_count = read_count;
                    tc.gene_id = "unknown";
                    sample_counts.push_back(tc);
                }
            }
        }

        counts_by_sample[sample.sample_id] = std::move(sample_counts);
        logging::info("  Loaded " + std::to_string(counts_by_sample[sample.sample_id].size()) +
                      " transcripts for " + sample.sample_id);
    }
}

std::vector<dtu_result> dtu::run_dtu_analysis(const dtu_contrast& contrast,
                                               const cxxopts::ParseResult& args) {
    std::vector<dtu_result> results;

    uint32_t min_samples_per_condition = args["min-samples"].as<uint32_t>();
    uint32_t min_reads_per_gene = args["min-reads-gene"].as<uint32_t>();

    std::vector<std::string> cond1_samples, cond2_samples;
    for (const auto& s : samples) {
        if (s.condition == contrast.condition1) cond1_samples.push_back(s.sample_id);
        else if (s.condition == contrast.condition2) cond2_samples.push_back(s.sample_id);
    }

    if (cond1_samples.size() < min_samples_per_condition ||
        cond2_samples.size() < min_samples_per_condition) {
        logging::warning("Insufficient samples for contrast: " +
                         contrast.condition1 + " (" + std::to_string(cond1_samples.size()) +
                         ") vs " + contrast.condition2 + " (" +
                         std::to_string(cond2_samples.size()) + ")");
        return results;
    }

    // Aggregate counts by gene and transcript
    std::map<std::string, std::map<std::string, std::map<std::string, std::vector<uint32_t>>>> gene_data;

    for (const auto& [sample_id, counts] : counts_by_sample) {
        std::string condition;
        for (const auto& s : samples) {
            if (s.sample_id == sample_id) {
                condition = s.condition;
                break;
            }
        }

        if (condition != contrast.condition1 && condition != contrast.condition2) {
            continue;
        }

        for (const auto& tc : counts) {
            gene_data[tc.gene_id][tc.transcript_id][condition].push_back(tc.read_count);
        }
    }

    // Calculate DTU for each gene/transcript
    for (const auto& [gene_id, transcripts] : gene_data) {
        std::map<std::string, double> gene_totals;
        for (const auto& [tx_id, cond_counts] : transcripts) {
            for (const auto& [cond, counts] : cond_counts) {
                gene_totals[cond] += std::accumulate(counts.begin(), counts.end(), 0.0);
            }
        }

        if (gene_totals[contrast.condition1] < min_reads_per_gene ||
            gene_totals[contrast.condition2] < min_reads_per_gene) {
            continue;
        }

        for (const auto& [tx_id, cond_counts] : transcripts) {
            dtu_result r;
            r.gene_id = gene_id;
            r.transcript_id = tx_id;

            double mean1 = 0, mean2 = 0;
            if (cond_counts.count(contrast.condition1)) {
                const auto& counts = cond_counts.at(contrast.condition1);
                mean1 = std::accumulate(counts.begin(), counts.end(), 0.0) / counts.size();
            }
            if (cond_counts.count(contrast.condition2)) {
                const auto& counts = cond_counts.at(contrast.condition2);
                mean2 = std::accumulate(counts.begin(), counts.end(), 0.0) / counts.size();
            }

            r.prop_condition1 = (gene_totals[contrast.condition1] > 0) ?
                                mean1 / gene_totals[contrast.condition1] : 0;
            r.prop_condition2 = (gene_totals[contrast.condition2] > 0) ?
                                mean2 / gene_totals[contrast.condition2] : 0;

            r.delta_proportion = r.prop_condition1 - r.prop_condition2;
            r.p_value = 1.0;  // Placeholder for proper statistical test
            r.fdr = 1.0;
            r.significant = false;

            results.push_back(r);
        }
    }

    return results;
}

void dtu::apply_fdr_correction(std::vector<dtu_result>& results, double fdr_threshold) {
    if (results.empty()) return;

    std::vector<size_t> indices(results.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
              [&results](size_t a, size_t b) {
                  return results[a].p_value < results[b].p_value;
              });

    size_t n = results.size();
    double prev_fdr = 1.0;

    for (size_t i = n; i > 0; --i) {
        size_t idx = indices[i - 1];
        double fdr = results[idx].p_value * n / i;
        fdr = std::min(fdr, prev_fdr);
        fdr = std::min(fdr, 1.0);
        results[idx].fdr = fdr;
        results[idx].significant = (fdr <= fdr_threshold);
        prev_fdr = fdr;
    }
}

void dtu::write_results(const std::string& path, const dtu_contrast& contrast,
                        const std::vector<dtu_result>& results) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open output file: " + path);
        return;
    }

    out << "gene_id\tgene_name\ttranscript_id\t"
        << "prop_" << contrast.condition1 << "\t"
        << "prop_" << contrast.condition2 << "\t"
        << "delta_proportion\tp_value\tfdr\tsignificant\n";

    std::vector<const dtu_result*> sorted;
    for (const auto& r : results) {
        sorted.push_back(&r);
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const dtu_result* a, const dtu_result* b) {
                  if (a->significant != b->significant) return a->significant > b->significant;
                  return std::abs(a->delta_proportion) > std::abs(b->delta_proportion);
              });

    for (const auto* r : sorted) {
        out << r->gene_id << "\t"
            << r->gene_name << "\t"
            << r->transcript_id << "\t"
            << r->prop_condition1 << "\t"
            << r->prop_condition2 << "\t"
            << r->delta_proportion << "\t"
            << r->p_value << "\t"
            << r->fdr << "\t"
            << (r->significant ? "yes" : "no") << "\n";
    }

    logging::info("  Results written to: " + path);
}

void dtu::write_proportions(const std::string& path) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open output file: " + path);
        return;
    }

    out << "gene_id\ttranscript_id";
    for (const auto& s : samples) {
        out << "\t" << s.sample_id;
    }
    out << "\n";

    std::map<std::pair<std::string, std::string>, std::map<std::string, double>> proportions;

    for (const auto& [sample_id, counts] : counts_by_sample) {
        std::map<std::string, double> gene_totals;
        for (const auto& tc : counts) {
            gene_totals[tc.gene_id] += tc.read_count;
        }

        for (const auto& tc : counts) {
            auto key = std::make_pair(tc.gene_id, tc.transcript_id);
            double prop = (gene_totals[tc.gene_id] > 0) ?
                          tc.read_count / gene_totals[tc.gene_id] : 0;
            proportions[key][sample_id] = prop;
        }
    }

    for (const auto& [key, sample_props] : proportions) {
        out << key.first << "\t" << key.second;
        for (const auto& s : samples) {
            auto it = sample_props.find(s.sample_id);
            out << "\t" << (it != sample_props.end() ? it->second : 0.0);
        }
        out << "\n";
    }

    logging::info("Proportions written to: " + path);
}

void dtu::write_summary(const std::string& path,
                        const std::vector<std::pair<dtu_contrast, std::vector<dtu_result>>>& all_results,
                        const cxxopts::ParseResult& args) {
    std::ofstream out(path);
    if (!out.is_open()) {
        logging::error("Cannot open output file: " + path);
        return;
    }

    double fdr_threshold = args["fdr"].as<double>();

    out << "=== DTU Analysis Summary ===\n\n";

    out << "Samples: " << samples.size() << "\n";
    out << "Contrasts analyzed: " << contrasts.size() << "\n\n";

    std::map<std::string, size_t> cond_counts;
    for (const auto& s : samples) {
        cond_counts[s.condition]++;
    }
    out << "Samples per condition:\n";
    for (const auto& [cond, count] : cond_counts) {
        out << "  " << cond << ": " << count << "\n";
    }
    out << "\n";

    out << "Results per contrast:\n";
    for (const auto& [contrast, results] : all_results) {
        size_t sig = std::count_if(results.begin(), results.end(),
                                   [](const dtu_result& r) { return r.significant; });
        out << "  " << contrast.condition1 << " vs " << contrast.condition2 << ":\n";
        out << "    Tested transcripts: " << results.size() << "\n";
        out << "    Significant (FDR < " << fdr_threshold << "): " << sig << "\n";
    }
    out << "\n";

    out << "Parameters:\n";
    out << "  Method: " << args["method"].as<std::string>() << "\n";
    out << "  FDR threshold: " << fdr_threshold << "\n";
    out << "  Min reads per transcript: " << args["min-reads-transcript"].as<uint32_t>() << "\n";
    out << "  Min reads per gene: " << args["min-reads-gene"].as<uint32_t>() << "\n";
    out << "  Min samples per condition: " << args["min-samples"].as<uint32_t>() << "\n";

    logging::info("Summary written to: " + path);
}

} // namespace subcall