#include "base.hpp"

// standard
#include <fstream>
#include <iostream>

// class
#include "utility.hpp"
#include "genogrove_builder.hpp"

// genogrove
#include <genogrove/io/gff_reader.hpp>
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/filetype_detector.hpp>

namespace gdt = genogrove::data_type;
namespace ggs = genogrove::structure;
namespace gio = genogrove::io;

base::base(const cxxopts::ParseResult& params)
    : params{params}, ftype{gio::filetype::UNKNOWN}, compression{false} {

    // Check if a genogrove file is provided - if so load genogrove structure
    if(params.count("genogrove")) {
        std::string gg_path = params["genogrove"].as<std::string>();
        if(std::filesystem::exists(gg_path)) {
            logging::info("Loading existing genogrove from: " + gg_path);
            load_genogrove(gg_path);
        } else {
            logging::warning("Genogrove file not found: " + gg_path);
        }
    }

    if(params.count("build-from")) {
        auto build_files = params["build-from"].as<std::vector<std::string>>();
        create_genogrove(build_files);
    }
}

base::~base() {
    if(structures) {
        delete structures;
    }
}

void base::validate(const cxxopts::ParseResult& args){
    // check if the input file exists
    if(!std::filesystem::exists(args["input"].as<std::string>())) {
        logging::error("Input file does not exist: " + args["input"].as<std::string>());
        throw std::runtime_error("Input file not found");
    }
}

void base::process() {
    logging::info("Starting atroplex pipeline...");
    start();
    // TODO: Add further processing steps
}

void base::start() {
    logging::info("Initializing genogrove structure...");






    // If no structures loaded and build-from files are provided, create new structures
    if(!structures && params.count("build-from")) {
        auto build_files = params["build-from"].as<std::vector<std::string>>();

        logging::info("Creating genomic structures from " + std::to_string(build_files.size()) + " file(s)");
        create_genogrove(build_files);
    }

    if(!structures) {
        logging::warning("No genomic structures loaded or created");
    } else {
        logging::info("Genomic structures ready (grove + transcript graph)");
    }
}

void base::create_genogrove(const std::vector<std::string>& build_files) {
    int order = params["order"].as<int>();

    // Use the genogrove_builder to handle multiple files and file types
    structures = genogrove_builder::build_from_files(build_files, order);

    if (!structures) {
        logging::error("Failed to create genomic structures from provided files");
    }
}

void base::load_genogrove(const std::string& gg_path) {
    // TODO: Implement genogrove deserialization
    logging::error("Genogrove loading not yet implemented");
}

void base::detect_input_filetype(){
    logging::info("Detecting input file type...");

    std::string input_path = params["input"].as<std::string>();
    gio::filetype_detector detector;
    auto [detected_ftype, is_gzipped] = detector.detect_filetype(input_path);

    ftype = detected_ftype;
    gzipped = is_gzipped;

    std::string ftype_str;
    switch(detected_ftype) {
        case gio::filetype::FASTQ:
            ftype_str = "FASTQ";
            break;
        case gio::filetype::BAM:
            ftype_str = "BAM";
            break;
        case gio::filetype::GFF:
            ftype_str = "GFF";
            break;
        case gio::filetype::GTF:
            ftype_str = "GTF";
            break;
        case gio::filetype::GG:
            ftype_str = "Genogrove";
            break;
        default:
            ftype_str = "UNKNOWN";
            logging::error("Unknown file type for input file");
            throw std::runtime_error("Unsupported file type");
    }

    logging::info("Detected file type: " + ftype_str +
                 (is_gzipped ? " (gzipped)" : ""));
}

void base::align_reads() {
    // TODO: Implement alignment for FASTQ input
    logging::warning("Alignment not yet implemented");
}

