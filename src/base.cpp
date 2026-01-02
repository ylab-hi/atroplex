#include "base.hpp"

// standard
#include <fstream>
#include <iostream>

// class
#include "utility.hpp"
#include "builder.hpp"

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
    if(grove) {
        delete grove;
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

    // If no grove loaded and build-from files are provided, create new grove
    if(!grove && params.count("build-from")) {
        auto build_files = params["build-from"].as<std::vector<std::string>>();

        logging::info("Creating grove from " + std::to_string(build_files.size()) + " file(s)");
        create_genogrove(build_files);
    }

    if(!grove) {
        logging::warning("No grove loaded or created");
    } else {
        logging::info("Grove ready with spatial index and graph structure");
    }
}

void base::create_genogrove(const std::vector<std::string>& build_files) {
    int order = params["order"].as<int>();

    // Use the genogrove_builder to handle multiple files and file types
    grove = builder::build_from_files(build_files, order);

    if (!grove) {
        logging::error("Failed to create grove from provided files");
    }
}

void base::load_genogrove(const std::string& gg_path) {
    // TODO: Implement genogrove deserialization
    logging::error("Genogrove loading not yet implemented");
}

void base::align_reads() {
    // TODO: Implement alignment for FASTQ input
    logging::warning("Alignment not yet implemented");
}