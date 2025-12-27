#include "base.hpp"

// standard
#include <fstream>
#include <iostream>

// class
#include "utility.hpp"

// genogrove
#include <genogrove/io/gff_reader.hpp>
#include <genogrove/structure/grove/grove.hpp>
#include <genogrove/data_type/interval.hpp>
#include <genogrove/io/filetype_detector.hpp>

namespace gdt = genogrove::data_type;
namespace ggs = genogrove::structure;
namespace gio = genogrove::io;

base::base(const cxxopts::ParseResult& params)
    : params{params}, ftype{filetype::UNKNOWN}, gzipped{false} {
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

    // Check if a genogrove file is provided
    if(params.count("genogrove")) {
        std::string gg_path = params["genogrove"].as<std::string>();
        if(std::filesystem::exists(gg_path)) {
            logging::info("Loading existing genogrove from: " + gg_path);
            load_genogrove(gg_path);
        } else {
            logging::warning("Genogrove file not found: " + gg_path);
        }
    }

    // If no grove loaded and annotation is provided, create new grove
    if(!grove && params.count("annotation")) {
        std::string annotation_path = params["annotation"].as<std::string>();
        if(!std::filesystem::exists(annotation_path)) {
            logging::error("Annotation file not found: " + annotation_path);
            throw std::runtime_error("Annotation file not found");
        }

        logging::info("Creating genogrove from annotation: " + annotation_path);
        create_genogrove(annotation_path);
    }

    if(!grove) {
        logging::warning("No genogrove structure loaded or created");
    } else {
        logging::info("Genogrove structure ready");
    }
}

void base::create_genogrove(const std::string& annotation_path) {
    int order = params["order"].as<int>();
    logging::info("Creating grove with order: " + std::to_string(order));

    // Create grove with interval as key and gff_entry as data
    grove = new ggs::grove<gdt::interval, gff_entry>(order);

    // Read and insert GFF entries
    gff_reader reader(annotation_path);
    gff_entry entry;
    size_t count = 0;

    while(reader.has_next()) {
        if(!reader.read_next(entry)) {
            if(reader.get_error_message().empty()) {
                break; // EOF
            }
            logging::error("Error reading annotation at line " +
                         std::to_string(reader.get_current_line()) +
                         ": " + reader.get_error_message());
            continue;
        }

        // Insert entry: seqid as index, interval as key, gff_entry as data
        grove->insert_data(entry.seqid, entry.interval, entry);
        count++;
    }

    logging::info("Inserted " + std::to_string(count) + " entries into genogrove");
}

void base::load_genogrove(const std::string& gg_path) {
    // TODO: Implement genogrove deserialization
    logging::error("Genogrove loading not yet implemented");
}

void base::detect_input_filetype(){
    logging::info("Detecting input file type...");

    std::string input_path = params["input"].as<std::string>();
    filetype_detector detector;
    auto [detected_ftype, is_gzipped] = detector.detect_filetype(input_path);

    ftype = detected_ftype;
    gzipped = is_gzipped;

    std::string ftype_str;
    switch(detected_ftype) {
        case filetype::FASTQ:
            ftype_str = "FASTQ";
            break;
        case filetype::BAM:
            ftype_str = "BAM";
            break;
        case filetype::GFF:
            ftype_str = "GFF";
            break;
        case filetype::GTF:
            ftype_str = "GTF";
            break;
        case filetype::GG:
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

