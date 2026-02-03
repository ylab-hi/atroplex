/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include <iostream>
#include <memory>
#include <string>

#include <cxxopts.hpp>

#include "config.hpp"
#include "utility.hpp"

#include "subcall/subcall.hpp"
#include "subcall/discover.hpp"
#include "subcall/build.hpp"
#include "subcall/dtu.hpp"

std::unique_ptr<subcall::subcall> create_subcall(const std::string& name) {
    if (name == "discover") {
        return std::make_unique<subcall::discover>();
    } else if (name == "build") {
        return std::make_unique<subcall::build>();
    } else if (name == "dtu") {
        return std::make_unique<subcall::dtu>();
    }
    return nullptr;
}

void showVersion(std::ostream& out) {
    out << "atroplex v" << atroplex_VERSION_MAJOR;
    out << "." << atroplex_VERSION_MINOR << ".";
    out << atroplex_VERSION_PATCH << " - ";
    out << "Pan-transcriptome indexing and analysis toolkit" << std::endl;
}

void print_general_help(cxxopts::Options& options) {
    std::cout << options.help() << "\n";
    std::cout << "Subcommands:\n";
    std::cout << "\tbuild\t\tBuild pan-transcriptome index from annotation sources\n";
    std::cout << "\tdtu\t\tDifferential transcript usage analysis\n";
    std::cout << "\tdiscover\tDiscover novel transcripts from long-read data\n";
    std::cout << "\nFor subcommand help: atroplex <subcommand> --help\n";
}

int main(int argc, char** argv) {
    cxxopts::Options options("atroplex",
            "Pan-transcriptome indexing and analysis toolkit");

    options.add_options()
        ("subcall", "The subcommand to run", cxxopts::value<std::string>())
        ("h,help", "Print help")
        ("v,version", "Print version")
        ;
    options.parse_positional({"subcall"});
    options.allow_unrecognised_options();

    // Parse command-line arguments
    auto args = options.parse(argc, argv);

    // Handle version flag
    if (args.count("version")) {
        showVersion(std::cout);
        return 0;
    }

    // Handle help without subcommand
    if (args.count("help") && !args.count("subcall")) {
        print_general_help(options);
        return 0;
    }

    // Require subcommand
    if (!args.count("subcall")) {
        logging::error("No subcommand specified.");
        print_general_help(options);
        return 1;
    }

    std::string subcall_name = args["subcall"].as<std::string>();
    std::unique_ptr<subcall::subcall> command = create_subcall(subcall_name);

    if (!command) {
        logging::error("Unknown subcommand: " + subcall_name);
        print_general_help(options);
        return 1;
    }

    // Parse subcommand-specific options
    try {
        cxxopts::Options subcall_options = command->parse_args(argc - 1, argv + 1);
        cxxopts::ParseResult subcall_args = subcall_options.parse(argc - 1, argv + 1);

        // Handle help for subcommand
        if (subcall_args.count("help")) {
            std::cout << subcall_options.help() << "\n";
            return 0;
        }

        command->run(subcall_args);

    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Run 'atroplex " << subcall_name << " --help' for usage.\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}