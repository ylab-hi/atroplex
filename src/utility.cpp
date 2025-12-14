#include "utility.hpp"

// standard
#include <iostream>

namespace logging {
    // ANSI color codes
    const std::string RESET = "\033[0m";
    const std::string YELLOW = "\033[33m";
    const std::string RED = "\033[31m";

    // Internal helper to get timestamp
    static std::string get_timestamp() {
        auto now = std::chrono::system_clock::now();
        std::time_t current_time = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&current_time), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }

    void info(const std::string& message) {
        std::cout << "[MORRIGAN] " << get_timestamp() << " - " << message << std::endl;
    }

    void warning(const std::string& message) {
        std::cout << YELLOW << "[MORRIGAN] " << get_timestamp() << " - WARNING: " << message << RESET << std::endl;
    }

    void error(const std::string& message) {
        std::cerr << RED << "[MORRIGAN] " << get_timestamp() << " - ERROR: " << message << RESET << std::endl;
    }
}