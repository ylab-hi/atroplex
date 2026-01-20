#include "utility.hpp"

// standard
#include <iostream>
#include <iomanip>
#include <sstream>
#include <locale>

namespace logging {
    // ANSI color codes
    const std::string RESET = "\033[0m";
    const std::string YELLOW = "\033[33m";
    const std::string RED = "\033[31m";
    const std::string CYAN = "\033[36m";

    // Progress tracking
    static std::chrono::steady_clock::time_point progress_start_time;

    // Internal helper to get timestamp
    static std::string get_timestamp() {
        auto now = std::chrono::system_clock::now();
        std::time_t current_time = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&current_time), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }

    // Format number with commas (e.g., 1234567 -> "1,234,567")
    static std::string format_number(size_t n) {
        std::string s = std::to_string(n);
        int insert_pos = static_cast<int>(s.length()) - 3;
        while (insert_pos > 0) {
            s.insert(insert_pos, ",");
            insert_pos -= 3;
        }
        return s;
    }

    void info(const std::string& message) {
        std::cout << "[atroplex] " << get_timestamp() << " - " << message << std::endl;
    }

    void warning(const std::string& message) {
        std::cout << YELLOW << "[atroplex] " << get_timestamp() << " - WARNING: " << message << RESET << std::endl;
    }

    void error(const std::string& message) {
        std::cerr << RED << "[atroplex] " << get_timestamp() << " - ERROR: " << message << RESET << std::endl;
    }

    void progress_start() {
        progress_start_time = std::chrono::steady_clock::now();
    }

    void progress(size_t lines, const std::string& prefix) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - progress_start_time).count();

        // Calculate rate (lines per second)
        double rate = (elapsed > 0) ? (static_cast<double>(lines) / elapsed * 1000.0) : 0.0;

        std::cout << "\r" << CYAN << "[atroplex] " << prefix << "... "
                  << format_number(lines) << " lines"
                  << " (" << format_number(static_cast<size_t>(rate)) << "/sec)"
                  << RESET << "          " << std::flush;  // Extra spaces to clear previous content
    }

    void progress_done(size_t segments, const std::string& prefix) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - progress_start_time).count();
        double seconds = static_cast<double>(elapsed) / 1000.0;

        // Clear the progress line and print final message
        std::cout << "\r" << "[atroplex] " << get_timestamp() << " - "
                  << prefix << " " << format_number(segments) << " segments"
                  << " in " << std::fixed << std::setprecision(1) << seconds << "s"
                  << "                    " << std::endl;  // Extra spaces to clear, then newline
    }
}