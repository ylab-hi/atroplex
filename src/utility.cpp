#include "utility.hpp"

// standard
#include <iostream>
#include <iomanip>
#include <sstream>
#include <locale>
#include <mutex>

std::string normalize_chromosome(const std::string& seqid) {
    // Already has chr prefix - check for MT variant
    if (seqid.size() >= 3 && seqid.substr(0, 3) == "chr") {
        return seqid;
    }

    // Ensembl mitochondria
    if (seqid == "MT") {
        return "chrM";
    }

    // Numeric chromosomes (1-22) or sex chromosomes (X, Y) or M
    if (!seqid.empty()) {
        // Check if it's a valid chromosome to prefix
        // Numeric: 1, 2, ..., 22
        // Letters: X, Y, M
        bool is_numeric = true;
        for (char c : seqid) {
            if (!std::isdigit(c)) {
                is_numeric = false;
                break;
            }
        }

        if (is_numeric || seqid == "X" || seqid == "Y" || seqid == "M") {
            return "chr" + seqid;
        }
    }

    // Unknown format - return as-is
    return seqid;
}

namespace logging {
    // ANSI color codes
    const std::string RESET = "\033[0m";
    const std::string YELLOW = "\033[33m";
    const std::string RED = "\033[31m";
    const std::string CYAN = "\033[36m";

    // Thread safety
    static std::mutex log_mutex;

    // Progress tracking
    static std::chrono::steady_clock::time_point progress_start_time;
    static bool progress_enabled = false;  // Disabled by default (for SLURM, CI, etc.)

    void set_progress_enabled(bool enabled) {
        std::lock_guard<std::mutex> lock(log_mutex);
        progress_enabled = enabled;
    }

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
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << "[atroplex] " << get_timestamp() << " - " << message << std::endl;
    }

    void warning(const std::string& message) {
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cout << YELLOW << "[atroplex] " << get_timestamp() << " - WARNING: " << message << RESET << std::endl;
    }

    void error(const std::string& message) {
        std::lock_guard<std::mutex> lock(log_mutex);
        std::cerr << RED << "[atroplex] " << get_timestamp() << " - ERROR: " << message << RESET << std::endl;
    }

    void progress_start() {
        std::lock_guard<std::mutex> lock(log_mutex);
        progress_start_time = std::chrono::steady_clock::now();
    }

    void progress(size_t lines, const std::string& prefix) {
        std::lock_guard<std::mutex> lock(log_mutex);
        if (!progress_enabled) return;

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
        std::lock_guard<std::mutex> lock(log_mutex);
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - progress_start_time).count();
        double seconds = static_cast<double>(elapsed) / 1000.0;

        // Clear the progress line (if enabled) and print final message
        if (progress_enabled) {
            std::cout << "\r";  // Clear progress line
        }
        std::cout << "[atroplex] " << get_timestamp() << " - "
                  << prefix << " " << format_number(segments) << " segments"
                  << " in " << std::fixed << std::setprecision(1) << seconds << "s"
                  << (progress_enabled ? "                    " : "")  // Extra spaces only if clearing
                  << std::endl;
    }
}