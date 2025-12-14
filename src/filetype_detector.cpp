#include "filetype_detector.hpp"

#include <cstring>
#include <fstream>


#include <zlib.h>

std::tuple<filetype, bool> filetype_detector::detect_filetype(
    const std::filesystem::path& filepath) {

    // Check for .gg extension (genogrove format)
    std::string extension = filepath.extension().string();
    if (!extension.empty() && extension[0] == '.') {
        extension = extension.substr(1);
    }
    if (extension == "gg") {
        return std::make_tuple(filetype::GG, false);
    }

    std::ifstream file(filepath, std::ios::binary);
    if(!file) {
        throw std::runtime_error("Failed to open file: " + filepath.string());
    }

    // Read first few bytes to check magic numbers
    char buffer[256];
    file.read(buffer, sizeof(buffer));
    std::streamsize bytes_read = file.gcount();
    file.close();

    if (bytes_read < 2) {
        return std::make_tuple(filetype::UNKNOWN, false);
    }

    // Check if the file is gzipped (magic bytes: 0x1f 0x8b)
    bool is_gzipped = (static_cast<unsigned char>(buffer[0]) == 0x1f &&
                       static_cast<unsigned char>(buffer[1]) == 0x8b);

    if (is_gzipped) {
        // For gzipped files, we need to decompress and check content
        return detect_gzipped_filetype(filepath);
    } else {
        // For plain files, check content directly
        return detect_plain_filetype(buffer, bytes_read);
    }
}

std::tuple<filetype, bool> filetype_detector::detect_plain_filetype(const char* buffer, std::streamsize size) {
    if (size < 4) {
        return std::make_tuple(filetype::UNKNOWN, false);
    }

    // Check for GFF/GTF version header
    if (size >= 14 && strncmp(buffer, "##gff-version", 13) == 0) {
        // Check if it's GFF3 or GTF
        if (size >= 16 && buffer[14] == '3') {
            return std::make_tuple(filetype::GFF, false);
        } else if (size >= 16 && buffer[14] == '2') {
            return std::make_tuple(filetype::GTF, false);
        }
        return std::make_tuple(filetype::GFF, false);  // Default to GFF
    }

    // Check for GTF format (often starts with chromosome line without version header)
    // GTF format: chr \t source \t feature \t start \t end \t score \t strand \t frame \t attributes
    if (buffer[0] != '#' && buffer[0] != '@') {
        // Look for tab-separated fields typical of GFF/GTF
        int tab_count = 0;
        for (std::streamsize i = 0; i < size && buffer[i] != '\n'; ++i) {
            if (buffer[i] == '\t') {
                tab_count++;
            }
        }
        // GFF/GTF has 8 tabs (9 fields)
        if (tab_count >= 7) {
            // Look for typical GTF/GFF keywords in the line
            std::string line(buffer, std::min(size, static_cast<std::streamsize>(200)));
            if (line.find("gene_id") != std::string::npos ||
                line.find("transcript_id") != std::string::npos) {
                return std::make_tuple(filetype::GTF, false);
            } else if (line.find("ID=") != std::string::npos ||
                       line.find("Parent=") != std::string::npos) {
                return std::make_tuple(filetype::GFF, false);
            }
        }
    }

    // Check for SAM header (@HD, @SQ, @RG, @PG, @CO)
    if (buffer[0] == '@' && (
        (size >= 3 && buffer[1] == 'H' && buffer[2] == 'D') ||
        (size >= 3 && buffer[1] == 'S' && buffer[2] == 'Q') ||
        (size >= 3 && buffer[1] == 'R' && buffer[2] == 'G') ||
        (size >= 3 && buffer[1] == 'P' && buffer[2] == 'G') ||
        (size >= 3 && buffer[1] == 'C' && buffer[2] == 'O'))) {
        return std::make_tuple(filetype::SAM, false);
    }

    // Check for FASTQ format (starts with @ but different structure than SAM)
    // FASTQ has 4-line structure: @name, sequence, +, quality
    if (buffer[0] == '@') {
        // Look for newline, then sequence line, then '+' marker
        bool could_be_fastq = false;
        for (std::streamsize i = 1; i < size - 1; ++i) {
            if (buffer[i] == '\n') {
                // After first line, should have DNA/RNA sequence
                if (i + 1 < size && (buffer[i+1] == 'A' || buffer[i+1] == 'C' ||
                                      buffer[i+1] == 'G' || buffer[i+1] == 'T' ||
                                      buffer[i+1] == 'N')) {
                    could_be_fastq = true;
                }
                break;
            }
        }
        if (could_be_fastq) {
            return std::make_tuple(filetype::FASTQ, false);
        }
    }

    return std::make_tuple(filetype::UNKNOWN, false);
}

std::tuple<filetype, bool> filetype_detector::detect_gzipped_filetype(const std::filesystem::path& filepath) {
    // Open gzipped file
    gzFile gzfile = gzopen(filepath.string().c_str(), "rb");
    if (!gzfile) {
        throw std::runtime_error("Failed to open gzipped file: " + filepath.string());
    }

    // Read decompressed header
    char buffer[256];
    int bytes_read = gzread(gzfile, buffer, sizeof(buffer));
    gzclose(gzfile);

    if (bytes_read < 4) {
        return std::make_tuple(filetype::UNKNOWN, true);
    }

    // Check for BAM magic bytes: "BAM\1"
    if (bytes_read >= 4 &&
        buffer[0] == 'B' && buffer[1] == 'A' && buffer[2] == 'M' && buffer[3] == '\1') {
        return std::make_tuple(filetype::BAM, true);
    }

    // Check for FASTQ or SAM in gzipped form
    auto [ftype, _] = detect_plain_filetype(buffer, bytes_read);
    return std::make_tuple(ftype, true);
}