#include "gff_reader.hpp"

// standard
#include <sstream>
#include <stdexcept>


gff_reader::gff_reader(const std::filesystem::path& filepath)
    : line_num(0), eof_reached(false), is_gtf(false) {

    file.open(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open GFF/GTF file: " + filepath.string());
    }

    // Detect format from first non-comment line or header
    std::string line;
    std::streampos start_pos = file.tellg();
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line.substr(0, 14) == "##gff-version ") {
            // Explicit GFF version
            if (line.find("2") != std::string::npos) {
                is_gtf = true;  // GFF2 is essentially GTF
            }
            break;
        } else if (line[0] != '#') {
            // Check attributes format to distinguish GTF from GFF3
            size_t attr_pos = 0;
            int tab_count = 0;
            for (size_t i = 0; i < line.length(); ++i) {
                if (line[i] == '\t') {
                    tab_count++;
                    if (tab_count == 8) {
                        attr_pos = i + 1;
                        break;
                    }
                }
            }
            if (attr_pos > 0 && attr_pos < line.length()) {
                std::string attrs = line.substr(attr_pos);
                // GTF uses space-separated "key value;" format
                // GFF3 uses semicolon-separated "key=value" format
                is_gtf = (attrs.find("gene_id") != std::string::npos ||
                          attrs.find("transcript_id") != std::string::npos) &&
                         attrs.find('=') == std::string::npos;
            }
            break;
        }
    }

    // Reset to beginning
    file.clear();
    file.seekg(start_pos);
}

gff_reader::~gff_reader() {
    if (file.is_open()) {
        file.close();
    }
}

bool gff_reader::read_next(gff_entry& entry) {
    std::string line;

    while (std::getline(file, line)) {
        line_num++;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Try to parse the line
        if (parse_line(line, entry)) {
            return true;
        }
    }

    eof_reached = true;
    return false;
}

bool gff_reader::parse_line(const std::string& line, gff_entry& entry) {
    std::istringstream iss(line);
    std::string field;
    int field_num = 0;

    // Clear previous entry
    entry = gff_entry();

    // Parse tab-separated fields
    while (std::getline(iss, field, '\t')) {
        switch (field_num) {
            case 0: // seqid
                entry.seqid = field;
                break;
            case 1: // source
                entry.source = field;
                break;
            case 2: // type
                entry.type = field;
                break;
            case 3: // start
                try {
                    entry.start = std::stoull(field);
                } catch (...) {
                    return false;
                }
                break;
            case 4: // end
                try {
                    entry.end = std::stoull(field);
                } catch (...) {
                    return false;
                }
                break;
            case 5: // score
                if (field != ".") {
                    try {
                        entry.score = std::stod(field);
                    } catch (...) {
                        entry.score = -1.0;
                    }
                }
                break;
            case 6: // strand
                if (field.length() > 0) {
                    entry.strand = field[0];
                }
                break;
            case 7: // phase
                if (field != ".") {
                    try {
                        entry.phase = std::stoi(field);
                    } catch (...) {
                        entry.phase = -1;
                    }
                }
                break;
            case 8: // attributes
                parse_attributes(field, entry);
                break;
            default:
                // Extra fields ignored
                break;
        }
        field_num++;
    }

    // Must have at least 8 fields
    return field_num >= 8;
}

void gff_reader::parse_attributes(const std::string& attr_str, gff_entry& entry) {
    if (is_gtf) {
        parse_gtf_attributes(attr_str, entry);
    } else {
        parse_gff3_attributes(attr_str, entry);
    }
}

void gff_reader::parse_gff3_attributes(const std::string& attr_str, gff_entry& entry) {
    // GFF3 format: key1=value1;key2=value2;key3=value3
    std::istringstream iss(attr_str);
    std::string pair;

    while (std::getline(iss, pair, ';')) {
        if (pair.empty()) continue;

        size_t eq_pos = pair.find('=');
        if (eq_pos != std::string::npos) {
            std::string key = pair.substr(0, eq_pos);
            std::string value = pair.substr(eq_pos + 1);

            // Trim whitespace
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            entry.attributes[key] = value;
        }
    }
}

void gff_reader::parse_gtf_attributes(const std::string& attr_str, gff_entry& entry) {
    // GTF format: key1 "value1"; key2 "value2"; key3 "value3";
    std::istringstream iss(attr_str);
    std::string token;
    std::string key;
    bool reading_key = true;

    while (iss >> token) {
        if (reading_key) {
            key = token;
            reading_key = false;
        } else {
            // Remove quotes and semicolon
            std::string value = token;
            if (value.front() == '"') value.erase(0, 1);
            if (value.back() == ';') value.pop_back();
            if (value.back() == '"') value.pop_back();

            entry.attributes[key] = value;
            reading_key = true;
        }
    }
}