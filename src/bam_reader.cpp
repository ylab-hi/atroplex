#include "bam_reader.hpp"

bam_reader::bam_reader(const std::filesystem::path& path)
    : line_num(0), eof_reached(false) {
    file = sam_open(path.c_str(), "r");
    if (!file) {
        throw std::runtime_error("Failed to open BAM/SAM file: " + path.string());
    }

    header = sam_hdr_read(file);
    if (!header) {
        sam_close(file);
        throw std::runtime_error("Failed to read BAM/SAM header: " + path.string());
    }

    record = bam_init1();
    if (!record) {
        bam_hdr_destroy(header);
        sam_close(file);
        throw std::runtime_error("Failed to initialize BAM record");
    }
}

bool bam_reader::read_next(alignment_entry& entry) {
    int result = sam_read1(file, header, record);

    if (result < 0) {
        eof_reached = true;
        if (result == -1) {
            // End of file
            return false;
        } else {
            error_message = "Error reading BAM/SAM file";
            return false;
        }
    }

    line_num++;

    // Extract fields from BAM record
    entry.qname = bam_get_qname(record);
    entry.rname = header->target_name[record->core.tid];
    entry.pos = record->core.pos + 1;  // Convert to 1-based
    entry.mapq = record->core.qual;
    entry.flag = record->core.flag;

    // Get CIGAR string
    uint32_t* cigar = bam_get_cigar(record);
    entry.cigar.clear();
    for (uint32_t i = 0; i < record->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        entry.cigar += std::to_string(len);
        entry.cigar += bam_cigar_opchr(op);
    }

    // Get sequence
    uint8_t* seq = bam_get_seq(record);
    entry.seq.clear();
    for (int i = 0; i < record->core.l_qseq; ++i) {
        entry.seq += seq_nt16_str[bam_seqi(seq, i)];
    }

    // Get quality
    uint8_t* qual = bam_get_qual(record);
    entry.qual.clear();
    for (int i = 0; i < record->core.l_qseq; ++i) {
        entry.qual += static_cast<char>(qual[i] + 33);
    }

    // Template length
    entry.tlen = record->core.isize;

    return true;
}

bool bam_reader::has_next() {
    return !eof_reached;
}

std::string bam_reader::get_error_message() {
    return error_message;
}

size_t bam_reader::get_current_line() {
    return line_num;
}

bam_reader::~bam_reader() {
    if (record) {
        bam_destroy1(record);
    }
    if (header) {
        bam_hdr_destroy(header);
    }
    if (file) {
        sam_close(file);
    }
}