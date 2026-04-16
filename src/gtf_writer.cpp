/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 */

#include "gtf_writer.hpp"

namespace gtf_writer {

static void write_attributes(std::ostream& os,
                              const std::vector<std::pair<std::string, std::string>>& attrs) {
    for (size_t i = 0; i < attrs.size(); ++i) {
        if (i > 0) os << ' ';
        os << attrs[i].first << " \"" << attrs[i].second << "\";";
    }
}

void write_header(std::ostream& os, const sample_info& info) {
    os << "##provider atroplex\n";
    os << "##format GTF\n";
    os << "##id: " << info.id << "\n";
    os << "##type: " << info.type << "\n";
    if (!info.description.empty())
        os << "##description: " << info.description << "\n";
    if (!info.species.empty())
        os << "##species: " << info.species << "\n";
    if (!info.annotation_source.empty())
        os << "##annotation_source: " << info.annotation_source << "\n";
    if (!info.annotation_version.empty())
        os << "##annotation_version: " << info.annotation_version << "\n";
    if (!info.assay.empty())
        os << "##assay: " << info.assay << "\n";
    if (!info.biosample.empty())
        os << "##biosample: " << info.biosample << "\n";
    if (!info.platform.empty())
        os << "##platform: " << info.platform << "\n";
}

void write_gene(std::ostream& os, const gtf_gene& g) {
    os << g.seqid << '\t'
       << g.source << '\t'
       << "gene" << '\t'
       << g.start << '\t'
       << g.end << '\t'
       << '.' << '\t'
       << g.strand << '\t'
       << '.' << '\t';

    std::vector<std::pair<std::string, std::string>> attrs;
    attrs.emplace_back("gene_id", g.gene_id);
    if (!g.gene_name.empty())
        attrs.emplace_back("gene_name", g.gene_name);
    if (!g.gene_biotype.empty())
        attrs.emplace_back("gene_biotype", g.gene_biotype);
    write_attributes(os, attrs);
    os << '\n';
}

void write_transcript(std::ostream& os, const gtf_transcript& t) {
    os << t.seqid << '\t'
       << t.source << '\t'
       << "transcript" << '\t'
       << t.start << '\t'
       << t.end << '\t'
       << '.' << '\t'
       << t.strand << '\t'
       << '.' << '\t';

    std::vector<std::pair<std::string, std::string>> attrs;
    attrs.emplace_back("gene_id", t.gene_id);
    if (!t.gene_name.empty())
        attrs.emplace_back("gene_name", t.gene_name);
    attrs.emplace_back("transcript_id", t.transcript_id);
    if (!t.transcript_biotype.empty())
        attrs.emplace_back("transcript_type", t.transcript_biotype);

    write_attributes(os, attrs);
    os << '\n';
}

void write_exon(std::ostream& os, const gtf_exon& e) {
    os << e.seqid << '\t'
       << e.source << '\t'
       << "exon" << '\t'
       << e.start << '\t'
       << e.end << '\t'
       << '.' << '\t'
       << e.strand << '\t'
       << '.' << '\t';

    std::vector<std::pair<std::string, std::string>> attrs;
    attrs.emplace_back("gene_id", e.gene_id);
    if (!e.gene_name.empty())
        attrs.emplace_back("gene_name", e.gene_name);
    attrs.emplace_back("transcript_id", e.transcript_id);
    attrs.emplace_back("exon_number", std::to_string(e.exon_number));
    write_attributes(os, attrs);
    os << '\n';
}

} // namespace gtf_writer