#!/usr/bin/env python3
"""
Inject expression counts from TALON quantification TSV into GTF transcript attributes.
Matches on talon_transcript ID and adds counts "N" attribute to transcript lines.
Also prepends ##property: value headers for atroplex sample_info parsing.

Usage:
    python inject_expression.py <gtf> <tsv> <output> [--header key=value ...]

Example:
    python inject_expression.py \
        rep01_ENCFF994XIB.gtf \
        rep01_ENCFF275JBU.tsv \
        rep01_ENCFF994XIB.counts.gtf \
        --header id=ENCSR583KAF_rep01 \
        --header assay="long-read RNA-seq" \
        --header biosample=GM12878 \
        --header pipeline=TALON
"""

import argparse
import re
import sys


def load_counts(tsv_path):
    """Load transcript_ID -> count mapping from TALON quantification TSV."""
    counts = {}
    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        count_col = len(header) - 1
        tid_col = header.index('transcript_ID')

        for line in f:
            fields = line.strip().split('\t')
            talon_tid = fields[tid_col]
            count = int(fields[count_col])
            counts[talon_tid] = count

    return counts


def extract_talon_transcript(attrs_str):
    """Extract talon_transcript value from GTF attribute string."""
    match = re.search(r'talon_transcript "([^"]+)"', attrs_str)
    return match.group(1) if match else None


def merge_gtf(gtf_path, counts, output_path, headers):
    """Add counts attribute to transcript lines in GTF."""
    matched = 0
    unmatched = 0

    with open(gtf_path) as fin, open(output_path, 'w') as fout:
        # Write atroplex headers
        for key, value in headers.items():
            fout.write(f"##{key}: {value}\n")

        for line in fin:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                fout.write(line)
                continue

            feature_type = fields[2]

            if feature_type == 'transcript':
                talon_tid = extract_talon_transcript(fields[8])
                if talon_tid and talon_tid in counts:
                    count = counts[talon_tid]
                    attrs = fields[8].rstrip()
                    if not attrs.endswith(';'):
                        attrs += ';'
                    attrs += f' counts "{count}";'
                    fields[8] = attrs
                    matched += 1
                else:
                    unmatched += 1

            fout.write('\t'.join(fields) + '\n')

    return matched, unmatched


def main():
    parser = argparse.ArgumentParser(
        description='Inject TALON expression counts into GTF transcript attributes')
    parser.add_argument('gtf', help='Input GTF file (TALON output)')
    parser.add_argument('tsv', help='TALON quantification TSV')
    parser.add_argument('output', help='Output GTF with counts injected')
    parser.add_argument('--header', action='append', default=[],
                        help='Add ##key: value header (can specify multiple times)')
    args = parser.parse_args()

    # Parse --header key=value pairs
    headers = {}
    for h in args.header:
        if '=' not in h:
            print(f"Error: header must be key=value, got: {h}", file=sys.stderr)
            sys.exit(1)
        key, value = h.split('=', 1)
        headers[key] = value

    print(f"Loading counts from: {args.tsv}")
    counts = load_counts(args.tsv)
    print(f"Loaded {len(counts)} transcript counts")

    print(f"Merging into: {args.output}")
    matched, unmatched = merge_gtf(args.gtf, counts, args.output, headers)
    print(f"Matched: {matched}, Unmatched: {unmatched}")


if __name__ == '__main__':
    main()