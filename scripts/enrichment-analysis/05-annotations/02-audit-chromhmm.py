#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    contig_lengths = {}
    contig_order = {}
    with open(args.fai) as handle:
        for index, line in enumerate(handle):
            name, length, *_ = line.rstrip("\n").split("\t")
            contig_lengths[name] = int(length)
            contig_order[name] = index

    counts = Counter()
    state_intervals = Counter()
    state_bases = Counter()
    previous = None
    union_chrom = None
    union_end = 0
    with gzip.open(args.input, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                counts["malformed"] += 1
                continue
            chrom, start_text, end_text, state = fields[:4]
            start, end = int(start_text), int(end_text)
            counts["intervals"] += 1
            counts["bases"] += max(0, end - start)
            state_intervals[state] += 1
            state_bases[state] += max(0, end - start)
            if chrom not in contig_lengths:
                counts["unknown_contig"] += 1
                continue
            if start < 0 or end <= start or end > contig_lengths[chrom]:
                counts["invalid_interval"] += 1
            current = (contig_order[chrom], start)
            if previous is not None and current < previous[:2]:
                counts["sort_errors"] += 1
            if previous is not None and chrom == previous[2] and start < previous[4]:
                counts["overlap_pairs"] += 1
                counts["overlap_bases"] += min(end, previous[4]) - start
                if state != previous[5]:
                    counts["conflicting_state_overlaps"] += 1
            previous = (*current, chrom, start, end, state)
            if chrom != union_chrom:
                union_chrom = chrom
                union_end = end
                counts["covered_bases"] += end - start
            elif start >= union_end:
                counts["covered_bases"] += end - start
                union_end = end
            elif end > union_end:
                counts["covered_bases"] += end - union_end
                union_end = end

    with pysam.TabixFile(str(args.input)) as tabix:
        counts["tabix_contigs"] = len(tabix.contigs)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for metric in sorted(counts):
            writer.writerow([metric, counts[metric]])
        for state in sorted(state_intervals):
            writer.writerow([f"state_intervals:{state}", state_intervals[state]])
            writer.writerow([f"state_bases:{state}", state_bases[state]])

    critical = ["malformed", "unknown_contig", "invalid_interval", "sort_errors", "overlap_pairs"]
    failures = [f"{key}={counts[key]}" for key in critical if counts[key]]
    if failures:
        raise SystemExit("ChromHMM audit failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
