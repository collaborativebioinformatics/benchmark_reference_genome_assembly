#!/usr/bin/env python3

import argparse
import bisect
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


BUILDS = ("hg15", "hg16", "hg17", "hg18", "hg19", "hg38", "hs1")


def read_mask(path):
    intervals = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            chrom, start, end = line.rstrip().split("\t")[:3]
            intervals[chrom].append((int(start), int(end)))
    return {
        chrom: ([start for start, _ in values], values)
        for chrom, values in intervals.items()
    }


def contained(index, chrom, start, end):
    if chrom not in index:
        return False
    starts, intervals = index[chrom]
    position = bisect.bisect_right(starts, start) - 1
    return position >= 0 and intervals[position][0] <= start and end <= intervals[position][1]


def size_bin(length):
    if length < 100:
        return "50-99"
    if length < 500:
        return "100-499"
    if length < 1000:
        return "500-999"
    if length < 5000:
        return "1k-4.9k"
    if length < 10000:
        return "5k-9.9k"
    if length < 50000:
        return "10k-49.9k"
    return ">=50k"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    masks = {}
    common = read_mask(root / "results/masks/common-7way.hs1.bed.gz")
    for build in BUILDS:
        masks[("common_7way", build)] = common
        masks[("source_specific", build)] = read_mask(root / f"results/masks/eligible/{build}.hs1.bed.gz")

    all_counts = Counter()
    kept_counts = Counter()
    with gzip.open(root / "results/annotations/sv-annotations.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length_class = size_bin(int(row["SV_length"]))
            strata = (("SV_type", row["SV_type"]), ("size_bin", length_class))
            for mode in ("common_7way", "source_specific"):
                for stratum_type, stratum in strata:
                    key = (
                        mode, row["arm"], row["reference"], row["benchmark_status"],
                        stratum_type, stratum,
                    )
                    all_counts[key] += 1
                    if contained(
                        masks[(mode, row["reference"])], row["chrom"],
                        int(row["start"]), int(row["end"]),
                    ):
                        kept_counts[key] += 1

    output = root / "results/enrichment/eligibility-by-svtype-size.tsv"
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "analysis_mask", "arm", "reference", "benchmark_status", "stratum_type",
            "stratum", "all_variants", "eligible_variants", "retained_fraction",
        ])
        for key in sorted(all_counts):
            total, kept = all_counts[key], kept_counts[key]
            writer.writerow([*key, total, kept, f"{kept / total:.8f}"])


if __name__ == "__main__":
    main()
