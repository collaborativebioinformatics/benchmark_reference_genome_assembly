#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--callsets", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    counts = Counter()
    fieldnames = None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.callsets, newline="") as callsets, gzip.open(args.output, "wt") as output:
        for callset in csv.DictReader(callsets, delimiter="\t"):
            arm, build = callset["arm"], callset["build"]
            path = args.root / "results/annotations/variants" / arm / build / "annotations.tsv.gz"
            with gzip.open(path, "rt") as source:
                reader = csv.DictReader(source, delimiter="\t")
                if fieldnames is None:
                    fieldnames = reader.fieldnames
                    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t")
                    writer.writeheader()
                elif reader.fieldnames != fieldnames:
                    raise SystemExit(f"Header mismatch: {path}")
                for row in reader:
                    row["mappability"] = (
                        "high" if float(row["mappability_unique_k100_fraction"]) >= 0.5 else "low"
                    )
                    writer.writerow(row)
                    counts[(arm, build, row["benchmark_status"])] += 1

    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["arm", "reference", "benchmark_status", "rows"])
        for key in sorted(counts):
            writer.writerow([*key, counts[key]])


if __name__ == "__main__":
    main()
