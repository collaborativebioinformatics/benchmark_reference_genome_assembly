#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path

import pysam


CATEGORY_ORDER = {"UTR": 0, "exon": 1, "promoter": 2, "intron": 3, "intergenic": 4}


def fetch(tabix, chrom, start, end):
    try:
        return [line.split("\t") for line in tabix.fetch(chrom, start, end)]
    except ValueError:
        return []


def union_overlap(rows, start, end):
    intervals = sorted((max(start, int(row[1])), min(end, int(row[2]))) for row in rows)
    total = 0
    current_start = current_end = None
    for left, right in intervals:
        if left >= right:
            continue
        if current_start is None:
            current_start, current_end = left, right
        elif left > current_end:
            total += current_end - current_start
            current_start, current_end = left, right
        else:
            current_end = max(current_end, right)
    if current_start is not None:
        total += current_end - current_start
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--mappability", type=Path, required=True)
    parser.add_argument("--annotation-dir", type=Path, required=True)
    parser.add_argument("--chromhmm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    multi_mappability = {}
    with open(args.mappability) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            multi_mappability[fields[0]] = float(fields[4])

    tracks = {
        "genes": pysam.TabixFile(str(args.annotation_dir / "genes.bed.gz")),
        "features": pysam.TabixFile(str(args.annotation_dir / "features.bed.gz")),
        "repeatmasker": pysam.TabixFile(str(args.annotation_dir / "repeatmasker.bed.gz")),
        "segdup": pysam.TabixFile(str(args.annotation_dir / "segdup.bed.gz")),
        "unique": pysam.TabixFile(str(args.annotation_dir / "mappability-unique-k100.bed.gz")),
        "chromhmm": pysam.TabixFile(str(args.chromhmm)),
    }

    with open(args.metadata, newline="") as source, gzip.open(args.output, "wt", newline="") as output:
        reader = csv.DictReader(source, delimiter="\t")
        extra = [
            "gene", "GENCODE_category", "ChromHMM_state", "repeat_type",
            "segmental_duplication", "mappability", "mappability_multi_k100_mean",
            "mappability_unique_k100_fraction",
        ]
        writer = csv.DictWriter(output, fieldnames=[*reader.fieldnames, *extra], delimiter="\t")
        writer.writeheader()
        for row in reader:
            chrom, start, end = row["chrom"], int(row["start"]), int(row["end"])
            gene_rows = fetch(tracks["genes"], chrom, start, end)
            feature_rows = fetch(tracks["features"], chrom, start, end)
            genes = {item[3] for item in gene_rows if len(item) > 3 and item[3] != "."}
            genes.update(item[4] for item in feature_rows if len(item) > 4 and item[4] != ".")
            categories = {item[3] for item in feature_rows if len(item) > 3}
            if not categories:
                categories = {"intron"} if gene_rows else {"intergenic"}
            category_text = ";".join(sorted(categories, key=lambda value: CATEGORY_ORDER.get(value, 99)))

            repeat_rows = fetch(tracks["repeatmasker"], chrom, start, end)
            repeats = {
                f"{item[3]}/{item[4]}" for item in repeat_rows
                if len(item) > 4 and item[3] not in {"", "."}
            }
            chromhmm_rows = fetch(tracks["chromhmm"], chrom, start, end)
            states = {item[3] for item in chromhmm_rows if len(item) > 3}
            segdup = bool(fetch(tracks["segdup"], chrom, start, end))
            unique_rows = fetch(tracks["unique"], chrom, start, end)
            unique_fraction = union_overlap(unique_rows, start, end) / (end - start)
            multi_mean = multi_mappability.get(row["row_id"], 0.0)
            row.update({
                "gene": ";".join(sorted(genes)) or ".",
                "GENCODE_category": category_text,
                "ChromHMM_state": ";".join(sorted(states)) or "unavailable",
                "repeat_type": ";".join(sorted(repeats)) or "none",
                "segmental_duplication": "yes" if segdup else "no",
                "mappability": "high" if unique_fraction >= 0.5 else "low",
                "mappability_multi_k100_mean": f"{multi_mean:.6f}",
                "mappability_unique_k100_fraction": f"{unique_fraction:.6f}",
            })
            writer.writerow(row)

    for tabix in tracks.values():
        tabix.close()


if __name__ == "__main__":
    main()
