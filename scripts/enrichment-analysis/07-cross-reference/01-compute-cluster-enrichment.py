#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
from collections import Counter, defaultdict
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
    current = None
    for left, right in intervals:
        if left >= right:
            continue
        if current is None:
            current = [left, right]
        elif left > current[1]:
            total += current[1] - current[0]
            current = [left, right]
        else:
            current[1] = max(current[1], right)
    if current is not None:
        total += current[1] - current[0]
    return total

def repeat_category(value):
    normalized = value.lower().replace("-", "_").replace(" ", "_")
    if normalized == "line":
        return "LINE"
    if normalized == "sine":
        return "SINE"
    if normalized == "ltr":
        return "LTR"
    if "satellite" in normalized:
        return "satellite"
    if normalized == "simple_repeat":
        return "simple_repeat"
    if normalized == "low_complexity":
        return "low_complexity"
    return "other_repeat"


def categories(row):
    repeats = (
        {"none"} if row["repeat_type"] == "none" else
        {repeat_category(item.split("/", 1)[0]) for item in row["repeat_type"].split(";")}
    )
    return {
        "GENCODE": set(row["GENCODE_category"].split(";")),
        "RepeatMasker": repeats,
        "ChromHMM": set(row["ChromHMM_state"].split(";")),
        "segmental_duplication": {row["segmental_duplication"]},
        "mappability_k100": {row["mappability"]},
    }


def wilson(count, total, z=1.959963984540054):
    if total == 0:
        return math.nan, math.nan
    proportion = count / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    margin = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, center - margin), min(1.0, center + margin)


def fmt(value):
    return "." if math.isnan(value) else f"{value:.10g}"


def sharing_class(count):
    if count == 1:
        return "reference_specific"
    if count == 7:
        return "shared_all_7"
    return "shared_2_to_6"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument(
        "--analysis", choices=("strict_full", "common_endpoints"), default="strict_full"
    )
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = root / "results/cross-reference"
    cluster_dir = (
        output_dir if args.analysis == "strict_full" else output_dir / "common-endpoints"
    )

    expected = {}
    expected_bases = {}
    mask_bases = None
    with open(root / "results/enrichment/oe-enrichment.tsv", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["analysis_mask"] != "common_7way":
                continue
            key = (row["annotation_layer"], row["category"])
            expected[key] = float(row["expected_fraction"])
            expected_bases[key] = int(row["expected_bases"])
            mask_bases = int(row["eligible_mask_bases"])

    annotation_dir = root / "results/annotations/hs1"
    tracks = {
        "genes": pysam.TabixFile(str(annotation_dir / "genes.bed.gz")),
        "features": pysam.TabixFile(str(annotation_dir / "features.bed.gz")),
        "repeatmasker": pysam.TabixFile(str(annotation_dir / "repeatmasker.bed.gz")),
        "segdup": pysam.TabixFile(str(annotation_dir / "segdup.bed.gz")),
        "unique": pysam.TabixFile(str(annotation_dir / "mappability-unique-k100.bed.gz")),
        "chromhmm": pysam.TabixFile(str(
            root / "results/annotations/chromhmm/hs1/E116_GM12878_15state.hs1.bed.gz"
        )),
    }

    annotated_path = cluster_dir / "cluster-annotations.tsv.gz"
    counts = Counter()
    totals = Counter()
    annotation_fields = [
        "arm", "cluster_id", "representative_source", "SV_ID", "chrom", "start", "end",
        "SV_type", "SV_length", "reference_count", "call_count", "references",
        "sharing_class", "sharing_count_class", "gene", "GENCODE_category", "ChromHMM_state", "repeat_type",
        "segmental_duplication", "mappability", "mappability_unique_k100_fraction",
    ]
    with gzip.open(annotated_path, "wt", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=annotation_fields, delimiter="\t")
        writer.writeheader()
        for arm in ("minimap2", "winnowmap"):
            with open(cluster_dir / arm / "clusters.tsv", newline="") as handle:
                for cluster in csv.DictReader(handle, delimiter="\t"):
                    representative = cluster["member_ids"].split(";", 1)[0]
                    build, original_id = representative.split("::", 1)
                    chrom = cluster["chrom"]
                    start, end = int(cluster["start"]), int(cluster["end"])
                    gene_rows = fetch(tracks["genes"], chrom, start, end)
                    feature_rows = fetch(tracks["features"], chrom, start, end)
                    genes = {item[3] for item in gene_rows if len(item) > 3 and item[3] != "."}
                    genes.update(item[4] for item in feature_rows if len(item) > 4 and item[4] != ".")
                    feature_categories = {item[3] for item in feature_rows if len(item) > 3}
                    if not feature_categories:
                        feature_categories = {"intron"} if gene_rows else {"intergenic"}
                    repeat_rows = fetch(tracks["repeatmasker"], chrom, start, end)
                    repeats = {
                        f"{item[3]}/{item[4]}" for item in repeat_rows
                        if len(item) > 4 and item[3] not in {"", "."}
                    }
                    chromhmm_rows = fetch(tracks["chromhmm"], chrom, start, end)
                    states = {item[3] for item in chromhmm_rows if len(item) > 3}
                    unique_rows = fetch(tracks["unique"], chrom, start, end)
                    unique_fraction = union_overlap(unique_rows, start, end) / (end - start)
                    count = int(cluster["reference_count"])
                    label = sharing_class(count)
                    exact_label = f"detected_{count}_of_7"
                    row = {
                        "arm": arm,
                        "cluster_id": cluster["cluster_id"],
                        "representative_source": build,
                        "SV_ID": original_id,
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "SV_type": cluster["SV_type"],
                        "SV_length": cluster["SV_length"],
                        "reference_count": count,
                        "call_count": cluster["call_count"],
                        "references": cluster["references"],
                        "sharing_class": label,
                        "sharing_count_class": exact_label,
                        "gene": ";".join(sorted(genes)) or ".",
                        "GENCODE_category": ";".join(sorted(feature_categories, key=lambda value: CATEGORY_ORDER.get(value, 99))),
                        "ChromHMM_state": ";".join(sorted(states)) or "unavailable",
                        "repeat_type": ";".join(sorted(repeats)) or "none",
                        "segmental_duplication": "yes" if fetch(tracks["segdup"], chrom, start, end) else "no",
                        "mappability": "high" if unique_fraction >= 0.5 else "low",
                        "mappability_unique_k100_fraction": f"{unique_fraction:.6f}",
                    }
                    writer.writerow(row)
                    for group in (label, exact_label, "all_clusters"):
                        totals[(arm, group)] += 1
                        for layer, values in categories(row).items():
                            for category in values:
                                counts[(arm, group, layer, category)] += 1

    for tabix in tracks.values():
        tabix.close()

    fields = [
        "analysis_mask", "arm", "sharing_class", "annotation_layer", "category",
        "observed_count", "cluster_count", "observed_fraction", "expected_bases",
        "eligible_mask_bases", "expected_fraction", "oe_enrichment", "oe_ci95_low",
        "oe_ci95_high",
    ]
    with open(cluster_dir / "cluster-oe-enrichment.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        layers = sorted({layer for layer, _ in expected})
        for arm in ("minimap2", "winnowmap"):
            groups = [f"detected_{count}_of_7" for count in range(1, 8)]
            groups.extend(("reference_specific", "shared_2_to_6", "shared_all_7", "all_clusters"))
            for group in groups:
                total = totals[(arm, group)]
                for layer in layers:
                    for category in sorted(value for key, value in expected if key == layer):
                        count = counts[(arm, group, layer, category)]
                        observed = count / total if total else math.nan
                        expected_fraction = expected[(layer, category)]
                        enrichment = observed / expected_fraction if expected_fraction > 0 and total else math.nan
                        low, high = wilson(count, total)
                        writer.writerow({
                            "analysis_mask": args.analysis,
                            "arm": arm,
                            "sharing_class": group,
                            "annotation_layer": layer,
                            "category": category,
                            "observed_count": count,
                            "cluster_count": total,
                            "observed_fraction": fmt(observed),
                            "expected_bases": expected_bases[(layer, category)],
                            "eligible_mask_bases": mask_bases,
                            "expected_fraction": fmt(expected_fraction),
                            "oe_enrichment": fmt(enrichment),
                            "oe_ci95_low": fmt(low / expected_fraction if expected_fraction > 0 else math.nan),
                            "oe_ci95_high": fmt(high / expected_fraction if expected_fraction > 0 else math.nan),
                        })


if __name__ == "__main__":
    main()
