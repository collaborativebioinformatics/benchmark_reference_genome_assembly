#!/usr/bin/env python3

import argparse
import bisect
import csv
import gzip
import math
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


def point_contained(index, chrom, position):
    if chrom not in index:
        return False
    starts, intervals = index[chrom]
    offset = bisect.bisect_right(starts, position) - 1
    return offset >= 0 and intervals[offset][0] <= position < intervals[offset][1]


def full_contained(index, chrom, start, end):
    if chrom not in index:
        return False
    starts, intervals = index[chrom]
    offset = bisect.bisect_right(starts, start) - 1
    return offset >= 0 and intervals[offset][0] <= start and end <= intervals[offset][1]


def endpoints_contained(index, chrom, start, end):
    return point_contained(index, chrom, start) and point_contained(index, chrom, end - 1)


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


def observed_categories(row):
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output_dir = root / "results/enrichment"

    masks = {}
    common = read_mask(root / "results/masks/common-7way.hs1.bed.gz")
    truth = read_mask(root / "results/masks/eligible/hs1.hs1.bed.gz")
    for build in BUILDS:
        masks[("common_endpoints", build)] = common
        masks[("source_endpoints", build)] = read_mask(root / f"results/masks/eligible/{build}.hs1.bed.gz")
        masks[("giab_full", build)] = truth

    expected = {}
    with open(output_dir / "oe-enrichment.tsv", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source_mode = row["analysis_mask"]
            key = (
                source_mode, row["reference"], row["annotation_layer"], row["category"]
            )
            expected[key] = (
                int(row["expected_bases"]), int(row["eligible_mask_bases"]),
                float(row["expected_fraction"]),
            )

    counts = Counter()
    totals = Counter()
    raw_totals = Counter()
    with gzip.open(root / "results/annotations/sv-annotations.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            arm, build, status = row["arm"], row["reference"], row["benchmark_status"]
            raw_totals[(arm, build, status)] += 1
            start, end = int(row["start"]), int(row["end"])
            for mode in ("common_endpoints", "source_endpoints", "giab_full"):
                index = masks[(mode, build)]
                eligible = (
                    full_contained(index, row["chrom"], start, end)
                    if mode == "giab_full" else
                    endpoints_contained(index, row["chrom"], start, end)
                )
                if not eligible:
                    continue
                group = (mode, arm, build, status)
                totals[group] += 1
                for layer, categories in observed_categories(row).items():
                    for category in categories:
                        counts[(*group, layer, category)] += 1


    with open(output_dir / "relaxed-eligibility-summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "analysis_mask", "arm", "reference", "benchmark_status", "all_variants",
            "eligible_variants", "retained_fraction",
        ])
        for mode in ("common_endpoints", "source_endpoints", "giab_full"):
            for arm in ("minimap2", "winnowmap"):
                for build in BUILDS:
                    for status in ("TP", "FP", "FN"):
                        raw = raw_totals[(arm, build, status)]
                        kept = totals[(mode, arm, build, status)]
                        writer.writerow([
                            mode, arm, build, status, raw, kept,
                            f"{kept / raw:.8f}" if raw else ".",
                        ])

    fields = [
        "analysis_mask", "arm", "reference", "benchmark_status", "annotation_layer",
        "category", "observed_count", "eligible_variant_count", "observed_fraction",
        "expected_bases", "eligible_mask_bases", "expected_fraction", "oe_enrichment",
        "oe_ci95_low", "oe_ci95_high",
    ]
    with open(output_dir / "relaxed-oe-enrichment.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for mode in ("common_endpoints", "source_endpoints", "giab_full"):
            for arm in ("minimap2", "winnowmap"):
                for build in BUILDS:
                    source_mode = "common_7way" if mode == "common_endpoints" else "source_specific"
                    expected_build = "hs1" if mode == "giab_full" else build
                    keys = sorted(
                        key for key in expected
                        if key[0] == source_mode and key[1] == expected_build
                    )
                    for status in ("TP", "FP", "FN"):
                        group = (mode, arm, build, status)
                        total = totals[group]
                        for _, _, layer, category in keys:
                            expected_count, mask_count, expected_fraction = expected[
                                (source_mode, expected_build, layer, category)
                            ]
                            count = counts[(*group, layer, category)]
                            observed = count / total if total else math.nan
                            enrichment = observed / expected_fraction if expected_fraction > 0 and total else math.nan
                            low, high = wilson(count, total)
                            writer.writerow({
                                "analysis_mask": mode,
                                "arm": arm,
                                "reference": build,
                                "benchmark_status": status,
                                "annotation_layer": layer,
                                "category": category,
                                "observed_count": count,
                                "eligible_variant_count": total,
                                "observed_fraction": fmt(observed),
                                "expected_bases": expected_count,
                                "eligible_mask_bases": mask_count,
                                "expected_fraction": fmt(expected_fraction),
                                "oe_enrichment": fmt(enrichment),
                                "oe_ci95_low": fmt(low / expected_fraction if expected_fraction > 0 else math.nan),
                                "oe_ci95_high": fmt(high / expected_fraction if expected_fraction > 0 else math.nan),
                            })

    with open(output_dir / "relaxed-definitions.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["analysis_mask", "variant_eligibility", "expected_base_denominator"])
        writer.writerows([
            ("common_endpoints", "Both represented event endpoints covered by the common 7-way mask; insertion/BND uses represented one-base anchor", "Common 7-way mask"),
            ("source_endpoints", "Both represented event endpoints covered by the source-specific mask; insertion/BND uses represented one-base anchor", "Named source-specific mask"),
            ("giab_full", "Entire represented event interval contained in one GIAB hs1 high-confidence interval", "Full GIAB hs1 benchmark BED"),
        ])


if __name__ == "__main__":
    main()
