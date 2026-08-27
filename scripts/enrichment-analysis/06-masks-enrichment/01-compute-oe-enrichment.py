#!/usr/bin/env python3

import argparse
import bisect
import csv
import gzip
import math
from collections import Counter, defaultdict
from pathlib import Path


BUILDS = ("hg15", "hg16", "hg17", "hg18", "hg19", "hg38", "hs1")


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def merge(intervals):
    result = {}
    for chrom, values in intervals.items():
        out = []
        for start, end in sorted(values):
            if start >= end:
                continue
            if out and start <= out[-1][1]:
                out[-1][1] = max(out[-1][1], end)
            else:
                out.append([start, end])
        result[chrom] = [(start, end) for start, end in out]
    return result


def read_bed(path, category_column=None, category_mapper=None):
    grouped = defaultdict(lambda: defaultdict(list))
    with opener(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            category = "covered" if category_column is None else fields[category_column]
            if category_mapper is not None:
                category = category_mapper(category)
            grouped[category][fields[0]].append((int(fields[1]), int(fields[2])))
    return {category: merge(values) for category, values in grouped.items()}


def total_bases(intervals):
    return sum(end - start for values in intervals.values() for start, end in values)


def intersection_bases(left, right):
    total = 0
    for chrom in left.keys() & right.keys():
        a, b = left[chrom], right[chrom]
        i = j = 0
        while i < len(a) and j < len(b):
            total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
    return total


def union_categories(categories):
    combined = defaultdict(list)
    for intervals in categories.values():
        for chrom, values in intervals.items():
            combined[chrom].extend(values)
    return merge(combined)


def containment_index(intervals):
    return {
        chrom: ([start for start, _ in values], values)
        for chrom, values in intervals.items()
    }


def contained(index, chrom, start, end):
    if chrom not in index:
        return False
    starts, values = index[chrom]
    position = bisect.bisect_right(starts, start) - 1
    return position >= 0 and values[position][0] <= start and end <= values[position][1]


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
    gencode = set(row["GENCODE_category"].split(";"))
    chromhmm = set(row["ChromHMM_state"].split(";"))
    if row["repeat_type"] == "none":
        repeats = {"none"}
    else:
        repeats = {repeat_category(item.split("/", 1)[0]) for item in row["repeat_type"].split(";")}
    return {
        "GENCODE": gencode,
        "ChromHMM": chromhmm,
        "RepeatMasker": repeats,
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


def format_number(value):
    if isinstance(value, float) and math.isnan(value):
        return "."
    return f"{value:.10g}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    annotation_dir = root / "results/annotations/hs1"
    mask_dir = root / "results/masks"
    output_dir = root / "results/enrichment"
    output_dir.mkdir(parents=True, exist_ok=True)

    masks = {}
    common = read_bed(mask_dir / "common-7way.hs1.bed.gz")["covered"]
    for build in BUILDS:
        masks[("common_7way", build)] = common
        masks[("source_specific", build)] = read_bed(
            mask_dir / f"eligible/{build}.hs1.bed.gz"
        )["covered"]

    features = read_bed(annotation_dir / "features.bed.gz", 3)
    repeatmasker = read_bed(annotation_dir / "repeatmasker.bed.gz", 3, repeat_category)
    chromhmm = read_bed(
        root / "results/annotations/chromhmm/hs1/E116_GM12878_15state.hs1.bed.gz", 3
    )
    segdup = read_bed(annotation_dir / "segdup.bed.gz")
    mappability = read_bed(annotation_dir / "mappability-unique-k100.bed.gz")

    tracks = {
        "GENCODE": features,
        "RepeatMasker": repeatmasker,
        "ChromHMM": chromhmm,
        "segmental_duplication": {"yes": segdup["covered"]},
        "mappability_k100": {"high": mappability["covered"]},
    }
    complement_names = {
        "GENCODE": "intergenic",
        "RepeatMasker": "none",
        "ChromHMM": "unavailable",
        "segmental_duplication": "no",
        "mappability_k100": "low",
    }

    expected = {}
    mask_sizes = {}
    for key, mask in masks.items():
        eligible_bases = total_bases(mask)
        mask_sizes[key] = eligible_bases
        for layer, categories in tracks.items():
            for category, intervals in categories.items():
                expected[(key, layer, category)] = intersection_bases(mask, intervals)
            covered = intersection_bases(mask, union_categories(categories))
            expected[(key, layer, complement_names[layer])] = eligible_bases - covered

    counts = Counter()
    totals = Counter()
    raw_totals = Counter()
    indices = {key: containment_index(mask) for key, mask in masks.items()}
    table = root / "results/annotations/sv-annotations.tsv.gz"
    with gzip.open(table, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            arm, build, status = row["arm"], row["reference"], row["benchmark_status"]
            base_group = (arm, build, status)
            raw_totals[base_group] += 1
            start, end = int(row["start"]), int(row["end"])
            for mode in ("common_7way", "source_specific"):
                mask_key = (mode, build)
                if not contained(indices[mask_key], row["chrom"], start, end):
                    continue
                group = (mode, arm, build, status)
                totals[group] += 1
                for layer, categories in observed_categories(row).items():
                    for category in categories:
                        counts[(*group, layer, category)] += 1

    with open(output_dir / "eligibility-summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "analysis_mask", "arm", "reference", "benchmark_status", "mask_bases",
            "all_variants", "eligible_variants", "retained_fraction",
        ])
        for mode in ("common_7way", "source_specific"):
            for arm in ("minimap2", "winnowmap"):
                for build in BUILDS:
                    for status in ("TP", "FP", "FN"):
                        raw = raw_totals[(arm, build, status)]
                        kept = totals[(mode, arm, build, status)]
                        writer.writerow([
                            mode, arm, build, status, mask_sizes[(mode, build)], raw, kept,
                            f"{kept / raw:.8f}" if raw else ".",
                        ])

    fieldnames = [
        "analysis_mask", "arm", "reference", "benchmark_status", "annotation_layer",
        "category", "observed_count", "eligible_variant_count", "observed_fraction",
        "expected_bases", "eligible_mask_bases", "expected_fraction", "oe_enrichment",
        "oe_ci95_low", "oe_ci95_high",
    ]
    with open(output_dir / "oe-enrichment.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for mode in ("common_7way", "source_specific"):
            for arm in ("minimap2", "winnowmap"):
                for build in BUILDS:
                    mask_key = (mode, build)
                    eligible_bases = mask_sizes[mask_key]
                    for status in ("TP", "FP", "FN"):
                        group = (mode, arm, build, status)
                        total = totals[group]
                        for layer, categories in tracks.items():
                            names = set(categories) | {complement_names[layer]}
                            names |= {
                                key[-1] for key in counts
                                if key[:4] == group and key[4] == layer
                            }
                            for category in sorted(names):
                                count = counts[(*group, layer, category)]
                                observed_fraction = count / total if total else math.nan
                                expected_bases = expected.get((mask_key, layer, category), 0)
                                expected_fraction = expected_bases / eligible_bases if eligible_bases else math.nan
                                enrichment = (
                                    observed_fraction / expected_fraction
                                    if total and expected_fraction > 0 else math.nan
                                )
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
                                    "observed_fraction": format_number(observed_fraction),
                                    "expected_bases": expected_bases,
                                    "eligible_mask_bases": eligible_bases,
                                    "expected_fraction": format_number(expected_fraction),
                                    "oe_enrichment": format_number(enrichment),
                                    "oe_ci95_low": format_number(low / expected_fraction if expected_fraction > 0 else math.nan),
                                    "oe_ci95_high": format_number(high / expected_fraction if expected_fraction > 0 else math.nan),
                                })

    with open(output_dir / "definitions.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["item", "definition"])
        writer.writerows([
            ("common_7way", "GIAB hs1 high-confidence BED intersected with hs1 chain blocks liftable from hg15, hg16, hg17, hg18, hg19, and hg38"),
            ("source_specific", "GIAB hs1 high-confidence BED intersected with hs1 chain blocks liftable from the named source; hs1 uses the full GIAB BED"),
            ("variant_eligibility", "Entire event interval must be contained in one mask interval; insertions use their one-base anchor"),
            ("oe_enrichment", "fraction of eligible variants overlapping a category divided by fraction of eligible mask bases covered by that category"),
            ("confidence_interval", "95% Wilson interval for the observed variant fraction divided by the fixed expected base fraction"),
            ("overlapping_categories", "GENCODE, RepeatMasker, and ChromHMM categories are tested independently and may overlap for one SV"),
            ("mappability", "High when at least half of the event interval is covered by the hs1 Umap k100 unique track; expected high is its base coverage"),
            ("ChromHMM_unavailable", "Eligible hs1 bases not covered by accepted lifted E116 GM12878 GRCh38 15-state bins"),
        ])


if __name__ == "__main__":
    main()
