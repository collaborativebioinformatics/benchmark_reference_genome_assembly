#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def count_rows(path, compressed=False):
    open_fn = gzip.open if compressed else open
    with open_fn(path, "rt") as handle:
        return sum(1 for _ in handle) - 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    checks = []

    def check(name, passed, detail):
        checks.append((name, "PASS" if passed else "FAIL", detail))

    with open(root / "state/vcf-input-audit.tsv") as handle:
        callsets = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["role"] == "callset"
        ]
    versions = sorted({row["source"] for row in callsets})
    check(
        "uniform_sniffles_version",
        len(callsets) == 14 and versions == ["Sniffles2_2.8.0"],
        f"callsets={len(callsets)};sources={','.join(versions)}",
    )

    annotation = root / "results/annotations/sv-annotations.tsv.gz"
    annotation_rows = count_rows(annotation, compressed=True)
    with open(root / "results/annotations/annotation-row-counts.tsv") as handle:
        expected_annotation_rows = sum(int(row["rows"]) for row in csv.DictReader(handle, delimiter="\t"))
    check("annotation_row_accounting", annotation_rows == expected_annotation_rows,
          f"combined={annotation_rows};summary={expected_annotation_rows}")

    with open(root / "results/masks/mask-summary.tsv") as handle:
        masks = {row["mask"]: row for row in csv.DictReader(handle, delimiter="\t")}
    common_fraction = float(masks["common_7way"]["eligible_truth_fraction"])
    check("common_mask_fraction", 0 < common_fraction <= 1,
          f"fraction={common_fraction:.8f};bases={masks['common_7way']['eligible_truth_bases']}")
    required_masks = [root / "results/masks/common-7way.hs1.bed.gz"] + [
        root / f"results/masks/eligible/{build}.hs1.bed.gz"
        for build in ("hg15", "hg16", "hg17", "hg18", "hg19", "hg38", "hs1")
    ]
    missing_masks = [str(path) for path in required_masks if not path.is_file() or not Path(f"{path}.tbi").is_file()]
    check("mask_files_and_indexes", not missing_masks,
          "all present" if not missing_masks else ";".join(missing_masks))

    for name in ("oe-enrichment.tsv", "relaxed-oe-enrichment.tsv"):
        rows = 0
        invalid = 0
        with open(root / "results/enrichment" / name) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows += 1
                observed = float(row["observed_fraction"])
                expected = float(row["expected_fraction"])
                enrichment = float(row["oe_enrichment"])
                if not (0 <= observed <= 1 and 0 <= expected <= 1 and enrichment >= 0):
                    invalid += 1
        check(f"{name}_numeric_ranges", invalid == 0, f"rows={rows};invalid={invalid}")

    total_cluster_annotations = 0
    for mode, base in (
        ("strict", root / "results/cross-reference"),
        ("common_endpoints", root / "results/cross-reference/common-endpoints"),
    ):
        mode_clusters = 0
        for arm in ("minimap2", "winnowmap"):
            cluster_rows = count_rows(base / arm / "clusters.tsv")
            with open(base / arm / "reference-count-distribution.tsv") as handle:
                distribution_rows = sum(int(row["clusters"]) for row in csv.DictReader(handle, delimiter="\t"))
            check(f"{mode}_{arm}_cluster_distribution", cluster_rows == distribution_rows,
                  f"clusters={cluster_rows};distribution={distribution_rows}")
            mode_clusters += cluster_rows
        annotation_rows = count_rows(base / "cluster-annotations.tsv.gz", compressed=True)
        check(f"{mode}_cluster_annotations", mode_clusters == annotation_rows,
              f"clusters={mode_clusters};annotations={annotation_rows}")
        exact_groups = defaultdict(set)
        with open(base / "cluster-oe-enrichment.tsv") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["sharing_class"].startswith("detected_"):
                    exact_groups[row["arm"]].add(row["sharing_class"])
        required_groups = {f"detected_{count}_of_7" for count in range(1, 8)}
        exact_ok = all(exact_groups[arm] == required_groups for arm in ("minimap2", "winnowmap"))
        check(f"{mode}_exact_sharing_enrichment", exact_ok,
              ";".join(f"{arm}={len(exact_groups[arm])}" for arm in ("minimap2", "winnowmap")))
        total_cluster_annotations += annotation_rows

    figures = root / "results/figures"
    sets = {extension: {path.stem for path in figures.glob(f"*.{extension}")} for extension in ("png", "pdf", "svg")}
    figure_sets_equal = sets["png"] == sets["pdf"] == sets["svg"]
    check("figure_triplets", figure_sets_equal and len(sets["png"]) == 16,
          f"png={len(sets['png'])};pdf={len(sets['pdf'])};svg={len(sets['svg'])}")

    output = root / "results/AUDIT.tsv"
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["check", "status", "detail"])
        writer.writerows(checks)
    failures = [name for name, status, _ in checks if status != "PASS"]
    if failures:
        raise SystemExit("Audit failures: " + ", ".join(failures))


if __name__ == "__main__":
    main()
