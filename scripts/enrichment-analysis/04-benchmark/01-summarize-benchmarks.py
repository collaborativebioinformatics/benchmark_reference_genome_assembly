#!/usr/bin/env python3

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import pysam


SIZE_BINS = [
    (50, 100, "50-99"),
    (100, 500, "100-499"),
    (500, 1000, "500-999"),
    (1000, 5000, "1k-4.9k"),
    (5000, 10000, "5k-9.9k"),
    (10000, 50000, "10k-49.9k"),
    (50000, None, ">=50k"),
]


def variant_size(record):
    value = record.info.get("SVLEN")
    if isinstance(value, tuple):
        value = value[0] if value else None
    try:
        return abs(int(value))
    except (TypeError, ValueError):
        pass
    if record.alts and len(record.alts) == 1 and not record.alts[0].startswith("<"):
        return abs(len(record.alts[0]) - len(record.ref))
    return None


def size_bin(size):
    if size is None:
        return "unknown"
    for lower, upper, label in SIZE_BINS:
        if size >= lower and (upper is None or size < upper):
            return label
    return "<50"


def count_vcf(path):
    by_type = Counter()
    by_size = Counter()
    by_type_size = Counter()
    with pysam.VariantFile(str(path)) as vcf:
        for record in vcf:
            value = record.info.get("SVTYPE", "MISSING")
            if isinstance(value, tuple):
                value = value[0] if value else "MISSING"
            svtype = str(value)
            label = size_bin(variant_size(record))
            by_type[svtype] += 1
            by_size[label] += 1
            by_type_size[(svtype, label)] += 1
    return by_type, by_size, by_type_size


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def metrics(tp_comp, fp, tp_base, fn):
    precision = ratio(tp_comp, tp_comp + fp)
    recall = ratio(tp_base, tp_base + fn)
    if precision is None or recall is None or precision + recall == 0:
        f1 = None
    else:
        f1 = 2 * precision * recall / (precision + recall)
    return precision, recall, f1


def write_rows(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--callsets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    overall = []
    stratified = {"svtype": [], "size": [], "svtype_size": []}
    with open(args.callsets, newline="") as handle:
        callsets = list(csv.DictReader(handle, delimiter="\t"))

    for row in callsets:
        arm, build = row["arm"], row["build"]
        result = args.root / "results/TRUVARI" / arm / build
        summary_path = result / "summary.json"
        if not summary_path.exists():
            continue
        summary = json.loads(summary_path.read_text())
        overall.append({
            "arm": arm,
            "build": build,
            "sniffles_version": row["sniffles_version"],
            "TP_base": summary["TP-base"],
            "TP_comp": summary["TP-comp"],
            "FP": summary["FP"],
            "FN": summary["FN"],
            "precision": summary["precision"],
            "recall": summary["recall"],
            "f1": summary["f1"],
            "gt_concordance": summary.get("gt_concordance"),
            "base_count": summary["base cnt"],
            "comp_count": summary["comp cnt"],
        })

        files = {
            "tp_comp": result / "tp-comp.vcf.gz",
            "fp": result / "fp.vcf.gz",
            "tp_base": result / "tp-base.vcf.gz",
            "fn": result / "fn.vcf.gz",
        }
        counts = {name: count_vcf(path) for name, path in files.items()}
        dimensions = [
            ("svtype", 0),
            ("size", 1),
            ("svtype_size", 2),
        ]
        for dimension, index in dimensions:
            keys = set().union(*(value[index].keys() for value in counts.values()))
            for key in sorted(keys):
                tp_comp = counts["tp_comp"][index][key]
                fp = counts["fp"][index][key]
                tp_base = counts["tp_base"][index][key]
                fn = counts["fn"][index][key]
                precision, recall, f1 = metrics(tp_comp, fp, tp_base, fn)
                result_row = {
                    "arm": arm,
                    "build": build,
                    "TP_base": tp_base,
                    "TP_comp": tp_comp,
                    "FP": fp,
                    "FN": fn,
                    "precision": precision,
                    "recall": recall,
                    "f1": f1,
                }
                if dimension == "svtype":
                    result_row["SVTYPE"] = key
                elif dimension == "size":
                    result_row["size_bin"] = key
                else:
                    result_row["SVTYPE"], result_row["size_bin"] = key
                stratified[dimension].append(result_row)

    common = ["arm", "build", "TP_base", "TP_comp", "FP", "FN", "precision", "recall", "f1"]
    write_rows(
        args.output_dir / "overall.tsv",
        overall,
        ["arm", "build", "sniffles_version", "TP_base", "TP_comp", "FP", "FN",
         "precision", "recall", "f1", "gt_concordance", "base_count", "comp_count"],
    )
    write_rows(args.output_dir / "by-svtype.tsv", stratified["svtype"], ["arm", "build", "SVTYPE", *common[2:]])
    write_rows(args.output_dir / "by-size.tsv", stratified["size"], ["arm", "build", "size_bin", *common[2:]])
    write_rows(
        args.output_dir / "by-svtype-size.tsv",
        stratified["svtype_size"],
        ["arm", "build", "SVTYPE", "size_bin", *common[2:]],
    )


if __name__ == "__main__":
    main()

