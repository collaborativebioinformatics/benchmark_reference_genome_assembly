#!/usr/bin/env python3

import argparse
import bisect
import gzip
from collections import defaultdict
from pathlib import Path

import pysam


BUILDS = ("hg15", "hg16", "hg17", "hg18", "hg19", "hg38", "hs1")


def scalar(value):
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def variant_length(record):
    value = scalar(record.info.get("SVLEN"))
    try:
        return abs(int(value))
    except (TypeError, ValueError):
        return max(1, record.stop - record.start)


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


def point_contained(index, chrom, position):
    if chrom not in index:
        return False
    starts, intervals = index[chrom]
    offset = bisect.bisect_right(starts, position) - 1
    return offset >= 0 and intervals[offset][0] <= position < intervals[offset][1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--arm", required=True, choices=("minimap2", "winnowmap"))
    parser.add_argument("--eligibility", choices=("full", "endpoints"), default="full")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    mask = read_mask(root / "results/masks/common-7way.hs1.bed.gz")

    first = root / f"results/lifted/{args.arm}/hg15/{args.arm}.hg15.to-hs1.vcf.gz"
    with pysam.VariantFile(str(first)) as source:
        header = source.header.copy()
    if "ID" not in header.info:
        header.add_line(
            '##INFO=<ID=ID,Number=.,Type=String,Description="Input call identifier retained by source VCF">'
        )
    if "SOURCE_BUILD" not in header.info:
        header.add_line(
            '##INFO=<ID=SOURCE_BUILD,Number=1,Type=String,Description="Source reference assembly before hs1 projection">'
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with pysam.VariantFile(str(args.output), "w", header=header) as output:
        for build in BUILDS:
            path = root / f"results/lifted/{args.arm}/{build}/{args.arm}.{build}.to-hs1.vcf.gz"
            with pysam.VariantFile(str(path)) as source:
                for index, record in enumerate(source, 1):
                    start, end = record.start, max(record.stop, record.start + 1)
                    svtype = str(scalar(record.info.get("SVTYPE", "MISSING")))
                    if svtype != "BND" and variant_length(record) < 50:
                        continue
                    eligible = (
                        contained(mask, record.contig, start, end)
                        if args.eligibility == "full" else
                        point_contained(mask, record.contig, start)
                        and point_contained(mask, record.contig, end - 1)
                    )
                    if not eligible:
                        continue
                    if "ID" in record.info:
                        del record.info["ID"]
                    record.translate(header)
                    original_id = record.id if record.id not in {None, "."} else str(index)
                    record.id = f"{build}::{original_id}"
                    record.info["SOURCE_BUILD"] = build
                    output.write(record)


if __name__ == "__main__":
    main()
