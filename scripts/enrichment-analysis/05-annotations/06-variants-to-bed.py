#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import pysam


def scalar(value):
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def variant_length(record):
    value = scalar(record.info.get("SVLEN"))
    try:
        return abs(int(value))
    except (TypeError, ValueError):
        return max(1, record.stop - record.pos)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--build", required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()

    fields = [
        "row_id", "arm", "reference", "coordinate_system", "benchmark_status",
        "SV_ID", "chrom", "start", "end", "POS", "VCF_END", "SV_type", "SV_length",
    ]
    with open(args.bed, "w") as bed, open(args.metadata, "w", newline="") as metadata:
        writer = csv.DictWriter(metadata, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for status in ("TP", "FP", "FN"):
            path = args.result_dir / f"{status}.vcf.gz"
            with pysam.VariantFile(str(path)) as vcf:
                for index, record in enumerate(vcf, 1):
                    row_id = f"{args.arm}.{args.build}.{status}.{index}"
                    start = record.start
                    end = max(record.stop, start + 1)
                    svtype = str(scalar(record.info.get("SVTYPE", "MISSING")))
                    bed.write(f"{record.contig}\t{start}\t{end}\t{row_id}\n")
                    writer.writerow({
                        "row_id": row_id,
                        "arm": args.arm,
                        "reference": args.build,
                        "coordinate_system": "hs1",
                        "benchmark_status": status,
                        "SV_ID": record.id or ".",
                        "chrom": record.contig,
                        "start": start,
                        "end": end,
                        "POS": record.pos,
                        "VCF_END": record.stop,
                        "SV_type": svtype,
                        "SV_length": variant_length(record),
                    })


if __name__ == "__main__":
    main()
