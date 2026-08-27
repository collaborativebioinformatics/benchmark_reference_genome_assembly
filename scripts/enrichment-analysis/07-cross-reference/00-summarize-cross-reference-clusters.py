#!/usr/bin/env python3

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def scalar(value):
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def sv_length(record):
    value = scalar(record.info.get("SVLEN"))
    try:
        return abs(int(value))
    except (TypeError, ValueError):
        return max(1, record.stop - record.start)


def collapse_id(record):
    value = record.info.get("MatchId")
    if isinstance(value, tuple):
        return str(value[0])
    return str(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", required=True)
    parser.add_argument("--kept", type=Path, required=True)
    parser.add_argument("--removed", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--distribution", type=Path, required=True)
    args = parser.parse_args()

    members = defaultdict(list)
    representatives = {}
    with pysam.VariantFile(str(args.kept)) as vcf:
        for index, record in enumerate(vcf, 1):
            group = str(record.info.get("CollapseId", f"singleton.{index}"))
            build = str(record.info["SOURCE_BUILD"])
            members[group].append((build, record.id or "."))
            representatives[group] = record.copy()
    with pysam.VariantFile(str(args.removed)) as vcf:
        for record in vcf:
            group = collapse_id(record)
            members[group].append((str(record.info["SOURCE_BUILD"]), record.id or "."))

    args.clusters.parent.mkdir(parents=True, exist_ok=True)
    distribution = Counter()
    with open(args.clusters, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "arm", "cluster_id", "chrom", "start", "end", "SV_type", "SV_length",
            "reference_count", "call_count", "references", "member_ids",
        ])
        for group, record in sorted(
            representatives.items(), key=lambda item: (item[1].contig, item[1].start, item[1].stop)
        ):
            builds = sorted({build for build, _ in members[group]})
            distribution[len(builds)] += 1
            writer.writerow([
                args.arm, group, record.contig, record.start, max(record.stop, record.start + 1),
                scalar(record.info.get("SVTYPE", "MISSING")), sv_length(record), len(builds),
                len(members[group]), ";".join(builds),
                ";".join(identifier for _, identifier in members[group]),
            ])

    with open(args.distribution, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["arm", "references_detected", "clusters", "fraction"])
        total = sum(distribution.values())
        for count in range(1, 8):
            writer.writerow([
                args.arm, count, distribution[count],
                f"{distribution[count] / total:.8f}" if total else ".",
            ])


if __name__ == "__main__":
    main()
