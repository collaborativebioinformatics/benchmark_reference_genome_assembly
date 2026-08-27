#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import pysam


BUILDS = ("hg15", "hg16", "hg17", "hg18", "hg19", "hg38")


def merge(intervals):
    merged = {}
    for chrom, values in intervals.items():
        out = []
        for start, end in sorted(values):
            if start >= end:
                continue
            if out and start <= out[-1][1]:
                out[-1][1] = max(out[-1][1], end)
            else:
                out.append([start, end])
        merged[chrom] = [(start, end) for start, end in out]
    return merged


def read_bed(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    intervals = defaultdict(list)
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            intervals[fields[0]].append((int(fields[1]), int(fields[2])))
    return merge(intervals)


def chain_query_blocks(path):
    intervals = defaultdict(list)
    query_name = None
    query_size = None
    query_strand = None
    query_position = None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "chain":
                if len(fields) != 13:
                    raise ValueError(f"Malformed chain header in {path}: {line.rstrip()}")
                query_name = fields[7]
                query_size = int(fields[8])
                query_strand = fields[9]
                query_position = int(fields[10])
                continue
            if query_name is None or len(fields) not in {1, 3}:
                raise ValueError(f"Malformed chain block in {path}: {line.rstrip()}")
            size = int(fields[0])
            if query_strand == "+":
                start, end = query_position, query_position + size
            elif query_strand == "-":
                start = query_size - (query_position + size)
                end = query_size - query_position
            else:
                raise ValueError(f"Unknown query strand {query_strand} in {path}")
            intervals[query_name].append((start, end))
            if len(fields) == 3:
                query_position += size + int(fields[2])
            else:
                query_name = None
                query_size = query_strand = query_position = None
    return merge(intervals)


def intersect(left, right):
    result = defaultdict(list)
    for chrom in left.keys() & right.keys():
        a, b = left[chrom], right[chrom]
        i = j = 0
        while i < len(a) and j < len(b):
            start = max(a[i][0], b[j][0])
            end = min(a[i][1], b[j][1])
            if start < end:
                result[chrom].append((start, end))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
    return merge(result)


def bases(intervals):
    return sum(end - start for values in intervals.values() for start, end in values)


def write_bgzip(path, intervals, order):
    path.parent.mkdir(parents=True, exist_ok=True)
    plain = path.with_suffix("")
    chroms = sorted(intervals, key=lambda chrom: (order.get(chrom, 10**9), chrom))
    with open(plain, "w") as handle:
        for chrom in chroms:
            for start, end in intervals[chrom]:
                handle.write(f"{chrom}\t{start}\t{end}\n")
    pysam.tabix_compress(str(plain), str(path), force=True)
    pysam.tabix_index(str(path), preset="bed", force=True)
    plain.unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root / "results/masks"
    truth = read_bed(root / "inputs/truth/hs1/HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed")
    fai = root.parent / "references/hs1/hs1.fa.fai"
    order = {}
    with open(fai) as handle:
        for index, line in enumerate(handle):
            order[line.split("\t", 1)[0]] = index

    liftable = {}
    eligible = {"hs1": truth}
    rows = []
    truth_bases = bases(truth)
    for build in BUILDS:
        chain = root.parent / f"liftover/chains/{build}/hs1/{build}ToHs1.over.chain.gz"
        blocks = chain_query_blocks(chain)
        liftable[build] = blocks
        eligible[build] = intersect(truth, blocks)
        write_bgzip(output / f"liftable/{build}.hs1.bed.gz", blocks, order)
        write_bgzip(output / f"eligible/{build}.hs1.bed.gz", eligible[build], order)
        rows.append((build, bases(blocks), bases(eligible[build])))

    write_bgzip(output / "eligible/hs1.hs1.bed.gz", truth, order)
    common = truth
    for build in BUILDS:
        common = intersect(common, liftable[build])
    write_bgzip(output / "common-7way.hs1.bed.gz", common, order)

    with open(output / "mask-summary.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "mask", "liftable_hs1_bases", "eligible_truth_bases",
            "truth_bases", "eligible_truth_fraction",
        ])
        for build, lift_bases, eligible_bases in rows:
            writer.writerow([
                build, lift_bases, eligible_bases, truth_bases,
                f"{eligible_bases / truth_bases:.8f}",
            ])
        writer.writerow(["hs1", ".", truth_bases, truth_bases, "1.00000000"])
        common_bases = bases(common)
        writer.writerow([
            "common_7way", ".", common_bases, truth_bases,
            f"{common_bases / truth_bases:.8f}",
        ])


if __name__ == "__main__":
    main()
