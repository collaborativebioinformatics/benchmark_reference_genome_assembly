#!/usr/bin/env python3
"""
Cross-reference severe read-trace outcomes against reference assembly gaps.

Consumes read_trace_best.tsv from trace_sv_reads_part1_samechr.sh (no pandas needed).
For each hs1-supporting-read that is MATCH in at least one build but SEVERE
(OTHER_CHROM / NO_SV_CALL) in another, checks whether a query-build coordinate
falls in or near a known assembly gap for that build.

OTHER_CHROM rows carry a real query-build coordinate (a candidate SV was found,
just in the wrong place) and are checked directly. NO_SV_CALL rows have no
coordinate of their own in this build (no SV was called anywhere for that read),
so the hs1 anchor's own chrom:pos is used as an approximate stand-in -- it is
NOT a lifted-over or observed position, just the nearest thing available without
a BAM lookup, and is labeled "hs1_anchor_approx" in the output accordingly.
"""
import argparse
import bisect
import csv
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path

SEVERE = {"OTHER_CHROM", "NO_SV_CALL"}


def find_latest_trace_dir():
    candidates = sorted(glob.glob(str(Path.home() / "read_trace_*")))
    return Path(candidates[-1]) if candidates else None


def load_gap_index(bed_path: Path):
    by_chrom = defaultdict(list)
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return by_chrom
    with open(bed_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            gap_type = parts[3] if len(parts) > 3 else "gap"
            by_chrom[chrom].append((start, end, gap_type))
    for chrom, intervals in by_chrom.items():
        intervals.sort()
    return by_chrom


def nearest_gap(index, chrom, pos):
    intervals = index.get(chrom)
    if not intervals:
        return None
    starts = [iv[0] for iv in intervals]
    i = bisect.bisect_right(starts, pos)
    best_dist, best_iv = None, None
    for iv in (intervals[i] if i < len(intervals) else None,
               intervals[i - 1] if i > 0 else None):
        if iv is None:
            continue
        start, end, _ = iv
        if start <= pos <= end:
            dist = 0
        elif pos < start:
            dist = start - pos
        else:
            dist = pos - end
        if best_dist is None or dist < best_dist:
            best_dist, best_iv = dist, iv
    return best_dist, best_iv


def anchors_with_match(trace_tsv: Path):
    matched = set()
    with open(trace_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["classification"] == "MATCH":
                matched.add(row["anchor_id"])
    return matched


def main():
    ap = argparse.ArgumentParser(
        description="Flag whether reference-sensitive SV read-trace misses land near assembly gaps."
    )
    ap.add_argument("--trace-dir", default=None,
                     help="Directory containing read_trace_best.tsv (default: latest ~/read_trace_* dir)")
    ap.add_argument("--gaps-dir", default=str(Path.home() / "Desktop"),
                     help="Directory containing <build>.gaps.bed files (default: ~/Desktop)")
    ap.add_argument("--window", type=int, default=1000,
                     help="Max distance (bp) from a gap to still call it 'near_gap' (default: 1000)")
    ap.add_argument("--out", default=None,
                     help="Output TSV path (default: <trace-dir>/gap_check/gap_candidate_summary.tsv)")
    args = ap.parse_args()

    trace_dir = Path(args.trace_dir) if args.trace_dir else find_latest_trace_dir()
    if trace_dir is None or not trace_dir.exists():
        sys.exit("FATAL: no trace directory found/given (expected read_trace_best.tsv inside it)")
    trace_tsv = trace_dir / "read_trace_best.tsv"
    if not trace_tsv.exists():
        sys.exit(f"FATAL: {trace_tsv} not found")

    gaps_dir = Path(args.gaps_dir)
    out_path = Path(args.out) if args.out else trace_dir / "gap_check" / "gap_candidate_summary.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Trace dir: {trace_dir}")
    print(f"Gaps dir:  {gaps_dir}")
    print("Pass 1/2: finding anchors with at least one MATCH...")
    matched_somewhere = anchors_with_match(trace_tsv)
    print(f"  {len(matched_somewhere)} anchor SVs matched in >=1 build")

    gap_index_cache = {}

    def gap_index_for(build):
        if build not in gap_index_cache:
            gap_index_cache[build] = load_gap_index(gaps_dir / f"{build}.gaps.bed")
        return gap_index_cache[build]

    out_header = [
        "anchor_id", "anchor_chrom", "anchor_pos", "anchor_end", "anchor_svtype", "anchor_svlen",
        "read_id", "build", "classification",
        "query_chrom", "query_pos", "coordinate_source",
        "near_gap", "gap_distance_bp", "gap_start", "gap_end", "gap_type",
    ]

    counts = defaultdict(lambda: defaultdict(int))  # counts[build][bucket]

    print("Pass 2/2: scanning severe rows...")
    n_severe = 0
    with open(trace_tsv, newline="") as f_in, open(out_path, "w", newline="") as f_out:
        reader = csv.DictReader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(out_header)

        for row in reader:
            if row["classification"] not in SEVERE:
                continue
            if row["anchor_id"] not in matched_somewhere:
                continue
            n_severe += 1
            build = row["build"]

            if row["classification"] == "OTHER_CHROM":
                q_chrom, q_pos = row["candidate_chrom"], int(row["candidate_pos"])
                coord_source = "candidate"
                result = nearest_gap(gap_index_for(build), q_chrom, q_pos)
                if result is None:
                    near, dist, gstart, gend, gtype = False, "NA", "NA", "NA", "NA"
                    counts[build]["other_chrom_no_gap_data"] += 1
                else:
                    dist, (gstart, gend, gtype) = result
                    near = dist <= args.window
                    counts[build]["other_chrom_near_gap" if near else "other_chrom_far_from_gap"] += 1
            else:  # NO_SV_CALL: no real coordinate in this build, approximate with the hs1 anchor position
                q_chrom, q_pos = row["anchor_chrom"], int(row["anchor_pos"])
                coord_source = "hs1_anchor_approx"
                result = nearest_gap(gap_index_for(build), q_chrom, q_pos)
                if result is None:
                    near, dist, gstart, gend, gtype = False, "NA", "NA", "NA", "NA"
                    counts[build]["no_sv_call_no_gap_data"] += 1
                else:
                    dist, (gstart, gend, gtype) = result
                    near = dist <= args.window
                    counts[build]["no_sv_call_near_gap" if near else "no_sv_call_far_from_gap"] += 1

            writer.writerow([
                row["anchor_id"], row["anchor_chrom"], row["anchor_pos"], row["anchor_end"],
                row["anchor_svtype"], row["anchor_svlen"], row["read_id"], build, row["classification"],
                q_chrom, q_pos, coord_source, near, dist, gstart, gend, gtype,
            ])

    print(f"\n{n_severe} severe (reference-sensitive) rows written to {out_path}\n")
    header = ["build", "OTHER_CHROM near_gap", "OTHER_CHROM far", "OTHER_CHROM no_gap_data",
              "NO_SV_CALL near_gap(approx)", "NO_SV_CALL far(approx)", "NO_SV_CALL no_gap_data"]
    widths = [8, 22, 18, 24, 28, 24, 24]
    print(" ".join(h.ljust(w) for h, w in zip(header, widths)))
    for build in sorted(counts):
        c = counts[build]
        vals = [build, c["other_chrom_near_gap"], c["other_chrom_far_from_gap"], c["other_chrom_no_gap_data"],
                c["no_sv_call_near_gap"], c["no_sv_call_far_from_gap"], c["no_sv_call_no_gap_data"]]
        print(" ".join(str(v).ljust(w) for v, w in zip(vals, widths)))


if __name__ == "__main__":
    main()
