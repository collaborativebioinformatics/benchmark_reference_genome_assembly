#!/usr/bin/env python3
"""Refine Part I classifications using the Part Ib alignment scan.

NO_SV_CALL in Part I conflates three mechanisms. With per-read alignment
coordinates they separate into:

    READ_ABSENT           read name not in the build's BAM at all
    READ_UNMAPPED         present, no mapped record
    MAPPED_SPLIT_NO_CALL  mapped with supplementary alignments, no SV called
    MAPPED_NO_CALL        mapped contiguously, no SV called

Every row also gets the gap context of the read's primary alignment in that
build (N-run overlap, distance to nearest gap, distance to contig end) and an
approximate read length, which Part I could not report.
"""

import argparse
import gzip
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

import pandas as pd

REFINED_NO_CALL = ["READ_ABSENT", "READ_UNMAPPED", "MAPPED_SPLIT_NO_CALL", "MAPPED_NO_CALL"]
CLASSES = ["MATCH", "SIZE_CHANGED", "TYPE_CHANGED", "BELOW_50", "OTHER_CHROM"] + REFINED_NO_CALL
NEAR_GAP_BP = 10_000


def norm_chrom(c):
    """chr1 / 1 / CHR1 -> 1. Build FASTAs do not agree on this."""
    if not isinstance(c, str):
        return c
    c = c.strip()
    return c[3:].upper() if c[:3].lower() == "chr" else c.upper()


class Intervals:
    """Disjoint sorted intervals per chromosome; overlap and nearest distance."""

    def __init__(self, bed_path=None):
        self.starts, self.ends = defaultdict(list), defaultdict(list)
        if bed_path and Path(bed_path).exists():
            with open(bed_path) as fh:
                for line in fh:
                    f = line.split()
                    if len(f) < 3:
                        continue
                    c = norm_chrom(f[0])
                    self.starts[c].append(int(f[1]))
                    self.ends[c].append(int(f[2]))
        for c in self.starts:
            order = sorted(range(len(self.starts[c])), key=lambda i: self.starts[c][i])
            self.starts[c] = [self.starts[c][i] for i in order]
            self.ends[c] = [self.ends[c][i] for i in order]

    def distance(self, chrom, s0, e0):
        """0 if the half-open [s0,e0) overlaps an interval; else bp to nearest.
        None if the chromosome has no intervals."""
        c = norm_chrom(chrom)
        st, en = self.starts.get(c), self.ends.get(c)
        if not st:
            return None
        i = bisect_right(st, e0 - 1)
        best = None
        if i > 0:
            if en[i - 1] > s0:
                return 0
            best = s0 - en[i - 1]
        if i < len(st):
            d = st[i] - e0
            best = d if best is None else min(best, d)
        return max(best, 0) if best is not None else None


def load_alignments(path, build):
    """Per (read_id) rollup of one build's alignment scan."""
    rows = {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r = f[idx["read_id"]]
            t = f[idx["aln_type"]]
            d = rows.get(r)
            if d is None:
                d = rows[r] = {
                    "n_primary": 0, "n_supplementary": 0, "n_secondary": 0, "n_unmapped": 0,
                    "primary_chrom": None, "primary_pos": None, "primary_end": None,
                    "primary_mapq": None, "primary_clip5": 0, "primary_clip3": 0,
                    "primary_span": 0, "_spans": {},
                }
            d["n_" + t] += 1
            if t in ("primary", "supplementary"):
                chrom = f[idx["chrom"]]
                pos, end = int(f[idx["pos"]]), int(f[idx["end"]])
                if t == "primary":
                    d["primary_chrom"] = chrom
                    d["primary_pos"], d["primary_end"] = pos, end
                    d["primary_mapq"] = int(f[idx["mapq"]])
                    d["primary_clip5"] = int(f[idx["clip5"]])
                    d["primary_clip3"] = int(f[idx["clip3"]])
                    d["primary_span"] = int(f[idx["ref_span"]])
                # keep spans per chromosome: record order in the BAM is arbitrary,
                # so the primary may arrive after its supplementaries
                mm = d["_spans"].get(chrom)
                d["_spans"][chrom] = [pos, end] if mm is None else [min(mm[0], pos), max(mm[1], end)]

    for d in rows.values():
        mm = d["_spans"].get(d["primary_chrom"])
        d["locus_span"] = (mm[1] - mm[0] + 1) if mm else None
        if d["n_primary"] == 0:
            d["read_len_approx"] = None
        else:
            # ref_span + clipping on the primary alignment approximates read length
            d["read_len_approx"] = d["primary_span"] + d["primary_clip5"] + d["primary_clip3"]
        del d["_spans"]

    out = pd.DataFrame.from_dict(rows, orient="index")
    if out.empty:
        return out
    out.index.name = "read_id"
    out = out.reset_index()
    out["build"] = build
    out["read_mapped"] = out["n_primary"] > 0
    return out


def refine(row):
    c = row["classification"]
    if c != "NO_SV_CALL":
        return c
    if not row["read_in_bam"]:
        return "READ_ABSENT"
    if not row["read_mapped"]:
        return "READ_UNMAPPED"
    if row["n_supplementary"] > 0:
        return "MAPPED_SPLIT_NO_CALL"
    return "MAPPED_NO_CALL"


def gap_context(d):
    if d is None or pd.isna(d):
        return "NA"
    if d == 0:
        return "IN_GAP"
    return "NEAR_GAP" if d <= NEAR_GAP_BP else "CLEAR"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trace-dir", required=True, help="Part I output directory")
    ap.add_argument("--bam-dir", required=True, help="Part Ib output directory")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--anchor-annotation", action="append", default=[], metavar="NAME=BED",
                    help="hs1 BED to intersect against anchor intervals; repeatable "
                         "(e.g. segdup=chm13v2.0_SD.bed)")
    args = ap.parse_args()

    trace, bamd = Path(args.trace_dir), Path(args.bam_dir)
    out = Path(args.out_dir) if args.out_dir else trace / "refined"
    out.mkdir(parents=True, exist_ok=True)

    best = pd.read_csv(trace / "read_trace_best.tsv", sep="\t", low_memory=False)
    builds = sorted(best["build"].unique())

    # ---- contig naming check ------------------------------------------------
    naming = {}
    for b in builds:
        s = pd.read_csv(bamd / f"sizes.{b}.tsv", sep="\t", header=None, names=["chrom", "len"]) \
            if (bamd / f"sizes.{b}.tsv").exists() else None
        if s is not None and len(s):
            naming[b] = "chr-prefixed" if str(s["chrom"].iloc[0]).lower().startswith("chr") else "bare"
    if len(set(naming.values())) > 1:
        print(f"WARNING: contig naming differs across builds: {naming}\n"
              "  Part I compares chromosome names as raw strings, so its same_chrom /\n"
              "  OTHER_CHROM calls are not trustworthy for the mismatched builds.")

    # ---- alignments ---------------------------------------------------------
    aln = []
    for b in builds:
        p = bamd / f"aln.{b}.tsv.gz"
        if not p.exists():
            print(f"WARNING: missing {p}, build {b} keeps its Part I classifications")
            continue
        aln.append(load_alignments(p, b))
    aln = pd.concat(aln, ignore_index=True) if aln else pd.DataFrame(columns=["read_id", "build"])

    df = best.merge(aln, on=["read_id", "build"], how="left")
    df["read_in_bam"] = df["n_primary"].notna()
    for c in ["n_primary", "n_supplementary", "n_secondary"]:
        df[c] = df[c].fillna(0).astype(int)
    df["read_mapped"] = df["read_mapped"].fillna(False)
    df["classification_refined"] = df.apply(refine, axis=1)

    # OTHER_CHROM is a statement about the SV record, not about the read.
    df["primary_chrom_matches_anchor"] = [
        None if not isinstance(pc, str) else int(norm_chrom(pc) == norm_chrom(ac))
        for pc, ac in zip(df["primary_chrom"], df["anchor_chrom"])
    ]

    # ---- gap context of the primary alignment in each build -----------------
    gaps = {b: Intervals(bamd / f"gaps.{b}.bed") for b in builds}
    sizes = {}
    for b in builds:
        p = bamd / f"sizes.{b}.tsv"
        if p.exists():
            s = pd.read_csv(p, sep="\t", header=None, names=["chrom", "len"])
            sizes[b] = {norm_chrom(c): int(l) for c, l in zip(s["chrom"], s["len"])}

    dist, to_end = [], []
    for b, c, pos, end in zip(df["build"], df["primary_chrom"], df["primary_pos"], df["primary_end"]):
        if not isinstance(c, str) or c == "*" or pd.isna(pos):
            dist.append(None); to_end.append(None); continue
        dist.append(gaps[b].distance(c, int(pos) - 1, int(end)))
        L = sizes.get(b, {}).get(norm_chrom(c))
        to_end.append(min(int(pos) - 1, L - int(end)) if L else None)
    df["dist_to_gap"] = dist
    df["dist_to_contig_end"] = to_end
    df["gap_context"] = [gap_context(d) for d in dist]

    df.to_csv(out / "read_trace_best_refined.tsv", sep="\t", index=False)

    # ---- per-build summary --------------------------------------------------
    tab = (df.groupby(["build", "classification_refined"]).size()
             .unstack(fill_value=0).reindex(columns=CLASSES, fill_value=0))
    tab["n"] = tab.sum(axis=1)
    pctcols = {}
    for c in CLASSES:
        pctcols[f"{c.lower()}_pct"] = 100.0 * tab[c] / tab["n"]
    summary = pd.concat([tab, pd.DataFrame(pctcols)], axis=1).reset_index()
    summary.to_csv(out / "build_summary_refined.tsv", sep="\t", index=False, float_format="%.2f")

    # ---- gap attribution ----------------------------------------------------
    nocall = df[df["classification_refined"].isin(REFINED_NO_CALL)]
    gap_attr = (nocall.groupby(["build", "classification_refined", "gap_context"]).size()
                .rename("n").reset_index())
    gap_attr.to_csv(out / "no_call_gap_attribution.tsv", sep="\t", index=False)

    # ---- per anchor x build failure mode ------------------------------------
    def mode(s):
        vc = s.value_counts()
        return max(CLASSES, key=lambda c: (vc.get(c, 0), -CLASSES.index(c)))

    per_anchor = (df.groupby(["anchor_id", "build"])
                    .agg(anchor_chrom=("anchor_chrom", "first"),
                         anchor_pos=("anchor_pos", "first"),
                         anchor_end=("anchor_end", "first"),
                         anchor_svtype=("anchor_svtype", "first"),
                         anchor_svlen=("anchor_svlen", "first"),
                         n_reads=("read_id", "size"),
                         dominant_class=("classification_refined", mode),
                         frac_match=("classification_refined", lambda s: (s == "MATCH").mean()),
                         frac_in_gap=("gap_context", lambda s: (s == "IN_GAP").mean()),
                         frac_near_gap=("gap_context", lambda s: (s == "NEAR_GAP").mean()),
                         median_gap_dist=("dist_to_gap", "median"),
                         median_read_len=("read_len_approx", "median"))
                    .reset_index())
    per_anchor.to_csv(out / "anchor_build_failure_modes.tsv", sep="\t",
                      index=False, float_format="%.4f")

    # ---- problem regions on the hs1 side ------------------------------------
    piv = per_anchor.pivot(index="anchor_id", columns="build", values="dominant_class")
    piv = piv.reindex(sorted(piv.columns), axis=1).add_prefix("cls_")
    gapf = per_anchor.pivot(index="anchor_id", columns="build", values="frac_in_gap")
    gapf = gapf.reindex(sorted(gapf.columns), axis=1).add_prefix("ingap_")
    meta = per_anchor.drop_duplicates("anchor_id").set_index("anchor_id")[
        ["anchor_chrom", "anchor_pos", "anchor_end", "anchor_svtype", "anchor_svlen"]]
    prob = meta.join(piv).join(gapf)
    cls_cols = [c for c in prob.columns if c.startswith("cls_")]
    prob["n_builds_match"] = (prob[cls_cols] == "MATCH").sum(axis=1)
    prob["n_builds_lost"] = prob[cls_cols].isin(REFINED_NO_CALL).sum(axis=1)
    prob["n_builds_gap_attributed"] = (prob[[c for c in prob.columns
                                             if c.startswith("ingap_")]] > 0.5).sum(axis=1)
    prob["reference_sensitive"] = (prob["n_builds_match"] > 0) & \
                                  ((prob[cls_cols] != "MATCH").sum(axis=1) > 0)

    for spec in args.anchor_annotation:
        name, _, path = spec.partition("=")
        iv = Intervals(path)
        prob[f"anno_{name}"] = [
            1 if (iv.distance(c, int(s) - 1, int(e)) == 0) else 0
            for c, s, e in zip(prob["anchor_chrom"], prob["anchor_pos"], prob["anchor_end"])
        ]

    prob = prob.sort_values(["n_builds_lost", "n_builds_match"], ascending=[False, False])
    prob.to_csv(out / "problem_regions.tsv", sep="\t", index=True,
                index_label="anchor_id", float_format="%.4f")

    # ---- terminal summary ---------------------------------------------------
    show = ["build", "n"] + [c for c in CLASSES if tab[c].sum() > 0]
    print("\nRefined classification per supporting read:\n")
    print(summary[show].to_string(index=False))
    print("\nNO_SV_CALL breakdown by gap context:\n")
    print(gap_attr.pivot_table(index=["build", "classification_refined"],
                               columns="gap_context", values="n",
                               fill_value=0).to_string())
    print(f"\nReference-sensitive anchors: {int(prob['reference_sensitive'].sum())} / {len(prob)}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
