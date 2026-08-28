#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

CLASSES = ["MATCH", "OTHER_CHROM", "TYPE_CHANGED", "SIZE_CHANGED", "BELOW_50", "NO_SV_CALL"]
SEVERE = {"OTHER_CHROM", "NO_SV_CALL"}


def pct(x, n):
    return 100.0 * x / n if n else float("nan")


def build_summary(df, label):
    rows = []
    for build, sub in df.groupby("build", sort=False):
        counts = sub["classification"].value_counts()
        n = len(sub)
        row = {"dataset": label, "build": build, "n": n}
        for c in CLASSES:
            k = int(counts.get(c, 0))
            row[c.lower()] = k
            row[f"{c.lower()}_pct"] = pct(k, n)
        rows.append(row)
    return pd.DataFrame(rows)


def dominant_class(row):
    counts = {
        "MATCH": row["n_match"],
        "OTHER_CHROM": row["n_other_chrom"],
        "TYPE_CHANGED": row["n_type_changed"],
        "SIZE_CHANGED": row["n_size_changed"],
        "BELOW_50": row["n_below_50"],
        "NO_SV_CALL": row["n_no_sv_call"],
    }
    # Stable tie break: prefer MATCH, then increasingly severe alternatives in CLASSES order.
    return max(CLASSES, key=lambda c: (counts[c], -CLASSES.index(c)))


def main():
    ap = argparse.ArgumentParser(description="Summarize HS1 supporting-read tracing across reference builds.")
    ap.add_argument("--input-dir", default=".", help="Directory containing tracing TSV outputs")
    ap.add_argument("--output-dir", default=None, help="Output directory; default: <input-dir>/summary")
    args = ap.parse_args()

    inp = Path(args.input_dir)
    out = Path(args.output_dir) if args.output_dir else inp / "summary"
    out.mkdir(parents=True, exist_ok=True)

    best = pd.read_csv(inp / "read_trace_best.tsv", sep="\t")
    sv = pd.read_csv(inp / "sv_trace_summary.tsv", sep="\t")

    # 1) Global build summary across all supporting reads.
    # (No per-SV "longest read" view: read length isn't tracked anywhere upstream --
    # Part I traces reads from VCF RNAMES only, and it's the same underlying reads
    # across all seven builds' BAMs, so getting it would mean a separate full BAM
    # scan. Section 2 below covers the "one call per SV per build" role instead,
    # via majority vote across all supporting reads rather than picking one read.)
    all_sum = build_summary(best, "all_supporting_reads")
    all_sum.to_csv(out / "build_summary_all_reads.tsv", sep="\t", index=False, float_format="%.2f")

    # 2) Equal-weight-per-SV consensus from all supporting reads.
    sv = sv.copy()
    sv["dominant_class"] = sv.apply(dominant_class, axis=1)
    sv["dominant_fraction"] = sv.apply(
        lambda r: max(r["n_match"], r["n_other_chrom"], r["n_type_changed"], r["n_size_changed"], r["n_below_50"], r["n_no_sv_call"]) / r["n_hs1_supporting_reads"],
        axis=1,
    )
    consensus_sum = build_summary(
        sv.rename(columns={"dominant_class": "classification"}),
        "per_sv_majority_of_all_supporting_reads",
    )
    consensus_sum.to_csv(out / "build_summary_per_sv_consensus.tsv", sep="\t", index=False, float_format="%.2f")

    # 3) One row per HS1 SV: majority-vote classification pattern across builds.
    meta_cols = ["anchor_id", "anchor_chrom", "anchor_pos", "anchor_end", "anchor_svtype", "anchor_svlen"]
    meta = sv[meta_cols].drop_duplicates("anchor_id").set_index("anchor_id")
    class_mat = sv.pivot(index="anchor_id", columns="build", values="dominant_class")
    class_mat = class_mat.reindex(sorted(class_mat.columns), axis=1)
    class_mat.to_csv(out / "consensus_classification_matrix.tsv", sep="\t")

    frac_mat = sv.pivot(index="anchor_id", columns="build", values="match_fraction")
    frac_mat = frac_mat.reindex(sorted(frac_mat.columns), axis=1)
    frac_mat.to_csv(out / "all_reads_match_fraction_matrix.tsv", sep="\t", float_format="%.4f")

    pattern = meta.join(class_mat.add_prefix("cons_"), how="left").join(frac_mat.add_prefix("matchfrac_"), how="left")
    class_cols = [c for c in pattern.columns if c.startswith("cons_")]
    pattern["n_builds_match_cons"] = (pattern[class_cols] == "MATCH").sum(axis=1)
    pattern["n_builds_disrupted_cons"] = (pattern[class_cols] != "MATCH").sum(axis=1)
    pattern["n_builds_severe_cons"] = pattern[class_cols].isin(SEVERE).sum(axis=1)
    pattern["cons_pattern"] = pattern[class_cols].astype(str).agg(" | ".join, axis=1)
    pattern["reference_sensitive"] = (pattern["n_builds_match_cons"] > 0) & (pattern["n_builds_disrupted_cons"] > 0)
    pattern["gap_followup_priority"] = pattern.apply(
        lambda r: "HIGH" if r["reference_sensitive"] and r["n_builds_severe_cons"] > 0
        else ("MEDIUM" if r["reference_sensitive"] else "LOW"), axis=1
    )
    priority_rank = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    pattern["_priority_rank"] = pattern["gap_followup_priority"].map(priority_rank)
    pattern = pattern.sort_values(["_priority_rank", "n_builds_severe_cons", "n_builds_disrupted_cons"], ascending=[True, False, False]).drop(columns="_priority_rank")
    pattern.to_csv(out / "sv_reference_patterns.tsv", sep="\t", index=True, index_label="anchor_id", float_format="%.4f")

    # 4) Candidate rows for systematic gap/reference-structure follow-up.
    # Require a severe outcome for a read in one build and a MATCH for some read
    # (possibly a different read of the same SV) in at least one other build.
    # Built from every supporting read (not just one representative read per SV),
    # so this is the full candidate set -- feed straight into check_gap_followup.py.
    matched_somewhere = set(best.loc[best["classification"] == "MATCH", "anchor_id"])
    gap_candidates = best[
        best["anchor_id"].isin(matched_somewhere) & best["classification"].isin(SEVERE)
    ].copy()
    keep = [
        "anchor_id", "anchor_chrom", "anchor_pos", "anchor_end", "anchor_svtype", "anchor_svlen",
        "read_id", "build", "classification", "candidate_id", "candidate_chrom",
        "candidate_pos", "candidate_end", "candidate_svtype", "candidate_svlen", "candidate_filter",
        "n_candidates_for_read", "ambiguous_multi_candidate"
    ]
    gap_candidates[keep].sort_values(["anchor_id", "build"]).to_csv(
        out / "gap_followup_candidates.tsv", sep="\t", index=False
    )

    # 5) Small terminal summary.
    print("\nMajority-vote classification per SV, across all supporting reads:\n")
    display_cols = ["build", "n"] + sum(([c.lower(), f"{c.lower()}_pct"] for c in CLASSES), [])
    print(consensus_sum[display_cols].to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    print(f"\nReference-sensitive SVs: {int(pattern['reference_sensitive'].sum())} / {len(pattern)}")
    print(f"High-priority gap/reference-structure candidates: {(pattern['gap_followup_priority'] == 'HIGH').sum()}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
