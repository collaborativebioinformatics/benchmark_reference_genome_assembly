#!/usr/bin/env python3
"""Relate SV recovery to reference assembly properties.

Two levels, deliberately separated:

  build level  n = 7, and every Table 1 metric is close to monotone in release
               year, so nothing here identifies a mechanism. Reported with the
               collinearity in view.

  locus level  n = anchors x builds. This is the level that can support a model
               and the level that generalises to a reference not in the panel,
               because its predictors are computed from the reference rather
               than looked up per build.

Table 1 statistics come from NCBI assembly reports. Contig counts and gap
content are recomputed here from the FASTAs that were actually indexed; where
the two disagree the recomputed value is the relevant one.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Table 1 (NCBI Datasets assembly reports). hg38 = GRCh38 initial; the p14 row
# has identical N50s so the choice does not affect any rank statistic.
TABLE1 = pd.DataFrame([
    ("hs1",  2022, 3117275501, 150617247, 150617247,   24, "GCF_009914755.1"),
    ("hg38", 2013, 3099734149,  57879411,  67794873,  999, "GCF_000001405.26"),
    ("hg19", 2009, 3101788170,  38508932,  46395641,  350, "GCF_000001405.13"),
    ("hg18", 2006, 3093104542,  38509590,  38509590, 1006, "GCF_000001405.12"),
    ("hg17", 2004, 3091360260,  37760040,  38509590, 1215, "GCF_000001405.11"),
    ("hg16", 2003, 3091959510,  28857747,  29104798, 1756, "GCF_000001405.10"),
    ("hg15", 2003, 3095784245,  23437594,  25443670, 2240, "GCF_000001405.8"),
], columns=["build", "year", "total_len_ncbi", "contig_n50", "scaffold_n50",
            "n_contigs_ncbi", "accession"])


def measured_stats(bam_dir, builds):
    """Contig count and N content of the FASTA each build was indexed from."""
    rows = []
    for b in builds:
        sp, gp = Path(bam_dir) / f"sizes.{b}.tsv", Path(bam_dir) / f"gaps.{b}.bed"
        if not sp.exists():
            continue
        s = pd.read_csv(sp, sep="\t", header=None, names=["chrom", "len"])
        gap_bp = n_gaps = 0
        if gp.exists():
            g = pd.read_csv(gp, sep="\t", header=None, names=["chrom", "start", "end"])
            gap_bp, n_gaps = int((g["end"] - g["start"]).sum()), len(g)
        rows.append({"build": b, "n_contigs_indexed": len(s),
                     "total_len_indexed": int(s["len"].sum()),
                     "gap_bp": gap_bp, "n_gaps": n_gaps,
                     "gap_fraction": gap_bp / max(int(s["len"].sum()), 1)})
    return pd.DataFrame(rows)


def partial_spearman(x, y, z):
    """Spearman of x,y after removing a linear fit on ranked z."""
    rx, ry, rz = (pd.Series(v).rank().to_numpy(float) for v in (x, y, z))
    A = np.column_stack([np.ones_like(rz), rz])
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return spearmanr(ex, ey).statistic


def auc(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float)
    pos, neg = y == 1, y == 0
    if pos.sum() == 0 or neg.sum() == 0:
        return float("nan")
    r = pd.Series(p).rank().to_numpy()
    return (r[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refined-dir", required=True,
                    help="output directory of reclassify_with_alignments.py")
    ap.add_argument("--bam-dir", required=True, help="Part Ib output directory")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    ref = Path(args.refined_dir)
    out = Path(args.out_dir) if args.out_dir else ref
    out.mkdir(parents=True, exist_ok=True)

    per_anchor = pd.read_csv(ref / "anchor_build_failure_modes.tsv", sep="\t")
    builds = sorted(per_anchor["build"].unique())

    stats = TABLE1.merge(measured_stats(args.bam_dir, builds), on="build", how="right")

    bad = stats[stats["n_contigs_ncbi"].notna() &
                (abs(stats["n_contigs_ncbi"] - stats["n_contigs_indexed"]) >
                 0.25 * stats["n_contigs_ncbi"])]
    if len(bad):
        print("WARNING: Table 1 contig counts disagree with the indexed FASTAs by >25% for: "
              + ", ".join(f"{r.build} ({int(r.n_contigs_ncbi)} vs {int(r.n_contigs_indexed)})"
                          for r in bad.itertuples()))
        print("  Reads were mapped to the indexed FASTA, so n_contigs_indexed is the\n"
              "  covariate that corresponds to the experiment.\n")

    # ---- build level --------------------------------------------------------
    recov = (per_anchor.assign(is_match=lambda d: d["dominant_class"] == "MATCH")
             .groupby("build")
             .agg(n_anchors=("anchor_id", "size"),
                  match_rate=("is_match", "mean"),
                  mean_frac_in_gap=("frac_in_gap", "mean"))
             .reset_index())
    stats = stats.merge(recov, on="build", how="inner")
    stats.to_csv(out / "build_assembly_stats.tsv", sep="\t", index=False, float_format="%.6g")

    metrics = ["year", "contig_n50", "scaffold_n50", "n_contigs_ncbi",
               "n_contigs_indexed", "total_len_indexed", "gap_fraction", "n_gaps"]
    metrics = [m for m in metrics if m in stats and stats[m].notna().sum() >= 4
               and stats[m].nunique() > 1]
    if len(stats) < 4 or not metrics:
        print(f"\nOnly {len(stats)} build(s) with stats; skipping build-level correlation.")
        metrics = []

    rows = []
    for m in metrics:
        # dict.fromkeys keeps the selection unique when m is "year" itself
        sub = stats[list(dict.fromkeys(["match_rate", m, "year"]))].dropna()
        rho = spearmanr(sub[m], sub["match_rate"]).statistic
        rows.append({"metric": m, "n": len(sub), "spearman_rho": rho,
                     "rho_vs_year": spearmanr(sub[m], sub["year"]).statistic,
                     "partial_rho_given_year": (float("nan") if m == "year"
                                                else partial_spearman(sub[m], sub["match_rate"],
                                                                      sub["year"]))})
    if rows:
        corr = pd.DataFrame(rows).sort_values("spearman_rho", key=abs, ascending=False)
        corr.to_csv(out / "build_metric_correlations.tsv", sep="\t",
                    index=False, float_format="%.4f")
        print(f"\nBuild-level association with per-anchor match rate "
              f"(n = {len(stats)} builds):\n")
        print(corr.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        print("\n  rho_vs_year is the diagnostic, not a nuisance column. Any metric with\n"
              "  |rho_vs_year| near 1 is a restatement of release date at this n, and the\n"
              "  partial correlation is fitted on ~4 residual df. hg15-hg19 are successive\n"
              "  revisions of one assembly lineage, so the effective number of independent\n"
              "  references is closer to three (NCBI33-36, GRCh37/38, T2T) than to seven.")

    # ---- locus level --------------------------------------------------------
    d = per_anchor.copy()
    d["y"] = (d["dominant_class"] == "MATCH").astype(int)
    d["log_gap_dist"] = np.log10(d["median_gap_dist"].fillna(1e7).clip(lower=1) + 1)
    d["in_gap"] = (d["frac_in_gap"] > 0.5).astype(int)
    d["log_svlen"] = np.log10(d["anchor_svlen"].clip(lower=1))
    d["log_readlen"] = np.log10(d["median_read_len"].fillna(d["median_read_len"].median())
                                .clip(lower=1))
    feats = ["in_gap", "log_gap_dist", "log_svlen", "log_readlen"]

    try:
        import statsmodels.api as sm
    except ImportError:
        print("\nstatsmodels not installed; skipping the locus-level model.")
        return

    dropped = [f for f in feats if d[f].nunique() < 2]
    if dropped:
        print(f"\nDropping zero-variance predictors: {', '.join(dropped)}")
        feats = [f for f in feats if f not in dropped]
    if not feats or d["y"].nunique() < 2:
        print("\nNot enough variation for the locus-level model.")
        return

    X = pd.concat([d[feats], pd.get_dummies(d["build"], prefix="build", drop_first=True)], axis=1)
    X = sm.add_constant(X.astype(float), has_constant="add")
    fit = sm.GLM(d["y"], X, family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": d["anchor_id"]})
    res = pd.DataFrame({"coef": fit.params, "se": fit.bse, "z": fit.tvalues, "p": fit.pvalues,
                        "odds_ratio": np.exp(fit.params)})
    res.to_csv(out / "locus_model_coefficients.tsv", sep="\t", float_format="%.4g")
    print("\nLocus-level logistic model, P(MATCH), build fixed effects, "
          "SEs clustered by anchor:\n")
    print(res.reindex(["const"] + feats).dropna(how="all")
             .to_string(float_format=lambda x: f"{x:.3f}"))

    # Leave-one-build-out. hs1 is the only genuinely independent held-out
    # reference; the GRCh-lineage folds leak into each other.
    lobo = []
    for b in builds:
        tr, te = d["build"] != b, d["build"] == b
        if d.loc[tr, "y"].nunique() < 2 or d.loc[te, "y"].nunique() < 2:
            lobo.append({"held_out_build": b, "n": int(te.sum()), "auc": float("nan")})
            continue
        Xt = sm.add_constant(d.loc[tr, feats].astype(float), has_constant="add")
        Xe = sm.add_constant(d.loc[te, feats].astype(float), has_constant="add")
        m = sm.GLM(d.loc[tr, "y"], Xt, family=sm.families.Binomial()).fit()
        lobo.append({"held_out_build": b, "n": int(te.sum()),
                     "auc": auc(d.loc[te, "y"], m.predict(Xe))})
    lobo = pd.DataFrame(lobo)
    lobo.to_csv(out / "locus_model_lobo.tsv", sep="\t", index=False, float_format="%.4f")
    print("\nLeave-one-build-out AUC (reference-only predictors, no build effects):\n")
    print(lobo.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nOutputs: {out}")


if __name__ == "__main__":
    main()
