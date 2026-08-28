# Results — HS1 read-trace (minimap2 arm)

Traces every hs1-anchor SV's Sniffles2-supporting reads (by read name, no liftover)
across the other six minimap2/Sniffles2 callsets (`scripts/trace_sv_reads_part1_samechr.sh`,
DNAnexus applet under `scripts/read-trace-minimap2/`), summarizes the outcomes
(`scripts/summarize_sv_trace.py`), and cross-references misses against per-build
assembly gap BEDs (`scripts/check_gap_followup.py`).

22,704 hs1 PASS anchor SVs (size ≥50 bp, ≥1 supporting read; mean 8.56 supporting
reads/SV, min 3, max 108) traced against hg15, hg16, hg17, hg18, hg19, hg38.

Full per-read/per-SV TSVs (`read_trace_best.tsv`, `gap_followup_candidates.tsv`,
`sv_reference_patterns.tsv`, `gap_candidate_summary.tsv`, classification/match-fraction
matrices) are in DNAnexus under `Group4_2026:/results/read_trace/minimap2/20260828-104421/`.
Only the two small per-build summary tables are kept here.

| File | Contents |
|---|---|
| `build_summary_all_reads.tsv` | Classification counts/percentages per build, one row per supporting read |
| `build_summary_per_sv_consensus.tsv` | Same, but one row per (SV, build) using the majority-vote classification across that SV's supporting reads |

## Per-read outcome, by build

| build | n reads traced | MATCH | OTHER_CHROM | TYPE_CHANGED | SIZE_CHANGED | BELOW_50 | NO_SV_CALL |
|---|---:|---:|---:|---:|---:|---:|---:|
| hg15 | 194,301 | 31.9% | 2.25% | 8.4% | 9.8% | 0.2% | 47.4% |
| hg16 | 194,301 | 32.2% | 3.09% | 8.5% | 9.8% | 0.3% | 46.2% |
| hg17 | 194,301 | 32.5% | 2.53% | 8.6% | 10.1% | 0.2% | 46.1% |
| hg18 | 194,301 | 32.5% | 2.10% | 8.6% | 10.1% | 0.3% | 46.5% |
| hg19 | 194,301 | 32.6% | 1.31% | 8.7% | 10.1% | 0.2% | 47.2% |
| hg38 | 194,301 | 31.5% | 0.87% | 8.5% | 9.3% | 0.2% | 49.7% |

`NO_SV_CALL` (no variant record in that build's VCF lists the read as support at all) is
the single largest outcome in every build — larger than `MATCH`.

## Per-SV majority-vote outcome

1,299 / 22,704 SVs (5.7%) are **reference-sensitive**: majority-MATCH in at least one
build and majority-disrupted in at least one other. 1,089 of those are **HIGH-priority**
(a severe outcome — `OTHER_CHROM`/`NO_SV_CALL` — in at least one build).

Read support count is not a clean predictor of outcome — it's non-monotonic:

| supporting reads (in hs1) | read-frac MATCH | read-frac NO_SV_CALL |
|---|---:|---:|
| 3–4 | 24.3% | 47.2% |
| 5–8 | 35.1% | 40.4% |
| 9–15 | 33.9% | 47.1% |
| 16+ | 27.1% | **58.0%** |

The worst outcomes are at the *highest* support, not the lowest — inconsistent with a
simple "too few reads to clear Sniffles' auto min-support threshold" explanation (all
three calling scripts use identical Sniffles2 settings, `--minsvlen 50 --output-rnames`,
no `--minsupport` override, so the threshold is the same everywhere). More likely: high
supporting-read counts mark repetitive/complex loci (segmental duplications, VNTRs)
where reads multi-map, and those same loci are where reference assemblies diverge most.

## Gap-proximity check

Of 38,814 severe (`OTHER_CHROM`/`NO_SV_CALL`) read-level misses, very few `OTHER_CHROM`
candidates (the only classification with a real query-build coordinate) land within 1 kb
of a known assembly gap:

| build | OTHER_CHROM near gap | OTHER_CHROM far from gap |
|---|---:|---:|
| hg15 | 29 | 857 |
| hg16 | 16 | 1,032 |
| hg17 | 0 | 367 |
| hg18 | 0 | 339 |
| hg19 | 2 | 145 |
| hg38 | 1 | 54 |

Assembly gaps explain almost none of the severe misses.

## Correlation with assembly quality (Table 1, `data/README.md`)

Pearson r across the 6 non-anchor builds:

| comparison | r |
|---|---:|
| Contig N50 vs `OTHER_CHROM`% | −0.76 |
| Release year vs `OTHER_CHROM`% | −0.86 |
| Contig N50 vs `NO_SV_CALL`% | +0.71 |
| Release year vs `NO_SV_CALL`% | **+0.91** |
| Contig N50 vs `MATCH`% | −0.38 |
| Contig N50 vs gap-proximity fraction | −0.29 |

Two opposite trends:

- **Cross-chromosome misplacement tracks assembly fragmentation as expected** — older,
  more fragmented builds (hg15–hg17: 1,000–2,240 contigs) misplace reads onto the wrong
  chromosome far more than hg38 (0.87% vs up to 3.09%).
- **The dominant failure mode goes the other way.** `NO_SV_CALL` (~46–50% of all reads)
  *increases* with assembly quality/recency (r=+0.91 with release year). Plausible
  mechanism: better repeat-resolution in newer assemblies normalizes/collapses
  repetitive loci that older, less-refined assemblies left more "open," letting a read
  cluster into a callable (if geometrically wrong) signal there instead of splitting
  across secondary alignments below Sniffles' support threshold in the better assembly.
- Gap-proximity doesn't explain either trend (r=−0.29) — it's assembly-wide repeat
  handling, not literal annotated gaps, driving this.

**Caveat:** n=6 data points — treat these correlations as suggestive, not statistically
conclusive.

## Reproducibility

```bash
# Part I: read-name tracing (DNAnexus applet, scripts/read-trace-minimap2/)
scripts/trace_sv_reads_part1_samechr.sh

# Part II: summarize
python3 scripts/summarize_sv_trace.py --input-dir <Part I output dir>

# Part III: gap follow-up (needs <build>.gaps.bed per query build)
python3 scripts/check_gap_followup.py --trace-dir <Part I output dir> --gaps-dir <dir with gap BEDs>
```
