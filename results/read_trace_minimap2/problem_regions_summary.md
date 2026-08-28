# Summary: `problem_regions.tsv`

One row per hs1-anchor structural variant (n = 22,704), classified against each
of the 6 non-anchor reference builds (hg15, hg16, hg17, hg18, hg19, hg38) using
the refined per-read alignment/gap classification from
`reclassify_with_alignments.py`.

## Three broad outcome groups

| Group | n | % | Meaning |
|---|---:|---:|---|
| Universally recovered | 6,740 | 29.7% | MATCH in all 6 non-anchor builds — a reference-agnostic "core" SV set |
| **Reference-sensitive** | **1,302** | **5.7%** | MATCH in ≥1 build, lost/failed in ≥1 other — the group that actually answers "does reference choice matter" |
| Never recovered elsewhere | 14,662 | 64.6% | Zero matches in any of the 6 other builds — either hs1-private true calls only resolvable at T2T resolution, or false positives |

Of the "never recovered elsewhere" group, 8,503 are fully lost — a hard
failure (read absent, unmapped, or mapped with no call) in *all six* builds.
The remaining ~6,159 show some other non-match outcome in at least one build
(wrong SV type, size drift, wrong chromosome) without ever being a clean loss
everywhere.

## The reference-sensitive 1,302, broken down further

- **SV type**: INS (658) and DEL (629) dominate; DUP/INV are negligible (9/6).
  INS is modestly over-represented here relative to its 45.7% share of all
  anchors — insertions appear somewhat more reference-dependent to detect
  than deletions.
- **Size**: skews larger on average than non-sensitive anchors (mean SVLEN
  ≈ 20.9 kb vs ≈ 8.9 kb). That mean is pulled by a few very large outliers,
  so treat it as "skews larger" rather than a precise point estimate without
  pulling the full size distribution.

## Headline finding: gap attribution is surprisingly small

Of the 11,613 anchors that lose in at least one build, only **1,012 (≈8.7%)**
have that loss actually landing inside a real assembly gap in that build.

**Assembly gaps are *not* the dominant explanation for reference-driven SV
loss in this dataset.** The bulk of losses trace to other failure modes
(see `05_build_summary_refined.tsv`): unmapped reads, split/supplementary
alignments, or a clean mapping with simply no SV call. This is a fairly
direct, quotable counter to the a priori assumption that older/gappier
assemblies mostly fail *because* of their gaps.

---

*Source: `06_problem_regions.tsv`, derived from
`anchor_build_failure_modes.tsv` (`reclassify_with_alignments.py`), itself
built from Part I (`trace_sv_reads_part1_samechr.sh`) and Part Ib
(`trace_sv_reads_part1b_bamscan.sh`, run on DNAnexus,
job-JB8p9pQ0Z7Q3KpvjqKV57pb4) outputs.*
