# BCM Hackathon 2026 - Benchmarking the Impact of Reference Genome Quality on SV Discovery

## Intro statement

## Problem

Structural variant detection is always relative to a reference assembly, and while humans now have a complete telomere-to-telomere reference, most species are stuck with fragmented draft assemblies. Researchers calling SVs against such drafts have no principled way to know what their callset is missing. Therefore, we would like to systematically investigate how assembly quality influences SV calling.

## Why should we solve it

## Workflow Diagram [in-works]

![Workflow Illustration](group4_flowchart_diagram.png)

## Methodology

Two arms — **winnowmap→sniffles** and **minimap2→sniffles** — are run on identical reads against each reference build (e.g. `hs1`, `hg38`, `hg19`, `hg18`, `hg17`, `hg16`, `hg15`), so SV calls can be compared across references/aligners.

### 0. Read subset (fixed across both arms)

20 of 39 SMRT cells (every second cell after sorting movie names) = 86.4 of 166.2 GB, 52.0% of bases, ≈15x coverage. Selected once, before alignment:

```bash
dx download Group4_2026:/config/subset_cells_15x.txt -o ~/cells.txt
```

### 1. Winnowmap arm

```bash
dx run app-cloud_workstation --ssh --instance-type mem2_ssd1_v2_x32 -imax_session_length=24h
# inside the workstation:
unset DX_WORKSPACE_ID; dx cd $DX_PROJECT_CONTEXT_ID:
dx download Group4_2026:/scripts/run_winnowmap.sh
bash run_winnowmap.sh hs1 2>&1 | tee run_hs1.log
```

`run_winnowmap.sh` picks up `~/cells.txt` automatically if present, and refuses to run if `COVERAGE` is also set — the two subsampling methods can't be combined by accident.

**Tools:** Winnowmap 2.03 (`-ax map-pb`, `-W` repetitive k-mers from meryl 1.4.2 at k=15, distinct 0.9998, `-Y`, 30 threads) → samtools 1.24 sort/index → Sniffles2 2.3.2 (`--minsvlen 50 --output-rnames`). Alignment is whole-genome; restricted to chr6/chr22 afterward.

Output: `Group4_2026:/results/winnowmap/<build>/`

If Sniffles needs to be re-run against an already-aligned winnowmap BAM (e.g. after a tool-version fix) without redoing alignment, use `rerun_sniffles_winnowmap.sh` instead — it finds the existing BAM/reference on the workstation, rebuilds the tool env (pinned to Python 3.10), and re-calls:

```bash
dx download Group4_2026:/scripts/rerun_sniffles_winnowmap.sh
bash rerun_sniffles_winnowmap.sh hg17 2>&1 | tee rerun_sniffles_hg17.log
```

### 2. Minimap2 arm

Alignment runs as the `minimap2_hg002` DNAnexus applet (BAM lands at `Group4_2026:/results/minimap2/<build>/`). SV calling on that BAM:

```bash
dx download Group4_2026:/scripts/sniffles_mm2.sh
BUILD=hs1 bash sniffles_mm2.sh
```

Output: `Group4_2026:/results/sniffles/minimap2/<build>/`

## Results
