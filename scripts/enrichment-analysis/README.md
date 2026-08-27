# HG002 SV analysis scripts

This is the comment-cleaned export of the workflow used to download, preprocess,
project, benchmark, annotate, compare, and plot the 14 HG002 Sniffles callsets.
Scripts are grouped by pipeline stage and numbered locally within each directory.

Run SLURM submissions from the workspace root containing `benchmarking/`,
`references/`, and `liftover/`. Partition directives are intentionally omitted.
Required input and parameter tables are in `config/`.

Prepare the analysis environment and annotation tracks with `00-common/`,
`05-annotations/`, and `08-figures/00-install-plot-tools.slurm`. The main
submitter includes the DNANexus input stage. Then run:

```bash
bash 09-orchestration/01-submit-analysis.sh
```

Results are written to `benchmarking/results/`.
