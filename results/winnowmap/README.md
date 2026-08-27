# Results — Winnowmap arm

Outputs from the Winnowmap → Sniffles2 arm, run on the fixed 20-cell subset
(`config/subset_cells_15x.txt`) plus one full 39-cell run on hs1 for the coverage control.

| File | Contents |
|---|---|
| `sv_results_summary.csv` | Per-run alignment statistics and SV counts (5 runs) |
| `sv_results_summary.md` | Same table with interpretation notes |
| `reference_build_index_stats.csv` | Per-reference contig counts and repetitive k-mer set sizes, measured from the downloaded FASTA files and the 2026-08-25 reference-asset release |
| `reference_build_index_stats.md` | Same table with the caption used in the manuscript |

Raw BAMs, VCFs and SNFs are in the DNAnexus project under
`Group4_2026:/results/winnowmap/<build>/`.

## Runs included

hs1 at 30x (39 cells) and 15x (20 cells); hg38, hg19 and hg18 at 15x.
hg17, hg16 and hg15 are not yet run on this arm.

## Reproducibility notes

Two environment constraints are required, neither of which is enforced by the
conda solve:

1. **Pin `python=3.11`** (or 3.10). The default solve resolves to Python 3.14
   free-threaded, under which Sniffles2's multiprocessing fails immediately with
   `TypeError: cannot pickle '_thread._ThreadHandle' object`.
2. **Install `psutil`.** It is not a declared dependency of the Sniffles2 conda
   package, but without it Sniffles2 starts its workers, exits with status 0 and
   writes an empty VCF containing headers only. This failure produces no error
   message — always check the SV count in the final log line before uploading.

Working invocation used for the re-called runs:

```bash
micromamba create -y -p ~/snf -c conda-forge -c bioconda "python=3.11" sniffles=2.3.2
~/snf/bin/pip install psutil
~/snf/bin/sniffles --input HG002.<build>.winnowmap.bam --reference <build>.fa \
  --vcf HG002.<build>.winnowmap.sniffles.vcf.gz \
  --snf HG002.<build>.winnowmap.snf \
  --minsvlen 50 --output-rnames --threads 30 --allow-overwrite
```
