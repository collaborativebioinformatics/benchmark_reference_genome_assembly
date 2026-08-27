# Results — Winnowmap arm

Alignment statistics and structural variant counts for the Winnowmap → Sniffles2 arm,
run on the fixed 20-cell subset (`config/subset_cells_15x.txt`) plus one full 39-cell
run on hs1 as a coverage control.

| File | Contents |
|---|---|
| `sv_results_by_caller.csv` | Per-run alignment statistics and SV counts under both Sniffles2 versions |
| `reference_build_index_stats.csv` | Per-reference contig counts and repetitive k-mer set sizes, measured from the downloaded FASTA files and the 2026-08-25 reference-asset release |

BAMs, VCFs and SNFs are in the DNAnexus project under
`Group4_2026:/results/winnowmap/<build>/`. Both callers' VCFs are retained; the
2.8.0 outputs carry a `sniffles280` infix.

## Runs

hs1 at 30x (39 cells) and 15x (20 cells); hg38, hg19 and hg18 at 15x.
hg17, hg16 and hg15 are not yet run on this arm.

## Reference composition

| Build | Assembly | Contigs | Repetitive k-mers (k=15) |
|---|---|---:|---:|
| hs1  | T2T-CHM13 | 25  | 66,872 |
| hg38 | GRCh38    | 455 | 66,836 |
| hg19 | GRCh37    | 93  | 66,667 |
| hg18 | NCBI36    | 49  | 66,717 |

Repetitive k-mer counts vary by only 0.3% across builds, confirming equivalent k-mer
downweighting — differences in the callsets are not attributable to index parameters.
Contig counts are as measured from the UCSC FASTA files and differ from NCBI assembly
reports, which count scaffolds differently.

## SV counts under both callers

| Run | Sniffles2 2.3.2 | Sniffles2 2.8.0 | Change |
|---|---:|---:|---:|
| hs1 30x  | 36,473 | 29,444 | −19.3% |
| hs1 15x  | 32,910 | 25,931 | −21.2% |
| hg38 15x | 28,108 | 23,294 | −17.1% |
| hg18 15x | 25,939 | 23,484 | −9.5% |
| hg19 15x | 25,840 | 23,100 | −10.6% |

Mapping rate is ~99.88% for every build and does not discriminate between references.
Supplementary alignment counts do: hg18 (193,283) and hg19 (167,871) fragment reads
2–5× more than hs1 (35,598) or hg38 (87,980) at matched coverage. Contig count alone
does not predict this — hg38 has 455 contigs but fewer supplementary alignments than
hg19 (93) or hg18 (49).

## The caller version changes the conclusions

Sniffles2 2.8.0 calls 10–21% fewer SVs than 2.3.2 on identical alignments, and the
SNF candidate counts roughly halve, so the difference arises upstream of final
reporting. This is not a uniform offset, and it materially changes three findings:

| Effect (matched 15x) | Under 2.3.2 | Under 2.8.0 |
|---|---:|---:|
| hs1 vs hg19 | +27.4% | +12.3% |
| hg38 vs hg19 | +8.8% | +0.8% |
| Coverage (hs1 30x → 15x) | −9.8% | −11.9% |

Consequences:

- The reference effect roughly halves under 2.8.0.
- "Reference quality matters more than sequencing depth" holds under 2.3.2 (27% vs 10%)
  but **not** under 2.8.0, where the two effects are comparable (12% vs 12%).
- hg38's advantage over the older builds collapses. Under 2.8.0, hg38, hg18 and hg19
  fall within 1.7% of each other and are effectively indistinguishable.

What survives both versions: **hs1 leads clearly**, and mapping rate does not
discriminate between references.

Any comparison that mixes caller versions is confounded. Note that in the wider study
the minimap2 arm and the Winnowmap hg15–hg17 runs used 2.8.0, so cross-arm comparisons
should use the 2.8.0 columns here.

## Reproducibility notes

Two environment constraints are required, neither enforced by the conda solve:

1. **Pin `python=3.11`** (or 3.10). The default solve resolves to Python 3.14
   free-threaded, under which Sniffles2 2.3.2's multiprocessing fails immediately with
   `TypeError: cannot pickle '_thread._ThreadHandle' object`.
2. **Install `psutil`.** It is not a declared dependency of the Sniffles2 conda package,
   but without it Sniffles2 2.3.2 starts its workers, exits with status 0 and writes an
   empty VCF containing headers only. This failure produces no error message — always
   check the SV count in the final log line before uploading.

Working invocation:

```bash
micromamba create -y -p ~/snf -c conda-forge -c bioconda "python=3.11" sniffles=2.8.0
~/snf/bin/pip install psutil
~/snf/bin/sniffles --input HG002.<build>.winnowmap.bam --reference <build>.fa \
  --vcf HG002.<build>.winnowmap.sniffles280.vcf.gz \
  --snf HG002.<build>.winnowmap.sniffles280.snf \
  --minsvlen 50 --output-rnames --threads 30
```
