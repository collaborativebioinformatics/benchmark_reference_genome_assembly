**Table.** HG002 PacBio CCS 15 kb aligned with Winnowmap v2.03 (SV-aware mode, k=15 Meryl
databases from the 2026-08-25 reference-asset release) and called with Sniffles v2.3.2
(--minsvlen 50 --output-rnames, 30 threads).

| Build | Assembly | Coverage | Contigs | Repetitive k-mers | Suppl. aligns | Mapped | SVs called |
|---|---|---|---:|---:|---:|---:|---:|
| hs1  | T2T-CHM13 | 30x (39 cells) | 25  | 66,872 | 68,629  | 99.88% | 36,473 |
| hs1  | T2T-CHM13 | 15x (20 cells) | 25  | 66,872 | 35,598  | 99.88% | 32,910 |
| hg38 | GRCh38    | 15x (20 cells) | 455 | 66,836 | 87,980  | 99.89% | 28,108 |
| hg18 | NCBI36    | 15x (20 cells) | 49  | 66,717 | 193,283 | 99.87% | 25,939 |
| hg19 | GRCh37    | 15x (20 cells) | 93  | 66,667 | 167,871 | 99.88% | 25,840 |

Key observations:
- At matched 15x coverage, SV recovery orders hs1 > hg38 > hg18 ~ hg19. hs1 calls 27% more
  SVs than hg19; hg38 sits between at +9%.
- The coverage effect is smaller than the reference effect: halving depth on hs1 costs ~10%
  of calls (36,473 -> 32,910), while switching hs1 -> hg19 at fixed depth costs 21%.
- Repetitive k-mer counts vary by only 0.3% across builds, confirming equivalent k-mer
  downweighting; differences are not attributable to index parameters.
- Mapping rate is ~99.88% for every build and does not discriminate between references.
- hg18 and hg19 produce 2-5x more supplementary alignments than hs1 or hg38, indicating more
  read fragmentation on the older assemblies. Contig count alone does not predict this:
  hg38 has 455 contigs but fewer supplementary alignments than hg19 (93) or hg18 (49).

Known issue: bioconda resolved Python to 3.14 free-threaded, which breaks Sniffles'
multiprocessing (cannot pickle '_thread._ThreadHandle'); psutil was also absent, causing
Sniffles to start workers and silently write 0 SVs. Both fixed by pinning python=3.11 and
installing psutil. Pin these in run_winnowmap.sh for reproducibility.
