**Table 1.** Reference build composition and Winnowmap index parameters.
Alignments used Winnowmap v2.03 in SV-aware mode with pre-computed Meryl v1.4.2
k-mer databases (k = 15) from the project reference-asset release (2026-08-25).

| Build | Assembly | Contigs | Repetitive k-mers (k=15) |
|---|---|---:|---:|
| hs1  | T2T-CHM13 | 25  | 66,872 |
| hg38 | GRCh38    | 455 | 66,836 |
| hg19 | GRCh37    | 93  | 66,667 |
| hg18 | NCBI36    | 49  | 66,717 |

Repetitive k-mer counts vary by only 0.3% across builds, confirming equivalent
k-mer downweighting; reference composition differs substantially (hg38 includes
unplaced scaffolds and alt loci; hs1 is chromosome-level with none).
hs1 was aligned at both 30x (39 SMRT cells) and 15x (20 cells) using the same index.
