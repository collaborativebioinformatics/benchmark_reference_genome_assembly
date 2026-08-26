# Raw data

## PacBio CCS reads (HG002)
- PacBio CCS 15 kb (~30x coverage)
- Data link: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/
- Location: `Group4_2026:/Data/CCS_15kb_30x/` — 39 SMRT-cell FASTQs (`*.Q20.fastq`), ~166.2 GB total, ~30x coverage.
- Subset actually used: 20 of 39 cells (every second cell after sorting movie names) = 86.4 GB, 52.0% of bases, ≈15x coverage. Exact list: [`config/subset_cells_15x.txt`](../config/subset_cells_15x.txt).
- Source: not recorded in DNAnexus metadata. The files were uploaded manually via a JupyterLab session (job `JupyterLab_data_upload_test`, run by `user-yunjialiusv`) rather than an automated download, so there's no stored source URL to restore. The movie-ID filenames (e.g. `m54238_180901_011437.Q20.fastq`) follow the naming convention of the public GIAB HG002 PacBio CCS 15kb dataset, but that's an inference from naming, not a confirmed source — worth checking with whoever ran the upload if exact provenance is needed.
- Data transfer: Due to the large number and size of the files, the datasets were first downloaded within a DNAnexus JupyterLab/workstation session and then uploaded to the shared project directory. 

## Reference genome assemblies

Seven UCSC builds spanning ~20 years of assembly quality, downloaded from UCSC goldenPath and repackaged (with per-build aligner indices) at `Group4_2026:/reference-assets/releases/2026-08-25/`:

| build | assembly | source |
|---|---|---|
| hs1  | T2T-CHM13 v2.0 | https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz |
| hg38 | GRCh38 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz |
| hg19 | GRCh37 | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz |
| hg18 | NCBI36 | https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/hg18.fa.gz |
| hg17 | NCBI35 | https://hgdownload.soe.ucsc.edu/goldenPath/hg17/bigZips/hg17.fa.gz |
| hg16 | NCBI34 | https://hgdownload.soe.ucsc.edu/goldenPath/hg16/bigZips/hg16.fa.gz |
| hg15 | NCBI33 (April 2003) | https://hgdownload.soe.ucsc.edu/goldenPath/hg15/bigZips/chromFa.zip |

Source table: [`config/genomes.tsv`](../config/genomes.tsv) (pulled from the release's provenance archive, `Group4_2026:/reference-assets/releases/2026-08-25/provenance/build-workflow-scripts.tar.zst`). Re-download any build directly with [`scripts/download_reference.sh`](../scripts/download_reference.sh):

```bash
scripts/download_reference.sh hs1
```

This is a standalone version of the original `download_reference.sh` from that provenance archive, with the SLURM/array-job plumbing stripped out. Each reference tarball in DNAnexus also has a `.sha256` sidecar for integrity verification.

Note: minimap2-arm results exist for all seven builds; the winnowmap arm is currently missing an `hg15` result (has `hg16`–`hg19`, `hg38`, `hs1`, plus an `hs1_15x` variant) — worth confirming whether that run is still pending.


## Reference Genome Assemblies

**Table 1. Reference genome assemblies and their summary statistics.**

| Reference (UCSC / build) | Release date | Total seq length (bp) | Contig N50 (bp) | Scaffold N50 (bp) | No. contigs | Accession (GenBank / RefSeq) |
|---|---|---:|---:|---:|---:|---|
| **T2T-CHM13 v2.0 (hs1)** | Jan 2022 | 3,117,275,501 | 150,617,247 | 150,617,247 | 24 | GCF_009914755.1 |
| **GRCh38.p14 (hg38)** | Feb 2022 (p14) | 3,099,441,038 | 57,879,411 | 67,794,873 | 996 | GCF_000001405.40 |
| **GRCh38 initial (hg38)** | Dec 2013 | 3,099,734,149 | 57,879,411 | 67,794,873 | 999 | GCF_000001405.26 |
| **GRCh37 (hg19)** | Feb 2009 | 3,101,788,170 | 38,508,932 | 46,395,641 | 350 | GCF_000001405.13 |
| **NCBI36 (hg18)** | Mar 2006 | 3,093,104,542 | 38,509,590 | 38,509,590 | 1,006 | GCF_000001405.12 |
| **NCBI35 (hg17)** | May 2004 | 3,091,360,260 | 37,760,040 | 38,509,590 | 1,215 | GCF_000001405.11 |
| **NCBI34 (hg16)** | Jul 2003 | 3,091,959,510 | 28,857,747 | 29,104,798 | 1,756 | GCF_000001405.10 |
| **NCBI33 (hg15)** | Apr 2003 | 3,095,784,245 | 23,437,594 | 25,443,670 | 2,240 | GCF_000001405.8 |

### Notes

- All assembly statistics reported in Table 1 were obtained from **NCBI Datasets assembly reports** using the same set of assembly metrics across all reference builds. This allows contig number, total sequence length, contig N50, and scaffold N50 to be compared using consistent definitions.
- For **T2T-CHM13 v2.0**, contig N50 and scaffold N50 are identical because the assembly is gap-free at chromosome scale.
- RefSeq (`GCF_`) accessions are reported in the table for consistency. Corresponding GenBank (`GCA_`) accessions can also be included if both identifiers are required.
- UCSC build names and release dates were cross-referenced with the **UCSC Genome Browser release history**.
- Assembly statistics may differ depending on the exact sequence set used (e.g., primary assembly versus additional alternate/unplaced sequences). The accession listed in the table should therefore be used to identify the exact assembly version.

### Sources and References

**Assembly statistics**  
NCBI Datasets, Genome. Per-accession assembly reports:  
https://www.ncbi.nlm.nih.gov/datasets/genome/

**UCSC build names and release dates**  
UCSC Genome Browser, Release Log:  
https://genome.ucsc.edu/goldenPath/releaseLog.html

**T2T-CHM13 v2.0**  
Nurk S, Koren S, Rhie A, Rautiainen M, et al. The complete sequence of a human genome. *Science*. 2022;376(6588):44–53. doi:10.1126/science.abj6987. Assembly: `GCA_009914755.4`.

**GRCh38**  
Schneider VA, Graves-Lindsay T, Howe K, et al. Evaluation of GRCh38 and de novo haploid genome assemblies demonstrates the enduring quality of the reference assembly. *Genome Research*. 2017;27(5):849–864. doi:10.1101/gr.213611.116. Assembly: `GCA_000001405.15` (RefSeq `GCF_000001405.26` / `GCF_000001405.40`).

**GRCh37 and the GRC assembly model**  
Church DM, Schneider VA, Graves T, et al. Modernizing reference genome assemblies. *PLoS Biology*. 2011;9(7):e1001091. doi:10.1371/journal.pbio.1001091. Assembly: `GCA_000001405.1`.

**Earlier human reference builds (NCBI33–36)**  
International Human Genome Sequencing Consortium. Finishing the euchromatic sequence of the human genome. *Nature*. 2004;431(7011):931–945. doi:10.1038/nature03001.

**Definition of N50 and related assembly metrics**  
Earl D, Bradnam K, St John J, et al. Assemblathon 1: a competitive assessment of de novo short read assembly methods. *Genome Research*. 2011;21(12):2224–2241. doi:10.1101/gr.126599.111.
