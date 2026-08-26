# Benchmarking the Impact of Reference Genome Quality on SV Discovery

---

## Research Focus
> How much of the SV callset changes as the reference assembly improves from an early draft (hg15) to a complete telomere-to-telomere assembly (T2T-CHM13) and what fraction of that change can be measured directly?

---

## Background

Structural variant detection is always relative to a reference assembly, and while humans now have a complete telomere-to-telomere reference, most species are stuck with fragmented draft assemblies. Researchers calling SVs against such drafts have no principled way to know what their callset is missing. Therefore, we would like to systematically investigate how assembly quality influences SV calling.

---

## Methods  
**Reference Genome Assemblies**  
7 reference genome assemblies, spanning the history of the human genome reference from its first public draft to the current telomere-to-telomere assembly, were used:
- hg15 (2003)
- hg16 (NCBI34, 2003)
- hg17 (NCBI35, 2004)
- hg18 (NCBI36, 2006)
- hg19 (GRCh37, 2009)
- hg38 (GRCh38, 2013)
- hs1 (T2T-CHM13 v2.0, 2022)

**Benchmark Sample**   
Our current workflow only consists of the benchmark sample HG002 from the PacBio CCS long-reads 15kb (~30x coverage). Which was then subsampled deterministically prior to alignment, to ~15x rather than the full 30x, primarily to reduce runtime across 7 reference builds x 2 aligners.

**Read Alignment**  
Reads were aligned to each reference build independently using 2 long-read aligners, run in parallel: **Winnowmap v2.03** and **minimap2**.  

- **Winnowmap**: a repetitive k-mer list was pre-computed for each reference build with Meryl v1.4.2 and supplied via `-W` flag and alignment was run with `-ax map-pb` preset
- **minimap2**: alignment was run with the `-ax map-hifi` preset, intended for CCS/HiFi data  

Both aligners were run genome-wide to avoid displacing reads that align more accurately elsewhere in the genome onto the target regions, which would display false structural-variant signal. Alignments were coordinate-sorted and indexed with samtools v1.24.  

**Structural Variant Calling - Sniffles2 v2.3.2**  
Structural variants were called from each whole-genome alignment using Sniffles2 v2.3.2 using the following parameter: `--minsvlen 50`, and produced 1 VCF and 1 SNF file per reference build x aligner combination.  

**Benchmarking Against Truth Sets**  
[In works]

## Solution

## Workflow Diagram [in-works]

![Workflow Illustration](group4_flowchart_diagram.png)

## Results

## Members
| Team member | Role |
| :--- | :--- |
| **Susanne Pfeifer** | Group Leader |
| Alisa Iakupova | xxx |
| Daniil Khlebnikov | xxx |
| Gerald McCollam | xxx |
| German Demidov | xxx |
| Miguel Angel Trejo Acosta | xxx |
| Nirajan Bhattarai | xxx |
| Steven Chen | xxx |
| Tammy Sisodiya | xxx |
| Vy Dang | xxx |
| Yunkia Liu | xxx |

## References
