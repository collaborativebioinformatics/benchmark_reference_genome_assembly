# Supplementary tables

Group 4, HG002 SV reference-benchmarking project. Sample HG002/NA24385,
PacBio HiFi CCS 15 kb, GIAB v5.0q benchmarks. All SV calling used Sniffles2 2.8.0.

## S1_reference_composition.csv
Build composition and Winnowmap index parameters. `contigs` is sequence-record
count; GRCh38's 455 include unplaced scaffolds and alt loci, CHM13v2.0 is
chromosome-level. Repetitive k-mer set sizes vary by only 0.3% across builds,
showing equivalent k-mer downweighting.

## S2_alignment_calling_summary.csv
Alignment and calling summary for the Winnowmap arms, including the CHM13v2.0 30x
versus 15x comparison. Mapping rate is ~99.88% for every build and does not
discriminate between references. Supplementary alignments do: NCBI36 and GRCh37
produce 2-5x more than CHM13v2.0 or GRCh38, indicating greater read fragmentation
on older assemblies, not predicted by contig count.

## S3_armB_inherited_tp_ladder.csv
Arm B read-name tracing, minimap2. Non-anchor calls inherited a TP label where
supporting read names overlapped a CHM13v2.0 Truvari-confirmed true positive
(minimum overlap 0.5, matching SVTYPE). Flat at 37.4-38.5% across 2003-2013
builds. Not a recall; no false-positive estimate for non-anchor builds.

## S4_armC_native_benchmarking.csv
Arm C native benchmarking: each build with a native GIAB v5.0q benchmark scored
against its own benchmark, no projection. Each pair appears unrestricted and
restricted to primary chromosomes; TP/FP/FN identical within every pair, so the
restriction is inert because --includebed already confines scoring to benchmark
regions. Unrefined (truvari refine not applied) - not comparable with Arm A.
benchmark_n differs by build (24,087 / 27,642 / 28,145): three region sets, not
one variant set.

## S5_projection_composition.csv
Class composition before and after projecting the GIAB benchmark VCF between
GRCh38 and CHM13v2.0, for two implementations and four liftOver -minMatch values.
distortion_fold is the fold change in INS:DEL ratio. Insertion survivors are
identical (33,249) in all five CHM13-to-GRCh38 configurations. Deletion survival
varies by implementation but not by -minMatch.

## S6_projection_stage_attribution.csv
Every record assigned one outcome: unmapped, fragmented (>1 target interval),
span_filtered (span changed >10%), or survived. By SVTYPE and length, both
directions. Unmapped is class-neutral; fragmentation and span filtering are not,
and fragmentation rises steeply with deletion length. Records <50 bp excluded,
matching Truvari -s 50.

## S7_tolerance_sweep_DEL_survival.csv
Deletion survival (%) against span tolerance, by length bin. Flat from 10% to 75%,
rising sharply only at 100% where the criterion is inactive: failing records are
grossly rather than marginally mis-projected. Averaged unweighted across both
directions.

## S8_armA_projected_PROVISIONAL.csv
Arm A precision, recall and F1 per source build after projection into CHM13v2.0
coordinates. PROVISIONAL: the six projected builds' values were read from a
plotted figure and are accurate to about +/-0.005. The CHM13v2.0 rows are exact,
matching Arm C and Arm B's anchor. Replace projected rows from the Arm A Truvari
summary, then drop the PROVISIONAL suffix.
