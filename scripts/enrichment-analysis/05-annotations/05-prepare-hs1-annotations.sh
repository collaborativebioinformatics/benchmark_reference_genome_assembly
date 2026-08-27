#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
root=$BENCHMARK_ROOT
output="$root/results/annotations/hs1"
tmp="$root/state/tmp/annotation-prepare.${SLURM_JOB_ID:-local}"
bedtools="$root/.mamba-env/bin/bedtools"
bgzip="$root/.mamba-env/bin/bgzip"
tabix="$root/.mamba-env/bin/tabix"
fai="$root/../references/hs1/hs1.fa.fai"
mkdir -p "$output" "$tmp"

"$root/.mamba-env/bin/python" "$ANALYSIS_EXPORT_ROOT/05-annotations/05-prepare-hs1-annotations.py" \
  --gff "$root/inputs/annotations/gencode/hs1/chm13.draft_v2.0.gene_annotation.gff3" \
  --repeatmasker "$root/inputs/annotations/repeatmasker/hs1/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed" \
  --segdup "$root/inputs/annotations/segdup/hs1/chm13v2.0_SD.bed" \
  --fai "$fai" \
  --output-dir "$tmp"

for track in genes features repeatmasker segdup; do
  "$bedtools" sort -faidx "$fai" -i "$tmp/$track.raw.bed" \
    | uniq > "$tmp/$track.bed"
  "$bgzip" -c "$tmp/$track.bed" > "$output/$track.bed.gz"
  "$tabix" -f -p bed "$output/$track.bed.gz"
done

"$root/tools/ucsc/bigBedToBed" \
  "$root/inputs/annotations/mappability/hs1/k100.Unique.Mappability.bb" \
  "$tmp/mappability-unique-k100.raw.bed"
"$bedtools" sort -faidx "$fai" -i "$tmp/mappability-unique-k100.raw.bed" \
  | cut -f1-3 | uniq > "$tmp/mappability-unique-k100.bed"
"$bgzip" -c "$tmp/mappability-unique-k100.bed" > "$output/mappability-unique-k100.bed.gz"
"$tabix" -f -p bed "$output/mappability-unique-k100.bed.gz"

for track in genes features repeatmasker segdup mappability-unique-k100; do
  records=$(gzip -cd "$output/$track.bed.gz" | wc -l)
  bases=$(gzip -cd "$output/$track.bed.gz" | awk '{sum += $3-$2} END {print sum+0}')
  printf '%s\t%s\t%s\n' "$track" "$records" "$bases"
done > "$output/track-summary.tsv"
sha256sum "$output"/*.bed.gz > "$output/sha256.tsv"
