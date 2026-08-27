#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
root=$BENCHMARK_ROOT
source_bed_gz="$root/inputs/annotations/chromhmm/grch38/E116_15_coreMarks_hg38lift_mnemonics.bed.gz"
chain="$root/../liftover/chains/hg38/hs1/hg38ToHs1.over.chain.gz"
output_dir="$root/results/annotations/chromhmm/hs1"
output="$output_dir/E116_GM12878_15state.hs1.bed.gz"
failed="$output_dir/E116_GM12878_15state.hg38-to-hs1.failed-bins.bed.gz"
distorted="$output_dir/E116_GM12878_15state.hg38-to-hs1.distorted-bins.bed.gz"
ambiguous="$output_dir/E116_GM12878_15state.hg38-to-hs1.ambiguous-overlaps.bed.gz"
tmp="$root/state/tmp/chromhmm-liftover.${SLURM_JOB_ID:-local}"
transanno="$root/../tools/liftover/bin/transanno"
bedtools="$root/.mamba-env/bin/bedtools"
bgzip="$root/.mamba-env/bin/bgzip"
tabix="$root/.mamba-env/bin/tabix"
mkdir -p "$output_dir" "$tmp"

"$root/.mamba-env/bin/python" "$ANALYSIS_EXPORT_ROOT/05-annotations/01-split-chromhmm-bins.py" \
  --input "$source_bed_gz" \
  --output "$tmp/source.200bp.bed" \
  --bin-size 200
"$transanno" liftbed \
  --chain "$chain" \
  --output "$tmp/mapped.unsorted.bed" \
  --failed "$tmp/failed.bed" \
  "$tmp/source.200bp.bed"
"$root/.mamba-env/bin/python" "$ANALYSIS_EXPORT_ROOT/05-annotations/01-filter-chromhmm-lift.py" \
  --input "$tmp/mapped.unsorted.bed" \
  --accepted "$tmp/mapped.accepted.unsorted.bed" \
  --rejected "$tmp/mapped.distorted.bed" \
  --minimum-ratio 0.80 \
  --maximum-ratio 1.20
"$bedtools" sort \
  -faidx "$root/../references/hs1/hs1.fa.fai" \
  -i "$tmp/mapped.accepted.unsorted.bed" > "$tmp/mapped.accepted.sorted.bed"
"$root/.mamba-env/bin/python" "$ANALYSIS_EXPORT_ROOT/05-annotations/01-merge-chromhmm-state-bins.py" \
  --input "$tmp/mapped.accepted.sorted.bed" \
  --output "$tmp/mapped.merged.bed" \
  --ambiguous "$tmp/mapped.ambiguous.bed"
"$bgzip" -c "$tmp/mapped.merged.bed" > "$tmp/output.bed.gz"
mv "$tmp/output.bed.gz" "$output"
"$bgzip" -c "$tmp/failed.bed" > "$failed"
"$bgzip" -c "$tmp/mapped.distorted.bed" > "$distorted"
"$bgzip" -c "$tmp/mapped.ambiguous.bed" > "$ambiguous"
"$tabix" -f -p bed "$output"

source_records=$(wc -l < "$tmp/source.200bp.bed")
mapped_bin_records=$(wc -l < "$tmp/mapped.accepted.unsorted.bed")
mapped_records=$(wc -l < "$tmp/mapped.merged.bed")
failed_records=$(gzip -cd "$failed" | wc -l)
distorted_records=$(gzip -cd "$distorted" | wc -l)
ambiguous_records=$(gzip -cd "$ambiguous" | wc -l)
source_bases=$(awk '{sum += $3-$2} END {print sum+0}' "$tmp/source.200bp.bed")
mapped_bases=$(gzip -cd "$output" | awk '{sum += $3-$2} END {print sum+0}')
failed_bases=$(gzip -cd "$failed" | awk '{sum += $3-$2} END {print sum+0}')
distorted_bases=$(gzip -cd "$distorted" | awk '{sum += $3-$2} END {print sum+0}')
ambiguous_bases=$(gzip -cd "$ambiguous" | awk '{sum += $3-$2} END {print sum+0}')
{
  printf 'source_bins\tmapped_bins\tmerged_intervals\tfailed_bins\tdistorted_bins\tambiguous_intervals\tsource_bases\tmapped_bases\tfailed_bases\tdistorted_target_bases\tambiguous_target_bases\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$source_records" "$mapped_bin_records" "$mapped_records" "$failed_records" "$distorted_records" \
    "$ambiguous_records" "$source_bases" "$mapped_bases" "$failed_bases" "$distorted_bases" "$ambiguous_bases"
} > "$output_dir/liftover-summary.tsv"
