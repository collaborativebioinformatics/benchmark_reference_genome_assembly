#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
load_task "${1:?task_id is required}"

result="$BENCHMARK_ROOT/results/TRUVARI/$ARM/$BUILD"
annotation_dir="$BENCHMARK_ROOT/results/annotations/hs1"
output_dir="$BENCHMARK_ROOT/results/annotations/variants/$ARM/$BUILD"
tmp="$BENCHMARK_ROOT/state/tmp/variant-annotation/$ARM.$BUILD.${SLURM_JOB_ID:-local}"
mkdir -p "$output_dir" "$tmp"

"$BENCHMARK_ENV/bin/python" "$ANALYSIS_EXPORT_ROOT/05-annotations/06-variants-to-bed.py" \
  --result-dir "$result" \
  --arm "$ARM" \
  --build "$BUILD" \
  --bed "$tmp/variants.bed" \
  --metadata "$tmp/metadata.tsv"
"$BENCHMARK_ROOT/tools/ucsc/bigWigAverageOverBed" \
  "$BENCHMARK_ROOT/inputs/annotations/mappability/hs1/k100.Umap.MultiTrackMappability.bw" \
  "$tmp/variants.bed" \
  "$tmp/mappability.tsv"
"$BENCHMARK_ENV/bin/python" "$ANALYSIS_EXPORT_ROOT/05-annotations/06-annotate-variants.py" \
  --metadata "$tmp/metadata.tsv" \
  --mappability "$tmp/mappability.tsv" \
  --annotation-dir "$annotation_dir" \
  --chromhmm "$BENCHMARK_ROOT/results/annotations/chromhmm/hs1/E116_GM12878_15state.hs1.bed.gz" \
  --output "$tmp/annotations.tsv.gz"
mv "$tmp/annotations.tsv.gz" "$output_dir/annotations.tsv.gz"
