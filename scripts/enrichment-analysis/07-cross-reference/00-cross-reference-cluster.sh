#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/../00-common/00-library.sh"

ARM=$1
ANALYSIS=${2:-strict_full}
case "$ANALYSIS" in
  strict_full)
    ELIGIBILITY=full
    OUTPUT_DIR="$BENCHMARK_ROOT/results/cross-reference/$ARM"
    ;;
  common_endpoints)
    ELIGIBILITY=endpoints
    OUTPUT_DIR="$BENCHMARK_ROOT/results/cross-reference/common-endpoints/$ARM"
    ;;
  *)
    echo "Unknown analysis: $ANALYSIS" >&2
    exit 2
    ;;
esac
TMP_DIR="$BENCHMARK_ROOT/state/tmp/cross-reference/$ANALYSIS.$ARM.${SLURM_JOB_ID:-local}"
mkdir -p "$OUTPUT_DIR" "$TMP_DIR"

TRUVARI=$(find "$BENCHMARK_ENV/bin" -maxdepth 1 -type f -name 'truv*' -print -quit)
[[ -x "$TRUVARI" ]] || { echo "Truvari executable not found" >&2; exit 1; }

"$BENCHMARK_ENV/bin/python" "$SCRIPT_DIR/00-prepare-cross-reference-vcf.py" \
  --root "$BENCHMARK_ROOT" \
  --arm "$ARM" \
  --eligibility "$ELIGIBILITY" \
  --output "$TMP_DIR/combined.unsorted.vcf"

"$BENCHMARK_ENV/bin/bcftools" sort \
  -Oz -o "$OUTPUT_DIR/common-7way.calls.vcf.gz" \
  "$TMP_DIR/combined.unsorted.vcf"
"$BENCHMARK_ENV/bin/tabix" -f -p vcf "$OUTPUT_DIR/common-7way.calls.vcf.gz"

"$TRUVARI" collapse \
  --input "$OUTPUT_DIR/common-7way.calls.vcf.gz" \
  --output "$TMP_DIR/clusters.kept.raw.vcf.gz" \
  --removed-output "$TMP_DIR/clusters.members.raw.vcf.gz" \
  --reference "$HS1_FASTA" \
  --refdist 500 \
  --pctseq 0.70 \
  --pctsize 0.70 \
  --pctovl 0.00 \
  --sizemin 50 \
  --sizemax -1 \
  --passonly

"$BENCHMARK_ENV/bin/bcftools" sort \
  -Oz -o "$OUTPUT_DIR/clusters.kept.vcf.gz" \
  "$TMP_DIR/clusters.kept.raw.vcf.gz"
"$BENCHMARK_ENV/bin/bcftools" sort \
  -Oz -o "$OUTPUT_DIR/clusters.members.vcf.gz" \
  "$TMP_DIR/clusters.members.raw.vcf.gz"
"$BENCHMARK_ENV/bin/tabix" -f -p vcf "$OUTPUT_DIR/clusters.kept.vcf.gz"
"$BENCHMARK_ENV/bin/tabix" -f -p vcf "$OUTPUT_DIR/clusters.members.vcf.gz"
"$BENCHMARK_ENV/bin/python" "$SCRIPT_DIR/00-summarize-cross-reference-clusters.py" \
  --arm "$ARM" \
  --kept "$OUTPUT_DIR/clusters.kept.vcf.gz" \
  --removed "$OUTPUT_DIR/clusters.members.vcf.gz" \
  --clusters "$OUTPUT_DIR/clusters.tsv" \
  --distribution "$OUTPUT_DIR/reference-count-distribution.tsv"
