#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
load_task "${1:?task_id is required}"

python="$BENCHMARK_ENV/bin/python"
input="$BENCHMARK_ROOT/results/lifted/$ARM/$BUILD/$ARM.$BUILD.to-hs1.vcf.gz"
output="$BENCHMARK_ROOT/results/TRUVARI/$ARM/$BUILD"
tmp_root="$BENCHMARK_ROOT/state/tmp/TRUVARI/$ARM.$BUILD.${SLURM_JOB_ID:-local}"
tmp_output="$tmp_root/bench"

require_file "$input"
require_file "$input.tbi"
require_file "$TRUTH_VCF"
require_file "$TRUTH_VCF.tbi"
require_file "$TRUTH_BED"
require_file "$HS1_FASTA"

if [[ -s "$output/summary.json" ]]; then
  echo "EXISTS $output/summary.json"
  exit 0
fi
if [[ -e "$output" ]]; then
  echo "Incomplete output directory already exists: $output" >&2
  exit 1
fi
mkdir -p "$tmp_root" "${output%/*}"

"$python" -m truvari bench \
  --base "$TRUTH_VCF" \
  --comp "$input" \
  --output "$tmp_output" \
  --reference "$HS1_FASTA" \
  --write-resolved \
  --refdist 500 \
  --pctseq 0.70 \
  --pctsize 0.70 \
  --pctovl 0.00 \
  --pick single \
  --sizemin 50 \
  --sizefilt 30 \
  --sizemax -1 \
  --includebed "$TRUTH_BED"

mv "$tmp_output" "$output"
ln "$output/tp-comp.vcf.gz" "$output/TP.vcf.gz"
ln "$output/tp-comp.vcf.gz.tbi" "$output/TP.vcf.gz.tbi"
ln "$output/fp.vcf.gz" "$output/FP.vcf.gz"
ln "$output/fp.vcf.gz.tbi" "$output/FP.vcf.gz.tbi"
ln "$output/fn.vcf.gz" "$output/FN.vcf.gz"
ln "$output/fn.vcf.gz.tbi" "$output/FN.vcf.gz.tbi"
ln "$output/tp-base.vcf.gz" "$output/TP.truth.vcf.gz"
ln "$output/tp-base.vcf.gz.tbi" "$output/TP.truth.vcf.gz.tbi"
