#!/usr/bin/env bash

set -euo pipefail

ANALYSIS_EXPORT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BENCHMARK_ROOT=$(cd "$ANALYSIS_EXPORT_ROOT/.." && pwd)
BENCHMARK_ENV="$BENCHMARK_ROOT/.mamba-env"
CALLSETS_TSV="$ANALYSIS_EXPORT_ROOT/config/callsets.tsv"
HS1_FASTA="$BENCHMARK_ROOT/../references/hs1/hs1.fa"
TRUTH_VCF="$BENCHMARK_ROOT/inputs/truth/hs1/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz"
TRUTH_BED="$BENCHMARK_ROOT/inputs/truth/hs1/HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed"

require_file() {
  [[ -s "$1" ]] || { echo "Missing or empty file: $1" >&2; exit 1; }
}

load_task() {
  local task_id=$1
  local row
  row=$(awk -F '\t' -v id="$task_id" 'NR > 1 && $1 == id {print; exit}' "$CALLSETS_TSV")
  [[ -n "$row" ]] || { echo "Unknown task_id: $task_id" >&2; exit 1; }
  IFS=$'\t' read -r TASK_ID ARM BUILD SNIFFLES_VERSION INPUT_VCF SOURCE_FASTA CHAIN_TO_HS1 <<< "$row"
  INPUT_VCF="$BENCHMARK_ROOT/$INPUT_VCF"
  SOURCE_FASTA="$BENCHMARK_ROOT/$SOURCE_FASTA"
  if [[ "$CHAIN_TO_HS1" != "." ]]; then
    CHAIN_TO_HS1="$BENCHMARK_ROOT/$CHAIN_TO_HS1"
  fi
}
