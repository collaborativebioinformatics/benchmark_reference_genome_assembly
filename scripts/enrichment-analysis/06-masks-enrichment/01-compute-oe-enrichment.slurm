#!/usr/bin/env bash
#SBATCH --job-name=svbench-enrichment
#SBATCH --output=logs/svbench-enrichment_%j.log
#SBATCH --error=logs/svbench-enrichment_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

cd "${ANALYSIS_WORKSPACE_ROOT:-$SLURM_SUBMIT_DIR}"
source benchmarking/analysis-export/00-common/00-library.sh
"$BENCHMARK_ENV/bin/python" benchmarking/analysis-export/06-masks-enrichment/01-compute-oe-enrichment.py \
  --root "$BENCHMARK_ROOT"
