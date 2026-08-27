#!/usr/bin/env bash
#SBATCH --job-name=svbench-relaxed-oe
#SBATCH --output=logs/svbench-relaxed-oe_%j.log
#SBATCH --error=logs/svbench-relaxed-oe_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G

set -euo pipefail

cd "${ANALYSIS_WORKSPACE_ROOT:-$SLURM_SUBMIT_DIR}"
source benchmarking/analysis-export/00-common/00-library.sh
"$BENCHMARK_ENV/bin/python" benchmarking/analysis-export/06-masks-enrichment/03-compute-relaxed-enrichment.py \
  --root "$BENCHMARK_ROOT"
