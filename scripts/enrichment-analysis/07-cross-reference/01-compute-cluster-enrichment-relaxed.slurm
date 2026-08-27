#!/usr/bin/env bash
#SBATCH --job-name=svbench-cluster-oe-relaxed
#SBATCH --output=logs/svbench-cluster-oe-relaxed_%j.log
#SBATCH --error=logs/svbench-cluster-oe-relaxed_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

cd "${ANALYSIS_WORKSPACE_ROOT:-$SLURM_SUBMIT_DIR}"
source benchmarking/analysis-export/00-common/00-library.sh
"$BENCHMARK_ENV/bin/python" benchmarking/analysis-export/07-cross-reference/01-compute-cluster-enrichment.py \
  --root "$BENCHMARK_ROOT" \
  --analysis common_endpoints
