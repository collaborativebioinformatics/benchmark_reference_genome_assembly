#!/usr/bin/env bash
#SBATCH --job-name=svbench-final-audit
#SBATCH --output=logs/svbench-final-audit_%j.log
#SBATCH --error=logs/svbench-final-audit_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

cd "${ANALYSIS_WORKSPACE_ROOT:-$SLURM_SUBMIT_DIR}"
source benchmarking/analysis-export/00-common/00-library.sh
"$BENCHMARK_ENV/bin/python" benchmarking/analysis-export/09-orchestration/00-audit-results.py \
  --root "$BENCHMARK_ROOT"
