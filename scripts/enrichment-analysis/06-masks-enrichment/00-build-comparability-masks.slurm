#!/usr/bin/env bash
#SBATCH --job-name=svbench-masks
#SBATCH --output=logs/svbench-masks_%j.log
#SBATCH --error=logs/svbench-masks_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

cd "${ANALYSIS_WORKSPACE_ROOT:-$SLURM_SUBMIT_DIR}"
bash benchmarking/analysis-export/06-masks-enrichment/00-build-comparability-masks.sh
