#!/usr/bin/env bash
#SBATCH --job-name=svbench-cluster-relaxed
#SBATCH --output=logs/svbench-cluster-relaxed_%A_%a.log
#SBATCH --error=logs/svbench-cluster-relaxed_%A_%a.err
#SBATCH --array=1-2
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd "${ANALYSIS_WORKSPACE_ROOT:-$SLURM_SUBMIT_DIR}"
case "$SLURM_ARRAY_TASK_ID" in
  1) ARM=minimap2 ;;
  2) ARM=winnowmap ;;
  *) exit 2 ;;
esac
bash benchmarking/analysis-export/07-cross-reference/00-cross-reference-cluster.sh "$ARM" common_endpoints
