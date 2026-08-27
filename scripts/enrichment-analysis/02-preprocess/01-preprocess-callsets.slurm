#!/usr/bin/env bash
#SBATCH --job-name=svbench-preprocess
#SBATCH --array=1-14%7
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/preprocess_%A_%a.log

set -euo pipefail
bash benchmarking/analysis-export/02-preprocess/01-preprocess-callset.sh "$SLURM_ARRAY_TASK_ID"

