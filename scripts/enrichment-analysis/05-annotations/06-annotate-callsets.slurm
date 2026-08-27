#!/usr/bin/env bash
#SBATCH --job-name=svbench-annotate
#SBATCH --array=1-14%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarking/logs/annotate_%A_%a.log

set -euo pipefail
bash benchmarking/analysis-export/05-annotations/06-annotate-callset.sh "$SLURM_ARRAY_TASK_ID"
