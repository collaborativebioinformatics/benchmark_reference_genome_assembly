#!/usr/bin/env bash
#SBATCH --job-name=svbench-liftover
#SBATCH --array=1-14%7
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=benchmarking/logs/liftover_%A_%a.log

set -euo pipefail
bash benchmarking/analysis-export/03-liftover/01-lift-callset-to-hs1.sh "$SLURM_ARRAY_TASK_ID"

