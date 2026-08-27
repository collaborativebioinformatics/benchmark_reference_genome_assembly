#!/usr/bin/env bash
#SBATCH --job-name=svbench-truvari
#SBATCH --array=1-14%4
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarking/logs/TRUVARI_%A_%a.log

set -euo pipefail
bash benchmarking/analysis-export/04-benchmark/00-benchmark-callset.sh "$SLURM_ARRAY_TASK_ID"

