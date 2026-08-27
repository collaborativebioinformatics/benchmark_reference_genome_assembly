#!/usr/bin/env bash
#SBATCH --job-name=svbench-summary
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/summary_%j.log

set -euo pipefail
benchmarking/.mamba-env/bin/python benchmarking/analysis-export/04-benchmark/01-summarize-benchmarks.py \
  --root benchmarking \
  --callsets benchmarking/analysis-export/config/callsets.tsv \
  --output-dir benchmarking/results/summary

