#!/usr/bin/env bash
#SBATCH --job-name=svbench-annotation-finalize
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/annotation-finalize_%j.log

set -euo pipefail
benchmarking/.mamba-env/bin/python benchmarking/analysis-export/05-annotations/07-finalize-annotations.py \
  --root benchmarking \
  --callsets benchmarking/analysis-export/config/callsets.tsv \
  --output benchmarking/results/annotations/sv-annotations.tsv.gz \
  --summary benchmarking/results/annotations/annotation-row-counts.tsv
