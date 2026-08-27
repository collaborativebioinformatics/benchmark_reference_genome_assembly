#!/usr/bin/env bash
#SBATCH --job-name=svbench-chromhmm-audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/chromhmm-audit_%j.log

set -euo pipefail
benchmarking/.mamba-env/bin/python benchmarking/analysis-export/05-annotations/02-audit-chromhmm.py \
  --input benchmarking/results/annotations/chromhmm/hs1/E116_GM12878_15state.hs1.bed.gz \
  --fai references/hs1/hs1.fa.fai \
  --output benchmarking/results/annotations/chromhmm/hs1/audit.tsv
