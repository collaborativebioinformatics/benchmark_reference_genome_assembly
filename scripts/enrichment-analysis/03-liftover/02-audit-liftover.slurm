#!/usr/bin/env bash
#SBATCH --job-name=svbench-lift-audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/liftover-audit_%j.log

set -euo pipefail
benchmarking/.mamba-env/bin/python benchmarking/analysis-export/03-liftover/02-audit-liftover.py \
  --root benchmarking \
  --callsets benchmarking/analysis-export/config/callsets.tsv \
  --output-prefix benchmarking/state/liftover-audit

