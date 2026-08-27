#!/usr/bin/env bash
#SBATCH --job-name=svbench-chromhmm-lift
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=benchmarking/logs/chromhmm-liftover_%j.log

set -euo pipefail
bash benchmarking/analysis-export/05-annotations/01-lift-chromhmm-to-hs1.sh
