#!/usr/bin/env bash
#SBATCH --job-name=svbench-chromhmm-dl
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/chromhmm-download_%j.log

set -euo pipefail
bash benchmarking/analysis-export/05-annotations/00-download-chromhmm.sh

