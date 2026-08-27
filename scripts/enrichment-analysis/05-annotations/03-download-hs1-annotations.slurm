#!/usr/bin/env bash
#SBATCH --job-name=svbench-annotation-download
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarking/logs/annotation-download_%j.log

set -euo pipefail
bash benchmarking/analysis-export/05-annotations/03-download-hs1-annotations.sh
