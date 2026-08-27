#!/usr/bin/env bash
#SBATCH --job-name=svbench-ucsc-tools
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/ucsc-tools_%j.log

set -euo pipefail
bash benchmarking/analysis-export/05-annotations/04-install-ucsc-tools.sh
