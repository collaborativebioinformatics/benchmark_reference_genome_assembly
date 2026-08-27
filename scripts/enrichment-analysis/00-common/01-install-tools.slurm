#!/usr/bin/env bash
#SBATCH --job-name=svbench-install
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/install-tools_%j.log

set -euo pipefail
bash benchmarking/analysis-export/00-common/01-install-tools.sh

