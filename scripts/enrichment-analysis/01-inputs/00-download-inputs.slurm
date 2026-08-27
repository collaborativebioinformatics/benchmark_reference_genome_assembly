#!/usr/bin/env bash
#SBATCH --job-name=svbench-download
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=benchmarking/logs/download-inputs_%j.log

set -euo pipefail
bash ./01-inputs/00-download-inputs.sh

