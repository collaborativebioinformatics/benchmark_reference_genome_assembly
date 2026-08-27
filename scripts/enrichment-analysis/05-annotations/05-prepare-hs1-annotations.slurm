#!/usr/bin/env bash
#SBATCH --job-name=svbench-annotation-prepare
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=benchmarking/logs/annotation-prepare_%j.log

set -euo pipefail
bash benchmarking/analysis-export/05-annotations/05-prepare-hs1-annotations.sh
