#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/00-library.sh"
root=$BENCHMARK_ROOT
env_dir="$root/.mamba-env"
mamba_bin=/home/dkhlebnikov/miniforge3/bin/mamba
liftover_commit=e62dcbf636b4370c6ec33b04836ea1de091eafdb
liftover_dir="$root/tools/liftoverSV"

if [[ -d "$env_dir/conda-meta" ]]; then
  "$mamba_bin" env update -p "$env_dir" -f "$ANALYSIS_EXPORT_ROOT/00-common/environment.yml" --prune
else
  "$mamba_bin" env create -p "$env_dir" -f "$ANALYSIS_EXPORT_ROOT/00-common/environment.yml"
fi

mkdir -p "$root/tools"
if [[ -d "$liftover_dir/.git" ]]; then
  git -C "$liftover_dir" fetch origin
else
  git clone https://github.com/lgmgeo/liftoverSV.git "$liftover_dir"
fi
git -C "$liftover_dir" checkout --detach "$liftover_commit"

{
  "$env_dir/bin/python" --version
  "$env_dir/bin/bcftools" --version | head -1
  "$env_dir/bin/samtools" --version | head -1
  "$env_dir/bin/python" -c 'import truvari; print("truvari", truvari.__version__)'
  /home/dkhlebnikov/miniforge3/envs/jasmine/bin/jasmine --version 2>&1 | head -1
  "$env_dir/bin/python" "$liftover_dir/bin/liftoverSV.py" --version
} > "$root/state/tool-versions.txt"
