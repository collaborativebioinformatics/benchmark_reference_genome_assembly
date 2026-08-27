#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
root=$BENCHMARK_ROOT
manifest="$ANALYSIS_EXPORT_ROOT/config/dnanexus-inputs.tsv"
dx_bin=/path/to/dx

while IFS=$'\t' read -r role arm build project remote id local; do
  [[ "$role" == role ]] && continue
  output="$root/$local"
  mkdir -p "${output%/*}"
  if [[ -s "$output" ]]; then
    echo "EXISTS $local"
  else
    echo "DOWNLOAD $project:$id $local"
    "$dx_bin" download "$project:$id" -o "$output"
  fi
done < "$manifest"

checksum_tmp="$root/state/input-sha256.tsv.${SLURM_JOB_ID:-local}.tmp"
while IFS=$'\t' read -r role arm build project remote id local; do
  [[ "$role" == role ]] && continue
  sha256sum "$root/$local"
done < "$manifest" \
  | sed "s#  $root/#  benchmarking/#" > "$checksum_tmp"
mv "$checksum_tmp" "$root/state/input-sha256.tsv"

"$root/.mamba-env/bin/python" "$ANALYSIS_EXPORT_ROOT/01-inputs/01-audit-input-vcfs.py" \
  --root "$root" \
  --manifest "$manifest" \
  --output-prefix "$root/state/vcf-input-audit"
