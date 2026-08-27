#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
root=$BENCHMARK_ROOT
tool_dir="$root/tools/ucsc"
mkdir -p "$tool_dir"

for tool in bigBedToBed bigWigAverageOverBed; do
  curl -fL --retry 5 --retry-delay 5 \
    -o "$tool_dir/$tool.part" \
    "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$tool"
  chmod 0755 "$tool_dir/$tool.part"
  mv "$tool_dir/$tool.part" "$tool_dir/$tool"
done

sha256sum "$tool_dir/bigBedToBed" "$tool_dir/bigWigAverageOverBed" \
  > "$root/state/ucsc-annotation-tools.sha256"
