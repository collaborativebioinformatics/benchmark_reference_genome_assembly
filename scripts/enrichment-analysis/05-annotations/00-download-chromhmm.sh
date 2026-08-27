#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
root=$BENCHMARK_ROOT
url=https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E116_15_coreMarks_hg38lift_mnemonics.bed.gz
output="$root/inputs/annotations/chromhmm/grch38/E116_15_coreMarks_hg38lift_mnemonics.bed.gz"
mkdir -p "${output%/*}"
if [[ ! -s "$output" ]]; then
  curl --fail --location --retry 5 --output "$output.part" "$url"
  mv "$output.part" "$output"
fi
gzip -t "$output"
sha256sum "$output" > "$output.sha256"
