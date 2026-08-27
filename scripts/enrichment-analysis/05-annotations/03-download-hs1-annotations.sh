#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
root=$BENCHMARK_ROOT

download() {
  local url=$1
  local output=$2
  mkdir -p "$(dirname "$output")"
  curl -fL --retry 5 --retry-delay 5 -o "$output.part" "$url"
  mv "$output.part" "$output"
}

download \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3 \
  "$root/inputs/annotations/gencode/hs1/chm13.draft_v2.0.gene_annotation.gff3"
download \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed \
  "$root/inputs/annotations/repeatmasker/hs1/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
download \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed \
  "$root/inputs/annotations/segdup/hs1/chm13v2.0_SD.bed"
download \
  http://hgdownload.soe.ucsc.edu/gbdb/hs1/hoffmanMappability/k100.Unique.Mappability.bb \
  "$root/inputs/annotations/mappability/hs1/k100.Unique.Mappability.bb"
download \
  http://hgdownload.soe.ucsc.edu/gbdb/hs1/hoffmanMappability/k100.Umap.MultiTrackMappability.bw \
  "$root/inputs/annotations/mappability/hs1/k100.Umap.MultiTrackMappability.bw"

find "$root/inputs/annotations" -type f \
  \( -path '*/gencode/hs1/*' -o -path '*/repeatmasker/hs1/*' -o \
     -path '*/segdup/hs1/*' -o -path '*/mappability/hs1/*' \) \
  -print0 | sort -z | xargs -0 sha256sum > "$root/state/hs1-annotation-inputs.sha256"
