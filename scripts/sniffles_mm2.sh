#!/usr/bin/env bash

set -euo pipefail

PROJ="Group4_2026"
REL="/reference-assets/releases/2026-08-25"
BUILD="${BUILD:?set BUILD, e.g. BUILD=hg15}"
THREADS="${THREADS:-30}"
PRESET="map-hifi"

RUN_DIR="${PROJ}:/results/minimap2/${BUILD}"
PREFIX="HG002.${BUILD}.minimap2.${PRESET}"
BAM_REMOTE="${RUN_DIR}/${PREFIX}.bam"
BAI_REMOTE="${RUN_DIR}/${PREFIX}.bam.bai"
CHECKSUM_REMOTE="${RUN_DIR}/${PREFIX}.sha256"
OUTDIR="${PROJ}:/results/sniffles/minimap2/${BUILD}"

# -----------------------------------------------------------------------------
# 0. Point dx at the PROJECT, not the temporary job workspace
# -----------------------------------------------------------------------------
unset DX_WORKSPACE_ID
dx select "$PROJ"
dx cd /

# -----------------------------------------------------------------------------
# 1. Tools
# -----------------------------------------------------------------------------
export PATH="$HOME/env/bin:$PATH"
if ! command -v sniffles >/dev/null 2>&1; then
  echo "installing tools"
  if ! command -v micromamba >/dev/null 2>&1; then
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
    export PATH="$PWD/bin:$PATH"
  fi
  micromamba create -y -p "$HOME/env" -c conda-forge -c bioconda \
  python=3.11 sniffles samtools rasusa zstd
  export PATH="$HOME/env/bin:$PATH"
fi

# -----------------------------------------------------------------------------
# 2. Reference (a .tar.zst under genomes/<build>/reference/)
# -----------------------------------------------------------------------------
dx download "${PROJ}:${REL}/genomes/${BUILD}/reference/${BUILD}.reference.tar.zst" -f -o ref.tar.zst
unzstd -f ref.tar.zst
tar -xf ref.tar
REF="references/${BUILD}/${BUILD}.fa"

# -----------------------------------------------------------------------------
# 3. Obtain BAM files
# -----------------------------------------------------------------------------

BAM="${PREFIX}.bam"
dx download "$BAM_REMOTE" -o "$BAM"
dx download "$BAI_REMOTE" -o "${BAM}.bai"



# -----------------------------------------------------------------------------
# 4. Run Sniffles
# -----------------------------------------------------------------------------
VCF="${PREFIX}.sniffles.vcf"
SNF="${PREFIX}.sniffles.snf"
sniffles --input "$BAM" --reference "$REF" \
--vcf "$VCF" --snf "$SNF" \
--minsvlen 50 --output-rnames --threads "$THREADS"

{
  echo "sniffles: $(sniffles --version)"
  echo "source_alignment: ${BAM_REMOTE} (map-hifi)"
} > "versions.${BUILD}.mm2.txt"

dx upload "$VCF" "$SNF" "versions.${BUILD}.mm2.txt" --destination "${OUTDIR}/" --wait
