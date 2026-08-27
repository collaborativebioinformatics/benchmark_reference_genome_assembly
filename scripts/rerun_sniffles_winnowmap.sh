#!/usr/bin/env bash

## Follow the run_winnowmap.sh
## bash rerun_sniffles.sh hg17 2>&1 | tee rerun_sniffles_hg17.log
set -euo pipefail

PROJ="Group4_2026"
BUILD="${1:-hs1}"
THREADS="${THREADS:-32}"
COVERAGE="${COVERAGE:-0}"
SEED="${SEED:-42}"
OUTDIR="/results/winnowmap/${BUILD}"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

unset DX_WORKSPACE_ID || true
dx select "$PROJ"
dx cd "${PROJ}:"

# -----------------------------------------------------------------------------
# 1. Tools — rebuild with Python 3.10
# -----------------------------------------------------------------------------

if ! command -v micromamba >/dev/null 2>&1; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
  export PATH="$PWD/bin:$PATH"
fi

micromamba remove -y -p "$HOME/env" --all || true

micromamba create -y -p "$HOME/env" -c conda-forge -c bioconda \
    python=3.10 \
    winnowmap=2.03 meryl=1.4.2 sniffles samtools rasusa zstd

export PATH="$HOME/env/bin:$PATH"

python --version
sniffles --version

# -----------------------------------------------------------------------------
# 2. Find existing alignment and reference
# -----------------------------------------------------------------------------

BAM=$(find /home/dnanexus -type f -name "HG002.${BUILD}.winnowmap.bam" 2>/dev/null | head -1)
REF=$(find /home/dnanexus -type f -name "${BUILD}.fa" 2>/dev/null | head -1)

[[ -n "$BAM" ]] || { echo "ERROR: HG002.${BUILD}.winnowmap.bam not found"; exit 1; }
[[ -n "$REF" ]] || { echo "ERROR: ${BUILD}.fa not found"; exit 1; }

WORK=$(dirname "$BAM")
cd "$WORK"

BAM=$(basename "$BAM")
REF=$(realpath "$REF")

echo "[check] BAM: $WORK/$BAM"
echo "[check] REF: $REF"

# -----------------------------------------------------------------------------
# 6. Call
# -----------------------------------------------------------------------------

VCF="HG002.${BUILD}.winnowmap.sniffles.vcf.gz"

log "calling with sniffles"

sniffles --input "$BAM" --reference "$REF" \
         --vcf "$VCF" --snf "HG002.${BUILD}.winnowmap.sniffles.snf" \
         --minsvlen 50 --output-rnames --threads "$THREADS" --allow-overwrite

# -----------------------------------------------------------------------------
# 7. Versions
# -----------------------------------------------------------------------------

{
  echo "build: ${BUILD}"
  echo "date_utc: $(date -u +%FT%TZ)"
  echo "coverage: ${COVERAGE} (0 = all reads); seed: ${SEED}"
  echo "python: $(python --version 2>&1)"
  winnowmap --version 2>&1 | sed 's/^/winnowmap: /'
  samtools --version 2>&1 | head -1 | sed 's/^/samtools: /'
  sniffles --version 2>&1 | head -1 | sed 's/^/sniffles: /'
} > "versions.${BUILD}.txt"

cat "versions.${BUILD}.txt"

# -----------------------------------------------------------------------------
# 8. Upload and verify
# -----------------------------------------------------------------------------

dx mkdir -p "$OUTDIR"

dx upload --path "${OUTDIR}/" \
    "$VCF" \
    "HG002.${BUILD}.winnowmap.sniffles.snf" \
    "versions.${BUILD}.txt"

echo "== ${OUTDIR} =="
dx ls -l "$OUTDIR"

log "done: ${BUILD}"
echo "Anything showing 'closing' is still uploading — do not terminate yet."
