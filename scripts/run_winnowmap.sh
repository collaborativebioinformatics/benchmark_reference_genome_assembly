#!/usr/bin/env bash
# =============================================================================
# Group 4 — Winnowmap arm.  v2: real DNAnexus paths, tarball layout, .fastq reads.
#
# Launch (on your Mac):
#   dx run app-cloud_workstation --ssh --instance-type mem2_ssd1_v2_x32 \
#       -imax_session_length=24h
# Inside the workstation:
#   unset DX_WORKSPACE_ID; dx cd $DX_PROJECT_CONTEXT_ID:
#   dx download Group4_2026:/scripts/run_winnowmap.sh
#   bash run_winnowmap.sh hs1 2>&1 | tee run_hs1.log
#
# Optional coverage reduction (MUST match whatever Alisa uses for minimap2):
#   COVERAGE=10 bash run_winnowmap.sh hs1 2>&1 | tee run_hs1.log
#   COVERAGE=0 (default) uses all reads.
# =============================================================================

set -euo pipefail

PROJ="Group4_2026"
REL="/reference-assets/releases/2026-08-25"
BUILD="${1:-hs1}"
THREADS="${THREADS:-30}"
COVERAGE="${COVERAGE:-0}"          # 0 = use all reads; else target coverage for rasusa
SEED="${SEED:-42}"                 # fixed seed so every reference sees the same reads

if [[ -f "$HOME/cells.txt" && "$COVERAGE" != "0" ]]; then
  echo "ERROR: cells.txt and COVERAGE are both set. Use one subsampling method, not both." >&2
  exit 1
fi
OUTDIR="${PROJ}:/results/winnowmap/${BUILD}"

WORK="/home/dnanexus/work/${BUILD}"
mkdir -p "$WORK" && cd "$WORK"
log() { echo "[$(date -u +%H:%M:%S)] $*"; }
log "build=${BUILD} threads=${THREADS} coverage=${COVERAGE:-all} work=${WORK}"

# -----------------------------------------------------------------------------
# 0. Point dx at the PROJECT, not the temporary job workspace
# -----------------------------------------------------------------------------
unset DX_WORKSPACE_ID || true
dx select "$PROJ"
dx cd "${PROJ}:"
dx mkdir -p "$OUTDIR"
df -h /home/dnanexus | tail -1 | awk '{print "[disk] "$4" free"}'

# -----------------------------------------------------------------------------
# 1. Tools
# -----------------------------------------------------------------------------
export PATH="$HOME/env/bin:$PATH"
if ! command -v winnowmap >/dev/null 2>&1; then
  log "installing tools"
  if ! command -v micromamba >/dev/null 2>&1; then
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
    export PATH="$PWD/bin:$PATH"
  fi
  micromamba create -y -p "$HOME/env" -c conda-forge -c bioconda \
      winnowmap=2.03 meryl=1.4.2 sniffles samtools rasusa zstd
  export PATH="$HOME/env/bin:$PATH"
fi
command -v rasusa >/dev/null 2>&1 || micromamba install -y -p "$HOME/env" -c bioconda rasusa

# -----------------------------------------------------------------------------
# 2. Reference (a .tar.zst under genomes/<build>/reference/)
# -----------------------------------------------------------------------------
if [[ ! -s "${BUILD}.fa" ]]; then
  log "fetching reference"
  dx download -f "${PROJ}:${REL}/genomes/${BUILD}/reference/${BUILD}.reference.tar.zst"
  tar --use-compress-program=unzstd -xf "${BUILD}.reference.tar.zst"
  FA=$(find . -name "${BUILD}.fa" -o -name "${BUILD}.fa.gz" | head -1)
  [[ -n "$FA" ]] || { echo "ERROR: no ${BUILD}.fa in the reference tarball"; \
      tar --use-compress-program=unzstd -tf "${BUILD}.reference.tar.zst"; exit 1; }
  [[ "$FA" == *.gz ]] && { zcat "$FA" > "${BUILD}.fa"; } || { [[ "$FA" != "./${BUILD}.fa" ]] && ln -sf "$FA" "${BUILD}.fa"; }
fi
[[ -s "${BUILD}.fa.fai" ]] || samtools faidx "${BUILD}.fa"
echo "[check] contigs: $(wc -l < "${BUILD}.fa.fai")"
awk '$1=="chr6"||$1=="chr22"{printf "[check] %s %s bp\n",$1,$2}' "${BUILD}.fa.fai"

# -----------------------------------------------------------------------------
# 3. Winnowmap assets (meryl DB, .wmi indices, repetitive k-mer list) — all in
#    one tarball. The k-mer list is the file that actually matters at run time.
# -----------------------------------------------------------------------------
KMERS="indexes/${BUILD}/winnowmap/${BUILD}.repetitive-k15.txt"
WMI="indexes/${BUILD}/winnowmap/${BUILD}.map-pb.wmi"
if [[ ! -s "$KMERS" ]]; then
  log "fetching winnowmap assets (~3 GB)"
  dx download -f "${PROJ}:${REL}/indexes/${BUILD}/winnowmap/${BUILD}.winnowmap.tar.zst"
  tar --use-compress-program=unzstd -xf "${BUILD}.winnowmap.tar.zst"
fi
[[ -s "$KMERS" ]] || { echo "ERROR: repetitive k-mer list missing"; exit 1; }
echo "[check] repetitive k-mers: $(wc -l < "$KMERS")"

# -----------------------------------------------------------------------------
# 4. Reads — 39 uncompressed SMRT-cell FASTQs, ~160 GB total
# -----------------------------------------------------------------------------
mkdir -p reads && cd reads
mapfile -t READ_FILES < <(dx ls "${PROJ}:/Data/CCS_15kb_30x/" | grep -E '\.fastq$')
# Fixed subset: if ~/cells.txt exists, use exactly those files (see config/subset_cells_15x.txt)
if [[ -f "$HOME/cells.txt" ]]; then
  mapfile -t READ_FILES < "$HOME/cells.txt"
  log "using fixed cell subset: ${#READ_FILES[@]} files from ~/cells.txt"
fi
[[ ${#READ_FILES[@]} -gt 0 ]] || { echo "ERROR: no FASTQ found"; exit 1; }
log "downloading ${#READ_FILES[@]} FASTQ files"
for f in "${READ_FILES[@]}"; do
  [[ -s "$f" ]] || dx download -f "${PROJ}:/Data/CCS_15kb_30x/${f}"
done
cd "$WORK"
df -h /home/dnanexus | tail -1 | awk '{print "[disk] "$4" free after reads"}'

if [[ "$COVERAGE" != "0" ]]; then
  SUB="reads.${COVERAGE}x.seed${SEED}.fastq.gz"
  if [[ ! -s "$SUB" ]]; then
    log "subsampling to ${COVERAGE}x (seed ${SEED})"
    cat reads/*.fastq | rasusa reads -c "${COVERAGE}" -g 3.1gb -s "${SEED}" -o "$SUB"
    # keep the exact subset so every reference and both aligners use identical reads
    dx upload "$SUB" --path "${PROJ}:/Data/subsampled/" || true
  fi
  QUERY=("$SUB")
else
  QUERY=(reads/*.fastq)
fi

# -----------------------------------------------------------------------------
# 5. Align — WHOLE genome. Do not subset the reference: reads from elsewhere
#    would pile onto chr6/chr22 and manufacture false SVs. Restrict later.
# -----------------------------------------------------------------------------
BAM="HG002.${BUILD}.winnowmap.bam"
if [[ ! -s "${BAM}.bai" ]]; then
  log "aligning with winnowmap"
  winnowmap -W "$KMERS" -ax map-pb -t "$THREADS" -Y \
      -R "@RG\tID:${BUILD}_wm\tSM:HG002\tPL:PACBIO\tLB:CCS_15kb" \
      "${BUILD}.fa" "${QUERY[@]}" \
    | samtools sort -@ 4 -m 3G -T sort_tmp -o "$BAM" -
  samtools index -@ 8 "$BAM"
fi
samtools flagstat -@ 8 "$BAM" > "HG002.${BUILD}.winnowmap.flagstat"
cat "HG002.${BUILD}.winnowmap.flagstat"

# upload the expensive artefact immediately, before variant calling
dx upload --path "${OUTDIR}/" "$BAM" "${BAM}.bai" "HG002.${BUILD}.winnowmap.flagstat"
log "BAM uploaded"

# -----------------------------------------------------------------------------
# 6. Call. --output-rnames is how variants are matched across references.
#    It cannot be added retrospectively.
# -----------------------------------------------------------------------------
VCF="HG002.${BUILD}.winnowmap.sniffles.vcf.gz"
log "calling with sniffles"
sniffles --input "$BAM" --reference "${BUILD}.fa" \
         --vcf "$VCF" --snf "HG002.${BUILD}.winnowmap.snf" \
         --minsvlen 50 --output-rnames --threads "$THREADS"

# -----------------------------------------------------------------------------
# 7. Versions — these fill the Methods placeholders
# -----------------------------------------------------------------------------
{
  echo "build: ${BUILD}"
  echo "date_utc: $(date -u +%FT%TZ)"
  echo "coverage: ${COVERAGE} (0 = all reads); seed: ${SEED}"
  echo "n_read_files: ${#READ_FILES[@]}"
  winnowmap --version 2>&1 | sed 's/^/winnowmap: /'
  samtools --version 2>&1 | head -1 | sed 's/^/samtools: /'
  sniffles --version 2>&1 | head -1 | sed 's/^/sniffles: /'
} > "versions.${BUILD}.txt"
cat "versions.${BUILD}.txt"

# -----------------------------------------------------------------------------
# 8. Upload and verify
# -----------------------------------------------------------------------------
dx upload --path "${OUTDIR}/" "$VCF" "HG002.${BUILD}.winnowmap.sniffles.snf" "versions.${BUILD}.txt"
echo "== ${OUTDIR} =="
dx ls -l "$OUTDIR"
log "done: ${BUILD}"
echo "Anything showing 'closing' is still uploading — do not terminate yet."
