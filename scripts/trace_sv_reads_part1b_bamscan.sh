#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# PART Ib: recover per-read alignment coordinates in every build
#
# Part I can only see reads that Sniffles put into an RNAMES field. A read that
# aligned but produced no SV call, a read that failed to align, and a read that
# is absent from the BAM all collapse into NO_SV_CALL. This step separates them.
#
# For each build:
#   1. scan the BAM for the anchor supporting reads -> per-alignment TSV
#   2. scan the reference FASTA for N-runs        -> gap BED + contig sizes
#
# Both are streamed; nothing large is kept on disk between builds.
#
# Usage:
#   ./trace_sv_reads_part1b_bamscan.sh <dx-path-of-part-I-output-dir>
# e.g.
#   ./trace_sv_reads_part1b_bamscan.sh results/read_trace/minimap2/20260828-101500
# =============================================================================

TRACE_DIR="${1:?usage: $0 <dx-path-of-part-I-output-dir>}"

PROJECT="project-JB6zPY00Z7Q7KBpj1P5yKYv4"
BUILDS=(hg15 hg16 hg17 hg18 hg19 hg38 hs1)
BAM_ROOT="results/minimap2"          # expects $BAM_ROOT/<build>/*.bam
REF_ROOT="reference-assets/releases/2026-08-25"
OUT="results/read_trace_bam/minimap2/$(date -u +%Y%m%d-%H%M%S)"
THREADS="${THREADS:-8}"

unset DX_WORKSPACE_ID || true
dx select "$PROJECT"

# samtools -N (read-name filter file) requires >= 1.12
sv=$(samtools --version | head -1 | awk '{print $2}')
awk -v v="$sv" 'BEGIN{split(v,a,".");
    if (a[1]<1 || (a[1]==1 && a[2]<12)) {print "FATAL: samtools >= 1.12 required, found " v > "/dev/stderr"; exit 1}}'

# References are bundled as <build>.reference.tar.zst, not flat <build>.fa.gz
command -v zstd >/dev/null 2>&1 || { apt-get update -qq && apt-get install -y -qq zstd; }

WORK="$HOME/bamscan_$(date -u +%Y%m%d-%H%M%S)"
mkdir -p "$WORK"; cd "$WORK"
echo "Working directory: $WORK"

# -----------------------------------------------------------------------------
# Anchor read names, taken from Part I so the two stages cannot drift apart
# -----------------------------------------------------------------------------
dx download -f "$PROJECT:$TRACE_DIR/read_trace_best.tsv" -o read_trace_best.tsv
awk -F'\t' 'NR>1 && $7!="NA" {print $7}' read_trace_best.tsv | sort -u > anchor_reads.txt
n_reads=$(wc -l < anchor_reads.txt)
[[ "$n_reads" -gt 0 ]] || { echo "FATAL: no read names in read_trace_best.tsv" >&2; exit 1; }
echo "Anchor supporting reads: $n_reads"

# -----------------------------------------------------------------------------
# SAM -> compact per-alignment TSV
# -----------------------------------------------------------------------------
cat > parse_sam.py <<'PY'
import re, sys

BUILD = sys.argv[1]
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
out = sys.stdout
out.write("read_id\tbuild\tchrom\tpos\tend\tmapq\tflag\taln_type\tclip5\tclip3\tref_span\tnm\n")

for line in sys.stdin:
    if line.startswith("@"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 11:
        continue
    qname, flag, rname, pos, mapq, cigar = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5]

    if flag & 0x100:
        aln_type = "secondary"
    elif flag & 0x800:
        aln_type = "supplementary"
    elif flag & 0x4:
        aln_type = "unmapped"
    else:
        aln_type = "primary"

    if aln_type == "unmapped" or cigar == "*":
        out.write(f"{qname}\t{BUILD}\t*\t0\t0\t{mapq}\t{flag}\t{aln_type}\t0\t0\t0\tNA\n")
        continue

    ops = CIG.findall(cigar)
    span = sum(int(n) for n, op in ops if op in "MDN=X")
    clip5 = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
    clip3 = int(ops[-1][0]) if len(ops) > 1 and ops[-1][1] in "SH" else 0

    nm = "NA"
    for t in f[11:]:
        if t.startswith("NM:i:"):
            nm = t[5:]
            break

    out.write(f"{qname}\t{BUILD}\t{rname}\t{pos}\t{pos + span - 1}\t{mapq}\t{flag}\t"
              f"{aln_type}\t{clip5}\t{clip3}\t{span}\t{nm}\n")
PY

# -----------------------------------------------------------------------------
# FASTA -> N-run BED (0-based half-open) + contig sizes
# Fast path: most lines contain no N, so the membership test short-circuits.
# -----------------------------------------------------------------------------
cat > nruns.py <<'PY'
import re, sys

RUN = re.compile(r"[Nn]+")
bed = open(sys.argv[1], "w")
sizes = open(sys.argv[2], "w")

chrom = None
offset = 0          # 0-based offset of the current line's first base
open_start = None   # start of an N-run continuing across lines

def close_run(end):
    global open_start
    if open_start is not None:
        bed.write(f"{chrom}\t{open_start}\t{end}\n")
        open_start = None

for line in sys.stdin:
    if line[0] == ">":
        close_run(offset)
        if chrom is not None:
            sizes.write(f"{chrom}\t{offset}\n")
        chrom = line[1:].split()[0]
        offset = 0
        continue
    s = line.rstrip()
    if not s:
        continue
    if "N" not in s and "n" not in s:
        close_run(offset)
        offset += len(s)
        continue
    prev_end = None
    for m in RUN.finditer(s):
        a, b = offset + m.start(), offset + m.end()
        if open_start is not None and m.start() == 0:
            pass                      # run continues from the previous line
        else:
            close_run(prev_end if prev_end is not None else a)
            open_start = a
        prev_end = b
        if m.end() != len(s):         # run ends inside this line
            close_run(b)
    offset += len(s)

close_run(offset)
if chrom is not None:
    sizes.write(f"{chrom}\t{offset}\n")
bed.close(); sizes.close()
PY

# -----------------------------------------------------------------------------
# Per build: BAM scan, then FASTA scan
# -----------------------------------------------------------------------------
for b in "${BUILDS[@]}"; do
    echo
    echo "=== $b ==="

    bam=$(dx find data --path "$PROJECT:$BAM_ROOT/$b" \
              --name "*.bam" --brief 2>&1 | head -1 || true)
    if [[ -z "$bam" || "$bam" == *Error* ]]; then
        echo "WARNING: no BAM found under $BAM_ROOT/$b, skipping ($bam)" >&2
    else
        echo "BAM: $bam"
        dx download -f "$bam" -o "$b.bam"
        samtools view -@ "$THREADS" -N anchor_reads.txt "$b.bam" \
            | python3 parse_sam.py "$b" \
            | gzip -c > "aln.$b.tsv.gz"
        rm -f "$b.bam"
        echo "  alignments: $(( $(zcat "aln.$b.tsv.gz" | wc -l) - 1 ))"
    fi

    fa_tar=$(dx find data --path "$PROJECT:$REF_ROOT" \
                 --name "$b.reference.tar.zst" --brief 2>&1 | head -1 || true)
    if [[ -z "$fa_tar" || "$fa_tar" == *Error* ]]; then
        echo "WARNING: no reference tarball found for $b under $REF_ROOT, skipping gap scan ($fa_tar)" >&2
    else
        dx download -f "$fa_tar" -o "$b.reference.tar.zst"
        tar --zstd -xOf "$b.reference.tar.zst" "references/$b/$b.fa" \
            | python3 nruns.py "gaps.$b.bed" "sizes.$b.tsv"
        rm -f "$b.reference.tar.zst"
        echo "  N-runs: $(wc -l < "gaps.$b.bed")  gap bp: $(awk '{s+=$3-$2} END{print s+0}' "gaps.$b.bed")"
    fi
done

# -----------------------------------------------------------------------------
# Provenance and upload
# -----------------------------------------------------------------------------
{
    echo "date_utc: $(date -u +%FT%TZ)"
    echo "project: $PROJECT"
    echo "part1_trace_dir: $TRACE_DIR"
    echo "bam_root: $BAM_ROOT"
    echo "ref_root: $REF_ROOT"
    echo "builds: ${BUILDS[*]}"
    echo "n_anchor_reads: $n_reads"
    echo "samtools: $sv"
} > provenance_bamscan.txt

# swiss-army-knife jobs only ever get VIEW project access (see app spec), so a
# script-issued `dx upload`/`dx mkdir` 401s no matter how privileged the user
# who launched the job is. Instead, drop results in ~/out/out/ and let the
# app's own output-collection step (which runs with its own authority, not
# this job's token) upload them to whatever --destination `dx run` was given.
mkdir -p "$HOME/out/out"
for f in aln.*.tsv.gz gaps.*.bed sizes.*.tsv anchor_reads.txt provenance_bamscan.txt; do
    [[ -s "$f" ]] || continue
    echo "Staging $f for output"
    cp "$f" "$HOME/out/out/"
done

echo
echo "Done. Results will land at whatever --destination this job was run with."
