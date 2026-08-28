#!/usr/bin/env bash
# HG002 SV benchmark across the hg15..hs1 reference ladder — read-name-tracing arm.
# Sniffles output filename conventions differ by aligner (confirmed against
# the real project layout):
#   minimap2:  results/sniffles/minimap2/<build>/HG002.<build>.minimap2.map-hifi.sniffles.vcf
#   winnowmap: results/sniffles/winnowmap/<build>/HG002.<build>.winnowmap.sniffles.vcf.gz

# =============================================================================
# Launch (on your Mac):
#   dx run app-cloud_workstation --ssh --instance-type mem2_ssd1_v2_x32 -imax_session_length=24h
# Inside the workstation:
#   unset DX_WORKSPACE_ID; dx cd $DX_PROJECT_CONTEXT_ID:
#   dx download Group4_2026:/scripts/run_benchmark_rname.sh
#   ANCHOR=hg38 ALIGNERS=minimap2 bash run_benchmark_rname.sh   # env vars, not script args; ANCHOR defaults to hs1
# =============================================================================
set -euo pipefail
 
PROJ="Group4_2026"
ANCHOR="${ANCHOR:-hs1}" # must be a build with a native GIAB truth set
 
BUILDS=(hg15 hg16 hg17 hg18 hg19 hg38 hs1)
IFS=' ' read -r -a ALIGNERS <<< "${ALIGNERS:-minimap2 winnowmap}"
case "$ANCHOR" in
  hs1)
    TRUTH="Data/truth_set/HG002_CHM13v2.0_v5.0q_stvar"
    ;;
  hg19)
    TRUTH="Data/truth_set/HG002_GRCh37_v5.0q_stvar"
    ;;
  hg38)
    TRUTH="Data/truth_set/HG002_GRCh38_v5.0q_stvar"
    ;;
  *)
    echo "FATAL: ANCHOR=$ANCHOR has no known native GIAB truth set (must be hg19, hg38, or hs1)." >&2
    exit 1
    ;;
esac
 
REF="reference-assets/releases/2026-08-25/genomes/${ANCHOR}/reference/${ANCHOR}.reference.tar.zst"
QUERIES_ROOT=results/sniffles
OUT=results/benchmark_rname/benchmark_$ANCHOR
OVERLAP_MIN="${OVERLAP_MIN:-0.5}" # fraction of an anchor TP call's reads that must reappear in a build's call to inherit its label
REQUIRE_TYPE_MATCH="${REQUIRE_TYPE_MATCH:-true}"  # only match against anchor TP calls of the same SVTYPE as the build's call
 
sniffles_src_path() {
  local aligner="$1" b="$2"
  case "$aligner" in
    minimap2)
      echo "$QUERIES_ROOT/minimap2/$b/HG002.$b.minimap2.map-hifi.sniffles.vcf"
      ;;
    winnowmap)
      echo "$QUERIES_ROOT/winnowmap/$b/HG002.$b.winnowmap.sniffles.vcf.gz"
      ;;
    *)
      echo "FATAL: no known Sniffles filename convention for aligner '$aligner'" >&2
      exit 1
      ;;
  esac
}
 
cd ~ && mkdir -p results/benchmark_rname && cd results/benchmark_rname
 
# -----------------------------------------------------------------------------
# 0. POINT TO PROJECT
# -----------------------------------------------------------------------------
unset DX_WORKSPACE_ID || true
dx select "$PROJ"
dx cd "${PROJ}:"
dx mkdir -p "$OUT"
df -h /home/dnanexus | tail -1 | awk '{print "[disk] "$4" free"}'
# -----------------------------------------------------------------------------
# 1. INSTALLATION OF TOOLS
# -----------------------------------------------------------------------------
MAMBA_ENV="${MAMBA_ENV:-benchmark_rname}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

if ! command -v micromamba >/dev/null 2>&1; then
  echo "micromamba not found — installing a static binary to \$HOME/bin ..." >&2
  mkdir -p "$HOME/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "$HOME" bin/micromamba 2>&1
  export PATH="$HOME/bin:$PATH"
fi

eval "$(micromamba shell hook --shell bash)"

if micromamba env list | grep -q "^${MAMBA_ENV}[[:space:]]"; then
  echo "Reusing existing micromamba environment '$MAMBA_ENV'." >&2
else
  echo "Creating micromamba environment '$MAMBA_ENV' (bcftools, samtools, htslib, truvari, python, matplotlib)..." >&2
  micromamba create -n "$MAMBA_ENV" -y -c bioconda -c conda-forge \
    bcftools samtools htslib truvari "python=3.10" matplotlib numpy
fi

micromamba activate "$MAMBA_ENV"

# Checks for installation of tool
for tool in truvari bcftools samtools tabix python3; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "FATAL: '$tool' not found after activating micromamba env '$MAMBA_ENV'." >&2
    exit 1
  }
done
 
{ echo "truvari_version: $(truvari version)"; \
  echo "bcftools_version: $(bcftools --version | head -1)"; \
  echo "samtools_version: $(samtools --version | head -1)"; \
  echo "tabix_version: $(tabix --version 2>&1 | head -1)"; \
  echo "python3_version: $(python3 --version 2>&1)"; \
  echo "mamba_env: $MAMBA_ENV"; } > packages_version.txt
cat packages_version.txt
# -----------------------------------------------------------------------------
# 2. REF + TRUTH for ANCHOR BUILD
# -----------------------------------------------------------------------------
dx download -f "$PROJ:$TRUTH.vcf.gz"{,.tbi} "$PROJ:$TRUTH.benchmark.bed"
mv "$(basename "$TRUTH").vcf.gz" truth.vcf.gz
mv "$(basename "$TRUTH").vcf.gz.tbi" truth.vcf.gz.tbi
sort -k1,1 -k2,2n "$(basename "$TRUTH").benchmark.bed" > truth.bed

dx download -f "$PROJ:$REF" -o anchor.tar.zst
mkdir -p ref
tar --use-compress-program=unzstd -xf anchor.tar.zst -C ref
ANCHOR_FA=$(find ref -maxdepth 4 -name '*.fa' | head -1)
[[ -n "$ANCHOR_FA" ]] || { echo "FATAL: no FASTA in anchor reference tarball" >&2; exit 1; }
[[ -s "${ANCHOR_FA}.fai" ]] || samtools faidx "$ANCHOR_FA"
 
{ echo "date: $(date -u +%FT%TZ)"; echo "anchor: $ANCHOR"; echo "aligners: ${ALIGNERS[*]}"; \
  echo "overlap_min: $OVERLAP_MIN"; echo "require_type_match: $REQUIRE_TYPE_MATCH"; \
  truvari version; bcftools --version | head -1; } > provenance.txt
 
cat > match_by_rnames.py <<'PY'
import sys, gzip, json
from collections import defaultdict
 
def op(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
 
def norm_type(t):
    """truvari bench ran with dup_to_ins=true (see params.json), so the
    anchor's own TP calls may already have DUP folded into INS by truvari
    before this script ever sees them. Treat DUP == INS for the type check;
    keep everything else strict."""
    return 'INS' if t == 'DUP' else t
 
def rname_sets(vcf_path):
    """SV id -> {'reads': frozenset(supporting read names), 'svtype': str},
    for every record with RNAMES."""
    out = {}
    with op(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            info = dict(
                kv.split('=', 1) if '=' in kv else (kv, True)
                for kv in cols[7].split(';')
            )
            rn = info.get('RNAMES')
            if not rn:
                continue
            call_id = cols[2] if cols[2] != '.' else f'{cols[0]}:{cols[1]}'
            out[call_id] = {
                'reads': frozenset(rn.split(',')),
                'svtype': info.get('SVTYPE', 'NA'),
            }
    return out
 
def best_overlap(query_reads, candidates):
    """Which candidate anchor TP call shares the largest fraction of ITS OWN
    reads with this query call — i.e. how much of the anchor's evidence
    does this build's call reproduce. `candidates` is a dict of
    call_id -> read set, already restricted to the matching SV type
    (or the full anchor set, if type-matching is off)."""
    best_id, best_frac = None, 0.0
    for aid, aset in candidates.items():
        if not aset:
            continue
        shared = len(query_reads & aset)
        frac = shared / len(aset)
        if frac > best_frac:
            best_id, best_frac = aid, frac
    return best_id, best_frac
 
if __name__ == '__main__':
    anchor_tp_vcf, build_vcf, min_overlap, out_json, require_type_match = sys.argv[1:6]
    require_type_match = require_type_match.strip().lower() in ('1', 'true', 'yes')
 
    anchor_tp = rname_sets(anchor_tp_vcf)
    build_calls = rname_sets(build_vcf)
 
    # bucket anchor TP calls by normalized SV type so a build call is only
    # ever compared against anchor calls of the same type
    anchor_by_type = defaultdict(dict)
    for aid, a in anchor_tp.items():
        anchor_by_type[norm_type(a['svtype'])][aid] = a['reads']
    anchor_all = {aid: a['reads'] for aid, a in anchor_tp.items()}
 
    results, matched, type_blocked = [], 0, 0
    for call_id, c in build_calls.items():
        svtype = c['svtype']
        if require_type_match:
            candidates = anchor_by_type.get(norm_type(svtype), {})
        else:
            candidates = anchor_all
        aid, frac = best_overlap(c['reads'], candidates)
 
        # for visibility: would this call have inherited a TP label under the
        # old, type-blind matching? only meaningful when type-matching is on
        # and this call didn't clear the (type-restricted) threshold
        if require_type_match and frac < float(min_overlap):
            alt_aid, alt_frac = best_overlap(c['reads'], anchor_all)
            if alt_frac >= float(min_overlap) and norm_type(anchor_tp[alt_aid]['svtype']) != norm_type(svtype):
                type_blocked += 1
 
        label = 'inherited_TP' if frac >= float(min_overlap) else 'unmatched'
        if label == 'inherited_TP':
            matched += 1
        results.append({
            'call': call_id,
            'call_svtype': svtype,
            'best_anchor_match': aid,
            'anchor_svtype': anchor_tp[aid]['svtype'] if aid else None,
            'read_overlap_frac': round(frac, 3),
            'label': label,
        })
 
    n = len(build_calls)
    json.dump({
        'n_calls': n,
        'n_inherited_tp': matched,
        'inherited_tp_rate': round(matched / n, 3) if n else 0,
        'require_type_match': require_type_match,
        'n_excluded_by_type_check': type_blocked,
        'calls': results,
    }, open(out_json, 'w'), indent=2)
PY
# -----------------------------------------------------------------------------
# 3. OCCURS ONCE PER ALIGNER - truvari bench on anchor, read name traced ladder
#     for the other 6 builds.
# -----------------------------------------------------------------------------
for ALIGNER in "${ALIGNERS[@]}"; do
  echo "==== aligner: $ALIGNER ====" >&2
 
  # 3a. every build's NATIVE-coordinate Sniffles VCF for this aligner
  for b in "${BUILDS[@]}"; do
    rel="$(sniffles_src_path "$ALIGNER" "$b")"
    src="$PROJ:$rel"
    dx describe "$src" >/dev/null 2>&1 || { echo "skip $ALIGNER/$b: no VCF at $rel" >&2; continue; }
    raw="$ALIGNER.$b.raw.vcf"
    [[ "$rel" == *.gz ]] && raw="$raw.gz"
    dx download -f "$src" -o "$raw"
    bcftools view -f PASS "$raw" \
      | bcftools filter -e 'INFO/SVTYPE=="BND" || INFO/SVTYPE=="INV"' \
      | bcftools sort -Oz -o "$ALIGNER.$b.norm.vcf.gz"
    tabix -f -p vcf "$ALIGNER.$b.norm.vcf.gz"
  done
 
  if [[ ! -s "$ALIGNER.$ANCHOR.norm.vcf.gz" ]]; then
    echo "SKIP aligner=$ALIGNER entirely: no $ANCHOR VCF found for it — nothing to bench or match against" >&2
    continue
  fi
 
  # 3b. truvari on the anchor only
  rm -rf "bench_${ALIGNER}_${ANCHOR}"
  truvari bench -b truth.vcf.gz -c "$ALIGNER.$ANCHOR.norm.vcf.gz" -f "$ANCHOR_FA" --includebed truth.bed \
    -r 500 -p 0.7 -P 0.7 -s 50 -S 30 --sizemax 1000000 \
    --passonly --dup-to-ins -o "bench_${ALIGNER}_${ANCHOR}"
  truvari refine --reference "$ANCHOR_FA" --regions truth.bed \
    --use-original-vcfs "bench_${ALIGNER}_${ANCHOR}" >/dev/null 2>&1 || true
 
  ANCHOR_TP="bench_${ALIGNER}_${ANCHOR}/tp-comp.vcf.gz"
  [[ -s "$ANCHOR_TP" ]] || ANCHOR_TP="bench_${ALIGNER}_${ANCHOR}/tp-call.vcf.gz"
  [[ -s "$ANCHOR_TP" ]] || { echo "FATAL: no TP vcf found under bench_${ALIGNER}_${ANCHOR}/" >&2; exit 1; }
 
  # 3c. read-name matcher: extend this aligner's anchor TP labels to every other build aligned with the SAME aligner
  printf 'build\tn_calls\tn_inherited_tp\tinherited_tp_rate\n' > "summary_rname.${ALIGNER}.tsv"
  for b in "${BUILDS[@]}"; do
    [[ "$b" == "$ANCHOR" ]] && continue
    [[ -s "$ALIGNER.$b.norm.vcf.gz" ]] || continue
    out_json="$ALIGNER.$b.rname_match.json"
    python3 match_by_rnames.py "$ANCHOR_TP" "$ALIGNER.$b.norm.vcf.gz" "$OVERLAP_MIN" "$out_json" "$REQUIRE_TYPE_MATCH"
    read -r n tp rate < <(python3 -c "
import json; d = json.load(open('$out_json'))
print(d['n_calls'], d['n_inherited_tp'], d['inherited_tp_rate'])")
    printf '%s\t%s\t%s\t%s\n' "$b" "$n" "$tp" "$rate" >> "summary_rname.${ALIGNER}.tsv"
  done
 
  echo "-- $ALIGNER --"
  column -t "summary_rname.${ALIGNER}.tsv"
done
# -----------------------------------------------------------------------------
# 4. VISUALS — one set of figures covering every aligner that actually ran
# -----------------------------------------------------------------------------
cat > make_plots.py <<'PY'
import glob, json, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
 
ANCHOR = os.environ['ANCHOR']
BUILD_ORDER = ['hg15', 'hg16', 'hg17', 'hg18', 'hg19', 'hg38', 'hs1']
COLORS = {'minimap2': '#2E86AB', 'winnowmap': '#C1440E'}
 
aligners = sorted({f.split('.')[1] for f in glob.glob('summary_rname.*.tsv')})
if not aligners:
    print('No summary_rname.*.tsv files found — nothing to plot.')
    raise SystemExit(0)
 
def color(a, i):
    return COLORS.get(a, plt.cm.tab10(i))
 
# ---- Figure 1: inherited_tp_rate ladder, one line per aligner ----
fig, ax = plt.subplots(figsize=(8, 5))
for i, aligner in enumerate(aligners):
    rates = {}
    with open(f'summary_rname.{aligner}.tsv') as f:
        next(f)
        for line in f:
            build, n_calls, n_tp, rate = line.rstrip('\n').split('\t')
            rates[build] = float(rate)
    xs = [b for b in BUILD_ORDER if b in rates]
    ys = [rates[b] for b in xs]
    ax.plot(xs, ys, marker='o', label=aligner, color=color(aligner, i))
ax.set_xlabel('reference build')
ax.set_ylabel('inherited_tp_rate')
ax.set_title(f'Read-name-inherited TP rate across builds (anchor={ANCHOR})')
ax.set_ylim(0, 1)
ax.legend(title='aligner')
fig.tight_layout()
fig.savefig('fig1_ladder_inherited_tp_rate.png', dpi=150)
plt.close(fig)
 
# ---- Figure 2: anchor's real precision/recall/F1, grouped by aligner ----
metrics = ['precision', 'recall', 'f1']
fig, ax = plt.subplots(figsize=(7, 5))
width = 0.8 / max(len(aligners), 1)
x = np.arange(len(metrics))
plotted_any = False
for i, aligner in enumerate(aligners):
    summary_path = f'bench_{aligner}_{ANCHOR}/summary.json'
    if not os.path.exists(summary_path):
        summary_path = f'bench_{aligner}_{ANCHOR}/refine.variant_summary.json'
    if not os.path.exists(summary_path):
        continue
    with open(summary_path) as f:
        s = json.load(f)
    vals = [s.get(m, 0) for m in metrics]
    ax.bar(x + i * width, vals, width, label=aligner, color=color(aligner, i))
    plotted_any = True
if plotted_any:
    ax.set_xticks(x + width * (len(aligners) - 1) / 2)
    ax.set_xticklabels(metrics)
    ax.set_ylabel('score')
    ax.set_ylim(0, 1)
    ax.set_title(f'Anchor ({ANCHOR}) real truvari accuracy, by aligner')
    ax.legend(title='aligner')
    fig.tight_layout()
    fig.savefig('fig2_anchor_accuracy.png', dpi=150)
plt.close(fig)
 
# ---- Figure 3: per-build call outcome breakdown, one panel per aligner ----
fig, axes = plt.subplots(1, len(aligners), figsize=(6 * len(aligners), 5), squeeze=False)
axes = axes[0]
for ax, aligner in zip(axes, aligners):
    builds, inherited, blocked, other = [], [], [], []
    for b in BUILD_ORDER:
        if b == ANCHOR:
            continue
        p = f'{aligner}.{b}.rname_match.json'
        if not os.path.exists(p):
            continue
        with open(p) as f:
            d = json.load(f)
        builds.append(b)
        inherited.append(d['n_inherited_tp'])
        blocked.append(d.get('n_excluded_by_type_check', 0))
        other.append(d['n_calls'] - d['n_inherited_tp'] - d.get('n_excluded_by_type_check', 0))
    x = np.arange(len(builds))
    ax.bar(x, inherited, label='inherited_TP', color='#2E7D32')
    ax.bar(x, blocked, bottom=inherited, label='excluded (type mismatch)', color='#C1440E')
    bottom2 = [a + b for a, b in zip(inherited, blocked)]
    ax.bar(x, other, bottom=bottom2, label='unmatched (no overlap)', color='#B0B0B0')
    ax.set_xticks(x)
    ax.set_xticklabels(builds, rotation=45, ha='right')
    ax.set_title(aligner)
    ax.set_ylabel('# calls')
fig.suptitle(f'Call outcome breakdown per build (anchor={ANCHOR})', y=0.99)
handles, labels = axes[0].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.93))
fig.tight_layout(rect=[0, 0, 1, 0.86])
fig.savefig('fig3_type_check_breakdown.png', dpi=150)
plt.close(fig)
 
# ---- Figure 4: precision/recall/F1 across all 7 builds, one panel per aligner ----
# NOTE: only the anchor build has real ground truth, so only it gets real
# precision/recall/F1 from truvari. The other six builds have no FP/FN
# concept under this method — there is nothing to plot for their precision
# or F1. Their point on this chart is inherited_tp_rate, a recall-like
# proxy, marked with a different style so it isn't mistaken for real recall.
from matplotlib.lines import Line2D
fig, axes = plt.subplots(1, len(aligners), figsize=(7 * len(aligners), 5.5), squeeze=False)
axes = axes[0]
for ax, aligner in zip(axes, aligners):
    rates = {}
    tsv_path = f'summary_rname.{aligner}.tsv'
    if os.path.exists(tsv_path):
        with open(tsv_path) as f:
            next(f)
            for line in f:
                build, n_calls, n_tp, rate = line.rstrip('\n').split('\t')
                rates[build] = float(rate)
 
    summary_path = f'bench_{aligner}_{ANCHOR}/summary.json'
    if not os.path.exists(summary_path):
        summary_path = f'bench_{aligner}_{ANCHOR}/refine.variant_summary.json'
    anchor_metrics = None
    if os.path.exists(summary_path):
        with open(summary_path) as f:
            anchor_metrics = json.load(f)
        rates[ANCHOR] = anchor_metrics.get('recall')
 
    xs = [b for b in BUILD_ORDER if b in rates]
    ys = [rates[b] for b in xs]
    ax.plot(xs, ys, '-', color='#1B5E20', alpha=0.5, zorder=1)
    for xi, b in enumerate(xs):
        is_anchor = (b == ANCHOR)
        ax.scatter(xi, rates[b], marker='*' if is_anchor else 'o',
                   s=260 if is_anchor else 70,
                   color='#F9A825' if is_anchor else '#1B5E20', zorder=3)
 
    if anchor_metrics and ANCHOR in xs:
        ai = xs.index(ANCHOR)
        ax.scatter(ai, anchor_metrics.get('precision'), marker='s', s=140, color='#1565C0', zorder=3)
        ax.scatter(ai, anchor_metrics.get('f1'), marker='D', s=140, color='#6A1B9A', zorder=3)
 
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels(xs, rotation=45, ha='right')
    ax.set_ylim(0, 1)
    ax.set_ylabel('score')
    ax.set_title(aligner)
 
legend_elems = [
    Line2D([0], [0], marker='o', color='#1B5E20', linestyle='-', label='recall proxy (inherited_tp_rate)'),
    Line2D([0], [0], marker='*', color='#F9A825', linestyle='None', markersize=14, label=f'real recall (anchor={ANCHOR})'),
    Line2D([0], [0], marker='s', color='#1565C0', linestyle='None', markersize=10, label=f'real precision (anchor={ANCHOR})'),
    Line2D([0], [0], marker='D', color='#6A1B9A', linestyle='None', markersize=10, label=f'real F1 (anchor={ANCHOR})'),
]
fig.suptitle('Recall proxy across all builds — real precision/recall/F1 only exist for the anchor', y=0.99, fontsize=11)
fig.legend(handles=legend_elems, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.90), fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.78])
fig.savefig('fig4_metrics_by_build.png', dpi=150)
plt.close(fig)
 
print('Wrote fig1_ladder_inherited_tp_rate.png, fig2_anchor_accuracy.png, fig3_type_check_breakdown.png, fig4_metrics_by_build.png')
PY
ANCHOR="$ANCHOR" python3 make_plots.py
# -----------------------------------------------------------------------------
# 5. UPLOADS
# -----------------------------------------------------------------------------
dx mkdir -p "$PROJ:$OUT"
dx upload summary_rname.*.tsv provenance.txt packages_version.txt fig*.png \
  ./*.rname_match.json bench_*/*.json \
  --destination "$PROJ:$OUT/" --brief
echo "-> $PROJ:$OUT"