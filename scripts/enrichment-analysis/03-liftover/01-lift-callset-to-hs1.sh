#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
load_task "${1:?task_id is required}"

bcftools="$BENCHMARK_ENV/bin/bcftools"
bgzip="$BENCHMARK_ENV/bin/bgzip"
liftover_sv="$BENCHMARK_ROOT/tools/liftoverSV/bin/liftoverSV.py"
input="$BENCHMARK_ROOT/results/preprocessed/$ARM/$BUILD/$ARM.$BUILD.pass.vcf.gz"
output_dir="$BENCHMARK_ROOT/results/lifted/$ARM/$BUILD"
output="$output_dir/$ARM.$BUILD.to-hs1.vcf.gz"
unmapped="$output_dir/$ARM.$BUILD.to-hs1.unmapped.tsv"
rejected="$output_dir/$ARM.$BUILD.to-hs1.canonicalization-rejected.tsv"
tmp_dir="$BENCHMARK_ROOT/state/tmp/liftover/$ARM.$BUILD.${SLURM_JOB_ID:-local}"
mkdir -p "$output_dir" "$tmp_dir"
mkdir -p "$tmp_dir/work"

require_file "$input"
require_file "$HS1_FASTA"

if [[ "$BUILD" == hs1 ]]; then
  candidate="$input"
  printf 'native_hs1\n' > "$unmapped"
else
  require_file "$CHAIN_TO_HS1"
  chain="$CHAIN_TO_HS1"
  if [[ "$CHAIN_TO_HS1" == *.gz ]]; then
    chain="$tmp_dir/input.chain"
    gzip -cd "$CHAIN_TO_HS1" > "$chain"
  fi
  "$BENCHMARK_ENV/bin/python" "$liftover_sv" \
    --chain "$chain" \
    --input-file "$input" \
    --ref-fasta-seq "$HS1_FASTA" \
    --output-base-name "$tmp_dir/lifted" \
    --n-workers "${SLURM_CPUS_PER_TASK:-4}" \
    --chunk-size 5000 \
    --percent 0.05 \
    --tmp-dir "$tmp_dir/work"
  require_file "$tmp_dir/lifted.sort.vcf.gz"
  "$bcftools" view -h "$input" | grep '^##INFO=' > "$tmp_dir/input-info-header.txt"
  gzip -cd "$tmp_dir/lifted.sort.vcf.gz" \
    | awk -v info="$tmp_dir/input-info-header.txt" '
        /^#CHROM/ {
          while ((getline line < info) > 0) print line
          close(info)
          print
          next
        }
        /^#/ {print}
      ' > "$tmp_dir/repaired-header.txt"
  {
    cat "$tmp_dir/repaired-header.txt"
    gzip -cd "$tmp_dir/lifted.sort.vcf.gz" | grep -v '^#'
  } | "$bgzip" -c > "$tmp_dir/lifted.reheadered.vcf.gz"
  candidate="$tmp_dir/lifted.reheadered.vcf.gz"
  if [[ -f "$tmp_dir/lifted.unmapped" ]]; then
    mv "$tmp_dir/lifted.unmapped" "$unmapped"
  else
    printf 'No unmapped report emitted\n' > "$unmapped"
  fi
fi

"$BENCHMARK_ENV/bin/python" "$ANALYSIS_EXPORT_ROOT/03-liftover/00-canonicalize-hs1-sv.py" \
  --input "$candidate" \
  --output "$tmp_dir/canonical.vcf.gz" \
  --reference "$HS1_FASTA" \
  --rejected "$rejected"
"$bcftools" sort -T "$tmp_dir/final-sort" -Oz -o "$output" "$tmp_dir/canonical.vcf.gz"
"$bcftools" index -f -t "$output"
"$bcftools" index -n "$output" > "$output.records.txt"
