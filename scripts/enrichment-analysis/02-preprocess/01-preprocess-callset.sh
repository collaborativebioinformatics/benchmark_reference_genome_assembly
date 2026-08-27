#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../00-common/00-library.sh"
load_task "${1:?task_id is required}"

bcftools="$BENCHMARK_ENV/bin/bcftools"
output_dir="$BENCHMARK_ROOT/results/preprocessed/$ARM/$BUILD"
output="$output_dir/$ARM.$BUILD.pass.vcf.gz"
tmp_dir="$BENCHMARK_ROOT/state/tmp/preprocess/$ARM.$BUILD.${SLURM_JOB_ID:-local}"
mkdir -p "$output_dir" "$tmp_dir"

require_file "$INPUT_VCF"
require_file "$SOURCE_FASTA"
[[ $("$bcftools" query -l "$INPUT_VCF" | wc -l) -eq 1 ]] || {
  echo "Expected exactly one sample in $INPUT_VCF" >&2
  exit 1
}

"$bcftools" view -f PASS -Ou "$INPUT_VCF" \
  | "$bcftools" sort -T "$tmp_dir/sort" -Oz -o "$tmp_dir/filtered.vcf.gz"
"$BENCHMARK_ENV/bin/python" "$ANALYSIS_EXPORT_ROOT/02-preprocess/00-deduplicate-sv-vcf.py" \
  --input "$tmp_dir/filtered.vcf.gz" \
  --output "$tmp_dir/deduplicated.vcf.gz" \
  --audit "$output_dir/$ARM.$BUILD.post-sniffles-duplicates.tsv"

printf '%s\n' "$ARM.$BUILD" > "$tmp_dir/sample-name.txt"
"$bcftools" reheader -s "$tmp_dir/sample-name.txt" \
  -o "$tmp_dir/reheadered.vcf.gz" "$tmp_dir/deduplicated.vcf.gz"
mv "$tmp_dir/reheadered.vcf.gz" "$output"
"$bcftools" index -f -t "$output"
"$bcftools" index -n "$output" > "$output.records.txt"
printf '%s\t%s\t%s\t%s\n' "$ARM" "$BUILD" "$SNIFFLES_VERSION" "$INPUT_VCF" \
  > "$output.provenance.tsv"
