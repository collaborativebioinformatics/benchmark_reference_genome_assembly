#!/usr/bin/env bash
set -euo pipefail

main() {
  echo "Value of reads: '$reads'"
  echo "Value of reference_build: '$reference_build'"
  echo "Value of aligner: '$aligner'"
  echo "Value of preset: '$preset'"
  echo "Value of min_sv_len: '${min_sv_len:-50}'"

  # -----------------------------------------------------------------------
  # 0. Fetch the input reads file (dx-download-all-inputs pulls named inputs
  #    into $HOME/in/<input_name>/<filename>)
  # -----------------------------------------------------------------------
  dx-download-all-inputs --parallel

  READS_PATH=$(find "$HOME/in/reads" -type f | head -n 1)
  echo "Reads downloaded to: $READS_PATH"

  # -----------------------------------------------------------------------
  # 1. Install tools (same pattern as sniffles_mm2.sh — micromamba env)
  # -----------------------------------------------------------------------
  if ! command -v sniffles >/dev/null 2>&1; then
    echo "installing tools"
    if ! command -v micromamba >/dev/null 2>&1; then
      curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
      export PATH="$PWD/bin:$PATH"
    fi
    micromamba create -y -p "$HOME/env" -c conda-forge -c bioconda \
      samtools minimap2 winnowmap sniffles rasusa
  fi
  export PATH="$HOME/env/bin:$PATH"

  # -----------------------------------------------------------------------
  # 2. Fetch the chosen reference genome from the project's shared assets
  # -----------------------------------------------------------------------
  PROJ="Group4_2026"
  REL="/reference-assets/releases/2026-08-25"
  REF_REMOTE="${PROJ}:${REL}/genomes/${reference_build}/reference/${reference_build}.reference.tar.zst"

  dx download "$REF_REMOTE" -o ref.tar.zst
  unzstd -f ref.tar.zst
  tar -xf ref.tar
  REF=$(find . -maxdepth 3 -name "*.fa" | head -n 1)
  echo "Using reference FASTA: $REF"

  # -----------------------------------------------------------------------
  # 3. Downsample reads to a target coverage (optional, on by default)
  #
  #    This normalizes coverage across different input datasets/technologies
  #    so SV-calling comparisons across reference builds are fair, per the
  #    project's benchmarking design.
  # -----------------------------------------------------------------------
  ALIGN_INPUT="$READS_PATH"

  if [ "${downsample:-true}" == "true" ]; then
    echo "Downsampling to ${target_coverage:-15}x (genome size ${genome_size:-3.1g})"
    DOWNSAMPLED="downsampled.$(basename "$READS_PATH")"

    rasusa reads \
      --coverage "${target_coverage:-15}" \
      --genome-size "${genome_size:-3.1g}" \
      -o "$DOWNSAMPLED" \
      "$READS_PATH"

    ALIGN_INPUT="$DOWNSAMPLED"
    echo "Downsampled reads written to: $ALIGN_INPUT"
  else
    echo "Downsampling disabled; using full input read set."
  fi

  # -----------------------------------------------------------------------
  # 4. Align
  # -----------------------------------------------------------------------
  PREFIX="input.${reference_build}.${aligner}.${preset}"
  BAM="${PREFIX}.bam"

  if [ "$aligner" == "minimap2" ]; then
    minimap2 -ax "$preset" -t "$(nproc)" "$REF" "$ALIGN_INPUT" \
      | samtools sort -@ "$(nproc)" -o "$BAM"
  else
    winnowmap -ax "$preset" -t "$(nproc)" "$REF" "$ALIGN_INPUT" \
      | samtools sort -@ "$(nproc)" -o "$BAM"
  fi
  samtools index "$BAM"

  # -----------------------------------------------------------------------
  # 6. Call structural variants with Sniffles2
  # -----------------------------------------------------------------------
  VCF="${PREFIX}.sniffles.vcf"
  SNF="${PREFIX}.sniffles.snf"

  sniffles --input "$BAM" --reference "$REF" \
    --vcf "$VCF" --snf "$SNF" \
    --minsvlen "${min_sv_len:-50}" --output-rnames --threads "$(nproc)"

  # -----------------------------------------------------------------------
  # 7. Upload outputs back through the DNAnexus job output system
  # -----------------------------------------------------------------------
  sv_vcf=$(dx upload "$VCF" --brief)
  sv_snf=$(dx upload "$SNF" --brief)
  alignment_bam=$(dx upload "$BAM" --brief)

  dx-jobutil-add-output sv_vcf "$sv_vcf" --class=file
  dx-jobutil-add-output sv_snf "$sv_snf" --class=file
  dx-jobutil-add-output alignment_bam "$alignment_bam" --class=file
}
