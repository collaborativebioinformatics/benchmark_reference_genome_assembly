import os
import tempfile

import pandas as pd
import streamlit as st

# DNAnexus backend functions: submit a job, wait for it, and fetch results
from src.sv_app.dnanexus_client import (
    submit_job,
    wait_for_job,
    get_output_files,
    download_file,
)

# VCF parsing function: turns raw Sniffles2 output into a tidy DataFrame
from src.sv_app.vcf_parser import parse_vcf

st.set_page_config(
    page_title="SV Discovery",
    page_icon="🧬",
    layout="wide"
)

st.header("SV Discovery 🧬")

st.markdown("""This tool aligns your long-read sequencing data against a
**reference genome build of your choice** and calls structural variants
with **Sniffles2**. Upload a FASTQ file to get started!""")
st.divider()

# reference builds available in the shared project's reference-assets folder
REFERENCE_BUILDS = ["hg15", "hg16", "hg17", "hg18", "hg19", "hg38", "hs1"]
ALIGNERS = ["minimap2", "winnowmap"]
PRESETS = ["map-hifi", "map-pb", "map-ont"]


def main():
    st.title("Structural Variant Discovery")

    # ---------------------------------------------------------------
    # Inputs — mirrors the dxapp.json inputSpec exactly, so what the
    # user sees here maps 1:1 to what the DNAnexus applet expects
    # ---------------------------------------------------------------
    uploaded_file = st.file_uploader(
        "Upload a FASTQ file or drag and drop it here",
        type=["fastq", "fq", "gz"],
        key="fastq_uploader",
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        reference_build = st.selectbox("Reference genome build", REFERENCE_BUILDS, index=REFERENCE_BUILDS.index("hg38"))
    with col2:
        aligner = st.selectbox("Aligner", ALIGNERS)
    with col3:
        preset = st.selectbox("Sequencing preset", PRESETS)

    with st.expander("Advanced options"):
        downsample = st.checkbox("Downsample to a target coverage", value=True)
        target_coverage = st.number_input("Target coverage (x)", min_value=1.0, value=15.0, step=1.0)
        genome_size = st.text_input("Genome size (for coverage calculation)", value="3.1g")
        min_sv_len = st.number_input("Minimum SV length (bp)", min_value=1, value=50, step=1)

    # ---------------------------------------------------------------
    # Run button — uploads the file, submits the job, and waits
    # ---------------------------------------------------------------
    if uploaded_file and st.button("Run SV discovery"):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fastq") as tmp:
            tmp.write(uploaded_file.getvalue())
            tmp_path = tmp.name

        try:
            status_placeholder = st.empty()

            with st.spinner("Submitting job to DNAnexus..."):
                job = submit_job(
                    local_fastq_path=tmp_path,
                    reference_build=reference_build,
                    aligner=aligner,
                    preset=preset,
                    downsample=downsample,
                    target_coverage=target_coverage,
                    genome_size=genome_size,
                    min_sv_len=min_sv_len,
                )

            # poll job status, updating the placeholder each tick so the
            # user can see it move through queued -> running -> done
            def on_tick(state):
                status_placeholder.info(f"Job status: **{state}**")

            with st.spinner("Aligning reads and calling structural variants... this can take a while for a full genome."):
                final_state = wait_for_job(job, on_tick=on_tick)

            if final_state != "done":
                st.error(f"Job did not complete successfully (final state: {final_state}). "
                         f"Check the job log on DNAnexus for details.")
            else:
                st.success("Done! Results are ready below.")
                outputs = get_output_files(job)
                st.session_state["last_job_outputs"] = outputs

                # VCF is small (KB-MB range) — download it right away so we
                # can chart it immediately, without an extra click
                if "sv_vcf" in outputs:
                    vcf_info = outputs["sv_vcf"]
                    local_vcf_path = os.path.join(tempfile.gettempdir(), vcf_info["name"])
                    with st.spinner("Fetching variant calls for visualization..."):
                        download_file(vcf_info["file_id"], local_vcf_path)
                    st.session_state["last_vcf_df"] = parse_vcf(local_vcf_path)

        finally:
            os.unlink(tmp_path)  # clean up temp file

    # ---------------------------------------------------------------
    # Results — show once a job has completed in this session
    # ---------------------------------------------------------------
    if "last_job_outputs" in st.session_state:
        outputs = st.session_state["last_job_outputs"]

        st.subheader("Results")

        # -----------------------------------------------------------
        # Visualizations — only shown once we have a parsed VCF
        # -----------------------------------------------------------
        if "last_vcf_df" in st.session_state:
            df = st.session_state["last_vcf_df"]

            if df.empty:
                st.info("No structural variants were called for this sample.")
            else:
                st.markdown("### Structural variant summary")

                # --- top-line metrics ---
                m1, m2, m3 = st.columns(3)
                m1.metric("Total SVs called", len(df))
                m2.metric("Chromosomes affected", df["CHROM"].nunique())
                pass_count = (df["FILTER"] == "PASS").sum()
                m3.metric("Passing filter", f"{pass_count}/{len(df)}")

                chart_col1, chart_col2 = st.columns(2)

                # --- SV type breakdown (bar chart) ---
                with chart_col1:
                    st.markdown("**Variant types**")
                    type_counts = df["SVTYPE"].value_counts()
                    st.bar_chart(type_counts)

                # --- SV size distribution (bar chart on absolute length) ---
                with chart_col2:
                    st.markdown("**Variant sizes (bp, absolute value)**")
                    sizes = df["SVLEN"].dropna().abs()
                    if not sizes.empty:
                        st.bar_chart(sizes.reset_index(drop=True))
                    else:
                        st.write("No size information available for these variants.")

                # --- position across the genome (scatter: chromosome vs position) ---
                st.markdown("**Where each variant sits in the genome**")
                position_df = df[["CHROM", "POS", "SVTYPE"]].copy()
                st.scatter_chart(position_df, x="POS", y="CHROM", color="SVTYPE")

                # --- full table, sortable/searchable by Streamlit natively ---
                st.markdown("**All variant calls**")
                st.dataframe(
                    df[["CHROM", "POS", "SVTYPE", "SVLEN", "SUPPORT", "VAF", "FILTER", "GENOTYPE"]],
                    use_container_width=True,
                )

                # let the user grab the parsed table as CSV, separate from the raw VCF download below
                csv_bytes = df.to_csv(index=False).encode()
                st.download_button(
                    "⬇ Download variant table as CSV",
                    csv_bytes,
                    "structural_variants.csv",
                    "text/csv",
                )

            st.divider()

        for key, label in [
            ("sv_vcf", "Structural variant calls (VCF)"),
            ("sv_snf", "Sniffles snapshot (SNF)"),
            ("alignment_bam", "Alignment (BAM)"),
        ]:
            if key in outputs:
                info = outputs[key]
                size_mb = info["size_bytes"] / (1024 * 1024)
                st.write(f"**{label}**: {info['name']} ({size_mb:.1f} MB)")

                if st.button(f"Prepare download: {info['name']}", key=f"prep_{key}"):
                    local_path = os.path.join(tempfile.gettempdir(), info["name"])
                    with st.spinner(f"Downloading {info['name']} from DNAnexus..."):
                        download_file(info["file_id"], local_path)
                    with open(local_path, "rb") as f:
                        st.download_button(
                            f"⬇ Download {info['name']}",
                            f,
                            file_name=info["name"],
                            key=f"dl_{key}",
                        )
    else:
        st.info("Upload a FASTQ file and click 'Run SV discovery' to get started")


if __name__ == "__main__":
    main()
