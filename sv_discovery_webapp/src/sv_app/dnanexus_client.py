"""
Thin wrapper around dxpy so the Streamlit app never has to touch DNAnexus
concepts (projects, job states, file IDs) directly. Mirrors the pattern of
vcf_parser.py / clinvar_merge.py: pure functions in, pure data out.
"""

import time
import dxpy

# Name of the applet built from dxapp.json (see: dx build sv_discovery_app/)
APPLET_NAME = "sv_discovery"

# The shared team project where the applet and reference assets live
PROJECT_NAME = "Group4_2026"

# How often to check job status while waiting (seconds)
POLL_INTERVAL_SECONDS = 15


def get_project():
    """Resolve the project name to a dxpy project handle."""
    project_id = dxpy.find_one_project(name=PROJECT_NAME, zero_ok=False)["id"]
    return dxpy.DXProject(project_id)


def get_applet(project):
    """Locate the sv_discovery applet inside the project."""
    applet_id = dxpy.find_one_data_object(
        classname="applet",
        name=APPLET_NAME,
        project=project.get_id(),
        zero_ok=False,
    )["id"]
    return dxpy.DXApplet(applet_id)


def submit_job(local_fastq_path, reference_build, aligner, preset,
               downsample, target_coverage, genome_size, min_sv_len):
    """
    Upload the user's FASTQ and launch the sv_discovery applet with the
    options chosen in the UI. Returns the DXJob handle immediately —
    does not block.
    """
    project = get_project()
    applet = get_applet(project)

    # Upload the user's file into the project so the applet can read it
    uploaded_reads = dxpy.upload_local_file(
        filename=local_fastq_path,
        project=project.get_id(),
        folder="/webapp_uploads",
        parents=True,
    )

    job = applet.run(
        applet_input={
            "reads": dxpy.dxlink(uploaded_reads),
            "reference_build": reference_build,
            "aligner": aligner,
            "preset": preset,
            "downsample": downsample,
            "target_coverage": target_coverage,
            "genome_size": genome_size,
            "min_sv_len": min_sv_len,
        },
        project=project.get_id(),
        folder="/webapp_results",
    )
    return job


def wait_for_job(job, on_tick=None):
    """
    Block until the job finishes. Calls on_tick(state) each poll so the
    caller (Streamlit) can update a spinner/status message.
    """
    while True:
        state = job.describe()["state"]
        if on_tick:
            on_tick(state)

        if state in ("done", "failed", "terminated"):
            return state

        time.sleep(POLL_INTERVAL_SECONDS)


def get_output_files(job):
    """
    Once a job is done, return local download URLs (via dxpy file links)
    for each output defined in dxapp.json's outputSpec.
    """
    output = job.describe()["output"]
    files = {}
    for name in ("sv_vcf", "sv_snf", "alignment_bam"):
        if name in output:
            file_id = output[name]["$dnanexus_link"]
            dxfile = dxpy.DXFile(file_id)
            files[name] = {
                "name": dxfile.describe()["name"],
                "file_id": file_id,
                "size_bytes": dxfile.describe()["size"],
            }
    return files


def download_file(file_id, local_path):
    """Download a DNAnexus file object to local disk for the download button."""
    dxpy.download_dxfile(file_id, local_path)
    return local_path
