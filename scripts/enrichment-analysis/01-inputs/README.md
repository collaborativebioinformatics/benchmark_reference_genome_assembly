# Inputs

`00-download-inputs.*` downloads files listed in `../config/dnanexus-inputs.tsv`,
records checksums, and runs `01-audit-input-vcfs.py`. Existing nonempty files are
retained.
