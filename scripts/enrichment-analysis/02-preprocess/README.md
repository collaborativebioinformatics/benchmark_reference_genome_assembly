# Preprocessing

`01-preprocess-callset.*` keeps PASS calls, sorts them, applies exact SV
deduplication with `00-deduplicate-sv-vcf.py`, standardizes the sample name, and
indexes the result.
