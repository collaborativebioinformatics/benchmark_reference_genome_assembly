# Annotation

Stages 00-05 download, lift, audit, and prepare ChromHMM and native hs1 tracks.
Stage 06 annotates TP/FP/FN variants, and stage 07 combines the per-callset
tables. Helpers share the number of the stage that invokes them.
