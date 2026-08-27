#!/usr/bin/env bash
set -euo pipefail

ANALYSIS_EXPORT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
ANALYSIS_WORKSPACE_ROOT=$(cd "$ANALYSIS_EXPORT_ROOT/../.." && pwd)
export ANALYSIS_WORKSPACE_ROOT
cd "$ANALYSIS_WORKSPACE_ROOT"

input_audit=$(sbatch --parsable benchmarking/analysis-export/01-inputs/00-download-inputs.slurm)
preprocess=$(sbatch --parsable --dependency="afterok:$input_audit" benchmarking/analysis-export/02-preprocess/01-preprocess-callsets.slurm)
liftover=$(sbatch --parsable --dependency="afterok:$preprocess" benchmarking/analysis-export/03-liftover/01-lift-callsets-to-hs1.slurm)
liftover_audit=$(sbatch --parsable --dependency="afterok:$liftover" benchmarking/analysis-export/03-liftover/02-audit-liftover.slurm)
benchmark=$(sbatch --parsable --dependency="afterok:$liftover_audit" benchmarking/analysis-export/04-benchmark/00-benchmark-callsets.slurm)
summary=$(sbatch --parsable --dependency="afterok:$benchmark" benchmarking/analysis-export/04-benchmark/01-summarize-benchmarks.slurm)
annotation=$(sbatch --parsable --dependency="afterok:$benchmark" benchmarking/analysis-export/05-annotations/06-annotate-callsets.slurm)
annotation_finalize=$(sbatch --parsable --dependency="afterok:$annotation" benchmarking/analysis-export/05-annotations/07-finalize-annotations.slurm)
masks=$(sbatch --parsable --dependency="afterok:$liftover_audit" benchmarking/analysis-export/06-masks-enrichment/00-build-comparability-masks.slurm)
enrichment=$(sbatch --parsable --dependency="afterok:$annotation_finalize:$masks" benchmarking/analysis-export/06-masks-enrichment/01-compute-oe-enrichment.slurm)
retention=$(sbatch --parsable --dependency="afterok:$enrichment" benchmarking/analysis-export/06-masks-enrichment/02-summarize-mask-retention.slurm)
strict_clusters=$(sbatch --parsable --dependency="afterok:$enrichment" benchmarking/analysis-export/07-cross-reference/00-cross-reference-cluster.slurm)
strict_cluster_enrichment=$(sbatch --parsable --dependency="afterok:$strict_clusters" benchmarking/analysis-export/07-cross-reference/01-compute-cluster-enrichment.slurm)
relaxed_enrichment=$(sbatch --parsable --dependency="afterok:$enrichment" benchmarking/analysis-export/06-masks-enrichment/03-compute-relaxed-enrichment.slurm)
relaxed_clusters=$(sbatch --parsable --dependency="afterok:$relaxed_enrichment" benchmarking/analysis-export/07-cross-reference/00-cross-reference-cluster-relaxed.slurm)
relaxed_cluster_enrichment=$(sbatch --parsable --dependency="afterok:$relaxed_clusters:$strict_cluster_enrichment" benchmarking/analysis-export/07-cross-reference/01-compute-cluster-enrichment-relaxed.slurm)
figures=$(sbatch --parsable --dependency="afterok:$summary:$retention:$strict_cluster_enrichment:$relaxed_cluster_enrichment" benchmarking/analysis-export/08-figures/01-make-figures.slurm)
final_audit=$(sbatch --parsable --dependency="afterok:$figures" benchmarking/analysis-export/09-orchestration/00-audit-results.slurm)

printf 'stage\tjob_id\n'
printf 'input_audit\t%s\n' "$input_audit"
printf 'preprocess\t%s\n' "$preprocess"
printf 'liftover\t%s\n' "$liftover"
printf 'liftover_audit\t%s\n' "$liftover_audit"
printf 'benchmark\t%s\n' "$benchmark"
printf 'summary\t%s\n' "$summary"
printf 'annotation\t%s\n' "$annotation"
printf 'annotation_finalize\t%s\n' "$annotation_finalize"
printf 'masks\t%s\n' "$masks"
printf 'enrichment\t%s\n' "$enrichment"
printf 'retention\t%s\n' "$retention"
printf 'strict_clusters\t%s\n' "$strict_clusters"
printf 'strict_cluster_enrichment\t%s\n' "$strict_cluster_enrichment"
printf 'relaxed_enrichment\t%s\n' "$relaxed_enrichment"
printf 'relaxed_clusters\t%s\n' "$relaxed_clusters"
printf 'relaxed_cluster_enrichment\t%s\n' "$relaxed_cluster_enrichment"
printf 'figures\t%s\n' "$figures"
printf 'final_audit\t%s\n' "$final_audit"
