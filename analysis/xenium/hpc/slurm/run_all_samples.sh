#!/bin/bash
# Master submitter: chains ingest → preprocess → (phenotype → spatial, phenotype → integrate)
# Spatial and integration both depend only on phenotyping, so they run in parallel.
#
# Usage:
#   cd "${XENIUM_ANALYSIS_ROOT:?set XENIUM_ANALYSIS_ROOT on the cluster}"
#   bash scripts/slurm/run_all_samples.sh

set -euo pipefail

cd "$(dirname "$(readlink -f "$0")")/../.."
mkdir -p logs

ING=$(sbatch --parsable scripts/slurm/00_ingest_array.sh)
PRE=$(sbatch --parsable --dependency="afterok:$ING" scripts/slurm/01_preprocess_array.sh)
PHE=$(sbatch --parsable --dependency="afterok:$PRE" scripts/slurm/02_phenotype_array.sh)
SPA=$(sbatch --parsable --dependency="afterok:$PHE" scripts/slurm/03_spatial_array.sh)
INT=$(sbatch --parsable --dependency="afterok:$PHE" scripts/slurm/05_tissue_comparison.sh)

cat <<EOF
Submitted chained pipeline for samples 0041323, 0041326:
  ingest     $ING  (array 0-1, CPU)
  preprocess $PRE  (array 0-1, CPU)  depends on ingest
  phenotype  $PHE  (array 0-1, GPU)  depends on preprocess
  spatial    $SPA  (array 0-1, CPU)  depends on phenotype
  integrate  $INT  (single,    GPU)  depends on phenotype

Monitor:  squeue -u \$USER
Logs:     tail -f logs/prep_${PRE}_0.out   (or pheno/spatial/integrate)
EOF
