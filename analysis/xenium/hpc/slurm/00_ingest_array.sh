#!/bin/bash
#SBATCH --job-name=xen_ingest
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --mail-user=smith6jt@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=04:00:00
#SBATCH --output=logs/ingest_%A_%a.out
#SBATCH --error=logs/ingest_%A_%a.err

set -euo pipefail

SAMPLES=(0041323 0041326)
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
export XENIUM_SAMPLE_ID="$SAMPLE_ID"
export XENIUM_BUNDLE_DIR="${XENIUM_DATA_ROOT:?set XENIUM_DATA_ROOT}/output-XETG00298__${SAMPLE_ID}__Region_1__20260416__180001"
export XENIUM_RAW_DIR="../data/raw"

echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID: sample=$SAMPLE_ID on $(hostname) at $(date)"

module purge
module load conda
source activate xenium_analysis

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

cd "${XENIUM_ANALYSIS_ROOT:?set XENIUM_ANALYSIS_ROOT on the cluster}"
mkdir -p "notebooks/executed/${SAMPLE_ID}" data/raw

jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=14400 \
    --output-dir="notebooks/executed/${SAMPLE_ID}" \
    --output="00_ingest_${SAMPLE_ID}.ipynb" \
    notebooks/00_ingest_xenium.ipynb

echo "Ingest completed for $SAMPLE_ID at $(date)"
