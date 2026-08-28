#!/bin/bash
#SBATCH --job-name=xen_prep
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --mail-user=smith6jt@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256gb
#SBATCH --time=12:00:00
#SBATCH --output=logs/prep_%A_%a.out
#SBATCH --error=logs/prep_%A_%a.err

set -euo pipefail

SAMPLES=(0041323 0041326)
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
export XENIUM_SAMPLE_ID="$SAMPLE_ID"

echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID: sample=$SAMPLE_ID on $(hostname) at $(date)"

module purge
module load conda
source activate xenium_analysis

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

cd "${XENIUM_ANALYSIS_ROOT:?set XENIUM_ANALYSIS_ROOT on the cluster}"
mkdir -p "notebooks/executed/${SAMPLE_ID}" \
         "data/processed/${SAMPLE_ID}" \
         "figures/${SAMPLE_ID}/01_preprocessing"

jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=86400 \
    --output-dir="notebooks/executed/${SAMPLE_ID}" \
    --output="01_preprocessing_${SAMPLE_ID}.ipynb" \
    notebooks/01_preprocessing_v2.ipynb

echo "Preprocessing completed for $SAMPLE_ID at $(date)"
