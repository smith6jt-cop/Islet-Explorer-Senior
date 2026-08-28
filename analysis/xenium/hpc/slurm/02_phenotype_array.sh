#!/bin/bash
#SBATCH --job-name=xen_pheno
#SBATCH --account=maigan
#SBATCH --qos=maigan
#SBATCH --mail-user=smith6jt@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-1%1
#SBATCH --ntasks=1
#SBATCH --partition=hpg-b200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=08:00:00
#SBATCH --output=logs/pheno_%A_%a.out
#SBATCH --error=logs/pheno_%A_%a.err

set -euo pipefail

SAMPLES=(0041323 0041326)
SAMPLE_ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
export XENIUM_SAMPLE_ID="$SAMPLE_ID"

echo "Job $SLURM_JOB_ID array $SLURM_ARRAY_TASK_ID: sample=$SAMPLE_ID on $(hostname) at $(date)"
nvidia-smi || echo "WARNING: nvidia-smi not available"

module purge
module load conda
source activate xenium_analysis

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

cd "${XENIUM_ANALYSIS_ROOT:?set XENIUM_ANALYSIS_ROOT on the cluster}"
mkdir -p "notebooks/executed/${SAMPLE_ID}" \
         "figures/${SAMPLE_ID}/02_phenotyping"

jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=86400 \
    --output-dir="notebooks/executed/${SAMPLE_ID}" \
    --output="02_phenotyping_${SAMPLE_ID}.ipynb" \
    notebooks/02_phenotyping.ipynb

echo "Phenotyping completed for $SAMPLE_ID at $(date)"
