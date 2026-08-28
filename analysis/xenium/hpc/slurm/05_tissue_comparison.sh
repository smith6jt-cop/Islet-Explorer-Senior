#!/bin/bash
#SBATCH --job-name=xen_integrate
#SBATCH --account=maigan
#SBATCH --qos=maigan
#SBATCH --mail-user=smith6jt@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --ntasks=1
#SBATCH --partition=hpg-b200
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=06:00:00
#SBATCH --output=logs/integrate_%j.out
#SBATCH --error=logs/integrate_%j.err

set -euo pipefail

export XENIUM_SAMPLE_IDS="${XENIUM_SAMPLE_IDS:-0041323,0041326}"

echo "Job $SLURM_JOB_ID: integrate samples=$XENIUM_SAMPLE_IDS on $(hostname) at $(date)"
nvidia-smi || echo "WARNING: nvidia-smi not available"

module purge
module load conda
source activate xenium_analysis

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_MAX_THREADS=$SLURM_CPUS_PER_TASK

cd "${XENIUM_ANALYSIS_ROOT:?set XENIUM_ANALYSIS_ROOT on the cluster}"
TAG=$(echo "$XENIUM_SAMPLE_IDS" | tr ',' '_')
mkdir -p "notebooks/executed/integrated_${TAG}" \
         "data/processed/integrated" \
         "figures/integrated_${TAG}"

jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=86400 \
    --output-dir="notebooks/executed/integrated_${TAG}" \
    --output="05_tissue_comparisons_${TAG}.ipynb" \
    notebooks/05_tissue_comparisons.ipynb

echo "Integration completed for $XENIUM_SAMPLE_IDS at $(date)"
