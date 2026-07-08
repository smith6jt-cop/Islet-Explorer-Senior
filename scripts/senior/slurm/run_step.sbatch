#!/bin/bash
# Generic SLURM wrapper for one Senior-pipeline step, so long runs survive an
# interactive session expiring. Submit any pipeline command through this:
#
#   sbatch scripts/senior/slurm/run_step.sbatch python scripts/senior/restore_normalize.py --cells-dir data/cells_redsea ...
#
# GPU steps (scVI) get a B200; CPU steps ignore the unused GPU request. Override
# resources per-submit, e.g.:  sbatch --time=48:00:00 --job-name=scvi run_step.sbatch python scripts/senior/train_scvi_senior.py
#
#SBATCH --job-name=senior_step
#SBATCH --partition=hpg-b200
#SBATCH --qos=maigan
#SBATCH --account=maigan
#SBATCH --gres=gpu:b200:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
cd /blue/maigan/smith6jt/Islet-Explorer-Senior-main
mkdir -p logs

source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
# canonical env from scvi_env.yml; falls back to xenium_analysis if scvi-env absent
if conda env list | grep -q "/envs/scvi-env"; then conda activate scvi-env; else conda activate xenium_analysis; fi

echo "[run_step] host=$(hostname) job=${SLURM_JOB_ID:-none} env=${CONDA_DEFAULT_ENV} start=$(date -Is)"
echo "[run_step] cmd: $*"
"$@"
rc=$?
echo "[run_step] exit=$rc end=$(date -Is)"
exit $rc
