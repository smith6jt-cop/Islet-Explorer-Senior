#!/bin/bash
#SBATCH --job-name=ss_refine
#SBATCH --account=maigan
#SBATCH --partition=hpg-default
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=06:00:00
#SBATCH --output=/blue/maigan/smith6jt/Xenium_Analysis/logs/ss_refine_%j.out
#SBATCH --error=/blue/maigan/smith6jt/Xenium_Analysis/logs/ss_refine_%j.err
#
# DAPI-driven refinement of the 6533 same-slide seed. CPU/RAM bound (two DAPI planes,
# ~1.7 GB + ~3.9 GB, plus FFT working set) -- no GPU requested, so this does not idle a
# GPU node while it runs.
set -euo pipefail

# Submitted from inside a GPU allocation: these leak in and make srun request a
# GRES this CPU-only job does not hold.
unset SLURM_GPUS SLURM_JOB_GPUS SLURM_GPUS_ON_NODE SLURM_GPUS_PER_NODE CUDA_VISIBLE_DEVICES

REPO=/blue/maigan/smith6jt/Islet-Explorer-Senior-main/core/phenocycler
DATA=/blue/maigan/smith6jt/Xenium_Analysis/data
O="$REPO/data/integration/sameslide/donor_id=6533/roi=panc/serial=0080870"

source /blue/maigan/smith6jt/miniforge3/etc/profile.d/conda.sh
conda activate xenium_analysis
cd "$REPO"

hostname; date; nproc; free -g | head -2

python3 scripts/sameslide/refine_dapi.py \
    --qptiff "$DATA/phenocycler/0080870_6533_pLNpanc_Scan3_cells.er.qptiff" \
    --bundle "$DATA/output-XETG00298__0080870__panc__20260910__211917" \
    --out    "$O" \
    --seed-stem sameslide_v2

date
echo "=== report ==="
python3 -m json.tool "$O/sameslide_dapi_report.json" || true
