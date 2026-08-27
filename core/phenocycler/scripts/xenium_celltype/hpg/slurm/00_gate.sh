#!/bin/bash
# Release gate. Proves shard+merge == monolithic KODAMA for the KODAMA build installed in
# R_LIBS_USER, on a crop of a REAL donor, before a single production core-hour is spent.
# Everything downstream carries --dependency=afterok on this job.
#
#   sbatch scripts/xenium_celltype/hpg/slurm/00_gate.sh
#SBATCH --job-name=xct_gate
#SBATCH --account=maigan
#SBATCH --qos=maigan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
#SBATCH --output=logs/xct_gate_%j.out
#SBATCH --error=logs/xct_gate_%j.err
set -uo pipefail
source "$(dirname "$(dirname "${BASH_SOURCE[0]}")")/env.sh"
cd "${XCT_ROOT}"; mkdir -p logs
xct_R

# Smallest cohort donor (6450, 255,149 cells) already prepped, cropped to --n. If prep has
# not run, fall back to the synthetic matrix, which exercises every line but not the data.
GATE_INPUT="${XCT_RUN}/donor_id=6450/roi=panc/input"
ARGS=(--n "${GATE_N:-4000}" --m "${GATE_M:-3}" --shards "${GATE_SHARDS:-3}"
      --landmarks "${GATE_LANDMARKS:-500}" --n-cores "${SLURM_CPUS_PER_TASK:-4}")
[[ -f "${GATE_INPUT}/pcs.mtx" ]] && ARGS+=(--input "${GATE_INPUT}")

echo "[gate] KODAMA lock:"; cat "${XCT_HPG}/kodama.lock" 2>/dev/null || echo "  (no lock file)"
Rscript "${XCT_HPG}/check_shard_equivalence.R" "${ARGS[@]}"
rc=$?
if [[ $rc -ne 0 ]]; then
  echo "[gate] FAILED — the sharded path is not equivalent to a monolithic KODAMA run on"
  echo "[gate] this build. Do not submit the production array. Compare kodama.lock against"
  echo "[gate] the pin in setup_envs.sh before anything else."
fi
exit $rc
