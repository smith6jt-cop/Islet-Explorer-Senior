#!/bin/bash
# Cohort assembly + the verification that decides whether this run is a result.
# Runs on the investment QoS, not burst: it is small, and it is the job whose absence would
# leave 20 donors of output with nothing saying whether they are complete.
#
#   sbatch --dependency=afterok:$CALL scripts/xenium_celltype/hpg/slurm/50_assemble.sh
#SBATCH --job-name=xct_assemble
#SBATCH --account=maigan
#SBATCH --qos=maigan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --time=06:00:00
#SBATCH --output=logs/xct_assemble_%j.out
#SBATCH --error=logs/xct_assemble_%j.err
set -uo pipefail
source "$(dirname "$(dirname "${BASH_SOURCE[0]}")")/env.sh"
cd "${XCT_ROOT}"; mkdir -p logs
xct_py
python "${XCT_HPG}/verify_run.py" --run "${XCT_RUN}" --assemble
rc=$?
echo "[assemble] verify_run.py exit ${rc} (non-zero = the run is NOT complete; read the report)"
exit $rc
