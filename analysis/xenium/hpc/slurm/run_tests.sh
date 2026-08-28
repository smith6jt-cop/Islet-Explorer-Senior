#!/bin/bash
#SBATCH --job-name=xenium_tests
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
#SBATCH --output=logs/tests_%j.out
#SBATCH --error=logs/tests_%j.err
#SBATCH --qos=your_qos
#SBATCH --account=your_account

# ----------------------------------------------------------------------------
# HiPerGator integration test runner.
#
# Runs the full pytest suite (including tests marked ``slow`` / ``hipergator``
# that need the real conda env and Xenium data).  Complements the lightweight
# GitHub Actions workflow: that one checks lint + unit + repo-health without
# scanpy; this one verifies that the full pipeline works end-to-end.
#
# Submit with:
#     sbatch scripts/slurm/run_tests.sh
# ----------------------------------------------------------------------------

set -euo pipefail

echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID:-local}"

# Load required modules (HiPerGator)
if command -v module >/dev/null 2>&1; then
    module purge
    module load conda
fi

# Activate conda environment
# shellcheck disable=SC1091
source activate xenium_analysis

mkdir -p logs

export NUMEXPR_MAX_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
# Avoid pytest picking up a stale cache from GitHub Actions / a different host.
export PYTEST_ADDOPTS="--cache-clear"

# Verify core scientific imports before spending time on the slow tests.
python - <<'PY'
import importlib
required = ["numpy", "pandas", "scanpy", "squidpy", "anndata", "zarr"]
missing = [m for m in required if not importlib.util.find_spec(m)]
if missing:
    raise SystemExit(f"Missing required packages in conda env: {missing}")
print("Environment OK:", ", ".join(required))
PY

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "== Ruff =="
python -m ruff check .

echo "== pytest (fast) =="
python -m pytest tests/ -v

echo "== pytest (slow/integration) =="
python -m pytest tests/ -v --run-slow --run-hipergator

echo "Job finished at: $(date)"
