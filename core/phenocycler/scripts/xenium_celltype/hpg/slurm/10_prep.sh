#!/bin/bash
# One donor per array task: counts -> lognorm -> PCA -> the three files the KODAMA scripts
# read, plus the feature-space kNN that BOTH arms and ALL shards would otherwise each
# recompute. Caching it here is the single largest saving in the whole plan: it is a
# function of the PC matrix alone, so it is identical across arms and shards.
#
# Sized for the largest donor in the cohort (6476, 1,109,281 cells). The 1.1M x 501 double
# neighbour arrays are ~8.8 GB and Rnanoflann makes a transient copy, so the memory here is
# real and is not a guess: --mem is set from the measured per-donor cell counts in
# work_manifest.json, and 240gb covers the top of that distribution with headroom.
#
#   sbatch --array=0-19 scripts/xenium_celltype/hpg/slurm/10_prep.sh
#SBATCH --job-name=xct_prep
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240gb
#SBATCH --time=12:00:00
#SBATCH --array=0-19
#SBATCH --output=logs/xct_prep_%A_%a.out
#SBATCH --error=logs/xct_prep_%A_%a.err
set -uo pipefail
source "$(dirname "$(dirname "${BASH_SOURCE[0]}")")/env.sh"
cd "${XCT_ROOT}"; mkdir -p logs
xct_py                      # the manifest reads below need a python

DONOR="$(python - "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json, os, sys
m = json.load(open(f"{os.environ['XCT_RUN']}/work_manifest.json"))
print(m["donors"][int(sys.argv[1])]["donor_id"])
PY
)"
LAND="$(python -c "import json,os;print(json.load(open(os.environ['XCT_RUN']+'/work_manifest.json'))['kodama']['landmarks'])")"
echo "[prep] array task ${SLURM_ARRAY_TASK_ID} -> donor ${DONOR} (landmarks=${LAND})"

python "${XCT_HPG}/stage_prep.py" --run "${XCT_RUN}" --donor "${DONOR}" \
       --n-jobs "${SLURM_CPUS_PER_TASK:-16}" || exit 1

# The kNN must come from the same Rnanoflann KODAMA.matrix would have called, so it is built
# in R, not in Python. Cached as knn_<neighbors>.rds beside pcs.mtx; kodama_merge.R reuses it.
xct_R
Rscript "${XCT_HPG}/precompute_knn.R" \
  --input "${XCT_RUN}/donor_id=${DONOR}/roi=panc/input" \
  --landmarks "${LAND}" \
  --n-cores "${SLURM_CPUS_PER_TASK:-16}"
