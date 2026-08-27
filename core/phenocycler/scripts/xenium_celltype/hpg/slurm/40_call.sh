#!/bin/bash
# One (donor, arm) per array task: KMeans the embedding across the k sweep, run the
# null-gated annotation, and write group calls + composition with binomial CIs.
#
# This stage is cheap ONLY because nullgate's contrast is linear in the gene set (see
# nullgate_fast.py): the whole k sweep and every null draw reduce to one sparse
# group-indicator matmul per k. Without that rewrite a 200-group, 200-draw, 8-lineage gate
# over 1.1M cells is ~320,000 passes over the count matrix and is not runnable at all.
#
#   sbatch --array=0-39%8 --dependency=afterok:$MERGE scripts/xenium_celltype/hpg/slurm/40_call.sh
#SBATCH --job-name=xct_call
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120gb
#SBATCH --time=08:00:00
#SBATCH --output=logs/xct_call_%A_%a.out
#SBATCH --error=logs/xct_call_%A_%a.err
set -uo pipefail
source "$(dirname "$(dirname "${BASH_SOURCE[0]}")")/env.sh"
cd "${XCT_ROOT}"; mkdir -p logs
xct_py                      # the manifest reads below need a python

read -r DONOR ROI ARM < <(python - "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json, os, sys
m = json.load(open(f"{os.environ['XCT_RUN']}/work_manifest.json"))
pairs = [(d["donor_id"], d["roi"], a) for d in m["donors"] for a in m["arms"]]
print(*pairs[int(sys.argv[1])])
PY
)
echo "[call] task ${SLURM_ARRAY_TASK_ID}: donor=${DONOR} arm=${ARM}"

python "${XCT_HPG}/stage_call.py" --run "${XCT_RUN}" --donor "${DONOR}" --roi "${ROI}" \
       --arm "${ARM}" --n-jobs "${SLURM_CPUS_PER_TASK:-16}"
