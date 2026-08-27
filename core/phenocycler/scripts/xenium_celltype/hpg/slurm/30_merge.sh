#!/bin/bash
# One (donor, arm) per array task: rbind the shards, recompute the dissimilarity over the
# full M, and embed. This is the stage that DOES scale with cores -- the dissimilarity loop
# iterates over cells, not over M, so it chunks at 16 and uses everything it is given -- so
# it is the one place a fat request is honest.
#
# Checkpointed at row-block granularity (--block 50000 -> ~23 blocks for the largest donor).
# A task killed at hour 23 resumes at the next unfinished block; nothing is recomputed.
#
#   sbatch --array=0-39%4 --dependency=afterok:$SHARD scripts/xenium_celltype/hpg/slurm/30_merge.sh
#SBATCH --job-name=xct_merge
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/xct_merge_%A_%a.out
#SBATCH --error=logs/xct_merge_%A_%a.err
set -uo pipefail
source "$(dirname "$(dirname "${BASH_SOURCE[0]}")")/env.sh"
cd "${XCT_ROOT}"; mkdir -p logs
xct_py                      # the manifest reads below need a python

read -r DONOR ROI ARM M MSH LAND SEED < <(python - "${SLURM_ARRAY_TASK_ID}" <<'PY'
import json, os, sys
m = json.load(open(f"{os.environ['XCT_RUN']}/work_manifest.json"))
pairs = [(d["donor_id"], d["roi"], a) for d in m["donors"] for a in m["arms"]]
d, roi, arm = pairs[int(sys.argv[1])]
k = m["kodama"]
print(d, roi, arm, k["M"], k["m_shard"], k["landmarks"], k["seed"])
PY
)
DEST="${XCT_RUN}/donor_id=${DONOR}/roi=${ROI}"
echo "[merge] task ${SLURM_ARRAY_TASK_ID}: donor=${DONOR} arm=${ARM} M=${M}"

xct_R
Rscript "${XCT_HPG}/kodama_merge.R" \
  --input "${DEST}/input" --output "${DEST}" --arm "${ARM}" \
  --M "${M}" --m "${MSH}" --landmarks "${LAND}" --seed "${SEED}" \
  --n-cores "${SLURM_CPUS_PER_TASK:-32}" --block "${XCT_MERGE_BLOCK:-50000}"
