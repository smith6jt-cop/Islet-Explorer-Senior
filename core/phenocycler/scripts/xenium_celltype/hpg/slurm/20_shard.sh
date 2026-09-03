#!/bin/bash
# The expensive stage: one (donor, arm, M-shard) per array task.
#
# WHY --cpus-per-task IS SMALL, AND WHY THAT IS NOT A COMPROMISE. KODAMA 3.3 dispatches its
# Monte-Carlo loop in chunks of ceiling(M / (n.cores*4)), so a single process cannot use
# more than about sqrt(m)/2 cores on that loop no matter what it is given. MEASURED on this
# hardware at M=16: 2 cores 100.7 s, 4 cores 103.7 s, 16 cores 118.6 s -- more cores was
# slower. kodama_shard.R enforces the m >= 4*cores^2 rule and refuses to start otherwise, so
# a fat request fails loudly instead of idling. Parallelism comes from SLURM here, not from
# mclapply, which is exactly why M is sharded.
#
# THE BINDING CONSTRAINT IS MEMORY, NOT CPU. At ~48gb a task, maigan's 625 GB investment
# holds ~13 concurrent and clive's 312 GB ~6, using ~38 of 120 CPUs. Raising cpus-per-task
# would not help (see above); the lever that would is the res_only patch described in
# submit_all.sh, which removes the shard's redundant kNN and dissimilarity entirely.
#
# The array size is n_shard_tasks from work_manifest.json — preflight prints it. Submit as
#   sbatch --array=0-$((N-1))%${XCT_MAX_CONCURRENT:-13} .../20_shard.sh
# The %N throttle is what keeps a burst allocation from being requested all at once.
#SBATCH --job-name=xct_shard
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/xct_shard_%A_%a.out
#SBATCH --error=logs/xct_shard_%A_%a.err
set -uo pipefail
source "$(dirname "$(dirname "${BASH_SOURCE[0]}")")/env.sh"
cd "${XCT_ROOT}"; mkdir -p logs
xct_py                      # the manifest reads below need a python

T="$(xct_task shard_tasks "${SLURM_ARRAY_TASK_ID}")"
get() { python -c "import json,sys;print(json.loads(sys.argv[1])['$1'])" "${T}"; }
DONOR="$(get donor_id)"; ROI="$(get roi)"; ARM="$(get arm)"
SHARD="$(get shard)"; M="$(get m)"
LAND="$(python -c "import json,os;print(json.load(open(os.environ['XCT_RUN']+'/work_manifest.json'))['kodama']['landmarks'])")"
SEED="$(python -c "import json,os;print(json.load(open(os.environ['XCT_RUN']+'/work_manifest.json'))['kodama']['seed'])")"
DEST="${XCT_RUN}/donor_id=${DONOR}/roi=${ROI}"
echo "[shard] task ${SLURM_ARRAY_TASK_ID}: donor=${DONOR} arm=${ARM} shard=${SHARD} m=${M}"

xct_R
Rscript "${XCT_HPG}/kodama_shard.R" \
  --input "${DEST}/input" --output "${DEST}" --arm "${ARM}" \
  --m "${M}" --shard "${SHARD}" --base-seed "${SEED}" --landmarks "${LAND}" \
  --n-cores "${SLURM_CPUS_PER_TASK:-2}"
