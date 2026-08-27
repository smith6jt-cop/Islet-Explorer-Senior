#!/usr/bin/env bash
# =============================================================================
# submit_all.sh — the whole production run as one dependency chain.
#
#   gate -> preflight -> prep -> calibrate -> shard -> merge -> call -> assemble
#
# Every stage is `--dependency=afterok` on the one before, so a failure stops
# the chain instead of feeding half-written inputs to the next stage. Every
# stage is idempotent: re-running this script after a failure resubmits only
# what has not been written, because each script checks its own output first.
#
# THE TWO ALLOCATIONS. maigan (CPU=80, MEM=625GB) and clive (CPU=40, MEM=312GB)
# are separate SLURM accounts, so one array cannot span both. The shard array is
# therefore submitted TWICE against disjoint index ranges — SHARD_SPLIT sets the
# boundary — and both halves write into the same run directory. verify_run.py
# checks the union, so a half that dies is visible as missing shards rather than
# as a quietly smaller M.
#
# WHY THE SHARD ARRAY IS MEMORY-THROTTLED AND NOT CPU-THROTTLED. KODAMA's M loop
# runs min(chunk_size, n.cores) iterations at a time with
# chunk_size = ceiling(M / (n.cores*4)), so one process cannot use more than
# ~sqrt(m)/2 cores no matter what it is given (MEASURED at M=16: 2 cores 100.7 s,
# 4 cores 103.7 s, 16 cores 118.6 s — more cores was slower). Throughput comes
# from running many small tasks, and how many run at once is set by memory:
# ~13 on maigan and ~6 on clive at 48gb each. That leaves CPUs idle, and the
# only real fix is upstream — see THE ONE PATCH WORTH MAKING at the bottom.
#
# Usage:
#   bash submit_all.sh --dry-run          # print every sbatch, submit nothing
#   bash submit_all.sh                    # submit the chain
#   bash submit_all.sh --from shard       # resume at a stage
#   XCT_RUN_ID=2026-08-26 bash submit_all.sh
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export XCT_RUN_ID="${XCT_RUN_ID:-$(date -u +%Y-%m-%dT%H%M%SZ)}"
source "${HERE}/env.sh"

DRY=0; FROM="gate"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY=1; shift ;;
    --from) FROM="$2"; shift 2 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done

ACCT_A="${XCT_ACCOUNT_A:-maigan}";  QOS_A="${XCT_QOS_A:-maigan-b}"
ACCT_B="${XCT_ACCOUNT_B:-clive}";   QOS_B="${XCT_QOS_B:-clive-b}"
CONC_A="${XCT_CONC_A:-13}"          # 625 GB / 48 GB
CONC_B="${XCT_CONC_B:-6}"           # 312 GB / 48 GB

mkdir -p "${XCT_RUN}" "${XCT_ROOT}/logs"
S="${XCT_HPG}/slurm"

run() {  # run <label> <sbatch args...> ; echoes the job id (or a placeholder when --dry-run)
  local label="$1"; shift
  if [[ $DRY -eq 1 ]]; then
    echo "  sbatch $*" >&2
    echo "DRYRUN_${label}"
  else
    local id; id="$(sbatch --parsable "$@")" || { echo "[submit] ${label} FAILED" >&2; exit 1; }
    echo "[submit] ${label} -> ${id}" >&2
    echo "${id}"
  fi
}
stage_at_or_after() {  # is $1 at or after ${FROM} in the pipeline order?
  local order=(gate preflight prep calibrate shard merge call assemble) i j
  for i in "${!order[@]}"; do [[ "${order[$i]}" == "${FROM}" ]] && break; done
  for j in "${!order[@]}"; do [[ "${order[$j]}" == "$1" ]] && break; done
  [[ $j -ge $i ]]
}

echo "run id     : ${XCT_RUN_ID}"
echo "run dir    : ${XCT_RUN}"
echo "code       : ${XCT_CODE} @ $(git -C "${XCT_CODE}" rev-parse --short HEAD 2>/dev/null || echo '?')"
echo "KODAMA lock: $(python -c "import json;print(json.load(open('${HERE}/kodama.lock'))['kodama_sha'][:12])" 2>/dev/null || echo 'MISSING — run setup_envs.sh')"
echo

# --- 0. preflight runs on the LOGIN NODE: it is seconds, and every array size below --------
#        comes out of the manifest it writes, so it cannot be a queued dependency.
if stage_at_or_after preflight; then
  echo "[submit] preflight (login node)"
  if [[ $DRY -eq 0 ]]; then
    source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate xct_py
    python "${XCT_HPG}/preflight.py" --out "${XCT_RUN}" \
      --M "${XCT_M:-96}" --m-shard "${XCT_M_SHARD:-32}" \
      --landmarks "${XCT_LANDMARKS:-2000}" --n-comps "${XCT_NCOMPS:-50}" || exit 1
    conda deactivate
  fi
fi
if [[ -f "${XCT_RUN}/work_manifest.json" ]]; then
  N_DONORS=$(python -c "import json;print(json.load(open('${XCT_RUN}/work_manifest.json'))['n_donors'])")
  N_SHARD=$(python -c "import json;print(json.load(open('${XCT_RUN}/work_manifest.json'))['n_shard_tasks'])")
  N_ARMS=$(python -c "import json;print(len(json.load(open('${XCT_RUN}/work_manifest.json'))['arms']))")
else
  N_DONORS=20; N_SHARD=120; N_ARMS=2
fi
N_PAIR=$((N_DONORS * N_ARMS))
SPLIT="${SHARD_SPLIT:-$(( N_SHARD * 2 / 3 ))}"   # maigan gets the larger share

DEP=""
if stage_at_or_after gate; then
  GATE=$(run gate --account="${ACCT_A}" --qos="${XCT_QOS_INVEST:-${ACCT_A}}" "${S}/00_gate.sh")
  DEP="--dependency=afterok:${GATE}"
fi
if stage_at_or_after prep; then
  PREP=$(run prep ${DEP} --array=0-$((N_DONORS-1))%"${CONC_A}" "${S}/10_prep.sh")
  DEP="--dependency=afterok:${PREP}"
fi
if stage_at_or_after calibrate; then
  CAL=$(run calibrate ${DEP} --array=0-11 "${S}/05_calibrate.sh")
  DEP="--dependency=afterok:${CAL}"
fi
if stage_at_or_after shard; then
  SH_A=$(run shard.maigan ${DEP} --account="${ACCT_A}" --qos="${QOS_A}" \
           --array=0-$((SPLIT-1))%"${CONC_A}" "${S}/20_shard.sh")
  SH_B=$(run shard.clive  ${DEP} --account="${ACCT_B}" --qos="${QOS_B}" \
           --array=${SPLIT}-$((N_SHARD-1))%"${CONC_B}" "${S}/20_shard.sh")
  DEP="--dependency=afterok:${SH_A}:${SH_B}"
fi
if stage_at_or_after merge; then
  MERGE=$(run merge ${DEP} --array=0-$((N_PAIR-1))%4 "${S}/30_merge.sh")
  DEP="--dependency=afterok:${MERGE}"
fi
if stage_at_or_after call; then
  CALL=$(run call ${DEP} --array=0-$((N_PAIR-1))%8 "${S}/40_call.sh")
  DEP="--dependency=afterok:${CALL}"
fi
if stage_at_or_after assemble; then
  run assemble ${DEP} "${S}/50_assemble.sh" >/dev/null
fi

cat <<NOTES

--------------------------------------------------------------------------------
WATCHING IT
  squeue -u "\$USER" -o '%.10i %.12j %.4t %.10M %.6D %R'
  sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS,ReqMem,NodeList
  tail -f ${XCT_ROOT}/logs/xct_shard_*_0.out
  python ${XCT_HPG}/verify_run.py --run ${XCT_RUN}      # anytime; exit 2 = not done

CONFIRM THE QoS LIMITS BEFORE TRUSTING --time=72:00:00
  sacctmgr show qos ${ACCT_A},${QOS_A},${ACCT_B},${QOS_B} \\
    format=Name,MaxWall,MaxTRESPU,GrpTRES,Priority,Preempt,Flags -P
  A burst QoS at UF is the lower-priority, opportunistic tier and carries a
  SHORTER MaxWall than the investment QoS. If MaxWall on ${QOS_A} is below the
  --time in 20_shard.sh, sbatch rejects the job at submit time rather than at
  hour 23 — so submit the dry run first and read the error. If Preempt is set,
  a shard can be killed mid-flight; that costs one shard, because kodama_shard.R
  skips any shard whose .json already says complete with a matching digest.

RESUMING AFTER ANYTHING
  bash submit_all.sh --from shard        # same XCT_RUN_ID; finished work is skipped
  Nothing needs to be cleaned up first. A dead shard leaves no .json, so it is
  recomputed; a dead merge leaves its finished row blocks, so it resumes there.

THE ONE PATCH WORTH MAKING (optional, high value)
  Each shard recomputes the feature kNN and the dissimilarity reweighting and
  then throws both away — kodama_merge.R redoes them once over the full M. That
  is the shard's memory peak and a large part of its wall clock. KODAMA is your
  own fork (smith6jt-cop/KODAMA), so the fix is a few lines in it: an argument
  that returns right after the M loop,

      KODAMA.matrix(..., res_only = FALSE)
      ...
      if (res_only) return(list(acc = accu, v = vect_acc, res = res,
                                res_constrain = res_constrain, n.cores = n.cores))
      # (existing dissimilarity code follows unchanged)

  inserted immediately before `one_dissimilarity_row` is defined. It changes no
  result — check_shard_equivalence.R still compares the shards' res against a
  monolithic run — and it would cut shard memory from ~19 GB to ~3 GB, which is
  what currently limits the run to ~19 concurrent tasks across both accounts.
--------------------------------------------------------------------------------
NOTES
