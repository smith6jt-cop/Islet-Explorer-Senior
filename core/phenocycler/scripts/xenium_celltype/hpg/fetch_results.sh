#!/usr/bin/env bash
# =============================================================================
# fetch_results.sh — bring the run's RESULTS back to the workstation, and
# nothing else.
#
# WHAT COMES BACK, AND WHAT MUST NOT. The results are small: per (donor, arm) an
# embedding (~30 MB at 1.1M cells), group calls and a composition (a few hundred
# KB), plus the manifests and the verification. About 1.5 GB for the cohort.
# The inputs are not: pcs.mtx is ~500 MB per donor as text and the cached kNN is
# ~9 GB per donor. Those are reproducible from the manifest and stay on /blue.
# The exclude list below is what keeps a 200 GB accident from happening.
#
# RUN THIS FROM HIPERGATOR, pushing to the workstation, if the workstation is
# reachable; otherwise run it ON the workstation pulling, which is the normal
# direction and the one that works today (see the ssh notes in submit_all.sh).
# Either way it is rsync over ssh with --partial, so an interrupted transfer
# resumes rather than restarting.
#
#   # from the workstation, pulling:
#   bash fetch_results.sh --pull hpg:/blue/maigan/smith6jt/Xenium_Celltype/runs/<id>
#   # or via Globus, which is already set up for this data (see config.ini's
#   # xenium_root note): the same source and destination paths apply.
# =============================================================================
set -uo pipefail

MODE="pull"; SRC=""; DEST="${XCT_LOCAL_DEST:-$HOME/IO60panc2nd/Islet-Explorer-Senior/data/development/xenium_celltype}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pull) MODE="pull"; SRC="$2"; shift 2 ;;
    --push) MODE="push"; SRC="$2"; shift 2 ;;
    --dest) DEST="$2"; shift 2 ;;
    *) echo "unknown arg: $1"; exit 2 ;;
  esac
done
[[ -z "${SRC}" ]] && { echo "usage: fetch_results.sh --pull <remote run dir> [--dest DIR]"; exit 2; }

# Results live under the repo-root data/development/, a level above core/, per CLAUDE.md:
# mutable review artifacts never write a stage manifest and never touch data/runs/.
mkdir -p "${DEST}"

INCLUDE=(
  --include='work_manifest.json' --include='verification.json'
  --include='cohort_composition.csv' --include='cohort_arm_difference.csv'
  --include='*/' 
  --include='embedding_*.csv' --include='kodama_*.json'
  --include='group_calls_*.csv' --include='composition_*.csv' --include='call_*.json'
  --include='prep.json' --include='res_*_shard*.json'
  --include='calibration/**'
)
EXCLUDE=(
  --exclude='pcs.mtx' --exclude='knn_*.rds' --exclude='res_*_shard*.rds'
  --exclude='embedding_*.rds' --exclude='dissim_*/' --exclude='coords.csv'
  --exclude='cell_id.csv' --exclude='*'
)

if [[ "${MODE}" == "pull" ]]; then
  echo "[fetch] ${SRC} -> ${DEST}"
  rsync -avh --partial --info=progress2 "${INCLUDE[@]}" "${EXCLUDE[@]}" "${SRC}/" "${DEST}/"
else
  echo "[fetch] ${SRC} -> ${DEST} (push)"
  rsync -avh --partial --info=progress2 "${INCLUDE[@]}" "${EXCLUDE[@]}" "${SRC}/" "${DEST}/"
fi

echo
echo "[fetch] what arrived:"
du -sh "${DEST}" 2>/dev/null
find "${DEST}" -name 'composition_*.csv' | wc -l | xargs echo "  composition files:"
find "${DEST}" -name 'embedding_*.csv' | wc -l | xargs echo "  embeddings:"
if [[ -f "${DEST}/verification.json" ]]; then
  python - "${DEST}/verification.json" <<'PY'
import json, sys
v = json.load(open(sys.argv[1]))
print(f"  verification: {'COMPLETE' if v['complete'] else 'INCOMPLETE'}"
      f" ({len(v['problems'])} problem(s))")
for p in v["problems"][:10]:
    print(f"    {p}")
PY
else
  echo "  verification.json ABSENT — the assemble stage did not run; this is not a finished run"
fi
