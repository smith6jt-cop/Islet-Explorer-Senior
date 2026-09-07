#!/usr/bin/env bash
# Integration stages for the two representative donors once BOTH development releases exist
# (plan 2026-09-05, Phase 3.2).  Stage by stage, fail-fast, each stage logged; the verdict-gated stages
# (match, cellmatch, grid's registered half, crossmodal) are skipped by the pipeline itself where
# qc.ALLOWED_BY_VERDICT does not permit them.  Run from core/phenocycler.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=../../.venv/bin/python
LOG=../../data/development/integration_logs
mkdir -p "$LOG"
DONORS=(--donor 6436 --donor 6442)
run() {  # stage [extra flags]
  local stage=$1; shift
  echo "== $(date '+%F %T') stage $stage $*" | tee -a "$LOG/pair.log"
  $PY -m phenocycler.integration.pipeline --only "$stage" --tissue pancreas "${DONORS[@]}" "$@" 2>&1 | grep -v "^\[phenospace\] .* does not exist yet" | tee "$LOG/${stage}.log" | tail -3
}
$PY -c "from phenocycler.config import load_config; from phenocycler.integration.phenospace import identity_attached; import sys; sys.exit(0 if identity_attached(load_config('config.ini')) else 2)" 2>/dev/null || { echo "identity not attached -- both development releases must load through phenospace first"; exit 2; }
run import_xenium --force
run postxen --force
run vocab --force
run structures --force
run register --force
run qc --force
run donor --force
run grid --force
run match --force
run cellmatch --force
run crossmodal --force
run qc --force
run figures --force
$PY -m phenocycler.integration.pipeline --status 2>&1 | tail -20 | tee -a "$LOG/pair.log"
echo "== $(date '+%F %T') integration pair done" | tee -a "$LOG/pair.log"
