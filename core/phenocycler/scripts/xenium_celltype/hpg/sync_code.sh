#!/usr/bin/env bash
# =============================================================================
# sync_code.sh — put this repository's Xenium cell-typing lane on HiPerGator and
# keep it honest about WHICH commit is running.
#
# git clone, not rsync. rsync makes it impossible to answer "which version
# produced this run?" six months later, and the run manifest embeds the commit
# so that question has an answer. The repo is private, so the clone uses the
# maintainer's existing GitHub access from a HiPerGator login node -- that also
# removes any need for this workstation to reach HiPerGator at all (see the
# ssh notes at the bottom of submit_all.sh).
#
# The lane is small: core/phenocycler/scripts/xenium_celltype plus the
# phenocycler package it imports for config + panel taxonomy. A sparse checkout
# keeps 158 GB of gitignored run data and the Shiny app off /blue.
#
# Run this ON HiPerGator:
#   bash sync_code.sh                 # clone or fast-forward
#   BRANCH=feat/xenium-celltype bash sync_code.sh
# =============================================================================
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/env.sh"

REPO="${REPO:-git@github.com:smith6jt-cop/Islet-Explorer-Senior.git}"
BRANCH="${BRANCH:-main}"

mkdir -p "${XCT_ROOT}" "${XCT_RUNS}" "${XCT_ROOT}/logs"

if [[ ! -d "${XCT_CODE}/.git" ]]; then
  echo "[sync] cloning ${REPO} (${BRANCH}) -> ${XCT_CODE}"
  git clone --filter=blob:none --no-checkout --branch "${BRANCH}" "${REPO}" "${XCT_CODE}" \
    || { echo "[sync] clone FAILED — is the SSH key for github registered on this node?"; exit 1; }
  git -C "${XCT_CODE}" sparse-checkout init --cone
  git -C "${XCT_CODE}" sparse-checkout set core/phenocycler
  git -C "${XCT_CODE}" checkout "${BRANCH}"
else
  echo "[sync] updating ${XCT_CODE}"
  git -C "${XCT_CODE}" fetch --prune origin "${BRANCH}" || exit 1
  if [[ -n "$(git -C "${XCT_CODE}" status --porcelain)" ]]; then
    echo "[sync] REFUSING to update: the checkout has local modifications."
    git -C "${XCT_CODE}" status --short
    echo "[sync] A run whose code differs from its recorded commit is unreproducible."
    exit 1
  fi
  git -C "${XCT_CODE}" checkout "${BRANCH}" && git -C "${XCT_CODE}" merge --ff-only "origin/${BRANCH}" || exit 1
fi

SHA="$(git -C "${XCT_CODE}" rev-parse HEAD)"
echo "[sync] ${BRANCH} @ ${SHA}"
echo "[sync] scripts: ${XCT_LANE}"
ls -1 "${XCT_LANE}" "${XCT_HPG}" 2>/dev/null | sed 's/^/    /'

# The pinned R package is NOT in this repo; it comes from setup_envs.sh and is recorded in
# kodama.lock, which the gate prints and verify_run.py copies into the run directory.
echo
echo "[sync] next:"
echo "    bash ${XCT_HPG}/setup_envs.sh"
echo "    bash ${XCT_HPG}/submit_all.sh --dry-run"
