#!/usr/bin/env bash
# Thin convenience wrapper for the manifest-driven PhenoCycler source workflow.
#
#   QuPath cohort manifest -> ingest -> geometry QC -> REDSEA -> expression
#
# This wrapper runs ONLY that source workflow. The layer and typing/app artifacts built on top of
# it have their own driver:  scripts/senior/build_layers.py --check | --plan | --run
#
# Examples:
#   bash scripts/senior/run_pipeline.sh
#   bash scripts/senior/run_pipeline.sh run --pipelined --through expression --jobs 8
#   bash scripts/senior/run_pipeline.sh status
#
# Two things this wrapper exists to get right, both of which it got wrong until 2026-09-18:
#
#   1. THE CONFIG. core/phenocycler/ owns the only config.ini. The default here pointed at a
#      repo-root ${REPO}/config.ini, which does not exist; root CLAUDE.md records that one DID
#      exist until 2026-08-10 and silently supplied stale defaults, so a repo-root config.ini is
#      a bug rather than a fallback. Set PHENOCYCLER_CONFIG to use a different configuration.
#   2. THE INTERPRETER. `python` on PATH is not the repo environment. The .venv is used when it
#      is present, so `phenocycler` and its dependencies resolve to the ones this checkout pins.
#
# `export` and the stages after `expression` are blocked by NO_VERIFIED_TYPING; the latch has not
# moved and this wrapper does not touch it. `export` is still accepted so the refusal comes from
# the pipeline, which names the latch, rather than from an argument check here.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PACKAGE_DIR="${REPO}/core/phenocycler"
CONFIG="${PHENOCYCLER_CONFIG:-${PACKAGE_DIR}/config.ini}"
PYTHON="${REPO}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || PYTHON="python"

if [[ ! -f "${CONFIG}" ]]; then
    echo "$0: no configuration at ${CONFIG}" >&2
    echo "  core/phenocycler/ owns the only config.ini; set PHENOCYCLER_CONFIG to use another." >&2
    exit 2
fi

if [[ $# -eq 0 ]]; then
    set -- run
fi

COMMAND="$1"
shift
case "$COMMAND" in
    run|status|export)
        ;;
    *)
        echo "usage: $0 [run|status|export] [pipeline arguments...]" >&2
        exit 2
        ;;
esac

# The workflow runs FROM core/phenocycler, which owns the config and the data roots it names.
cd "${PACKAGE_DIR}"
exec "${PYTHON}" -m phenocycler.pipeline "$COMMAND" --config "$CONFIG" "$@"
