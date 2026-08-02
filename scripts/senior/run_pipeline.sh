#!/usr/bin/env bash
# Thin convenience wrapper for the manifest-driven Phenocycler workflow.
#
# Examples:
#   bash scripts/senior/run_pipeline.sh
#   bash scripts/senior/run_pipeline.sh run --jobs 8
#   bash scripts/senior/run_pipeline.sh status
#   bash scripts/senior/run_pipeline.sh export
#
# Set PHENOCYCLER_CONFIG to use a configuration other than the repository
# config.ini. Raw source paths belong in the cohort manifest, not environment
# variables consumed by this wrapper.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONFIG="${PHENOCYCLER_CONFIG:-${REPO}/config.ini}"

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

cd "$REPO"
exec python -m phenocycler.pipeline "$COMMAND" --config "$CONFIG" "$@"
