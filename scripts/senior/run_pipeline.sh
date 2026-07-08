#!/usr/bin/env bash
# Run the Phenocycler_Analysis phenotyping pipeline against Islet-Explorer-Senior's data.
#
# The pipeline engine lives in the Phenocycler_Analysis git submodule (canonical steps 1-5:
# cells -> REDSEA -> RESTORE(+extra) -> hormone_floor -> broad lineage -> qupath -> figures).
# Install it once into the pipeline env:
#     pip install -e Phenocycler_Analysis        # inside `conda activate scvi-env`
# then invoke stages, e.g.:
#     bash scripts/senior/run_pipeline.sh --status
#     bash scripts/senior/run_pipeline.sh --only hormone_floor lineage figures --jobs 8
#
# Machine-specific inputs are passed as env-var overrides so config.ini stays portable
# (e.g. for the HiPerGator /blue/maigan paths just export the three vars before calling).
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export PHENOCYCLER_IMAGES_DIR="${PHENOCYCLER_IMAGES_DIR:-/home/smith6jt/IO60panc2nd/Images}"
export PHENOCYCLER_DONOR_METADATA="${PHENOCYCLER_DONOR_METADATA:-/home/smith6jt/IO60panc2nd/donor_metadata_panc.xlsx}"
export PHENOCYCLER_CELLS_CSV="${PHENOCYCLER_CELLS_CSV:-/home/smith6jt/IO60panc2nd/Cellmeasurements.csv}"
cd "$REPO"
exec python -m phenocycler.pipeline --config "$REPO/config.ini" "$@"
