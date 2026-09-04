#!/usr/bin/env bash
# Execute the hypothesis notebooks in dependency order.
#
# The tiers are not cosmetic. H10 builds the feature matrix H8/H9/H19 read;
# H15/H16 write the cell-type tables H20 reads; H18/H19 pool everything. Run a
# later tier before an earlier one and it will quietly consume the previous
# round's numbers.
#
#   scripts/run_notebooks.sh            # every tier, in order
#   scripts/run_notebooks.sh 1 2        # only tiers 1 and 2
#
# Each notebook is executed in place, so outputs land in the file. Failures are
# recorded and the run continues, because one broken notebook should not hide
# the state of the other nineteen; the exit status is non-zero if any failed.
set -uo pipefail

PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$PROJECT/../.." && pwd)"
DEFAULT_PY="$REPO_ROOT/.venv/bin/python"
if [ ! -x "$DEFAULT_PY" ]; then
    DEFAULT_PY="$(command -v python)"
fi
PY="${FOLLICLE_PYTHON:-$DEFAULT_PY}"
cd "$PROJECT/analysis" || exit 1

# This shell exports MPLBACKEND=Agg. Under Agg, plt.show() is a no-op inside a
# Jupyter kernel: every notebook still writes its PNGs through save_figure, but
# not one figure is embedded in the .ipynb — which is silent, since nothing
# errors and the files on disk look right. The notebooks are tracked in git and
# are meant to be readable there, so force the inline backend.
export MPLBACKEND='module://matplotlib_inline.backend_inline'

TIER1="H1_vessel_density H2_vessel_morphology H3_follicle_structure
       H4_tissue_proportions H5_follicle_vascularization
       H6_follicle_reconstruction H7_snp_screen"
TIER2="H10_multivariate H8_clinical_covariates H9_polygenic_risk"
TIER3="H12_spatial_bootstrap H13_cropped_replication"
TIER4="H11_lda_deconfounded H14_cell_boundaries H15_single_cell_analysis
       H16_cell_phenotyping H17_perivascular_phenotyping
       H20_bcell_density_feasibility"
TIER5="H18_HE_follicles H19_meta_analysis"

# Args are tier numbers, or notebook stems to run exactly (in the order given).
# The latter is for resuming after a fix without repeating the expensive
# notebooks that already succeeded — H11 alone is 33 minutes.
WANT=("$@")
[ ${#WANT[@]} -eq 0 ] && WANT=(1 2 3 4 5)

if [[ "${WANT[0]}" == H* ]]; then
    TIER_EXPLICIT="${WANT[*]}"
    WANT=(_EXPLICIT)
fi

failed=()
for tier in "${WANT[@]}"; do
    eval "nbs=\$TIER$tier"
    echo "=============== TIER $tier ==============="
    for nb in $nbs; do
        printf '%-42s ' "$nb"
        start=$SECONDS
        if "$PY" -m jupyter nbconvert --to notebook --execute --inplace \
               --ExecutePreprocessor.timeout=-1 \
               "$nb.ipynb" >"/tmp/nb_$nb.log" 2>&1; then
            echo "ok   ($((SECONDS-start))s)"
        else
            echo "FAIL ($((SECONDS-start))s)  -> /tmp/nb_$nb.log"
            failed+=("$nb")
        fi
    done
done

echo
if [ ${#failed[@]} -eq 0 ]; then
    echo "all notebooks executed cleanly"
else
    echo "FAILED (${#failed[@]}): ${failed[*]}"
    exit 1
fi
