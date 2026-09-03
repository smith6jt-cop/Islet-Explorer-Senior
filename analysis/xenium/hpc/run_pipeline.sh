#!/bin/bash
# Consolidated Xenium pipeline driver — raw → final, one TISSUE-aware entry point.
#
# Supersedes the three older drivers (run_local_pipeline.sh, run_identity_core_redesign.sh,
# run_thymus_post_processing.sh). Runs inside the allocated HiPerGator VSCode tunnel
# session (a single B200 GPU). DO NOT sbatch this script.
#
# Stage vocabulary (STAGES env var or positional $1):
#   01         per-sample  nbconvert 01_preprocessing_v2.ipynb   → {s}_preprocessed.h5ad
#   02         per-sample  nbconvert 02_phenotyping.ipynb        → {s}_phenotyped.h5ad (raw labels)
#   phenotype  whole-tissue [backup] finalize_panel → hierarchical_phenotype (--output inplace)
#                          → validate_redesign (GATE) → render_stage02_phenotyping
#                          ⇒ {s}_phenotyped.h5ad gets the identity_core labels
#   03         per-sample  nbconvert 03_spatial_analysis.ipynb   → {s}_spatial_analysis.h5ad + CSVs
#   final      tissue-disp pancreas: nb07 (07_immune_islet_proximity) → insulitis_analysis.py
#                          thymus:   thymic_zonation_analysis.py
#   explorer   per-sample  export_xenium_explorer_groups.py                        (opt-in)
#   05         batch       nbconvert 05_tissue_comparisons.ipynb                    (opt-in)
#
# IDEMPOTENCY: stage 02 writes its own labels into {s}_phenotyped.h5ad, so any `02`
# is auto-followed by `phenotype` (which re-stamps the identity_core labels). Set
# XENIUM_SKIP_MERGE=1 to suppress that coupling (e.g. to inspect raw notebook labels).
#
# Usage:
#   bash scripts/run_pipeline.sh                                   # default stages, pancreas
#   STAGES="03 final" TISSUE=thymus bash scripts/run_pipeline.sh   # regen thymus finals
#   STAGES="01 02 03 final" TISSUE=pancreas bash scripts/run_pipeline.sh   # full from-scratch
#   bash scripts/run_pipeline.sh 02 0041323                        # one stage, one sample
#
# Logs:  logs/local_<stage>_<sample>_<ts>.{out,err}; executed nbs under notebooks/executed/.

set -uo pipefail

PROJECT_ROOT="${XENIUM_ANALYSIS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$PROJECT_ROOT"

source "${CONDA_BASE:-/blue/maigan/smith6jt/miniforge3}/etc/profile.d/conda.sh"
conda activate xenium_analysis

NCPU=$(nproc)
export OMP_NUM_THREADS="$NCPU"
export MKL_NUM_THREADS="$NCPU"
export NUMEXPR_MAX_THREADS="$NCPU"
export OPENBLAS_NUM_THREADS="$NCPU"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}"

# CUDA toolkit env for rapids_singlecell / cuml / cupy / cudf (CLAUDE.md §7).
source "${PROJECT_ROOT}/scripts/env_rapids.sh"

# ---- TISSUE switch: per-tissue sample defaults, QC floors, finalize audit CSV ----
TISSUE="${TISSUE:-pancreas}"
case "$TISSUE" in
    pancreas)
        SAMPLES_DEFAULT="0041323 0041326"
        export XENIUM_MIN_COUNTS="${XENIUM_MIN_COUNTS:-50}"
        export XENIUM_MIN_GENES="${XENIUM_MIN_GENES:-20}"
        AUDIT_ARG=()                         # finalize_panel.py default audit CSV
        ;;
    thymus)
        SAMPLES_DEFAULT="0041435 0041534"
        # Thymus median ~107-124 counts/cell — pancreas MIN_COUNTS=50 would cut
        # ~20-25% of legitimate small thymocytes. Lower floors retain ≥90% of cells.
        export XENIUM_MIN_COUNTS="${XENIUM_MIN_COUNTS:-20}"
        export XENIUM_MIN_GENES="${XENIUM_MIN_GENES:-10}"
        AUDIT_ARG=(--audit-csv data/processed/panel_audit_thymus.csv)
        ;;
    *) echo "unknown TISSUE=$TISSUE (expected pancreas|thymus)" >&2; exit 2 ;;
esac
SAMPLES="${SAMPLES:-$SAMPLES_DEFAULT}"

# Positional override: $1 = single stage, $2 = single sample.
if [ $# -ge 1 ]; then STAGES="$1"; fi
if [ $# -ge 2 ]; then SAMPLES="$2"; fi
STAGES="${STAGES:-01 02 phenotype 03}"

mkdir -p logs

# ---- Auto-couple 02 → phenotype unless explicitly skipped ----
EFFECTIVE_STAGES=""
for s in $STAGES; do
    EFFECTIVE_STAGES="$EFFECTIVE_STAGES $s"
    if [ "$s" = "02" ] && [ "${XENIUM_SKIP_MERGE:-0}" != "1" ]; then
        case " $STAGES " in
            *" phenotype "*) : ;;                                  # user already asked for it
            *) EFFECTIVE_STAGES="$EFFECTIVE_STAGES phenotype"
               echo "[driver] note: auto-inserting 'phenotype' after '02' (set XENIUM_SKIP_MERGE=1 to disable)";;
        esac
    fi
done

echo "[$(date)] driver config: TISSUE=$TISSUE SAMPLES=[$SAMPLES] STAGES=[${EFFECTIVE_STAGES# }]"

# ---- helpers ----------------------------------------------------------------
run () {  # run a python/CLI step with a banner + rc echo; returns the step's rc
    echo; echo "######## $* ########"
    "$@"; local rc=$?
    echo "## rc=$rc"
    [ $rc -ne 0 ] && echo "!! FAILED: $*"
    return $rc
}

run_per_sample () {  # nbconvert a per-sample notebook (01/02/03)
    local stage="$1" sample="$2" nb
    case "$stage" in
        01) nb="notebooks/01_preprocessing_v2.ipynb" ;;
        02) nb="notebooks/02_phenotyping.ipynb" ;;
        03) nb="notebooks/03_spatial_analysis.ipynb" ;;
        *)  echo "unknown per-sample stage: $stage" >&2; return 2 ;;
    esac
    local ts; ts=$(date +%Y%m%d_%H%M%S)
    local outdir="notebooks/executed/${sample}"
    mkdir -p "$outdir" "figures/${sample}/01_preprocessing" \
             "figures/${sample}/02_phenotyping" "figures/${sample}/03_spatial_analysis"
    local out="logs/local_${stage}_${sample}_${ts}.out"
    local err="logs/local_${stage}_${sample}_${ts}.err"
    echo "[$(date)] START stage=$stage sample=$sample  nb=$nb"
    echo "  logs: $out / $err"
    XENIUM_SAMPLE_ID="$sample" XENIUM_TISSUE="$TISSUE" \
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=86400 \
        --output-dir="$outdir" --output="${stage}_${sample}.ipynb" \
        "$nb" >"$out" 2>"$err"
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "[$(date)] FAIL stage=$stage sample=$sample (rc=$rc); tail $err:"; tail -30 "$err" || true
        return "$rc"
    fi
    echo "[$(date)] OK   stage=$stage sample=$sample"
}

run_batch_notebook () {  # nbconvert a whole-tissue batch notebook (nb07, 05)
    local nb="$1" tag="$2"
    local ts; ts=$(date +%Y%m%d_%H%M%S)
    local outdir="notebooks/executed/${tag}"; mkdir -p "$outdir"
    local out="logs/local_${tag}_${ts}.out" err="logs/local_${tag}_${ts}.err"
    echo "[$(date)] START batch nb=$nb tag=$tag"
    echo "  logs: $out / $err"
    XENIUM_TISSUE="$TISSUE" \
    XENIUM_DATA_ROOT="$PROJECT_ROOT/data/processed" \
    XENIUM_FIG_ROOT="$PROJECT_ROOT/figures" \
    XENIUM_SAMPLE_IDS="$(echo "$SAMPLES" | tr ' ' ',')" \
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=86400 \
        --output-dir="$outdir" --output="$(basename "${nb%.ipynb}")_${tag}.ipynb" \
        "$nb" >"$out" 2>"$err"
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "[$(date)] FAIL batch nb=$nb (rc=$rc); tail $err:"; tail -30 "$err" || true
        return "$rc"
    fi
    echo "[$(date)] OK   batch nb=$nb"
}

stage_phenotype () {  # identity_core merge: backup → finalize → hierarchical → validate → figs
    local ts; ts=$(date +%Y%m%d_%H%M%S)
    local bk="data/processed/_archive_${ts}_premerge"
    echo "[$(date)] PHENOTYPE: mandatory pre-merge backup → $bk"
    mkdir -p "$bk" || { echo "!! backup mkdir failed"; return 1; }
    for s in $SAMPLES; do
        local src="data/processed/$s/${s}_phenotyped.h5ad"
        [ -f "$src" ] || { echo "!! missing $src (run stage 02 first)"; return 1; }
        cp -f "$src" "$bk/${s}_phenotyped.h5ad" || { echo "!! BACKUP FAILED for $s — aborting before in-place merge"; return 1; }
    done
    run python scripts/finalize_panel.py --tissue "$TISSUE" --samples $SAMPLES "${AUDIT_ARG[@]+"${AUDIT_ARG[@]}"}" || return 1
    run python scripts/hierarchical_phenotype.py --tissue "$TISSUE" --samples $SAMPLES \
        --output inplace --dl-workers 0 || return 1
    run python scripts/validate_redesign.py --tissue "$TISSUE" --samples $SAMPLES || return 1   # GATE
    run python scripts/render_stage02_phenotyping.py --tissue "$TISSUE" $SAMPLES || return 1
}

stage_final () {  # tissue-dispatched terminal analysis
    case "$TISSUE" in
        pancreas)
            run_batch_notebook notebooks/07_immune_islet_proximity.ipynb "07_pancreas" || return 1
            for s in $SAMPLES; do
                [ -f "data/processed/$s/${s}_immune_phenotyped.h5ad" ] || {
                    echo "!! nb07 did not produce ${s}_immune_phenotyped.h5ad — aborting insulitis"; return 1; }
            done
            run python scripts/insulitis_analysis.py --samples $SAMPLES || return 1
            ;;
        thymus)
            run python scripts/thymic_zonation_analysis.py --samples $SAMPLES || return 1
            ;;
    esac
}

stage_explorer () {  # opt-in Xenium Explorer export
    for s in $SAMPLES; do
        run python scripts/export_xenium_explorer_groups.py \
            -i "data/processed/$s/${s}_phenotyped.h5ad" || return 1
    done
}

run_integration () {  # opt-in cross-tissue (stage 05)
    local csv; csv=$(echo "$SAMPLES" | tr ' ' ',')
    local tag; tag=$(echo "$csv" | tr ',' '_')
    local ts; ts=$(date +%Y%m%d_%H%M%S)
    local outdir="notebooks/executed/integrated_${tag}"
    mkdir -p "$outdir" "figures/integrated_${tag}" "data/processed/integrated"
    local out="logs/local_05_${tag}_${ts}.out" err="logs/local_05_${tag}_${ts}.err"
    echo "[$(date)] START stage=05 samples=$csv"
    XENIUM_SAMPLE_IDS="$csv" XENIUM_FIGURES_DIR="figures/integrated_${tag}" XENIUM_TISSUE="$TISSUE" \
    jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=86400 \
        --output-dir="$outdir" --output="05_tissue_comparisons_${tag}.ipynb" \
        notebooks/05_tissue_comparisons.ipynb >"$out" 2>"$err"
    local rc=$?
    [ "$rc" -ne 0 ] && { echo "[$(date)] FAIL stage=05 (rc=$rc)"; tail -30 "$err" || true; return "$rc"; }
    echo "[$(date)] OK   stage=05"
}

# ---- dispatch ---------------------------------------------------------------
for stage in $EFFECTIVE_STAGES; do
    case "$stage" in
        01|02|03)  for sample in $SAMPLES; do run_per_sample "$stage" "$sample" || exit 1; done ;;
        phenotype) stage_phenotype || exit 1 ;;
        final)     stage_final     || exit 1 ;;
        explorer)  stage_explorer  || exit 1 ;;
        05)        run_integration || exit 1 ;;
        *)         echo "Unknown stage: $stage (expected 01 02 phenotype 03 final explorer 05)" >&2; exit 2 ;;
    esac
done

echo "[$(date)] Pipeline finished. TISSUE=$TISSUE stages=[${EFFECTIVE_STAGES# }] samples=[$SAMPLES]"
