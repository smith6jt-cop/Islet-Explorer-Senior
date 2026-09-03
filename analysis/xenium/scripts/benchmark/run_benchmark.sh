#!/usr/bin/env bash
# =============================================================================
# run_benchmark.sh — driver for the spatial-method benchmark (mirrors
# scripts/run_pipeline.sh conventions). Each method runs as its OWN
# process in its OWN env, wrapped by profile_wrap.sh for wall-time + peak RSS.
#
#   prep      [xenium_analysis]  build shared inputs (h5ad + MTX bridge)
#   banksy    [bench_py]         BANKSY domains
#   utag      [bench_py]         UTAG domains
#   baseline  [xenium_analysis]  stage-03 niche (5th domain method)
#   kodama    [xenium R+R_libs]  KODAMA domains
#   misty     [xenium R+R_libs]  MISTy interaction model
#   eval      [xenium_analysis]  scorecard + report + figures
#
# Env control (all optional):
#   SAMPLES="0041323 0041435"   VARIANTS="crop50k"   STEPS="prep methods eval"
#   GENE_PANEL=embedding_panel  NCORES=8  BENCH_TIMEOUT=21600
#   SWEEP_VARIANTS="crop25k crop50k crop100k"   (eval scalability curves)
# STEPS tokens: prep | banksy | utag | baseline | kodama | misty | methods | eval
#   ("methods" = banksy utag baseline kodama misty)
#
# Examples:
#   STEPS="prep methods eval" bash scripts/benchmark/run_benchmark.sh
#   SAMPLES=0041435 VARIANTS=crop25k STEPS=eval bash scripts/benchmark/run_benchmark.sh
# =============================================================================
set -uo pipefail

ROOT="${XENIUM_ANALYSIS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "${ROOT}"
CONDA_SH="${CONDA_SH:-${CONDA_BASE:-/blue/maigan/smith6jt/miniforge3}/etc/profile.d/conda.sh}"

SAMPLES="${SAMPLES:-0041323 0041435}"
VARIANTS="${VARIANTS:-crop50k}"
STEPS="${STEPS:-prep methods eval}"
GENE_PANEL="${GENE_PANEL:-embedding_panel}"
NCORES="${NCORES:-${SLURM_CPUS_PER_TASK:-8}}"
export R_LIBS_USER="${R_LIBS_USER:-/blue/maigan/smith6jt/R_libs}"
SWEEP_VARIANTS="${SWEEP_VARIANTS:-}"
export BENCH_TIMEOUT="${BENCH_TIMEOUT:-}"

PROF="bash ${ROOT}/scripts/benchmark/profile_wrap.sh"
M="${ROOT}/scripts/benchmark/methods"
want() { [[ " ${STEPS} " == *" $1 "* ]] || [[ " ${STEPS} " == *" methods "* && " banksy utag baseline kodama misty " == *" $1 "* ]]; }

echo "benchmark: SAMPLES=[${SAMPLES}] VARIANTS=[${VARIANTS}] STEPS=[${STEPS}] NCORES=${NCORES} PANEL=${GENE_PANEL}"

# ---- prep (xenium_analysis) ----
if want prep; then
  ( source "${CONDA_SH}"; conda activate xenium_analysis
    for s in ${SAMPLES}; do
      python "${ROOT}/scripts/benchmark/prep_benchmark_input.py" \
        --sample "${s}" --variants ${VARIANTS} --gene-panel "${GENE_PANEL}"
    done )
fi

# ---- Python domain methods: BANKSY + UTAG (bench_py + LD_PRELOAD) ----
if want banksy || want utag; then
  ( source "${ROOT}/scripts/benchmark/env_bench_py.sh"
    for s in ${SAMPLES}; do for v in ${VARIANTS}; do
      want banksy && ${PROF} "${s}" "${v}" banksy -- python "${M}/run_banksy.py" --sample "${s}" --variant "${v}"
      want utag   && ${PROF} "${s}" "${v}" utag   -- python "${M}/run_utag.py"   --sample "${s}" --variant "${v}"
    done; done )
fi

# ---- baseline niche (xenium_analysis) ----
if want baseline; then
  ( source "${CONDA_SH}"; conda activate xenium_analysis
    for s in ${SAMPLES}; do for v in ${VARIANTS}; do
      ${PROF} "${s}" "${v}" baseline -- python "${M}/run_baseline_niche.py" --sample "${s}" --variant "${v}"
    done; done )
fi

# ---- R methods: KODAMA + MISTy (xenium_analysis R + R_libs) ----
if want kodama || want misty; then
  ( source "${CONDA_SH}"; conda activate xenium_analysis
    for s in ${SAMPLES}; do for v in ${VARIANTS}; do
      want kodama && ${PROF} "${s}" "${v}" kodama -- Rscript "${M}/run_kodama.R" --sample "${s}" --variant "${v}" --n-cores "${NCORES}"
      want misty  && ${PROF} "${s}" "${v}" misty  -- Rscript "${M}/run_misty.R"  --sample "${s}" --variant "${v}" --n-cores "${NCORES}"
    done; done )
fi

# ---- eval (xenium_analysis) ----
if want eval; then
  ( source "${CONDA_SH}"; conda activate xenium_analysis
    for s in ${SAMPLES}; do for v in ${VARIANTS}; do
      sweep_arg=""; [[ -n "${SWEEP_VARIANTS}" ]] && sweep_arg="--sweep-variants ${SWEEP_VARIANTS}"
      python "${ROOT}/scripts/benchmark/evaluate.py" --sample "${s}" --variant "${v}" ${sweep_arg}
    done; done )
fi

echo "run_benchmark.sh complete."
