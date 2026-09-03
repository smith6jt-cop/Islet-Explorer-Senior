#!/usr/bin/env bash
# =============================================================================
# setup_envs.sh — build the ISOLATED environments for the spatial-method
# benchmark (BANKSY, UTAG, KODAMA, MISTy). The production `xenium_analysis`
# env is NEVER modified.
#
#   * Python methods (BANKSY, UTAG) -> new conda env `bench_py`.
#       UTAG hard-pins scanpy>=1.9.3,<1.12, which is INCOMPATIBLE with
#       xenium_analysis (scanpy 1.12 + rapids_singlecell 0.14.1, decision #9),
#       so they must live in their own env.
#   * R methods (KODAMA, MISTy) -> reuse the in-env R 4.5.3 but install ONLY
#       the missing packages into the project-local, already-writable R library
#       `/blue/maigan/smith6jt/R_libs`. An R user-library physically cannot
#       touch Python site-packages, so the production pins stay safe.
#
# Idempotent: re-running skips work already done. Prints versions at the end.
#
# Usage:
#   bash scripts/benchmark/setup_envs.sh            # both
#   bash scripts/benchmark/setup_envs.sh --py-only  # just bench_py
#   bash scripts/benchmark/setup_envs.sh --r-only   # just R packages
# =============================================================================
set -uo pipefail

CONDA_BASE="${CONDA_BASE:-/blue/maigan/smith6jt/miniforge3}"  # HiPerGator default
BENCH_ENV="bench_py"
R_LIBS_USER="${R_LIBS_USER:-/blue/maigan/smith6jt/R_libs}"
export R_LIBS_USER

DO_PY=1; DO_R=1
case "${1:-}" in
  --py-only) DO_R=0 ;;
  --r-only)  DO_PY=0 ;;
  "") ;;
  *) echo "unknown arg: $1"; exit 2 ;;
esac

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# -----------------------------------------------------------------------------
# 1. Python env: bench_py
# -----------------------------------------------------------------------------
if [[ $DO_PY -eq 1 ]]; then
  echo "=============================================================="
  echo "[bench_py] building Python env for BANKSY + UTAG"
  echo "=============================================================="
  if conda env list | grep -qE "^${BENCH_ENV}\s"; then
    echo "[bench_py] env already exists — skipping create"
  else
    # Base from conda-forge (C-extension pkgs are more reliable via conda);
    # scanpy/squidpy/pybanksy/utag installed together via pip so the resolver
    # produces one consistent scanpy<1.12 stack.
    mamba create -y -n "${BENCH_ENV}" -c conda-forge \
        python=3.11 numpy scipy pandas scikit-learn h5py matplotlib \
        leidenalg python-igraph numba || { echo "[bench_py] create FAILED"; exit 1; }
  fi

  conda activate "${BENCH_ENV}"
  echo "[bench_py] python: $(python --version 2>&1) @ $(which python)"

  # pip-resolve the analysis stack in one pass (utag pulls scanpy<1.12 + squidpy).
  # setuptools<81 is REQUIRED: squidpy->spatialdata->xarray_schema imports the
  # deprecated pkg_resources, which setuptools>=81 no longer ships.
  pip install --no-input \
      "setuptools<81" "scanpy>=1.9.3,<1.12" "squidpy>=1.4" pybanksy \
      "git+https://github.com/ElementoLab/utag.git@main" \
      || { echo "[bench_py] pip install FAILED"; exit 1; }

  echo "[bench_py] verifying imports ..."
  python - <<'PY'
import importlib
mods = ["banksy", "utag", "scanpy", "squidpy", "anndata", "leidenalg", "igraph", "sklearn"]
for m in mods:
    mod = importlib.import_module(m)
    print(f"  {m:10s} {getattr(mod,'__version__','?')}")
print("bench_py OK")
PY
  conda deactivate
fi

# -----------------------------------------------------------------------------
# 2. R packages into the project R library (reusing xenium_analysis R 4.5.3)
# -----------------------------------------------------------------------------
if [[ $DO_R -eq 1 ]]; then
  echo "=============================================================="
  echo "[R] installing KODAMA + mistyR into ${R_LIBS_USER}"
  echo "=============================================================="
  conda activate xenium_analysis
  mkdir -p "${R_LIBS_USER}"
  echo "[R] Rscript: $(which Rscript)"

  Rscript --no-save - <<'RS'
options(repos = c(CRAN = "https://cloud.r-project.org"))
lib <- Sys.getenv("R_LIBS_USER")
.libPaths(c(lib, .libPaths()))
cat("R lib target:", lib, "\n")
cat(R.version.string, "\n")

have <- function(p) requireNamespace(p, quietly = TRUE)
need_cran <- c("BiocManager","remotes","ranger","distances","furrr","future",
               "Rtsne","uwot","R6","rlist","deldir","filelock","digest")
to_get <- need_cran[!vapply(need_cran, have, logical(1))]
if (length(to_get)) {
  cat("installing CRAN:", paste(to_get, collapse=", "), "\n")
  install.packages(to_get, lib = lib, Ncpus = max(1, parallel::detectCores() - 1))
} else cat("CRAN deps already present\n")

if (!have("mistyR")) {
  cat("installing mistyR (Bioconductor) ...\n")
  BiocManager::install("mistyR", lib = lib, update = FALSE, ask = FALSE)
} else cat("mistyR present:", as.character(packageVersion("mistyR")), "\n")

if (!have("KODAMA")) {
  cat("installing KODAMA (GitHub tkcaccia/KODAMA) ...\n")
  remotes::install_github("tkcaccia/KODAMA", lib = lib, upgrade = "never")
} else cat("KODAMA present:", as.character(packageVersion("KODAMA")), "\n")
RS

  echo "[R] verifying library loads ..."
  Rscript --no-save -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); \
    suppressMessages({library(mistyR); library(KODAMA); library(ranger); library(future)}); \
    cat("mistyR", as.character(packageVersion("mistyR")), "| KODAMA", as.character(packageVersion("KODAMA")), "\n"); \
    cat("R bench OK\n")' \
    || { echo "[R] verification FAILED"; exit 1; }
  conda deactivate
fi

echo "=============================================================="
echo "setup_envs.sh complete."
echo "=============================================================="
