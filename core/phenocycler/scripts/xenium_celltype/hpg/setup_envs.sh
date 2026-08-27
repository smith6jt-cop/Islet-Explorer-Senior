#!/usr/bin/env bash
# =============================================================================
# setup_envs.sh — build the ISOLATED environments for the Xenium cell-typing
# production run on HiPerGator. Mirrors the discipline of
# /blue/maigan/smith6jt/Xenium_Analysis/scripts/benchmark/setup_envs.sh:
# the production `xenium_analysis` conda env is NEVER modified, and R packages
# go into the project-local, already-writable library that
# scripts/benchmark/run_benchmark.sh:37 established as
# R_LIBS_USER=/blue/maigan/smith6jt/R_libs.
#
#   * Python  -> new conda env `xct_py` (anndata / scipy / sklearn / h5py /
#                pyarrow). It is deliberately small: this lane does not need
#                scanpy or squidpy, and pinning fewer things is what keeps it
#                from colliding with `xenium_analysis` or with `bench_py`.
#   * R       -> KODAMA at a PINNED COMMIT into R_LIBS_USER, plus its deps.
#
# WHY THE PIN IS NOT OPTIONAL HERE. The benchmark's setup_envs.sh:115 runs
#     remotes::install_github("tkcaccia/KODAMA", lib = lib, upgrade = "never")
# — unpinned, so it takes upstream HEAD on whatever day it runs. Every KODAMA
# timing this deployment is sized from, and the shard-equivalence proof the run
# depends on, were made against smith6jt-cop/KODAMA@d2dbe30ee6 (reports itself
# as 3.3, dated 2026-03-17). Upstream tkcaccia/KODAMA has moved since
# (8be47354, 2026-08-22). An unpinned install would put a DIFFERENT KODAMA
# behind the same version string, and `kodama_merge.R` re-implements ~20 lines
# of KODAMA.matrix — so a silent upgrade is exactly the failure that would not
# announce itself. The pin is recorded in kodama.lock and re-proved by
# check_shard_equivalence.R before every production submission.
#
# Idempotent. Usage:
#   bash setup_envs.sh              # both
#   bash setup_envs.sh --py-only
#   bash setup_envs.sh --r-only
# =============================================================================
set -uo pipefail

CONDA_BASE="${CONDA_BASE:-/blue/maigan/smith6jt/miniforge3}"
PY_ENV="${PY_ENV:-xct_py}"
export R_LIBS_USER="${R_LIBS_USER:-/blue/maigan/smith6jt/R_libs}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The one line that decides which KODAMA the whole run is.
KODAMA_REPO="${KODAMA_REPO:-smith6jt-cop/KODAMA}"
KODAMA_SHA="${KODAMA_SHA:-d2dbe30ee66509b82c084616d6961c5b292cd059}"

DO_PY=1; DO_R=1
case "${1:-}" in
  --py-only) DO_R=0 ;;
  --r-only)  DO_PY=0 ;;
  "") ;;
  *) echo "unknown arg: $1"; exit 2 ;;
esac

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# -----------------------------------------------------------------------------
# 1. Python env: xct_py
# -----------------------------------------------------------------------------
if [[ $DO_PY -eq 1 ]]; then
  echo "=============================================================="
  echo "[${PY_ENV}] Python env for prep / calling / verification"
  echo "=============================================================="
  if conda env list | grep -qE "^${PY_ENV}\s"; then
    echo "[${PY_ENV}] exists — skipping create"
  else
    mamba create -y -n "${PY_ENV}" -c conda-forge \
        python=3.12 numpy scipy pandas scikit-learn h5py pyarrow anndata \
        matplotlib pytest \
      || { echo "[${PY_ENV}] create FAILED"; exit 1; }
  fi
  conda activate "${PY_ENV}"
  echo "[${PY_ENV}] python: $(python --version 2>&1) @ $(which python)"
  python - <<'PY' || { echo "[xct_py] import check FAILED"; exit 1; }
import importlib
for m in ["numpy","scipy","pandas","sklearn","h5py","pyarrow","anndata"]:
    mod = importlib.import_module(m)
    print(f"  {m:10s} {getattr(mod,'__version__','?')}")
print("xct_py OK")
PY
  conda deactivate
fi

# -----------------------------------------------------------------------------
# 2. R: KODAMA at the pinned commit, into the project R library
# -----------------------------------------------------------------------------
if [[ $DO_R -eq 1 ]]; then
  echo "=============================================================="
  echo "[R] KODAMA ${KODAMA_REPO}@${KODAMA_SHA:0:10} -> ${R_LIBS_USER}"
  echo "=============================================================="
  conda activate xenium_analysis        # reuse its R only; install nothing into it
  mkdir -p "${R_LIBS_USER}"
  echo "[R] Rscript: $(which Rscript)"

  KODAMA_REPO="${KODAMA_REPO}" KODAMA_SHA="${KODAMA_SHA}" \
  Rscript --no-save - <<'RS' || { echo "[R] install FAILED"; exit 1; }
options(repos = c(CRAN = "https://cloud.r-project.org"))
lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
repo <- Sys.getenv("KODAMA_REPO"); sha <- Sys.getenv("KODAMA_SHA")
cat("R lib target:", lib, "\n"); cat(R.version.string, "\n")

have <- function(p) requireNamespace(p, quietly = TRUE)
# Rnanoflann is KODAMA's kNN backend and kodama_merge.R calls it directly;
# uwot/Rtsne back KODAMA.visualization; digest and jsonlite are ours.
need <- c("remotes", "Matrix", "Rcpp", "RcppArmadillo", "Rnanoflann",
          "uwot", "Rtsne", "irlba", "jsonlite", "digest")
todo <- need[!vapply(need, have, logical(1))]
if (length(todo)) {
  cat("installing CRAN:", paste(todo, collapse = ", "), "\n")
  install.packages(todo, lib = lib, Ncpus = max(1, parallel::detectCores() - 1))
} else cat("CRAN deps already present\n")

want_sha <- sha
cur <- if (have("KODAMA")) {
  d <- utils::packageDescription("KODAMA")
  if (!is.null(d$GithubSHA1)) d$GithubSHA1 else if (!is.null(d$RemoteSha)) d$RemoteSha else ""
} else ""
if (!identical(cur, want_sha)) {
  cat(sprintf("installing KODAMA %s@%s (had: %s)\n", repo, want_sha,
              if (nzchar(cur)) cur else "<none>"))
  remotes::install_github(paste0(repo, "@", want_sha), lib = lib, upgrade = "never")
} else cat("KODAMA already at the pinned commit\n")
RS

  echo "[R] verifying and writing the lock ..."
  Rscript --no-save - > "${HERE}/kodama.lock" <<'RS' || { echo "[R] verification FAILED"; exit 1; }
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(KODAMA); library(Rnanoflann); library(uwot); library(jsonlite)})
d <- utils::packageDescription("KODAMA")
sha <- if (!is.null(d$GithubSHA1)) d$GithubSHA1 else d$RemoteSha
want <- Sys.getenv("KODAMA_SHA", sha)
if (nzchar(want) && !identical(sha, want))
  stop(sprintf("installed KODAMA sha %s != pinned %s", sha, want))
jsonlite::write_json(list(
  kodama_repo = paste(d$RemoteUsername, d$RemoteRepo, sep = "/"),
  kodama_sha = sha, kodama_version = as.character(packageVersion("KODAMA")),
  r_version = R.version.string, r_libs_user = Sys.getenv("R_LIBS_USER"),
  Rnanoflann = as.character(packageVersion("Rnanoflann")),
  uwot = as.character(packageVersion("uwot")),
  host = Sys.info()[["nodename"]],
  written_at_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")),
  stdout(), auto_unbox = TRUE, pretty = TRUE)
RS
  cat "${HERE}/kodama.lock"
  conda deactivate
fi

echo "=============================================================="
echo "setup_envs.sh complete.  Next: bash ${HERE}/../hpg/submit_all.sh --dry-run"
echo "=============================================================="
