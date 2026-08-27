#!/usr/bin/env bash
# Source this, never execute it. One place that says where everything is on HiPerGator, so
# no SLURM script carries its own copy of a path. Mirrors the layout the sibling benchmark
# established: work under /blue/maigan/smith6jt, conda at .../miniforge3, R packages in
# .../R_libs (see Xenium_Analysis/scripts/benchmark/run_benchmark.sh:30,37).
#
#   XCT_ROOT   work area (code + runs + logs); /blue is the fast scratch, /orange is data
#   XCT_RUN    the run directory; every stage reads work_manifest.json from here
#   XENIUM_SITE  selects the bundle/h5ad layout in sites.json
export XCT_ROOT="${XCT_ROOT:-/blue/maigan/smith6jt/Xenium_Celltype}"
export XCT_CODE="${XCT_CODE:-${XCT_ROOT}/code}"
export XCT_RUNS="${XCT_RUNS:-${XCT_ROOT}/runs}"
export XCT_RUN="${XCT_RUN:-${XCT_RUNS}/${XCT_RUN_ID:-current}}"
export XCT_HPG="${XCT_CODE}/core/phenocycler/scripts/xenium_celltype/hpg"
export XCT_LANE="${XCT_CODE}/core/phenocycler/scripts/xenium_celltype"

export CONDA_BASE="${CONDA_BASE:-/blue/maigan/smith6jt/miniforge3}"
export R_LIBS_USER="${R_LIBS_USER:-/blue/maigan/smith6jt/R_libs}"
export XENIUM_SITE="${XENIUM_SITE:-hipergator}"
export XENIUM_SITES="${XENIUM_SITES:-${XCT_HPG}/sites.json}"
export PYTHONPATH="${XCT_CODE}/core/phenocycler:${PYTHONPATH:-}"

# BLAS threads are set explicitly because R's kmeans and Rnanoflann already parallelise; an
# OpenBLAS that also grabs every core turns a 2-cpu task into a 128-thread thrash and SLURM
# charges for it either way.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

xct_py()  { source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate xct_py; }
xct_R()   { source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate xenium_analysis; }
# Manifest reads need a python; every SLURM script therefore calls `xct_py` BEFORE it calls
# xct_task or any inline python, and only then switches to `xct_R`. conda activate replaces
# rather than stacks, so the order is: xct_py (read the manifest) -> xct_R (run KODAMA).
xct_task() {  # $1 = manifest key ("shard_tasks"), $2 = index -> one JSON object on stdout
  python - "$1" "$2" <<'PYX'
import json, os, sys
m = json.load(open(os.path.join(os.environ["XCT_RUN"], "work_manifest.json")))
print(json.dumps(m[sys.argv[1]][int(sys.argv[2])]))
PYX
}
