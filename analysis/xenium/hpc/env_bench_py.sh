#!/usr/bin/env bash
# Source this to enter the bench_py env with the libstdc++ preload that BANKSY
# and UTAG need (their matplotlib needs a newer CXXABI than the system
# /lib64/libstdc++.so.6 — same class of fix as scripts/env_rapids.sh).
#   source scripts/benchmark/env_bench_py.sh
source "${CONDA_BASE:-/blue/maigan/smith6jt/miniforge3}/etc/profile.d/conda.sh"
conda activate bench_py
export LD_PRELOAD="${CONDA_PREFIX}/lib/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}"
