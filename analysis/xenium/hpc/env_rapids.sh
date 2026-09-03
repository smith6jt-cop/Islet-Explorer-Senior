# Source-able env block for the rapids/cudf-aware xenium_analysis environment.
# Required so that `import squidpy` / `import rapids_singlecell` / `import cudf`
# do not crash with the pynvjitlink/numba_cuda incompatibility documented in
# CLAUDE.md §7. Idempotent — safe to source repeatedly.
#
# Usage (interactive shell):
#   source /blue/maigan/smith6jt/Xenium_Analysis/scripts/env_rapids.sh
#
# Usage (driver scripts): source after `conda activate xenium_analysis` so that
# $CONDA_PREFIX resolves to the env's prefix.

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "env_rapids.sh: CONDA_PREFIX is unset — activate xenium_analysis first" >&2
    return 1 2>/dev/null || exit 1
fi

export CUDA_PATH="${CUDA_PATH:-/apps/compilers/cuda/12.8.1}"
export CUDA_HOME="${CUDA_HOME:-/apps/compilers/cuda/12.8.1}"

case ":${PATH}:" in
    *":${CUDA_PATH}/bin:"*) : ;;
    *) export PATH="${CUDA_PATH}/bin:${PATH}" ;;
esac

case ":${LD_LIBRARY_PATH:-}:" in
    *":${CUDA_PATH}/lib64:"*) : ;;
    *) export LD_LIBRARY_PATH="${CUDA_PATH}/lib64:${CUDA_HOME}/nvvm/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}" ;;
esac

_xenium_libcudart="${CUDA_HOME}/lib64/libcudart.so.12"
_xenium_libstdcxx="${CONDA_PREFIX}/lib/libstdc++.so.6"
case ":${LD_PRELOAD:-}:" in
    *":${_xenium_libcudart}:"*) : ;;
    *) export LD_PRELOAD="${_xenium_libstdcxx}:${_xenium_libcudart}${LD_PRELOAD:+:${LD_PRELOAD}}" ;;
esac
unset _xenium_libcudart _xenium_libstdcxx

export NUMBA_CUDA_ENABLE_PYNVJITLINK=1
