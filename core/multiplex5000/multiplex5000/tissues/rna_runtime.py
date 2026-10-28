"""Numerical runtime for reproducible pancreas and pLN Xenium clustering.

The seeded clustering methods change calls when BLAS, OpenMP or Numba thread
counts change. Load the libraries they use before setting threadpool limits;
limiting only the libraries loaded at Stage 06 entry misses later imports.
"""

from __future__ import annotations

from contextlib import contextmanager
from importlib.metadata import version

NUMERICAL_PACKAGES = ("numpy", "scipy", "scikit-learn", "scanpy", "anndata",
                      "igraph", "numba", "pynndescent", "threadpoolctl")


def versions() -> dict[str, str]:
    """Record the installed numerical libraries that can change a partition."""
    return {name: version(name) for name in NUMERICAL_PACKAGES}


@contextmanager
def fixed_threads(count: int):
    """Use the project's thread count for every numerical pool in an RNA worker."""
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("RNA compute threads must be a positive integer")

    import igraph  # noqa: F401
    import numba
    import numpy  # noqa: F401
    import pynndescent  # noqa: F401
    import scanpy  # noqa: F401
    import scipy  # noqa: F401
    import sklearn.decomposition  # noqa: F401
    import sklearn.neighbors  # noqa: F401
    from threadpoolctl import threadpool_info, threadpool_limits

    previous = numba.get_num_threads()
    try:
        numba.set_num_threads(count)
    except ValueError as exc:
        raise ValueError(f"RNA compute threads {count} exceed this Numba runtime's maximum") from exc
    try:
        with threadpool_limits(limits=count):
            pools = [{"api": item["internal_api"], "version": item.get("version"),
                      "threads": item["num_threads"]} for item in threadpool_info()]
            if numba.get_num_threads() != count or any(
                    item["threads"] != count for item in pools):
                raise RuntimeError("RNA numerical libraries did not adopt the project thread count")
            yield {"numba_threads": count, "pools": pools}
    finally:
        numba.set_num_threads(previous)
