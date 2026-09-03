"""The KINTSUGI clustering, imported from a pinned checkout and run headless.

KINTSUGI's analysis end is ``notebooks/4_Segmentation_Analysis.ipynb`` plus the module it was
extracted into, ``notebooks/Kanalysis.py``.  Its clustering is a **self-organising map**:
``Kanalysis.visualize_som`` scales the marker block with ``StandardScaler`` and calls
``pyFlowSOM.som`` / ``map_data_to_nodes``.  The notebook then takes the SOM label into scanpy,
runs Leiden beside it, crosses the two into ``combined_clusters``, and collapses that back
down by correlation-distance Ward linkage.

Three things about that pipeline are recorded here rather than silently reproduced or
silently dropped.

**1. The progressive loop is a figure, not a method.**  ``visualize_som`` loops
``for rlen_iter in range(1, num_passes + 1)`` and calls ``pyFlowSOM.som(..., rlen_iter, ...)``
from scratch each time, rendering a frame per iteration into a ``stackview`` stack.  Only the
last iteration is the result, and it is exactly ``som(..., rlen=num_passes, ...)``.  Running
the loop headless would cost ``num_passes`` full SOM fits to reach an identical answer, so
:func:`som_cluster` calls it once.  :func:`som_cluster_progressive` reproduces the loop for
anyone who wants the frames.

**2. pyFlowSOM is not a declared KINTSUGI dependency.**  It is a lazy import inside
``Kanalysis`` and appears in no extra, in ``requirements.txt``, in ``constraints.txt`` or in
any ``envs/*.yml`` at the pinned commit.  So ``kintsugi install analysis`` does **not** supply
it and it has to be installed explicitly.  This module says so in the error rather than
letting an ImportError surface from three frames down.

**3. The notebook fits PCA before it normalises, and this lane does not.**  In
``4_Segmentation_Analysis.ipynb`` cell 97 runs ``sc.pp.pca`` and cells 99 runs
``sc.pp.normalize_total`` / ``log1p`` / ``scale`` *after* it.  ``sc.pp.neighbors`` at cell 101
then reuses the ``X_pca`` computed from unnormalised values, so the neighbour graph, PAGA,
UMAP and Leiden all descend from the pre-normalisation embedding.  That is an ordering
mistake, not a design choice -- nothing downstream reads the normalised PCA.  Reproducing it
would put a known defect underneath every number this lane produces, so :func:`embed` runs the
documented order (normalise -> log1p -> scale -> PCA -> neighbours) and
``--notebook-pca-order`` restores the notebook's for comparison.  The deviation is declared
here, in the README, and in the output manifest.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
import panel                                        # noqa: E402

#: The KINTSUGI checkout this lane was built against -- a gitignored clone under /external/,
#: following the convention ``uniform_adapter.py`` already uses.  KINTSUGI is pip-installable,
#: but the analysis end lives in ``notebooks/`` and is not part of the installed package, so
#: pinning the distribution version would not pin the code this lane actually runs.
KINTSUGI_COMMIT = "716d1ebd87db916ab04910c6374e04d2d9ea508e"
KINTSUGI_VERSION = "1.2.3"
KINTSUGI_PATH = panel.REPO / "external" / "KINTSUGI"

#: ``visualize_som``'s defaults are ``num_nodes=10, num_passes=5, lr_start=0.05, lr_end=0.01``;
#: the values below are what notebook 4 actually passes at every call site (a 5x5 grid).  They
#: are the starting point for the pre-registered sweep, not the answer: 25 nodes was chosen
#: for a 26-marker lymph node, and this panel carries 58 markers against 37 candidate
#: phenotypes.  The sweep is reported in full; the grid is never chosen to make a cell-type
#: frequency come out.
NOTEBOOK_SOM = {
    "num_nodes": 25,          # xdim = ydim = 5
    "num_passes": 10,
    "lr_start": 0.06,
    "lr_end": 0.001,
    "seed": 627,
}

#: Notebook 4 cell 106 / 113: Leiden beside the SOM, crossed with it, then collapsed by
#: correlation-distance Ward linkage to a fixed count.  ``n_final_clusters = 23`` is hardcoded
#: there.  This lane keeps the *mechanism* and refuses the constant -- see ``meta_cluster``.
NOTEBOOK_LEIDEN_RESOLUTION = 0.6
NOTEBOOK_N_FINAL_CLUSTERS = 23


def _assert_pinned_checkout() -> pathlib.Path:
    if not KINTSUGI_PATH.exists():
        raise RuntimeError(
            f"KINTSUGI checkout missing at {KINTSUGI_PATH}.\n"
            f"  git clone https://github.com/smith6jt-cop/KINTSUGI.git external/KINTSUGI\n"
            f"  git -C external/KINTSUGI checkout {KINTSUGI_COMMIT}")
    head = subprocess.run(["git", "-C", str(KINTSUGI_PATH), "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    if head != KINTSUGI_COMMIT:
        raise RuntimeError(
            f"KINTSUGI checkout is at {head}, this lane is pinned to {KINTSUGI_COMMIT}. "
            "A different commit can change the clustering without anything in the output "
            "being able to show it.")
    return KINTSUGI_PATH


def import_kanalysis():
    """Import KINTSUGI's own ``Kanalysis`` from the pinned checkout."""
    path = _assert_pinned_checkout()
    sys.path.insert(0, str(path / "notebooks"))
    import Kanalysis                                # noqa: F401
    return Kanalysis


def import_pyflowsom():
    """Import pyFlowSOM, with the reason it is not already installed."""
    _assert_pinned_checkout()
    try:
        import pyFlowSOM
    except ImportError as exc:                      # pragma: no cover - environment issue
        raise RuntimeError(
            "pyFlowSOM is not installed. KINTSUGI imports it lazily inside "
            "Kanalysis.visualize_som and declares it in no extra, in requirements.txt, in "
            f"constraints.txt or in any envs/*.yml at {KINTSUGI_COMMIT[:7]}, so "
            "`kintsugi install analysis` does not supply it.\n"
            "  pip install -c external/KINTSUGI/constraints.txt pyFlowSOM"
        ) from exc
    return pyFlowSOM


def som_cluster(values: np.ndarray, channels: list[str], *,
                num_nodes: int = NOTEBOOK_SOM["num_nodes"],
                num_passes: int = NOTEBOOK_SOM["num_passes"],
                lr_start: float = NOTEBOOK_SOM["lr_start"],
                lr_end: float = NOTEBOOK_SOM["lr_end"],
                seed: int = NOTEBOOK_SOM["seed"]) -> tuple[np.ndarray, pd.DataFrame]:
    """``Kanalysis.visualize_som``'s statistical core, without the figure.

    Identical to the last frame of the notebook's progressive loop.  ``values`` is
    ``(n_cells, n_channels)`` on whatever scale the caller established -- for this lane, the
    UniFORM-aligned scale -- and the returned per-cluster profile is the mean of **those**
    values, not of the z-scaled ones.  That matters: the scaling exists so Euclidean distance
    is not dominated by bright markers, and it would make the profile unreadable against a
    gate stated on the aligned scale.  ``visualize_som`` does the same thing, by building
    ``train_data`` from the raw columns and scaling only ``train_array``.

    Returns ``(cluster_per_cell, profile)`` where ``profile`` is indexed ``Cluster_1..N``.
    """
    pyFlowSOM = import_pyflowsom()
    from sklearn.preprocessing import StandardScaler

    values = np.asarray(values, dtype=np.float64)
    if values.shape[1] != len(channels):
        raise ValueError(f"values has {values.shape[1]} columns, channels has {len(channels)}")
    if not np.isfinite(values).all():
        raise ValueError(
            "non-finite value in the SOM feature block. A missing measurement is "
            "unavailable, not zero -- drop the cell or the marker, never impute.")

    grid_size = int(np.sqrt(num_nodes))
    if grid_size * grid_size != num_nodes:
        raise ValueError(
            f"num_nodes={num_nodes} is not square; visualize_som takes grid_size = "
            "int(sqrt(num_nodes)) and would silently fit a smaller map")

    train_array = StandardScaler().fit_transform(values)
    som = pyFlowSOM.som(train_array, grid_size, grid_size, num_passes,
                        alpha_range=(lr_start, lr_end), seed=seed)
    clusters, _dist = pyFlowSOM.map_data_to_nodes(som, train_array)

    # visualize_som relabels to consecutive integers so empty nodes do not leave gaps.
    unique = np.unique(clusters)
    remap = {old: new for new, old in enumerate(unique)}
    clusters = np.array([remap[c] for c in clusters])

    frame = pd.DataFrame(values, columns=channels)
    frame["cluster"] = clusters
    profile = frame.groupby("cluster").mean()
    profile.index = [f"Cluster_{i + 1}" for i in range(len(profile))]
    return clusters, profile


def meta_cluster(profile: pd.DataFrame, n_final: int) -> dict[str, int]:
    """Notebook cell 113: correlation distance, Ward linkage, ``fcluster`` at ``maxclust``.

    The mechanism is KINTSUGI's and is kept.  The notebook's ``n_final_clusters = 23`` is a
    hardcoded constant and is **not** kept: a cluster count chosen to land near an expected
    number of cell types is the frequency reasoning ``CLAUDE.md``'s first non-negotiable
    forbids.  ``n_final`` is supplied by the caller from the pre-registered stability sweep,
    and the whole sweep is reported.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import pdist

    if n_final < 2 or n_final > len(profile):
        raise ValueError(f"n_final={n_final} outside [2, {len(profile)}]")
    linkage_matrix = linkage(pdist(profile, metric="correlation"), method="ward")
    assignments = fcluster(linkage_matrix, n_final, criterion="maxclust")
    return dict(zip(profile.index, (int(a) for a in assignments)))


def provenance() -> dict:
    """What produced these clusters, for the output manifest."""
    import importlib.metadata as md
    try:
        pyflowsom_version = md.version("pyFlowSOM")
    except Exception:                                # pragma: no cover
        pyflowsom_version = "not-installed"
    return {
        "kintsugi_commit": KINTSUGI_COMMIT,
        "kintsugi_version": KINTSUGI_VERSION,
        "pyflowsom_version": pyflowsom_version,
        "som_parameters": dict(NOTEBOOK_SOM),
        "deviation_from_notebook": [
            "PCA is fitted after normalise/log1p/scale, not before as in notebook cell 97; "
            "the notebook's neighbour graph descends from an unnormalised embedding.",
            "n_final_clusters comes from the stability sweep, not the notebook's hardcoded 23.",
        ],
    }
