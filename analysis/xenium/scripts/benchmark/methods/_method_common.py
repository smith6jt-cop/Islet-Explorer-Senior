"""Shared helpers for the Python method runners (BANKSY, UTAG, baseline).

Importable in both the `bench_py` env (BANKSY/UTAG) and `xenium_analysis`
(baseline) — only imports scanpy/sklearn lazily inside functions so importing
this module is cheap.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_BENCH = _HERE.parent.parent            # scripts/benchmark
for p in (str(_BENCH.parent.parent), str(_BENCH)):   # repo root + scripts/benchmark
    if p not in sys.path:
        sys.path.insert(0, p)

import bench_common as bc  # noqa: E402  (resolved via scripts/benchmark on path)


def load_input(sample: str, variant: str):
    import scanpy as sc
    p = bc.input_paths(sample, variant)["h5ad"]
    if not p.exists():
        raise FileNotFoundError(f"missing input — run prep first: {p}")
    return sc.read_h5ad(p)


def cell_ids(adata) -> np.ndarray:
    if "cell_id" in adata.obs.columns:
        return adata.obs["cell_id"].astype(str).to_numpy()
    return adata.obs_names.astype(str).to_numpy()


def write_labels(sample: str, variant: str, method: str,
                 ids: np.ndarray, labels: np.ndarray) -> Path:
    bc.results_dir(sample).mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({"cell_id": np.asarray(ids).astype(str),
                       "domain": np.asarray(labels).astype(str)})
    out = bc.labels_path(sample, variant, method)
    df.to_csv(out, index=False)
    print(f"  [{method}] wrote {df['domain'].nunique()} domains x {len(df):,} cells -> {out.name}")
    return out


def pick_closest_label_set(label_sets: dict, target: int):
    """Given {key: 1d-label-array}, return (key, labels, n) whose unique-count is
    nearest `target` (ties -> more domains, then smaller key)."""
    best = None
    for key, lab in label_sets.items():
        n = int(pd.unique(np.asarray(lab)).size)
        score = (abs(n - target), -n, str(key))
        if best is None or score < best[0]:
            best = (score, key, np.asarray(lab), n)
    return best[1], best[2], best[3]


def calibrate_leiden(adata, use_rep: str, target: int, seed: int = 0,
                     neighbors_k: int = 15, fixed_res: float | None = None,
                     res_grid=(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.3, 1.6, 2.0)):
    """Neighbours on `use_rep`, then leiden. If fixed_res: one run. Else sweep
    res_grid and return the partition whose cluster count is closest to target.
    Returns (labels, resolution, n_clusters)."""
    import scanpy as sc
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=neighbors_k, random_state=seed)
    if fixed_res is not None:
        sc.tl.leiden(adata, resolution=fixed_res, random_state=seed,
                     flavor="igraph", n_iterations=2, directed=False)
        lab = adata.obs["leiden"].to_numpy()
        return lab, fixed_res, int(pd.unique(lab).size)
    sets = {}
    for r in res_grid:
        sc.tl.leiden(adata, resolution=r, random_state=seed, flavor="igraph",
                     n_iterations=2, directed=False, key_added="_leiden_cal")
        sets[r] = adata.obs["_leiden_cal"].to_numpy().copy()
    res, lab, n = pick_closest_label_set(sets, target)
    return lab, res, n
