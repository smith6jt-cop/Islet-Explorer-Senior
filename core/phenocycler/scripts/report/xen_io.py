"""Xenium inputs for the QC figures: the bundle (every cell) and the processed base (retained cells).

The bundle is the denominator.  ``cells.parquet`` lists every segmented cell with its per-type
count columns; ``cell_feature_matrix.h5`` holds the six feature types the panel decodes.  The base
(``/data/xenium/processed/<serial>/panc/<serial>_phenotyped.h5ad``) holds the cells that survived
``build_phenotyped_base.qc_mask`` and records, in ``uns['prep']['qc']``, the thresholds and the
per-mask drop counts.  Those counts are recomputed here by calling the SAME ``qc_mask`` on the
same gene matrix and compared with the recorded ones, so X01/X02 draw a reproduction rather than
a transcript.

The base is read with ``h5py`` rather than ``anndata`` on purpose: the figures need ``obs`` and
the per-cell non-zero count of ``layers['counts']`` (genes per cell), which is ``diff(indptr)`` --
no reason to materialise 641k x 5101 floats twice.
"""
from __future__ import annotations

import pathlib
import sys
from dataclasses import dataclass, field

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402
from build_phenotyped_base import qc_mask, MIN_GENES, MIN_COUNTS, DENSITY_QUANTILE  # noqa: E402  (scripts/xenium_celltype)

GENE = "Gene Expression"


def _decode(v):
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray) and v.dtype.kind in "SO":
        return [x.decode() if isinstance(x, bytes) else x for x in v.tolist()]
    if isinstance(v, np.generic):
        return v.item()
    return v


def read_prep(base: pathlib.Path) -> dict:
    """``uns['prep']`` as plain Python."""
    def walk(g):
        out = {}
        for k, v in g.items():
            out[k] = walk(v) if isinstance(v, h5py.Group) else _decode(v[()])
        return out
    with h5py.File(base, "r") as f:
        return walk(f["uns/prep"])


def read_base_obs(base: pathlib.Path) -> pd.DataFrame:
    """cell_id-indexed obs with ``n_genes_base`` = non-zeros per row of ``layers['counts']``."""
    with h5py.File(base, "r") as f:
        idx = _decode(f["obs/_index"][:])
        obs = pd.DataFrame(index=pd.Index(idx, name="cell_id", dtype=object))
        for c in ("cell_area", "nucleus_area", "transcript_counts"):
            obs[c] = f["obs"][c][:]
        ip = f["layers/counts/indptr"][:]
        obs["n_genes_base"] = np.diff(ip).astype(int)
    return obs


def read_matrix(bundle: pathlib.Path) -> tuple:
    """(barcodes, feature_type[], names[], X csr cells x features) -- every feature type, every cell."""
    with h5py.File(bundle / "cell_feature_matrix.h5", "r") as f:
        g = f["matrix"]
        X = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:])[::-1])
        ftype = np.array(_decode(g["features/feature_type"][:]))
        names = np.array(_decode(g["features/name"][:]))
        bar = np.array(_decode(g["barcodes"][:]))
    return bar, ftype, names, X


def assert_feature_types(ftype: np.ndarray) -> dict:
    got = pd.Series(ftype).value_counts().to_dict()
    if got != pins.XENIUM_FEATURE_TYPES:
        raise RuntimeError(f"cell_feature_matrix.h5 feature types {got} != pinned {pins.XENIUM_FEATURE_TYPES}")
    return got


@dataclass
class XenTables:
    donor: str
    serial: str
    status: str
    bundle: pathlib.Path
    base: pathlib.Path
    cells: pd.DataFrame          # every bundle cell (matrix barcode order): counts, areas, n_genes, n_counts, masks
    prep: dict                   # uns['prep'] of the base
    qc_recomputed: dict          # qc_mask() on the bundle gene matrix
    feature_types: dict          # type -> n features
    per_type_total: dict         # type -> total decoded counts over all cells
    per_type_mean_per_cell: dict
    column_check: dict           # cells.parquet count columns vs matrix per-type totals
    metrics_summary: pd.DataFrame | None
    n_base: int
    inputs: list = field(default_factory=list)

    @property
    def n(self) -> int:
        return int(len(self.cells))


def load_tables(donor: str) -> XenTables:
    bundle = pins.xenium_bundle(donor)
    base = pins.xenium_base(donor)
    serial = pins.XENIUM_SERIAL[str(donor)]
    prep = read_prep(base)
    if str(prep.get("serial")) != serial or str(prep.get("donor_id")) != str(donor):
        raise RuntimeError(f"base {base} is stamped serial={prep.get('serial')} donor={prep.get('donor_id')}, "
                           f"expected {serial}/{donor}")
    bar, ftype, names, X = read_matrix(bundle)
    feature_types = assert_feature_types(ftype)
    cells = pd.read_parquet(bundle / "cells.parquet")
    cells["cell_id"] = cells["cell_id"].astype(str)
    pos = pd.Series(np.arange(len(cells)), index=cells["cell_id"].to_numpy())
    hit = pos.reindex(bar)
    if hit.isna().any() or len(cells) != len(bar):
        raise RuntimeError(f"donor {donor}: cells.parquet and cell_feature_matrix.h5 do not list the same cells")
    cells = cells.iloc[hit.to_numpy().astype(int)].reset_index(drop=True)      # matrix barcode order

    gene = ftype == GENE
    Xg = X[:, gene].tocsr()
    n_counts = np.asarray(Xg.sum(axis=1)).ravel().astype(int)
    n_genes = np.asarray((Xg > 0).sum(axis=1)).ravel().astype(int)
    cells["n_counts_matrix"] = n_counts
    cells["n_genes"] = n_genes
    keep, qc = qc_mask(Xg, cells["cell_area"].to_numpy(float))
    density = np.divide(n_counts, cells["cell_area"].to_numpy(float), out=np.zeros(len(cells)),
                        where=cells["cell_area"].to_numpy(float) > 0)
    cells["density"] = density
    cells["fail_min_genes"] = n_genes < MIN_GENES
    cells["fail_min_counts"] = n_counts < MIN_COUNTS
    cells["fail_density"] = density > qc["density_cap_counts_per_um2"]
    cells["retained_recomputed"] = keep
    obs = read_base_obs(base)
    cells["retained_base"] = cells["cell_id"].isin(obs.index).to_numpy()
    cells["n_genes_base"] = obs["n_genes_base"].reindex(cells["cell_id"]).to_numpy()
    qc["n_retained_recomputed"] = int(keep.sum())
    qc["n_retained_base"] = int(len(obs))
    qc["retained_sets_identical"] = bool((cells["retained_recomputed"] == cells["retained_base"]).all())
    qc["matches_uns_prep"] = all(int(qc[k]) == int(prep["qc"][k]) for k in
                                 ("dropped_min_genes", "dropped_min_counts", "dropped_density")) \
        and abs(float(qc["density_cap_counts_per_um2"]) - float(prep["qc"]["density_cap_counts_per_um2"])) < 1e-9
    qc["n_genes_base_equals_matrix"] = bool(
        np.array_equal(cells.loc[cells["retained_base"], "n_genes_base"].to_numpy().astype(int),
                       cells.loc[cells["retained_base"], "n_genes"].to_numpy().astype(int)))
    qc["transcript_counts_equals_gene_sum"] = bool(np.array_equal(cells["transcript_counts"].to_numpy(int), n_counts))

    totals = np.bincount(X.indices, weights=X.data, minlength=X.shape[1])
    per_type_total = {t: float(totals[ftype == t].sum()) for t in pins.XENIUM_FEATURE_TYPES}
    per_type_mean = {t: v / len(cells) for t, v in per_type_total.items()}
    column_check = {}
    for col, t in pins.XENIUM_COUNT_COLUMNS.items():
        if col in cells.columns:
            column_check[t] = {"column": col, "column_total": float(cells[col].sum()),
                               "matrix_total": per_type_total[t],
                               "equal": bool(abs(float(cells[col].sum()) - per_type_total[t]) < 0.5)}
    ms_path = pins.metrics_summary_path(donor)
    ms = pd.read_csv(ms_path) if ms_path else None
    inputs = [bundle / "cells.parquet", bundle / "cell_feature_matrix.h5", base, pins.XENIUM_PATHS_CSV] \
        + ([ms_path] if ms_path else [])
    return XenTables(donor=str(donor), serial=serial, status=pins.donor_status(donor), bundle=bundle, base=base,
                     cells=cells, prep=prep, qc_recomputed=qc, feature_types=feature_types,
                     per_type_total=per_type_total, per_type_mean_per_cell=per_type_mean, column_check=column_check,
                     metrics_summary=ms, n_base=int(len(obs)), inputs=inputs)


QC_CONSTANTS = {"min_genes": MIN_GENES, "min_counts": MIN_COUNTS, "density_quantile": DENSITY_QUANTILE}
