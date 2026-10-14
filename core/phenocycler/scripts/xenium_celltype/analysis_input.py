"""Explicit, nucleus-aware section inputs shared by workstation and HPG consumers."""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import xenio


def section_entry(manifest, donor, roi, serial=None):
    rows = [d for d in manifest["donors"] if str(d["donor_id"]) == str(donor)
            and d.get("roi", "panc") == roi
            and (serial is None or str(d.get("serial", "")) == str(serial))]
    if len(rows) != 1:
        raise ValueError(f"expected one manifest section for {donor}/{roi}/{serial}; found {len(rows)}")
    return rows[0]


def section_dir(root, donor, roi, serial=None):
    dest = Path(root) / f"donor_id={donor}" / f"roi={roi}"
    return dest / f"serial={serial}" if serial else dest


def load_base(path, *, donor, roi, serial=None, detection_support=None):
    """Load precisely the common QC population; optional support never changes the base file.

    A support-masked matrix is a protein-guided comparison, not independent RNA evidence.
    """
    import anndata as ad

    path = Path(path).resolve()
    a = ad.read_h5ad(path)
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError("QC base cell and gene IDs must be unique")
    prep = a.uns.get("prep", {})
    if not prep.get("qc", {}).get("nucleus_required", False):
        raise ValueError("base lacks nucleus-aware QC provenance; rebuild it")
    for column, expected in (("donor_id", str(donor)), ("roi", roi), ("sample", serial)):
        if column not in a.obs or (expected is not None and
                                  set(a.obs[column].astype(str)) != {str(expected)}):
            raise ValueError(f"base {column} does not match requested section {expected}")
    for column in ("nucleus_area", "cell_area", "gdna_lambda"):
        if column not in a.obs or not np.isfinite(a.obs[column]).all():
            raise ValueError(f"base lacks finite {column}")
    if (a.obs.nucleus_area <= 0).any() or (a.obs.cell_area <= 0).any():
        raise ValueError("base includes cells without positive nucleus/cell area")
    if "counts" not in a.layers or "spatial" not in a.obsm:
        raise ValueError("base requires raw counts layer and spatial coordinates")
    counts = sp.csr_matrix(a.layers["counts"])
    xy = np.asarray(a.obsm["spatial"], dtype=float)
    if not np.isfinite(xy).all() or not np.isfinite(counts.data).all() or (counts.data < 0).any():
        raise ValueError("base contains invalid coordinates or counts")
    method = "rna_only"
    if detection_support is not None:
        support = Path(detection_support)
        ids = pd.read_csv(support / "cells.csv", dtype=str)["cell_id"].to_numpy()
        genes = pd.read_csv(support / "genes.csv", dtype=str)["gene"].to_numpy()
        if not np.array_equal(ids, a.obs_names.astype(str)) or not np.array_equal(genes, a.var_names.astype(str)):
            raise ValueError("detection support axes do not exactly match the QC base")
        mask = sp.load_npz(support / "support.npz").tocsr()
        if mask.shape != counts.shape or mask.dtype != np.dtype(bool):
            raise ValueError("detection support must be a boolean matrix matching the QC base")
        counts = counts.multiply(mask).tocsr()
        method = "protein_guided"
    resolved_serial = str(prep.get("serial", serial or ""))
    dd = xenio.DonorData(str(donor), a.obs_names.astype(str).to_numpy(),
                         list(a.var_names.astype(str)), counts, xy, "qc_base", 0,
                         serial=resolved_serial, bundle=Path(prep["bundle"]) if prep.get("bundle") else None)
    return dd, {"base": str(path), "method": method,
                "detection_support": str(Path(detection_support).resolve()) if detection_support else "",
                "serial": resolved_serial, "roi": roi, "donor_id": str(donor)}
