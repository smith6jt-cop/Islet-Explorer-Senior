"""Loading Xenium expression + spatial for one donor, from either of the two sources.

There are two, and they are NOT interchangeable by donor:

* ``h5ad_out/<donor>_PancXenium.h5ad`` — 16 of the 20 run-cohort donors. Raw counts and
  donor metadata; **no** cell types, no ``obsm['spatial']``, no layers.
* the bundle's ``cell_feature_matrix.h5`` — the remaining four (6450, 6505, 6534, 6623,
  i.e. the whole of Xenium runs MM046 and MM047), which have no h5ad at all.

Dropping the bundle path would silently delete two of the eleven Xenium runs, and those two
carry two of the cohort's three Aab+ donors. Cross-validated on donor 6539, which has both:
identical gene order, scores agreeing to 6e-7, 100% argmax agreement over 20,000 cells. So
the choice of source is a fact about which files exist, not a methodological fork.

Spatial always comes from the bundle's ``cells.parquet`` (``x_centroid``/``y_centroid``, in
microns) — neither expression source carries coordinates. The h5ad is a strict subset of the
bundle's cells (317,610 of 319,005 for 6539; the difference is the Xenium QC drop), so the
join is an inner join on ``cell_id`` and the dropped count is reported, never assumed.
"""

from __future__ import annotations

import glob
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

XENIUM_ROOT = Path("/data/xenium")
H5AD_DIR = XENIUM_ROOT / "h5ad_out"

#: Feature rows in cell_feature_matrix.h5 that are real genes. The file also carries
#: negative-control probes/codewords and deprecated codewords — 5,101 genes out of 9,638
#: or 9,644 features depending on the run's control-set version.
GENE_FEATURE = "Gene Expression"


class XeniumIOError(RuntimeError):
    """A donor's expression or spatial data cannot be loaded — never silently skipped."""


@dataclass(frozen=True)
class DonorData:
    donor_id: str
    cell_id: np.ndarray
    genes: list[str]
    counts: sp.csr_matrix     # cells x genes, RAW
    xy: np.ndarray            # cells x 2, microns
    source: str               # "h5ad" | "bundle"
    n_dropped_no_spatial: int


def _bundle_dir(donor_id: str, xenium_paths_csv: Path, roi: str = "panc") -> Path:
    xen = pd.read_csv(xenium_paths_csv, dtype=str)
    hit = xen[(xen.donor_id == donor_id) & (xen.roi == roi)]
    hit = hit[hit.source_path.notna() & (hit.source_path != "NA")]
    if hit.empty:
        raise XeniumIOError(f"donor {donor_id} roi={roi} has no bundle path recorded")
    return XENIUM_ROOT / str(hit.iloc[0].source_path).rstrip("/")


def _read_bundle_matrix(bundle: Path):
    import h5py
    p = bundle / "cell_feature_matrix.h5"
    if not p.exists():
        raise XeniumIOError(f"{p} is absent; transfer the bundle's cell_feature_matrix.h5")
    with h5py.File(p, "r") as f:
        g = f["matrix"]
        # stored CSC over (features, cells) -> the same buffers are CSR over (cells, features)
        X = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]),
                          shape=tuple(g["shape"][:])[::-1])
        name = np.array([x.decode() for x in g["features"]["name"][:]])
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
        bar = np.array([x.decode() for x in g["barcodes"][:]])
    keep = ftype == GENE_FEATURE
    if not keep.any():
        raise XeniumIOError(f"{p} carries no '{GENE_FEATURE}' features")
    return bar, list(name[keep]), X[:, keep].tocsr()


def _h5ad_path(donor_id: str, roi: str = "panc") -> Path | None:
    tag = "Panc" if roi == "panc" else "pLN"
    for p in sorted(H5AD_DIR.glob(f"*_{tag}Xenium.h5ad")):
        stem = p.name.split("_")[0]
        if stem.lstrip("HDL").rstrip("abc") == donor_id:
            return p
    return None


def load_donor(donor_id: str, xenium_paths_csv: Path, *, roi: str = "panc",
               prefer: str = "h5ad") -> DonorData:
    """Counts + microns for one donor, from whichever source exists.

    ``prefer`` picks the source when both exist. It changes which cells are present (the
    h5ad is the upstream QC subset) but not the values — see the module docstring.
    """
    bundle = _bundle_dir(donor_id, xenium_paths_csv, roi)
    cells_p = bundle / "cells.parquet"
    if not cells_p.exists():
        raise XeniumIOError(f"{cells_p} is absent; spatial coordinates live only there")
    cells = pd.read_parquet(cells_p, columns=["cell_id", "x_centroid", "y_centroid"])
    cells["cell_id"] = cells.cell_id.astype(str)

    h5 = _h5ad_path(donor_id, roi) if prefer == "h5ad" else None
    if h5 is not None:
        import anndata as ad
        a = ad.read_h5ad(h5)
        ids = np.array(list(map(str, a.obs_names)))
        genes = list(map(str, a.var_names))
        counts = sp.csr_matrix(a.X)
        source = "h5ad"
    else:
        ids, genes, counts = _read_bundle_matrix(bundle)
        source = "bundle"

    pos = pd.Series(np.arange(len(cells)), index=cells.cell_id.values)
    hit = pos.reindex(ids)
    ok = hit.notna().to_numpy()
    if not ok.any():
        raise XeniumIOError(f"donor {donor_id}: no cell_id from {source} matched cells.parquet")
    rows = hit[ok].to_numpy().astype(int)
    xy = cells.loc[:, ["x_centroid", "y_centroid"]].to_numpy()[rows]

    return DonorData(donor_id=donor_id, cell_id=ids[ok], genes=genes,
                     counts=counts[ok].tocsr(), xy=xy, source=source,
                     n_dropped_no_spatial=int((~ok).sum()))


def lognorm(counts: sp.csr_matrix, *, target: float | None = None) -> sp.csr_matrix:
    """Per-cell depth normalisation to `target` (default: the median depth) then log1p.

    For an argmax across lineages within a cell this is nearly a no-op — measured on 6539,
    raw and lognorm agree to within 0.2 percentage points on every lineage share, because a
    per-cell scalar cancels in a within-cell comparison. It matters here because the
    embedding compares cells to each other, where depth is a genuine confound.
    """
    depth = np.asarray(counts.sum(axis=1)).ravel()
    if target is None:
        target = float(np.median(depth[depth > 0])) if (depth > 0).any() else 1.0
    scale = np.divide(target, depth, out=np.zeros_like(depth, dtype=float), where=depth > 0)
    return sp.csr_matrix(counts.multiply(scale[:, None]).log1p())
