"""Read one accepted Xenium region without donor-level source substitution.

The matrix and cell geometry are separate instrument files. Their native cell IDs
are only unique within a region; downstream tables must use ``cell_uid``. This
module reads raw gene counts and records nucleus eligibility, but makes no cell
type call and applies no expression cutoff.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pyarrow.parquet as pq
from scipy import sparse

from .xenium_sources import XeniumSource

GENE_FEATURE = "Gene Expression"
CELL_COLUMNS = ("cell_id", "x_centroid", "y_centroid", "nucleus_count")


class XeniumRegionError(ValueError):
    """A region's matrix and cell table cannot be joined exactly."""


@dataclass(frozen=True)
class XeniumRegion:
    source: XeniumSource
    cell_id: np.ndarray
    cell_uid: np.ndarray
    genes: tuple[str, ...]
    counts: sparse.csr_matrix  # raw, cells x Gene Expression features
    x_um: np.ndarray
    y_um: np.ndarray
    nucleus_count: np.ndarray
    nucleus_missing: np.ndarray
    eligible: np.ndarray

    @property
    def n_cells(self) -> int:
        return len(self.cell_id)


def _decode(values) -> np.ndarray:
    return np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v)
                       for v in values], dtype=str)


def load_region(source: XeniumSource) -> XeniumRegion:
    """Read exactly the indexed bundle, preserving native and acquisition-qualified IDs.

    Every matrix barcode must occur exactly once in ``cells.parquet``; the table
    may have a different row order, so geometry is joined by ID. Zero-nucleus
    cells remain in the returned matrix for accounting and are marked ineligible.
    """
    if source.status != "ready":
        raise XeniumRegionError(f"{source.acquisition_id}: source is {source.status}, not ready")
    folder = Path(source.source_path)
    with h5py.File(folder / "cell_feature_matrix.h5", "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(v) for v in matrix["shape"][:])
        barcodes = _decode(matrix["barcodes"][:])
        names = _decode(matrix["features"]["name"][:])
        types = _decode(matrix["features"]["feature_type"][:])
        if len(shape) != 2 or shape != (len(names), len(barcodes)):
            raise XeniumRegionError(f"{source.acquisition_id}: inconsistent matrix shape")
        if len(set(barcodes)) != len(barcodes):
            raise XeniumRegionError(f"{source.acquisition_id}: duplicate matrix cell_id")
        genes = tuple(names[types == GENE_FEATURE])
        if not genes or len(set(genes)) != len(genes):
            raise XeniumRegionError(f"{source.acquisition_id}: no unique gene-expression names")
        # Xenium stores CSC(features x cells); these same buffers are CSR(cells x features).
        raw = sparse.csr_matrix((matrix["data"][:], matrix["indices"][:],
                                 matrix["indptr"][:]), shape=(shape[1], shape[0]))
        if raw.nnz and (np.any(raw.data < 0) or not np.issubdtype(raw.dtype, np.integer)):
            raise XeniumRegionError(f"{source.acquisition_id}: raw counts are not nonnegative integers")
        counts = raw[:, types == GENE_FEATURE].tocsr()

    table = pq.read_table(folder / "cells.parquet", columns=list(CELL_COLUMNS))
    cells = table.to_pydict()
    cell_ids = _decode(cells["cell_id"])
    if len(set(cell_ids)) != len(cell_ids):
        raise XeniumRegionError(f"{source.acquisition_id}: duplicate cells.parquet cell_id")
    row_of = {cell_id: n for n, cell_id in enumerate(cell_ids)}
    if len(row_of) != len(barcodes) or set(row_of) != set(barcodes):
        raise XeniumRegionError(f"{source.acquisition_id}: matrix and cells.parquet cell IDs differ")
    order = np.fromiter((row_of[cell_id] for cell_id in barcodes), dtype=np.int64,
                        count=len(barcodes))
    x = np.asarray(cells["x_centroid"], dtype=float)[order]
    y = np.asarray(cells["y_centroid"], dtype=float)[order]
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise XeniumRegionError(f"{source.acquisition_id}: non-finite cell coordinate")
    missing = np.asarray([v is None for v in cells["nucleus_count"]], dtype=bool)[order]
    nuclei = np.asarray([v if v is not None else -1 for v in cells["nucleus_count"]],
                        dtype=np.int64)[order]
    if np.any(nuclei[~missing] < 0):
        raise XeniumRegionError(f"{source.acquisition_id}: negative nucleus count")
    eligible = nuclei > 0
    cell_uid = np.asarray([f"{source.acquisition_id}:cells:{cell_id}" for cell_id in barcodes])
    return XeniumRegion(source, barcodes, cell_uid, genes, counts, x, y, nuclei, missing,
                        eligible)
