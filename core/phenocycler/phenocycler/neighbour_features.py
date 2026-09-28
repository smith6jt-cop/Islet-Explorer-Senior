"""Per-cell spatial-neighbour features, from centroids, for the v3 typing lanes.

Three things this module computes, all from a k-nearest-neighbour query built once per donor with
``scipy.spatial.cKDTree`` and evaluated in row chunks so a donor of millions of cells never
materialises an all-pairs distance matrix:

* :func:`max_neighbour_intensity` -- the highest value of a requested marker among a cell's own
  neighbours (for example, to flag hormone bleed-through onto a non-endocrine cell sitting next to
  a bright endocrine one);
* :func:`neighbour_fraction` -- the fraction of a cell's neighbours that are "positive" under a
  caller-supplied boolean/nullable-boolean mask (a marker call, a class, a contact criterion --
  this module does not define positivity, it only aggregates it over the neighbourhood);
* :func:`distance_to_nearest_exceeding` -- the distance from a cell to the nearest OTHER cell
  whose value on a marker exceeds a caller-supplied threshold (the bleed-through distance itself,
  e.g. distance to the nearest GCG+ cell). The threshold is never computed here from a target
  frequency or any other data-derived shortcut -- it is the caller's reference level, passed in.

Every function follows the same missing-data rule as the rest of this workflow: a cell with no
neighbour within ``radius_um`` (or, for the exceeding-threshold distance, a donor with no cell at
all above the threshold) reports ``NaN`` -- unavailable, never zero and never treated as infinity.

**Polygons**: this module works from centroids only. A caller with real cell polygons and a true
touching/adjacency graph can still use :func:`neighbour_fraction` -- pass an adjacency-derived
positive mask (1.0 for a touching neighbour, 0.0 otherwise) alongside a query built with a radius
generous enough to include every geometrically adjacent centroid; this module does not itself
rasterize or intersect polygons.

Not wired into any stage, ``build_layers.py`` node, or the Makefile. Additive capability only.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .artifacts import ContractError


@dataclass(frozen=True)
class NeighbourQuery:
    """The k-nearest-neighbour result for one donor's centroids, self excluded.

    ``indices[i, j]`` is the global row index of cell ``i``'s ``j``-th nearest neighbour within
    ``radius_um``, or ``-1`` if cell ``i`` has fewer than ``j + 1`` such neighbours.
    ``distances[i, j]`` is the matching Euclidean distance in the same units as the input
    centroids, or ``NaN`` where ``indices`` is ``-1``.
    """

    indices: np.ndarray
    distances: np.ndarray
    k: int
    radius_um: float


def _validate_coords(coords: np.ndarray) -> np.ndarray:
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ContractError("coords must be an (n, 2) array of centroid x, y micrometres")
    if coords.shape[0] == 0:
        raise ContractError("coords must contain at least one cell")
    return coords


def query_neighbours(
    coords: np.ndarray,
    *,
    k: int,
    radius_um: float,
    chunk_size: int = 100_000,
) -> NeighbourQuery:
    """The k nearest OTHER cells within ``radius_um`` of every cell, chunked for memory.

    ``k`` bounds how many neighbours are kept per cell; ``radius_um`` bounds how far a neighbour
    may be. A cell can have fewer than ``k`` valid neighbours (padded with ``-1`` / ``NaN``); it
    never has more.
    """
    coords = _validate_coords(coords)
    n = coords.shape[0]
    if k < 1:
        raise ContractError("k must be >= 1")
    if not (radius_um > 0 and np.isfinite(radius_um)):
        raise ContractError("radius_um must be finite and > 0")
    if chunk_size < 1:
        raise ContractError("chunk_size must be >= 1")

    tree = cKDTree(coords)
    k_query = min(k + 1, n)  # +1 to make room for self, which is always found at distance 0.
    width = max(k_query - 1, 0)
    out_idx = np.full((n, k), -1, dtype=np.int64)
    out_dist = np.full((n, k), np.nan, dtype=np.float64)

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        d, idx = tree.query(
            coords[start:end], k=k_query, distance_upper_bound=radius_um, workers=-1
        )
        if k_query == 1:
            d = d[:, None]
            idx = idx[:, None]
        row_ids = np.arange(start, end)[:, None]
        self_mask = idx == row_ids
        # Stable sort on "is this column self" (False before True) drops exactly the one self
        # occurrence per row while preserving the existing distance order among the rest.
        order = np.argsort(self_mask, axis=1, kind="stable")
        idx_sorted = np.take_along_axis(idx, order, axis=1)[:, :width]
        d_sorted = np.take_along_axis(d, order, axis=1)[:, :width]
        invalid = (idx_sorted == n) | ~np.isfinite(d_sorted)
        idx_sorted = np.where(invalid, -1, idx_sorted)
        d_sorted = np.where(invalid, np.nan, d_sorted)
        out_idx[start:end, :width] = idx_sorted
        out_dist[start:end, :width] = d_sorted

    return NeighbourQuery(indices=out_idx, distances=out_dist, k=k, radius_um=radius_um)


def max_neighbour_intensity(query: NeighbourQuery, values: np.ndarray) -> np.ndarray:
    """Per cell, the highest ``values`` among its own neighbours; ``NaN`` with none available."""
    values = np.asarray(values, dtype=float)
    n = query.indices.shape[0]
    if len(values) != n:
        raise ContractError("values must have exactly one entry per cell in the query")

    valid = query.indices >= 0
    gathered = np.full(query.indices.shape, np.nan, dtype=float)
    gathered[valid] = values[query.indices[valid]]
    out = np.full(n, np.nan, dtype=float)
    any_finite = np.isfinite(gathered).any(axis=1)
    out[any_finite] = np.nanmax(gathered[any_finite], axis=1)
    return out


def neighbour_fraction(query: NeighbourQuery, positive: np.ndarray | pd.Series) -> np.ndarray:
    """Per cell, the fraction of its neighbours flagged ``True`` in ``positive``.

    ``positive`` may carry ``NA``/``NaN`` for a neighbour whose own measurement was unavailable;
    such a neighbour is excluded from both the numerator and the denominator rather than counted
    as negative. A cell with zero neighbours whose ``positive`` value is itself known reports
    ``NaN`` -- the fraction is undefined, not zero.
    """
    n = query.indices.shape[0]
    as_float = pd.array(positive, dtype="boolean").to_numpy(dtype=float, na_value=np.nan)
    if len(as_float) != n:
        raise ContractError("positive must have exactly one entry per cell in the query")

    valid = query.indices >= 0
    gathered = np.full(query.indices.shape, np.nan, dtype=float)
    gathered[valid] = as_float[query.indices[valid]]
    known = np.isfinite(gathered)
    denom = known.sum(axis=1)
    numer = np.where(known, gathered, 0.0).sum(axis=1)
    out = np.full(n, np.nan, dtype=float)
    has_denom = denom > 0
    out[has_denom] = numer[has_denom] / denom[has_denom]
    return out


def distance_to_nearest_exceeding(
    coords: np.ndarray,
    values: np.ndarray,
    threshold: float,
    *,
    exclude_self: bool = True,
    chunk_size: int = 100_000,
) -> np.ndarray:
    """Per cell, the distance to the nearest cell whose ``values`` exceed ``threshold``.

    ``threshold`` is taken exactly as given -- this function fits nothing and infers nothing from
    the data; it is the caller's reference level (e.g. a hormone's donor reference level). A
    donor in which no cell exceeds the threshold reports ``NaN`` for every cell, and with
    ``exclude_self=True`` (the default) a cell that itself exceeds the threshold reports its
    distance to the nearest OTHER such cell rather than 0.
    """
    coords = _validate_coords(coords)
    values = np.asarray(values, dtype=float)
    n = coords.shape[0]
    if len(values) != n:
        raise ContractError("values must have exactly one entry per cell in coords")
    if not np.isfinite(threshold):
        raise ContractError("threshold must be finite")
    if chunk_size < 1:
        raise ContractError("chunk_size must be >= 1")

    out = np.full(n, np.nan, dtype=float)
    subset_global = np.flatnonzero(np.isfinite(values) & (values > threshold))
    if subset_global.size == 0:
        return out

    subset_tree = cKDTree(coords[subset_global])
    k_query = 2 if (exclude_self and subset_global.size > 1) else 1

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        d, idx = subset_tree.query(coords[start:end], k=k_query, workers=-1)
        if k_query == 1:
            d = d[:, None]
            idx = idx[:, None]
        row_global = subset_global[idx]
        if exclude_self:
            row_ids = np.arange(start, end)[:, None]
            self_mask = row_global == row_ids
            order = np.argsort(self_mask, axis=1, kind="stable")
            row_global = np.take_along_axis(row_global, order, axis=1)
            d = np.take_along_axis(d, order, axis=1)
            is_still_self = row_global[:, 0] == np.arange(start, end)
            out[start:end] = np.where(is_still_self, np.nan, d[:, 0])
        else:
            out[start:end] = d[:, 0]

    return out


__all__ = [
    "NeighbourQuery",
    "distance_to_nearest_exceeding",
    "max_neighbour_intensity",
    "neighbour_fraction",
    "query_neighbours",
]
