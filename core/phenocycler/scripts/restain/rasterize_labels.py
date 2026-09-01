"""Shapely polygons in microns -> an int32 label image on a :class:`GridSpec`.

``geo.rasterize`` already puts shapely-in-microns onto a grid, but through PIL mode ``"1"``,
so it is binary: adjacent islets merge irrecoverably and per-cell identity is impossible.
``scripts/dct/window_mask.rasterize_windows`` keeps labels but reads a QuPath GeoJSON path in
qptiff pixels.  Neither is the thing a per-cell measurement needs, so this is the third
rasteriser -- and the only one that is both labelled and in micron space.

**Pixel centres, not edges.**  A pixel is inside a polygon when its CENTRE is, so micron
``x`` maps to column ``(x - origin) / pixel_um - 0.5``.  ``geo.rasterize`` (PIL, used for the
registration tissue mask) and ``redsea._draw_geometry`` use the edge convention and are
deliberately left alone: their half-pixel is a constant translation that the fitted affine
absorbs, and changing them would move numbers that have nothing to do with this lane.

**A segmentation tiles without gaps**, so at any finite resolution two neighbours compete for
the pixel on their shared boundary.  Later draws win, and the order is smallest-last
(``area_desc``) so the cell with less to lose keeps the pixel.  ``n_contested`` counts how
often it happened rather than asserting it is negligible.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phenocycler.integration.rasterize import GridSpec


class RasterError(RuntimeError):
    """Geometry cannot be rasterised onto the grid it was given."""


@dataclass(frozen=True)
class LabelRaster:
    """``labels`` holds ``i + 1`` for ``geoms[i]`` and 0 for background."""

    labels: np.ndarray
    grid: GridSpec
    drawn_px: np.ndarray          # pixels the polygon covered, before anything overdrew it
    owned_px: np.ndarray          # pixels it still holds in `labels`
    n_contested: int              # pixels claimed by more than one polygon

    @property
    def n_labels(self) -> int:
        return int(self.drawn_px.size)


def _rings_in_pixel_centres(poly, grid: GridSpec):
    """Exterior and interior rings as (rows, cols) in pixel-centre coordinates."""
    p, ox, oy = grid.pixel_um, grid.origin_x_um, grid.origin_y_um

    def conv(coords):
        a = np.asarray(coords, dtype=np.float64)
        return (a[:, 1] - oy) / p - 0.5, (a[:, 0] - ox) / p - 0.5

    parts = [poly] if poly.geom_type == "Polygon" else list(poly.geoms)
    for part in parts:
        if part.is_empty:
            continue
        yield conv(part.exterior.coords), [conv(r.coords) for r in part.interiors]


def rasterize_labels(geoms, grid: GridSpec, *, order: str = "area_desc",
                     dtype=np.int32) -> LabelRaster:
    """Label image for ``geoms`` (shapely, microns) on ``grid``.

    ``geoms`` may contain ``Polygon`` or ``MultiPolygon``; every part of a multipart geometry
    draws at the same label and its pixels sum into one entry.
    """
    from skimage.draw import polygon as draw_polygon

    geoms = list(geoms)
    n = len(geoms)
    if n and np.iinfo(dtype).max < n:
        raise RasterError(f"{n} polygons will not fit in {np.dtype(dtype).name}")

    labels = np.zeros(grid.shape, dtype=dtype)
    drawn = np.zeros(n, dtype=np.int64)
    contested = 0

    if order == "area_desc":
        seq = sorted(range(n), key=lambda i: -float(getattr(geoms[i], "area", 0.0)))
    elif order == "input":
        seq = list(range(n))
    else:
        raise RasterError(f"unknown order {order!r}")

    # Flat indices rather than a full-grid boolean per polygon. The mask-per-polygon form is
    # O(n_polygons * grid_size): correct, but on a whole section (958,744 nuclei over a
    # 10260x5040 grid) that is ~50 MB of allocation per polygon and never finishes.
    W = grid.width
    flat_labels = labels.ravel()
    for i in seq:
        g = geoms[i]
        if g is None or getattr(g, "is_empty", True):
            continue
        parts, holes = [], []
        for (rr_ext, cc_ext), rings in _rings_in_pixel_centres(g, grid):
            rr, cc = draw_polygon(rr_ext, cc_ext, shape=grid.shape)
            if rr.size:
                parts.append(rr.astype(np.int64) * W + cc.astype(np.int64))
            for hole_rr, hole_cc in rings:
                hr, hc = draw_polygon(hole_rr, hole_cc, shape=grid.shape)
                if hr.size:
                    holes.append(hr.astype(np.int64) * W + hc.astype(np.int64))
        if not parts:
            continue
        flat = np.unique(np.concatenate(parts))       # multipart parts may touch
        if holes:
            flat = np.setdiff1d(flat, np.concatenate(holes), assume_unique=True)
        if flat.size == 0:
            continue
        drawn[i] = int(flat.size)
        contested += int(np.count_nonzero(flat_labels[flat]))
        flat_labels[flat] = i + 1

    owned = np.bincount(labels.ravel(), minlength=n + 1)[1:].astype(np.int64)
    return LabelRaster(labels=labels, grid=grid, drawn_px=drawn, owned_px=owned,
                       n_contested=contested)


def rasterize_binary(geoms, grid: GridSpec) -> np.ndarray:
    """Union of ``geoms`` as a boolean mask -- the alignment gate's template."""
    return rasterize_labels(geoms, grid, order="input").labels > 0
