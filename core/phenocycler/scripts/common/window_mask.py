"""Windowed rasterization: label masks for the seven donors that have none.

Thirteen of the twenty cohort donors have an archived whole-slide ``int32`` mask from the
2026-06-19 REDSEA run.  The other seven -- every batch-2 donor, which is to say every panel
variant except cd38-59ch -- have none, because their GeoJSONs were regenerated after those
masks were made.  A gate that ran only on the thirteen would be a gate on one panel.

``redsea.rasterize_mask`` can build them, but it allocates a whole-slide ``int32`` raster:
~2.4 GB for a 23,760 x 24,960 slide, and we need eight 1024-px cores per donor.  So this
module streams each GeoJSON **once** and draws only the polygons that intersect one of the
requested windows, reusing ``redsea._draw_geometry`` unchanged so the polygon fill, the
interior-ring carving and the MultiPolygon handling are the same code that built the
archived masks.

**The label invariant is preserved exactly**: ``label = feature_index + 1`` over the whole
stream, so a window's labels are global and sparse, ``object_ids[label - 1]`` is the cell's
QuPath UUID, and nothing renumbers.  That is what lets a window join to the run's expression
parquet on ``object_id``.

Before this is trusted on a donor that has no archived mask, ``verify_against_archive``
must reproduce an archived one **bit for bit** over many windows.  A mismatch stops the
lane; it is not a tolerance to widen.
"""
from __future__ import annotations

import pathlib
import sys
from dataclasses import dataclass

import numpy as np

_PKG = pathlib.Path(__file__).resolve().parents[2]
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))

from phenocycler.artifacts import geometry_present, iter_qupath_geojson_features  # noqa: E402
from phenocycler.redsea import _draw_geometry  # noqa: E402

__all__ = ["Window", "rasterize_windows", "verify_against_archive"]


@dataclass(frozen=True)
class Window:
    y0: int
    x0: int
    size: int

    @property
    def y1(self) -> int:
        return self.y0 + self.size

    @property
    def x1(self) -> int:
        return self.x0 + self.size


def _bbox(geom):
    """(x0, y0, x1, y1) over every ring, without materializing the whole coordinate list."""
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    xs0 = ys0 = np.inf
    xs1 = ys1 = -np.inf
    for poly in polys:
        ext = np.asarray(poly[0], dtype=np.float64)
        if ext.shape[0] < 3:
            continue
        xs0 = min(xs0, ext[:, 0].min()); xs1 = max(xs1, ext[:, 0].max())
        ys0 = min(ys0, ext[:, 1].min()); ys1 = max(ys1, ext[:, 1].max())
    return xs0, ys0, xs1, ys1


def _shift(geom, dx, dy):
    """A copy translated into window coordinates, so _draw_geometry needs no change."""
    def ring(r):
        a = np.asarray(r, dtype=np.float64).copy()
        a[:, 0] -= dx
        a[:, 1] -= dy
        return a.tolist()
    if geom["type"] == "MultiPolygon":
        coords = [[ring(r) for r in poly] for poly in geom["coordinates"]]
    else:
        coords = [ring(r) for r in geom["coordinates"]]
    return {"type": geom["type"], "coordinates": coords}


def rasterize_windows(geojson, windows, *, with_nucleus: bool = False,
                      downsample: float = 1.0, progress=None):
    """One pass over the GeoJSON; one int32 raster per window, with GLOBAL labels.

    Returns ``(masks, object_ids[, nucleus_masks])``.  ``object_ids`` is the full feature
    list, so ``object_ids[label - 1]`` holds for every window.
    """
    windows = [w if isinstance(w, Window) else Window(*w) for w in windows]
    if not windows:
        raise ValueError("no windows requested -- an empty rasterization is not a result")
    if any(w.size <= 0 for w in windows):
        raise ValueError("a window has non-positive size")
    s = 1.0 / downsample
    masks = [np.zeros((w.size, w.size), dtype=np.int32) for w in windows]
    nucs = [np.zeros((w.size, w.size), dtype=np.int32) for w in windows] if with_nucleus else None
    oids = []

    for i, ft in enumerate(iter_qupath_geojson_features(pathlib.Path(geojson))):
        oids.append(ft.get("id"))
        label = i + 1
        geom = ft.get("geometry")
        if not geom:
            continue
        bx0, by0, bx1, by1 = _bbox(geom)
        if not np.isfinite(bx0):
            continue
        for k, w in enumerate(windows):
            # bbox test in the same scaled space the drawing uses
            if bx1 * s < w.x0 or bx0 * s > w.x1 or by1 * s < w.y0 or by0 * s > w.y1:
                continue
            shape = (w.size, w.size)
            _draw_geometry(masks[k], _shift(geom, w.x0 / s, w.y0 / s), label, shape, s)
            if with_nucleus:
                ng = ft.get("nucleusGeometry")
                if geometry_present(ng):
                    _draw_geometry(nucs[k], _shift(ng, w.x0 / s, w.y0 / s), label, shape, s)
        if progress and (i + 1) % 200000 == 0:
            progress(i + 1)

    if with_nucleus:
        # Same rule as rasterize_mask: a nucleus pixel the cell raster no longer claims is
        # dropped, so nucleus is a strict subset of cell per label.
        for m, nz in zip(masks, nucs):
            np.putmask(nz, nz != m, 0)
        return masks, oids, nucs
    return masks, oids


def verify_against_archive(geojson, archive_tif, windows, *, downsample: float = 1.0):
    """Bit-identity against an archived whole-slide mask.  Returns (all_ok, per-window rows).

    This is the gate on the rasterizer itself.  It runs on donors that HAVE an archived mask
    so the four that do not can be trusted; a mismatch stops the lane rather than being
    absorbed into a tolerance.
    """
    import tifffile

    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from tiff_window import read_window

    windows = [w if isinstance(w, Window) else Window(*w) for w in windows]
    fresh, _ = rasterize_windows(geojson, windows, downsample=downsample)
    rows, ok_all = [], True
    with tifffile.TiffFile(archive_tif) as tf:
        page = tf.pages[0]
        for w, got in zip(windows, fresh):
            want = read_window(page, w.y0, w.y1, w.x0, w.x1).astype(np.int32)
            same = np.array_equal(want, got)
            ok_all &= same
            diff = int((want != got).sum())
            rows.append({"y0": w.y0, "x0": w.x0, "size": w.size, "identical": same,
                         "pixels_differing": diff,
                         "labels_archive": int(np.unique(want[want > 0]).size),
                         "labels_fresh": int(np.unique(got[got > 0]).size)})
    return ok_all, rows
