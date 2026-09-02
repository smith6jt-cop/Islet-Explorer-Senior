#!/usr/bin/env python
"""The registered re-stain, in the file formats Xenium Explorer imports.

This lane fitted a transform that maps **Xenium microns to IF microns** and measured it at a
median per-cell distance of 0.546 um (6436) and 0.630 um (6442). Xenium Explorer wants the
same correspondence expressed the other way round and in a different unit -- so the one thing
this module must get exactly right is the change of frame, and every step of it is derived
here rather than assumed.

Three sources fix the target convention, and they agree:

* 10x's file-format documentation: the morphology images are **0.2125 um/px** at full
  resolution, "the origin of the coordinate system is the upper left of the TIFF image", and
  micron coordinates "can be converted to pixel coordinates by dividing by the pixel size".
  So Xenium pixel = Xenium micron / 0.2125, with no half-pixel offset.
* 10x's image-alignment documentation: the alignment file is a **3x3 transformation matrix
  (row-major), last row [0,0,1], no header**, and "the transformation value units are in
  pixels".
* ``spatialdata_io.xenium_aligned_image`` reads that CSV with ``header=None`` and applies it
  as ``Affine(alignment, input_axes=("x","y"), output_axes=("x","y"))`` into the same global
  frame it puts Xenium geometry in by scaling microns by ``1/pixel_size``. Which settles the
  direction that prose alone leaves open: the matrix maps **imported-image pixels to Xenium
  morphology pixels**, not the reverse.

The IF side is where a half pixel can hide. ``export_dapi.py`` measured this qptiff's pyramid
to be a **box mean** (r=0.99992 against 0.9891 for the best decimation), so level-``L`` pixel
``k`` is the mean of level-0 pixels ``f*k .. f*k+f-1`` and its centre sits at ``(k + 0.5) * s``
microns for ``s = f * 0.25``. That ``s/2`` is the frame the affine was fitted in, so it belongs
in the matrix. What is *not* settled by any document is whether Explorer indexes an imported
image by pixel corner or pixel centre; the gap that leaves is ``s/2`` = **0.125 um at full
resolution**, which is a fifth of the residual this lane reports and is recorded in the
bundle's manifest rather than quietly assumed away.

Nothing here re-fits anything. It reads the transform ``run_refine.py`` chose and re-expresses
it, so a defect in this module shows up as a constant offset in Explorer rather than as a
number that looks plausible.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

#: 10x Genomics, Xenium file-format documentation: the morphology images are 0.2125 um/px at
#: full resolution and micron coordinates convert to pixels by dividing by it. Not a fit and
#: not negotiable -- the alignment file lands in Explorer's frame or it lands nowhere.
XENIUM_UM_PER_PX = 0.2125

#: 10x's conversion guidance: "tile size should be 1024 px (important for performance in
#: Xenium Explorer; smaller tile sizes may cause laggy performance)". The source qptiff is
#: tiled at 512, so re-tiling is the point of rewriting it rather than an incidental change.
EXPLORER_TILE = 1024

#: The bar ``hormone.py`` was pre-registered at, before any of it was run. Quoted here because
#: this module decides what may be called a hormone-positive cell, and that decision is not
#: allowed to move the bar in either direction.
PREREGISTERED_AUROC_BAR = 0.75


class ExplorerExportError(RuntimeError):
    """A bundle cannot be written as the thing it claims to be -- never silently degraded."""


# --------------------------------------------------------------------------- #
# The change of frame
# --------------------------------------------------------------------------- #

def alignment_matrix(affine: np.ndarray, *, if_um_per_px: float, if_origin_um: float,
                     xenium_um_per_px: float = XENIUM_UM_PER_PX) -> np.ndarray:
    """Xenium Explorer's 3x3 alignment matrix: exported IF pixels -> Xenium morphology pixels.

    ``affine`` is this lane's own 2x3, mapping **Xenium microns to IF microns**. Explorer
    needs the inverse of that, sandwiched between the two pixel/micron scalings::

        S_if   = [[s, 0, o], [0, s, o], [0, 0, 1]]   exported IF pixel -> IF micron
        A      = [[a, b, tx], [c, d, ty], [0, 0, 1]] Xenium micron     -> IF micron
        S_xen  = [[k, 0, 0], [0, k, 0], [0, 0, 1]]   Xenium pixel      -> Xenium micron

        M = S_xen^-1 @ A^-1 @ S_if

    ``o`` (``if_origin_um``) is not cosmetic. It is where exported pixel coordinate 0 sits in
    the micron frame the affine was fitted in, and ``export_dapi.origin_um`` derives it from
    the *measured* pyramid rule. Dropping it puts every cell a constant ``o`` off in both axes.
    """
    A = np.asarray(affine, dtype=np.float64)
    if A.shape != (2, 3):
        raise ExplorerExportError(f"expected a 2x3 micron affine, got shape {A.shape}")
    if not np.isfinite(A).all():
        raise ExplorerExportError("the affine carries non-finite entries")
    if abs(float(np.linalg.det(A[:, :2]))) < 1e-12:
        raise ExplorerExportError("the affine's linear part is singular and cannot be inverted")
    for name, v in (("if_um_per_px", if_um_per_px), ("xenium_um_per_px", xenium_um_per_px)):
        if not float(v) > 0.0:
            raise ExplorerExportError(f"{name} must be positive, got {v!r}")

    A_h = np.vstack([A, [0.0, 0.0, 1.0]])
    s, o, k = float(if_um_per_px), float(if_origin_um), float(xenium_um_per_px)
    S_if = np.array([[s, 0.0, o], [0.0, s, o], [0.0, 0.0, 1.0]])
    S_xen_inv = np.array([[1.0 / k, 0.0, 0.0], [0.0, 1.0 / k, 0.0], [0.0, 0.0, 1.0]])

    M = S_xen_inv @ np.linalg.inv(A_h) @ S_if
    M[2] = [0.0, 0.0, 1.0]          # exact, rather than 1e-17 away from it
    return M


def matrix_round_trip_um(M: np.ndarray, affine: np.ndarray, *, if_um_per_px: float,
                         if_origin_um: float, xenium_um_per_px: float = XENIUM_UM_PER_PX,
                         xy_xenium_um: np.ndarray) -> float:
    """Worst-case microns by which ``M`` disagrees with ``affine`` on the given points.

    The matrix is a re-expression, so this has no fitting freedom in it: it is a check that
    the composition above was assembled correctly, and it is written into the manifest so a
    reader does not have to take the algebra on trust. Points are pushed Xenium micron ->
    IF micron -> IF pixel -> (via M) Xenium pixel -> Xenium micron and compared with where
    they started.
    """
    xy = np.asarray(xy_xenium_um, dtype=np.float64).reshape(-1, 2)
    A = np.asarray(affine, dtype=np.float64)
    if_um = xy @ A[:, :2].T + A[:, 2]
    if_px = (if_um - float(if_origin_um)) / float(if_um_per_px)
    xen_px = np.column_stack([if_px, np.ones(len(if_px))]) @ np.asarray(M).T
    back = xen_px[:, :2] * float(xenium_um_per_px)
    return float(np.max(np.linalg.norm(back - xy, axis=1)))


def write_alignment_csv(path: Path, M: np.ndarray) -> Path:
    """The matrix as Explorer reads it: three rows, no header, last row exactly ``0,0,1``.

    Written with ``repr`` rather than a fixed number of decimals. The linear part here is
    ~1.18 and the translation ~1e4 pixels, so rounding the matrix to six places would round
    the translation to a tenth of a pixel while leaving the scale exact -- a formatting
    choice that would show up as a real displacement in the viewer.
    """
    M = np.asarray(M, dtype=np.float64)
    if M.shape != (3, 3):
        raise ExplorerExportError(f"an alignment matrix is 3x3, got {M.shape}")
    if not np.allclose(M[2], [0.0, 0.0, 1.0]):
        raise ExplorerExportError(f"the last row must be [0,0,1], got {M[2].tolist()}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        for row in M:
            w.writerow([repr(float(v)) for v in row])
    return path


# --------------------------------------------------------------------------- #
# The image
# --------------------------------------------------------------------------- #

def qptiff_level_shapes(qptiff: Path) -> list[tuple[int, int]]:
    """``(height, width)`` of every pyramid level of the qptiff's first series."""
    import tifffile

    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        return [(int(lv.shape[-2]), int(lv.shape[-1]))
                for lv in getattr(series, "levels", [series])]


def integer_factors(shapes: list[tuple[int, int]]) -> list[int]:
    """Downsample factor of each level, refused unless every one is an exact integer.

    ``shapes[0][0] / shapes[lvl][0]`` is not the factor -- the levels are floor-halved, so the
    ratio can be 4.000176 where the factor is 4, a scale error worth 0.6 um at the far corner
    of a section. ``export_dapi.pyramid_factor`` makes the same point for the DAPI export and
    this is the same rule applied to every level at once, because the exported OME-TIFF
    carries all of them and Explorer will navigate between them.
    """
    import export_dapi as ED

    return [ED.pyramid_factor(shapes, lvl) for lvl in range(len(shapes))]


def _tile_stream(zstore, channels: list[int], height: int, width: int, dtype,
                 tile: int = EXPLORER_TILE):
    """Tiles in the order ``TiffWriter`` consumes them: per channel, per row, per column.

    Streamed rather than materialised. A full-resolution channel of these sections is 2.9 GB,
    and the point of writing tiles is that the peak is one tile rather than the slide.
    """
    for c in channels:
        for y0 in range(0, height, tile):
            for x0 in range(0, width, tile):
                patch = np.asarray(zstore[c, y0:y0 + tile, x0:x0 + tile])
                if patch.shape != (tile, tile):
                    pad = np.zeros((tile, tile), dtype=dtype)
                    pad[:patch.shape[0], :patch.shape[1]] = patch
                    patch = pad
                yield patch


def write_explorer_ome_tiff(qptiff: Path, out_path: Path, channels: list[int],
                            channel_names: list[str], *, base_um_per_px: float,
                            compression: str = "zlib", tile: int = EXPLORER_TILE,
                            max_levels: int | None = None, workers: int = 8) -> dict:
    """The selected qptiff channels as a pyramidal, tiled OME-TIFF Explorer can import.

    10x's requirement is "pyramidal, tiled, OME-TIFF"; this writes exactly that, and takes
    each pyramid level **from the qptiff's own levels** rather than downsampling. That is not
    an optimisation. ``export_dapi.py`` measured this pyramid to be a box mean, and the
    alignment matrix is derived from that same rule -- so re-deriving the levels with some
    other filter would put the matrix and the image on subtly different grids.

    ``compression`` defaults to ``zlib``: lossless, universally readable, and it does not put
    a lossy codec between the colleague and a hormone channel they are being asked to judge.
    10x's own morphology files are JPEG-2000, which ``imagecodecs`` can also write; it is a
    flag rather than the default for that reason.
    """
    import tifffile

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if len(channels) != len(channel_names):
        raise ExplorerExportError("every exported channel needs a name")
    if not channels:
        raise ExplorerExportError("no channels selected; refusing to write an empty image")

    shapes = qptiff_level_shapes(qptiff)
    factors = integer_factors(shapes)
    n_levels = len(shapes) if max_levels is None else min(len(shapes), int(max_levels))

    written = []
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        dtype = series.dtype
        with tifffile.TiffWriter(str(out_path), bigtiff=True) as tw:
            for lvl in range(n_levels):
                h, w = shapes[lvl]
                store = series.aszarr(level=lvl)
                try:
                    import zarr
                    z = zarr.open(store, mode="r")
                    shape = (len(channels), h, w)
                    # Compression is the wall clock here -- a full-resolution channel is
                    # 2.9 GB and there are three of them -- and it is per-tile, so it
                    # parallelises exactly. The tile ORDER is still serial, so the file is
                    # byte-identical whatever this is set to.
                    opts = dict(tile=(tile, tile), compression=compression,
                                photometric="minisblack", maxworkers=int(workers))
                    if lvl == 0:
                        um = float(base_um_per_px)
                        tw.write(_tile_stream(z, channels, h, w, dtype, tile=tile),
                                 shape=shape, dtype=dtype, subifds=n_levels - 1,
                                 metadata={"axes": "CYX",
                                           "Channel": {"Name": list(channel_names)},
                                           "PhysicalSizeX": um, "PhysicalSizeY": um,
                                           "PhysicalSizeXUnit": "µm",
                                           "PhysicalSizeYUnit": "µm"},
                                 **opts)
                    else:
                        tw.write(_tile_stream(z, channels, h, w, dtype, tile=tile),
                                 shape=shape, dtype=dtype, subfiletype=1, **opts)
                finally:
                    store.close()
                written.append({"level": lvl, "height": h, "width": w,
                                "factor": int(factors[lvl]),
                                "um_per_px": float(base_um_per_px) * int(factors[lvl])})
                print(f"[explorer]   level {lvl}: {h} x {w} x {len(channels)}", flush=True)

    return {"path": str(out_path), "dtype": str(dtype), "tile": int(tile),
            "compression": compression, "channels": list(channel_names),
            "channel_indices": [int(c) for c in channels], "levels": written,
            "bytes": int(out_path.stat().st_size)}


# --------------------------------------------------------------------------- #
# The annotations
# --------------------------------------------------------------------------- #

#: QuPath's own packed-ARGB integers, so the classes arrive in Explorer coloured the way the
#: maintainer drew them rather than in whatever order an importer assigns.
QUPATH_COLOR = {"Tissue": -16744448, "Islet": -3670016}


def warp_geoms_to_xenium(geoms, transform):
    """IF-micron shapely geometries -> Xenium-micron shapely geometries.

    ``inverse_apply`` rather than ``apply``: this lane's transform is fitted Xenium -> IF, and
    the annotations were drawn on the IF. Getting this backwards does not error -- it produces
    polygons roughly 2 * (714, 282) um from where they belong on 6436, which looks like a
    registration failure rather than like a sign error.
    """
    import shapely

    return [shapely.transform(g, lambda a: transform.inverse_apply(a)) for g in geoms]


def _rings(geom) -> list[list[list[float]]]:
    """GeoJSON coordinate rings for a Polygon, exterior first."""
    import shapely

    if geom.geom_type != "Polygon":
        raise ExplorerExportError(f"expected a Polygon, got {geom.geom_type}")
    out = [[[float(x), float(y)] for x, y in shapely.get_coordinates(geom.exterior)]]
    for hole in geom.interiors:
        out.append([[float(x), float(y)] for x, y in shapely.get_coordinates(hole)])
    return out


def annotations_geojson(named_geoms: dict[str, list], path: Path) -> dict:
    """Class name -> geometries, written as a QuPath-shaped GeoJSON FeatureCollection.

    Coordinates are **Xenium microns**. Explorer's annotation import asks which of microns or
    pixels a file is in, so the frame is stated in the filename as well as here; a file that
    is silently in the wrong one imports without complaint and lands 4.7x off.

    ``MultiPolygon`` is split into one feature per part. Explorer's importer takes polygons,
    and the ``Tissue`` class arrives from QuPath as a single MultiPolygon of every fragment on
    the slide, which as one feature would be one annotation covering the whole section.
    """
    features = []
    for cls, geoms in named_geoms.items():
        parts = []
        for g in geoms:
            if g.is_empty:
                continue
            parts.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
        for i, g in enumerate(parts):
            features.append({
                "type": "Feature",
                "id": f"{cls}-{i + 1}",
                "geometry": {"type": "Polygon", "coordinates": _rings(g)},
                "properties": {"objectType": "annotation",
                               "name": f"{cls} {i + 1}",
                               "classification": {"name": cls,
                                                  "colorRGB": QUPATH_COLOR.get(cls, -1)}},
            })
    doc = {"type": "FeatureCollection", "features": features}
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))
    return {"path": str(path), "n_features": len(features),
            "per_class": {c: sum(1 for f in features
                                 if f["properties"]["classification"]["name"] == c)
                          for c in named_geoms}}


def annotations_csv(named_geoms: dict[str, list], path: Path) -> dict:
    """The same polygons as ``Selection,X,Y``, for Explorer versions older than the GeoJSON path.

    One row per vertex, exteriors only -- the CSV form has no way to express a hole, and
    silently dropping one is better stated than discovered. Coordinates are Xenium microns.
    """
    import shapely

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = n_sel = n_holes = 0
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Selection", "X", "Y"])
        for cls, geoms in named_geoms.items():
            parts = []
            for g in geoms:
                if g.is_empty:
                    continue
                parts.extend(list(g.geoms) if g.geom_type == "MultiPolygon" else [g])
            for i, g in enumerate(parts):
                n_sel += 1
                n_holes += len(g.interiors)
                for x, y in shapely.get_coordinates(g.exterior):
                    w.writerow([f"{cls} {i + 1}", repr(float(x)), repr(float(y))])
                    n_rows += 1
    return {"path": str(path), "n_selections": n_sel, "n_vertices": n_rows,
            "n_holes_dropped": n_holes}


# --------------------------------------------------------------------------- #
# The cell groups
# --------------------------------------------------------------------------- #

def write_cell_groups(path: Path, cell_ids, groups, colors: dict[str, str]) -> dict:
    """Explorer's cell-group CSV: ``cell_id,group,color``, header in that exact order.

    ``color`` is optional from Explorer v4.1 and is written here because the alternative is a
    palette assigned by import order, which would give the same label a different colour in
    two donors. A cell whose group has no colour is refused rather than written blank.
    """
    import pandas as pd

    cell_ids = np.asarray(cell_ids)
    groups = np.asarray(groups)
    if len(cell_ids) != len(groups):
        raise ExplorerExportError(f"{len(cell_ids)} cells against {len(groups)} labels")
    missing = sorted(set(map(str, np.unique(groups))) - set(colors))
    if missing:
        raise ExplorerExportError(f"no colour recorded for group(s) {missing}")

    df = pd.DataFrame({"cell_id": cell_ids, "group": groups,
                       "color": [colors[str(g)] for g in groups]})
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    counts = df.group.value_counts().to_dict()
    return {"path": str(path), "n_cells": int(len(df)),
            "groups": {str(k): int(v) for k, v in sorted(counts.items())}}


def two_component_threshold(values: np.ndarray, *, max_iter: int = 200,
                            tol: float = 1e-8, seed: int = 0) -> dict:
    """Where one donor's own log intensity stops looking like its lower component.

    A two-component Gaussian mixture fitted by EM on ``log1p(value)``, with the threshold at
    the crossing where the posterior for the upper component reaches 0.5.

    **What makes this admissible under the repository's own rule.** The rule forbids tuning
    toward an expected marker frequency, cross-donor uniformity or disease ordering. Nothing
    here is told what fraction should come out positive: the mixing weight is *estimated from
    this donor's own distribution* and then reported, so the positive fraction is an output to
    be read, never an input to be hit. Every fitted parameter is written into the manifest for
    the same reason -- a threshold nobody can see the derivation of is a threshold nobody can
    disagree with.

    It is still a model, and it is only ever applied to a channel whose AUROC has already
    cleared the pre-registered bar. On a channel that has not, this number would be a
    confident split of noise.
    """
    v = np.asarray(values, dtype=np.float64)
    x = np.log1p(v[np.isfinite(v)])
    if x.size < 100:
        raise ExplorerExportError(f"only {x.size} finite values; refusing to fit a mixture")

    lo, hi = np.percentile(x, [25.0, 75.0])
    mu = np.array([lo, hi], dtype=np.float64)
    if not mu[1] > mu[0]:
        raise ExplorerExportError("the channel has no spread between its quartiles")
    sd = np.full(2, max(float(x.std()), 1e-6))
    w = np.array([0.5, 0.5])

    ll_prev = -np.inf
    for it in range(int(max_iter)):
        p = w * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
        tot = p.sum(axis=1)
        tot[tot <= 0] = np.finfo(float).tiny
        r = p / tot[:, None]
        nk = r.sum(axis=0)
        nk[nk <= 0] = np.finfo(float).tiny
        w = nk / len(x)
        mu = (r * x[:, None]).sum(axis=0) / nk
        sd = np.sqrt((r * (x[:, None] - mu) ** 2).sum(axis=0) / nk)
        sd = np.maximum(sd, 1e-6)
        ll = float(np.log(tot).sum())
        if abs(ll - ll_prev) < tol * max(1.0, abs(ll)):
            break
        ll_prev = ll

    order = np.argsort(mu)
    mu, sd, w = mu[order], sd[order], w[order]

    grid = np.linspace(mu[0], mu[1], 20_001)
    pu = w[1] * np.exp(-0.5 * ((grid - mu[1]) / sd[1]) ** 2) / sd[1]
    pl = w[0] * np.exp(-0.5 * ((grid - mu[0]) / sd[0]) ** 2) / sd[0]
    cross = np.nonzero(pu >= pl)[0]
    if cross.size == 0:
        raise ExplorerExportError(
            "the upper component never outweighs the lower one between their means; "
            "this channel is not two-component and no threshold is defined for it")
    thr_log = float(grid[cross[0]])

    return {"threshold": float(np.expm1(thr_log)), "threshold_log1p": thr_log,
            "means_log1p": mu.tolist(), "sds_log1p": sd.tolist(), "weights": w.tolist(),
            "n_fitted": int(x.size), "iterations": int(it + 1),
            "log_likelihood": float(ll),
            "separation_sds": float((mu[1] - mu[0]) / np.sqrt((sd ** 2).mean()))}


def decile_bands(values: np.ndarray) -> tuple[np.ndarray, dict]:
    """Per-donor deciles of a channel, as ordered labels -- a display, not a call.

    Cells whose value is unavailable (off the imported image) are labelled ``Unavailable``
    with the reason recorded beside them, rather than pooled into a bottom band where they
    would read as "measured, and negative".
    """
    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    labels = np.full(v.shape, "Unavailable", dtype=object)
    if ok.sum() < 10:
        return labels, {"n_available": int(ok.sum()), "edges": []}
    edges = np.percentile(v[ok], np.arange(10, 100, 10))
    idx = np.searchsorted(edges, v[ok], side="right")
    labels[ok] = [f"d{i + 1:02d}" for i in idx]
    return labels, {"n_available": int(ok.sum()),
                    "edges": [float(e) for e in edges]}
