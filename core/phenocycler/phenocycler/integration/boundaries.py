"""
Cell polygons, in microns, from both platforms — the loader ``match_by_iou`` was missing.

:func:`sameslide.match_by_iou` has been implemented and correct since the integration layer
was written, and unused, because nothing could hand it polygons. This is that piece.

Why bother when centroids already work: **centroid distance and cell contact are different
questions.** Two acinar cells whose centroids sit 18 um apart can be touching; two lymphocytes
12 um apart are not. And for matching two segmentations of *one* cell, IoU separates "the same
cell outlined twice" from "two adjacent cells" in a way distance cannot — in packed islet
tissue a neighbouring nucleus is often closer than the partner outline's centroid. That is
precisely the regime the post-Xenium re-stain pairing (S1c) works in.

The two formats disagree about almost everything, which is most of the work here:

============  ==========================================  ===========================
              PhenoCycler                                 Xenium
============  ==========================================  ===========================
file          ``redsea_scratch/geojson/cells__<tag>``     ``cell_boundaries.parquet``
layout        GeoJSON features, one per cell              long table, one row/vertex
id            feature ``id`` = QuPath object UUID         ``cell_id``
units         **pixels** (full-res qptiff)                **microns**
============  ==========================================  ===========================

The units mismatch is the trap. The GeoJSON is rasterized straight against the qptiff by
``redsea.rasterize_mask``, so its coordinates are pixels, while every contract table is in
microns. Rather than read the qptiff's resolution tags — which needs the image mounted, and
silently returns ``None`` when the tags are missing — the scale is **fitted from the data**:
the contract table and the GeoJSON describe the same cells and are joined by the same id, so
a per-axis linear fit of polygon centroid (px) against ``x_um``/``y_um`` recovers the
conversion exactly and self-validates. A fit whose residual is large means the two files are
not the same section, which is worth an error rather than a silently wrong overlay.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig
from .contract import read_cell_table
from .transform import Transform

#: A PhenoCycler pixel is ~0.3-0.5 um. Anything outside this is not a unit conversion, it is
#: two files that do not belong together — see `fit_pixel_scale`.
PLAUSIBLE_UM_PER_PX = (0.05, 2.0)

#: Residual of the px->um fit above which the join is not believable, in microns.
MAX_FIT_RESIDUAL_UM = 5.0


class BoundaryError(RuntimeError):
    """Polygons cannot be loaded, or do not correspond to the cell table."""


def _require_geo():
    try:
        import geopandas as gpd
        import shapely  # noqa: F401
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise BoundaryError(
            "polygon loading needs geopandas + shapely, which are in the integration extra:\n"
            "  pip install -e '.[integration]'\n"
            "or use environment-integration.yml."
        ) from exc
    return gpd


def read_geojson_polygons(path: Path):
    """QuPath GeoJSON -> ``GeoSeries`` indexed by object id, in the file's own units (pixels).

    Handles Polygon and MultiPolygon, and keeps interior rings: a cell with a hole punched in
    it has a smaller area than its outline, and IoU is an area ratio.
    """
    gpd = _require_geo()
    from shapely.geometry import MultiPolygon, Polygon

    feats = json.loads(Path(path).read_text()).get("features", [])
    ids, geoms = [], []
    for ft in feats:
        g = ft.get("geometry") or {}
        coords = g.get("coordinates")
        if not coords:
            continue
        rings = coords if g.get("type") == "MultiPolygon" else [coords]
        parts = []
        for poly in rings:
            shell = poly[0]
            if len(shell) < 4:
                continue
            holes = [h for h in poly[1:] if len(h) >= 4]
            parts.append(Polygon(shell, holes))
        if not parts:
            continue
        ids.append(str(ft.get("id")))
        geoms.append(parts[0] if len(parts) == 1 else MultiPolygon(parts))
    return gpd.GeoSeries(geoms, index=pd.Index(ids, name="cell_id"))


def fit_pixel_scale(polys, cells: pd.DataFrame) -> tuple[float, float, float, float]:
    """Fit ``um = scale * px + offset`` per axis from cells present in both. Returns
    ``(sx, bx, sy, by)``.

    Self-calibrating on purpose. The alternative — reading ``XResolution`` from the qptiff —
    needs the image mounted and returns ``None`` when the tags are absent or unusual, and a
    missing scale would silently become 1.0 and put the polygons a thousand microns from
    their own cells. Fitting uses the correspondence that already exists: the same cells,
    the same ids, two coordinate systems.
    """
    cent = polys.centroid
    px = pd.DataFrame({"px": cent.x.to_numpy(), "py": cent.y.to_numpy()}, index=polys.index)
    joined = px.join(cells.set_index("cell_id").loc[:, ["x_um", "y_um"]], how="inner")
    if len(joined) < 10:
        raise BoundaryError(
            f"only {len(joined)} cells are in both the polygons and the cell table — the "
            f"GeoJSON feature ids must be the QuPath object ids that `cell_id` carries")

    out = []
    for p, u in (("px", "x_um"), ("py", "y_um")):
        a, b = np.polyfit(joined[p].to_numpy(float), joined[u].to_numpy(float), 1)
        out.extend([float(a), float(b)])
    sx, bx, sy, by = out

    for s in (sx, sy):
        if not (PLAUSIBLE_UM_PER_PX[0] <= abs(s) <= PLAUSIBLE_UM_PER_PX[1]):
            raise BoundaryError(
                f"fitted scale {s:.4f} um/px is outside {PLAUSIBLE_UM_PER_PX} — the GeoJSON "
                f"and the cell table do not look like the same section")

    resid = np.hypot(joined["x_um"] - (sx * joined["px"] + bx),
                     joined["y_um"] - (sy * joined["py"] + by))
    med = float(np.median(resid))
    if med > MAX_FIT_RESIDUAL_UM:
        raise BoundaryError(
            f"px->um fit leaves a {med:.1f} um median residual (max {MAX_FIT_RESIDUAL_UM}) — "
            f"the polygons do not correspond to these cells")
    return sx, bx, sy, by


def load_pheno_polygons(cfg: PipelineConfig, donor: str, roi: str, *,
                        geojson: Optional[Path] = None, base=None):
    """PhenoCycler cell polygons for one section, converted to microns.

    ``base`` overrides the contract directory, so a post-Xenium re-stain
    (``cells_pheno_xif``) can supply its own cells for the scale fit.
    """
    from shapely.affinity import affine_transform

    table = read_cell_table(base if base is not None else cfg.cells_pheno_dir,
                            donor, roi, modality="phenocycler")
    path = Path(geojson) if geojson else _find_geojson(cfg, table.df)
    if path is None or not path.exists():
        raise BoundaryError(
            f"no GeoJSON for donor {donor}/{roi} under {cfg.geojson_dir}. REDSEA rebuilds the "
            f"filename tag from the section's `image` column; run "
            f"scripts/groovy/export_cells_geojson.groovy for it.")

    polys = read_geojson_polygons(path)
    sx, bx, sy, by = fit_pixel_scale(polys, table.df)
    # shapely affine: [a, b, d, e, xoff, yoff] applied as (a*x + b*y + xoff, d*x + e*y + yoff)
    um = polys.apply(lambda g: affine_transform(g, [sx, 0.0, 0.0, sy, bx, by]))
    return um.loc[um.index.isin(set(table.df["cell_id"]))]


def _find_geojson(cfg: PipelineConfig, cells: pd.DataFrame) -> Optional[Path]:
    """Rebuild the GeoJSON filename from the section's image, exactly as REDSEA does."""
    import re

    if "image" not in cells.columns or cells.empty:
        return None
    tag = re.sub(r"[^A-Za-z0-9._-]", "_", str(cells["image"].iloc[0]))
    path = Path(cfg.geojson_dir) / f"cells__{tag}.geojson"
    return path if path.exists() else None


#: 10x writes both boundary tables in the same long format, one row per vertex, with a
#: ``label_id`` that is the RING key -- not the cell key.  In ``cell_boundaries.parquet`` the
#: two coincide (958,744 labels over 958,744 cells on donor 6442).  In
#: ``nucleus_boundaries.parquet`` they do NOT: 958,918 labels over 907,352 cells, up to 6
#: labels on one cell, because a cell may hold more than one nucleus and a cell with none is
#: simply absent.  Splitting rings on ``label_id`` is therefore exact; splitting them by
#: looking for a vertex that repeats the group's first vertex silently FUSES the 2nd..nth
#: nucleus of every multi-nucleus cell -- measured, 47,120 cells on donor 6442.
_BOUNDARY_COLUMNS = ("cell_id", "vertex_x", "vertex_y")


def _read_boundary_rings(path: Path, *, label_ids=None):
    """Long vertex table -> (cell_id per ring, label_id per ring, ring geometries).

    Vectorised on purpose.  The previous body read the whole table with ``pd.read_parquet``
    and built one ``Polygon`` per cell in a Python ``groupby`` loop; on donor 6442 that is
    23,951,541 rows and ~3.7 minutes, and passing ``cell_ids`` did not help because the filter
    was an ``.isin`` on the long table, applied AFTER the full read.

    ``shapely.linearrings(coords, indices=)`` then ``shapely.polygons(rings)`` is the correct
    pair.  ``shapely.polygons(coords, indices=)`` raises ("geometries must be a
    one-dimensional array"), and ``shapely.polygons(rings, indices=)`` turns the second ring
    of a group into a HOLE -- an area-zero polygon rather than an error.
    """
    import shapely

    if path.suffix == ".parquet":
        import pyarrow.parquet as pq

        pf = pq.ParquetFile(path)
        available = set(pf.schema_arrow.names)
        missing = [c for c in _BOUNDARY_COLUMNS if c not in available]
        if missing:
            raise BoundaryError(f"{path} lacks {missing}")
        cols = list(_BOUNDARY_COLUMNS) + (["label_id"] if "label_id" in available else [])
        tbl = pq.read_table(path, columns=cols)          # column pushdown: 4 of N columns
        cell_id = tbl.column("cell_id").to_numpy(zero_copy_only=False).astype(str)
        vx = tbl.column("vertex_x").to_numpy(zero_copy_only=False).astype(np.float64)
        vy = tbl.column("vertex_y").to_numpy(zero_copy_only=False).astype(np.float64)
        lab = (tbl.column("label_id").to_numpy(zero_copy_only=False)
               if "label_id" in cols else None)
    else:
        df = pd.read_csv(path)
        missing = [c for c in _BOUNDARY_COLUMNS if c not in df.columns]
        if missing:
            raise BoundaryError(f"{path} lacks {missing}")
        cell_id = df["cell_id"].astype(str).to_numpy()
        vx = df["vertex_x"].to_numpy(np.float64)
        vy = df["vertex_y"].to_numpy(np.float64)
        lab = df["label_id"].to_numpy() if "label_id" in df.columns else None

    if lab is None:
        # No ring key: fall back to runs of equal cell_id, which is what the long format
        # guarantees. This cannot separate two nuclei of one cell, so it is only ever
        # reached for a boundary table that is 1:1 with cells.
        key = cell_id
    else:
        key = np.asarray(lab)

    if label_ids is not None:
        if lab is None:
            raise BoundaryError(f"{path} has no label_id column to select on")
        keep = np.isin(key, np.asarray(list(label_ids)))
        key, cell_id, vx, vy = key[keep], cell_id[keep], vx[keep], vy[keep]

    if key.size == 0:
        return np.empty(0, dtype=object), np.empty(0, dtype=np.int64), np.empty(0, dtype=object)

    # Runs, not groupby: the table is written ring-major, so equal keys are contiguous.
    # Verified rather than assumed -- a re-sorted file would otherwise fragment every ring.
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    lengths = np.diff(np.r_[starts, key.size])
    if len(np.unique(key)) != len(starts):
        raise BoundaryError(
            f"{path} is not grouped by ring: {len(starts)} runs over {len(np.unique(key))} "
            "distinct keys. Vertices of one ring must be contiguous.")

    long_enough = lengths >= 3            # a ring needs 3 distinct vertices to have an area
    ring_index = np.repeat(np.arange(starts.size), lengths)
    vertex_keep = np.repeat(long_enough, lengths)
    compact = np.cumsum(long_enough) - 1

    coords = np.column_stack([vx[vertex_keep], vy[vertex_keep]])
    rings = shapely.linearrings(coords, indices=compact[ring_index[vertex_keep]])
    geoms = shapely.polygons(rings)
    return cell_id[starts][long_enough], key[starts][long_enough], geoms


def load_xenium_polygons(source: Path, *, cell_ids=None, label_ids=None):
    """Xenium ``cell_boundaries.parquet`` -> ``GeoSeries`` indexed by ``cell_id``, in microns.

    Long format upstream — one row per vertex — so vertices are grouped per cell and closed
    into polygons. Already in microns, so no scale fit is needed or wanted.

    ``label_ids`` selects by the table's own 1-based ring key and is pushed down; ``cell_ids``
    still filters after the read, because the cell key is a string column and selecting on it
    cheaply would mean carrying a second index this lane does not have.
    """
    gpd = _require_geo()

    src = Path(source)
    path = src if src.is_file() else _find_xenium_boundaries(src)
    if path is None:
        raise BoundaryError(
            f"no cell_boundaries.parquet / .csv.gz under {src} — it ships in the Xenium "
            f"output bundle; the zarr and the phenotyped h5ad do not carry polygons")

    ids, _labels, geoms = _read_boundary_rings(path, label_ids=label_ids)
    if cell_ids is not None:
        want = set(map(str, cell_ids))
        keep = np.array([i in want for i in ids], dtype=bool)
        ids, geoms = ids[keep], geoms[keep]
    return gpd.GeoSeries(list(geoms), index=pd.Index(list(ids), name="cell_id"))


def load_xenium_nucleus_polygons(source: Path, *, label_ids=None):
    """Xenium ``nucleus_boundaries.parquet`` -> ``(GeoSeries, cell_id array)``, in microns.

    One entry per NUCLEUS, not per cell: the GeoSeries is indexed by the nucleus ``label_id``
    and the returned array says which cell each belongs to.  A cell with two nuclei appears
    twice; a cell with none does not appear at all, which is why the caller is handed the
    cell ids rather than being left to infer them.
    """
    gpd = _require_geo()

    src = Path(source)
    path = src if src.is_file() else _find_xenium_nucleus_boundaries(src)
    if path is None:
        raise BoundaryError(
            f"no nucleus_boundaries.parquet / .csv.gz under {src} — it ships in the Xenium "
            f"output bundle alongside cell_boundaries.parquet")
    ids, labels, geoms = _read_boundary_rings(path, label_ids=label_ids)
    series = gpd.GeoSeries(list(geoms), index=pd.Index(list(labels), name="label_id"))
    return series, np.asarray(ids, dtype=object)


def _find_xenium_boundaries(bundle: Path) -> Optional[Path]:
    for name in ("cell_boundaries.parquet", "cell_boundaries.csv.gz", "cell_boundaries.csv"):
        p = Path(bundle) / name
        if p.exists():
            return p
    return None


def _find_xenium_nucleus_boundaries(bundle: Path) -> Optional[Path]:
    for name in ("nucleus_boundaries.parquet", "nucleus_boundaries.csv.gz",
                 "nucleus_boundaries.csv"):
        p = Path(bundle) / name
        if p.exists():
            return p
    return None


def warp_polygons(polys, tf: Transform):
    """Push polygons through a :class:`Transform`, non-rigid field included.

    Vertex-wise rather than affine-only: ``Transform`` may carry a displacement field, and
    applying just its affine part would leave exactly the local error the B-spline stage was
    added to remove.
    """
    from shapely.geometry import MultiPolygon, Polygon
    from shapely.ops import transform as shapely_transform

    def _apply(x, y, z=None):
        xy = tf.apply(np.column_stack([np.asarray(x, float), np.asarray(y, float)]))
        return xy[:, 0], xy[:, 1]

    return polys.apply(lambda g: shapely_transform(_apply, g)
                       if isinstance(g, (Polygon, MultiPolygon)) else g)
