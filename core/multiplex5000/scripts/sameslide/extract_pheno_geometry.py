#!/usr/bin/env python3
"""
QuPath cell GeoJSON -> nucleus centroids + cell polygons, in microns.

The export is in **qptiff pixel space** while every table in this lane is in microns, and
the conversion is FITTED rather than read from a tag. That follows ``boundaries.fit_pixel_scale``
and for the same reason: the qptiff's TIFF resolution rational and its XML PixelSizeMicrons
disagree (0.2504699 vs 0.25 on the post-Xenium cohort), so a tag-derived scale can be wrong by
tens of microns across a slide. Regressing polygon centroids against the measurement table's
own ``Centroid X/Y um`` -- the same objects, both spaces -- is self-checking: the residual says
whether the conversion is right.

Nuclei are why this file matters. Registration matches Xenium nucleus centroids, and without
a PhenoCycler nucleus the cell centroid has to stand in, which leaves a systematic offset of
roughly the nucleus-to-cell-centroid distance -- comparable to the ~5.4 um nearest-neighbour
spacing at this density, so it swamps the correspondence signal it is meant to provide.

Not every cell has one: this export carries nucleusGeometry on 940,640 of 1,028,140 cells.
Cells without a nucleus keep a NaN nucleus centroid and are simply absent from the nucleus
point cloud, rather than silently falling back to their cell centroid -- mixing the two
sources is exactly the bias this file exists to remove.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import shapely

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

#: Sanity band for the fitted scale, matching boundaries.PLAUSIBLE_UM_PER_PX.
PLAUSIBLE_UM_PER_PX = (0.05, 2.0)
#: A fit worse than this means the two spaces are not related the way we assume.
MAX_FIT_RESIDUAL_UM = 5.0


def _rings(geom: dict | None) -> list | None:
    """Exterior ring of a Polygon, or of the largest part of a MultiPolygon."""
    if not geom:
        return None
    t, c = geom.get("type"), geom.get("coordinates")
    if not c:
        return None
    if t == "Polygon":
        return c[0]
    if t == "MultiPolygon":
        return max((p[0] for p in c if p), key=len, default=None)
    return None


def vector_centroids(rings: list[list | None]) -> np.ndarray:
    """True polygon centroids, built vectorised. NaN where the ring is missing/degenerate."""
    out = np.full((len(rings), 2), np.nan)
    keep, flat, idx = [], [], []
    for i, r in enumerate(rings):
        if r is None or len(r) < 4:
            continue
        keep.append(i)
        flat.append(np.asarray(r, dtype=float))
        idx.append(np.full(len(r), len(keep) - 1))
    if not keep:
        return out
    coords = np.concatenate(flat)
    # shapely.polygons(indices=...) takes LINEARRINGS, not raw coordinates -- coordinates
    # have to become rings first, which is what linearrings(indices=...) is for.
    rings_arr = shapely.linearrings(coords, indices=np.concatenate(idx))
    polys = shapely.polygons(rings_arr)
    cent = shapely.centroid(polys)
    out[np.asarray(keep)] = np.column_stack([shapely.get_x(cent), shapely.get_y(cent)])
    return out


def fit_scale(px: np.ndarray, um: np.ndarray, axis: str) -> tuple[float, float, float]:
    """Least-squares ``um = s*px + b`` on one axis, with the median absolute residual."""
    ok = np.isfinite(px) & np.isfinite(um)
    if ok.sum() < 100:
        raise SystemExit(f"{axis}: only {int(ok.sum())} joined points; cannot fit a scale")
    s, b = np.polyfit(px[ok], um[ok], 1)
    resid = float(np.median(np.abs(um[ok] - (s * px[ok] + b))))
    if not (PLAUSIBLE_UM_PER_PX[0] <= abs(s) <= PLAUSIBLE_UM_PER_PX[1]):
        raise SystemExit(f"{axis}: fitted {s:.6f} um/px is outside {PLAUSIBLE_UM_PER_PX}")
    if resid > MAX_FIT_RESIDUAL_UM:
        raise SystemExit(f"{axis}: median residual {resid:.2f} um exceeds {MAX_FIT_RESIDUAL_UM}")
    return float(s), float(b), resid


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geojson", required=True, type=Path)
    ap.add_argument("--pheno", required=True, type=Path, help="cells_pheno.parquet (micron truth)")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--no-polygons", action="store_true", help="centroids only, skip WKB")
    args = ap.parse_args(argv)

    t0 = time.time()
    doc = json.loads(args.geojson.read_bytes())
    feats = doc["features"] if isinstance(doc, dict) else doc
    print(f"[extract] parsed {len(feats):,} features in {time.time()-t0:.0f}s", flush=True)

    ids = [str(f.get("id", "")) for f in feats]
    cell_rings = [_rings(f.get("geometry")) for f in feats]
    nuc_rings = [_rings(f.get("nucleusGeometry")) for f in feats]
    n_nuc = sum(r is not None for r in nuc_rings)
    print(f"[extract] nucleusGeometry present on {n_nuc:,} / {len(feats):,} "
          f"({n_nuc/len(feats):.1%})", flush=True)

    t1 = time.time()
    cell_px = vector_centroids(cell_rings)
    nuc_px = vector_centroids(nuc_rings)
    print(f"[extract] centroids in {time.time()-t1:.0f}s", flush=True)

    truth = pd.read_parquet(args.pheno, columns=["cell_id", "x_um", "y_um"])
    geo = pd.DataFrame({"cell_id": ids,
                        "cell_px_x": cell_px[:, 0], "cell_px_y": cell_px[:, 1],
                        "nuc_px_x": nuc_px[:, 0], "nuc_px_y": nuc_px[:, 1]})
    j = geo.merge(truth, on="cell_id", how="inner")
    print(f"[extract] joined {len(j):,} / {len(geo):,} to the measurement table", flush=True)

    sx, bx, rx = fit_scale(j["cell_px_x"].to_numpy(), j["x_um"].to_numpy(), "x")
    sy, by, ry = fit_scale(j["cell_px_y"].to_numpy(), j["y_um"].to_numpy(), "y")
    print(f"[extract] scale x {sx:.6f} um/px (resid {rx:.3f} um)   "
          f"y {sy:.6f} um/px (resid {ry:.3f} um)", flush=True)

    j["nuc_x_um"] = sx * j["nuc_px_x"] + bx
    j["nuc_y_um"] = sy * j["nuc_px_y"] + by
    j["cell_x_um_geo"] = sx * j["cell_px_x"] + bx
    j["cell_y_um_geo"] = sy * j["cell_px_y"] + by

    args.out.mkdir(parents=True, exist_ok=True)
    cols = ["cell_id", "nuc_x_um", "nuc_y_um", "cell_x_um_geo", "cell_y_um_geo", "x_um", "y_um"]
    j[cols].to_parquet(args.out / "pheno_centroids.parquet", index=False)

    meta = {"n_features": len(feats), "n_nucleus": int(n_nuc), "n_joined": int(len(j)),
            "scale_um_per_px": {"x": sx, "y": sy},
            "intercept_um": {"x": bx, "y": by},
            "fit_residual_um": {"x": rx, "y": ry},
            "geojson": str(args.geojson)}

    if not args.no_polygons:
        t2 = time.time()
        keep, flat, idx = [], [], []
        for i, r in enumerate(cell_rings):
            if r is None or len(r) < 4:
                continue
            a = np.asarray(r, dtype=float)
            a[:, 0] = sx * a[:, 0] + bx
            a[:, 1] = sy * a[:, 1] + by
            keep.append(i); flat.append(a); idx.append(np.full(len(a), len(keep) - 1))
        polys = shapely.polygons(
            shapely.linearrings(np.concatenate(flat), indices=np.concatenate(idx)))
        pd.DataFrame({"cell_id": [ids[i] for i in keep],
                      "wkb": shapely.to_wkb(polys)}).to_parquet(
            args.out / "pheno_cell_polygons.parquet", index=False)
        meta["n_polygons"] = len(keep)
        print(f"[extract] {len(keep):,} cell polygons -> micron WKB in {time.time()-t2:.0f}s",
              flush=True)

    (args.out / "extract_pheno_geometry.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"[extract] wrote {args.out}  (total {time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
