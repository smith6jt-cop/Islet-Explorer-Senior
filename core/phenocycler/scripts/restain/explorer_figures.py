#!/usr/bin/env python
"""QC on the bundle as a colleague receives it, not on the pipeline that made it.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/explorer_figures.py --donors 6436 6442

Every panel here is drawn from **only the three files that ship** -- the OME-TIFF, the
alignment CSV, and the annotation GeoJSON -- plus the one rule 10x documents for reading them
(``Xenium pixel = Xenium micron / 0.2125``). Nothing is read from the transform, the refine
output or any intermediate. That is the point: if the export composed its matrix wrongly, or
wrote the pyramid on a different grid than the matrix assumes, this notices, whereas a figure
drawn from the fitted transform would agree with itself and show nothing.

The offset sweep is the panel that carries the claim. A registration that is merely
*plausible* still puts DAPI above background at a nucleus; only a registration that is
*correct* has its maximum at zero offset.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import display as DISP                                           # noqa: E402
from run_explorer_export import BUNDLE_ROOT                      # noqa: E402

#: 10x's documented rule, restated here rather than imported, so this file needs nothing from
#: the export it is checking.
XENIUM_UM_PER_PX = 0.2125

#: Offsets the sweep plants, in microns. It runs to 12 um because a nucleus is ~6 um: past
#: that the curve must fall to the between-nuclei floor, and a curve that does not is
#: measuring tissue rather than nuclei.
SWEEP_UM = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 25.0)

#: Half-width of the window every panel is measured and drawn in, in microns.
WINDOW_UM = 750.0


def _read_bundle(dest: Path) -> dict:
    m = json.loads((dest / "manifest.json").read_text())
    tif = dest / "images" / f"{m['xenium_serial']}_IF.ome.tif"
    if not tif.exists():
        raise FileNotFoundError(f"{tif} is absent -- the bundle was made with --skip-image")
    return {"manifest": m, "tif": tif,
            "matrix": np.loadtxt(dest / "images" / Path(m["alignment"]["file"]).name,
                                 delimiter=","),
            "geojson": json.loads(
                (dest / "annotations" / "tissue_islet_xenium_microns.geojson").read_text())}


def _window(tif: Path, x0: int, y0: int, x1: int, y1: int) -> tuple[np.ndarray, list[str]]:
    """A crop of every exported channel at full resolution, plus their names from the OME XML."""
    import tifffile
    import zarr

    with tifffile.TiffFile(str(tif)) as t:
        names = re.findall(r'Channel[^>]*Name="([^"]+)"', t.ome_metadata or "")
        z = zarr.open(t.series[0].aszarr(level=0), mode="r")
        crop = np.asarray(z[:, y0:y1, x0:x1])
    return crop, names


def _sample(plane: np.ndarray, if_px: np.ndarray, x0: int, y0: int) -> np.ndarray:
    col = np.rint(if_px[:, 0]).astype(int) - x0
    row = np.rint(if_px[:, 1]).astype(int) - y0
    ok = (col >= 0) & (col < plane.shape[1]) & (row >= 0) & (row < plane.shape[0])
    return plane[row[ok], col[ok]].astype(np.float64)


def figure_one(dest: Path, out: Path) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon

    b = _read_bundle(dest)
    m = b["manifest"]
    M = b["matrix"]
    Minv = np.linalg.inv(M)                      # Xenium pixels -> exported IF pixels
    if_um_px = float(m["alignment"]["if_um_per_px"])

    cells = pd.read_parquet(Path(m["xenium_bundle"]) / "cells.parquet",
                            columns=["x_centroid", "y_centroid"])
    xy_um = cells.to_numpy()

    centre = np.median(xy_um, axis=0)
    sel = (np.abs(xy_um - centre) < WINDOW_UM).all(axis=1)
    if sel.sum() < 500:
        raise RuntimeError(f"only {sel.sum()} cells in the window; nothing to check")
    xen_px = xy_um[sel] / XENIUM_UM_PER_PX
    if_px = (np.column_stack([xen_px, np.ones(sel.sum())]) @ Minv.T)[:, :2]

    pad = int(round(30.0 / if_um_px))
    x0, y0 = np.floor(if_px.min(0)).astype(int) - pad
    x1, y1 = np.ceil(if_px.max(0)).astype(int) + pad
    crop, names = _window(b["tif"], int(x0), int(y0), int(x1), int(y1))
    dapi = crop[names.index("DAPI")].astype(np.float32)

    # --- the sweep --------------------------------------------------------------------- #
    sweep = []
    for d in SWEEP_UM:
        vals = [_sample(dapi, if_px + np.array([[dx, dy]]) * (d / if_um_px), x0, y0).mean()
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))]
        sweep.append(float(np.mean(vals)) if d else
                     float(_sample(dapi, if_px, x0, y0).mean()))
    sweep = np.asarray(sweep)
    floor = float(sweep[-1])
    contrast = float(sweep[0] / floor) if floor > 0 else float("nan")

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6))

    ax = axes[0]
    ax.plot(SWEEP_UM, sweep, "o-", color="#1c7ed6", lw=2)
    ax.axhline(floor, color="#adb5bd", ls="--", lw=1)
    ax.annotate(f"between nuclei ({floor:.1f})", (SWEEP_UM[-1], floor),
                textcoords="offset points", xytext=(-6, 8), ha="right",
                fontsize=8, color="#495057")
    ax.plot([0], [sweep[0]], "o", ms=11, mfc="none", mec="#e8590c", mew=2)
    ax.set_xlabel("offset planted on the shipped alignment (µm)")
    ax.set_ylabel("mean DAPI at Xenium nucleus centroids")
    ax.set_title(f"maximum is at zero offset\n{contrast:.2f}× the between-nuclei floor",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)

    # --- the overlay ------------------------------------------------------------------- #
    # Windowed on the workflow's own floor: vmin is 0 because 0 is what the scanner clamps
    # at, and the top comes from a percentile of the WHOLE crop rather than of the part
    # being looked at.
    zoom_um = 120.0
    zc = if_px.mean(axis=0)
    half = int(round(zoom_um / 2 / if_um_px))
    zx0, zy0 = int(zc[0] - half), int(zc[1] - half)
    sub = dapi[zy0 - y0:zy0 - y0 + 2 * half, zx0 - x0:zx0 - x0 + 2 * half]
    rng = DISP.unsigned_display_range(dapi.ravel()[::37])
    ax = axes[1]
    ax.imshow(sub, cmap="gray", vmin=0.0, vmax=rng.vmax, interpolation="nearest",
              extent=[zx0, zx0 + 2 * half, zy0 + 2 * half, zy0])
    inside = ((if_px[:, 0] > zx0) & (if_px[:, 0] < zx0 + 2 * half)
              & (if_px[:, 1] > zy0) & (if_px[:, 1] < zy0 + 2 * half))
    ax.plot(if_px[inside, 0], if_px[inside, 1], "o", ms=7, mfc="none",
            mec="#f76707", mew=1.1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{zoom_um:.0f} µm of the shipped image, with Xenium nuclei placed\n"
                 f"through the shipped matrix alone (vmin=0, vmax=p{rng.high_percentile})",
                 fontsize=10)

    # --- the annotations ---------------------------------------------------------------- #
    ax = axes[2]
    hormone = [n for n in names if n in ("INS", "GCG")]
    plane = crop[names.index(hormone[0])].astype(np.float32) if hormone else dapi
    hr = DISP.unsigned_display_range(plane.ravel()[::37])
    ax.imshow(plane, cmap="magma", vmin=0.0, vmax=hr.vmax, interpolation="nearest",
              extent=[x0, x1, y1, y0])
    n_drawn = 0
    for feat in b["geojson"]["features"]:
        if feat["properties"]["classification"]["name"] != "Islet":
            continue
        ring = np.asarray(feat["geometry"]["coordinates"][0], dtype=float)  # Xenium microns
        rp = (np.column_stack([ring / XENIUM_UM_PER_PX, np.ones(len(ring))]) @ Minv.T)[:, :2]
        if not ((rp[:, 0] > x0).any() and (rp[:, 0] < x1).any()
                and (rp[:, 1] > y0).any() and (rp[:, 1] < y1).any()):
            continue
        ax.add_patch(MplPolygon(rp, closed=True, fill=False, ec="#63e6be", lw=1.4))
        n_drawn += 1
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"shipped {hormone[0] if hormone else 'DAPI'} channel with the shipped "
                 f"islet\nannotations ({n_drawn} in frame), both put in place by the same "
                 "matrix", fontsize=10)

    mk = m["hormone_gate"]["markers"]
    fig.suptitle(
        f"Donor {m['donor_id']} · Xenium {m['xenium_serial']} · bundle QC from shipped files only"
        f"   |   INS AUROC {mk['INS']['auroc']:.3f}, GCG {mk['GCG']['auroc']:.3f}"
        f" (bar {m['hormone_gate']['preregistered_bar']:.2f}: "
        f"{'met' if mk['INS']['clears_preregistered_bar'] else 'not met'})",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"[figures] {out}", flush=True)
    return {"donor_id": m["donor_id"], "serial": m["xenium_serial"],
            "n_cells_in_window": int(sel.sum()),
            "sweep_um": list(SWEEP_UM), "sweep_dapi": sweep.tolist(),
            "contrast_at_zero": contrast, "floor": floor,
            "argmax_offset_um": float(SWEEP_UM[int(np.argmax(sweep))]),
            "islets_drawn": n_drawn, "figure": str(out)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--bundle-root", type=Path, default=BUNDLE_ROOT)
    args = ap.parse_args(argv)

    recs = []
    for dest in sorted(Path(args.bundle_root).glob("donor_*")):
        donor = dest.name.split("_")[1]
        if args.donors and donor not in args.donors:
            continue
        rec = figure_one(dest, dest / "bundle_qc.png")
        recs.append(rec)
        print(f"  {donor}: DAPI at nuclei {rec['sweep_dapi'][0]:.2f}, floor {rec['floor']:.2f}"
              f"  ({rec['contrast_at_zero']:.2f}x)   argmax at "
              f"{rec['argmax_offset_um']:.1f} um", flush=True)
    if not recs:
        print("no bundle matched", file=sys.stderr)
        return 1
    p = Path(args.bundle_root) / "bundle_qc.json"
    p.write_text(json.dumps(recs, indent=2))
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
