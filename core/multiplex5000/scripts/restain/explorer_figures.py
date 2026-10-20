#!/usr/bin/env python
"""QC on the bundle as a colleague receives it, not on the pipeline that made it.

    cd core/multiplex5000
    ../../.venv/bin/python scripts/restain/explorer_figures.py --donors 6436 6442

Every number here is drawn from **only the files that ship** -- the native OME-TIFF, the
matrix CSV, the similarity CSV, the keypoints CSV, the warped OME-TIFF, and the manifest --
plus Xenium's own cell table and the one rule 10x documents for reading them (``Xenium pixel =
Xenium micron / 0.2125``). Nothing is read from the transform, the refine output or any
intermediate. That is the point: if the export composed its matrix wrongly, wrote the pyramid
on a different grid than the matrix assumes, or warped with the wrong half-pixel, this notices,
whereas a figure drawn from the fitted transform would agree with itself and show nothing.

Three ways into the same overlay are checked the way Explorer would read each of them:

* the native image through the shipped **matrix** (Xenium nuclei -> IF pixels via ``M^-1``);
* the native image through the shipped **similarity** (the same, via ``S^-1``);
* the **warped** image with the **identity** (pixel = micron / 0.2125, corner origin) -- exactly
  what "Yes, skip alignment" does.

The offset sweep is the panel that carries the claim. A registration that is merely
*plausible* still puts DAPI above background at a nucleus; only a registration that is
*correct* has its maximum at zero offset. The warp is additionally pinned by resampling the
native image with an implementation that shares no code with the writer
(``scipy.ndimage.map_coordinates``) and differencing it against the shipped warped pixels.
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
import explorer as EX                                            # noqa: E402
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

#: Side of the Xenium-pixel window on which the shipped warp is differenced against an
#: independent resampling of the native image.
CROSSCHECK_PX = 512


def _read_bundle(dest: Path) -> dict:
    m = json.loads((dest / "manifest.json").read_text())
    serial = m["xenium_serial"]
    tif = dest / m["alignment"].get("applies_to", f"images/{serial}_IF.ome.tif")
    if not tif.exists():
        raise FileNotFoundError(f"{tif} is absent -- the bundle was made with --skip-image")
    out = {"manifest": m, "dest": dest, "tif": tif,
           "matrix": EX.read_alignment_csv(dest / m["alignment"]["file"]),
           "similarity": None, "keypoints": None, "warped": None}
    if "similarity_diagnostic" in m:
        out["similarity"] = EX.read_alignment_csv(dest / m["similarity_diagnostic"]["file"])
    if "keypoints" in m:
        out["keypoints"] = EX.read_keypoints_csv(dest / m["keypoints"]["file"])
    if "image_on_xenium_grid" in m:
        w = dest / m["image_on_xenium_grid"]["file"]
        if not w.exists():
            raise FileNotFoundError(f"{w} is declared by the manifest and told to the maintainer "
                                    "as attempt 4, but is absent")
        out["warped"] = w
    return out


def _window(tif: Path, x0: int, y0: int, x1: int, y1: int) -> tuple[np.ndarray, list[str]]:
    """A crop of every channel at full resolution, plus their names from the OME XML.

    Clamped to the image: a request that overhangs comes back zero-padded to the requested
    shape, so callers index it with the coordinates they asked for.
    """
    import tifffile
    import zarr

    with tifffile.TiffFile(str(tif)) as t:
        names = re.findall(r'Channel[^>]*Name="([^"]+)"', t.ome_metadata or "")
        z = zarr.open(t.series[0].aszarr(level=0), mode="r")
        C, H, W = z.shape
        cx0, cy0, cx1, cy1 = max(x0, 0), max(y0, 0), min(x1, W), min(y1, H)
        crop = np.zeros((C, y1 - y0, x1 - x0), dtype=z.dtype)
        if cx1 > cx0 and cy1 > cy0:
            crop[:, cy0 - y0:cy1 - y0, cx0 - x0:cx1 - x0] = np.asarray(z[:, cy0:cy1, cx0:cx1])
    return crop, names


def _sample(plane: np.ndarray, if_px: np.ndarray, x0: int, y0: int) -> np.ndarray:
    """IF pixels are centre-indexed (the matrix carries the measured half pixel): round."""
    col = np.rint(if_px[:, 0]).astype(int) - x0
    row = np.rint(if_px[:, 1]).astype(int) - y0
    ok = (col >= 0) & (col < plane.shape[1]) & (row >= 0) & (row < plane.shape[0])
    return plane[row[ok], col[ok]].astype(np.float64)


def _sample_on_xenium_grid(plane: np.ndarray, xy_um: np.ndarray, x0: int, y0: int) -> np.ndarray:
    """Xenium pixels are corner-origin: pixel k is [k, k+1) um/0.2125, so FLOOR, never round.

    This is how Explorer reads a file added with "Yes, skip alignment": pixel = micron / 0.2125.
    """
    col = np.floor(xy_um[:, 0] / XENIUM_UM_PER_PX).astype(int) - x0
    row = np.floor(xy_um[:, 1] / XENIUM_UM_PER_PX).astype(int) - y0
    ok = (col >= 0) & (col < plane.shape[1]) & (row >= 0) & (row < plane.shape[0])
    return plane[row[ok], col[ok]].astype(np.float64)


def sweep_native(dapi: np.ndarray, if_px: np.ndarray, x0: int, y0: int, if_um_px: float,
                 offsets_um=SWEEP_UM) -> np.ndarray:
    """Mean DAPI at the mapped nuclei, then at each planted offset (mean of the 4 directions)."""
    out = []
    for d in offsets_um:
        if d == 0:
            out.append(float(_sample(dapi, if_px, x0, y0).mean()))
            continue
        vals = [_sample(dapi, if_px + np.array([[dx, dy]]) * (d / if_um_px), x0, y0).mean()
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))]
        out.append(float(np.mean(vals)))
    return np.asarray(out)


def sweep_on_xenium_grid(plane: np.ndarray, xy_um: np.ndarray, x0: int, y0: int,
                         offsets_um=SWEEP_UM) -> np.ndarray:
    """The same sweep on an image that is already on the Xenium grid (identity placement)."""
    out = []
    for d in offsets_um:
        if d == 0:
            out.append(float(_sample_on_xenium_grid(plane, xy_um, x0, y0).mean()))
            continue
        vals = [_sample_on_xenium_grid(plane, xy_um + np.array([[dx, dy]]) * d, x0, y0).mean()
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))]
        out.append(float(np.mean(vals)))
    return np.asarray(out)


def _sig_digits(token: str) -> int:
    return len(token.replace("-", "").replace(".", "").lstrip("0").rstrip("0") or "0")


def csv_file_facts(path: Path) -> dict:
    """What the bytes of an alignment CSV are, checked against Explorer's own export layout."""
    b = Path(path).read_bytes()
    rows = b.split(b"\r\n")
    toks = b.decode("ascii", errors="replace").replace("\r\n", ",").split(",")
    return {"bytes": len(b), "crlf_count": b.count(b"\r\n"), "bare_lf": b.count(b"\n") - b.count(b"\r\n"),
            "trailing_newline": b.endswith((b"\n", b"\r")), "last_row": rows[-1].decode("ascii", "replace"),
            "max_significant_digits": max(_sig_digits(t) for t in toks), "ascii": all(c < 128 for c in b),
            "layout_ok": (b.count(b"\r\n") == 2 and b.count(b"\n") == 2 and not b.endswith((b"\n", b"\r"))
                          and rows[-1] == b"0,0,1" and all(c < 128 for c in b)
                          and max(_sig_digits(t) for t in toks) <= EX.ALIGNMENT_SIGNIFICANT_DIGITS)}


def file_facts(b: dict) -> dict:
    """Everything that can be checked without reading a pixel."""
    import tifffile

    m, dest = b["manifest"], b["dest"]
    a = m["alignment"]
    facts = {"matrix_csv": csv_file_facts(dest / a["file"])}
    facts["matrix_equals_manifest_as_written"] = bool(
        "matrix_as_written" in a and np.array_equal(b["matrix"], np.asarray(a["matrix_as_written"])))
    facts["matrix_close_to_manifest_matrix"] = bool(
        np.allclose(b["matrix"], np.asarray(a["matrix"]), rtol=1e-8, atol=0.0))
    if b["similarity"] is not None:
        s = m["similarity_diagnostic"]
        facts["similarity_csv"] = csv_file_facts(dest / s["file"])
        facts["similarity_is_exact_similarity"] = EX.is_exact_similarity(b["similarity"])
        facts["similarity_equals_manifest_as_written"] = bool(
            "matrix_as_written" in s and np.array_equal(b["similarity"], np.asarray(s["matrix_as_written"])))
    if b["keypoints"] is not None:
        kp = b["keypoints"]
        facts["keypoints_n"] = int(len(kp))
        facts["keypoints_residual_px"] = EX.keypoints_residual_px(kp, b["matrix"])
        if b["similarity"] is not None:
            refit = EX.similarity_from_keypoints(kp)
            facts["keypoints_refit_matches_shipped_similarity"] = bool(
                np.allclose(EX.rounded_alignment_matrix(refit), b["similarity"], rtol=0, atol=1e-9))
    if b["warped"] is not None:
        with tifffile.TiffFile(str(b["warped"])) as t:
            shapes = [[int(lv.shape[-2]), int(lv.shape[-1])] for lv in t.series[0].levels]
            pg = t.pages[0]
            ome = t.ome_metadata or ""
            facts["warped"] = {
                "is_ome": bool(t.is_ome), "shape": list(t.series[0].shape), "n_levels": len(shapes),
                "levels_equal_morphology": shapes == [list(s) for s in m["morphology"]["level_shapes"]],
                "physical_size_x": (re.search(r'PhysicalSizeX="([^"]+)"', ome) or [None, None])[1],
                "resolution_unit": int(pg.tags["ResolutionUnit"].value) if "ResolutionUnit" in pg.tags else None,
                "tile": [int(pg.tilelength), int(pg.tilewidth)], "compression": int(pg.compression),
                "channels": re.findall(r'Channel[^>]*Name="([^"]+)"', ome)}
        with tifffile.TiffFile(str(b["tif"])) as t:
            facts["warped"]["channels_equal_native"] = (
                facts["warped"]["channels"] == re.findall(r'Channel[^>]*Name="([^"]+)"', t.ome_metadata or ""))
    return facts


def warp_crosscheck(native_tif: Path, warped_tif: Path, M: np.ndarray, wx0: int, wy0: int,
                    size: int = CROSSCHECK_PX) -> dict:
    """The shipped warped DAPI against the native DAPI resampled independently at M^-1 (k+.5, l+.5).

    ``scipy.ndimage.map_coordinates`` shares nothing with ``cv2.warpAffine``; a half-pixel
    convention error would show as a systematic difference of tens of counts at nucleus edges,
    whereas the writer's 1/32-px fixed-point coordinates leave a p99 of about one count.
    """
    from scipy.ndimage import map_coordinates

    wcrop, wnames = _window(warped_tif, wx0, wy0, wx0 + size, wy0 + size)
    warped = wcrop[wnames.index("DAPI")].astype(np.float64)
    ks, ls = np.meshgrid(np.arange(wx0, wx0 + size), np.arange(wy0, wy0 + size))
    centres = np.column_stack([ks.ravel() + 0.5, ls.ravel() + 0.5, np.ones(ks.size)])
    src = (centres @ np.linalg.inv(M).T)[:, :2]
    x0, y0 = int(np.floor(src[:, 0].min())) - 3, int(np.floor(src[:, 1].min())) - 3
    x1, y1 = int(np.ceil(src[:, 0].max())) + 4, int(np.ceil(src[:, 1].max())) + 4
    ncrop, nnames = _window(native_tif, x0, y0, x1, y1)
    native = ncrop[nnames.index("DAPI")].astype(np.float64)
    resampled = map_coordinates(native, [src[:, 1] - y0, src[:, 0] - x0], order=1,
                                mode="constant", cval=0.0).reshape(size, size)
    # `_window` zero-pads a crop that overhangs the image, so scipy and cv2 both blend the
    # edge pixel toward zero there and would agree anyway; the border band is still excluded
    # so the reported numbers describe the warp of real pixels only, and the count is recorded.
    import tifffile
    with tifffile.TiffFile(str(native_tif)) as t:
        h_if, w_if = t.series[0].levels[0].shape[-2:]
    supported = ((src[:, 0] >= 0) & (src[:, 0] <= w_if - 1)
                 & (src[:, 1] >= 0) & (src[:, 1] <= h_if - 1)).reshape(size, size)
    diff = warped - np.rint(resampled)
    diff[~supported] = np.nan
    ad = np.abs(diff[supported])
    return {"window_xenium_px": [wx0, wy0, size], "n_compared": int(supported.sum()),
            "n_excluded_at_if_border": int((~supported).sum()),
            "max_abs_diff": float(ad.max()), "p99_abs_diff": float(np.percentile(ad, 99)),
            "mean_diff": float(diff[supported].mean()), "frac_within_1": float((ad <= 1).mean()),
            "diff": diff}


def figure_one(dest: Path, out: Path) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Polygon as MplPolygon

    b = _read_bundle(dest)
    m = b["manifest"]
    M, S = b["matrix"], b["similarity"]
    Minv = np.linalg.inv(M)                      # Xenium pixels -> exported IF pixels
    if_um_px = float(m["alignment"]["if_um_per_px"])

    cells = pd.read_parquet(Path(m["xenium_bundle"]) / "cells.parquet",
                            columns=["x_centroid", "y_centroid"])
    xy_um = cells.to_numpy()
    centre = np.median(xy_um, axis=0)
    sel = (np.abs(xy_um - centre) < WINDOW_UM).all(axis=1)
    if sel.sum() < 500:
        raise RuntimeError(f"only {sel.sum()} cells in the window; nothing to check")
    xy_sel = xy_um[sel]
    xen_px = xy_sel / XENIUM_UM_PER_PX
    hom = np.column_stack([xen_px, np.ones(len(xen_px))])
    if_px = (hom @ Minv.T)[:, :2]

    pad = int(round(30.0 / if_um_px))
    x0, y0 = np.floor(if_px.min(0)).astype(int) - pad
    x1, y1 = np.ceil(if_px.max(0)).astype(int) + pad
    crop, names = _window(b["tif"], int(x0), int(y0), int(x1), int(y1))
    dapi = crop[names.index("DAPI")].astype(np.float32)

    rec = {"donor_id": m["donor_id"], "serial": m["xenium_serial"],
           "n_cells_in_window": int(sel.sum()), "sweep_um": list(SWEEP_UM)}

    # --- the sweeps ------------------------------------------------------------------- #
    curves = {}
    sw = sweep_native(dapi, if_px, x0, y0, if_um_px)
    curves["native_through_matrix"] = sw
    if S is not None:
        if_px_S = (hom @ np.linalg.inv(S).T)[:, :2]
        curves["native_through_similarity"] = sweep_native(dapi, if_px_S, x0, y0, if_um_px)
    wdapi = None
    if b["warped"] is not None:
        wpad = int(round(30.0 / XENIUM_UM_PER_PX))
        wx0, wy0 = np.floor(xen_px.min(0)).astype(int) - wpad
        wx1, wy1 = np.ceil(xen_px.max(0)).astype(int) + wpad
        wcrop, wnames = _window(b["warped"], int(wx0), int(wy0), int(wx1), int(wy1))
        wdapi = wcrop[wnames.index("DAPI")].astype(np.float32)
        curves["warped_through_identity"] = sweep_on_xenium_grid(wdapi, xy_sel, int(wx0), int(wy0))
    for key, c in curves.items():
        floor = float(c[-1])
        rec[key] = {"sweep_dapi": c.tolist(), "floor": floor,
                    "contrast_at_zero": float(c[0] / floor) if floor > 0 else None,
                    "argmax_offset_um": float(SWEEP_UM[int(np.argmax(c))])}
    if "warped_through_identity" in rec:
        rec["warped_over_native_peak_ratio"] = float(
            curves["warped_through_identity"][0] / curves["native_through_matrix"][0])

    # --- the warp cross-check --------------------------------------------------------- #
    diff_img = None
    if b["warped"] is not None and int(m["alignment"].get("if_pyramid_level", 0)) == 0:
        cx, cy = np.floor(centre / XENIUM_UM_PER_PX).astype(int)
        cc = warp_crosscheck(b["tif"], b["warped"], M, int(cx - CROSSCHECK_PX // 2),
                             int(cy - CROSSCHECK_PX // 2))
        diff_img = cc.pop("diff")
        rec["warp_crosscheck"] = cc
    elif b["warped"] is not None:
        rec["warp_crosscheck"] = {"skipped": "the native image ships at pyramid level "
                                  f"{m['alignment']['if_pyramid_level']} but the warp is from "
                                  "level 0; a difference here would measure the pyramid, not "
                                  "the warp"}

    # --- file facts ------------------------------------------------------------------- #
    facts = file_facts(b)
    rec["file_facts"] = facts

    checks = {k: v["argmax_offset_um"] == 0.0 for k, v in rec.items()
              if isinstance(v, dict) and "argmax_offset_um" in v}
    checks["matrix_csv_layout"] = facts["matrix_csv"]["layout_ok"]
    checks["matrix_equals_manifest_as_written"] = facts["matrix_equals_manifest_as_written"]
    checks["matrix_close_to_manifest_matrix"] = facts["matrix_close_to_manifest_matrix"]
    if "similarity_csv" in facts:
        checks["similarity_csv_layout"] = facts["similarity_csv"]["layout_ok"]
        checks["similarity_is_exact_similarity"] = facts["similarity_is_exact_similarity"]
    if "keypoints_residual_px" in facts:
        checks["keypoints_residual_under_1e-6_px"] = facts["keypoints_residual_px"] < 1e-6
    if "keypoints_refit_matches_shipped_similarity" in facts:
        checks["keypoints_refit_matches_shipped_similarity"] = facts["keypoints_refit_matches_shipped_similarity"]
    if "warped" in facts:
        checks["warped_levels_equal_morphology"] = facts["warped"]["levels_equal_morphology"]
        checks["warped_channels_equal_native"] = facts["warped"]["channels_equal_native"]
        checks["warped_physical_size_0.2125"] = facts["warped"]["physical_size_x"] == "0.2125"
    if "warp_crosscheck" in rec and "p99_abs_diff" in rec["warp_crosscheck"]:
        checks["warp_crosscheck_p99_under_1_count"] = rec["warp_crosscheck"]["p99_abs_diff"] <= 1.0
    rec["checks"] = checks
    rec["pass"] = all(checks.values())

    # --- the figure --------------------------------------------------------------------- #
    fig, axes = plt.subplots(2, 3, figsize=(17.5, 11))
    ax = axes[0, 0]
    colours = {"native_through_matrix": "#1c7ed6", "native_through_similarity": "#f76707",
               "warped_through_identity": "#2b8a3e"}
    for key, c in curves.items():
        ax.plot(SWEEP_UM, c, "o-", color=colours[key], lw=2, ms=4,
                label=f"{key.replace('_', ' ')} (argmax {rec[key]['argmax_offset_um']:.1f} µm, "
                      f"{rec[key]['contrast_at_zero']:.2f}× floor)")
    ax.set_xlabel("offset planted on the shipped placement (µm)")
    ax.set_ylabel("mean DAPI at Xenium nucleus centroids")
    ax.set_title("every route must peak at zero offset", fontsize=10)
    ax.legend(fontsize=7, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)

    # Windowed on the workflow's own floor: vmin is 0 because 0 is what the scanner clamps
    # at, and the top comes from a percentile of the WHOLE crop rather than of the part
    # being looked at.
    zoom_um = 120.0
    zc = if_px.mean(axis=0)
    half = int(round(zoom_um / 2 / if_um_px))
    zx0, zy0 = int(zc[0] - half), int(zc[1] - half)
    sub = dapi[zy0 - y0:zy0 - y0 + 2 * half, zx0 - x0:zx0 - x0 + 2 * half]
    rng = DISP.unsigned_display_range(dapi.ravel()[::37])
    ax = axes[0, 1]
    # IF pixels are centre-indexed (column c sits AT coordinate c), and matplotlib's extent
    # gives pixel EDGES, so the edges are at -0.5; the warped panel below is corner-origin and
    # needs no such shift.
    ax.imshow(sub, cmap="gray", vmin=0.0, vmax=rng.vmax, interpolation="nearest",
              extent=[zx0 - 0.5, zx0 + 2 * half - 0.5, zy0 + 2 * half - 0.5, zy0 - 0.5])
    inside = ((if_px[:, 0] > zx0) & (if_px[:, 0] < zx0 + 2 * half)
              & (if_px[:, 1] > zy0) & (if_px[:, 1] < zy0 + 2 * half))
    ax.plot(if_px[inside, 0], if_px[inside, 1], "o", ms=7, mfc="none", mec="#f76707", mew=1.1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{zoom_um:.0f} µm of the native image, Xenium nuclei placed through the\n"
                 f"shipped matrix alone (vmin=0, vmax=p{rng.high_percentile})", fontsize=10)

    ax = axes[0, 2]
    if wdapi is not None:
        whalf = int(round(zoom_um / 2 / XENIUM_UM_PER_PX))
        wzc = np.floor(centre / XENIUM_UM_PER_PX).astype(int)
        wzx0, wzy0 = int(wzc[0] - whalf), int(wzc[1] - whalf)
        wsub = wdapi[wzy0 - wy0:wzy0 - wy0 + 2 * whalf, wzx0 - wx0:wzx0 - wx0 + 2 * whalf]
        wrng = DISP.unsigned_display_range(wdapi.ravel()[::37])
        ax.imshow(wsub, cmap="gray", vmin=0.0, vmax=wrng.vmax, interpolation="nearest",
                  extent=[wzx0, wzx0 + 2 * whalf, wzy0 + 2 * whalf, wzy0])
        winside = ((xen_px[:, 0] > wzx0) & (xen_px[:, 0] < wzx0 + 2 * whalf)
                   & (xen_px[:, 1] > wzy0) & (xen_px[:, 1] < wzy0 + 2 * whalf))
        ax.plot(xen_px[winside, 0], xen_px[winside, 1], "o", ms=7, mfc="none", mec="#2b8a3e",
                mew=1.1)
        ax.set_title(f"{zoom_um:.0f} µm of the WARPED image with nuclei at µm / 0.2125 —\n"
                     "identity placement, how Explorer reads it under 'skip alignment'",
                     fontsize=10)
    else:
        ax.text(0.5, 0.5, "no warped image in this bundle", ha="center", va="center")
    ax.set_xticks([]); ax.set_yticks([])

    # --- footprint --------------------------------------------------------------------- #
    ax = axes[1, 0]
    Wm, Hm = m["morphology"]["width"], m["morphology"]["height"]
    ax.add_patch(Rectangle((0, 0), Wm, Hm, fill=False, ec="#1c7ed6", lw=2, label="morphology FOV"))
    fp = np.asarray(m["alignment"]["footprint_xenium_px"]["if_corners_xenium_px"])
    order = [0, 1, 3, 2]
    ax.add_patch(MplPolygon(fp[order], closed=True, fill=False, ec="#f76707", lw=1.5,
                            label="native IF corners through the matrix"))
    if b["keypoints"] is not None:
        kp = b["keypoints"]
        ax.plot(kp[:, 0], kp[:, 1], "x", color="#2b8a3e", ms=6, label=f"{len(kp)} keypoints (fixed)")
    ax.set_xlim(min(fp[:, 0].min(), 0) - 2000, max(fp[:, 0].max(), Wm) + 2000)
    ax.set_ylim(max(fp[:, 1].max(), Hm) + 2000, min(fp[:, 1].min(), 0) - 2000)
    ax.set_aspect("equal")
    ax.set_xlabel("Xenium px"); ax.set_ylabel("Xenium px")
    ax.legend(fontsize=7, loc="lower right")
    ax.set_title("where the native image lands on the Xenium grid\n"
                 "(the picture of 'unaligned' under skip-alignment)", fontsize=10)

    # --- the cross-check --------------------------------------------------------------- #
    ax = axes[1, 1]
    if diff_img is not None:
        cc = rec["warp_crosscheck"]
        ax.hist(diff_img[np.isfinite(diff_img)].ravel(), bins=np.arange(-6.5, 7.5, 1.0),
                color="#2b8a3e")
        ax.set_yscale("log")
        ax.set_xlabel("shipped warped DAPI − native resampled with scipy at M⁻¹(k+½, l+½)")
        ax.set_title(f"{CROSSCHECK_PX}² Xenium px ({cc['n_excluded_at_if_border']} at the IF border "
                     f"excluded): p99 |diff| {cc['p99_abs_diff']:.1f}, max {cc['max_abs_diff']:.0f}, "
                     f"{100 * cc['frac_within_1']:.1f}% within 1 count", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    else:
        ax.text(0.5, 0.5, "no warped image in this bundle", ha="center", va="center")
        ax.set_xticks([]); ax.set_yticks([])

    # --- the facts --------------------------------------------------------------------- #
    ax = axes[1, 2]
    ax.axis("off")
    lines = [f"{'PASS' if rec['pass'] else 'FAIL'} — {sum(checks.values())}/{len(checks)} checks"]
    lines += [f"{'✓' if v else '✗'} {k}" for k, v in checks.items()]
    mc = facts["matrix_csv"]
    lines += ["", f"matrix csv: {mc['bytes']} B, {mc['crlf_count']} CRLF, trailing newline "
                  f"{mc['trailing_newline']}, last row '{mc['last_row']}', ≤{mc['max_significant_digits']} digits"]
    if "keypoints_residual_px" in facts:
        lines.append(f"keypoints residual through the shipped matrix: {facts['keypoints_residual_px']:.1e} px")
    if "warped" in facts:
        w = facts["warped"]
        lines.append(f"warped: {w['shape']}, {w['n_levels']} levels, PhysicalSizeX {w['physical_size_x']}, "
                     f"tile {w['tile']}, channels {w['channels']}")
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=8, family="monospace")

    fig.suptitle(f"Donor {m['donor_id']} · Xenium {m['xenium_serial']} · bundle QC from shipped "
                 f"files only   |   positivity ≥ {m['positivity']['threshold_counts']:.0f}/255: "
                 + ", ".join(f"{k} {v:,}" for k, v in sorted(m["positivity"]["counts"].items())),
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[figures] {out}", flush=True)
    rec["figure"] = str(out)
    return rec


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
        parts = [f"{k.replace('_', ' ')}: {v['sweep_dapi'][0]:.1f} vs floor {v['floor']:.1f} "
                 f"({v['contrast_at_zero']:.2f}x, argmax {v['argmax_offset_um']:.1f} um)"
                 for k, v in rec.items() if isinstance(v, dict) and "argmax_offset_um" in v]
        print(f"  {donor}: {'PASS' if rec['pass'] else 'FAIL'} | " + " | ".join(parts), flush=True)
        cc = rec.get("warp_crosscheck", {})
        if "p99_abs_diff" in cc:
            print(f"  {donor}: warp cross-check p99 |diff| {cc['p99_abs_diff']:.2f}, max "
                  f"{cc['max_abs_diff']:.0f} counts, {cc['n_excluded_at_if_border']} border px "
                  "excluded", flush=True)
        elif cc:
            print(f"  {donor}: warp cross-check skipped: {cc['skipped']}", flush=True)
        for k, v in rec["checks"].items():
            if not v:
                print(f"  {donor}: FAILED {k}", flush=True)
    if not recs:
        print("no bundle matched", file=sys.stderr)
        return 1
    p = Path(args.bundle_root) / "bundle_qc.json"
    p.write_text(json.dumps(recs, indent=2, allow_nan=False))
    print(f"wrote {p}")
    return 0 if all(r["pass"] for r in recs) else 2


if __name__ == "__main__":
    raise SystemExit(main())
