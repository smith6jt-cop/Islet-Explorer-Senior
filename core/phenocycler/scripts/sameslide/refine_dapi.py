#!/usr/bin/env python3
"""
DAPI-driven refinement of the same-slide seed. This is the step point matching cannot do.

Why not centroids. Measured on 6533: nearest-neighbour spacing is ~5.4 um at ~10,000
nuclei/mm2, so a grid search over translations scored by "fraction of Xenium nuclei with a
PhenoCycler nucleus within 2 um" peaks at 0.1658 against a surface median of 0.1573 -- a
1.05x bump, z=3.0. There is no peak to optimise against, because ~15.7% of cells have SOME
neighbour within 2 um at any offset. ICP inherits that: at a 35 um residual roughly 40 cells
are closer than the true partner, so almost every mutual pair it averages is wrong, which is
how the moment-seeded fit settled 710 um from the truth while its statistics looked ordinary.

What the seed already did and did not fix. The global rigid fit is correct and is NOT the
problem: rotation -179.9998 deg, scale 1.0000019, no mirror, and the best remaining global
translation is (+2, +6) um. Yet per-window offsets are ~35 um. A residual that is ~0 globally
and +-35 um locally is non-rigid by definition, so no rigid refinement can remove it -- it
needs a field.

So refinement is driven by DAPI texture, which is what the restain lane uses to reach
sub-micron on the same problem. Windows are measured by normalised cross-correlation of the
two DAPI planes, a global affine correction is fitted to those offsets, and the structured
part that survives becomes a displacement field. The gate is the same geometric one
align.py uses: a bar derived from the donor's own median cell diameter, not a chosen number.

Coordinate care. The working frame is the qptiff's own level-0 lattice, so the FIXED image is
never resampled -- only the moving side is interpolated, once. Working at level 0 also sidesteps
the pyramid-origin correction export_dapi.origin_um exists for: that half-pixel offset is a
property of how a level was downsampled (Xenium decimates, the PhenoImager box-means), and at
factor 1 it vanishes. Whether this 70-plex Fusion pyramid box-means is an empirical question
this script does not need to answer.

Usage (after a seed exists):
    python scripts/sameslide/refine_dapi.py \
        --qptiff  .../0080870_6533_pLNpanc_Scan3_cells.er.qptiff \
        --bundle  .../output-XETG00298__0080870__panc__20260910__211917 \
        --out     .../serial=0080870 --seed-stem sameslide_v2
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "restain"))

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402

#: Verified independently for this pair: the geojson->measurement-table fit lands at
#: 0.507800 um/px on both axes with a 24-25 nm median residual, and 19200 px * 0.5078 =
#: 9750 um matches the cell extent. Xenium's 0.2125 is the vendor constant, confirmed the
#: same way (31691 * 0.2125 = 6734 um).
PHENO_UM_PER_PX = 0.5078
XENIUM_UM_PER_PX = 0.2125

WINDOW_UM = 500.0          #: match align.py's probe size
SEARCH_UM = 60.0           #: half-width of the offset search
HIGHPASS_SIGMA_UM = 15.0   #: kill illumination gradients, keep nuclear texture
MIN_PEAK_Z = 6.0           #: a window must beat noise; the flat-surface null here sits at ~3.8
MIN_WINDOWS = 30           #: below this a correction is not worth fitting
FIELD_PIXEL_UM = 250.0     #: displacement-field bin, as run_refine.LOCAL_FIELD_PIXEL_UM
FIELD_MIN_UM = 0.30        #: parsimony: no field unless this much structure survives


def load_plane(path: Path, level: int = 0, channel: int = 0) -> np.ndarray:
    """One channel of one pyramid level, as a 2-D array."""
    with tifffile.TiffFile(path) as tf:
        lv = tf.series[0].levels[level]
        arr = lv.asarray()
    if arr.ndim == 3:
        arr = arr[channel]
    return np.asarray(arr)


def highpass(a: np.ndarray, sigma_px: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter
    a = a.astype(np.float32)
    return a - gaussian_filter(a, sigma_px)


def ncc_offset(a: np.ndarray, b: np.ndarray, pix_um: float,
               search_um: float = SEARCH_UM) -> dict | None:
    """Offset that best aligns ``a`` onto ``b``, by FFT cross-correlation.

    Returns the peak's sub-window offset in microns and its z against the rest of the
    search region -- the statistic that distinguishes a real match from a flat surface.
    """
    az, bz = a - a.mean(), b - b.mean()
    if az.std() == 0 or bz.std() == 0:
        return None
    cc = np.fft.irfft2(np.fft.rfft2(az) * np.conj(np.fft.rfft2(bz)), s=a.shape)
    cc /= az.std() * bz.std() * az.size
    cc = np.fft.fftshift(cc)
    r = int(round(search_um / pix_um))
    m0, m1 = a.shape[0] // 2, a.shape[1] // 2
    sub = cc[m0 - r:m0 + r + 1, m1 - r:m1 + r + 1]
    if sub.size == 0 or sub.std() == 0:
        return None
    k = np.unravel_index(np.argmax(sub), sub.shape)
    return {"dx_um": float((k[0] - r) * pix_um), "dy_um": float((k[1] - r) * pix_um),
            "peak": float(sub[k]), "z": float((sub[k] - sub.mean()) / sub.std())}


def sample_moving(moving: np.ndarray, tf_inv: Transform, grid: GridSpec,
                  x0_um: float, y0_um: float, n: int, pix_um: float) -> np.ndarray:
    """Fill a fixed-frame tile by pulling from the moving image through ``tf_inv``."""
    from scipy.ndimage import map_coordinates
    ys, xs = np.mgrid[0:n, 0:n]
    fx = x0_um + (xs + 0.5) * pix_um
    fy = y0_um + (ys + 0.5) * pix_um
    mv = tf_inv.apply(np.column_stack([fx.ravel(), fy.ravel()]))
    r = mv[:, 1] / XENIUM_UM_PER_PX - 0.5
    c = mv[:, 0] / XENIUM_UM_PER_PX - 0.5
    out = map_coordinates(moving, np.vstack([r, c]), order=1, mode="constant", cval=0.0)
    return out.reshape(n, n).astype(np.float32)


def measure(fixed: np.ndarray, moving: np.ndarray, tf: Transform, *,
            pix_um: float, step_um: float, min_signal: float) -> list[dict]:
    """One NCC offset per occupied window over the overlap."""
    tf_inv = tf.invert()
    n = int(round(WINDOW_UM / pix_um))
    sig = HIGHPASS_SIGMA_UM / pix_um
    H, W = fixed.shape
    grid = GridSpec(0.0, 0.0, pix_um, H, W)
    out = []
    ny = int((H * pix_um) // step_um)
    nx = int((W * pix_um) // step_um)
    for iy in range(ny):
        for ix in range(nx):
            x0, y0 = ix * step_um, iy * step_um
            r0, c0 = int(y0 / pix_um), int(x0 / pix_um)
            if r0 + n > H or c0 + n > W:
                continue
            ftile = fixed[r0:r0 + n, c0:c0 + n]
            if ftile.mean() < min_signal:
                continue
            mtile = sample_moving(moving, tf_inv, grid, x0, y0, n, pix_um)
            if mtile.mean() < min_signal:
                continue
            res = ncc_offset(highpass(mtile, sig), highpass(ftile, sig), pix_um)
            if res and res["z"] >= MIN_PEAK_Z:
                res.update(x_um=x0 + WINDOW_UM / 2, y_um=y0 + WINDOW_UM / 2)
                out.append(res)
    return out


def fit_affine_correction(wins: list[dict]) -> np.ndarray:
    """Least-squares affine taking measured window centres onto their corrected positions."""
    P = np.array([[w["x_um"], w["y_um"], 1.0] for w in wins])
    D = np.array([[w["dx_um"], w["dy_um"]] for w in wins])
    sol, *_ = np.linalg.lstsq(P, D, rcond=None)      # (3,2)
    A = np.eye(3)[:2]
    A[:, :2] += sol[:2].T
    A[:, 2] += sol[2]
    return A


def build_field(wins: list[dict], residual: np.ndarray, pix_um: float,
                shape_um: tuple[float, float]) -> tuple[np.ndarray, GridSpec]:
    """Binned, gap-filled, smoothed displacement field from what the affine left behind.

    A bin with no window is UNKNOWN, not zero: leaving it at zero would ramp a real
    displacement down to nothing at the tissue edge and fabricate a gradient exactly where
    the geometry is least constrained. Empty bins are filled from the nearest populated one,
    as refine.refine_local does, then the whole field is smoothed.
    """
    from scipy.ndimage import distance_transform_edt, gaussian_filter
    from scipy.stats import binned_statistic_2d
    xs = np.array([w["x_um"] for w in wins])
    ys = np.array([w["y_um"] for w in wins])
    w_px = int(np.ceil(shape_um[0] / pix_um)) + 1
    h_px = int(np.ceil(shape_um[1] / pix_um)) + 1
    edges_x = np.arange(w_px + 1) * pix_um
    edges_y = np.arange(h_px + 1) * pix_um
    field = np.zeros((2, h_px, w_px), np.float32)
    for k in (0, 1):
        stat, *_ = binned_statistic_2d(ys, xs, residual[:, k], statistic="median",
                                       bins=[edges_y, edges_x])
        known = np.isfinite(stat)
        if not known.any():
            continue
        idx = distance_transform_edt(~known, return_distances=False, return_indices=True)
        filled = stat[tuple(idx)]
        field[k] = gaussian_filter(np.nan_to_num(filled), 1.0)
    return field, GridSpec(0.0, 0.0, pix_um, h_px, w_px)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qptiff", required=True, type=Path)
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--seed-stem", default="sameslide_v2")
    ap.add_argument("--stem", default="sameslide_dapi")
    ap.add_argument("--level", type=int, default=0, help="qptiff pyramid level to work on")
    ap.add_argument("--xenium-level", type=int, default=0)
    ap.add_argument("--step-um", type=float, default=WINDOW_UM)
    ap.add_argument("--min-signal", type=float, default=1.0)
    ap.add_argument("--iterations", type=int, default=2)
    args = ap.parse_args(argv)

    t0 = time.time()
    seed = Transform.load(args.out, stem=args.seed_stem)
    pix_um = PHENO_UM_PER_PX * (2 ** args.level)

    print(f"[refine] loading qptiff DAPI level {args.level} ...", flush=True)
    fixed = load_plane(args.qptiff, level=args.level, channel=0)
    print(f"[refine] fixed  {fixed.shape} @ {pix_um:.4f} um/px", flush=True)
    dapi = args.bundle / "morphology_focus" / "ch0000_dapi.ome.tif"
    moving = load_plane(dapi, level=args.xenium_level, channel=0)
    print(f"[refine] moving {moving.shape} @ {XENIUM_UM_PER_PX} um/px "
          f"({time.time()-t0:.0f}s)", flush=True)

    tf = seed
    history = []
    for it in range(args.iterations):
        wins = measure(fixed, moving, tf, pix_um=pix_um, step_um=args.step_um,
                       min_signal=args.min_signal)
        off = np.hypot([w["dx_um"] for w in wins], [w["dy_um"] for w in wins]) if wins else np.array([])
        rec = {"iteration": it, "n_windows": len(wins),
               "median_offset_um": float(np.median(off)) if len(off) else None,
               "median_peak_z": float(np.median([w["z"] for w in wins])) if wins else None}
        print(f"[refine] iter {it}: {rec['n_windows']} windows, "
              f"median |offset| {rec['median_offset_um']}, z {rec['median_peak_z']}", flush=True)
        if len(wins) < MIN_WINDOWS:
            rec["stopped"] = f"only {len(wins)} windows above z>={MIN_PEAK_Z}"
            history.append(rec)
            break
        A = fit_affine_correction(wins)
        corr = Transform(affine=A)
        # compose(other) is self o other -- other first. The correction is measured in the
        # FIXED frame, so it must be applied AFTER the seed: corr o tf, not tf o corr.
        # The reversed order would silently apply a fixed-frame correction in moving space.
        tf = corr.compose(tf)
        history.append(rec)

    # Whatever structure survives the affine becomes a field -- only if it is worth it.
    wins = measure(fixed, moving, tf, pix_um=pix_um, step_um=args.step_um,
                   min_signal=args.min_signal)
    field_applied = False
    if len(wins) >= MIN_WINDOWS:
        resid = np.array([[w["dx_um"], w["dy_um"]] for w in wins])
        structured = float(np.median(np.linalg.norm(resid - np.median(resid, axis=0), axis=1)))
        if structured >= FIELD_MIN_UM:
            shape_um = (fixed.shape[1] * pix_um, fixed.shape[0] * pix_um)
            f, g = build_field(wins, resid, FIELD_PIXEL_UM, shape_um)
            tf = Transform(affine=tf.affine, displacement=f, grid=g)
            field_applied = True
        history.append({"stage": "field", "structured_residual_um": structured,
                        "applied": field_applied})

    final = measure(fixed, moving, tf, pix_um=pix_um, step_um=args.step_um,
                    min_signal=args.min_signal)
    off = np.hypot([w["dx_um"] for w in final], [w["dy_um"] for w in final]) if final else np.array([])

    # Gate: the same geometric bar align.py uses -- derived from the donor's own cells, not
    # picked. self_overlap_offset_um(d, 0.75) is the displacement at which two discs of
    # diameter d still share 75% of their area.
    bar = None
    try:
        import align as align_mod
        import pandas as pd
        area = pd.read_parquet(args.bundle / "cells.parquet", columns=["cell_area"])["cell_area"]
        diam = float(2 * np.sqrt(area.median() / np.pi))
        bar = {"median_cell_diameter_um": diam,
               "pass_offset_um": float(align_mod.self_overlap_offset_um(diam, 0.75)),
               "fail_offset_um": float(align_mod.self_overlap_offset_um(diam, 0.50))}
    except Exception as exc:                                            # noqa: BLE001
        bar = {"error": f"{type(exc).__name__}: {exc}"}

    args.out.mkdir(parents=True, exist_ok=True)
    tf.save(args.out, stem=args.stem)
    report = {
        "direction": "xenium_moving -> phenocycler_fixed",
        "seed_stem": args.seed_stem,
        "working_pixel_um": pix_um,
        "window_um": WINDOW_UM, "search_um": SEARCH_UM, "min_peak_z": MIN_PEAK_Z,
        "history": history,
        "final": {"n_windows": len(final),
                  "median_offset_um": float(np.median(off)) if len(off) else None,
                  "p90_offset_um": float(np.percentile(off, 90)) if len(off) else None,
                  "median_peak_z": float(np.median([w["z"] for w in final])) if final else None},
        "field_applied": field_applied,
        "gate": bar,
        "note": ("matched_within from refine.score_transform is NOT a valid gate for this pair: "
                 "its null displaces by 25 um and at ~5.4 um nearest-neighbour spacing that "
                 "leaves the distribution unchanged, so it reads 1.00x for a good transform "
                 "and a bad one alike. Judge on window offset against the geometric bar."),
        "seconds": round(time.time() - t0, 1),
    }
    (args.out / f"{args.stem}_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"[refine] final: {report['final']}")
    print(f"[refine] gate : {bar}")
    print(f"[refine] wrote {args.out / (args.stem + '.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
