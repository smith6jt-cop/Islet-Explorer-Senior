#!/usr/bin/env python
"""Which residual rotation actually aligns the DAPI? Judged on septa, at the septum scale.

The earlier rotation diagnostic scored islet centroids and a tissue mask. For donor 6436 the
islet curves were flat at every angle, which was read as "no correspondence exists". That
reading was wrong: the maintainer's own side-by-side shows the Xenium capture area sitting in
the PhenoCycler section at a clear angle, so the sections plainly do correspond and it is the
*islet structures* for that donor that are uninformative.

DAPI is the signal the registration is actually given, and interlobular septa are its sharpest
landmark -- high contrast, unambiguous, and distributed across the section. So this sweeps a
residual rotation and shows the two DAPI rasters superimposed at each one, zoomed to septa. The
correct angle is the one where the magenta and green septa become a single grey septum.

Run from ``core/phenocycler``::

    ../../.venv/bin/python scripts/restain/pheno_rotation_dapi_figure.py --donor 6436
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from phenocycler.config import load_config  # noqa: E402
from phenocycler.integration.contract import read_cell_table  # noqa: E402
from phenocycler.integration.rasterize import grid_for  # noqa: E402
from phenocycler.integration.register import _registration_rasters  # noqa: E402
from phenocycler.integration.transform import Transform, from_rigid  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
OUT_DIR = REPO / "data" / "development" / "restain" / "figures"


def _rgb(fixed, moving):
    rgb = np.zeros(fixed.shape + (3,), np.float32)
    rgb[..., 0] = np.clip(fixed, 0, 1)
    rgb[..., 2] = np.clip(fixed, 0, 1)
    rgb[..., 1] = np.clip(moving, 0, 1)
    return rgb


def _dapi_ncc(fixed, warped):
    """NCC where BOTH rasters have tissue -- the only place the comparison means anything.

    Score the *structure*, not the speckle. At 2 um a DAPI raster is mostly individual nuclei,
    and serial sections do not share nuclei, so a pixelwise NCC on raw DAPI measures
    uncorrelated texture and comes out flat and near zero whatever the angle. Blurring to the
    lobule/septum scale first is what makes the metric answer the question being asked.
    """
    both = (warped > 0.02) & (fixed > 0.02)
    if int(both.sum()) < 5000:
        return float("nan"), 0
    a = fixed[both].astype(np.float64)
    b = warped[both].astype(np.float64)
    a = a - a.mean()
    b = b - b.mean()
    d = float(np.linalg.norm(a) * np.linalg.norm(b))
    return (float((a * b).sum() / d) if d > 0 else float("nan")), int(both.sum())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Residual rotation judged on DAPI septa")
    ap.add_argument("--donor", required=True)
    ap.add_argument("--angles", type=float, nargs="+",
                    default=[0.0, -3.0, -5.0, -6.5, -8.0, -11.0])
    ap.add_argument("--zoom-um", type=float, default=2200.0)
    ap.add_argument("--details", type=int, default=3)
    ap.add_argument("--blur-um", type=float, default=60.0,
                    help="smooth to this scale before scoring, so septa and lobules dominate "
                         "the correlation instead of uncorrelated nuclear speckle")
    args = ap.parse_args(argv)

    cfg = load_config(CORE / "config.ini")
    donor = args.donor

    ph = read_cell_table(CORE / "data/integration/cells_pheno", donor, "panc",
                         modality="phenocycler")
    xe = read_cell_table(CORE / "data/integration/cells_xen", donor, "panc", modality="xenium")
    grid = grid_for(ph.xy(), cfg.reg_pixel_um, other=xe.xy())
    moving, fixed = _registration_rasters(cfg, xe, ph, grid, "xenium", "phenocycler", donor)

    base = replace(Transform.load(
        CORE / f"data/integration/registration/donor_id={donor}/roi=panc", stem="transform"),
        displacement=None)
    centre = base.apply(xe.xy()).mean(axis=0)
    cxy = (float(centre[0]), float(centre[1]))

    from scipy import ndimage

    sigma = max(1.0, args.blur_um / grid.pixel_um)
    fixed_s = ndimage.gaussian_filter(fixed, sigma)

    warps, scores = {}, {}
    for deg in args.angles:
        tf = base if deg == 0.0 else from_rigid(deg, (0.0, 0.0), center_um=cxy).compose(base)
        w = tf.warp_mask(moving, grid)
        warps[deg] = w
        scores[deg] = _dapi_ncc(fixed_s, ndimage.gaussian_filter(w, sigma))
        print(f"  residual {deg:+6.2f} deg -> DAPI NCC {scores[deg][0]:.4f} "
              f"over {scores[deg][1]:,} shared-tissue px", flush=True)

    # Detail windows FAR from the rotation centre. Rotation about a point moves nothing at that
    # point and grows with radius, so a central window cannot distinguish any of these angles:
    # at 5 mm out, 7 deg is ~610 um and unmissable; at 0.5 mm it is 61 um and invisible.
    ref = warps[0.0]
    half = max(6, int(round(args.zoom_um / grid.pixel_um / 2)))
    centre_px = grid.to_pixels(np.array([cxy]))[0]
    rows, cols = np.nonzero((ref > 0.02) & (fixed > 0.02))
    rng = np.random.default_rng(0)
    picks = rng.choice(rows.size, size=min(20000, rows.size), replace=False)
    cand = []
    for i in picks:
        r, c = int(rows[i]), int(cols[i])
        if r - half < 0 or c - half < 0 or r + half >= grid.height or c + half >= grid.width:
            continue
        radius = float(np.hypot(c - centre_px[0], r - centre_px[1])) * grid.pixel_um
        win = fixed[r - half:r + half, c - half:c + half]
        if float(win.std()) <= 0:
            continue
        cand.append((radius, r, c))
    cand.sort(reverse=True)                      # farthest from the rotation centre first
    chosen, seen = [], []
    for radius, r, c in cand:
        if all(abs(r - rr) > 2 * half or abs(c - cc) > 2 * half for rr, cc in seen):
            seen.append((r, c))
            chosen.append((r, c, radius))
        if len(chosen) >= args.details:
            break

    n = len(args.angles)
    fig, axes = plt.subplots(len(chosen), n, figsize=(4.1 * n, 4.3 * len(chosen)),
                             constrained_layout=True, squeeze=False)
    for col, deg in enumerate(args.angles):
        ncc, _npx = scores[deg]
        for row, (r, c, radius) in enumerate(chosen):
            sl = (slice(r - half, r + half), slice(c - half, c + half))
            ax = axes[row][col]
            ax.imshow(_rgb(fixed[sl], warps[deg][sl]), interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            if row == 0:
                mark = "  (stored)" if deg == 0.0 else ""
                ax.set_title(f"residual {deg:+.1f}°{mark}\nDAPI NCC {ncc:.4f}", fontsize=10)
            if col == 0:
                ax.set_ylabel(f"detail {row + 1}\n{2 * half * grid.pixel_um:.0f} µm across\n"
                              f"{radius / 1000:.1f} mm from rotation centre\n"
                              f"({np.radians(7) * radius:.0f} µm per 7°)", fontsize=8)

    best = max((d for d in args.angles if np.isfinite(scores[d][0])),
               key=lambda d: scores[d][0])
    fig.suptitle(
        f"donor {donor} — which residual rotation aligns the DAPI?  "
        f"PhenoCycler magenta · Xenium green · agreement grey.\n"
        f"The right angle merges each magenta/green septum pair into one grey septum.  "
        f"Scored on DAPI blurred to {args.blur_um:.0f} µm.  Best: {best:+.1f}° ({scores[best][0]:.4f} vs "
        f"{scores[0.0][0]:.4f} stored).",
        fontsize=12,
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"pheno_rotation_dapi_{donor}.png"
    fig.savefig(out, dpi=135)
    plt.close(fig)
    print(f"[figure] {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
