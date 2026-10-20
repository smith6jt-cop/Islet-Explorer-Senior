#!/usr/bin/env python
"""What the PhenoCycler->Xenium registration actually consumes, and how well it lands.

The annotation overlay sheets show the *output* placed on the Xenium cell field. This shows the
*input*: the two DAPI rasters `register.py` is given, on their shared micron grid, and the two
of them superimposed after the stored transform. Misalignment reads directly as colour fringing
rather than having to be inferred from a summary statistic.

Both rasters come from ``_registration_rasters``, so this is literally the array pair the
registration saw — including its own 1-99.5 percentile normalisation, which is part of the
registration contract rather than a display choice. Panels are drawn at a stated vmin/vmax of
0-1 on those normalised arrays; no extra stretch and no gamma is applied on top.

A caution the Xenium side needs: ``morphology_focus_0000.ome.tif`` declares a four-file OME
series, and this cohort's local copies carry only the first file, so tifffile zero-fills the
other three and says so on stderr. Channel 0 is the DAPI file that exists and is what gets
indexed, but the zero-filled planes are reported here so that is checkable rather than assumed.

Run from ``core/multiplex5000``::

    ../../.venv/bin/python scripts/restain/pheno_registration_inputs_figure.py --donors 6436 6442
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
from phenocycler.integration.transform import Transform  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
OUT_DIR = REPO / "data" / "development" / "restain" / "figures"


def _rgb(fixed: np.ndarray, moving: np.ndarray) -> np.ndarray:
    """PhenoCycler magenta, warped Xenium green. Agreement reads grey/white."""
    rgb = np.zeros(fixed.shape + (3,), dtype=np.float32)
    rgb[..., 0] = np.clip(fixed, 0, 1)
    rgb[..., 2] = np.clip(fixed, 0, 1)
    rgb[..., 1] = np.clip(moving, 0, 1)
    return rgb


def _panels(donor: str, cfg, n_zoom: int, zoom_um: float) -> dict:
    ph_tab = read_cell_table(CORE / "data/integration/cells_pheno", donor, "panc",
                             modality="phenocycler")
    xe_tab = read_cell_table(CORE / "data/integration/cells_xen", donor, "panc",
                             modality="xenium")
    grid = grid_for(ph_tab.xy(), cfg.reg_pixel_um, other=xe_tab.xy())

    # moving = xenium, fixed = phenocycler, matching config's fixed_modality
    moving, fixed = _registration_rasters(cfg, xe_tab, ph_tab, grid, "xenium", "phenocycler",
                                          donor)

    tf = Transform.load(CORE / f"data/integration/registration/donor_id={donor}/roi=panc",
                        stem="transform")
    warped_affine = replace(tf, displacement=None).warp_mask(moving, grid)
    warped_full = tf.warp_mask(moving, grid)

    rows, cols = np.nonzero(warped_affine > 0)
    zooms = []
    if rows.size:
        half = max(4, int(round(zoom_um / grid.pixel_um / 2)))
        rng = np.random.default_rng(0)
        picks = rng.choice(rows.size, size=min(n_zoom * 40, rows.size), replace=False)
        scored = []
        for i in picks:
            r, c = int(rows[i]), int(cols[i])
            if r - half < 0 or c - half < 0 or r + half >= grid.height or c + half >= grid.width:
                continue
            fw = fixed[r - half:r + half, c - half:c + half]
            mw = warped_affine[r - half:r + half, c - half:c + half]
            scored.append((float(fw.std() * mw.std()), r, c))     # busiest windows first
        scored.sort(reverse=True)
        seen: list[tuple[int, int]] = []
        for _s, r, c in scored:
            if all(abs(r - rr) > 4 * half or abs(c - cc) > 4 * half for rr, cc in seen):
                seen.append((r, c))
            if len(seen) >= n_zoom:
                break
        zooms = [(r, c, half) for r, c in seen]

    return {
        "grid": grid, "fixed": fixed, "moving": moving,
        "warped_affine": warped_affine, "warped_full": warped_full,
        "zooms": zooms,
        "moving_zero_fraction": float((moving == 0).mean()),
        "fixed_zero_fraction": float((fixed == 0).mean()),
        "rotation_deg": tf.rotation_deg,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Show the actual registration input rasters")
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--zooms", type=int, default=3)
    ap.add_argument("--zoom-um", type=float, default=1500.0)
    args = ap.parse_args(argv)

    cfg = load_config(CORE / "config.ini")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for donor in args.donors:
        d = _panels(donor, cfg, args.zooms, args.zoom_um)
        grid = d["grid"]
        n_zoom = len(d["zooms"])
        fig, axes = plt.subplots(2, max(3, n_zoom), figsize=(5.2 * max(3, n_zoom), 11),
                                 constrained_layout=True, squeeze=False)
        for ax in axes.ravel():
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_visible(False)

        for i, (img, title) in enumerate((
            (d["fixed"], "FIXED input: PhenoCycler DAPI\n(qptiff channel 0)"),
            (d["moving"], "MOVING input: Xenium DAPI\n(morphology_focus channel 0)"),
        )):
            ax = axes[0][i]
            ax.set_visible(True)
            ax.imshow(img, cmap="gray", vmin=0.0, vmax=1.0, interpolation="nearest")
            ax.set_title(f"{title}\non the shared {grid.pixel_um:g} µm grid, vmin/vmax 0–1",
                         fontsize=9)

        ax = axes[0][2]
        ax.set_visible(True)
        ax.imshow(_rgb(d["fixed"], d["warped_affine"]), interpolation="nearest")
        ax.set_title("after the stored transform (affine only)\n"
                     "PhenoCycler magenta · Xenium green · agreement grey", fontsize=9)

        for j, (r, c, half) in enumerate(d["zooms"]):
            ax = axes[1][j]
            ax.set_visible(True)
            sl = (slice(r - half, r + half), slice(c - half, c + half))
            ax.imshow(_rgb(d["fixed"][sl], d["warped_affine"][sl]), interpolation="nearest")
            ax.set_title(f"detail {j + 1} — {2 * half * grid.pixel_um:.0f} µm across",
                         fontsize=9)

        fig.suptitle(
            f"donor {donor} — what the registration actually consumes.  "
            f"DAPI→DAPI, stored rotation {d['rotation_deg']:+.2f}°.\n"
            f"Zero-valued fraction: PhenoCycler {d['fixed_zero_fraction']:.3f}, "
            f"Xenium {d['moving_zero_fraction']:.3f} "
            f"(a zero-filled morphology channel would read ~1.000).",
            fontsize=11,
        )
        out = OUT_DIR / f"pheno_registration_inputs_{donor}.png"
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f"[figure] donor {donor} -> {out}", flush=True)
        print(f"          fixed zero-fraction {d['fixed_zero_fraction']:.4f} | "
              f"moving zero-fraction {d['moving_zero_fraction']:.4f}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
