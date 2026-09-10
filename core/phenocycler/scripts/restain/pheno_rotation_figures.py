#!/usr/bin/env python
"""Rotation-response of the PhenoCycler->Xenium registration, per donor.

The question this answers is whether a donor's registration is merely mis-oriented. Sweep a
residual rotation on top of the stored transform and plot three curves against it:

* **islet RMSE** and **islet Dice** -- computed from islet centroids, which the *image* track
  never optimised against, so they are an out-of-sample check on it; and
* **footprint-in-tissue** -- the share of the warped Xenium footprint landing inside
  PhenoCycler tissue, which is the honest coverage measure when one section is a sub-field of
  the other.

A donor whose orientation is right shows a sharp spike in the islet curves at 0 deg. A donor
whose islet curves are flat is not mis-oriented by that amount: no rotation recovers a
correspondence that is not there, and a coverage-only gain is not evidence of one.

Run from ``core/phenocycler``::

    ../../.venv/bin/python scripts/restain/pheno_rotation_figures.py --donors 6436 6442
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
from phenocycler.integration.rasterize import (  # noqa: E402
    grid_for,
    rasterize_structures,
    tissue_mask,
)
from phenocycler.integration.register import (  # noqa: E402
    _load_structures,
    dice,
    nearest_neighbour_rmse,
)
from phenocycler.integration.transform import Transform, from_rigid  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
OUT_DIR = REPO / "data" / "development" / "restain" / "figures"
QC_UM = 10.0


def _curves(cfg, donor: str, angles: np.ndarray) -> dict:
    from scipy import ndimage

    ph = read_cell_table(CORE / "data/integration/cells_pheno", donor, "panc",
                         modality="phenocycler").xy()
    xe = read_cell_table(CORE / "data/integration/cells_xen", donor, "panc",
                         modality="xenium").xy()
    grid = grid_for(ph, QC_UM, other=xe)
    ph_mask = tissue_mask(ph, grid)
    xe_mask = tissue_mask(xe, grid)

    fx = _load_structures(cfg, "phenocycler", donor, "panc")
    mv = _load_structures(cfg, "xenium", donor, "panc")
    fx_xy = fx.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    mv_xy = mv.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    r = max(1, int(round(75.0 / grid.pixel_um)))
    mv_isl = ndimage.binary_dilation(rasterize_structures(mv, grid) > 0, iterations=r)
    fx_isl = ndimage.binary_dilation(rasterize_structures(fx, grid) > 0, iterations=r)

    base = replace(
        Transform.load(CORE / f"data/integration/registration/donor_id={donor}/roi=panc",
                       stem="transform"),
        displacement=None)
    centre = base.apply(xe).mean(axis=0)
    cxy = (float(centre[0]), float(centre[1]))

    rmse, idice, recall = [], [], []
    for deg in angles:
        tf = from_rigid(float(deg), (0.0, 0.0), center_um=cxy).compose(base)
        value, _n = nearest_neighbour_rmse(tf.apply(mv_xy), fx_xy, max_dist_um=300.0)
        rmse.append(value)
        idice.append(dice(tf.warp_mask(mv_isl, grid), fx_isl))
        warped = tf.warp_mask(xe_mask, grid)
        recall.append(float(np.logical_and(warped, ph_mask).sum() / max(warped.sum(), 1)))
    return {
        "angles": angles,
        "rmse": np.array(rmse, float),
        "islet_dice": np.array(idice, float),
        "footprint": np.array(recall, float),
        "n_islets": (len(mv_xy), len(fx_xy)),
        "rotation_deg": base.rotation_deg,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Rotation-response of the S4 registration")
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--span", type=float, default=20.0)
    ap.add_argument("--step", type=float, default=2.0)
    args = ap.parse_args(argv)

    cfg = load_config(CORE / "config.ini")
    angles = np.arange(-args.span, args.span + args.step / 2, args.step)

    data = {d: _curves(cfg, d, angles) for d in args.donors}
    n = len(args.donors)
    fig, axes = plt.subplots(3, n, figsize=(6.0 * n, 10.5), constrained_layout=True,
                             squeeze=False)

    for col, donor in enumerate(args.donors):
        c = data[donor]
        for row, (key, label, better) in enumerate((
            ("islet_dice", "islet Dice  (out-of-sample)", "higher is better"),
            ("rmse", "islet RMSE, µm  (out-of-sample)", "lower is better"),
            ("footprint", "Xenium footprint inside\nPhenoCycler tissue", "higher is better"),
        )):
            ax = axes[row][col]
            ax.plot(c["angles"], c[key], marker="o", markersize=3, linewidth=1.4,
                    color="#1f77b4")
            ax.axvline(0.0, color="#d62728", linestyle="--", linewidth=1.2)
            ax.set_xlabel("residual rotation applied on top of the stored transform (deg)",
                          fontsize=8)
            ax.set_ylabel(f"{label}\n({better})", fontsize=8)
            ax.tick_params(labelsize=8)
            ax.grid(alpha=0.25)
            if row == 0:
                ax.set_title(
                    f"donor {donor} — stored rotation {c['rotation_deg']:+.2f}°\n"
                    f"{c['n_islets'][0]} Xenium islets vs {c['n_islets'][1]} PhenoCycler islets",
                    fontsize=10)

    fig.suptitle(
        "Is the PhenoCycler→Xenium registration merely mis-oriented?  "
        "A right orientation spikes the out-of-sample islet curves at 0°.\n"
        "Red dashed line = the stored transform.",
        fontsize=11,
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"pheno_rotation_response_{'_'.join(args.donors)}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[figure] {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
