#!/usr/bin/env python
"""Where does the stored transform think the Xenium capture area sits on the PhenoCycler section?

This is the maintainer's own check, drawn from the registration instead of by hand: take the
four corners of the Xenium morphology image, push them through the stored transform, and draw
the resulting quadrilateral on the PhenoCycler DAPI. If the registration is right, that box
lands where the Xenium field actually is -- same position, same angle.

It is directly comparable to a box drawn by eye on the two images side by side, which makes it
the fastest way to confirm or refute a suspected rotation error. ``--residuals`` draws candidate
rotations alongside so the difference between them is visible rather than argued; it defaults to
the stored transform alone. It used to default to ``[0.0, -7.0]``, carrying a specific hypothesis
about donor 6436 into every figure anyone drew afterwards -- and that hypothesis was wrong. The
answer was -12.5 deg of rotation *and* 2.6 mm of translation, which no rotation-only candidate
could have expressed. A default is not the place to keep an open question.

Run from ``core/multiplex5000``::

    ../../.venv/bin/python scripts/restain/pheno_capture_box_figure.py --donors 6436 6442
"""

from __future__ import annotations

import argparse
import json
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
from phenocycler.integration.manifest import load_manifest  # noqa: E402
from phenocycler.integration.rasterize import (grid_for, normalize_image,  # noqa: E402
                                               image_to_grid, read_qptiff_channel)
from phenocycler.integration.transform import Transform, from_rigid  # noqa: E402
from phenocycler.integration.manifest import _pheno_paths  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
OUT_DIR = REPO / "data" / "development" / "restain" / "figures"
BUNDLE_ROOT = REPO / "data" / "development" / "restain" / "explorer_bundle"


def _xenium_corners_um(donor: str) -> np.ndarray:
    """The Xenium morphology image's four corners, in Xenium microns."""
    matches = sorted(BUNDLE_ROOT.glob(f"donor_{donor}_*"))
    if len(matches) != 1:
        raise SystemExit(f"donor {donor}: expected 1 Explorer bundle, found {len(matches)}")
    mo = json.loads((matches[0] / "manifest.json").read_text())["morphology"]
    w = int(mo["width"]) * float(mo["um_per_px"])
    h = int(mo["height"]) * float(mo["um_per_px"])
    return np.array([[0.0, 0.0], [w, 0.0], [w, h], [0.0, h], [0.0, 0.0]])


def _box_angle_deg(corners: np.ndarray) -> float:
    """Orientation of the box's top edge, in the PhenoCycler frame (y down, as displayed)."""
    dx = corners[1, 0] - corners[0, 0]
    dy = corners[1, 1] - corners[0, 1]
    return float(np.degrees(np.arctan2(dy, dx)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Draw the Xenium capture box on PhenoCycler DAPI")
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--residuals", type=float, nargs="+", default=[0.0],
                    help="extra rotations to draw for comparison, on top of the stored one")
    args = ap.parse_args(argv)

    cfg = load_config(CORE / "config.ini")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    colours = ["#ff2d2d", "#00e0ff", "#ffe000", "#7CFC00"]

    for donor in args.donors:
        ph = read_cell_table(CORE / "data/integration/cells_pheno", donor, "panc",
                             modality="phenocycler")
        xe = read_cell_table(CORE / "data/integration/cells_xen", donor, "panc",
                             modality="xenium")
        grid = grid_for(ph.xy(), cfg.reg_pixel_um, other=xe.xy())

        qptiff_path, _ = _pheno_paths(cfg, donor, "panc")
        img, um = read_qptiff_channel(Path(qptiff_path), channel=0, target_um=cfg.reg_pixel_um)
        fixed = normalize_image(image_to_grid(img, um, grid))

        base = replace(Transform.load(
            CORE / f"data/integration/registration/donor_id={donor}/roi=panc", stem="transform"),
            displacement=None)
        corners_um = _xenium_corners_um(donor)
        centre = base.apply(xe.xy()).mean(axis=0)
        cxy = (float(centre[0]), float(centre[1]))

        fig, ax = plt.subplots(figsize=(9, 13), constrained_layout=True)
        ax.imshow(fixed, cmap="gray", vmin=0.0, vmax=1.0, interpolation="nearest",
                  extent=(grid.origin_x_um,
                          grid.origin_x_um + grid.width * grid.pixel_um,
                          grid.origin_y_um + grid.height * grid.pixel_um,
                          grid.origin_y_um))
        ax.set_xticks([])
        ax.set_yticks([])

        lines = []
        for i, deg in enumerate(args.residuals):
            tf = base if deg == 0.0 else from_rigid(deg, (0.0, 0.0),
                                                    center_um=cxy).compose(base)
            box = tf.apply(corners_um)
            angle = _box_angle_deg(box)
            colour = colours[i % len(colours)]
            label = (f"stored transform ({tf.rotation_deg:+.2f}°)" if deg == 0.0
                     else f"stored {deg:+.1f}° residual ({tf.rotation_deg:+.2f}°)")
            ax.plot(box[:, 0], box[:, 1], color=colour, linewidth=3.0, label=label)
            lines.append(f"{label}: box top edge at {angle:+.2f}°")
            print(f"  donor {donor}: {label} -> transform rotation {tf.rotation_deg:+.2f}°, "
                  f"box top edge {angle:+.2f}°", flush=True)

        ax.legend(loc="lower left", fontsize=10, framealpha=0.85)
        ax.set_title(
            f"donor {donor} — where the registration puts the Xenium capture area\n"
            f"PhenoCycler DAPI, {grid.pixel_um:g} µm grid, vmin/vmax 0–1\n"
            + "\n".join(lines),
            fontsize=10)
        out = OUT_DIR / f"pheno_capture_box_{donor}.png"
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f"[figure] {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
