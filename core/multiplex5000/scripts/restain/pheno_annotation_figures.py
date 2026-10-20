#!/usr/bin/env python
"""Verification sheet for the PhenoCycler annotation sidecars, drawn from the shipped files.

Reads only what was published -- the sidecar GeoJSONs on the Xenium grid and the Xenium cell
table -- so the figure cannot agree with the export by construction. If the S4 transform placed
the annotations correctly, drawn islets land on Xenium tissue and the transformed Tissue outline
follows the Xenium cell field; if it did not, that is what the sheet shows.

Background is the Xenium cell field itself at a linear count scale with its own measured range
printed on the panel, not a percentile stretch.

Run from ``core/multiplex5000``::

    ../../.venv/bin/python scripts/restain/pheno_annotation_figures.py --donors 6436 6442
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

from phenocycler.integration.contract import read_cell_table  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
BUNDLE_ROOT = REPO / "data" / "development" / "restain" / "explorer_bundle"
CELLS_XEN = CORE / "data" / "integration" / "cells_xen"
OUT_DIR = REPO / "data" / "development" / "restain" / "figures"

OUTLINE_CLASSES = (("Tissue", "#1f77b4"), ("Islet", "#d62728"), ("Nerve", "#c8a800"))


def _rings(bundle: Path, serial: str, cls: str, limit: int | None = None) -> list[np.ndarray]:
    rings: list[np.ndarray] = []
    pattern = f"{serial}_pheno_{cls}_on_xenium_grid*.geojson"
    for path in sorted((bundle / "annotations").glob(pattern)):
        doc = json.loads(path.read_text())
        for feature in doc["features"]:
            rings.append(np.asarray(feature["geometry"]["coordinates"][0], dtype=float))
            if limit is not None and len(rings) >= limit:
                return rings
    return rings


def _figure(donor: str, bins: int) -> Path:
    matches = sorted(BUNDLE_ROOT.glob(f"donor_{donor}_*"))
    if len(matches) != 1:
        raise SystemExit(f"donor {donor}: expected 1 bundle, found {len(matches)}")
    bundle = matches[0]
    serial = bundle.name.rsplit("_", 1)[-1]
    qc = json.loads((bundle / "pheno_annotations_qc.json").read_text())
    upp = float(qc["morphology"]["um_per_px"])
    width = int(qc["morphology"]["width_px"])
    height = int(qc["morphology"]["height_px"])

    xen = read_cell_table(CELLS_XEN, donor, "panc", modality="xenium").xy()
    xen_px = xen / upp

    counts, xedges, yedges = np.histogram2d(
        xen_px[:, 0], xen_px[:, 1], bins=bins,
        range=[[0, width], [0, height]],
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 9), constrained_layout=True)
    for ax in axes:
        ax.set_xlim(0, width)
        ax.set_ylim(height, 0)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    vmax = float(counts.max())
    for ax in axes:
        ax.imshow(
            counts.T, origin="upper", extent=(0, width, height, 0),
            cmap="Greys", vmin=0.0, vmax=vmax, interpolation="nearest",
        )

    axes[0].set_title(
        f"donor {donor} / {serial} — Xenium cells only\n"
        f"{len(xen):,} cells, {bins}x{bins} bins, linear 0–{vmax:.0f} per bin",
        fontsize=10,
    )

    drawn: list[str] = []
    for cls, color in OUTLINE_CLASSES:
        rings = _rings(bundle, serial, cls)
        if not rings:
            continue
        for ring in rings:
            axes[1].plot(ring[:, 0], ring[:, 1], color=color, linewidth=0.5, alpha=0.9)
        drawn.append(f"{cls} n={len(rings):,}")

    reg = qc["registration"]
    axes[1].set_title(
        "+ PhenoCycler annotations through S4 (serial section)\n"
        + "  ".join(drawn)
        + f"\ntissue Dice {reg['tissue_dice']:.3f} · islet RMSE {reg['rmse_um']:.0f} µm · "
        f"displacement used: {reg['displacement_used']}",
        fontsize=10,
    )
    handles = [
        plt.Line2D([], [], color=color, linewidth=1.5, label=cls)
        for cls, color in OUTLINE_CLASSES
    ]
    axes[1].legend(handles=handles, loc="lower right", fontsize=8, framealpha=0.9)

    fig.suptitle(
        "PhenoCycler annotations placed on the Xenium morphology grid — "
        "serial section, not a cell-level correspondence",
        fontsize=11,
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"pheno_annotations_{donor}_{serial}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="QC sheet for the PhenoCycler annotation sidecars")
    ap.add_argument("--donors", nargs="+", required=True)
    ap.add_argument("--bins", type=int, default=700)
    args = ap.parse_args(argv)
    for donor in args.donors:
        print(f"[figure] donor {donor} -> {_figure(donor, args.bins)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
