#!/usr/bin/env python
"""Before/after: what changing the transform actually did to the placed annotations.

A whole-section overlay understates this badly. A change worth a few hundred microns is real and
larger than an islet, but rendered across a 14 mm section it is a handful of screen pixels on
objects a hundred microns wide, which reads as "no change". So this draws both placements
together, on the same axes, and zooms.

Every caption here is computed from the two transforms it was handed. An earlier version stated
the difference in prose -- "the scale fix", "an 11% single-axis stretch" -- and that prose was
still on the figure after the transform had changed for an entirely different reason, describing
a 12.5 degree rotation as a scale correction. Nothing here is written down that the figure
cannot re-derive.

The zoom panels are chosen at increasing distance from the rotation centre, because a transform
difference that is part rotation and part scale grows outward: the centre barely moves and the
edges move most. Showing only a central crop would understate it exactly as the whole-section
view does.

Run from ``core/multiplex5000``::

    ../../.venv/bin/python scripts/restain/pheno_transform_diff_figure.py --donor 6436 \\
        --before data/integration/registration.backup_before_coarse_widen_20260910
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

import pheno_annotations as PA  # noqa: E402

from phenocycler.integration.contract import read_cell_table  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402

REPO = HERE.parents[3]
CORE = HERE.parents[1]
OUT_DIR = REPO / "data" / "development" / "restain" / "figures"
BUNDLE_ROOT = REPO / "data" / "development" / "restain" / "explorer_bundle"

BEFORE_COLOUR = "#ff8c00"
AFTER_COLOUR = "#1e90ff"


def _rings_um(donor: str, cls_name: str, root: Path, stem_root: Path) -> list[np.ndarray]:
    manifest = json.loads((CORE / "data" / "qupath_manifest.json").read_text())
    stem = PA.image_stem_for(manifest, donor)
    scale = PA.pixel_size_for(manifest, donor)
    cls = PA.CLASSES_BY_NAME[cls_name]
    polygons, _audit = PA._exterior_polygons(PA.source_path(stem_root, cls, stem), cls.name)
    import shapely

    return [np.asarray(shapely.get_coordinates(p.exterior), float) * scale for p in polygons]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Before/after placement of the same annotations")
    ap.add_argument("--donor", required=True)
    ap.add_argument("--before", type=Path, required=True,
                    help="registration root holding the previous transform")
    ap.add_argument("--class", dest="cls", default="Islet")
    ap.add_argument("--annotation-root", type=Path,
                    default=Path.home() / "IO60panc2nd" / "results" / "Geojson")
    ap.add_argument("--zooms", type=float, nargs="+", default=[1500.0, 1500.0, 1500.0])
    args = ap.parse_args(argv)

    donor = args.donor
    bundle = sorted(BUNDLE_ROOT.glob(f"donor_{donor}_*"))[0]
    serial = bundle.name.rsplit("_", 1)[-1]
    mo = json.loads((bundle / "manifest.json").read_text())["morphology"]
    upp = float(mo["um_per_px"])
    width, height = int(mo["width"]), int(mo["height"])

    after = replace(Transform.load(
        CORE / f"data/integration/registration/donor_id={donor}/roi=panc", stem="transform"),
        displacement=None)
    before = replace(Transform.load(
        args.before / f"donor_id={donor}/roi=panc", stem="transform"), displacement=None)

    rings = _rings_um(donor, args.cls, args.before, args.annotation_root)
    a_rings = [after.inverse_apply(r, n_iter=10) / upp for r in rings]
    b_rings = [before.inverse_apply(r, n_iter=10) / upp for r in rings]

    disp = np.concatenate([
        np.linalg.norm((a - b) * upp, axis=1) for a, b in zip(a_rings, b_rings)])

    xen = read_cell_table(CORE / "data/integration/cells_xen", donor, "panc",
                          modality="xenium").xy() / upp
    counts, _xe, _ye = np.histogram2d(xen[:, 0], xen[:, 1], bins=600,
                                      range=[[0, width], [0, height]])

    centres = np.array([r.mean(axis=0) for r in a_rings])
    per_ring = np.array([np.linalg.norm((a - b) * upp, axis=1).mean()
                         for a, b in zip(a_rings, b_rings)])

    # Require the outline to sit on Xenium tissue. The largest movers are at the periphery,
    # outside the cell field entirely, so ranking on displacement alone picks crops that render
    # on blank background and show the shift against nothing.
    col = np.clip((centres[:, 0] / width * counts.shape[0]).astype(int), 0, counts.shape[0] - 1)
    row = np.clip((centres[:, 1] / height * counts.shape[1]).astype(int), 0, counts.shape[1] - 1)
    on_tissue = counts[col, row] > 0
    inside = ((centres[:, 0] > 0) & (centres[:, 0] < width)
              & (centres[:, 1] > 0) & (centres[:, 1] < height) & on_tissue)
    order = np.argsort(np.where(inside, per_ring, -np.inf))[::-1]

    # Spread the crops out, so three panels are not three views of one neighbourhood.
    picks: list[int] = []
    min_gap = max(args.zooms) / upp * 1.5
    for i in order:
        if not inside[i]:
            break
        if all(np.hypot(*(centres[i] - centres[j])) > min_gap for j in picks):
            picks.append(int(i))
        if len(picks) >= len(args.zooms):
            break

    fig = plt.figure(figsize=(6.2 * (1 + len(picks)), 11), constrained_layout=True)
    gs = fig.add_gridspec(2, 1 + len(picks), height_ratios=[3, 1])

    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(counts.T, origin="upper", extent=(0, width, height, 0), cmap="Greys",
              vmin=0.0, vmax=float(counts.max()), interpolation="nearest")
    for r in b_rings:
        ax.plot(r[:, 0], r[:, 1], color=BEFORE_COLOUR, linewidth=0.6)
    for r in a_rings:
        ax.plot(r[:, 0], r[:, 1], color=AFTER_COLOUR, linewidth=0.6)
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"donor {donor} — every {args.cls} outline, both placements\n"
                 f"whole section, {len(rings):,} outlines", fontsize=10)
    ax.legend(handles=[
        plt.Line2D([], [], color=BEFORE_COLOUR, lw=2,
                   label=f"before  (rot {before.rotation_deg:+.2f}°, "
                         f"scale {before.scale:.4f})"),
        plt.Line2D([], [], color=AFTER_COLOUR, lw=2,
                   label=f"after   (rot {after.rotation_deg:+.2f}°, "
                         f"scale {after.scale:.4f})")],
        loc="lower right", fontsize=9, framealpha=0.9)

    for j, (idx, span) in enumerate(zip(picks, args.zooms)):
        axz = fig.add_subplot(gs[0, 1 + j])
        cx, cy = centres[idx]
        half = span / upp / 2
        sl = (slice(max(0, int(cy - half)), int(cy + half)),
              slice(max(0, int(cx - half)), int(cx + half)))
        sub = counts.T[
            max(0, int((cy - half) / height * 600)):int((cy + half) / height * 600),
            max(0, int((cx - half) / width * 600)):int((cx + half) / width * 600)]
        axz.imshow(sub, origin="upper", extent=(cx - half, cx + half, cy + half, cy - half),
                   cmap="Greys", vmin=0.0, vmax=float(counts.max()), interpolation="nearest")
        for r in b_rings:
            axz.plot(r[:, 0], r[:, 1], color=BEFORE_COLOUR, linewidth=1.8)
        for r in a_rings:
            axz.plot(r[:, 0], r[:, 1], color=AFTER_COLOUR, linewidth=1.8)
        axz.set_xlim(cx - half, cx + half)
        axz.set_ylim(cy + half, cy - half)
        axz.set_aspect("equal")
        axz.set_xticks([])
        axz.set_yticks([])
        axz.set_title(f"detail {j + 1} — {span:.0f} µm across\n"
                      f"this outline moved {per_ring[idx]:.0f} µm", fontsize=10)
        del sl

    axh = fig.add_subplot(gs[1, 1:])
    axh.hist(disp, bins=80, color="#444444")
    for q, style in ((50, "-"), (90, "--")):
        v = np.percentile(disp, q)
        axh.axvline(v, color="#d62728", linestyle=style, linewidth=1.5,
                    label=f"p{q} = {v:.0f} µm")
    axh.set_xlabel("how far each annotation vertex moved between the two transforms (µm)",
                   fontsize=9)
    axh.set_ylabel("vertices", fontsize=9)
    axh.legend(fontsize=9)
    axh.tick_params(labelsize=8)

    # Say what actually changed, in the numbers the transforms carry, so the caption cannot
    # describe last month's edit.
    d_rot = after.rotation_deg - before.rotation_deg
    d_scale = after.scale / before.scale - 1.0
    d_trans = np.hypot(after.affine[0, 2] - before.affine[0, 2],
                       after.affine[1, 2] - before.affine[1, 2])
    parts = []
    if abs(d_rot) >= 0.05:
        parts.append(f"{d_rot:+.2f}° of rotation")
    if abs(d_scale) >= 5e-4:
        parts.append(f"{100 * d_scale:+.2f}% of scale")
    if d_trans >= 1.0:
        parts.append(f"{d_trans:.0f} µm of translation")
    what = ", ".join(parts) if parts else "no change in the affine"

    fig.suptitle(
        f"donor {donor} / {serial} — {what} moved every {args.cls} outline by a median of "
        f"{np.median(disp):.0f} µm (p90 {np.percentile(disp, 90):.0f}, "
        f"max {disp.max():.0f} µm).\n"
        f"A difference that is part rotation grows outward from the centre of rotation, which "
        f"is why the whole-section view understates it and the details do not.",
        fontsize=12)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"pheno_transform_diff_{donor}_{args.cls}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"[figure] {out}", flush=True)
    print(f"  median {np.median(disp):.0f} um | p90 {np.percentile(disp, 90):.0f} | "
          f"max {disp.max():.0f} um", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
