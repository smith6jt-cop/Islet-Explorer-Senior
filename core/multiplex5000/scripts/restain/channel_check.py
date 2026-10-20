#!/usr/bin/env python
"""Which qptiff channel actually holds insulin, shown rather than asserted.

    cd core/multiplex5000
    ../../.venv/bin/python scripts/restain/channel_check.py --donors 6436 6442

**Why this exists.** ``restain_sections.csv`` records ``Opal 690 = INS, Opal 520 = GCG`` for
every section, hand-recorded because the fact is in no file. The export lane believed it and
reported that donor 6436's INS channel does not distinguish islet from exocrine. Measured in
the IF's own frame, with no registration anywhere in the path:

    6436   Opal 520  36.58 inside islets / 1.21 outside = 30.2x     Opal 690  2.78 / 1.64 = 1.7x
    6442   Opal 690  39.04 / 2.87 = 13.6x                           Opal 520  23.18 / 0.40 = 57.8x

So 6436 has an excellent islet channel and it is not the one recorded as INS. That is either a
recording error on one section or a stain that landed in a different Opal, and **this module
does not decide which** -- the mapping is the maintainer's recorded fact and is resolved by
looking, not by a ratio. This draws what there is to look at.

Nothing here reads the transform, the AUROC, or the export. Every panel is the raw qptiff with
the maintainer's own QuPath polygons on top, in the qptiff's own frame.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))

import geo                                                       # noqa: E402
import hormone as HZ                                             # noqa: E402
from run_hormone import IMAGES, SECTIONS                         # noqa: E402
from run_register import GEOJSON_DIR, OUT, XENIUM_PATHS, discover, donor_of  # noqa: E402

#: Level to draw whole-section thumbnails from. Level 5 is ~8 um/px -- an islet is ~10 px, which
#: is enough to see WHERE the bright things are, which is the question this row answers.
SECTION_LEVEL = 5
#: Level to draw the close-ups from. Level 0 is 0.25 um/px: individual cells inside an islet.
ISLET_LEVEL = 0
#: Half-width of an islet close-up, in microns. An islet is 50-250 um across, so this frames one
#: with enough exocrine around it to judge the contrast honestly.
ISLET_HALF_UM = 110.0
#: How many islets to draw. Taken across the size range rather than the top N, so the sheet is
#: not three pictures of the same exceptional object.
N_ISLETS = 3


def _levels(qptiff: Path):
    import tifffile

    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        return [(int(lv.shape[-2]), int(lv.shape[-1]))
                for lv in getattr(series, "levels", [series])]


def _plane(qptiff: Path, channel: int, level: int, box=None) -> np.ndarray:
    """One channel at one level, optionally a ``(x0, y0, x1, y1)`` pixel box of that level."""
    import tifffile
    import zarr

    with tifffile.TiffFile(str(qptiff)) as tif:
        store = tif.series[0].aszarr(level=level)
        try:
            z = zarr.open(store, mode="r")
            if box is None:
                return np.asarray(z[channel]).astype(np.float32)
            x0, y0, x1, y1 = box
            return np.asarray(z[channel, y0:y1, x0:x1]).astype(np.float32)
        finally:
            store.close()


def _vmax(plane: np.ndarray, mask: np.ndarray | None, pct: float = 99.5) -> float:
    """Top of the display range, from the whole tissue rather than from what is drawn.

    Selection-independent on purpose: a vmax taken from the crop would rescale every panel to
    its own contents and make an empty channel look exactly like a full one, which is the
    error this whole sheet exists to expose. vmin is 0 everywhere and is never a percentile.
    """
    vals = plane[mask] if mask is not None else plane.ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 1.0
    return max(float(np.percentile(vals, pct)), 1.0)


def check_one(serial: str, gj: Path, sections, out_dir: Path) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shapely
    from matplotlib.patches import Polygon as MplPolygon

    from phenocycler.integration.rasterize import GridSpec

    donor = donor_of(serial, XENIUM_PATHS)
    row = sections[sections.xenium_serial == serial].iloc[0]
    qptiff = IMAGES / row.image_file
    names = HZ.channel_names(qptiff)
    roles = {}
    for col, val in ((k, str(v).strip()) for k, v in row.items()
                     if isinstance(k, str) and k.startswith("channel_")):
        page = HZ.channel_index_for_role(qptiff, row, val)
        if page is not None:
            roles[page] = val

    ann = geo.load_annotations(gj)
    islets, tissue = geo.islets(ann), geo.tissue_polygon(ann)
    shapes = _levels(qptiff)

    # --- the section row, and the enrichment each channel actually shows ----------------- #
    sh, sw = shapes[SECTION_LEVEL]
    s_um = geo.PIXEL_UM * (shapes[0][0] // sh)
    s_grid = GridSpec(origin_x_um=0.5 * s_um, origin_y_um=0.5 * s_um, pixel_um=s_um,
                      height=sh, width=sw)
    m_isl = geo.rasterize(islets, s_grid)
    m_tis = geo.rasterize([tissue], s_grid)
    m_out = m_tis & ~m_isl

    planes, stats = [], []
    for i, nm in enumerate(names):
        a = _plane(qptiff, i, SECTION_LEVEL)
        inside = float(a[m_isl].mean()) if m_isl.any() else float("nan")
        outside = float(a[m_out].mean()) if m_out.any() else float("nan")
        planes.append(a)
        stats.append({"channel": nm, "recorded_role": roles.get(i, ""),
                      "islet_mean": inside, "elsewhere_mean": outside,
                      "ratio": inside / outside if outside else float("nan"),
                      "vmax": _vmax(a, m_tis)})

    # --- the islets: across the size range, not the top three ---------------------------- #
    areas = np.array([g.area for g in islets])
    order = np.argsort(areas)[::-1]
    picks = [islets[order[k]] for k in
             np.unique(np.linspace(0, min(len(order), 40) - 1, N_ISLETS).astype(int))]

    n_ch = len(names)
    fig, axes = plt.subplots(1 + len(picks), n_ch, figsize=(3.1 * n_ch, 3.25 * (1 + len(picks))))
    axes = np.atleast_2d(axes)

    for j, (a, st) in enumerate(zip(planes, stats)):
        ax = axes[0, j]
        ax.imshow(a, cmap="gray", vmin=0.0, vmax=st["vmax"], interpolation="nearest")
        for g in islets:
            c = np.asarray(shapely.get_coordinates(g.exterior)) / s_um - 0.5
            ax.add_patch(MplPolygon(c, closed=True, fill=False, ec="#63e6be", lw=0.45))
        ax.set_xticks([]); ax.set_yticks([])
        role = st["recorded_role"] or "-"
        ax.set_title(f"{st['channel']}\nrecorded: {role}\n"
                     f"islet {st['islet_mean']:.2f} / elsewhere {st['elsewhere_mean']:.2f}"
                     f"  =  {st['ratio']:.1f}x", fontsize=9)
        if j == 0:
            ax.set_ylabel("whole section", fontsize=10)

    i_um = geo.PIXEL_UM * (shapes[0][0] // shapes[ISLET_LEVEL][0])
    for r, g in enumerate(picks, start=1):
        cx, cy = g.centroid.x, g.centroid.y
        x0 = int((cx - ISLET_HALF_UM) / i_um); x1 = int((cx + ISLET_HALF_UM) / i_um)
        y0 = int((cy - ISLET_HALF_UM) / i_um); y1 = int((cy + ISLET_HALF_UM) / i_um)
        ring = np.asarray(shapely.get_coordinates(g.exterior)) / i_um - 0.5
        for j, st in enumerate(stats):
            ax = axes[r, j]
            crop = _plane(qptiff, j, ISLET_LEVEL, box=(x0, y0, x1, y1))
            ax.imshow(crop, cmap="gray", vmin=0.0, vmax=st["vmax"],
                      interpolation="nearest", extent=[x0, x1, y1, y0])
            ax.add_patch(MplPolygon(ring, closed=True, fill=False, ec="#63e6be", lw=1.3))
            ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
            ax.set_xticks([]); ax.set_yticks([])
            if j == 0:
                ax.set_ylabel(f"islet {r}\n{g.area:.0f} µm²", fontsize=9)

    fig.suptitle(
        f"Donor {donor} · IF {serial} · every qptiff channel against the maintainer's own "
        f"QuPath islet polygons\nqptiff frame only — no registration, no Xenium, no transform. "
        f"vmin=0 and vmax=p99.5 over ALL tissue in that channel (shared down each column).",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"channel_check_{donor}_{serial}.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    print(f"[channel-check] {p}", flush=True)
    for st in stats:
        print(f"    {st['channel']:<10} recorded {st['recorded_role']:<8} "
              f"islet {st['islet_mean']:8.2f}  elsewhere {st['elsewhere_mean']:8.2f}  "
              f"{st['ratio']:6.2f}x", flush=True)
    return {"donor_id": donor, "serial": serial, "qptiff": str(qptiff),
            "n_islets": len(islets), "channels": stats, "figure": str(p)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT / "channel_check")
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    recs = []
    for serial, gj in discover(args.geojson_dir).items():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        recs.append(check_one(serial, gj, sections, args.out))
    if not recs:
        print("no donor matched", file=sys.stderr)
        return 1
    p = args.out / "channel_check.json"
    p.write_text(json.dumps(recs, indent=2))
    print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
