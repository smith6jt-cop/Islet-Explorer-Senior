#!/usr/bin/env python
"""Pictures for the alignment gate, so the verdict can be checked rather than believed.

    cd core/multiplex5000
    ../../.venv/bin/python scripts/restain/align_figures.py --donors 6442

Two sheets per donor, into ``data/development/restain/figures/donor_id=<donor>/``:

``alignment_field.png``   every islet's residual as an arrow on the section, with the
                          pass/fail geometry stated in microns.
``alignment_windows.png`` four islets -- median, p90, worst and one the gate could not
                          measure -- each as three panels: the IF DAPI with the warped Xenium
                          nuclei drawn AS REGISTERED, the same DAPI with the nuclei shifted by
                          the measured offset, and the correlation surface itself.

The middle panel is the one that settles whether the number is real. If the nuclei sit off the
DAPI on the left and land on it in the middle, the residual is a property of the registration.
If neither lines up, the measurement is wrong and the verdict should not be trusted.

Every image panel is windowed on the measured floor, 0, through ``display.py``; the window and
the fraction of pixels at the floor are printed on the panel.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib                                                 # noqa: E402
matplotlib.use("Agg")                       # headless: no display on this box

import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import align as AL                                                # noqa: E402
import display as DP                                              # noqa: E402
import geo                                                        # noqa: E402
import hormone as HZ                                              # noqa: E402
import tiles as TL                                                # noqa: E402
import xenio                                                      # noqa: E402
from rasterize_labels import rasterize_binary                     # noqa: E402
from run_hormone import IMAGES, SECTIONS                          # noqa: E402
from run_register import (GEOJSON_DIR, OUT, XENIUM_PATHS,         # noqa: E402
                          discover, donor_of, register_one)

from phenocycler.integration import boundaries as B               # noqa: E402

OUTLINE = "#00FFFF"          # nuclei, matching the repo's nucleus colour
SHIFTED = "#FFD700"          # the same nuclei after the measured shift
SCALE_BAR_UM = 50.0


def _scale_bar(ax, extent, um=SCALE_BAR_UM, color="white"):
    """A scale bar. Nothing in this repo draws one, so this sets the convention."""
    x0, x1, y1, y0 = extent
    span = x1 - x0
    xa = x0 + 0.06 * span
    ya = y1 - 0.08 * (y1 - y0)
    ax.plot([xa, xa + um], [ya, ya], color=color, lw=2.5, solid_capstyle="butt")
    ax.text(xa + um / 2, ya - 0.02 * (y1 - y0), f"{um:g} um", color=color,
            ha="center", va="bottom", fontsize=7)


def _rings(geom):
    parts = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
    for p in parts:
        yield np.asarray(p.exterior.coords)


def _records(a: dict) -> list[dict]:
    """Probe records, under whichever key this run wrote."""
    return a.get("windows", a.get("islets", []))


def field_figure(rec: dict, out_dir: Path) -> Path:
    a = rec["alignment"]
    isl = _records(a)
    meas = [o for o in isl if o["reason"] == "measured"]
    fig, ax = plt.subplots(figsize=(9, 11), dpi=150)

    if meas:
        x = np.array([o["cx_um"] for o in meas]); y = np.array([o["cy_um"] for o in meas])
        dx = np.array([o["dx_um"] for o in meas]); dy = np.array([o["dy_um"] for o in meas])
        mag = np.hypot(dx, dy)
        q = ax.quiver(x, y, dx, dy, mag, angles="xy", scale_units="xy",
                      scale=1.0 / 40.0, cmap="magma", width=0.004)
        cb = fig.colorbar(q, ax=ax, fraction=0.03, pad=0.02)
        cb.set_label("residual offset (um), arrows x40")
    un = [o for o in isl if o["reason"] != "measured"]
    if un and "cx_um" in (un[0] if un else {}):
        ax.scatter([o["cx_um"] for o in un], [o["cy_um"] for o in un], s=14,
                   c=DP.UNAVAILABLE, marker="x", label=f"not measurable (n={len(un)})")
        ax.legend(loc="lower right", fontsize=8)

    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xlabel("x (um, IF frame)"); ax.set_ylabel("y (um, IF frame)")
    source = {"ecc": "the tissue-outline registration",
              "valis": "the VALIS DAPI-to-DAPI registration",
              "refined": "the nucleus-matched registration"}.get(
                  a.get("transform_source") or "ecc", a.get("transform_source"))
    ax.set_title(
        f"donor {rec['donor_id']} - residual after {source}, "
        f"per {a.get('window_set', 'islet')} window\n"
        f"{a['verdict']}: median {a['median_offset_um']:.2f} um vs a "
        f"{a['median_cell_diameter_um']:.2f} um median cell "
        f"(75% self-overlap at {a['threshold_pass_um']:.2f}, 50% at {a['threshold_fail_um']:.2f})",
        fontsize=10)
    ax.text(0.01, -0.06,
            f"{a['n_measured']}/{a.get('n_windows', a.get('n_islets'))} "
            f"{a.get('window_set', 'islet')} windows measurable; "
            f"{a['n_truncated']} peaks at the {a['constants']['SEARCH_RADIUS_UM']:g} um search rim; "
            f"split-half precision {a['estimator_precision_um']:.3f} um; "
            f"DAPI is {a['contrast_median']:.1f}x brighter inside the warped nuclei than "
            f"between them; tissue Dice {rec['dice_before']:.3f} -> {rec['dice_after']:.3f}. "
            f"n=2 sections.",
            transform=ax.transAxes, fontsize=7.5, va="top")
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / "alignment_field.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[figures] {p}", flush=True)
    return p


def window_figure(rec: dict, picks: list[dict], windows: dict, out_dir: Path) -> Path:
    a = rec["alignment"]
    n = len(picks)
    fig, axes = plt.subplots(n, 3, figsize=(13.5, 4.6 * n), dpi=150, squeeze=False)
    for r, o in enumerate(picks):
        w = windows[o["islet_index"]]
        dapi, valid, warped, grid, cc, lvl_um = (
            w["dapi"], w["valid"], w["warped"], w["grid"], w["cc"], w["lvl_um"])
        rng = DP.unsigned_display_range(dapi[valid])
        extent = (grid.origin_x_um, grid.origin_x_um + grid.width * grid.pixel_um,
                  grid.origin_y_um + grid.height * grid.pixel_um, grid.origin_y_um)

        for c, (title, shift) in enumerate((
                ("as registered", (0.0, 0.0)),
                ("shifted by the measured offset",
                 (o.get("dx_um") or 0.0, o.get("dy_um") or 0.0)))):
            ax = axes[r][c]
            ax.imshow(dapi, cmap="gray", vmin=rng.vmin, vmax=rng.vmax, extent=extent,
                      interpolation="nearest")
            for g in warped:
                for ring in _rings(g):
                    ax.plot(ring[:, 0] + shift[0], ring[:, 1] + shift[1],
                            color=SHIFTED if c else OUTLINE, lw=0.45)
            ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
            ax.set_aspect("equal"); ax.set_axis_off()
            _scale_bar(ax, extent)
            ax.set_title(f"window {o['islet_index']} ({o['reason']}) - DAPI, nuclei {title}\n"
                         f"{rng.caption()}", fontsize=7.5)

        ax = axes[r][2]
        R = a["constants"]["SEARCH_RADIUS_UM"]
        im = ax.imshow(cc, cmap="magma", extent=(-R, R, R, -R), interpolation="nearest")
        ax.add_patch(plt.Circle((0, 0), R, fill=False, color="white", lw=0.8, ls="--"))
        if o.get("dx_um") is not None:
            ax.plot([o["dx_um"]], [o["dy_um"]], marker="+", color="#00FF00", ms=12, mew=2)
        ax.set_aspect("equal")
        ax.set_xlabel("dx (um)"); ax.set_ylabel("dy (um)")
        peak = "n/a" if o.get("ncc_peak") is None else f"{o['ncc_peak']:.3f}"
        mag = "n/a" if o.get("magnitude_um") is None else f"{o['magnitude_um']:.2f} um"
        ax.set_title(f"correlation surface - peak {peak}, |d| {mag}", fontsize=7.5)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label("normalised correlation")

    fig.suptitle(
        f"donor {rec['donor_id']} - does the Xenium nucleus segmentation land on the IF DAPI?\n"
        f"cyan = as registered, gold = after the measured shift. Black point is 0 on every panel.",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / "alignment_windows.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[figures] {p}", flush=True)
    return p


def _pick(islets: list[dict]) -> list[dict]:
    """Median, p90, worst, and one the gate could not measure.

    Deterministic and stated: the reader sees the typical case, the tail, the worst case, and
    what an exclusion looks like -- not four islets chosen after seeing them.
    """
    meas = sorted((o for o in islets if o["reason"] == "measured"),
                  key=lambda o: o["magnitude_um"])
    picks = []
    if meas:
        picks.append(dict(meas[len(meas) // 2], _label="median"))
        picks.append(dict(meas[min(len(meas) - 1, int(0.9 * (len(meas) - 1)))], _label="p90"))
        picks.append(dict(meas[-1], _label="worst"))
    un = [o for o in islets if o["reason"] != "measured"]
    if un:
        picks.append(dict(max(un, key=lambda o: o["n_nuclei"]), _label="not measurable"))
    seen, out = set(), []
    for p in picks:
        if p["islet_index"] not in seen:
            seen.add(p["islet_index"]); out.append(p)
    return out


def windows_for(rec: dict, picks: list[dict], sections, *, transform: str = "ecc",
                out_dir: Path = OUT) -> dict:
    """Rebuild exactly the crops the verdict was taken from, for the picked islets."""
    import shapely
    import tifffile

    serial = rec["serial"]
    gj = discover(GEOJSON_DIR)[serial]
    # Draw with the transform the verdict was taken from, or the picture is of a different
    # registration than the number beside it. Resolved first so the tissue-outline fit is not
    # re-run just to be discarded.
    supplied = None
    if transform == "refined":
        from run_align import _refined_transform
        supplied = _refined_transform(rec["donor_id"], rec["roi"], out_dir)
    elif transform == "valis":
        from run_align import _valis_transform
        supplied = _valis_transform(rec["donor_id"], rec["roi"], out_dir)
    _rec2, res, _dd, isl, _g = register_one(serial, gj, roi=rec["roi"], transform=supplied)
    row = sections[sections.xenium_serial == serial].iloc[0]
    qptiff = IMAGES / row.image_file
    pixel_um = float(row.pixel_size_um)
    bundle = xenio._bundle_dir(rec["donor_id"], XENIUM_PATHS, rec["roi"])

    polys, _cids = B.load_xenium_nucleus_polygons(bundle)
    arr = np.asarray(polys.values, dtype=object)
    cen_if = res.transform.apply(shapely.get_coordinates(shapely.centroid(arr)))

    out = {}
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series]))
        shapes = [(lv.shape[-2], lv.shape[-1]) for lv in levels]
        lvl = TL.level_for(shapes, pixel_um, AL.ALIGN_TARGET_UM)
        lvl_um = pixel_um * (shapes[0][0] / shapes[lvl][0])
        page = levels[lvl].pages[HZ.channel_index(qptiff, "NUCLEUS", row)]
        half = (rec["alignment"].get("window_um") or 0.0) / 2.0
        for o in picks:
            k = o["islet_index"]
            if half:                      # tissue tiles: the probe IS the window
                cx, cy = o["cx_um"], o["cy_um"]
                bounds = (cx - half, cy - half, cx + half, cy + half)
            else:                         # islet windows: bbox plus the recorded buffer
                x0, y0, x1, y1 = isl[k].bounds
                pad = AL.ISLET_BUFFER_UM
                bounds = (x0 - pad, y0 - pad, x1 + pad, y1 + pad)
            x0, y0, x1, y1 = bounds
            grid, _px = TL.window_grid(bounds, lvl_um)
            inside = ((cen_if[:, 0] >= x0) & (cen_if[:, 0] <= x1)
                      & (cen_if[:, 1] >= y0) & (cen_if[:, 1] <= y1))
            warped = shapely.transform(arr[inside], lambda a: res.transform.apply(a))
            dapi, valid = TL.read_window_on_grid(page, grid, lvl_um)
            template = rasterize_binary(list(warped), grid)
            surf, _why = AL.ncc_surface(dapi, valid, template, pixel_um=lvl_um)
            R = AL.SEARCH_RADIUS_UM
            if surf is not None:
                r = int(round(R / lvl_um))
                cy, cx = surf.shape[0] // 2, surf.shape[1] // 2
                surf = surf[max(0, cy - r):cy + r + 1, max(0, cx - r):cx + r + 1]
            out[k] = {"dapi": dapi, "valid": valid, "warped": list(warped), "grid": grid,
                      "cc": surf if surf is not None else np.zeros((3, 3)), "lvl_um": lvl_um}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--transform", default="ecc", choices=["ecc", "valis", "refined"])
    args = ap.parse_args(argv)

    suffix = "" if args.transform == "ecc" else f"_{args.transform}"
    src = args.out / f"alignment_{args.roi}{suffix}.json"
    if not src.exists():
        print(f"{src} is absent -- run run_align.py first", file=sys.stderr)
        return 1
    sections = HZ.load_section_map(SECTIONS)
    made = []
    for rec in json.loads(src.read_text()):
        if args.donors and rec["donor_id"] not in set(args.donors):
            continue
        d = args.out / "figures" / f"donor_id={rec['donor_id']}"
        if args.transform != "ecc":
            d = d / args.transform
        made.append(field_figure(rec, d))
        picks = _pick(_records(rec["alignment"]))
        if picks:
            made.append(window_figure(rec, picks,
                                      windows_for(rec, picks, sections,
                                                  transform=args.transform, out_dir=args.out),
                                      d))
    if not made:
        print("no donor matched", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
