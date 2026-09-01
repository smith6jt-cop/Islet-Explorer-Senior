#!/usr/bin/env python
"""Does the post-Xenium IF line up with its Xenium section at CELL scale?

Runs BEFORE any per-cell measurement, and its verdict gates that work:

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_align.py
    ../../.venv/bin/python scripts/restain/run_align.py --donors 6442 --windows islets

**Probes tile the whole TISSUE by default, not the islets.** Islets are where the hormone
question lives, but they are a biased and sparse sample of the section: on donor 6436 only 19
of 125 islet windows produced a peak, and every islet is by definition endocrine tissue in a
particular part of the section. A residual field estimated at islets alone is fitted where the
islets happen to be. A regular grid over the tissue answers "is this section aligned?" instead
of "are these islets aligned?", gives hundreds of probes instead of tens, and uses every cell
in the window rather than only the ones inside an annotation.

Exit 0 PASS, 3 INCONCLUSIVE, 4 FAIL, 1 I/O. See ``align.py`` for what is measured and why the
thresholds are geometry rather than preference.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import align as AL                                                # noqa: E402
import geo                                                        # noqa: E402
import hormone as HZ                                              # noqa: E402
import tiles as TL                                                # noqa: E402
import xenio                                                      # noqa: E402
from rasterize_labels import rasterize_binary                     # noqa: E402
from run_hormone import IMAGES, SECTIONS                          # noqa: E402
from run_register import (GEOJSON_DIR, OUT, XENIUM_PATHS,         # noqa: E402
                          discover, donor_of, register_one)

from phenocycler.integration import boundaries as B               # noqa: E402

#: Probe window side. Must comfortably exceed 2x the search radius so the correlation disk fits
#: inside the window; at 60 um search that is 120 um, and 500 um leaves room for the nuclei
#: themselves to carry the correlation rather than the window edge.
WINDOW_UM = 500.0
#: Step between probe centres. Equal to the window, so probes tile without overlap and each
#: measurement is independent.
WINDOW_STEP_UM = 500.0


def _median_cell_diameter_um(bundle: Path) -> float:
    """Equivalent diameter of the donor's own median cell -- the yardstick every offset is
    reported against, taken from this donor's data rather than from a literature figure."""
    import pyarrow.parquet as pq
    a = pq.read_table(bundle / "cells.parquet", columns=["cell_area"]).column(0).to_numpy()
    return float(2.0 * np.sqrt(np.nanmedian(a) / np.pi))


def tissue_windows(tissue, *, size_um: float = WINDOW_UM,
                   step_um: float = WINDOW_STEP_UM) -> list[dict]:
    """Probe windows tiling the tissue polygon, in IF microns.

    A tile is kept when its CENTRE is inside the tissue, so the probe set follows the section's
    real shape rather than its bounding box -- a pancreas section fills well under half of its
    own bbox, and probing the empty part would produce a page of unavailable reasons.
    """
    from shapely.geometry import Point

    x0, y0, x1, y1 = tissue.bounds
    half = size_um / 2.0
    out = []
    ys = np.arange(y0 + half, y1, step_um)
    xs = np.arange(x0 + half, x1, step_um)
    for cy in ys:
        for cx in xs:
            if not tissue.contains(Point(cx, cy)):
                continue
            out.append({"index": len(out), "kind": "tissue",
                        "cx_um": float(cx), "cy_um": float(cy),
                        "bounds": (cx - half, cy - half, cx + half, cy + half)})
    return out


def islet_windows(islets, *, pad_um: float = AL.ISLET_BUFFER_UM) -> list[dict]:
    """One window per kept islet -- the original probe set, kept for comparison."""
    out = []
    for k, g in enumerate(islets):
        x0, y0, x1, y1 = g.bounds
        out.append({"index": k, "kind": "islet",
                    "cx_um": float(g.centroid.x), "cy_um": float(g.centroid.y),
                    "bounds": (x0 - pad_um, y0 - pad_um, x1 + pad_um, y1 + pad_um)})
    return out


def _valis_transform(donor: str, roi: str, out_dir: Path):
    """The VALIS-derived affine (plus displacement field, if the mapping was not affine)."""
    from phenocycler.integration.rasterize import GridSpec
    from phenocycler.integration.transform import Transform

    p = out_dir / f"valis_registration_{roi}.json"
    if not p.exists():
        raise AL.AlignError(f"{p} is absent -- run valis_register.py then valis_transform.py")
    for rec in json.loads(p.read_text()):
        if rec["donor_id"] != donor:
            continue
        affine = np.asarray(rec["affine"], dtype=np.float64)
        d = rec.get("displacement")
        if not d:
            return Transform(affine=affine)
        z = np.load(d["file"])
        return Transform(affine=affine, displacement=z["field"],
                         grid=GridSpec.from_dict(d["grid"]))
    raise AL.AlignError(f"{p} carries no donor {donor}")


def _refined_transform(donor: str, roi: str, out_dir: Path):
    """The transform ``run_refine.py`` chose, loaded from the donor's own directory.

    Written by ``Transform.save``, so the affine, the grid and the displacement field come back
    exactly as they were fitted -- nothing here re-derives or re-scales any of them.
    """
    from phenocycler.integration.transform import Transform

    d = out_dir / "valis" / f"donor_id={donor}"
    if not (d / "refined.json").exists():
        raise AL.AlignError(f"{d/'refined.json'} is absent -- run run_refine.py first")
    return Transform.load(d, stem="refined")


def align_one(serial: str, gj: Path, sections: pd.DataFrame, *, roi: str = "panc",
              mode: str = "nucleus", windows: str = "tissue",
              transform: str = "ecc", out_dir: Path = OUT) -> dict:
    import shapely
    import tifffile

    # Resolve the transform BEFORE registering, so `register_one` can skip a fit it would only
    # discard and can measure Dice on the transform actually under test. Same gate, same probes,
    # same nuclei -- only the transform changes, so the comparison stays a measurement.
    donor = donor_of(serial, XENIUM_PATHS)
    if transform == "ecc":
        supplied = None
    elif transform == "valis":
        supplied = _valis_transform(donor, roi, out_dir)
    elif transform == "refined":
        supplied = _refined_transform(donor, roi, out_dir)
    else:
        raise AL.AlignError(f"unknown transform source {transform!r}")
    rec, res, dd, isl, _grid = register_one(serial, gj, roi=roi, transform=supplied)
    row = sections[sections.xenium_serial == serial]
    if row.empty:
        raise HZ.HormoneError(f"serial {serial} is not recorded in {SECTIONS}")
    row = row.iloc[0]
    qptiff = IMAGES / row.image_file
    if not qptiff.exists():
        raise HZ.HormoneError(f"{qptiff} is absent")

    pixel_um = float(row.pixel_size_um)
    if pixel_um != geo.PIXEL_UM:
        raise HZ.HormoneError(
            f"{row.image_file} records pixel_size_um={pixel_um} but annotations were scaled "
            f"by geo.PIXEL_UM={geo.PIXEL_UM}. These must be the same number.")

    bundle = xenio._bundle_dir(donor, XENIUM_PATHS, roi)
    diameter_um = _median_cell_diameter_um(bundle)

    if windows == "tissue":
        probes = tissue_windows(geo.tissue_polygon(geo.load_annotations(gj)))
    elif windows == "islets":
        probes = islet_windows(isl)
    else:
        raise AL.AlignError(f"unknown window set {windows!r}")

    if mode == "nucleus":
        polys, _cell_ids = B.load_xenium_nucleus_polygons(bundle)
    else:
        polys = B.load_xenium_polygons(bundle)
    arr = np.asarray(polys.values, dtype=object)

    # Warp FORWARD, Xenium -> IF. `Transform.invert()` drops any displacement field, so the
    # inverse direction is not merely slower here, it is a different transform.
    cen_if = res.transform.apply(shapely.get_coordinates(shapely.centroid(arr)))

    offsets: list[AL.Offset] = []
    contrasts: list[float] = []
    halves: list[float] = []
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series]))
        shapes = [(lv.shape[-2], lv.shape[-1]) for lv in levels]
        lvl = TL.level_for(shapes, pixel_um, AL.ALIGN_TARGET_UM)
        lvl_um = pixel_um * (shapes[0][0] / shapes[lvl][0])
        page = levels[lvl].pages[HZ.channel_index(qptiff, "NUCLEUS", row)]

        for w in probes:
            x0, y0, x1, y1 = w["bounds"]
            base = dict(islet_index=w["index"], cx_um=w["cx_um"], cy_um=w["cy_um"])
            grid, _px = TL.window_grid((x0, y0, x1, y1), lvl_um)
            inside = ((cen_if[:, 0] >= x0) & (cen_if[:, 0] <= x1)
                      & (cen_if[:, 1] >= y0) & (cen_if[:, 1] <= y1))
            n_nuc = int(inside.sum())
            if n_nuc == 0:
                offsets.append(AL.Offset(**base, n_nuclei=0,
                                         reason="no nuclei in this window"))
                continue
            warped = shapely.transform(arr[inside], lambda a: res.transform.apply(a))
            dapi, valid = TL.read_window_on_grid(page, grid, lvl_um)
            template = rasterize_binary(list(warped), grid)

            off = AL.offset_by_ncc(dapi, valid, template, pixel_um=lvl_um,
                                   islet_index=w["index"], n_nuclei=n_nuc)
            c = AL.contrast_inside_vs_gaps(dapi, valid, template)
            if c is not None:
                contrasts.append(c)

            # The check that decides whether the offset is REAL. Contrast at the as-registered
            # position says nothing on its own -- if the registration is off, of course DAPI is
            # not brighter inside the warped nuclei. But if the measured shift is the true
            # residual, applying it must make DAPI brighter inside them.
            c_shift = half = None
            if off.reason == "measured":
                c_shift = AL.contrast_inside_vs_gaps(
                    dapi, valid, template,
                    shift_px=(int(round(off.dy_um / lvl_um)),
                              int(round(off.dx_um / lvl_um))))
                if n_nuc >= 2 * AL.MIN_NUCLEI_PER_ISLET:
                    a_off = AL.offset_by_ncc(dapi, valid,
                                             rasterize_binary(list(warped[0::2]), grid),
                                             pixel_um=lvl_um, n_nuclei=n_nuc // 2)
                    b_off = AL.offset_by_ncc(dapi, valid,
                                             rasterize_binary(list(warped[1::2]), grid),
                                             pixel_um=lvl_um, n_nuclei=n_nuc - n_nuc // 2)
                    if a_off.reason == "measured" and b_off.reason == "measured":
                        half = float(np.hypot(a_off.dx_um - b_off.dx_um,
                                              a_off.dy_um - b_off.dy_um))
                        halves.append(half)
            offsets.append(AL.Offset(**{**off.__dict__, **base, "contrast": c,
                                        "contrast_shifted": c_shift,
                                        "half_split_um": half}))

    measured = [o for o in offsets if o.reason == "measured"]
    mags = np.array([o.magnitude_um for o in measured], dtype=float)
    n_total = len(offsets)
    median_um = float(np.median(mags)) if mags.size else None
    p90_um = float(np.percentile(mags, 90)) if mags.size else None

    # Bounding median: every unmeasured probe entered at the search ceiling. Unavailability
    # preferentially removes BADLY aligned windows, so a survivors-only median can only flatter.
    bounded = np.concatenate([mags, np.full(n_total - len(measured), AL.SEARCH_RADIUS_UM)])
    median_bounded = float(np.median(bounded)) if bounded.size else None

    shifted = [o.contrast_shifted for o in measured if o.contrast_shifted is not None]
    paired = [(o.contrast, o.contrast_shifted) for o in measured
              if o.contrast is not None and o.contrast_shifted is not None]
    contrast_fraction = (float(np.mean(np.asarray(contrasts) >= AL.MIN_CONTRAST))
                         if contrasts else 0.0)
    precision_um = float(np.median(halves)) if halves else None
    status, why = AL.verdict(median_um, diameter_um=diameter_um, n_measured=len(measured),
                             n_total=n_total, contrast_fraction=contrast_fraction,
                             mode=mode, precision_um=precision_um)

    reasons: dict[str, int] = {}
    for o in offsets:
        reasons[o.reason] = reasons.get(o.reason, 0) + 1

    out = dict(rec)
    out["alignment"] = {
        "mode": mode, "window_set": windows, "transform_source": transform,
        "affine": np.asarray(res.transform.affine).tolist(),
        "window_um": WINDOW_UM if windows == "tissue" else None,
        "qptiff": str(qptiff), "pyramid_level": int(lvl),
        "align_um_per_px": float(lvl_um),
        "certifies_mapping": {"pixel_um_full_res": pixel_um,
                              "source": "restain_sections.csv / qptiff XML PixelSizeMicrons"},
        "sign_convention": "shift to ADD to the warped Xenium mask to land it on the IF DAPI",
        "median_cell_diameter_um": diameter_um,
        "threshold_pass_um": AL.self_overlap_offset_um(diameter_um, 0.75),
        "threshold_fail_um": AL.self_overlap_offset_um(diameter_um, 0.50),
        "n_windows": n_total, "n_measured": len(measured),
        "n_truncated": int(sum(1 for o in measured if o.truncated)),
        "truncated_fraction": (float(np.mean([o.truncated for o in measured]))
                               if measured else 0.0),
        "measured_fraction": (len(measured) / n_total) if n_total else 0.0,
        "median_offset_um": median_um, "p90_offset_um": p90_um,
        "median_offset_um_bounded": median_bounded,
        "estimator_precision_um": precision_um,
        "contrast_median": float(np.median(contrasts)) if contrasts else None,
        "contrast_median_after_shift": float(np.median(shifted)) if shifted else None,
        "contrast_improved_fraction": (float(np.mean([b > a for a, b in paired]))
                                       if paired else None),
        "reasons": reasons,
        "verdict": status, "verdict_why": why,
        "constants": {"ALIGN_TARGET_UM": AL.ALIGN_TARGET_UM,
                      "SEARCH_RADIUS_UM": AL.SEARCH_RADIUS_UM,
                      "HIGHPASS_SIGMA_UM": AL.HIGHPASS_SIGMA_UM,
                      "MIN_NCC_PEAK": AL.MIN_NCC_PEAK, "MIN_PEAK_Z": AL.MIN_PEAK_Z,
                      "MIN_NUCLEI_PER_ISLET": AL.MIN_NUCLEI_PER_ISLET,
                      "ISLET_BUFFER_UM": AL.ISLET_BUFFER_UM,
                      "TRUNCATION_FRACTION_OF_R": AL.TRUNCATION_FRACTION_OF_R,
                      "WINDOW_UM": WINDOW_UM, "WINDOW_STEP_UM": WINDOW_STEP_UM,
                      "MIN_MEASURED_FRACTION": AL.MIN_MEASURED_FRACTION,
                      "MIN_CONTRAST": AL.MIN_CONTRAST,
                      "MIN_CONTRAST_ISLET_FRACTION": AL.MIN_CONTRAST_ISLET_FRACTION},
        "windows": [o.as_dict() for o in offsets],
    }
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--mode", default="nucleus", choices=["nucleus", "cell"])
    ap.add_argument("--windows", default="tissue", choices=["tissue", "islets"])
    ap.add_argument("--transform", default="ecc", choices=["ecc", "valis", "refined"],
                    help="which registration to measure: the tissue-outline ECC fit, the raw "
                         "VALIS DAPI-to-DAPI fit, or the nucleus-matched transform run_refine.py "
                         "chose")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    found = discover(args.geojson_dir)
    if not found:
        print(f"no geojson exports under {args.geojson_dir}", file=sys.stderr)
        return 1

    rows, worst = [], "PASS"
    for serial, gj in sorted(found.items()):
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in set(args.donors):
            continue
        rec = align_one(serial, gj, sections, roi=args.roi, mode=args.mode,
                        windows=args.windows, transform=args.transform, out_dir=args.out)
        a = rec["alignment"]
        rows.append(rec)
        med = "n/a" if a["median_offset_um"] is None else f"{a['median_offset_um']:.2f}"
        print(f"[{rec['donor_id']}] {a['verdict']}: median {med} um "
              f"(pass<={a['threshold_pass_um']:.2f}, fail>={a['threshold_fail_um']:.2f}), "
              f"[{a['transform_source']}] "
              f"{a['n_measured']}/{a['n_windows']} {a['window_set']} windows measured, "
              f"{a['n_truncated']} at the search rim | {a['verdict_why']}", flush=True)
        order = {"PASS": 0, "INCONCLUSIVE": 1, "FAIL": 2}
        if order[a["verdict"]] > order[worst]:
            worst = a["verdict"]

    if not rows:
        print("no donor matched", file=sys.stderr)
        return 1
    args.out.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.windows == "tissue" else f"_{args.windows}"
    if args.transform != "ecc":
        suffix += f"_{args.transform}"
    p = args.out / f"alignment_{args.roi}{suffix}.json"
    p.write_text(json.dumps(rows, indent=2, allow_nan=False))
    print(f"-> {p}", flush=True)
    return {"PASS": 0, "INCONCLUSIVE": 3, "FAIL": 4}[worst]


if __name__ == "__main__":
    raise SystemExit(main())
