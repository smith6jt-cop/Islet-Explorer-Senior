#!/usr/bin/env python
"""Does the post-Xenium IF line up with its Xenium section at CELL scale?

Runs BEFORE any per-cell measurement, and its verdict gates that work:

    cd core/multiplex5000
    ../../.venv/bin/python scripts/restain/run_align.py --transform refined
    ../../.venv/bin/python scripts/restain/run_align.py --donors 6442 --windows islets

**Probes tile the whole TISSUE, not the islets.** Islets are where the hormone
question lives, but they are a biased and sparse sample of the section: on donor 6436 only 19
of 125 islet windows produced a peak, and every islet is by definition endocrine tissue in a
particular part of the section. A residual field estimated at islets alone is fitted where the
islets happen to be. A regular grid over the tissue answers "is this section aligned?" instead
of "are these islets aligned?", gives hundreds of probes instead of tens, and uses every cell
in the window rather than only the ones inside an annotation.  When a section has no QuPath
GeoJSON, the same 500 um grid is anchored at the origin and only bins occupied by transformed
Xenium nuclei are retained.  That fallback uses no IF intensity and therefore cannot select
windows because they flatter the DAPI result it is about to measure.

Exit 0 PASS, 3 INCONCLUSIVE, 4 FAIL, 1 I/O. See ``align.py`` for what is measured and why the
thresholds are geometry rather than preference.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import align as AL  # noqa: E402
import geo  # noqa: E402
import registration_inputs as RI  # noqa: E402
import tiles as TL  # noqa: E402
from rasterize_labels import rasterize_binary  # noqa: E402
from run_register import OUT, register_one  # noqa: E402

from phenocycler.integration import boundaries as B  # noqa: E402
# The one file checksum (phenocycler/artifacts.py): every byte, hex digest. Bound to the
# local name because other modules and tests call it through this one.
from phenocycler.artifacts import sha256_file as _sha256  # noqa: E402

#: Probe window side. Must comfortably exceed 2x the search radius so the correlation disk fits
#: inside the window; at 60 um search that is 120 um, and 500 um leaves room for the nuclei
#: themselves to carry the correlation rather than the window edge.
WINDOW_UM = 500.0
#: Step between probe centres. Equal to the window, so probes tile without overlap and each
#: measurement is independent.
WINDOW_STEP_UM = 500.0


def _json_write(path: Path, value) -> Path:
    """Atomically write one strict JSON record.

    The per-donor file is a release gate consumed by the exporter.  A reader must see either
    the old complete verdict or the new complete verdict, never half of a JSON document.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False))
    tmp.replace(path)
    return path


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


def occupied_windows(xy_if_um: np.ndarray, *, size_um: float = WINDOW_UM,
                     step_um: float = WINDOW_STEP_UM) -> list[dict]:
    """Non-overlapping, occupied probe windows from transformed Xenium nuclei.

    The grid is anchored at micron zero, rather than at the observed point minimum, so its
    placement cannot change with point order or with a nucleus appearing at the tissue edge.
    Points are grouped once by integer bin -- O(n log n), with no ``n_windows * n_points``
    containment sweep in the measurement loop.  ``member_indices`` is deliberately an
    implementation detail and is omitted from the JSON alignment record.

    The registration gate's declared 500 um windows have ``size_um == step_um``.  Refusing a
    different combination keeps the precomputed membership exact: overlapping or gapped bins
    would need a different spatial-index contract rather than looking approximately right.
    """
    if not (float(size_um) > 0.0 and float(step_um) > 0.0):
        raise AL.AlignError(
            f"occupied-window size and step must be positive, got {size_um}, {step_um}")
    if not np.isclose(float(size_um), float(step_um), rtol=0.0, atol=1e-12):
        raise AL.AlignError(
            "occupied windows require size_um == step_um so every nucleus has exactly one bin")

    xy = np.asarray(xy_if_um, dtype=np.float64)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise AL.AlignError(f"nucleus coordinates must have shape (n, 2), got {xy.shape}")
    finite_idx = np.flatnonzero(np.isfinite(xy).all(axis=1))
    if finite_idx.size == 0:
        return []

    step = float(step_um)
    bx = np.floor(xy[finite_idx, 0] / step).astype(np.int64)
    by = np.floor(xy[finite_idx, 1] / step).astype(np.int64)
    # y, then x, then original point index: stable output and stable split-half membership.
    order = np.lexsort((finite_idx, bx, by))
    finite_idx, bx, by = finite_idx[order], bx[order], by[order]
    starts = np.r_[0, 1 + np.flatnonzero((bx[1:] != bx[:-1]) | (by[1:] != by[:-1]))]
    stops = np.r_[starts[1:], len(finite_idx)]

    half = float(size_um) / 2.0
    out = []
    for start, stop in zip(starts, stops, strict=True):
        ix, iy = int(bx[start]), int(by[start])
        cx, cy = (ix + 0.5) * step, (iy + 0.5) * step
        out.append({"index": len(out), "kind": "occupied",
                    "cx_um": float(cx), "cy_um": float(cy),
                    "bounds": (cx - half, cy - half, cx + half, cy + half),
                    "member_indices": finite_idx[start:stop]})
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
    return _saved_transform(donor, "refined", out_dir)


def _saved_transform(donor: str, stem: str, out_dir: Path):
    """Load one explicitly named, donor-local transform without any fallback.

    A gate and an export must name the same serialization.  In particular, an
    ``intensity_refined`` request must never quietly fall back to ``refined`` merely because
    the former is absent.
    """
    from phenocycler.integration.transform import Transform

    if not stem or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-"
                       for c in stem) or stem.startswith("."):
        raise AL.AlignError(f"unsafe saved-transform stem {stem!r}")
    d = out_dir / "valis" / f"donor_id={donor}"
    js = d / f"{stem}.json"
    if not js.exists():
        raise AL.AlignError(f"{js} is absent -- produce the {stem!r} transform first")
    payload = json.loads(js.read_text())
    npz = d / f"{stem}_displacement.npz"
    if payload.get("has_displacement") and not npz.exists():
        raise AL.AlignError(f"{js} declares a displacement field but {npz} is absent")
    return Transform.load(d, stem=stem)


def _select_window_fold(probes: list[dict], fold: str) -> list[dict]:
    """Select a deterministic checkerboard fold from the origin-anchored probe grid."""
    if fold == "all":
        return probes
    parity = 0 if fold == "even" else 1
    selected = []
    for w in probes:
        ix = int(np.floor(float(w["cx_um"]) / WINDOW_STEP_UM))
        iy = int(np.floor(float(w["cy_um"]) / WINDOW_STEP_UM))
        if (ix + iy) % 2 == parity:
            selected.append(w)
    return selected


def align_one(record: RI.RegistrationInput, *, mode: str = "nucleus", windows: str = "auto",
              transform: str = "refined", window_fold: str = "all", out_dir: Path = OUT,
              xenium_paths: Path = RI.XENIUM_PATHS) -> dict:
    """Measure one inventory record without requiring annotations or biomarker metadata.

    ``refined`` and ``valis`` load their transforms directly.  In particular, they never call
    :func:`run_register.register_one`, whose tissue-outline path requires GeoJSON and loads
    expression data irrelevant to this DAPI gate.  ``ecc`` remains available as the legacy
    annotation-derived transform and therefore still requires a GeoJSON.
    """
    import shapely
    import tifffile

    donor, serial, roi = record.donor_id, record.serial, record.roi
    qptiff, pixel_um, bundle = record.qptiff, float(record.pixel_um), record.bundle
    annotations = geo.load_annotations(record.geojson) if record.geojson is not None else None
    annotations_have_tissue = annotations is not None and geo.has_tissue(annotations)
    resolved_windows = ("tissue" if annotations_have_tissue else "occupied") \
        if windows == "auto" else windows
    if resolved_windows in {"tissue", "islets"} and annotations is None:
        raise AL.AlignError(
            f"donor {donor} has no GeoJSON; use --windows auto or --windows occupied")
    if resolved_windows == "tissue" and not annotations_have_tissue:
        raise AL.AlignError(
            f"donor {donor} GeoJSON has no {geo.TISSUE!r} annotation; classes present: "
            f"{sorted(annotations)}; use --windows auto or --windows occupied")
    if transform == "ecc" and not annotations_have_tissue:
        classes = [] if annotations is None else sorted(annotations)
        raise AL.AlignError(
            f"donor {donor} cannot use the annotation-derived ECC transform without a "
            f"{geo.TISSUE!r} annotation; classes present: {classes}")
    if (transform == "ecc" or resolved_windows in {"tissue", "islets"}) \
            and pixel_um != geo.PIXEL_UM:
        raise AL.AlignError(
            f"{qptiff.name} is {pixel_um} um/px but GeoJSON coordinates use "
            f"geo.PIXEL_UM={geo.PIXEL_UM}; annotation-derived registration/windows are refused")

    if transform == "ecc":
        if record.geojson is None:
            raise AL.AlignError(f"donor {donor} has no GeoJSON; the ECC transform requires one")
        rec, result, _dd, islets, _grid = register_one(
            serial, record.geojson, xenium_paths=xenium_paths, roi=roi)
        if rec["donor_id"] != donor:
            raise AL.AlignError(
                f"registration inventory says donor {donor}, ECC resolved {rec['donor_id']}")
        tf = result.transform
    elif transform == "valis":
        tf = _valis_transform(donor, roi, out_dir)
        islets = geo.islets(annotations) if annotations is not None else []
        rec = {"serial": serial, "donor_id": donor, "roi": roi,
               "geojson": None if record.geojson is None else str(record.geojson)}
    else:
        tf = (_refined_transform(donor, roi, out_dir) if transform == "refined"
              else _saved_transform(donor, transform, out_dir))
        islets = geo.islets(annotations) if annotations is not None else []
        rec = {"serial": serial, "donor_id": donor, "roi": roi,
               "geojson": None if record.geojson is None else str(record.geojson)}
    rec["registration_input"] = {
        "qptiff": str(qptiff), "bundle": str(bundle),
        "pixel_um": pixel_um, "pixel_size_source": "qptiff XML PixelSizeMicrons",
        "dapi_page": int(record.dapi_page), "dapi_channel_name": record.channel_names[record.dapi_page],
        "provenance_source": record.provenance_source,
        "geojson": None if record.geojson is None else str(record.geojson),
    }
    diameter_um = _median_cell_diameter_um(bundle)

    if mode == "nucleus":
        polys, _cell_ids = B.load_xenium_nucleus_polygons(bundle)
    else:
        polys = B.load_xenium_polygons(bundle)
    arr = np.asarray(polys.values, dtype=object)

    # Warp FORWARD, Xenium -> IF. `Transform.invert()` drops any displacement field, so the
    # inverse direction is not merely slower here, it is a different transform.
    cen_if = tf.apply(shapely.get_coordinates(shapely.centroid(arr)))

    if resolved_windows == "tissue":
        probes = tissue_windows(geo.tissue_polygon(annotations))
    elif resolved_windows == "islets":
        probes = islet_windows(islets)
    elif resolved_windows == "occupied":
        if mode != "nucleus":
            raise AL.AlignError(
                "annotation-free occupied windows require --mode nucleus; the cell mask cannot "
                "define a nucleus-occupancy grid and cannot certify alignment")
        probes = occupied_windows(cen_if)
    else:
        raise AL.AlignError(f"unknown window set {resolved_windows!r}")
    if window_fold not in {"all", "even", "odd"}:
        raise AL.AlignError(f"unknown window fold {window_fold!r}")
    if window_fold != "all" and resolved_windows != "occupied":
        raise AL.AlignError("checkerboard window folds are defined only for occupied windows")
    probes = _select_window_fold(probes, window_fold)

    offsets: list[AL.Offset] = []
    contrasts: list[float] = []
    halves: list[float] = []
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series]))
        shapes = [(lv.shape[-2], lv.shape[-1]) for lv in levels]
        if tuple(shapes) != tuple(record.pyramid_shapes):
            raise AL.AlignError(
                f"{qptiff.name} pyramid changed since inventory: {shapes} != "
                f"{list(record.pyramid_shapes)}")
        lvl = TL.level_for(shapes, pixel_um, AL.ALIGN_TARGET_UM)
        lvl_um = pixel_um * (2 ** lvl)
        page = levels[lvl].pages[record.dapi_page]

        for w in probes:
            x0, y0, x1, y1 = w["bounds"]
            base = dict(islet_index=w["index"], cx_um=w["cx_um"], cy_um=w["cy_um"])
            grid, _px = TL.window_grid((x0, y0, x1, y1), lvl_um)
            if "member_indices" in w:
                idx = np.asarray(w["member_indices"], dtype=np.int64)
            else:
                inside = ((cen_if[:, 0] >= x0) & (cen_if[:, 0] <= x1)
                          & (cen_if[:, 1] >= y0) & (cen_if[:, 1] <= y1))
                idx = np.flatnonzero(inside)
            n_nuc = int(len(idx))
            if n_nuc == 0:
                offsets.append(AL.Offset(**base, n_nuclei=0,
                                         reason="no nuclei in this window"))
                continue
            warped = shapely.transform(arr[idx], lambda a: tf.apply(a))
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
    saved_transform = transform not in {"ecc", "valis"}
    donor_dir = out_dir / "valis" / f"donor_id={donor}"
    transform_file = donor_dir / f"{transform}.json" if saved_transform else None
    displacement_file = (donor_dir / f"{transform}_displacement.npz"
                         if saved_transform and tf.displacement is not None else None)
    out["alignment"] = {
        "mode": mode, "window_set": resolved_windows,
        "window_request": windows, "window_fold": window_fold,
        "transform_source": transform,
        "transform_file": None if transform_file is None else str(transform_file),
        "transform_sha256": None if transform_file is None else _sha256(transform_file),
        "displacement_file": (None if displacement_file is None
                              else str(displacement_file)),
        "displacement_sha256": (None if displacement_file is None
                                else _sha256(displacement_file)),
        "has_displacement": tf.displacement is not None,
        "affine": np.asarray(tf.affine).tolist(),
        "window_um": WINDOW_UM if resolved_windows in {"tissue", "occupied"} else None,
        "window_source": ("QuPath annotation" if resolved_windows in {"tissue", "islets"}
                          else "transformed Xenium nucleus occupancy; no IF intensity"),
        "qptiff": str(qptiff), "pyramid_level": int(lvl),
        "align_um_per_px": float(lvl_um),
        "certifies_mapping": {"pixel_um_full_res": pixel_um,
                              "source": "qptiff XML PixelSizeMicrons",
                              "dapi_page": int(record.dapi_page),
                              "dapi_channel_name": record.channel_names[record.dapi_page]},
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
    ap.add_argument("--windows", default="auto",
                    choices=["auto", "tissue", "islets", "occupied"],
                    help="auto uses annotation tissue windows when GeoJSON exists and otherwise "
                         "uses deterministic occupied 500-um windows")
    ap.add_argument("--transform", default="refined",
                    help="which registration to measure: the tissue-outline ECC fit, the raw "
                         "VALIS DAPI-to-DAPI fit, or an explicitly named donor-local transform "
                         "such as refined or intensity_refined")
    ap.add_argument("--window-fold", default="all", choices=["all", "even", "odd"],
                    help="measure all occupied windows, or one origin-anchored checkerboard "
                         "fold (used to keep intensity-refinement validation held out)")
    ap.add_argument("--inventory", type=Path, default=RI.REGISTRATION_INPUTS)
    ap.add_argument("--images-dir", type=Path, default=RI.IMAGES)
    ap.add_argument("--geojson-dir", type=Path, default=RI.GEOJSON_DIR)
    ap.add_argument("--xenium-paths", type=Path, default=RI.XENIUM_PATHS)
    ap.add_argument("--xenium-root", type=Path, default=RI.XENIUM_ROOT)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    records = RI.load_registration_inputs(
        args.inventory, args.xenium_paths, args.images_dir, args.geojson_dir,
        args.xenium_root, roi=args.roi)
    selected = RI.select_registration_inputs(records, args.donors, context="run_align")

    rows, worst = [], "PASS"
    for record in selected:
        RI.require_registration_files(record)
        rec = align_one(record, mode=args.mode, windows=args.windows,
                        transform=args.transform, window_fold=args.window_fold,
                        out_dir=args.out,
                        xenium_paths=args.xenium_paths)
        a = rec["alignment"]
        if args.transform not in {"ecc", "valis"}:
            gate = (args.out / "valis" / f"donor_id={record.donor_id}"
                    / f"alignment_{args.transform}.json")
            a["gate_file"] = str(gate)
            _json_write(gate, rec)
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
    suffix = "" if args.windows in {"auto", "tissue"} else f"_{args.windows}"
    if args.transform != "ecc":
        suffix += f"_{args.transform}"
    if args.window_fold != "all":
        suffix += f"_{args.window_fold}"
    p = args.out / f"alignment_{args.roi}{suffix}.json"
    _json_write(p, rows)
    print(f"-> {p}", flush=True)
    return {"PASS": 0, "INCONCLUSIVE": 3, "FAIL": 4}[worst]


if __name__ == "__main__":
    raise SystemExit(main())
