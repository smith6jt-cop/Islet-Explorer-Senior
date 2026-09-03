#!/usr/bin/env python
"""Write each platform's DAPI plane as a plain single-channel OME-TIFF, for VALIS.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/export_dapi.py --donors 6442

VALIS registers FILES in a directory. Handing it the raw inputs means handing it a 5-channel
qptiff and a 4-file Xenium morphology set, and hoping it picks the DAPI plane out of each; this
writes exactly the two planes we mean, with the pixel size recorded in the OME-XML so VALIS can
reconcile the two resolutions itself (IF 0.25 um/px, Xenium 0.2125 um/px at full res, and
different pyramid levels chosen for each).

The pixel size written here comes from the qptiff's exact XML ``PixelSizeMicrons`` value,
not its approximate TIFF resolution rational and not a cohort-wide constant. The strict
registration inventory records image/serial/ROI provenance only; raw channel names,
calibration, DAPI page, dtype, and pyramid geometry are read and validated from each image.

**Each side also records where its pixel (0, 0) actually sits in microns**, because the two
pyramids are built differently and the naive ``micron = pixel * um_per_px`` is wrong for both,
in different amounts. Measured on donor 6436, three 384-px tiles each, level 0 against level 2:

    Xenium morphology   decimation phase (0,0)  r = 0.9936   box mean  r = 0.9766
    IF qptiff (DAPI)    box mean                r = 0.99992  best decimation r = 0.9891

So the Xenium pyramid KEEPS level-0 pixel ``f*k`` and the IF pyramid AVERAGES ``f*k .. f*k+f-1``.
Their pixel centres are therefore at ``(f*k + 0.5) * base`` and ``(f*k + f/2) * base``, i.e.
``origin_um`` of ``0.5 * base`` and ``0.5 * f * base``. Left uncorrected that is a fixed
+0.39 um per axis on the warped result -- one twenty-fifth of a cell, but a constant one, and
the residual this lane is chasing is itself a constant.
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
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import registration_inputs as RI                                  # noqa: E402
import tiles as TL                                                # noqa: E402
from run_register import OUT                                      # noqa: E402

#: Work near 1 um/px. Nuclei are ~6 um so they stay several pixels across, and a whole section
#: is then a few hundred megapixels rather than several billion.
TARGET_UM = 1.0
#: Xenium morphology full-resolution scale, per ``rasterize.read_xenium_morphology``.
XENIUM_BASE_UM = 0.2125


def pyramid_factor(shapes, lvl: int) -> int:
    """Integer downsample factor of pyramid level ``lvl``, refused if it is not one.

    ``shapes[0][0] / shapes[lvl][0]`` is NOT the factor. The levels are floor-halved, so on
    donor 6436 that ratio is 68119/17029 = 4.000176 where the factor is exactly 4 -- a 4.4e-5
    scale error, which is a displacement growing to 0.6 um at the far corner of the section.
    The factor is an integer by construction; anything else is a pyramid this code has not
    seen and must not guess at.
    """
    if not 0 <= lvl < len(shapes):
        raise ValueError(f"pyramid level {lvl} is outside 0..{len(shapes) - 1}")
    # These pyramids are built by repeated two-fold reduction.  Inferring the factor from a
    # rounded dimension ratio is subtly wrong at deep levels when level 0 has an odd length
    # (for example 82080/320 == 256.5).  The level index is the exact source of truth.
    f = 2 ** lvl
    for axis in (0, 1):
        # floor-halving: level `lvl` may be up to 1 px short of shapes[0]/f, never more.
        expect = shapes[0][axis] // f
        if not 0 <= expect - shapes[lvl][axis] <= 1:
            raise ValueError(
                f"level {lvl} shape {shapes[lvl]} is not level 0 {shapes[0]} downsampled by "
                f"{f} on axis {axis} (expected {expect})")
    return f


def origin_um(base_um: float, factor: int, *, averaged: bool) -> float:
    """Micron position of exported pixel coordinate 0, given how the pyramid was built.

    ``averaged`` -- level-`lvl` pixel k is the mean of level-0 pixels ``f*k .. f*k+f-1``, whose
    centre is ``(f*k + f/2) * base``; the exported pixel size already carries ``f*base``, so the
    leftover is ``0.5 * f * base``. Otherwise level-`lvl` pixel k IS level-0 pixel ``f*k``,
    centred at ``(f*k + 0.5) * base``, leaving ``0.5 * base``.
    """
    return 0.5 * base_um * (factor if averaged else 1)


def _write_ome(path: Path, arr: np.ndarray, um_per_px: float) -> Path:
    import tifffile

    path.parent.mkdir(parents=True, exist_ok=True)
    # The channel is NAMED "DAPI" because VALIS's fluorescence path selects a plane by channel
    # name (`ChannelGetter`, `if_processing_kwargs={"channel": "dapi"}`); an unnamed channel
    # would leave it guessing on a file that has exactly one.
    tifffile.imwrite(
        str(path), arr, photometric="minisblack", ome=True,
        metadata={"axes": "YX", "Channel": {"Name": "DAPI"},
                  "PhysicalSizeX": float(um_per_px),
                  "PhysicalSizeY": float(um_per_px),
                  "PhysicalSizeXUnit": "µm", "PhysicalSizeYUnit": "µm"},
    )
    print(f"[export] {path}  {arr.shape} {arr.dtype}  {um_per_px:.5f} um/px", flush=True)
    return path


def check_xenium_origin(bundle: Path, level0_shape, base_um: float = XENIUM_BASE_UM) -> dict:
    """The repo assumes morphology pixel (0,0) is micron (0,0). Nothing has ever tested it.

    A wrong origin does not raise -- it puts every Xenium cell a fixed distance from its own
    nucleus, which is exactly the failure this whole exercise is trying to remove. So it is
    checked against the cell table's own extent before anything is fitted on it.
    """
    import pyarrow.parquet as pq

    c = pq.read_table(bundle / "cells.parquet",
                      columns=["x_centroid", "y_centroid"]).to_pandas()
    x1, y1 = float(c.x_centroid.max()), float(c.y_centroid.max())
    img_x_um, img_y_um = level0_shape[1] * base_um, level0_shape[0] * base_um
    ok = (img_x_um >= x1) and (img_y_um >= y1) and (img_x_um < 1.5 * x1) and (img_y_um < 1.5 * y1)
    return {"image_extent_um": [img_x_um, img_y_um], "cell_extent_um": [x1, y1],
            "level0_shape": list(level0_shape), "base_um_per_px": base_um,
            "brackets_cell_table": bool(ok)}


def export_one(record: RI.RegistrationInput, *, out_root: Path = OUT / "valis") -> dict:
    import tifffile

    donor, serial, roi = record.donor_id, record.serial, record.roi
    qptiff, pixel_um = record.qptiff, record.pixel_um

    dst = out_root / f"donor_id={donor}"
    slides = dst / "slides"
    rec = {"donor_id": donor, "serial": serial, "roi": roi, "dir": str(dst)}

    # --- IF DAPI: one page of the qptiff, never the 5-channel decode ---------------------
    with tifffile.TiffFile(str(qptiff)) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series]))
        shapes = [(lv.shape[-2], lv.shape[-1]) for lv in levels]
        lvl = TL.level_for(shapes, pixel_um, TARGET_UM)
        f_if = pyramid_factor(shapes, lvl)
        if_um = pixel_um * f_if
        if_origin = origin_um(pixel_um, f_if, averaged=True)      # measured: box mean, r=0.99992
        page = levels[lvl].pages[record.dapi_page]
        if_img = page.asarray()
    rec["if"] = {"qptiff": str(qptiff), "level": int(lvl), "um_per_px": if_um,
                 "factor": int(f_if), "origin_um": if_origin, "pyramid": "box_mean",
                 "shape": list(if_img.shape), "dtype": str(if_img.dtype)}
    _write_ome(slides / "02_if_dapi.ome.tiff", if_img, if_um)
    del if_img

    # --- Xenium DAPI: morphology_focus channel 0 -----------------------------------------
    bundle = record.bundle
    focus = bundle / "morphology_focus" / "morphology_focus_0000.ome.tif"
    if not focus.exists():
        raise FileNotFoundError(f"{focus} is absent -- transfer it before exporting")
    with tifffile.TiffFile(str(focus)) as tif:
        series = tif.series[0]
        levels = list(getattr(series, "levels", [series]))
        shapes = [(lv.shape[-2], lv.shape[-1]) for lv in levels]
        rec["xenium_origin_check"] = check_xenium_origin(bundle, shapes[0])
        lvl = TL.level_for(shapes, XENIUM_BASE_UM, TARGET_UM)
        f_xen = pyramid_factor(shapes, lvl)
        xen_um = XENIUM_BASE_UM * f_xen
        xen_origin = origin_um(XENIUM_BASE_UM, f_xen, averaged=False)  # measured: decimated (0,0)
        # ONE page, not `levels[lvl].asarray()`. The OME series declares four channels across
        # four files and only ``_0000`` is on this machine, so the series form materialises
        # three zeroed planes to keep one, and silently depends on DAPI being channel 0.
        xen_img = levels[lvl].pages[0].asarray()
        if xen_img.ndim == 3:
            xen_img = xen_img[0]
    import pyarrow.parquet as pq
    cells = pq.read_table(bundle / "cells.parquet",
                          columns=["cell_id", "x_centroid", "y_centroid"]).to_pandas()
    np.save(dst / "xenium_centroids_um.npy",
            cells[["x_centroid", "y_centroid"]].to_numpy(np.float64))
    (dst / "xenium_cell_ids.npy").parent.mkdir(parents=True, exist_ok=True)
    np.save(dst / "xenium_cell_ids.npy", cells["cell_id"].astype(str).to_numpy())
    rec["n_centroids"] = int(len(cells))
    rec["xenium"] = {"file": str(focus), "level": int(lvl), "um_per_px": xen_um,
                     "factor": int(f_xen), "origin_um": xen_origin, "pyramid": "decimated",
                     "shape": list(xen_img.shape), "dtype": str(xen_img.dtype),
                     "levels": [list(s) for s in shapes]}
    _write_ome(slides / "01_xenium_dapi.ome.tiff", xen_img, xen_um)

    # The VALIS venv reads only this and the .npy arrays -- it imports nothing from the repo.
    (dst / "export.json").write_text(json.dumps(
        {"donor_id": donor, "serial": serial, "roi": roi,
         "xenium_um_per_px": xen_um, "if_um_per_px": if_um,
         "xenium_origin_um": xen_origin, "if_origin_um": if_origin,
         "xenium_bundle": str(bundle), "qptiff": str(qptiff),
         "if_pixel_um": pixel_um, "if_dapi_page": record.dapi_page,
         "if_channel_names": list(record.channel_names),
         "if_pyramid_shapes": [list(shape) for shape in record.pyramid_shapes],
         "pixel_size_source": "qptiff XML PixelSizeMicrons",
         "provenance_source": record.provenance_source,
         "geojson": None if record.geojson is None else str(record.geojson),
         "if": rec["if"], "xenium": rec["xenium"],
         "xenium_origin_check": rec["xenium_origin_check"]}, indent=2))
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
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
    rows = []
    for record in RI.select_registration_inputs(records, args.donors,
                                                context="export_dapi"):
        RI.require_registration_files(record)
        donor = record.donor_id
        rec = export_one(record, out_root=args.out / "valis")
        c = rec["xenium_origin_check"]
        print(f"[{donor}] xenium origin check: image {np.round(c['image_extent_um'],0).tolist()} um "
              f"vs cells {np.round(c['cell_extent_um'],0).tolist()} um -> "
              f"{'OK' if c['brackets_cell_table'] else 'MISMATCH'}", flush=True)
        rows.append(rec)
    if not rows:
        print("no donor matched", file=sys.stderr)
        return 1
    p = args.out / "valis" / f"export_{args.roi}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(rows, indent=2))
    print(f"-> {p}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
