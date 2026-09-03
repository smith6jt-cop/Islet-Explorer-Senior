#!/usr/bin/env python
"""The registered post-Xenium IF, plus Xenium's segmentations coloured by INS/GCG positivity.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_explorer_export.py --donors 6436 6442

Per donor:

    donor_<id>_<serial>/
      images/<serial>_IF.ome.tif                       DAPI + INS + GCG, native 0.25 um/px grid
      images/<serial>_IF_imagealignment.csv            the 3x3 for that image (IF px -> Xenium px)
      images/<serial>_IF_similarity_imagealignment.csv the nearest similarity, a diagnostic
      images/<serial>_IF_keypoints.csv                 20 exact nucleus pairs, Explorer's other input
      images/<serial>_IF_on_xenium_grid.ome.tif        the same IF resampled onto the Xenium
                                                       morphology grid, for "Yes, skip alignment"
      cell_groups/hormone_positivity.csv               cell_id, group, color
      README.md, manifest.json                         how to import it, and what the numbers are

Explorer 4.x refused the previous matrix file with "Failed to load file." 10x documents no
validation rule, so every way the file differed from Explorer's own export was removed (bytes),
the two other documented inputs were added (keypoints, which Explorer fits itself; a similarity,
which is the only shape Explorer ever exports), and the image is also shipped already on the
Xenium grid, the one route whose placement depends on nothing but documented behaviour.

Positivity is the mean 8-bit value of the channel inside **that cell's own Xenium polygon**,
against ``--threshold`` (default 20). One number, applied identically to every donor, recorded
in the manifest and in the file name of nothing -- change the flag and re-run to move it. It is
not fitted, not per-donor, and not tuned to a frequency; the resulting counts are outputs.

The evaluation happens in Xenium Explorer. This writes what is needed to do that.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import explorer as EX                                            # noqa: E402
import geo                                                       # noqa: E402
import hormone as HZ                                             # noqa: E402
import protein_rna as PR                                         # noqa: E402
import xenio                                                     # noqa: E402
from run_hormone import IMAGES, SECTIONS                         # noqa: E402
from run_register import OUT, XENIUM_PATHS, discover, donor_of   # noqa: E402

BUNDLE_ROOT = OUT / "explorer_bundle"

#: Channel roles written into the image. `restain_sections.csv` records Opal 570 FAILED and
#: Sample AF as autofluorescence; neither is exported, and the manifest records why.
EXPORT_ROLES = ("NUCLEUS", "GCG", "INS")

#: Mean 8-bit counts inside a cell's own polygon, at or above which that cell is called positive.
DEFAULT_THRESHOLD = 20.0

GROUP_COLORS = {"INS_positive": "#20c997", "GCG_positive": "#f76707",
                "INS_and_GCG_positive": "#ae3ec9", "Negative": "#adb5bd",
                "Unavailable": "#495057"}

#: The one convention the lane knows it cannot settle from any document; recorded in every
#: manifest rather than assumed away (it was, between 6435fe3 and this).
IF_PIXEL_CONVENTION_NOTE = (
    "no 10x document states whether Explorer indexes an imported image by pixel corner or pixel "
    "centre. This matrix assumes corner, matching 10x's own micron-to-pixel rule for the Xenium "
    "side. If Explorer assumes centre, every cell is off by if_origin_um, in both axes.")

#: The round trip Xenium um -> IF um -> IF px -> (M) Xenium px -> Xenium um over every cell; a
#: matrix that is a re-expression of the affine reproduces it to floating-point noise.
ROUND_TRIP_TOLERANCE_UM = 1e-6


class ExportError(RuntimeError):
    """A bundle cannot be assembled from what is on disk -- never partially, never guessed."""


def _reset(dest: Path, *, keep_images: bool) -> list[str]:
    """Clear the bundle so it is rebuilt whole, sparing the OME-TIFFs already in it if asked.

    A bundle is rewritten rather than updated, because a half-refreshed one -- new matrix,
    stale groups -- is the kind of thing that gets sent to a colleague without anyone noticing.
    ``keep_images`` is the exception, and a reset that deleted the multi-GB images would
    destroy exactly what that flag exists to preserve. Every CSV is deleted regardless, so no
    stale alignment file can ride alongside a new image. The stash is restored in ``finally``:
    a crash between the rename out and the rename back must not strand the images.
    """
    dest = Path(dest)
    if not dest.exists():
        (dest / "images").mkdir(parents=True)
        return []
    keep = sorted((dest / "images").glob("*.ome.tif")) if keep_images else []
    stash = dest.parent / f".{dest.name}.keep"
    if stash.exists():
        raise ExportError(f"{stash} is left over from an interrupted run; move or delete it "
                          "before rebuilding, so nothing from it is silently resurrected")
    try:
        if keep:
            stash.mkdir(parents=True)
            for p in keep:
                p.rename(stash / p.name)
        shutil.rmtree(dest)
        (dest / "images").mkdir(parents=True)
    finally:
        if keep and stash.exists():
            (dest / "images").mkdir(parents=True, exist_ok=True)
            for p in sorted(stash.iterdir()):
                p.rename(dest / "images" / p.name)
            stash.rmdir()
    return [p.name for p in keep]


def _describe_ome_tiff(path: Path) -> dict | None:
    """A kept image's manifest record, read from the file itself in the writer's own shape.

    A skipped rewrite must describe the file that is actually alongside, so the record comes
    from its bytes -- pixel size from the OME-XML, levels from the pyramid -- never from a
    previous manifest, which a failed run may have removed. ``None`` if the file is not there.
    """
    import re

    import tifffile

    path = Path(path)
    if not path.exists():
        return None
    with tifffile.TiffFile(str(path)) as t:
        ome = t.ome_metadata or ""
        m = re.search(r'PhysicalSizeX="([^"]+)"', ome)
        if m is None:
            raise ExportError(f"{path} carries no PhysicalSizeX; a kept image cannot be described")
        um = float(m.group(1))
        s = t.series[0]
        h0 = int(s.levels[0].shape[-2])
        levels = [{"level": L, "height": int(lv.shape[-2]), "width": int(lv.shape[-1]),
                   "factor": int(round(h0 / lv.shape[-2])),
                   "um_per_px": um * int(round(h0 / lv.shape[-2]))}
                  for L, lv in enumerate(s.levels)]
        pg = t.pages[0]
        rec = {"path": str(path), "dtype": str(s.dtype), "tile": int(pg.tilewidth),
               "compression": str(pg.compression.name).lower(),
               "channels": re.findall(r'Channel[^>]*Name="([^"]+)"', ome),
               "levels": levels, "bytes": int(path.stat().st_size),
               "skipped": True, "preserved_from_a_previous_run": [path.name],
               "described_from": "the file's own OME-XML and pyramid"}
    return rec


def _channel_plan(qptiff: Path, row) -> tuple[list[int], list[str], list[dict]]:
    """Which qptiff pages are exported, under which names, and why the rest are not.

    Indices come from the qptiff's own ``<Name>`` XML, never from page order: the cohort mixes
    channel orders between batches, and a hardcoded order hands over someone else's stain.
    """
    names = HZ.channel_names(qptiff)
    idx, out_names, dropped = [], [], []
    for col, role in ((k, str(v).strip()) for k, v in row.items()
                      if isinstance(k, str) and k.startswith("channel_")):
        page = HZ.channel_index_for_role(qptiff, row, role)
        if page is None:
            continue
        if role.upper() in {r.upper() for r in EXPORT_ROLES}:
            idx.append(page)
            out_names.append("DAPI" if role.upper() == "NUCLEUS" else role)
        else:
            dropped.append({"qptiff_channel": names[page], "recorded_role": role,
                            "reason": ("recorded FAILED; never to be used"
                                       if role.upper() == "FAILED"
                                       else "recorded AF; disregarded per the maintainer")})
    order = np.argsort(idx)
    return ([idx[i] for i in order], [out_names[i] for i in order], dropped)


def per_cell_values(qptiff: Path, bundle: Path, cell_ids, transform, row,
                    *, um_per_px: float) -> dict:
    """Mean INS and GCG inside each Xenium cell's own polygon, in the IF frame.

    Polygon-exact rather than a disc at the centroid: cells are ~9 um apart and a disc
    straddles the neighbour whose hormone this is not.
    """
    import tifffile
    import zarr

    from phenocycler.integration.rasterize import GridSpec
    from rasterize_labels import rasterize_labels

    base_um = float(row.pixel_size_um)
    lvl, factor = PR.pyramid_level_for(qptiff, um_per_px, base_um)
    with tifffile.TiffFile(str(qptiff)) as tif:
        h = int(tif.series[0].levels[lvl].shape[-2])
        w = int(tif.series[0].levels[lvl].shape[-1])
    # The grid IS the image: same origin, same pixel size, same extent, so a label index and a
    # pixel index are the same thing and no second alignment can creep in between them.
    grid = GridSpec(origin_x_um=0.5 * um_per_px, origin_y_um=0.5 * um_per_px,
                    pixel_um=um_per_px, height=h, width=w)

    geoms = PR.polygons_per_cell(bundle, "cell", cell_ids)
    keep = np.array([g is not None for g in geoms])
    idx = np.nonzero(keep)[0]
    lab = rasterize_labels(PR.warp_polygons([geoms[i] for i in idx], transform), grid,
                           order="input")

    out = {}
    with tifffile.TiffFile(str(qptiff)) as tif:
        store = tif.series[0].aszarr(level=lvl)
        try:
            z = zarr.open(store, mode="r")
            for marker in ("INS", "GCG"):
                page = HZ.channel_index(qptiff, marker, row)
                img = np.asarray(z[page]).astype(np.float32)
                sub = PR.mean_in_labels(img, lab.labels, len(idx))
                del img
                v = np.full(len(cell_ids), np.nan)
                v[idx] = sub
                out[marker] = v
        finally:
            store.close()
    return {"values": out, "n_without_polygon": int((~keep).sum()),
            "sampled_um_per_px": float(um_per_px), "pyramid_level": lvl,
            "pyramid_factor": int(factor)}


def _footprint(M: np.ndarray, if_hw: tuple[int, int]) -> dict:
    """Where the imported image's corners land in Xenium pixels -- the picture of 'unaligned'."""
    h, w = int(if_hw[0]), int(if_hw[1])
    # M's IF input is centre-indexed: pixel 0 is centred at 0, so the image's edge is at -0.5
    # and the far edge at w - 0.5. Half an IF pixel = 0.125 um; recorded rather than rounded.
    corners = np.array([[-0.5, -0.5, 1], [w - 0.5, -0.5, 1], [-0.5, h - 0.5, 1],
                        [w - 0.5, h - 0.5, 1]], dtype=np.float64)
    xy = (corners @ np.asarray(M).T)[:, :2]
    return {"if_corners_xenium_px": xy.tolist(),
            "x_range": [float(xy[:, 0].min()), float(xy[:, 0].max())],
            "y_range": [float(xy[:, 1].min()), float(xy[:, 1].max())],
            "convention": "IF image edges at -0.5 and w-0.5 in the matrix's centre-indexed "
                          "input frame; Xenium px corner-origin"}


def _warp_to_xenium_grid(qptiff: Path, pages: list[int], names: list[str], M0: np.ndarray,
                         grid: dict, out_path: Path, *, compression: str, workers: int) -> dict:
    """The selected qptiff channels resampled onto the morphology grid, then written.

    Always from qptiff level 0 through the level-0 matrix, whatever ``--level`` ships natively.
    One level-0 plane at a time in RAM (2.9-3.3 GB) into a preallocated (C, H, W) uint8 stack
    (7.6-10.2 GB), then the pyramid (its level-1 uint16 transient is another half of that),
    then the file. Peak about 16 GB (6436) / 22 GB (6442).
    """
    import time

    import cv2
    import tifffile

    H, W = int(grid["height"]), int(grid["width"])
    t0 = time.time()
    level0 = np.zeros((len(pages), H, W), dtype=np.uint8)
    with tifffile.TiffFile(str(qptiff)) as tif:
        lvl0 = tif.series[0].levels[0]
        for c, page in enumerate(pages):
            plane = lvl0.pages[page].asarray(maxworkers=int(workers))
            if plane.dtype != np.uint8:
                raise ExportError(f"qptiff page {page} is {plane.dtype}, not uint8")
            EX.warp_plane_to_xenium_grid(plane, M0, H, W, out=level0[c])
            del plane
            print(f"[explorer]   warped {names[c]} onto the {H} x {W} Xenium grid "
                  f"({time.time() - t0:.0f} s)", flush=True)
    # How much of the morphology grid has no IF under it: sampled on a 64-px lattice of pixel
    # centres, pushed through M0^-1 against the plane's extent.
    ys, xs = np.mgrid[32:H:64, 32:W:64]
    pts = np.column_stack([xs.ravel() + 0.5, ys.ravel() + 0.5, np.ones(xs.size)])
    src = pts @ np.linalg.inv(M0).T
    with tifffile.TiffFile(str(qptiff)) as tif:
        h_if, w_if = tif.series[0].levels[0].shape[-2:]
    outside = ((src[:, 0] < -0.5) | (src[:, 0] > w_if - 0.5)
               | (src[:, 1] < -0.5) | (src[:, 1] > h_if - 0.5))
    levels = EX.xenium_pyramid(level0, int(grid["n_levels"]))
    rec = EX.write_warped_ome_tiff(levels, out_path, names, compression=compression,
                                   workers=workers)
    del levels, level0
    rec["seconds"] = round(time.time() - t0, 1)
    rec["outside_if_fraction"] = float(outside.mean())
    rec["opencv_version"] = cv2.__version__
    for lv in rec["levels"]:
        print(f"[explorer]   level {lv['level']}: {lv['height']} x {lv['width']} x {len(names)}",
              flush=True)
    return rec


def positivity_groups(ins: np.ndarray, gcg: np.ndarray, threshold: float):
    """``INS_positive`` / ``GCG_positive`` / both / ``Negative`` / ``Unavailable``.

    A cell whose polygon fell outside the imaged area is ``Unavailable``, never ``Negative``:
    nothing was measured for it, and calling that negative is the one error that cannot be
    seen in the viewer.
    """
    ins = np.asarray(ins, dtype=float)
    gcg = np.asarray(gcg, dtype=float)
    ok = np.isfinite(ins) & np.isfinite(gcg)
    i_pos, g_pos = ok & (ins >= threshold), ok & (gcg >= threshold)
    labels = np.full(len(ins), "Unavailable", dtype=object)
    labels[ok] = "Negative"
    labels[i_pos] = "INS_positive"
    labels[g_pos] = "GCG_positive"
    labels[i_pos & g_pos] = "INS_and_GCG_positive"
    return labels


def export_one(serial: str, sections, *, roi: str = "panc", out_dir: Path = OUT,
               bundle_root: Path = BUNDLE_ROOT, level: int = 0, compression: str = "zlib",
               skip_image: bool = False, skip_warped: bool = False, workers: int = 8,
               threshold: float = DEFAULT_THRESHOLD,
               sample_um_per_px: float = PR.PROTEIN_UM_PER_PX) -> dict:
    from run_align import _refined_transform

    donor = donor_of(serial, XENIUM_PATHS)
    row = sections[sections.xenium_serial == serial]
    if row.empty:
        raise ExportError(f"serial {serial} is not recorded in {SECTIONS}")
    row = row.iloc[0]
    qptiff = IMAGES / row.image_file
    if not qptiff.exists():
        raise ExportError(f"{qptiff} is absent")
    pixel_um = float(row.pixel_size_um)
    if pixel_um != geo.PIXEL_UM:
        raise ExportError(f"{row.image_file} records pixel_size_um={pixel_um} but this lane's "
                          f"geometry is in geo.PIXEL_UM={geo.PIXEL_UM}")

    transform = _refined_transform(donor, roi, out_dir)
    if transform.displacement is not None:
        raise ExportError("the chosen transform carries a displacement field, and a 3x3 "
                          "alignment matrix cannot represent one; this refuses rather than "
                          "silently truncating it to the affine part")
    if transform.is_mirrored:
        raise ExportError("the chosen transform is mirrored; the similarity diagnostic is fitted "
                          "without reflection and would be wrong")

    bundle = xenio._bundle_dir(donor, XENIUM_PATHS, roi)
    import pandas as pd
    cells = pd.read_parquet(bundle / "cells.parquet",
                            columns=["cell_id", "x_centroid", "y_centroid"])
    cell_ids = cells.cell_id.astype(str).to_numpy()
    xy_um = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=np.float64)
    grid = EX.morphology_grid(bundle)
    H, W = int(grid["height"]), int(grid["width"])

    shapes = EX.qptiff_level_shapes(qptiff)
    factors = EX.integer_factors(shapes)
    if not 0 <= level < len(shapes):
        raise ExportError(f"level {level} is outside this qptiff's {len(shapes)} levels")
    f = factors[level]
    if_um_per_px = pixel_um * f
    # Measured, not assumed: export_dapi.py established this pyramid is a box mean, so the
    # exported pixel's centre sits half its own width in.
    if_origin_um = 0.5 * pixel_um * f

    M = EX.alignment_matrix(transform.affine, if_um_per_px=if_um_per_px,
                            if_origin_um=if_origin_um)
    M_written = EX.rounded_alignment_matrix(M)
    round_trip = EX.matrix_round_trip_um(M, transform.affine, if_um_per_px=if_um_per_px,
                                         if_origin_um=if_origin_um, xy_xenium_um=xy_um)
    if round_trip > ROUND_TRIP_TOLERANCE_UM:
        raise ExportError(f"the alignment matrix does not reproduce the affine: round trip "
                          f"{round_trip:.3e} um over {len(xy_um)} cells")
    round_trip_written = EX.matrix_round_trip_um(M_written, transform.affine,
                                                 if_um_per_px=if_um_per_px,
                                                 if_origin_um=if_origin_um, xy_xenium_um=xy_um)
    # Everything derived for Explorer is derived from the matrix AS WRITTEN, so the shipped
    # files agree with each other exactly; the 9-digit rounding is worth < 3e-4 px anywhere.
    Minv = np.linalg.inv(M_written)

    native_name = f"{serial}_IF.ome.tif"
    warped_name = f"{serial}_IF_on_xenium_grid.ome.tif"
    dest = Path(bundle_root) / f"donor_{donor}_{serial}"

    # --- everything that can refuse, BEFORE the bundle is cleared --------------------------- #
    # None of this needs the destination; a refusal here leaves the previous bundle whole.
    pages, ch_names, dropped = _channel_plan(qptiff, row)
    kp_rec = EX.keypoint_pairs(xy_um, M_written, morph_hw=(H, W), if_hw=shapes[level],
                               cell_ids=cell_ids)
    kp = kp_rec["keypoints"]
    S = EX.similarity_from_keypoints(kp)
    consistency = EX.similarity_consistency(S, transform.affine, if_um_per_px=if_um_per_px)
    if not consistency["consistent"]:
        raise ExportError(f"the similarity through the keypoints disagrees with the transform: "
                          f"{consistency}")
    ones = np.ones(4)
    fov_if = (np.column_stack([[0, W, 0, W], [0, 0, H, H], ones]) @ Minv.T)[:, :2]
    cells_if = (np.column_stack([xy_um / EX.XENIUM_UM_PER_PX, np.ones(len(xy_um))]) @ Minv.T)[:, :2]
    dev_kp = EX.similarity_deviation_um(M_written, S, kp[:, 2:4])
    dev_fov = EX.similarity_deviation_um(M_written, S, fov_if)
    dev_cells = EX.similarity_deviation_um(M_written, S, cells_if)
    M0 = EX.alignment_matrix(transform.affine, if_um_per_px=pixel_um,
                             if_origin_um=0.5 * pixel_um)
    image_rec = warped_rec = None
    if skip_image:
        image_rec = _describe_ome_tiff(dest / "images" / native_name)
        if image_rec is None:
            raise ExportError(f"--skip-image, but {dest / 'images' / native_name} is not there")
        if image_rec["channels"] != ch_names:
            raise ExportError(f"the kept native image carries channels {image_rec['channels']}, "
                              f"not the planned {ch_names}; re-export it without --skip-image")
    if skip_warped:
        warped_rec = _describe_ome_tiff(dest / "images" / warped_name)
        if warped_rec is None:
            raise ExportError(f"--skip-warped, but {dest / 'images' / warped_name} is not there")
        if warped_rec["channels"] != ch_names:
            raise ExportError(f"the kept warped image carries channels {warped_rec['channels']}, "
                              f"not the planned {ch_names}; re-export it without --skip-warped")

    # --- clear, then write ------------------------------------------------------------------ #
    _reset(dest, keep_images=skip_image or skip_warped)
    (dest / "cell_groups").mkdir(parents=True)

    align_csv = EX.write_alignment_csv(dest / "images" / f"{serial}_IF_imagealignment.csv", M)
    print(f"[explorer] {donor}: alignment matrix -> {align_csv.name}  "
          f"(round trip {round_trip:.2e} um, as written {round_trip_written:.2e} um)", flush=True)
    kp_csv = EX.write_keypoints_csv(dest / "images" / f"{serial}_IF_keypoints.csv", kp)
    sim_csv = EX.write_alignment_csv(
        dest / "images" / f"{serial}_IF_similarity_imagealignment.csv", S)
    print(f"[explorer] {donor}: {len(kp)} keypoints -> {kp_csv.name}; similarity -> "
          f"{sim_csv.name} (deviation over cells p50 {np.median(dev_cells):.2f} um, "
          f"max {dev_cells.max():.2f} um; FOV corners max {dev_fov.max():.2f} um)", flush=True)

    if not skip_image:
        image_rec = EX.write_explorer_ome_tiff(
            qptiff, dest / "images" / native_name, pages, ch_names,
            base_um_per_px=if_um_per_px, compression=compression,
            max_levels=len(shapes) - level, workers=workers)
    image_rec["dropped_channels"] = dropped

    if not skip_warped:
        warped_rec = _warp_to_xenium_grid(qptiff, pages, ch_names, M0, grid,
                                          dest / "images" / warped_name,
                                          compression=compression, workers=workers)
    if warped_rec.get("channels") != ch_names:
        raise ExportError("the warped image's channels differ from the native image's")

    # --- positivity ----------------------------------------------------------------------- #
    sampled = per_cell_values(qptiff, bundle, cell_ids, transform, row,
                              um_per_px=sample_um_per_px)
    ins, gcg = sampled["values"]["INS"], sampled["values"]["GCG"]
    labels = positivity_groups(ins, gcg, threshold)
    groups = EX.write_cell_groups(dest / "cell_groups" / "hormone_positivity.csv",
                                  cell_ids, labels, GROUP_COLORS)
    print(f"[explorer] {donor}: {groups['groups']}", flush=True)

    if_h0, if_w0 = shapes[0]
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "donor_id": donor, "xenium_serial": serial, "roi": roi,
        "xenium_bundle": str(bundle), "qptiff": str(qptiff),
        "transform": {"source": "refined", **transform.to_dict()},
        "morphology": {
            **grid,
            "if_level0_shape": [int(if_h0), int(if_w0)],
            "if_over_morphology_pixel_ratio": float(if_h0 * if_w0) / float(H * W),
        },
        "alignment": {
            "file": str(align_csv.relative_to(dest)),
            "matrix": np.asarray(M).tolist(),
            "matrix_as_written": M_written.tolist(),
            "maps": "exported IF image pixels -> Xenium morphology image pixels",
            "applies_to": f"images/{native_name}",
            "never_apply_to": f"images/{warped_name} (it is already on the Xenium grid; a "
                              "matrix would apply the transform twice)",
            "xenium_um_per_px": EX.XENIUM_UM_PER_PX,
            "if_um_per_px": if_um_per_px, "if_origin_um": if_origin_um,
            "if_pyramid_level": int(level), "if_pyramid_factor": int(f),
            "if_pyramid_rule": "box_mean (measured, export_dapi.py)",
            "round_trip_max_um": float(round_trip),
            "round_trip_as_written_max_um": float(round_trip_written),
            "round_trip_n_cells": int(len(xy_um)),
            "csv_format": EX.alignment_csv_format(),
            "structure": EX.matrix_anisotropy(M),
            "footprint_xenium_px": _footprint(M, shapes[level]),
            "translation_sign": [int(np.sign(M[0, 2])), int(np.sign(M[1, 2]))],
            "unresolved_convention_um": if_origin_um,
            "if_pixel_convention_note": IF_PIXEL_CONVENTION_NOTE,
        },
        "keypoints": {
            "file": str(kp_csv.relative_to(dest)),
            "header": list(EX.KEYPOINTS_HEADER),
            "n": int(len(kp)), "grid": kp_rec["grid"],
            "placement": "one Xenium cell centroid nearest each node of a grid over the cell "
                         "extent, so every keypoint sits on a nucleus visible in both images",
            "cell_ids": kp_rec["cell_ids"],
            "fixed_frame": "Xenium morphology image pixels (centroid um / 0.2125)",
            "alignment_frame": f"images/{native_name} pixels, exactly M^-1 @ fixed with M as "
                               "written in the shipped CSV",
            "residual_px_through_shipped_matrix": EX.keypoints_residual_px(kp, M_written),
            "residual_px_through_full_precision_matrix": EX.keypoints_residual_px(kp, M),
            "csv_format": EX.keypoints_csv_format(),
        },
        "similarity_diagnostic": {
            "file": str(sim_csv.relative_to(dest)),
            "matrix": np.asarray(S).tolist(),
            "matrix_as_written": EX.rounded_alignment_matrix(S).tolist(),
            "source": "least-squares similarity (Umeyama; Explorer documents Horn 1987) through "
                      "the shipped keypoints, canonicalised so a == d and b == -c exactly",
            "status": "diagnostic, second choice: a projection of the shipped matrix at chosen "
                      "points with the cost recorded below, not the transform and not a fit to "
                      "any data",
            "scale": float(np.sqrt(abs(np.linalg.det(S[:2, :2])))),
            "rotation_deg": float(np.degrees(np.arctan2(S[1, 0], S[0, 0]))),
            "is_exact_similarity": EX.is_exact_similarity(S),
            "consistency_with_transform": consistency,
            "deviation_from_affine_um": {
                "at_keypoints_max": float(dev_kp.max()),
                "at_fov_corners_max": float(dev_fov.max()),
                "over_cells": {"p50": float(np.percentile(dev_cells, 50)),
                               "p90": float(np.percentile(dev_cells, 90)),
                               "max": float(dev_cells.max()), "n": int(len(dev_cells))},
            },
        },
        "image": image_rec,
        "image_on_xenium_grid": {
            "file": f"images/{warped_name}",
            "import_route": "Images -> Add image -> 'Yes, skip alignment' (no alignment file)",
            "maps": "identity: imported pixel (k, l) is Xenium morphology pixel (k, l); "
                    "1 px = 0.2125 um",
            "width": W, "height": H, "n_levels": int(grid["n_levels"]),
            "level_shapes": grid["level_shapes"],
            "pyramid": "2x2 box mean of the warped level 0, floor-halved, round half up; level-L "
                       "pixel k covers level-0 pixels 2^L k .. 2^L (k+1) - 1, centre at "
                       "0.2125 * 2^L * (k + 0.5) um. 10x's own morphology pyramid is decimated, "
                       "which differs by under half a level-L pixel at coarse levels only",
            "resampling": "output pixel (k, l) = bilinear sample (cv2 INTER_LINEAR, 1/32-px "
                          "fixed point) of the native IF level 0 (0.25 um/px, pixel-centre "
                          "indexed) at M^-1 (k + 0.5, l + 0.5); 8-bit, no intensity rescaling; "
                          "outside the IF image = 0",
            "matrix_used": M0.tolist(),
            "same_transform_as": "alignment.matrix (identical for --level 0)",
            "tiff": warped_rec,
            "for_viewing_only": "positivity groups were measured on the native raster at "
                                f"{sample_um_per_px} um/px, not on this resampled one; hormone "
                                "values here are interpolated",
        },
        "positivity": {
            "file": "cell_groups/hormone_positivity.csv",
            "measure": "mean 8-bit channel value inside the cell's own Xenium cell polygon",
            "threshold_counts": float(threshold),
            "sampled_um_per_px": sampled["sampled_um_per_px"],
            "n_cells": int(len(cell_ids)),
            "n_without_polygon": sampled["n_without_polygon"],
            "counts": groups["groups"],
            "channels": {m: HZ.channel_names(qptiff)[HZ.channel_index(qptiff, m, row)]
                         for m in ("INS", "GCG")},
        },
        "uses_qupath_annotations": False,
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False))
    (dest / "README.md").write_text(_readme(manifest))
    return manifest


def _readme(m: dict) -> str:
    a, p, k, s, w = (m["alignment"], m["positivity"], m["keypoints"],
                     m["similarity_diagnostic"], m["image_on_xenium_grid"])
    serial = m["xenium_serial"]
    native = Path(a["applies_to"]).name
    warped = Path(w["file"]).name
    matrix_csv, sim_csv, kp_csv = (Path(a["file"]).name, Path(s["file"]).name,
                                   Path(k["file"]).name)
    bundle_name = Path(m["xenium_bundle"]).name
    W, H = m["morphology"]["width"], m["morphology"]["height"]
    tx, ty = m["transform"]["translation_um"]
    dev = s["deviation_from_affine_um"]
    ratio = m["morphology"]["if_over_morphology_pixel_ratio"]
    fp = a["footprint_xenium_px"]
    counts = "".join(f"| `{kk}` | {v:,} |\n" for kk, v in sorted(p["counts"].items()))
    return f"""# Donor {m['donor_id']} — post-Xenium IF registered to Xenium {serial}

Generated {m['generated_utc']}.

## Before you start

1. Note Xenium Explorer's version (**Help → About**). Every report below needs it.
2. **Remove any IF image layer already added to this dataset** from an earlier attempt. A layer
   placed with "Yes, skip alignment" on the native image sits about {abs(tx):.0f} × {abs(ty):.0f} µm
   off and will confuse every retry.
3. Open the Xenium output bundle **`{bundle_name}`** — the folder name contains `__Panc__`. A
   `__pLN__` bundle with the same serial exists and is the wrong one.
4. Explorer should show the morphology image as **{W} × {H} px** at 0.2125 µm/px
   (`experiment.xenium` → `pixel_size`). The IF is {ratio:.2f}× the morphology's pixel count,
   so it is not "much larger" than the DAPI image.

## Attempts, in this order

Keep only the first attempt that works; one image layer per attempt.

| # | Images → Add image → file | when asked about alignment | then upload |
|---|---|---|---|
| 1 | `images/{native}` | **"Yes, use an existing alignment file"** | `images/{matrix_csv}` |
| 2 | `images/{native}` | **"Yes, use an existing alignment file"** | `images/{sim_csv}` |
| 3 | `images/{native}` | **"Yes, use an existing alignment file"** | `images/{kp_csv}` |
| 4 | `images/{warped}` | **"Yes, skip alignment"** | nothing |

* Attempt 1 is the registered transform, an affine. Explorer 4.x refused the previous version of
  this file ("Failed to load file."); it is now byte-for-byte in the layout Explorer itself
  exports (CRLF, no trailing newline, last row `0,0,1`, at most nine significant digits).
* Attempt 2 is a **diagnostic**: the nearest similarity to the affine, the only shape Explorer's
  own aligner produces. It sits within {dev['over_cells']['p50']:.2f} µm (median) and
  {dev['over_cells']['max']:.2f} µm (max) of attempt 1 over every cell, and up to
  {dev['at_fov_corners_max']:.2f} µm at the far corners of the field.
* Attempt 3 is Explorer's other documented input: {k['n']} exact nucleus pairs. Explorer fits its
  own similarity through them (so it lands where attempt 2 does) and the keypoints stay
  adjustable afterwards; a matrix does not.
* On an older Explorer, or if the dialog differs: choose **"No, image is not aligned"**, place at
  least 3 keypoints, then use the image layer's dropdown → **"Import a New Alignment File"** with
  any of the three files.
* Attempt 4 needs no file: `{warped}` is the same IF resampled onto the morphology grid
  ({W} × {H} px, 0.2125 µm/px) through the same transform, which is exactly what 10x documents
  "Yes, skip alignment" for. **Never apply a CSV to this file**, and **never choose "Yes, skip
  alignment" for `{native}`**: that option places an image pixel-for-pixel on the morphology
  grid, and the native IF is on its own 0.25 µm/px grid, so it must appear unaligned.

## What each outcome means

| result | meaning |
|---|---|
| 1 works | the file's bytes were the problem; done |
| 1 fails, 2 works | Explorer requires a similarity; the keypoints file becomes the deliverable (Explorer's own adjustable fit), the affine stays the recorded transform |
| 1–2 fail, 3 works | the matrix path rejects for a reason not yet identified; keypoints deliver; the verbatim error is the next lead |
| 1–3 fail, 4 works | the transform is right and matrix import is unusable in this version; the warped image is the deliverable. Candidate: the native image's footprint reaches Xenium px x [{fp['x_range'][0]:.0f}, {fp['x_range'][1]:.0f}], y [{fp['y_range'][0]:.0f}, {fp['y_range'][1]:.0f}] against a {W} × {H} field, with a negative translation |
| 4 misaligned | contradicts the local QC (`bundle_qc.png`); stop and report |

## Report back

* Explorer version.
* For each attempt tried: OK, or the verbatim error text (a screenshot is fine) and whether an
  image layer appeared at all.
* One nucleus-scale screenshot of IF DAPI over Xenium DAPI from the first attempt that worked.
* The morphology dimensions Explorer shows for this dataset.

## Then the cell groups

**Cells → cell groups → import** → `cell_groups/hormone_positivity.csv`.

## The alignment files

`{matrix_csv}` maps **pixels of `{native}` to pixels of the Xenium morphology image**
(0.2125 µm/px), 3×3 row-major, last row `0,0,1`, no header — 10x's documented format and
direction. The IF ships at {a['if_um_per_px']} µm/px. It applies to `{native}` only. The same
file is valid input to `xeniumranger import-segmentation` and `spatialdata_io.xenium_aligned_image`.
The matrix reproduces the registration over all {a['round_trip_n_cells']:,} cells to
{a['round_trip_max_um']:.1e} µm; as written with nine digits, to
{a['round_trip_as_written_max_um']:.1e} µm.

`{sim_csv}` is the similarity described under attempt 2; `manifest.json` →
`similarity_diagnostic` records it, its cost, and that it is not the transform.

`{kp_csv}` is `fixedX,fixedY,alignmentX,alignmentY`: `fixed` in morphology pixels, `alignment` in
`{native}` pixels, every pair exact under the matrix (residual
{k['residual_px_through_shipped_matrix']:.1e} px).

One convention no 10x document settles is whether Explorer indexes an imported image by pixel
corner or pixel centre; the difference is ±{a['unresolved_convention_um']} µm and is recorded, not
assumed away.

## The positivity groups

{p['measure']}, at or above **{p['threshold_counts']:.0f} of 255**, sampled at
{p['sampled_um_per_px']} µm/px on the native image. INS is `{p['channels']['INS']}` and GCG is
`{p['channels']['GCG']}`.

| group | cells |
|---|---:|
{counts}
`Unavailable` is a cell whose polygon fell outside the imaged area — nothing was measured for
it, and it is not called negative. {p['n_without_polygon']:,} of {p['n_cells']:,} cells have no
cell polygon in the Xenium bundle.

The threshold is one number applied identically to both donors and is recorded here so it can
be changed: re-run with `--threshold`.

## Not in this bundle

* **Opal 570** is recorded FAILED for this section and is not in the image.
* **Sample AF** is disregarded per the maintainer and is not in the image.
* **No QuPath annotation of any class**, on either modality.
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--bundle-root", type=Path, default=BUNDLE_ROOT)
    ap.add_argument("--level", type=int, default=0,
                    help="qptiff pyramid level to ship. 0 is full resolution (0.25 um/px).")
    ap.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                    help="mean 8-bit counts inside a cell's polygon at or above which it is "
                         "called positive for that hormone.")
    ap.add_argument("--compression", default="zlib",
                    choices=["zlib", "lzw", "jpeg2000", "none"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--skip-image", action="store_true",
                    help="do not rewrite images/<serial>_IF.ome.tif; keep the one already there "
                         "(its manifest record is read from the file itself). Every CSV, the "
                         "cell groups, the manifest and the README are rewritten.")
    ap.add_argument("--skip-warped", action="store_true",
                    help="do not rewrite images/<serial>_IF_on_xenium_grid.ome.tif; same rules.")
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    made = []
    for serial in discover().keys():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        if sections[sections.xenium_serial == serial].empty:
            continue
        m = export_one(serial, sections, roi=args.roi, out_dir=args.out,
                       bundle_root=args.bundle_root, level=args.level,
                       compression=args.compression, skip_image=args.skip_image,
                       skip_warped=args.skip_warped, workers=args.workers,
                       threshold=args.threshold)
        made.append(m)
        print(f"[explorer] -> {args.bundle_root / f'donor_{donor}_{serial}'}\n", flush=True)

    if not made:
        print("no donor matched", file=sys.stderr)
        return 1
    (args.bundle_root / f"bundles_{args.roi}.json").write_text(
        json.dumps(made, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
