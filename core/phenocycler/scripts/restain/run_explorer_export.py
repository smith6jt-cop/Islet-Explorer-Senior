#!/usr/bin/env python
"""The registered post-Xenium IF, plus Xenium's segmentations coloured by INS/GCG positivity.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_explorer_export.py --donors 6436 6442

Three files per donor, and nothing else:

    donor_<id>_<serial>/
      images/<serial>_IF.ome.tif                 DAPI + INS + GCG, pyramidal, tiled 1024
      images/<serial>_IF_imagealignment.csv      the 3x3 Xenium Explorer imports with it
      cell_groups/hormone_positivity.csv         cell_id, group, color
      README.md, manifest.json                   how to import it, and what the numbers are

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


class ExportError(RuntimeError):
    """A bundle cannot be assembled from what is on disk -- never partially, never guessed."""


def _reset(dest: Path, *, keep_images: bool) -> list[str]:
    """Clear the bundle so it is rebuilt whole, sparing an OME-TIFF already in it if asked.

    A bundle is rewritten rather than updated, because a half-refreshed one -- new matrix,
    stale groups -- is the kind of thing that gets sent to a colleague without anyone noticing.
    ``--skip-image`` is the exception, and a reset that deleted the 2 GB image would destroy
    exactly what that flag exists to preserve.
    """
    dest = Path(dest)
    if not dest.exists():
        (dest / "images").mkdir(parents=True)
        return []
    keep = sorted((dest / "images").glob("*.ome.tif")) if keep_images else []
    stash = dest.parent / f".{dest.name}.keep"
    if keep:
        stash.mkdir(parents=True, exist_ok=True)
        for p in keep:
            p.rename(stash / p.name)
    shutil.rmtree(dest)
    (dest / "images").mkdir(parents=True)
    for p in sorted(stash.iterdir()) if stash.exists() else []:
        p.rename(dest / "images" / p.name)
    if stash.exists():
        stash.rmdir()
    return [p.name for p in keep]


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
               skip_image: bool = False, workers: int = 8,
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

    bundle = xenio._bundle_dir(donor, XENIUM_PATHS, roi)
    import pandas as pd
    cells = pd.read_parquet(bundle / "cells.parquet", columns=["cell_id"])
    cell_ids = cells.cell_id.astype(str).to_numpy()

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

    dest = Path(bundle_root) / f"donor_{donor}_{serial}"
    kept = _reset(dest, keep_images=skip_image)
    (dest / "cell_groups").mkdir(parents=True)

    align_csv = EX.write_alignment_csv(dest / "images" / f"{serial}_IF_imagealignment.csv", M)
    print(f"[explorer] {donor}: alignment matrix -> {align_csv.name}", flush=True)

    pages, ch_names, dropped = _channel_plan(qptiff, row)
    image_rec = {"skipped": True, "preserved_from_a_previous_run": kept,
                 "path": str(dest / "images" / kept[0]) if kept else None,
                 "channels": ch_names}
    if not skip_image:
        image_rec = EX.write_explorer_ome_tiff(
            qptiff, dest / "images" / f"{serial}_IF.ome.tif", pages, ch_names,
            base_um_per_px=if_um_per_px, compression=compression,
            max_levels=len(shapes) - level, workers=workers)
    image_rec["dropped_channels"] = dropped

    sampled = per_cell_values(qptiff, bundle, cell_ids, transform, row,
                              um_per_px=sample_um_per_px)
    ins, gcg = sampled["values"]["INS"], sampled["values"]["GCG"]
    labels = positivity_groups(ins, gcg, threshold)
    groups = EX.write_cell_groups(dest / "cell_groups" / "hormone_positivity.csv",
                                  cell_ids, labels, GROUP_COLORS)
    print(f"[explorer] {donor}: {groups['groups']}", flush=True)

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "donor_id": donor, "xenium_serial": serial, "roi": roi,
        "xenium_bundle": str(bundle), "qptiff": str(qptiff),
        "transform": {"source": "refined", **transform.to_dict()},
        "alignment": {
            "file": str(align_csv.relative_to(dest)),
            "matrix": np.asarray(M).tolist(),
            "maps": "exported IF image pixels -> Xenium morphology image pixels",
            "xenium_um_per_px": EX.XENIUM_UM_PER_PX,
            "if_um_per_px": if_um_per_px, "if_origin_um": if_origin_um,
            "if_pyramid_level": int(level), "if_pyramid_factor": int(f),
        },
        "image": image_rec,
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
    a, p, img = m["alignment"], m["positivity"], m["image"]
    serial = m["xenium_serial"]
    name = Path(img["path"]).name if img.get("path") else f"{serial}_IF.ome.tif"
    counts = "".join(f"| `{k}` | {v:,} |\n" for k, v in sorted(p["counts"].items()))
    return f"""# Donor {m['donor_id']} — post-Xenium IF registered to Xenium {serial}

Generated {m['generated_utc']}.

## Import

1. Open the Xenium output bundle for this section in Xenium Explorer.
2. **Images → add image** → `images/{name}`.
3. When it asks for an alignment, choose **Import a New Alignment File** (transformation
   matrix, not keypoints) → `images/{Path(a['file']).name}`.
4. **Cells → cell groups → import** → `cell_groups/hormone_positivity.csv`.

## The alignment file

3x3 row-major, no header, last row `0,0,1`, in pixels — 10x's documented format. It maps
**pixels of the image in `images/` to pixels of the Xenium morphology image**
({a['xenium_um_per_px']} µm/px). The IF ships at {a['if_um_per_px']} µm/px. The same file works
in `xeniumranger import-segmentation` and `spatialdata_io.xenium_aligned_image`.

## The positivity groups

{p['measure']}, at or above **{p['threshold_counts']:.0f} of 255**, sampled at
{p['sampled_um_per_px']} µm/px. INS is `{p['channels']['INS']}` and GCG is
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
                    help="write the matrix and the cell groups, keeping any OME-TIFF already "
                         "written. The image is the slow part; the rest is seconds.")
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
                       workers=args.workers, threshold=args.threshold)
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
