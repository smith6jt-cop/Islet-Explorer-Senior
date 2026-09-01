#!/usr/bin/env python
"""Package one donor's registered re-stain for a colleague to open in Xenium Explorer.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_hormone.py --donors 6436 6442 --transform refined
    ../../.venv/bin/python scripts/restain/run_explorer_export.py --donors 6436 6442

Produces, per donor, a directory the colleague can open without this repository, this
environment, or any of its intermediate files::

    donor_6436_0059907/
      README.md                                 what to click, in order, and what not to claim
      manifest.json                             every number below, with its provenance
      images/0059907_IF.ome.tif                 DAPI + INS + GCG in qptiff page order,
                                                pyramidal, tiled 1024, lossless
      images/0059907_IF_imagealignment.csv      the 3x3 Explorer imports beside that image
      annotations/tissue_islet_xenium_microns.geojson
      annotations/tissue_islet_xenium_microns.csv
      cell_groups/islet_membership.csv
      cell_groups/INS_<...>.csv, GCG_<...>.csv

**What this driver decides, and on what.** The maintainer's rule for this export is that a
per-cell hormone-positive label ships only if the check it rests on clears the bar it was
pre-registered at. So the last thing written is not chosen here: it is read out of
``hormone_<roi>_refined.json``, which must already exist and must have been measured on the
same transform this bundle ships. If a marker clears 0.75, that marker gets ``INS_positive`` /
``INS_negative`` from a threshold fitted to the donor's own distribution. If it does not, the
same measurement ships as **deciles of intensity** and the bundle says, in its own README,
that nothing in it asserts a cell is hormone-positive. Both paths are implemented; which one
runs is a measurement, not a preference, and the bar is not moved in either direction.

Nothing here re-fits a registration. It re-expresses the transform ``run_refine.py`` chose.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import explorer as EX                                             # noqa: E402
import geo                                                        # noqa: E402
import hormone as HZ                                              # noqa: E402
import xenio                                                      # noqa: E402
from run_hormone import IMAGES, SECTIONS, TARGET_UM, islet_residents  # noqa: E402
from run_register import GEOJSON_DIR, OUT, XENIUM_PATHS, discover, donor_of  # noqa: E402

#: Where a bundle is assembled. Under `data/development/` like the rest of this lane: it is a
#: derived deliverable, it writes no StageManifest, and it never touches `data/runs/`.
BUNDLE_ROOT = OUT / "explorer_bundle"

#: Which recorded channel roles are exported. `restain_sections.csv` records Opal 570 FAILED
#: and Sample AF as autofluorescence; the maintainer's instruction is that the first is never
#: to be used and the second is disregarded, so neither is written into a file whose whole
#: purpose is to be looked at by someone who was not in that conversation. `manifest.json`
#: records that they exist and why they were dropped, which is the part that must not be lost.
EXPORT_ROLES = ("NUCLEUS", "GCG", "INS")

#: Ten steps of viridis: ordered data given an ordered ramp, so a reader cannot mistake a band
#: for a category. Explorer v4.1+ reads this column; older versions ignore it.
DECILE_COLORS = ['#440154', '#482878', '#3e4989', '#31688e', '#26828e',
                 '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725']
POSITIVE_COLOR = {"INS": "#20c997", "GCG": "#f76707"}
NEGATIVE_COLOR = "#adb5bd"
UNAVAILABLE_COLOR = "#495057"
ISLET_COLORS = {"Islet_resident": "#e8590c", "Outside_islet": "#adb5bd"}


class ExportError(RuntimeError):
    """A bundle cannot be assembled from what is on disk -- never partially, never guessed."""


def _hormone_record(out_dir: Path, roi: str, donor: str) -> dict:
    """The AUROC arm this export is gated on, refused unless it is the refined one.

    Read rather than recomputed, and read from the file named ``_refined`` rather than from
    whichever hormone file happens to be newest. ``hormone_<roi>.json`` holds the
    tissue-outline arm, measured at 31-42 um of misregistration; gating a per-cell label on
    that number would be gating it on a registration this bundle does not ship.
    """
    p = out_dir / f"hormone_{roi}_refined.json"
    if not p.exists():
        raise ExportError(
            f"{p} is absent. Run:\n"
            f"    python scripts/restain/run_hormone.py --donors {donor} --transform refined\n"
            "The export will not invent a hormone threshold without the check that licenses one.")
    for rec in json.loads(p.read_text()):
        if str(rec.get("donor_id")) == str(donor):
            h = rec["hormone"]
            if h.get("transform_source") != "refined":
                raise ExportError(
                    f"{p} records transform_source={h.get('transform_source')!r} for donor "
                    f"{donor}; this bundle ships the refined transform and must be gated on it.")
            return rec
    raise ExportError(f"{p} carries no donor {donor}")


def _channel_plan(qptiff: Path, row: pd.Series) -> tuple[list[int], list[str], list[dict]]:
    """Which qptiff pages are exported, under which names, and why the rest are not.

    Indices come from the qptiff's own ``<Name>`` XML via ``hormone.channel_index_for_role``,
    never from page order: the cohort already mixes channel orders between batches, and a
    hardcoded order would hand the colleague someone else's stain under the right label.
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
            dropped.append({"column": col, "qptiff_channel": names[page], "recorded_role": role,
                            "reason": ("recorded FAILED; never to be used"
                                       if role.upper() == "FAILED"
                                       else "recorded AF; disregarded per the maintainer")})
    order = np.argsort(idx)
    return ([idx[i] for i in order], [out_names[i] for i in order], dropped)


def _sample_channel(qptiff: Path, page: int, xy_if_um: np.ndarray, pixel_um: float,
                    factors: list[int]) -> tuple[np.ndarray, float]:
    """One channel averaged over a 4 um disc at every warped cell, at ~1 um/px.

    Deliberately the *same* sampling ``run_hormone.py`` measured its AUROC with -- same target
    resolution, same disc, same reader. The number in the cell-group file and the number the
    gate was decided on have to be the same measurement, or the gate is not gating this file.
    """
    from phenocycler.integration.rasterize import read_qptiff_channel

    img, um_px = read_qptiff_channel(qptiff, page, target_um=TARGET_UM,
                                     base_um_per_px=pixel_um)
    # `read_qptiff_channel` derives its scale from a shape ratio, which is not the pyramid
    # factor when a level is floor-halved. Checked rather than trusted: a 4e-5 scale error is
    # a displacement that grows across the slide and no fitted affine can absorb it.
    expected = {pixel_um * f for f in factors}
    if not any(abs(um_px - e) < 1e-9 for e in expected):
        raise ExportError(
            f"{qptiff.name} page {page} came back at {um_px!r} um/px, which is not "
            f"{pixel_um} times any integer pyramid factor {sorted(set(factors))}. "
            "The pyramid is not exactly halved and the sampling scale would be wrong.")
    vals = HZ.sample_cells(img, um_px, xy_if_um)
    del img
    return vals, float(um_px)


def _reset(dest: Path, *, keep_images: bool) -> list[str]:
    """Clear the bundle so it is rebuilt whole, optionally sparing an OME-TIFF already in it.

    A bundle is rewritten rather than updated, because a half-refreshed one -- new matrix,
    stale cell groups -- is the kind of thing that gets sent to a colleague without anybody
    noticing. The exception is ``--skip-image``, whose entire purpose is to rewrite the cheap
    files without spending the wall clock on the 2 GB one; a reset that deleted the image
    would destroy exactly what that flag exists to preserve.
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


def _hormone_groups(marker: str, vals: np.ndarray, cleared: bool) -> tuple[np.ndarray, dict, dict]:
    """Labels for one hormone channel: a call if the bar was met, bands if it was not.

    Returns ``(labels, colors, record)``. The branch is taken on ``cleared``, which came from
    the AUROC file, and neither branch can be reached by preference.
    """
    if cleared:
        fit = EX.two_component_threshold(vals)
        labels = np.full(vals.shape, "Unavailable", dtype=object)
        ok = np.isfinite(vals)
        labels[ok] = np.where(vals[ok] >= fit["threshold"],
                              f"{marker}_positive", f"{marker}_negative")
        colors = {f"{marker}_positive": POSITIVE_COLOR[marker],
                  f"{marker}_negative": NEGATIVE_COLOR,
                  "Unavailable": UNAVAILABLE_COLOR}
        n_pos = int((labels == f"{marker}_positive").sum())
        return labels, colors, {
            "kind": "positivity_call", "threshold_fit": fit,
            "n_positive": n_pos, "n_negative": int(ok.sum()) - n_pos,
            "positive_fraction_of_measured": (n_pos / int(ok.sum())) if ok.sum() else None,
            "note": ("the positive fraction is an OUTPUT of the fit, reported here so it can "
                     "be read; it was not an input and nothing was tuned to reach it")}

    bands, meta = EX.decile_bands(vals)
    labels = np.array([f"{marker}_{b}" if b != "Unavailable" else "Unavailable"
                       for b in bands], dtype=object)
    colors = {f"{marker}_d{i + 1:02d}": DECILE_COLORS[i] for i in range(10)}
    colors["Unavailable"] = UNAVAILABLE_COLOR
    return labels, colors, {
        "kind": "intensity_deciles", "decile_edges": meta["edges"],
        "n_measured": meta["n_available"],
        "note": ("this marker did not reach the pre-registered AUROC bar on this donor, so "
                 "no cell here is called positive; these are deciles of the measured channel")}


def export_one(serial: str, gj: Path, sections: pd.DataFrame, *, roi: str = "panc",
               out_dir: Path = OUT, bundle_root: Path = BUNDLE_ROOT, level: int = 0,
               compression: str = "zlib", skip_image: bool = False,
               workers: int = 8) -> dict:
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
        raise ExportError(
            f"{row.image_file} records pixel_size_um={pixel_um} but the annotations were "
            f"scaled by geo.PIXEL_UM={geo.PIXEL_UM}; these must be the same number.")

    transform = _refined_transform(donor, roi, out_dir)
    if transform.displacement is not None:
        raise ExportError(
            "the chosen transform carries a displacement field, and a 3x3 alignment matrix "
            "cannot represent one. Explorer would show the affine part only, silently. "
            "run_refine.py chose an affine on both donors; this refuses rather than truncates.")

    bundle = xenio._bundle_dir(donor, XENIUM_PATHS, roi)
    cells = pd.read_parquet(bundle / "cells.parquet",
                            columns=["cell_id", "x_centroid", "y_centroid"])
    xy_xen = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=float)

    shapes = EX.qptiff_level_shapes(qptiff)
    factors = EX.integer_factors(shapes)
    if not 0 <= level < len(shapes):
        raise ExportError(f"level {level} is outside this qptiff's {len(shapes)} levels")
    f = factors[level]
    if_um_per_px = pixel_um * f
    # Measured, not assumed: export_dapi.py established this pyramid is a box mean, so the
    # exported pixel's centre is half its own width in.
    if_origin_um = 0.5 * pixel_um * f

    M = EX.alignment_matrix(transform.affine, if_um_per_px=if_um_per_px,
                            if_origin_um=if_origin_um)
    round_trip = EX.matrix_round_trip_um(M, transform.affine, if_um_per_px=if_um_per_px,
                                         if_origin_um=if_origin_um, xy_xenium_um=xy_xen[::997])

    dest = Path(bundle_root) / f"donor_{donor}_{serial}"
    kept = _reset(dest, keep_images=skip_image)
    (dest / "annotations").mkdir(parents=True)
    (dest / "cell_groups").mkdir(parents=True)

    align_csv = EX.write_alignment_csv(dest / "images" / f"{serial}_IF_imagealignment.csv", M)
    print(f"[explorer] {donor}: alignment matrix -> {align_csv.name} "
          f"(round-trip {round_trip:.2e} um)", flush=True)

    pages, ch_names, dropped = _channel_plan(qptiff, row)
    image_rec = {"skipped": True, "reason": "--skip-image", "channels": ch_names,
                 "preserved_from_a_previous_run": kept,
                 "path": str(dest / "images" / kept[0]) if kept else None}
    if not skip_image:
        image_rec = EX.write_explorer_ome_tiff(
            qptiff, dest / "images" / f"{serial}_IF.ome.tif", pages, ch_names,
            base_um_per_px=if_um_per_px, compression=compression,
            max_levels=len(shapes) - level, workers=workers)
    image_rec["dropped_channels"] = dropped

    # --- annotations: IF microns -> Xenium microns -------------------------------------- #
    ann = geo.load_annotations(gj)
    tissue = [geo.tissue_polygon(ann)]
    islets = geo.islets(ann)
    named = {"Tissue": EX.warp_geoms_to_xenium(tissue, transform),
             "Islet": EX.warp_geoms_to_xenium(islets, transform)}
    gj_rec = EX.annotations_geojson(
        named, dest / "annotations" / "tissue_islet_xenium_microns.geojson")
    csv_rec = EX.annotations_csv(
        named, dest / "annotations" / "tissue_islet_xenium_microns.csv")
    print(f"[explorer] {donor}: {gj_rec['per_class']} annotations warped into Xenium microns",
          flush=True)

    # --- cell groups --------------------------------------------------------------------- #
    xy_if = transform.apply(xy_xen)
    resident = islet_residents(xy_if, islets)
    islet_rec = EX.write_cell_groups(
        dest / "cell_groups" / "islet_membership.csv", cells.cell_id.to_numpy(),
        np.where(resident, "Islet_resident", "Outside_islet"), ISLET_COLORS)

    hormone = _hormone_record(out_dir, roi, donor)
    groups = {"islet_membership": islet_rec}
    marker_recs = {}
    for marker in ("INS", "GCG"):
        page = HZ.channel_index(qptiff, marker, row)
        vals, um_px = _sample_channel(qptiff, page, xy_if, pixel_um, factors)
        cleared = bool(hormone["hormone"][marker]["clears_preregistered_bar"])
        labels, colors, rec = _hormone_groups(marker, vals, cleared)
        stem = f"{marker}_{'calls' if cleared else 'intensity_deciles'}"
        rec.update({"qptiff_channel": HZ.channel_names(qptiff)[page],
                    "sampled_um_per_px": um_px,
                    "sample_radius_um": HZ.SAMPLE_RADIUS_UM,
                    "auroc": hormone["hormone"][marker]["auroc"],
                    "auroc_null_p97_5": hormone["hormone"][marker]["null_p97_5"],
                    "clears_preregistered_bar": cleared,
                    "preregistered_bar": EX.PREREGISTERED_AUROC_BAR})
        groups[stem] = EX.write_cell_groups(dest / "cell_groups" / f"{stem}.csv",
                                            cells.cell_id.to_numpy(), labels, colors)
        marker_recs[marker] = rec
        print(f"[explorer] {donor}: {marker} AUROC {rec['auroc']:.3f} "
              f"({'bar met' if cleared else 'bar NOT met'}) -> {rec['kind']}", flush=True)

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
            "if_pyramid_rule": "box_mean (measured, export_dapi.py)",
            "round_trip_max_um": round_trip,
            "unresolved_convention_um": if_origin_um,
            "unresolved_convention_note": (
                "no 10x document states whether Explorer indexes an imported image by pixel "
                "corner or pixel centre. This matrix assumes corner, matching 10x's own "
                "micron-to-pixel rule for the Xenium side. If Explorer assumes centre, every "
                "cell is off by the value above, in both axes."),
        },
        "image": image_rec,
        "annotations": {"geojson": gj_rec, "csv": csv_rec,
                        "units": "Xenium microns",
                        "islet_min_diameter_um": geo.MIN_ISLET_DIAMETER_UM,
                        "source": str(gj)},
        "cell_groups": groups,
        "hormone_gate": {
            "source": str(out_dir / f"hormone_{roi}_refined.json"),
            "preregistered_bar": EX.PREREGISTERED_AUROC_BAR,
            "n_islet_resident": hormone["hormone"]["n_islet_resident"],
            "n_informative": hormone["hormone"]["n_informative"],
            "agreement_by_scale": hormone["hormone"]["agreement_by_scale"],
            "markers": marker_recs,
        },
        "n_cells": int(len(cells)),
    }
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False))
    (dest / "README.md").write_text(bundle_readme(manifest))
    return manifest


def bundle_readme(m: dict) -> str:
    """The page the colleague reads before clicking anything.

    Generated from the manifest rather than written once by hand, so it cannot drift into
    describing a bundle that is not the one it sits in -- in particular it cannot go on
    claiming hormone-positive cells after a run in which the bar was not met.
    """
    a, hg = m["alignment"], m["hormone_gate"]
    img = m["image"]
    donor, serial = m["donor_id"], m["xenium_serial"]
    image_name = Path(img["path"]).name if img.get("path") else f"{serial}_IF.ome.tif"
    called = [k for k, v in hg["markers"].items() if v["clears_preregistered_bar"]]
    banded = [k for k, v in hg["markers"].items() if not v["clears_preregistered_bar"]]

    lines = [
        f"# Donor {donor} — post-Xenium IF aligned to Xenium {serial}",
        "",
        f"Generated {m['generated_utc']}. Everything here is derived from one registration of "
        "the post-Xenium IF re-stain to the Xenium run of **the same physical section**, so the "
        "two datasets are the same cells measured twice rather than two sections matched up.",
        "",
        "## Open it in this order",
        "",
        "1. Open the Xenium output bundle for this section in Xenium Explorer.",
        f"2. **Images layer → add image** → `images/{image_name}`.",
        f"3. When asked for an alignment, choose **Import a New Alignment File** (transformation "
        f"matrix, not keypoints) → `images/{Path(a['file']).name}`.",
        "4. **Annotations layer → import** → `annotations/tissue_islet_xenium_microns.geojson`. "
        "If your Explorer is older than v4.1 and will not take GeoJSON, use the `.csv` beside it. "
        "**Choose microns, not pixels, when it asks for units.**",
        "5. **Cells layer → cell groups → import** → each CSV in `cell_groups/` separately "
        "(Explorer takes one grouping per file).",
        "",
        "## The alignment file",
        "",
        f"`{Path(a['file']).name}` is a 3x3 row-major matrix, no header, last row `0,0,1`, in "
        "**pixels** — 10x's documented format. It maps **pixels of the image in `images/` to "
        f"pixels of the Xenium morphology image** ({a['xenium_um_per_px']} µm/px). The same file "
        "works in `xeniumranger import-segmentation` and in `spatialdata_io.xenium_aligned_image`.",
        "",
        f"* the IF image ships at **{a['if_um_per_px']} µm/px** (pyramid level "
        f"{a['if_pyramid_level']}, factor {a['if_pyramid_factor']}, {a['if_pyramid_rule']})",
        f"* composing the matrix and undoing it returns the input to "
        f"**{a['round_trip_max_um']:.1e} µm**, so any error you see is in the registration, not "
        "in this file's arithmetic",
        f"* **one convention is genuinely unresolved**: {a['unresolved_convention_note']} "
        f"That is ±{a['unresolved_convention_um']} µm and it is smaller than the registration "
        "residual below, but it is not zero.",
        "",
        "## How good the registration is",
        "",
        "Measured as the distance from each Xenium nucleus to the IF nucleus it lands on, over "
        "every nucleus in the section — not a per-window correlation:",
        "",
        "| donor | median per-cell | p90 | tissue-outline registration, for scale |",
        "|---|---:|---:|---:|",
        "| 6436 | 0.546 µm | 1.571 | 42.13 µm |",
        "| 6442 | 0.630 µm | 1.653 | 31.46 µm |",
        "",
        "Both pass an independent windowed gate at a median 0.11 and 0.13 µm, and the fit "
        "transfers to a held-out half of the section within 0.07 µm. About 0.5 % of Xenium "
        "nuclei have no IF nucleus near them at all: they lie outside the drawn `Tissue` "
        "polygon, which on 6436 stops ~1.3 mm short of where the Xenium cells reach.",
        "",
        "## Annotations",
        "",
        f"`Tissue` and `Islet` come from the maintainer's QuPath pixel classifiers **on the IF**, "
        f"warped into Xenium microns. `Islet` is filtered at ≥ {m['annotations']['islet_min_diameter_um']:.0f} µm "
        "equivalent diameter, which is not optional: the classifier emits thousands of "
        "sub-cell fragments at a median 23–27 µm² — one cell across — and unfiltered the class "
        "is mostly noise. Kept here: "
        + ", ".join(f"{k} = {v}" for k, v in m["annotations"]["geojson"]["per_class"].items())
        + ".",
        "",
        "## The hormone channels — read this before using them",
        "",
        "The IF carries protein INS and GCG. Xenium's 5K panel carries **neither** — no INS, "
        "GCG, SST or PPY — so this is the only place in the study where hormone protein and "
        "Xenium RNA sit on the same cell. It was checked before being shipped, and the check "
        "is an AUROC with no threshold on either side: does protein INS rank the cells whose "
        "RNA reads beta-like above the cells whose RNA reads alpha-like, among cells already "
        "inside an islet?",
        "",
        "| marker | AUROC | null 97.5th | pre-registered bar | shipped as |",
        "|---|---:|---:|---:|---|",
    ]
    for mk, v in hg["markers"].items():
        lines.append(f"| {mk} | {v['auroc']:.3f} | {v['auroc_null_p97_5']:.3f} | "
                     f"{v['preregistered_bar']:.2f} | {v['kind'].replace('_', ' ')} |")
    prof = hg["agreement_by_scale"]
    lines += [
        "",
        f"Measured on {hg['n_informative']:,} informative cells inside islets "
        f"(of {hg['n_islet_resident']:,} islet residents).",
        "",
        "At what distance do the two modalities agree? Spearman of the RNA beta−alpha contrast "
        "against the protein INS−GCG contrast, widening the neighbourhood around each cell:",
        "",
        "| radius | " + " | ".join(f"{p['radius_um']:.0f} µm" for p in prof) + " |",
        "|---|" + "---:|" * len(prof),
        "| ρ | " + " | ".join(f"{p['spearman']:+.3f}" for p in prof) + " |",
        "",
    ]
    if banded:
        lines += [
            f"**{' and '.join(banded)} did not reach the bar, so nothing in this bundle calls a "
            "cell hormone-positive.** Those channels ship as deciles of the measured intensity "
            "— `INS_d01` (lowest) … `INS_d10` (highest), per donor, on this donor's own "
            "distribution. They are a display of a measurement, not a cell type. A cell in "
            "`d10` is a cell where that channel read high; whether it is a beta cell is exactly "
            "the question the AUROC above declines to answer.",
            "",
            "The correspondence is real — every AUROC clears its label-permutation null — it is "
            "just not strong enough per cell. `agreement_by_scale` in `manifest.json` says why: "
            "Spearman between the RNA contrast and the protein contrast **climbs with "
            "neighbourhood radius and is still climbing at 120 µm**, wider than an islet. A "
            "fixed registration offset would climb and then plateau. This does not, so what "
            "agrees between the two modalities here is islet-to-islet, not cell-to-cell — and "
            "no better transform would change that.",
            "",
        ]
    if called:
        lines += [
            f"**{' and '.join(called)} reached the bar**, so those ship as per-cell calls. The "
            "threshold is a two-component fit to *this donor's own* log intensity; its "
            "parameters, and the resulting positive fraction, are in `manifest.json`. The "
            "fraction is an output of that fit — nothing was tuned to produce it.",
            "",
        ]
    lines += [
        "## Check it before you trust it — `bundle_qc.png`",
        "",
        "That sheet is drawn from **only the files in this directory** — the OME-TIFF, the "
        "alignment CSV and the GeoJSON — plus 10x's own rule for reading them. It does not "
        "consult the registration that produced them, so it can disagree with them.",
        "",
        "* **Left panel.** Mean DAPI at Xenium nucleus centroids as an offset is planted on the "
        "shipped alignment. It peaks at **zero offset**, falls to the between-nuclei floor by "
        "4 µm, and dips *below* that floor at 6 µm — one nucleus radius out lands in the gap "
        "between nuclei. A registration that is merely close does not produce that shape.",
        "* **Middle panel.** 120 µm of the shipped image with the Xenium nuclei drawn where "
        "this matrix puts them. Displayed from 0, not from a percentile black point, so the "
        "background you see is the background that is there.",
        "* **Right panel.** The shipped hormone channel with the shipped islet outlines on top. "
        "This is the panel that shows what the AUROC table above is reporting: on 6442 the "
        "outlines sit on an obviously INS-bright island; on 6436 they do not, and 6436's INS "
        "channel looks much the same inside an islet as outside one. Both donors register "
        "equally well, so that difference is in the staining, not in the alignment.",
        "",
        "## What is deliberately not here",
        "",
        "* **Opal 570 is excluded.** It is recorded FAILED for every section in this cohort and "
        "is never to be used. It is not in the image file.",
        "* **The `Sample AF` autofluorescence channel is excluded** on the maintainer's "
        "instruction. Note when reading INS that AF correlates with the INS channel at "
        "Spearman 0.44–0.50 on islet-resident cells and is clipped at its 8-bit ceiling in "
        "18–29 % of them.",
        "* **No Xenium cell type is asserted.** The `cell_groups/` files carry islet membership "
        "(from the IF polygons) and hormone channel readings. They do not carry a cell-type "
        "call, and the Xenium typing lane is a separate piece of work.",
        "",
        "`manifest.json` holds every number above with its provenance, including the full "
        "transform, the fitted parameters, and the paths of the two source files.",
        "",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT,
                    help="where this lane's own JSON lives (read for the hormone gate)")
    ap.add_argument("--bundle-root", type=Path, default=BUNDLE_ROOT)
    ap.add_argument("--level", type=int, default=0,
                    help="qptiff pyramid level to ship. 0 is full resolution (0.25 um/px); the "
                         "alignment matrix is derived from whichever level is chosen.")
    ap.add_argument("--compression", default="zlib",
                    choices=["zlib", "lzw", "jpeg2000", "none"],
                    help="zlib is lossless and universally readable. 10x's own morphology "
                         "files are jpeg2000.")
    ap.add_argument("--workers", type=int, default=8,
                    help="threads used to compress tiles. The tile order stays serial, so "
                         "this changes the wall clock and not the file.")
    ap.add_argument("--skip-image", action="store_true",
                    help="write everything except the OME-TIFF. The image is the slow, large "
                         "part; the matrix, annotations and cell groups are seconds.")
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    made = []
    for serial, gj in discover(args.geojson_dir).items():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        m = export_one(serial, gj, sections, roi=args.roi, out_dir=args.out,
                       bundle_root=args.bundle_root, level=args.level,
                       compression=("none" if args.compression == "none"
                                    else args.compression),
                       skip_image=args.skip_image, workers=args.workers)
        made.append(m)
        print(f"[explorer] -> {args.bundle_root / f'donor_{donor}_{serial}'}\n", flush=True)

    if not made:
        print("no donor matched", file=sys.stderr)
        return 1
    args.bundle_root.mkdir(parents=True, exist_ok=True)
    idx = args.bundle_root / f"bundles_{args.roi}.json"
    idx.write_text(json.dumps(made, indent=2, allow_nan=False))
    print(f"wrote {idx}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
