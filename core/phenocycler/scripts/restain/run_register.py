#!/usr/bin/env python
"""Register the post-Xenium IF to its own Xenium section, on tissue outlines alone.

The post-Xenium IF images the SAME PHYSICAL SECTION Xenium measured, so this is the one place
in the cohort where protein and RNA are on one cell.  What makes it reachable without
segmenting the re-stains is the maintainer's route: register the **Xenium segmentation's
tissue outline** to the IF's, and let the Xenium segmentation be the segmentation.

**No identity is used, on either side.**  The IF supplies a pixel-classifier tissue polygon;
Xenium supplies a cell footprint rasterised from ``cells.parquet``.  Neither is a cell type,
so this runs today, on all 20 donors in principle, and it is not gated on
``require_current_pheno``.  It also means the islet check below is genuinely independent: the
registration never saw an islet or a gene.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_register.py                # every export found
    ../../.venv/bin/python scripts/restain/run_register.py --donors 6442

Measured 2026-08-28 on the two exports available:

    donor 6436 (IF 0059907): Dice 0.753 -> 0.890,  125 islets >= 50 um
    donor 6442 (IF 0059947): Dice 0.922 -> 0.951,   95 islets >= 50 um

Both clear ``qc_tissue_dice_min = 0.80``.  Non-rigid refinement is skipped: ``SimpleITK`` is
not installed, and ``register_image(nonrigid=False)`` says so rather than silently producing
an affine-only result labelled as B-spline.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))                      # core/phenocycler
sys.path.insert(0, str(HERE.parents[1] / "scripts" / "xenium_celltype"))

import geo                                                     # noqa: E402
import validate as VAL                                          # noqa: E402
import xenio                                                   # noqa: E402
from phenocycler.integration import rasterize as RZ            # noqa: E402
from phenocycler.integration import register as RG             # noqa: E402

GEOJSON_DIR = Path.home() / "IO60panc2nd" / "PostXeniumIHC" / "geojson"
XENIUM_PATHS = HERE.parents[1] / "data" / "integration" / "xenium_paths.csv"
OUT = HERE.parents[3] / "data" / "development" / "restain"
SUBPANELS = HERE.parents[3] / "data" / "development" / "xenium_panels" / "subpanels"

#: Registration grid.  10 um is ~2 cell diameters -- fine enough that a tissue edge is sharp,
#: coarse enough that the ECC search is seconds rather than minutes.
PIXEL_UM = 10.0
#: The gate `integration/qc.py` already applies to a registration.  Reported, never relaxed.
DICE_MIN = 0.80

_SERIAL = re.compile(r"(\d{7})")


def discover(geojson_dir: Path = GEOJSON_DIR) -> dict[str, Path]:
    """IF serial -> geojson path, for every export present."""
    out = {}
    for p in sorted(Path(geojson_dir).glob("*.geojson*")):
        m = _SERIAL.search(p.name)
        if m:
            out[m.group(1)] = p
    return out


def donor_of(serial: str, xenium_paths: Path) -> str:
    """The donor a slide serial belongs to, from the manifest -- never from the filename.

    Two of the five local exports are named by donor (`HDL115`, `npod6591`) and three by
    serial, so a filename rule would need two cases and would be wrong the first time a
    fourth naming appeared.  The manifest already records serial -> donor.
    """
    xen = pd.read_csv(xenium_paths, dtype=str)
    hit = xen[xen.source_path.fillna("").str.contains(serial)]
    donors = sorted(set(hit.donor_id))
    if len(donors) != 1:
        raise geo.RestainGeoError(
            f"serial {serial} maps to donors {donors} in {xenium_paths}; a re-stain must "
            "resolve to exactly one donor")
    return donors[0]


def register_one(serial: str, gj: Path, *, xenium_paths: Path = XENIUM_PATHS,
                 roi: str = "panc") -> dict:
    donor = donor_of(serial, xenium_paths)
    ann = geo.load_annotations(gj)
    tissue = geo.tissue_polygon(ann)
    isl = geo.islets(ann)

    dd = xenio.load_donor(donor, xenium_paths, roi=roi, prefer="bundle")
    scale_check = geo.assert_scale_is_plausible(tissue, dd.xy)

    bx0, by0, bx1, by1 = tissue.bounds
    corners = np.array([[bx0, by0], [bx1, by1]])
    grid = RZ.grid_for(dd.xy, PIXEL_UM, other=corners, margin_um=300.0)

    fixed = geo.rasterize([tissue], grid)          # IF: the drawn tissue outline
    moving = RZ.tissue_mask(dd.xy, grid)           # Xenium: the cell footprint
    before = RG.dice(moving, fixed)
    res = RG.register_image(moving, fixed, grid, nonrigid=False)
    after = RG.dice(res.transform.warp_mask(moving, grid), fixed)

    return {
        "serial": serial, "donor_id": donor, "roi": roi,
        "geojson": str(gj), "n_islet_objects": len(ann.get(geo.ISLET, [])),
        "n_islets_kept": len(isl), "min_islet_diameter_um": geo.MIN_ISLET_DIAMETER_UM,
        "grid": grid.as_dict(), "dice_before": before, "dice_after": after,
        "dice_min": DICE_MIN, "passes_dice": bool(after >= DICE_MIN),
        "stages": list(res.stages), "notes": list(res.notes),
        "nonrigid": False, "nonrigid_reason": "SimpleITK not installed",
        "affine": np.asarray(res.transform.affine).tolist(),
        **scale_check,
    }, res, dd, isl, grid


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    found = discover(args.geojson_dir)
    if not found:
        print(f"no exports under {args.geojson_dir}")
        return 1

    rows = []
    for serial, gj in found.items():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        rec, res, dd, isl, grid = register_one(serial, gj, roi=args.roi)

        # The independent check: the registration used tissue outlines only, so islets and
        # genes are both unseen by it.  Reported for every donor, pass or fail.
        cen_xen = res.transform.inverse_apply(geo.centroids(isl))
        core, by_hint = VAL.endocrine_core(SUBPANELS)
        gidx = {g: i for i, g in enumerate(dd.genes)}
        L = xenio.lognorm(dd.counts)
        rec["endocrine"] = VAL.enrichment(VAL.gene_set_score(L, core, gidx), dd.xy, cen_xen)
        rec["endocrine_by_hint"] = {
            h: VAL.enrichment(VAL.gene_set_score(L, gs, gidx), dd.xy, cen_xen)
            for h, gs in sorted(by_hint.items()) if any(g in gidx for g in gs)}
        rows.append(rec)
        flag = "" if rec["passes_dice"] else f"   <-- BELOW {DICE_MIN}"
        e = rec["endocrine"]
        print(f"  {donor} (IF {serial}): Dice {rec['dice_before']:.3f} -> "
              f"{rec['dice_after']:.3f}  islets {rec['n_islets_kept']}/"
              f"{rec['n_islet_objects']}{flag}", flush=True)
        print(f"      endocrine core inside islets {e['mean_inside']:.4f} vs "
              f"{e['mean_outside']:.4f} elsewhere = {e['ratio']:.2f}x   "
              f"null {e['null_mean']:.4f} (p95 {e['null_p95']:.4f})   "
              f"{'clears null' if e['clears_null'] else 'INSIDE NULL'}", flush=True)
        print("      by hint: " + "  ".join(
            f"{h} {v['ratio']:.2f}x" for h, v in rec["endocrine_by_hint"].items()), flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / f"registration_{args.roi}.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {args.out / f'registration_{args.roi}.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
