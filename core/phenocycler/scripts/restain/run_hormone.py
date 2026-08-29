#!/usr/bin/env python
"""Does protein INS rank the RNA-beta cells above the RNA-alpha cells, on the same cell?

The only check that exists on Xenium's endocrine tier. Xenium has no INS, GCG, SST or PPY, so
every endocrine call it makes rests on surrogates; the post-Xenium IF stained the hormones on
that same section, and the tissue-outline registration puts the two in one frame.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_hormone.py --donors 6442

Neither side uses a threshold. The protein side is an AUROC. The RNA side's label is
``beta_surrogate > alpha_surrogate`` **within the same cell**, a contrast rather than a gate,
computed only over cells already inside an islet -- which is what makes this a test of subtype
instead of a test of whether the registration found islets at all.
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

import geo                                                        # noqa: E402
import hormone as HZ                                              # noqa: E402
import validate as VAL                                            # noqa: E402
import xenio                                                      # noqa: E402
from run_register import (GEOJSON_DIR, OUT, SUBPANELS, XENIUM_PATHS,  # noqa: E402
                          discover, donor_of, register_one)

IMAGES = Path.home() / "IO60panc2nd" / "PostXeniumIHC" / "Images"
SECTIONS = HERE.parents[1] / "data" / "integration" / "restain_sections.csv"
#: Channel resolution to read. ~1 um/px is 4x coarser than the acquisition and still ~8 px
#: across a Xenium cell, which is enough to average over the cell rather than one pixel.
TARGET_UM = 1.0


def islet_residents(xy_if_um: np.ndarray, islets: list) -> np.ndarray:
    """Which warped Xenium cells fall inside a kept IF islet polygon.

    Two shapely 2 conventions to get right, and both fail SILENTLY as an empty or wrong
    result rather than an error:

      * the predicate is evaluated INPUT-first -- `predicate="contains"` asks whether each
        POINT contains a polygon, which is never true and returns nothing at all;
        point-in-polygon is `"within"`.
      * `query` returns index PAIRS, with row 0 indexing the INPUT and row 1 the tree -- not
        a per-point list of hits.
    """
    from shapely import points as shp_points
    from shapely.strtree import STRtree
    pts = shp_points(np.asarray(xy_if_um, dtype=float))
    tree = STRtree(islets)
    pt_idx, _tree_idx = tree.query(pts, predicate="within")
    inside = np.zeros(len(xy_if_um), dtype=bool)
    inside[pt_idx] = True
    return inside


def run_one(serial: str, gj: Path, sections: pd.DataFrame, *, roi: str = "panc") -> dict:
    rec, res, dd, isl, grid = register_one(serial, gj, roi=roi)
    donor = rec["donor_id"]
    row = sections[sections.xenium_serial == serial]
    if row.empty:
        raise HZ.HormoneError(f"serial {serial} is not recorded in {SECTIONS}")
    row = row.iloc[0]
    qptiff = IMAGES / row.image_file
    if not qptiff.exists():
        raise HZ.HormoneError(f"{qptiff} is absent")

    # Xenium cells -> IF frame. The transform maps moving (Xenium) -> fixed (IF), so this is
    # the forward direction and no inverse is involved.
    xy_if = res.transform.apply(dd.xy)
    resident = islet_residents(xy_if, isl)

    core, by_hint = VAL.endocrine_core(SUBPANELS)
    gidx = {g: i for i, g in enumerate(dd.genes)}
    L = xenio.lognorm(dd.counts)
    beta = VAL.gene_set_score(L, by_hint["beta"], gidx)
    alpha = VAL.gene_set_score(L, by_hint["alpha"], gidx)
    # A CONTRAST inside one cell, not a gate: no cutoff, and no appeal to how many cells
    # should be beta. Cells where the two tie carry no evidence either way and are dropped.
    beta_like = beta > alpha
    informative = resident & (beta != alpha)

    out = dict(rec)
    out["hormone"] = {"qptiff": str(qptiff), "target_um": TARGET_UM,
                      "sample_radius_um": HZ.SAMPLE_RADIUS_UM,
                      "n_islet_resident": int(resident.sum()),
                      "n_informative": int(informative.sum()),
                      "n_beta_like": int(beta_like[informative].sum()),
                      "n_alpha_like": int((~beta_like[informative]).sum())}

    from phenocycler.integration.rasterize import read_qptiff_channel

    # Every channel, including the ones that are not markers: the FAILED one and the AF one
    # are what the hormone numbers have to be read against, so they are measured, not skipped.
    names = HZ.channel_names(qptiff)
    per_channel: dict[str, np.ndarray] = {}
    for i, nm in enumerate(names):
        img, um_px = read_qptiff_channel(qptiff, i, target_um=TARGET_UM)
        per_channel[nm] = HZ.sample_cells(img, um_px, xy_if)
        out["hormone"].setdefault("um_per_px", float(um_px))
        del img
    out["hormone"]["channels"] = names
    out["hormone"]["diagnostics"] = HZ.channel_diagnostics(
        {k: v[resident] for k, v in per_channel.items()})

    af_name = next((n for n, c in zip(names, ["channel_dapi", "channel_opal570",
                                              "channel_opal690", "channel_opal520",
                                              "channel_sampleaf"])
                    if str(row.get(c, "")).strip().upper() == "AF"), None)
    for marker, positive_is_beta in (("INS", True), ("GCG", False)):
        ch = HZ.channel_index(qptiff, marker, row)
        vals = per_channel[names[ch]][informative]
        pos = beta_like[informative] if positive_is_beta else ~beta_like[informative]
        rec_m = {"channel_index": ch, "channel_name": names[ch],
                 "n_off_image": int(np.isnan(vals).sum()),
                 "raw": HZ.auroc_with_null(vals, pos)}
        if af_name is not None:
            # Reported, never the default: the maintainer's instruction is to disregard the
            # AF channel, and an instruction is not silently overridden by a measurement.
            rec_m["af_regressed"] = HZ.auroc_with_null(
                HZ.regress_out(vals, per_channel[af_name][informative]), pos)
        rec_m.update(rec_m["raw"])          # headline stays raw
        out["hormone"][marker] = rec_m

    # How to read the AUROCs above: at what spatial scale do the two modalities agree?
    ins_r = pd.Series(per_channel[names[HZ.channel_index(qptiff, "INS", row)]]).rank(pct=True)
    gcg_r = pd.Series(per_channel[names[HZ.channel_index(qptiff, "GCG", row)]]).rank(pct=True)
    out["hormone"]["agreement_by_scale"] = HZ.agreement_by_scale(
        xy_if[resident], (beta - alpha)[resident],
        (ins_r - gcg_r).to_numpy()[resident])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--geojson-dir", type=Path, default=GEOJSON_DIR)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)

    sections = HZ.load_section_map(SECTIONS)
    rows = []
    for serial, gj in discover(args.geojson_dir).items():
        donor = donor_of(serial, XENIUM_PATHS)
        if args.donors and donor not in args.donors:
            continue
        rec = run_one(serial, gj, sections, roi=args.roi)
        rows.append(rec)
        h = rec["hormone"]
        print(f"\n  {donor} (IF {serial})  Dice {rec['dice_after']:.3f}   "
              f"islet residents {h['n_islet_resident']:,}   informative {h['n_informative']:,} "
              f"({h['n_beta_like']:,} beta-like / {h['n_alpha_like']:,} alpha-like)")
        for m in ("INS", "GCG"):
            a = h[m]
            tgt = "beta-like" if m == "INS" else "alpha-like"
            print(f"     {m:>3} ranks {tgt:>10}: AUROC {a['auroc']:.3f}   "
                  f"null {a['null_mean']:.3f} [{a['null_p2_5']:.3f}, {a['null_p97_5']:.3f}]   "
                  f"{'clears' if a['clears_null'] else 'INSIDE NULL'}"
                  f"   off-image {a['n_off_image']:,}", flush=True)
        prof = rec["hormone"]["agreement_by_scale"]
        print("     agreement by scale (rho): " + "  ".join(
            f"{p['radius_um']:.0f}um {p['spearman']:+.3f}" for p in prof), flush=True)

    args.out.mkdir(parents=True, exist_ok=True)
    p = args.out / f"hormone_{args.roi}.json"
    p.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
