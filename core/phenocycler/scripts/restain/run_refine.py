#!/usr/bin/env python
"""Fit every candidate registration on the two nucleus fields, score them alike, keep the best.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/run_refine.py --donors 6436 6442

Runs after ``export_dapi.py`` (the two DAPI planes and their micron frames), ``valis_register.py``
(the seed) and ``nuclei.py`` (the two point sets). Writes ``refined.json`` beside each donor's
export and one summary at ``data/development/restain/refine_<roi>.json``.

**Every candidate is scored by the same measurement and all of them are reported**, winner and
losers, so "this is the best of the packages available" is something a reader can check rather
than a claim. The score is the per-cell distance between a Xenium nucleus and the IF nucleus it
lands on, with a placement null beside it.

**The winner is chosen here; it is judged elsewhere.** ``run_align.py`` -- rasterised nucleus
masks against the IF DAPI, windowed normalised cross-correlation -- reads the chosen transform
afterwards and never contributes to choosing it.
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

import refine as RF                                               # noqa: E402
import registration_inputs as RI                                  # noqa: E402
from run_register import OUT                                      # noqa: E402

from phenocycler.integration.rasterize import GridSpec            # noqa: E402
from phenocycler.integration.transform import Transform           # noqa: E402

#: Fit a local field only when a global fit leaves this much structure behind. Below it the
#: residual is at the detector's own precision and a field would be fitting noise into geometry.
LOCAL_FIELD_MIN_UM = 0.30
#: Spacing of that field. Coarse on purpose: it can absorb a smooth stretch and cannot
#: shrink-wrap one nucleus field onto another.
LOCAL_FIELD_PIXEL_UM = 250.0
#: Candidates within this of the best are treated as tied, and the SIMPLEST of them is kept.
#: A tenth of the detector's own localisation error, so a difference smaller than this is not a
#: difference in registration -- and shipping extra freedom to buy it would be paying in
#: degrees of freedom for something that was never measured.
WINNER_TOLERANCE_UM = 0.05


def choose(cands):
    """The simplest candidate that is not measurably worse than the best one.

    Lowest-score-wins would always prefer more degrees of freedom, because more freedom cannot
    score worse on the data it was fitted to. Ranking by stage count inside a declared tolerance
    means a field has to EARN its place by a margin, rather than by a rounding difference.
    """
    def p50(c):
        return c.score["pair_distance_um"].get("p50", float("inf"))

    best = min(p50(c) for c in cands)
    tied = [c for c in cands if p50(c) <= best + WINNER_TOLERANCE_UM]
    return min(tied, key=lambda c: (len(c.stages), p50(c)))


def held_out(seed, xen, ifn, *, affine: bool = True) -> dict:
    """Fit on one spatial half of the section, measure on the other. Both directions.

    Mutual nearest-neighbour matching can be self-consistently WRONG: if the seed were out by
    more than half the spacing between nuclei, every point would pair with a stranger and the
    fit would converge confidently onto that. Nothing in the in-sample score would show it,
    because the score is computed on the pairs the fit chose.

    Splitting the section in two and scoring each fit on territory it never saw removes that.
    A transform that is a registration transfers; a transform that is a fit does not. The two
    halves' parameters are reported beside the scores for the same reason.
    """
    mid = float(np.median(xen[:, 0]))
    out = {"split_x_um": mid, "folds": []}
    for name, m in (("left", xen[:, 0] < mid), ("right", xen[:, 0] >= mid)):
        fit_xen, test_xen = xen[m], xen[~m]
        tf = RF.refine_icp(seed, fit_xen, ifn, affine=affine)
        out["folds"].append({
            "fitted_on": name, "n_fit": int(m.sum()), "n_test": int((~m).sum()),
            "in_sample": RF.score_transform(tf, fit_xen, ifn)["pair_distance_um"],
            "held_out": RF.score_transform(tf, test_xen, ifn)["pair_distance_um"],
            "transform": RF.describe(tf)})
    a, b = out["folds"][0]["transform"], out["folds"][1]["transform"]
    out["translation_disagreement_um"] = float(np.linalg.norm(
        np.array(a["translation_um"]) - np.array(b["translation_um"])))
    out["scale_disagreement_ppm"] = float(abs(a["scale_ppm"] - b["scale_ppm"]))
    return out


def _load(donor_dir: Path) -> dict:
    d = Path(donor_dir)
    need = ["export.json", "xenium_nuclei_um.npy", "if_nuclei_um.npy",
            "xenium_centroids_um.npy", "warped_if_um.npy"]
    missing = [n for n in need if not (d / n).exists()]
    if missing:
        raise RF.RefineError(
            f"{d.name} is missing {missing} -- run export_dapi.py, valis_register.py and "
            f"nuclei.py first")
    return {
        "meta": json.loads((d / "export.json").read_text()),
        "xen": np.load(d / "xenium_nuclei_um.npy"),
        "ifn": np.load(d / "if_nuclei_um.npy"),
        "cells": np.load(d / "xenium_centroids_um.npy"),
        "warped": np.load(d / "warped_if_um.npy"),
    }


def _valis_field_transform(donor_dir: Path, seed: Transform) -> Transform | None:
    """The transform the lane shipped on 2026-08-31: VALIS's affine plus its binned residual."""
    npz = Path(donor_dir) / "valis_displacement.npz"
    if not npz.exists():
        return None
    z = np.load(npz)
    grid = GridSpec(float(z["origin_x_um"]), float(z["origin_y_um"]), float(z["pixel_um"]),
                    int(z["height"]), int(z["width"]))
    return Transform(affine=np.array(seed.affine), displacement=z["field"], grid=grid,
                     meta={"source": "valis_transform.py 250um binned residual"})


def candidates_for(data: dict, donor_dir: Path, *, if_image=None) -> list[RF.Candidate]:
    xen, ifn = data["xen"], data["ifn"]
    seed = RF.seed_from_points(data["cells"], data["warped"])

    out = [RF.Candidate("valis_affine", seed, stages=["valis", "affine_fit"],
                        notes=["VALIS's own map reduced to an affine -- the seed, and the "
                               "baseline every other candidate has to beat"])]

    field_tf = _valis_field_transform(donor_dir, seed)
    if field_tf is not None:
        out.append(RF.Candidate("valis_field", field_tf, stages=["valis", "binned_residual"],
                                notes=["what the lane shipped before this work"]))

    try:
        out.append(RF.Candidate("nuclei_translation",
                                RF.refine_translation(seed, xen, ifn),
                                stages=["valis", "median_pair_shift"],
                                notes=["seed plus one constant: the smallest correction the "
                                       "measured residual asks for"]))
    except RF.RefineError as exc:
        out.append(RF.Candidate("nuclei_translation", seed, notes=[f"refused: {exc}"]))

    sim = RF.refine_icp(seed, xen, ifn, affine=False)
    out.append(RF.Candidate("nuclei_icp_similarity", sim, stages=["valis", "icp_similarity"],
                            notes=["rotation, uniform scale and translation on mutual "
                                   "nearest-neighbour nuclei"]))

    aff = RF.refine_icp(sim, xen, ifn, affine=True)
    out.append(RF.Candidate("nuclei_icp_affine", aff, stages=["valis", "icp_similarity",
                                                              "icp_affine"],
                            notes=["the similarity fit continued with shear and anisotropic "
                                   "scale allowed"]))

    if if_image is not None:
        m = data["meta"]
        dx, dy, peak = RF.global_translation_by_image(
            seed, xen, if_image, um_per_px=float(m["if_um_per_px"]),
            origin_um=float(m["if_origin_um"]))
        out.append(RF.Candidate(
            "image_translation", RF.translated(seed, dx, dy),
            stages=["valis", "whole_section_ncc"],
            notes=[f"intensity, not points: whole-section NCC of the warped nucleus raster "
                   f"against the high-passed IF DAPI, peak {peak:.4f}"]))

        px, py = RF.global_translation_by_phase(
            seed, xen, if_image, um_per_px=float(m["if_um_per_px"]),
            origin_um=float(m["if_origin_um"]))
        out.append(RF.Candidate(
            "phase_translation", RF.translated(seed, px, py),
            stages=["valis", "phase_correlation"],
            notes=["scikit-image phase correlation, upsample 100 -- a third mechanism for the "
                   "same translation, weighting every spatial frequency equally",
                   f"disagrees with the NCC estimate by "
                   f"{float(np.hypot(px - dx, py - dy)):.3f} um"]))
    return out


def refine_one(donor_dir: Path, *, with_image: bool = True) -> dict:
    donor_dir = Path(donor_dir)
    data = _load(donor_dir)
    meta = data["meta"]

    if_image = None
    if with_image:
        import tifffile
        if_image = tifffile.imread(str(donor_dir / "slides" / "02_if_dapi.ome.tiff"))

    cands = candidates_for(data, donor_dir, if_image=if_image)
    del if_image

    for c in cands:
        c.score = RF.score_transform(c.transform, data["xen"], data["ifn"])
        p50 = c.score["pair_distance_um"].get("p50", float("nan"))
        print(f"[{meta['donor_id']}] {c.name:24s} pair p50 {p50:6.3f} um  "
              f"p90 {c.score['pair_distance_um'].get('p90', float('nan')):6.3f}  "
              f"paired {100*c.score['paired_fraction']:5.1f}%  "
              f"within 2um {100*c.score['matched_within']['2um']:5.1f}% "
              f"(null {100*c.score['matched_within_null']['2um']:4.1f}%)", flush=True)

    # A global fit that still leaves structure gets one coarse field; one that does not, does not.
    best_global = min(
        (c for c in cands if c.transform.displacement is None),
        key=lambda c: c.score["pair_distance_um"].get("p50", float("inf")))
    left = best_global.score.get("residual_after_constant_um", 0.0)
    if left >= LOCAL_FIELD_MIN_UM:
        try:
            loc = RF.refine_local(best_global.transform, data["xen"], data["ifn"],
                                  pixel_um=LOCAL_FIELD_PIXEL_UM)
            c = RF.Candidate(f"{best_global.name}+local", loc,
                             stages=best_global.stages + ["local_field"],
                             notes=[f"a {LOCAL_FIELD_PIXEL_UM:g} um field over "
                                    f"{best_global.name}, fitted because its residual still "
                                    f"held {left:.3f} um of structure after a constant"])
            c.score = RF.score_transform(c.transform, data["xen"], data["ifn"])
            cands.append(c)
            print(f"[{meta['donor_id']}] {c.name:24s} pair p50 "
                  f"{c.score['pair_distance_um']['p50']:6.3f} um  "
                  f"p90 {c.score['pair_distance_um']['p90']:6.3f}", flush=True)
        except RF.RefineError as exc:
            print(f"[{meta['donor_id']}] local field refused: {exc}", flush=True)
    else:
        print(f"[{meta['donor_id']}] no local field: {best_global.name} leaves {left:.3f} um "
              f"after a constant, under the {LOCAL_FIELD_MIN_UM} um it would take to justify "
              f"fitting geometry to it", flush=True)

    winner = choose(cands)

    seed = RF.seed_from_points(data["cells"], data["warped"])
    cv = held_out(seed, data["xen"], data["ifn"], affine="icp_affine" in winner.stages)
    for f in cv["folds"]:
        print(f"[{meta['donor_id']}] hold-out: fitted on the {f['fitted_on']} half "
              f"({f['n_fit']:,} nuclei) -> in-sample p50 {f['in_sample']['p50']:.3f} um, "
              f"HELD OUT p50 {f['held_out']['p50']:.3f} um", flush=True)
    print(f"[{meta['donor_id']}] the two halves disagree by "
          f"{cv['translation_disagreement_um']:.3f} um of translation and "
          f"{cv['scale_disagreement_ppm']:.0f} ppm of scale", flush=True)

    winner.transform.meta.update({"candidate": winner.name, "stages": winner.stages})
    winner.transform.save(donor_dir, stem="refined")

    report = {
        "donor_id": meta["donor_id"], "serial": meta["serial"], "roi": meta["roi"],
        "n_xenium_nuclei": int(len(data["xen"])), "n_if_nuclei": int(len(data["ifn"])),
        # A Xenium nucleus can only be matched if the IF detector found one there, so the
        # matched fraction cannot exceed this however good the registration is. Recorded so the
        # two are never confused.
        "detection_ceiling": float(len(data["ifn"]) / max(1, len(data["xen"]))),
        "winner": winner.name,
        "held_out": cv,
        "transform_file": str(donor_dir / "refined.json"),
        "candidates": [{"name": c.name, "stages": c.stages, "notes": c.notes,
                        "transform": RF.describe(c.transform), "score": c.score}
                       for c in cands],
        "constants": {"PAIR_MAX_UM": RF.PAIR_MAX_UM, "MATCH_RADII_UM": list(RF.MATCH_RADII_UM),
                      "NULL_OFFSET_UM": RF.NULL_OFFSET_UM,
                      "LOCAL_FIELD_MIN_UM": LOCAL_FIELD_MIN_UM,
                      "LOCAL_FIELD_PIXEL_UM": LOCAL_FIELD_PIXEL_UM,
                      "WINNER_TOLERANCE_UM": WINNER_TOLERANCE_UM},
    }
    # The cohort summary is convenient, but a subset or resumed run necessarily replaces it.
    # Keep the complete candidate/hold-out record beside the transform that downstream stages
    # consume so parallel batches cannot erase one another's provenance.
    report_path = donor_dir / "refinement.json"
    tmp = report_path.with_name(f".{report_path.name}.tmp")
    tmp.write_text(json.dumps(report, indent=2))
    tmp.replace(report_path)
    return report


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
    ap.add_argument("--no-image", action="store_true",
                    help="skip the intensity candidate (it reads the whole IF plane)")
    args = ap.parse_args(argv)

    records = RI.load_registration_inputs(
        args.inventory, args.xenium_paths, args.images_dir, args.geojson_dir,
        args.xenium_root, roi=args.roi)
    selected = RI.select_registration_inputs(records, args.donors, context="run_refine")
    for record in selected:
        RI.require_registration_files(record)

    rows = []
    for record in selected:
        rows.append(refine_one(args.out / "valis" / f"donor_id={record.donor_id}",
                               with_image=not args.no_image))
    if not rows:
        print("no donor matched", file=sys.stderr)
        return 1
    p = args.out / f"refine_{args.roi}.json"
    p.write_text(json.dumps(rows, indent=2))
    print(f"-> {p}", flush=True)
    for r in rows:
        w = next(c for c in r["candidates"] if c["name"] == r["winner"])
        print(f"[{r['donor_id']}] WINNER {r['winner']}: pair p50 "
              f"{w['score']['pair_distance_um']['p50']:.3f} um, "
              f"p90 {w['score']['pair_distance_um']['p90']:.3f}, "
              f"{100*w['score']['matched_within']['2um']:.1f}% of Xenium nuclei within 2 um "
              f"(null {100*w['score']['matched_within_null']['2um']:.1f}%, detection ceiling "
              f"{100*r['detection_ceiling']:.1f}%)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
