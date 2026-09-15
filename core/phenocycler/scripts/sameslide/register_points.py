#!/usr/bin/env python3
"""
Point-set registration of a same-slide PhenoCycler/Xenium pair. No image reads.

Direction is **Xenium (moving) -> PhenoCycler (fixed)**, matching ``cfg.fixed_modality``
and ``sameslide.pair_donor``'s ``moving, fixed = (xen, pheno)``. Downstream work warps
Xenium geometry FORWARD into PhenoCycler space rather than inverting: ``Transform.invert()``
returns affine-only and drops any displacement field, so the inverse is a different
transform, not merely a slower one.

**Why a moment seed rather than ``register.register_point_set``.** That function pairs
STRUCTURES (islets) and needs both sides already callable. Xenium 0080870 cannot supply them
-- ``structures.hormone_seed_mask`` seeds islets from INS/GCG/SST, and none of the three is
on either Xenium panel. And the two frames differ by ~9.4 mm in y (the qptiff holds other
sections above the annotated pancreas), so ICP from identity has nothing to converge on.
Since this is one physical section, a rigid seed from the two clouds' own centroids and
principal axes is enough to get ICP into its basin, and it uses nothing but the coordinates
already in hand.

**Acceptance.** ``align.py``'s DAPI-window NCC gate needs image reads this lane skips, so it
does NOT apply. The gate here is the nucleus-pair statistic: ``refine.score_transform``'s
mutual-pair quantiles and ``matched_within`` radii, each against the 25 um placement null,
plus a spatial split-half hold-out so the number is not circular. That substitution is
recorded in the artifact.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "restain"))

from phenocycler.integration.register import icp  # noqa: E402
from phenocycler.integration.transform import Transform, from_rigid  # noqa: E402

import refine as refine_mod  # noqa: E402  (scripts/restain/refine.py)

#: ICP is annealed: a wide gate first so the seed's error is inside it, then tightened so
#: the fit is driven by cells that actually correspond rather than by whatever was nearest.
ICP_SCHEDULE_UM = (200.0, 50.0, 15.0, 5.0)

#: Fit a displacement field only when a global fit leaves this much structured residual.
#: Same parsimony rule as run_refine.LOCAL_FIELD_MIN_UM -- fewer stages wins on a tie.
LOCAL_FIELD_MIN_UM = 0.30
LOCAL_FIELD_PIXEL_UM = 250.0


def principal_frame(xy: np.ndarray) -> tuple[np.ndarray, float]:
    """Centroid and principal-axis angle (degrees) of a point cloud."""
    c = xy.mean(axis=0)
    cov = np.cov((xy - c).T)
    w, v = np.linalg.eigh(cov)
    major = v[:, int(np.argmax(w))]
    return c, float(np.degrees(np.arctan2(major[1], major[0])))


def seed_candidates(moving_xy: np.ndarray, fixed_xy: np.ndarray) -> list[Transform]:
    """Rigid seeds from the two clouds' moments.

    A principal axis is defined only up to sign, so the axis difference and its 180-degree
    flip are both plausible; a mirror is included so it can be *scored and rejected* rather
    than assumed away. A same-section pair should not be mirrored, so a mirrored winner is a
    finding worth surfacing, not a result to use quietly.
    """
    cm, am = principal_frame(moving_xy)
    cf, af = principal_frame(fixed_xy)
    out = []
    for mirror in (False, True):
        for extra in (0.0, 180.0):
            out.append(from_rigid(
                angle_deg=(af - am) + extra, translation_um=tuple(cf - cm),
                scale=1.0, mirror=mirror, center_um=tuple(cm)))
    return out


def quick_score(tf: Transform, moving_xy: np.ndarray, fixed_xy: np.ndarray,
                *, radius_um: float = 25.0, n: int = 20000, seed: int = 0) -> float:
    """Cheap seed ranking: fraction of a moving subsample landing within ``radius_um``."""
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(moving_xy), size=min(n, len(moving_xy)), replace=False)
    d, _ = cKDTree(fixed_xy).query(tf.apply(moving_xy[idx]), k=1, workers=-1)
    return float((d <= radius_um).mean())


def fit(moving_xy: np.ndarray, fixed_xy: np.ndarray, *, verbose: bool = True) -> tuple[Transform, dict]:
    """Moment seed -> annealed ICP. Returns the transform and a provenance dict."""
    cands = seed_candidates(moving_xy, fixed_xy)
    scores = [quick_score(t, moving_xy, fixed_xy) for t in cands]
    best = int(np.argmax(scores))
    tf = cands[best]
    prov = {
        "seed_scores_within_25um": [round(s, 4) for s in scores],
        "seed_chosen": best,
        "seed_mirrored": bool(best >= 2),
    }
    if verbose:
        print(f"[register] seed candidates {prov['seed_scores_within_25um']} -> #{best}"
              f"{' (MIRRORED)' if prov['seed_mirrored'] else ''}", flush=True)

    for gate in ICP_SCHEDULE_UM:
        t0 = time.time()
        tf = icp(moving_xy, fixed_xy, initial=tf, max_dist_um=gate)
        if verbose:
            print(f"[register] icp gate {gate:>6.1f} um  ({time.time()-t0:.0f}s)", flush=True)
    prov["icp_schedule_um"] = list(ICP_SCHEDULE_UM)
    a = tf.affine
    lin = a[:, :2]
    prov["scale"] = float(np.sqrt(abs(np.linalg.det(lin))))
    prov["rotation_deg"] = float(np.degrees(np.arctan2(lin[1, 0], lin[0, 0])))
    prov["translation_um"] = [float(a[0, 2]), float(a[1, 2])]
    prov["mirrored"] = bool(np.linalg.det(lin) < 0)
    return tf, prov


def held_out(moving_xy: np.ndarray, fixed_xy: np.ndarray, seed_tf: Transform) -> dict:
    """Spatial split-half: fit on each side of the median x, score on the other.

    An in-sample residual is the thing the fit minimised, so it cannot say whether the fit
    generalises. Splitting on geometry rather than at random keeps the two halves genuinely
    independent -- a random split leaves each half surrounded by the other's cells.
    """
    out = {}
    cut = float(np.median(moving_xy[:, 0]))
    halves = {"low": moving_xy[:, 0] <= cut, "high": moving_xy[:, 0] > cut}
    fits = {}
    for name, mask in halves.items():
        if mask.sum() < 1000:
            return {"error": f"half {name!r} has too few points ({int(mask.sum())})"}
        fits[name] = icp(moving_xy[mask], fixed_xy, initial=seed_tf, max_dist_um=15.0)
    for name, mask in halves.items():
        other = "high" if name == "low" else "low"
        s = refine_mod.score_transform(fits[other], moving_xy[mask], fixed_xy)
        out[f"fit_{other}_score_{name}"] = s.get("pair_distance_um", {}).get("p50")
    la, lb = fits["low"].affine, fits["high"].affine
    out["translation_disagreement_um"] = float(np.linalg.norm(la[:, 2] - lb[:, 2]))
    sa = np.sqrt(abs(np.linalg.det(la[:, :2])))
    sb = np.sqrt(abs(np.linalg.det(lb[:, :2])))
    out["scale_disagreement_ppm"] = float(abs(sa - sb) / max(sa, sb) * 1e6)
    return out


def xenium_nucleus_centroids(bundle: Path) -> np.ndarray:
    """One centroid per nucleus label from ``nucleus_boundaries.parquet`` (already microns)."""
    nb = pd.read_parquet(Path(bundle) / "nucleus_boundaries.parquet")
    key = "label_id" if "label_id" in nb.columns else "cell_id"
    g = nb.groupby(key)[["vertex_x", "vertex_y"]].mean()
    return g.to_numpy(dtype=float)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pheno", required=True, type=Path, help="cells_pheno.parquet from read_pheno_tsv")
    ap.add_argument("--bundle", required=True, type=Path, help="Xenium bundle directory")
    ap.add_argument("--out", required=True, type=Path, help="directory for transform + report")
    ap.add_argument("--stem", default="sameslide")
    ap.add_argument("--pheno-nuclei", type=Path, default=None,
                    help="optional cell GeoJSON; when given, PhenoCycler NUCLEUS centroids "
                         "are used instead of cell centroids")
    ap.add_argument("--subsample", type=int, default=0, help="0 = fit on all points")
    args = ap.parse_args(argv)

    pheno = pd.read_parquet(args.pheno, columns=["x_um", "y_um"])
    fixed_xy = pheno.to_numpy(dtype=float)
    moving_xy = xenium_nucleus_centroids(args.bundle)
    print(f"[register] moving (Xenium nuclei): {len(moving_xy):,}   "
          f"fixed (PhenoCycler): {len(fixed_xy):,}")

    nucleus_source = "xenium_nucleus_vs_phenocycler_CELL_centroid"
    if args.pheno_nuclei is None:
        print("[register] NOTE: PhenoCycler cell centroids are standing in for nucleus "
              "centroids. A nucleus sits off-centre in its cell, so this leaves a systematic "
              "offset the fit cannot remove. Pass --pheno-nuclei once the segmentations land.",
              flush=True)

    if args.subsample:
        rng = np.random.default_rng(0)
        moving_xy = moving_xy[rng.choice(len(moving_xy), min(args.subsample, len(moving_xy)), replace=False)]
        fixed_xy = fixed_xy[rng.choice(len(fixed_xy), min(args.subsample, len(fixed_xy)), replace=False)]

    tf, prov = fit(moving_xy, fixed_xy)
    score = refine_mod.score_transform(tf, moving_xy, fixed_xy)
    print(f"[register] pair p50 {score.get('pair_distance_um', {}).get('p50')} um   "
          f"matched_within {score.get('matched_within')}")

    field = None
    resid = score.get("residual_after_constant_um")
    if resid is not None and resid >= LOCAL_FIELD_MIN_UM:
        try:
            field = refine_mod.refine_local(tf, moving_xy, fixed_xy,
                                            pixel_um=LOCAL_FIELD_PIXEL_UM)
            fs = refine_mod.score_transform(field, moving_xy, fixed_xy)
            if (fs["pair_distance_um"]["p50"] < score["pair_distance_um"]["p50"] - 0.05):
                tf, score = field, fs
                prov["local_field"] = {"pixel_um": LOCAL_FIELD_PIXEL_UM, "applied": True}
            else:
                prov["local_field"] = {"applied": False, "reason": "no improvement beyond tolerance"}
        except Exception as exc:                                   # noqa: BLE001
            prov["local_field"] = {"applied": False, "reason": f"{type(exc).__name__}: {exc}"}
    else:
        prov["local_field"] = {"applied": False,
                               "reason": f"residual_after_constant_um {resid} < {LOCAL_FIELD_MIN_UM}"}

    args.out.mkdir(parents=True, exist_ok=True)
    tf.save(args.out, stem=args.stem)
    report = {
        "direction": "xenium_moving -> phenocycler_fixed",
        "nucleus_source": nucleus_source if args.pheno_nuclei is None else "both_nuclei",
        "gate": "nucleus-pair statistics vs a 25 um placement null; align.py's DAPI-window "
                "NCC gate does NOT apply because this lane reads no images",
        "provenance": prov,
        "score": score,
        "held_out": held_out(moving_xy, fixed_xy, tf),
    }
    (args.out / f"{args.stem}_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"[register] wrote {args.out / (args.stem + '.json')}")
    print(f"[register] held_out: {report['held_out']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
