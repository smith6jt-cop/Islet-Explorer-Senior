#!/usr/bin/env python
"""Turn VALIS's warped centroids into the lane's Xenium -> IF ``Transform``.

    cd core/phenocycler
    ../../.venv/bin/python scripts/restain/valis_transform.py --donors 6442

VALIS hands back where each Xenium cell centroid lands in the IF frame. That is a point
mapping, not a transform, so this fits one and -- more importantly -- **reports how well an
affine fits it**. The residual field measured independently by ``run_align.py`` was affine at
R^2 = 0.9997, so a sub-micron fit here is two methods agreeing; a large non-affine part would
mean they do not, and that is worth knowing before anything is built on either.

If the affine does not fit, the mapping is kept as a displacement field on a coarse grid rather
than forced through six parameters -- ``Transform`` already carries and applies one.
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

from run_register import OUT                                      # noqa: E402

from phenocycler.integration.rasterize import GridSpec            # noqa: E402
from phenocycler.integration.transform import (Transform,         # noqa: E402
                                               fit_affine_lstsq)

#: Above this the mapping is not affine and a field is kept instead. A cell is ~9 um, and the
#: pass threshold for per-cell work is ~1.9 um, so half a micron is comfortably inside both.
AFFINE_RESIDUAL_MAX_UM = 0.5
#: Spacing of the fallback displacement field. 250 um is far finer than any deformation a
#: fitted affine could have missed and still cheap to carry.
FIELD_PIXEL_UM = 250.0


def _field_from_residual(src, dst, affine, *, pixel_um=FIELD_PIXEL_UM):
    """Residual of the affine fit, binned onto a grid, as a (2, H, W) displacement in microns.

    Two things a naive binning gets wrong, both measured on donor 6442:

    * **A bin with no points is UNKNOWN, not zero.** 32 % of bins over the section's bounding
      box contain no cell, and leaving those at 0.0 makes the bilinear sampler ramp from a real
      displacement down to nothing across the tissue edge -- a fabricated gradient exactly where
      the geometry is least constrained. Empty bins are filled from the nearest populated one.
    * **A bin's value must be robust.** 0.33 % of VALIS's warped points land more than 50 um
      from the affine prediction (max 250 um), enough to drag a bin's MEAN hundreds of microns:
      the mean-binned field peaked at 245.9 um where the per-point residual median was 1.1 um.
      A per-bin median ignores them.
    """
    from scipy.ndimage import distance_transform_edt
    from scipy.stats import binned_statistic_2d

    pred = src @ np.asarray(affine)[:, :2].T + np.asarray(affine)[:, 2]
    res = dst - pred
    x0, y0 = src.min(axis=0) - pixel_um
    x1, y1 = src.max(axis=0) + pixel_um
    w = int(np.ceil((x1 - x0) / pixel_um)) + 1
    h = int(np.ceil((y1 - y0) / pixel_um)) + 1
    grid = GridSpec(float(x0), float(y0), float(pixel_um), h, w)
    # HALF A BIN. ``GridSpec`` maps micron ``origin + k * pixel_um`` to index k, so index k is
    # the bin's CENTRE, not its left edge. Starting the edges at ``x0`` puts every bin's median
    # half a bin (125 um) from the node the sampler reads it at, which no amount of refitting
    # upstream can correct because it happens after the transform.
    xedges = x0 - 0.5 * pixel_um + np.arange(w + 1) * pixel_um
    yedges = y0 - 0.5 * pixel_um + np.arange(h + 1) * pixel_um

    field = np.zeros((2, h, w), dtype=np.float64)
    filled = None
    for k in (0, 1):
        stat, *_ = binned_statistic_2d(src[:, 1], src[:, 0], res[:, k],
                                       statistic="median", bins=[yedges, xedges])
        empty = ~np.isfinite(stat)
        if filled is None:
            filled = ~empty
        if empty.all():
            continue
        # nearest populated bin, so the field never ramps to a fabricated zero
        _d, idx = distance_transform_edt(empty, return_distances=True, return_indices=True)
        field[k] = np.where(empty, stat[tuple(idx)], stat)
    return field, grid, res, filled


def transform_one(donor_dir: Path) -> dict:
    donor_dir = Path(donor_dir)
    src = np.load(donor_dir / "xenium_centroids_um.npy")
    dst_p = donor_dir / "warped_if_um.npy"
    if not dst_p.exists():
        raise FileNotFoundError(f"{dst_p} is absent -- run valis_register.py first")
    dst = np.load(dst_p)
    meta = json.loads((donor_dir / "export.json").read_text())

    ok = np.isfinite(src).all(axis=1) & np.isfinite(dst).all(axis=1)
    src, dst = src[ok], dst[ok]

    tf = fit_affine_lstsq(src, dst)
    pred = src @ np.asarray(tf.affine)[:, :2].T + np.asarray(tf.affine)[:, 2]
    resid = np.linalg.norm(dst - pred, axis=1)
    med, p90, p99 = (float(np.median(resid)), float(np.percentile(resid, 90)),
                     float(np.percentile(resid, 99)))

    rec = {"donor_id": meta["donor_id"], "serial": meta["serial"], "roi": meta["roi"],
           "n_points": int(len(src)), "n_dropped_nonfinite": int((~ok).sum()),
           "affine": np.asarray(tf.affine).tolist(),
           "affine_residual_um": {"median": med, "p90": p90, "p99": p99,
                                  "max_allowed": AFFINE_RESIDUAL_MAX_UM},
           "is_affine": bool(med <= AFFINE_RESIDUAL_MAX_UM)}

    if not rec["is_affine"]:
        field, grid, _res, filled = _field_from_residual(src, dst, tf.affine)
        np.savez_compressed(donor_dir / "valis_displacement.npz", field=field,
                            **{k: np.asarray(v) for k, v in grid.as_dict().items()})
        mag = np.hypot(field[0], field[1])
        rec["displacement"] = {"file": str(donor_dir / "valis_displacement.npz"),
                               "grid": grid.as_dict(), "pixel_um": FIELD_PIXEL_UM,
                               "median_um": float(np.median(mag)),
                               "p99_um": float(np.percentile(mag, 99)),
                               "max_um": float(mag.max()),
                               "bins_populated": int(filled.sum()),
                               "bins_total": int(filled.size),
                               "bins_filled_from_neighbour": int((~filled).sum())}
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=OUT / "valis")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--roi", default="panc")
    args = ap.parse_args(argv)

    dirs = sorted(d for d in args.root.glob("donor_id=*") if d.is_dir())
    if args.donors:
        want = {str(x) for x in args.donors}
        dirs = [d for d in dirs if d.name.split("=", 1)[1] in want]
    rows = []
    for d in dirs:
        if not (d / "warped_if_um.npy").exists():
            print(f"[{d.name}] no warped points yet -- skipping", flush=True)
            continue
        rec = transform_one(d)
        r = rec["affine_residual_um"]
        print(f"[{rec['donor_id']}] affine fit to VALIS on {rec['n_points']:,} centroids: "
              f"residual median {r['median']:.3f} um, p90 {r['p90']:.3f}, p99 {r['p99']:.3f} "
              f"-> {'AFFINE' if rec['is_affine'] else 'keeping a displacement field'}", flush=True)
        rows.append(rec)
    if not rows:
        print("nothing to do", file=sys.stderr)
        return 1
    args.out.mkdir(parents=True, exist_ok=True)
    p = args.out / f"valis_registration_{args.roi}.json"
    p.write_text(json.dumps(rows, indent=2))
    print(f"-> {p}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
