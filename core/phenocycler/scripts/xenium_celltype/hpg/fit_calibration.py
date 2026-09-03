#!/usr/bin/env python
"""Turn the calibration ladder into the two numbers that size the production array.

T(M) = fixed + M * per_iteration, at fixed n and fixed cores. `fixed` is the kNN plus the
dissimilarity reweighting, which every shard of a donor-arm pays again and which does use
its cores. `per_iteration` is the Monte-Carlo cost, which shards perfectly and which cannot
use more than about sqrt(m)/2 cores inside one process.

The two decide the shape of the whole run:

  * the shard wall clock is  fixed + (m / occupancy) * per_iteration
  * the redundant cost of sharding is  (n_shards - 1) * fixed  per donor-arm
  * so m is chosen as the smallest value whose wall fits comfortably inside the QoS limit,
    and cpus-per-task as floor(sqrt(m)/2), which is the most KODAMA's chunking can use.

This prints the table rather than choosing for you, because the QoS wall limit is a fact
about the cluster on the day and `sacctmgr show qos` is the authority on it, not this file.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def fit(points: list[tuple[int, float]]) -> tuple[float, float]:
    n = len(points)
    sx = sum(m for m, _ in points); sy = sum(t for _, t in points)
    sxx = sum(m * m for m, _ in points); sxy = sum(m * t for m, t in points)
    denom = n * sxx - sx * sx
    if denom == 0:
        raise SystemExit("calibration needs at least two distinct M values")
    slope = (n * sxy - sx * sy) / denom
    return (sy - slope * sx) / n, slope


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--wall-hours", type=float, default=72.0,
                    help="the QoS MaxWall you intend to submit under")
    ap.add_argument("--M", type=int, default=None, help="default: the manifest's M")
    a = ap.parse_args(argv)

    man = json.loads((a.run / "work_manifest.json").read_text())
    M = a.M or int(man["kodama"]["M"])
    biggest = max(d["n_cells_expression"] for d in man["donors"])

    print(f"{'donor@L':>14} {'n_cells':>10} {'fixed_s':>10} {'per_iter_s':>11} {'cores':>6}")
    fits = {}
    for cal in sorted((a.run / "calibration").glob("donor_id=*/landmarks=*")):
        donor = f"{cal.parent.name.split('=')[1]}@L{cal.name.split('=')[1]}"
        pts, cores = [], None
        for j in sorted(cal.glob("cal_M*.json")):
            d = json.loads(j.read_text())
            pts.append((int(d["m"]), float(d["seconds"])))
            cores = int(d["n_cores"]); n_cells = int(d["n_cells"])
        if len(pts) < 2:
            continue
        fixed, per = fit(sorted(pts))
        fits[donor] = (n_cells, fixed, per, cores)
        print(f"{donor:>14} {n_cells:>10,} {fixed:>10.0f} {per:>11.0f} {cores:>6}")

    if not fits:
        raise SystemExit("no calibration points found under --run/calibration")

    # Extrapolate to the largest donor with the measured n^1.44 law, anchored on whichever
    # calibration donor is closest in size -- stated, not hidden, because the exponent came
    # from the workstation ladder and this is the one place it is being trusted.
    anchor = max(fits.values(), key=lambda v: v[0])
    n0, fixed0, per0, cores0 = anchor
    scale = (biggest / n0) ** 1.44
    print(f"\nextrapolated to the largest donor ({biggest:,} cells) with n^1.44 from n={n0:,}:")
    print(f"  fixed {fixed0*scale/3600:.1f} h | per iteration {per0*scale/3600:.2f} h "
          f"(at {cores0} cores)\n")

    print(f"{'m':>5} {'shards':>7} {'cores':>6} {'occ':>4} {'shard_wall_h':>13} "
          f"{'redundant_h':>12} {'fits':>5}")
    for m in (1, 2, 4, 8, 16, 25, 32, 48, 64, 100):
        if M % m:
            continue
        cores = max(1, int(math.isqrt(m // 4))) if m >= 4 else 1
        occ = min(max(1, min(16, math.ceil(m / (cores * 4)))), cores) if cores > 1 else 1
        wall = (fixed0 * scale + (m / occ) * per0 * scale) / 3600
        redundant = (M // m - 1) * fixed0 * scale / 3600
        print(f"{m:>5} {M//m:>7} {cores:>6} {occ:>4} {wall:>13.1f} {redundant:>12.1f} "
              f"{'yes' if wall < a.wall_hours else 'NO':>5}")
    print(f"\nColumns: occ = iterations KODAMA can actually run at once in one process; "
          f"redundant_h = fixed cost paid again by every extra shard of one donor-arm.\n"
          f"Pick the smallest m whose shard_wall_h clears {a.wall_hours:.0f} h with margin.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
