#!/usr/bin/env python3
"""Check one re-exported QuPath cell table before anything downstream reads it.

The batch-2 export carried a second detection pass inside the islet annotations, so every
islet-core cell existed roughly twice (memory: batch2-export-defects, plan correction A). The
re-export is meant to fix that at source. This measures whether it did, and whether the
object ids survived -- which is what lets every join stay on the QuPath UUID with no centroid
matching anywhere.

Three measurements, and nothing is inferred from any of them about biology:

  1. SHAPE            rows and columns, against the other revised exports.
  2. DUPLICATE RATE   **DROPPED. The re-exports have ZERO duplicates and this check only ever
                      produced noise around that fact.**

                      QuPath removes them itself, in scripts/0b_duplicate_check_only.groovy:
                      centroids within 1.0 um AND cell areas within 15%, joined by union-find,
                      one survivor per group. That script IS the definition and it is what ran.
                      There is nothing for this file to second-guess.

                      Two failed attempts, recorded so neither is repeated:
                        * Exact string equality on centroid + areas. Returned 0 on the ORIGINAL
                          6534 export, which carries 36,596 known pairs -- it cannot see the
                          defect at all, so its "CLEAN" meant nothing.
                        * 0b's own tolerance, re-implemented here. Fires correctly on the
                          original (41.73% inside islets, 0.00% outside) but reports 0.03-0.11%
                          on the revised exports -- and those are NOT duplicates. They sit at a
                          median 0.990 um against the 1.0 um ceiling, where a real duplicate sits
                          at ~0.1 um, and their nuclei disagree outright. The export prints
                          centroids to 0.1 um, so a true 1.02 um separation rounds into the box.
                          A boundary artefact of print precision, reported as a defect count.

                      To verify a re-export, read 0b's own console output, or ingest and read the
                      geometry QC's duplicate_group -- rasterised overlap, on full-precision
                      geometry, which is what the detection actually rests on.

  3. ID CONTINUITY    share of the export's Object IDs that are already in the content-
                      addressed run's per-cell table for that donor. This does NOT need to be
                      100%: a re-run detection legitimately creates new ids. It is reported so
                      that a re-export which silently regenerated EVERY id -- and would
                      therefore break every existing join -- is visible immediately.

                      It also reports, for the ids the OLD run put in a duplicate group, which
                      member of each pair the re-export kept. This is the trap: QuPath's
                      de-duplication and the geometry QC's tie-break are unrelated, so the
                      kept member is about 50/50 (6534: 18,179 old survivors and 18,453 old
                      losers). The old run marked every loser `analysis_eligible = FALSE`, so
                      JOINING THE NEW EXPORT TO THE OLD RUN'S QC FLAGS WOULD SILENTLY DROP
                      HALF THE ISLET CELLS. The new export needs its own ingest and its own
                      geometry QC; the old flags do not transfer.

    python core/phenocycler/scripts/typing/check_revised_export.py data/RevisedQuPathExports/6534.tsv
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import subpanel_io as S

WANTED = ["Image", "Object ID", "Parent", "Centroid X µm", "Centroid Y µm"]
#: What makes two rows the same DETECTION rather than two cells at the same printed spot.
GEOMETRY = ["Nucleus: Area µm^2", "Cell: Area µm^2", "Cell: Max diameter µm"]


def duplicate_rows(records, radius_um: float = 1.0, area_tol: float = 0.15) -> int:
    """Rows sitting in a duplicate group, by 0b_duplicate_check_only.groovy's rule.

    ``records`` is an iterable of (x, y, area). Buckets on a ``radius_um`` grid and compares
    each cell only against its own bucket and the eight around it -- the same trick the Groovy
    script uses, and the reason this is linear rather than quadratic. Union-find then joins
    transitively, so a chain of three near-coincident detections counts as ONE group with two
    redundant members, exactly as ``removeObjects`` would leave it.
    """
    import math

    pts = [(x, y, a) for x, y, a in records if x is not None and y is not None and a is not None]
    buckets: dict[tuple[int, int], list[int]] = {}
    for i, (x, y, _a) in enumerate(pts):
        buckets.setdefault((int(x // radius_um), int(y // radius_um)), []).append(i)

    parent = list(range(len(pts)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    for (bx, by), members in buckets.items():
        near = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                near.extend(buckets.get((bx + dx, by + dy), ()))
        for i in members:
            xi, yi, ai = pts[i]
            for j in near:
                if j <= i:
                    continue
                xj, yj, aj = pts[j]
                if math.hypot(xi - xj, yi - yj) > radius_um:
                    continue
                if abs(ai - aj) / max(ai, aj, 1e-9) < area_tol:
                    union(i, j)

    sizes: dict[int, int] = {}
    for i in range(len(pts)):
        r = find(i)
        sizes[r] = sizes.get(r, 0) + 1
    return sum(n for n in sizes.values() if n > 1)


def run_duplicate_sets(donor: str):
    """(survivors, losers) as the OLD run's geometry QC labelled them, or (None, None)."""
    import pyarrow.parquet as pq
    part = S.RUN / f"geometry_qc/donor_id={donor}"
    files = sorted(part.glob("*.parquet")) if part.is_dir() else []
    if not files:
        return None, None
    import pandas as pd
    frames = [pq.read_table(f, columns=["object_id", "duplicate_group", "duplicate_survivor"]
                            ).to_pandas() for f in files]
    g = pd.concat(frames, ignore_index=True)
    g = g[g["duplicate_group"].notna()]
    if g.empty:
        return set(), set()
    keep = g["duplicate_survivor"].fillna(False).astype(bool)
    return set(g.loc[keep, "object_id"]), set(g.loc[~keep, "object_id"])


def run_ids(donor: str) -> set[str] | None:
    """Object ids already in the run for this donor, or None if it has no partition."""
    import pyarrow.parquet as pq
    # The run writes `data_0.parquet`, not `data.parquet`; glob rather than
    # name it, because a fixed name silently reported "no partition" for a
    # donor that has one, which reads as a missing donor rather than a typo.
    part = S.RUN / f"cells/donor_id={donor}"
    files = sorted(part.glob("*.parquet")) if part.is_dir() else []
    if not files:
        return None
    ids = set()
    for f in files:
        ids |= set(pq.read_table(f, columns=["object_id"]).column("object_id").to_pylist())
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export", type=pathlib.Path)
    ap.add_argument("--donor", default=None, help="default: the file stem")
    args = ap.parse_args()

    donor = args.donor or args.export.stem
    started = time.time()
    print(f"{args.export}  ({args.export.stat().st_size / 1e9:.2f} GB)  donor {donor}\n")

    delim = "\t" if args.export.suffix == ".tsv" else ","
    seen: dict[bool, Counter] = {True: Counter(), False: Counter()}
    points: dict[tuple[str, bool], list] = defaultdict(list)
    parents = Counter()
    ids: set[str] = set()
    n = 0
    with args.export.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader)
        index = {name: i for i, name in enumerate(header)}
        missing = [c for c in WANTED if c not in index]
        if missing:
            print(f"  MISSING COLUMNS: {missing}")
            return 2
        print(f"  {len(header)} columns")
        gx, gy = index["Centroid X µm"], index["Centroid Y µm"]
        gp, gi, gm = index["Parent"], index["Object ID"], index["Image"]
        ga = index["Cell: Area µm^2"]
        for row in reader:
            n += 1
            parent = row[gp]
            parents[parent] += 1
            in_islet = "islet" in parent.lower()
            seen[in_islet][(row[gm], row[gx], row[gy])] += 1
            # Grouped by IMAGE: two cells in different images cannot be duplicates
            # of each other however close their coordinates look.
            points[(row[gm], in_islet)].append((_f(row[gx]), _f(row[gy]), _f(row[ga])))
            ids.add(row[gi])
            if n % 5_000_000 == 0:
                print(f"    {n:,} rows  {time.time() - started:.0f}s", flush=True)

    print(f"  {n:,} rows, {len(ids):,} distinct object ids  ({time.time() - started:.0f}s)\n")
    if len(ids) != n:
        print(f"  !! {n - len(ids):,} rows share an object id with another row\n")

    print("\n  parents, by share")
    for parent, c in parents.most_common(6):
        print(f"    {100 * c / n:5.1f}%  {c:>10,}  {parent[:64] or '(blank)'}")

    known = run_ids(donor)
    if known is None:
        print(f"\n  no run partition for {donor}; id continuity not measurable here")
    else:
        both = len(ids & known)
        print(f"\n  id continuity against run {S.RUN.name}")
        print(f"    {both:>11,} of {len(ids):,} export ids are in the run  ({100 * both / len(ids):.1f}%)")
        print(f"    {len(known) - both:>11,} of {len(known):,} run ids are not in the export")
        surv, lose = run_duplicate_sets(donor)
        if surv:
            print(f"\n  of the old run's {len(surv):,} duplicate PAIRS, the re-export kept")
            print(f"    {len(surv & ids):>8,} old survivors  and  {len(lose & ids):>8,} old losers")
            print("    -- roughly one per pair, but WHICH one is arbitrary: QuPath's")
            print("    de-duplication and the geometry QC's raster tie-break are unrelated.")
            print("    The old run marks every loser analysis_eligible=FALSE, so joining this")
            print("    export to the OLD QC flags would drop about half the islet cells.")
            print("    Re-ingest and recompute geometry QC; do not carry the old flags over.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
