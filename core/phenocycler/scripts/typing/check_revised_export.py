#!/usr/bin/env python3
"""Check one re-exported QuPath cell table before anything downstream reads it.

The batch-2 export carried a second detection pass inside the islet annotations, so every
islet-core cell existed roughly twice (memory: batch2-export-defects, plan correction A). The
re-export is meant to fix that at source. This measures whether it did, and whether the
object ids survived -- which is what lets every join stay on the QuPath UUID with no centroid
matching anywhere.

Three measurements, and nothing is inferred from any of them about biology:

  1. SHAPE            rows and columns, against the other revised exports.
  2. DUPLICATE RATE   a duplicated detection is the SAME cell segmented twice, so the test
                      is same image AND same centroid AND same geometry -- nucleus area, cell
                      area and max diameter. Reported separately for islet parents and
                      everything else, because the defect was specific to the islet
                      annotations.

                      CENTROID ALONE IS NOT THE TEST, and an earlier version of this script
                      used it. The export prints centroids to 0.1 um, so in dense tissue two
                      genuinely different cells collide on the printed value now and then --
                      measured on 6534's re-export: 79 keys, 158 rows, and ZERO of them
                      agreed on any geometry column (nucleus areas 45.97 vs 55.24, 17.74 vs
                      18.74, and so on). Reporting those 158 as "0.06% duplicates" was
                      reporting the export's coordinate precision, not a defect. The
                      collision count is still printed, labelled for what it is.
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
    same_cell: dict[bool, Counter] = {True: Counter(), False: Counter()}
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
        gg = [index[c] for c in GEOMETRY]
        for row in reader:
            n += 1
            parent = row[gp]
            parents[parent] += 1
            in_islet = "islet" in parent.lower()
            place = (row[gm], row[gx], row[gy])
            seen[in_islet][place] += 1
            same_cell[in_islet][place + tuple(row[i] for i in gg)] += 1
            ids.add(row[gi])
            if n % 5_000_000 == 0:
                print(f"    {n:,} rows  {time.time() - started:.0f}s", flush=True)

    print(f"  {n:,} rows, {len(ids):,} distinct object ids  ({time.time() - started:.0f}s)\n")
    if len(ids) != n:
        print(f"  !! {n - len(ids):,} rows share an object id with another row\n")

    print("  duplicate detections (same image + centroid + nucleus area + cell area + max diameter)")
    for in_islet, label in ((True, "islet parents"), (False, "everything else")):
        counts = same_cell[in_islet]
        total = sum(counts.values())
        if not total:
            print(f"    {label:<18} no rows")
            continue
        dup_rows = sum(c for c in counts.values() if c > 1)
        verdict = "CLEAN" if dup_rows == 0 else "DEFECT PRESENT"
        print(f"    {label:<18} {total:>11,} rows  {dup_rows:>8,} duplicated  "
              f"({100 * dup_rows / total:5.2f}%)  {verdict}")
    print("    for reference, the FIRST batch-2 export measured, inside islets:")
    print("      6450 57.0%  6505 51.3%  6523 58.2%  6534 58.5%")
    print("      6566 50.6%  6591 58.7%  6623 47.8%   -- and 1.4-5.2% outside them")

    print("\n  printed-centroid collisions -- NOT duplicates, reported so the number is not")
    print("  mistaken for one. Centroids print to 0.1 um; distinct cells sometimes share a value.")
    for in_islet, label in ((True, "islet parents"), (False, "everything else")):
        counts = seen[in_islet]
        total = sum(counts.values())
        if not total:
            continue
        hits = sum(c for c in counts.values() if c > 1)
        print(f"    {label:<18} {hits:>8,} rows share a printed centroid with another")

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
