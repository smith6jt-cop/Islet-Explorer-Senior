#!/usr/bin/env python3
"""Bring the per-cell crossed fractions up to the islet level.

``build_crossed_fraction.py`` measures, for every islet-associated cell, how much of its own
segmented area a vessel, lymphatic or nerve object covers. Per islet and region that becomes two
different questions, and both are written because they answer different things:

    crossed_any_<structure>   the share of the region's cells that any of that structure touches
    crossed_mean_<structure>  the mean covered fraction over the cells it touches

The first is reach, the second is depth. Reporting only the first is what the old per-structure
aggregate did, and it made a structure that grazes half the cells look like one that engulfs them:
vessels touch ~25-54 % of islet-associated cells but cover a median of only ~8 % of a touched
cell's area.

The per-cell file stays the record -- it joins to anything by ``object_id``. This is the rollup the
app can key by ``(Case ID, islet_key)``.

    python core/phenocycler/scripts/typing/summarise_crossed_fraction.py
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"
CROSSED = S.OUT / "crossed"
STRUCTURES = ("vessel", "lymphatic", "nerve")

#: Below this many cells a regional mean is one or two cells' geometry, not a summary.
MIN_CELLS = 5


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--no-merge", action="store_true",
                        help="accepted for old commands; this script never merges any more")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    frames = []
    for donor in donors:
        path = CROSSED / f"crossed_{donor}.parquet"
        if not path.exists():
            print(f"  {donor}: no crossed file -- skipped")
            continue
        crossed = pd.read_parquet(path)
        cells = pq.read_table(
            S.OUT / f"cells/donor_id={donor}/data.parquet",
            columns=["object_id", "islet_key", "cell_region_20um"]).to_pandas()
        merged = crossed.merge(cells, on="object_id", how="left", validate="one_to_one")
        merged["region"] = merged["cell_region_20um"].map({"core": "core", "expansion": "peri"})
        merged = merged[merged["islet_key"].notna() & merged["region"].notna()]
        merged["case_id"] = int(donor)
        frames.append(merged)
        print(f"  {donor}: {len(merged):>7,} cells", flush=True)

    if not frames:
        print("nothing to summarise")
        return 1
    cells = pd.concat(frames, ignore_index=True)

    out = []
    for region, suffix in (("core", "core"), ("peri", "peri")):
        here = cells[cells["region"] == region]
        grouped = here.groupby(["case_id", "islet_key"])
        summary = grouped.size().rename(f"crossed_n_{suffix}").to_frame()
        for structure in STRUCTURES:
            column = f"crossed_frac_{structure}"
            touched = here[column] > 0
            summary[f"crossed_any_{structure}_{suffix}"] = (
                here.assign(_t=touched).groupby(["case_id", "islet_key"])["_t"].mean())
            # Depth is over the touched cells only; an islet none of them touch has no depth.
            depth = here[touched].groupby(["case_id", "islet_key"])[column].mean()
            summary[f"crossed_mean_{structure}_{suffix}"] = depth
        thin = summary[f"crossed_n_{suffix}"] < MIN_CELLS
        value_cols = [c for c in summary.columns if c.startswith("crossed_") and "_n_" not in c]
        summary.loc[thin, value_cols] = np.nan
        print(f"  {region:>5}: {len(summary):,} islets, "
              f"{int(thin.sum()):,} below the {MIN_CELLS}-cell floor (Unavailable)")
        out.append(summary)

    summary = out[0].join(out[1], how="outer").reset_index()
    summary.to_parquet(CROSSED / "crossed_islet_summary.parquet", index=False)
    print(f"\n{len(summary):,} islets -> {CROSSED / 'crossed_islet_summary.parquet'}")
    for structure in STRUCTURES:
        reach = summary[f"crossed_any_{structure}_peri"].median()
        depth = summary[f"crossed_mean_{structure}_peri"].median()
        print(f"  peri {structure:<10s} median reach {reach:.3f}  median depth {depth:.3f}")

    # The merge into islets.parquet is build_parquet_views.py's (its ROLLUPS table reads the
    # summary written above).  This script writes one file.
    return 0
    islets = pd.read_parquet(APP / "parquet/islets.parquet")
    stale = [c for c in islets.columns if c.startswith("crossed_")]
    if stale:
        islets = islets.drop(columns=stale)
    merged = islets.merge(summary, on=["case_id", "islet_key"], how="left", validate="one_to_one")
    if len(merged) != len(islets):
        raise ValueError("the merge changed the islet count")
    merged.to_parquet(APP / "parquet/islets.parquet", index=False)
    covered = merged["crossed_any_vessel_peri"].notna().sum()
    print(f"  merged into islets.parquet: {covered:,} of {len(merged):,} islets carry a peri "
          f"vessel reach")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
