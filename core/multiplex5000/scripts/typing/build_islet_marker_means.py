#!/usr/bin/env python3
"""Per-islet MEAN MARKER VALUES -- the thing Phase 5 promised and never built.

The Trajectory tab's Multi-Feature Heatmap is meant to show how MARKERS move along pseudotime.
It was showing the eight cell classes instead, and faithfully so: `adata_ins_root.h5ad` is
(7302 islets x 8 cell classes), its `X` IS the composition, and there is no marker anywhere in
the file -- no layer, and `obsm` holds only the ALR composition, the diffmap and the UMAP. The
heatmap read `var_names` and got what was there. `aggregate_islets.py` never computed marker
means, so nothing downstream could have shown them.

This builds them from the per-cell UniFORM-aligned values -- the REDSEA-corrected, aligned
intensities that are the authoritative one-value-per-marker -- averaged within each islet.

THREE CHOICES THAT ARE NOT ARBITRARY:

  * ELIGIBLE CELLS ONLY. `qc_analysis_eligible` excludes detections with no nucleus, which are
    not known to be cells at all, and averaging their intensities into an islet's mean would
    import measurements from objects the rest of the analysis refuses to count.

  * CORE AND PERI ARE SEPARATE COLUMNS, never pooled. An islet's core mean and the mean of
    the 20 um ring around it are different quantities -- that is the entire point of the region
    split -- and a pooled mean would move with how much ring a given islet happens to have.

  * ISLET MEMBERSHIP COMES FROM `typing/cells`, JOINED BY `object_id` -- never from the aligned
    table's own `islet_num` / `cell_region`. Measured 2026-09-17: on the seven re-exported
    batch-2 donors the aligned table still carries the ORIGINAL export's islet numbering, and it
    does not survive the re-export. On the cells the two tables share and both label, the two
    numberings agree on 0.00% (6623, 6450) and 0.15% (6505) -- against 99.94-99.97% on the
    batch-1 donors, where nothing was renumbered. Deriving `islet_key` from the aligned column
    therefore attached every batch-2 marker mean to a DIFFERENT islet than the one that key
    names in `islet_geometry`, `islet_celltype_proportions` and `neighborhood_metrics`: a
    right-looking islet served the wrong cells' intensities. It also silently dropped most of
    those donors' islets (6505: 646 of 3,117) because the stale numbers had nothing to match.
    The membership in `typing/cells` is the geometric one every other app artifact joins on, so
    it is the one this reads.

A marker absent from a donor's panel yields NaN, not zero: unavailable is not absence. A cell
with no aligned row has no measured intensity to contribute and is counted, never imputed.

    python core/multiplex5000/scripts/typing/build_islet_marker_means.py [--donors 115]
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time

import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

OUT = S.OUT / "app"
META = {"donor_id", "object_id", "qc_analysis_eligible", "X_centroid", "Y_centroid",
        "cell_region", "islet_num"}


#: The per-cell table that owns islet membership. `cell_region_20um` is the app's "expansion";
#: the column names this build writes keep the app's older spelling, `_peri`.
REGION_COLUMNS = (("core", ""), ("expansion", "_peri"))


def _membership(donor: str) -> pd.DataFrame | None:
    """object_id -> the islet and region every other app artifact puts that cell in.

    `Excluded` -- exocrine inside an islet border -- is dropped here for the same reason
    `aggregate_islets.py` drops it from the composition denominator: it is out of analysis, so
    its intensities do not belong in an islet's mean either. Not cosmetic. Measured 2026-09-17,
    it is 2.7% of 115's islet cells but 21.9% of 6505's and 22.7% of 6534's, so keeping it would
    have pulled those two donors' means toward exocrine while their compositions ignored the same
    cells -- two tables on one screen answering the same question over different denominators.
    It also produced an islet the composition has no row for at all: 6539_Islet_536 holds exactly
    one eligible cell and that cell is `Excluded`.
    """
    path = S.OUT / f"cells/donor_id={donor}/data.parquet"
    if not path.exists():
        return None
    cells = pq.read_table(path, columns=["object_id", "islet_key", "combined_islet_id",
                                         "cell_region_20um", "cell_class",
                                         "qc_analysis_eligible"]).to_pandas()
    cells = cells[cells["qc_analysis_eligible"].fillna(False).astype(bool)]
    cells = cells[cells["cell_class"] != "Excluded"]
    regions = [region for region, _ in REGION_COLUMNS]
    return cells[cells["islet_key"].notna() & cells["cell_region_20um"].isin(regions)]


def donor_means(donor: str) -> tuple[pd.DataFrame, dict[str, int]] | None:
    path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    if not path.exists():
        return None
    members = _membership(donor)
    if members is None or members.empty:
        return None

    df = pq.read_table(path).to_pandas()
    markers = [c for c in df.columns if c not in META]
    df = df[df["qc_analysis_eligible"].fillna(False).astype(bool)]

    # An inner join, so an islet cell with no aligned row contributes nothing rather than a zero.
    # Both sides are unique on object_id; validate says so rather than trusting it.
    merged = members.merge(df[["object_id"] + markers], on="object_id", how="inner",
                           validate="one_to_one")
    unmeasured = len(members) - len(merged)
    if merged.empty:
        return None

    out = []
    for region, suffix in REGION_COLUMNS:
        part = merged[merged["cell_region_20um"] == region]
        if part.empty:
            continue
        grouped = part.groupby("islet_key")[markers].mean()
        grouped.columns = [f"marker_mean{suffix}_{m}" for m in markers]
        out.append(grouped)
    if len(out) != 2:
        raise RuntimeError(
            f"{donor}: expected core AND peri means, built {len(out)}. A region label that "
            f"matches nothing yields a silently half-built table, which is how the first run "
            f"produced core-only columns without erroring.")
    frame = pd.concat(out, axis=1)
    # The combined id is carried from the per-cell table, not re-spelled here, so this artifact
    # cannot drift from the key the other app tables are joined on.
    keys = merged.drop_duplicates("islet_key").set_index("islet_key")["combined_islet_id"]
    frame.insert(0, "combined_islet_id", keys.reindex(frame.index).to_numpy())
    frame.insert(1, "case_id", donor)
    stats = {"islets": len(frame), "cells": len(merged), "unmeasured_cells": unmeasured,
             "markers": len(markers)}
    return frame.reset_index(), stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--include-marginal", action="store_true")
    args = ap.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    frames = []
    for donor in donors:
        started = time.time()
        built = donor_means(donor)
        if built is None:
            print(f"  {donor}: no aligned per-cell values"); continue
        frame, stats = built
        frames.append(frame)
        print(f"  {donor}: {stats['islets']:>6,} islets | {stats['markers']} markers x "
              f"core+expansion | {stats['cells']:>8,} islet cells measured, "
              f"{stats['unmeasured_cells']:,} with no aligned row (unavailable) | "
              f"{time.time()-started:.0f}s", flush=True)
    if not frames:
        return 1
    everything = pd.concat(frames, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "islet_marker_means.parquet"
    everything.to_parquet(target, index=False)
    print(f"\n{len(everything):,} islets, {everything.shape[1]} columns -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
