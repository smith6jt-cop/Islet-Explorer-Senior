#!/usr/bin/env python3
"""Bring the single-cell ordering up to the islet level, where the app already has views.

``trajectory_singlecell.R`` orders individual islet and peri-islet cells.  Every islet-level view
in the app -- the Plot tab's measures, Statistics, the Spatial cards -- is keyed by
``(Case ID, islet_key)``, so the ordering is only usable there once it is summarised per islet and
region.  This writes that summary and merges it into ``islets.parquet``.

Two things stated rather than assumed:

* **The summary covers only the cells that were ordered.**  The trajectory runs on a proportional
  per-donor subsample, so an islet's summary is over the sampled cells in it, and the count is
  carried beside every statistic (``sc_n_core``, ``sc_n_peri``).  An islet with no sampled cell in
  a region has no value there -- Unavailable, not zero.
* **The median, not the mean.**  Diffusion pseudotime is an ordering, not a quantity with a
  meaningful scale, so an islet's position on it is better described by the middle of its cells
  than by their average, which a single extreme cell moves.  The mean is carried too, for anyone
  who wants it.

    python core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import subpanel_io as S

APP = S.OUT / "app"
TRAJ = S.OUT / "trajectory"

#: Minimum sampled cells in a region before an islet gets a value there. Below this the median is
#: one or two cells' position, which is not a summary of anything.
MIN_SAMPLED_CELLS = 5


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-merge", action="store_true",
                        help="write the summary but leave islets.parquet alone")
    args = parser.parse_args()

    cells = pd.read_parquet(TRAJ / "singlecell_pseudotime.parquet")
    cells = cells[np.isfinite(cells["pseudotime"])]
    cells["case_id"] = cells["donor_id"].astype(int)
    print(f"{len(cells):,} ordered cells over {cells['case_id'].nunique()} donors")

    frames = []
    for region, suffix in (("core", "core"), ("peri", "peri")):
        here = cells[cells["region"] == region]
        grouped = here.groupby(["case_id", "islet_key"])["pseudotime"]
        summary = grouped.agg(["median", "mean", "size"]).rename(columns={
            "median": f"sc_pseudotime_{suffix}",
            "mean": f"sc_pseudotime_mean_{suffix}",
            "size": f"sc_n_{suffix}"})
        thin = summary[f"sc_n_{suffix}"] < MIN_SAMPLED_CELLS
        summary.loc[thin, [f"sc_pseudotime_{suffix}", f"sc_pseudotime_mean_{suffix}"]] = np.nan
        print(f"  {region:>5}: {len(summary):,} islets, "
              f"{int(thin.sum()):,} below the {MIN_SAMPLED_CELLS}-cell floor (Unavailable)")
        frames.append(summary)

    summary = frames[0].join(frames[1], how="outer").reset_index()
    # The peri-minus-core difference is the quantity the study question asks for: how the cells
    # around an islet sit on the ordering relative to the cells inside it.
    summary["sc_pseudotime_peri_minus_core"] = (summary["sc_pseudotime_peri"]
                                                - summary["sc_pseudotime_core"])
    summary.to_parquet(TRAJ / "singlecell_islet_summary.parquet", index=False)
    print(f"\n{len(summary):,} islets -> {TRAJ / 'singlecell_islet_summary.parquet'}")
    both = summary["sc_pseudotime_peri_minus_core"].notna()
    print(f"  {int(both.sum()):,} islets have both regions; median peri-core difference "
          f"{summary.loc[both, 'sc_pseudotime_peri_minus_core'].median():+.4g}")

    if args.no_merge:
        return 0

    islets = pd.read_parquet(APP / "parquet/islets.parquet")
    stale = [c for c in islets.columns if c.startswith("sc_pseudotime") or c.startswith("sc_n_")]
    if stale:
        islets = islets.drop(columns=stale)
    merged = islets.merge(summary, on=["case_id", "islet_key"], how="left", validate="one_to_one")
    if len(merged) != len(islets):
        raise ValueError("the merge changed the islet count")
    merged.to_parquet(APP / "parquet/islets.parquet", index=False)
    covered = merged["sc_pseudotime_core"].notna().sum()
    print(f"  merged into islets.parquet: {covered:,} of {len(merged):,} islets carry a core "
          f"single-cell pseudotime")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
