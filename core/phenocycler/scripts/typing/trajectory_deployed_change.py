#!/usr/bin/env python3
"""What the within-donor rebuild cost the ten donors that were already deployed.

The maintainer decided on 2026-09-17 to refit the trajectory within donor for all twenty donors,
**accepting that the ten donors already carrying a pseudotime would get different values**. This
measures how different, per donor, so the decision is seen rather than described.

Two comparisons, because they answer different questions:

*Cells.* The two builds sampled independently -- the pooled build took 49,999 cells by proportional
cohort quota, this one takes up to 25,000 per donor -- so only the cells in both samples can be
compared. Their count is reported per donor; it is an overlap, not a subsample of anything, and it
is not weighted.

*Islets.* This is the comparison that matters, because the islet summary is what the app reads and
the card draws. For every islet with a core value in both builds: the Spearman correlation between
the old and new within-donor ordering, and the median absolute change in an islet's **percentile
position within its own donor**. Percentile rather than raw rank, because the two builds order
different numbers of islets per donor and a raw-rank change would mostly measure that.

A correlation here is a comparison of two builds, not a result about tissue, and it gates nothing.

    .venv/bin/python core/phenocycler/scripts/typing/trajectory_deployed_change.py
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import subpanel_io as S

TRAJ = S.OUT / "trajectory"
OLD = TRAJ / "superseded_pooled_2026-09-14"


def percentile_within(frame: pd.DataFrame, key: str, value: str) -> pd.Series:
    """Rank within donor, scaled to (0, 1] -- comparable when the counts differ."""
    return frame.groupby(key, observed=True)[value].rank(method="average", pct=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
    if not OLD.is_dir():
        raise SystemExit(f"the superseded pooled build is not at {OLD}; nothing to compare")

    new_meta = json.loads((TRAJ / "singlecell_trajectory.json").read_text())
    old_meta = json.loads((OLD / "singlecell_trajectory.json").read_text())
    if new_meta.get("fit_scope") != "within_donor":
        raise SystemExit("the current build is not the within-donor one; refusing to compare")

    old_cells = pd.read_parquet(OLD / "singlecell_pseudotime.parquet",
                                columns=["donor_id", "object_id", "pseudotime"])
    new_cells = pd.read_parquet(TRAJ / "singlecell_pseudotime.parquet",
                                columns=["donor_id", "object_id", "donor_status", "pseudotime"])
    deployed = sorted(old_cells["donor_id"].astype(str).unique())
    print(f"the deployed ten: {', '.join(deployed)}")
    print(f"  pooled build   {old_meta['n_ordered']:>7,} cells of {old_meta['n_population']:>9,} "
          f"over {len(deployed)} donors, one shared embedding")
    print(f"  within-donor   {new_meta['n_ordered']:>7,} cells of {new_meta['n_population']:>9,} "
          f"over {new_meta['n_donors_ordered']} donors, one embedding each\n")

    old_summary = pd.read_parquet(OLD / "singlecell_islet_summary.parquet")
    new_summary = pd.read_parquet(TRAJ / "singlecell_islet_summary.parquet")
    for frame in (old_summary, new_summary):
        frame["donor_id"] = frame["case_id"].astype(str)

    merged_cells = old_cells.merge(
        new_cells, on=["donor_id", "object_id"], how="inner", suffixes=("_old", "_new"))

    rows, report = [], {}
    print(f"  {'donor':>6} {'status':>5} | {'cells both':>10} {'rho cells':>9} | "
          f"{'islets both':>11} {'rho islets':>10} {'med |dpct|':>10} {'max |dpct|':>10} | "
          f"{'rho delta':>9}")
    for donor in deployed:
        cell_pair = merged_cells[merged_cells["donor_id"] == donor]
        rho_cells = (float(spearmanr(cell_pair["pseudotime_old"],
                                     cell_pair["pseudotime_new"]).statistic)
                     if len(cell_pair) >= 10 else np.nan)

        old_d = old_summary[(old_summary.donor_id == donor)
                            & old_summary.sc_pseudotime_core.notna()].copy()
        new_d = new_summary[(new_summary.donor_id == donor)
                            & new_summary.sc_pseudotime_core.notna()].copy()
        old_d["pct"] = old_d["sc_pseudotime_core"].rank(method="average", pct=True)
        new_d["pct"] = new_d["sc_pseudotime_core"].rank(method="average", pct=True)
        pair = old_d.merge(new_d, on=["donor_id", "islet_key"], suffixes=("_old", "_new"))
        if len(pair) >= 10:
            rho_islets = float(spearmanr(pair.sc_pseudotime_core_old,
                                         pair.sc_pseudotime_core_new).statistic)
            shift = (pair["pct_new"] - pair["pct_old"]).abs()
            med_shift, max_shift = float(shift.median()), float(shift.max())
        else:
            rho_islets = med_shift = max_shift = np.nan

        both_delta = pair.dropna(subset=["sc_pseudotime_peri_minus_core_old",
                                         "sc_pseudotime_peri_minus_core_new"])
        rho_delta = (float(spearmanr(both_delta.sc_pseudotime_peri_minus_core_old,
                                     both_delta.sc_pseudotime_peri_minus_core_new).statistic)
                     if len(both_delta) >= 10 else np.nan)

        status = str(new_cells.loc[new_cells.donor_id == donor, "donor_status"].iloc[0])
        report[donor] = {
            "donor_status": status,
            "cells_in_both_samples": int(len(cell_pair)),
            "spearman_cells_old_vs_new": None if not np.isfinite(rho_cells) else round(rho_cells, 4),
            "islets_ordered_pooled": int(len(old_d)),
            "islets_ordered_within_donor": int(len(new_d)),
            "islets_in_both": int(len(pair)),
            "spearman_islets_old_vs_new": None if not np.isfinite(rho_islets) else round(rho_islets, 4),
            "median_abs_change_in_within_donor_percentile":
                None if not np.isfinite(med_shift) else round(med_shift, 4),
            "max_abs_change_in_within_donor_percentile":
                None if not np.isfinite(max_shift) else round(max_shift, 4),
            "islets_with_both_deltas": int(len(both_delta)),
            "spearman_peri_minus_core_old_vs_new":
                None if not np.isfinite(rho_delta) else round(rho_delta, 4),
        }
        rows.append((donor, status, len(cell_pair), rho_cells, len(pair), rho_islets,
                     med_shift, max_shift, rho_delta))
        print(f"  {donor:>6} {status:>5} | {len(cell_pair):>10,} {rho_cells:>9.3f} | "
              f"{len(pair):>11,} {rho_islets:>10.3f} {med_shift:>10.3f} {max_shift:>10.3f} | "
              f"{rho_delta:>9.3f}")

    shifts = np.array([r[6] for r in rows], dtype=float)
    rhos = np.array([r[5] for r in rows], dtype=float)
    print(f"\n  over the ten: median |percentile change| ranges "
          f"{np.nanmin(shifts):.3f}-{np.nanmax(shifts):.3f} (median of the ten "
          f"{np.nanmedian(shifts):.3f}); islet-level rho ranges "
          f"{np.nanmin(rhos):+.3f}-{np.nanmax(rhos):+.3f}")
    print("  A percentile change of 0.33 is what random re-ordering gives, so it is the yardstick:")
    print("  the expected |difference| of two independent uniforms is 1/3.")

    target = TRAJ / "deployed_ten_change.json"
    target.write_text(json.dumps({
        "what_this_is": ("How far the ten already-deployed donors moved when the trajectory was "
                         "refit within donor. The maintainer accepted this change on 2026-09-17; "
                         "this is its size."),
        "pooled_build": {"n_ordered": old_meta["n_ordered"],
                         "n_population": old_meta["n_population"],
                         "sampling": old_meta["sampling"], "fit_scope": "pooled_across_donors"},
        "within_donor_build": {"n_ordered": new_meta["n_ordered"],
                               "n_population": new_meta["n_population"],
                               "sampling": new_meta["sampling"], "fit_scope": "within_donor"},
        "random_reordering_yardstick": 1 / 3,
        "per_donor": report,
    }, indent=2))
    print(f"\n-> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
