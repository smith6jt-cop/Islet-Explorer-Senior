#!/usr/bin/env python3
"""Does the within-donor ordering reproduce itself? Measured per donor, reported, never a gate.

``trajectory_singlecell.R`` fits each donor on a subsample of that donor's cells. Run with a
replicate seed offset it fits each donor again on an INDEPENDENTLY DRAWN subsample of the same
cells and writes the result beside the first. This compares the two.

**Three numbers, because they fail separately and the distinction is the whole point.**

``rho_shared_cells``
    Spearman between the two passes over the cells that happen to be in both subsamples. This is
    the axis itself. If it is low, the two passes found different orderings of the same donor's
    cells, and the ordering is a property of the draw rather than of the donor.

``rho_marker_profile``
    Pearson between the two passes' Spearman(pseudotime, marker) profiles over the 38 markers.
    The same question asked without needing shared cells, and it separates "a different axis" from
    "the same axis with few cells in common".

``rho_islet_core_median``
    Spearman between the two passes' per-islet core medians -- the quantity the app actually keys
    on. This can be low while the axis is perfectly reproducible, because an islet's median is over
    however many of its cells the draw happened to take: with a donor of 3,000 islets and a 13 %
    sample, two passes summarise almost disjoint cells per islet. ``median_sampled_core_cells``
    is reported beside it so the two causes are told apart rather than guessed at.

**Nothing here gates and nothing here is a biological check.** A donor whose passes disagree keeps
its ordering; what changes is that the disagreement is written down next to it instead of being
invisible. This is internal coherence -- does the build reproduce itself -- which is legitimate in
a way that reproducing a previous biological finding is not.

    .venv/bin/python core/phenocycler/scripts/typing/trajectory_reproducibility.py
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import subpanel_io as S

TRAJ = S.OUT / "trajectory"
MIN_SHARED_CELLS = 200
MIN_ISLET_CELLS = 5          #: the islet floor the summary already uses
#: Declared before the numbers were looked at: below this the two passes are not ordering the same
#: way and the donor's ordering is a property of which cells were drawn. It labels, never excludes.
RHO_REPRODUCIBLE = 0.70


def islet_core_medians(frame: pd.DataFrame) -> pd.Series:
    core = frame[frame["region"] == "core"]
    grouped = core.groupby("islet_key")["pseudotime"].agg(["median", "size"])
    return grouped.loc[grouped["size"] >= MIN_ISLET_CELLS, "median"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()

    primary_path = TRAJ / "singlecell_pseudotime.parquet"
    replicate_path = TRAJ / "singlecell_pseudotime_replicate.parquet"
    if not replicate_path.is_file():
        raise SystemExit(
            f"no replicate at {replicate_path}. Produce one with:\n"
            "  Rscript core/phenocycler/scripts/typing/trajectory_singlecell.R 25000 99001")

    meta = json.loads((TRAJ / "singlecell_trajectory.json").read_text())
    markers = list(meta["markers"])
    columns = ["donor_id", "object_id", "islet_key", "region", "donor_status", "pseudotime"]
    primary = pd.read_parquet(primary_path, columns=columns)
    replicate = pd.read_parquet(replicate_path, columns=columns)
    values = pq.read_table(TRAJ / "singlecell_features.parquet",
                           columns=["donor_id", "object_id"] + markers).to_pandas()

    print(f"  {'donor':>6} {'status':>5} {'shared':>8} {'ovlp':>6} {'rho cells':>10} "
          f"{'rho profile':>12} {'islets':>7} {'rho islet':>10} {'core/islet':>11}  verdict")
    report, rows = {}, []
    for donor in sorted(primary["donor_id"].unique()):
        a = primary[primary.donor_id == donor]
        b = replicate[replicate.donor_id == donor]
        marker_frame = values[values.donor_id == donor][["object_id"] + markers]

        shared = a.merge(b[["object_id", "pseudotime"]], on="object_id", suffixes=("_a", "_b"))
        # How much of the two passes' input is literally the same cells. A donor ordered whole
        # (6380) has an overlap of 1.000, so any disagreement there is not about the draw at all.
        overlap = len(shared) / len(a) if len(a) else np.nan
        rho_cells = (float(spearmanr(shared.pseudotime_a, shared.pseudotime_b).statistic)
                     if len(shared) >= MIN_SHARED_CELLS else np.nan)

        profiles = []
        for frame in (a, b):
            joined = frame.merge(marker_frame, on="object_id")
            profiles.append(np.array([
                spearmanr(joined[marker], joined["pseudotime"]).statistic for marker in markers]))
        rho_profile = float(np.corrcoef(profiles[0], profiles[1])[0, 1])

        medians = pd.concat([islet_core_medians(a).rename("a"),
                             islet_core_medians(b).rename("b")], axis=1).dropna()
        rho_islet = (float(spearmanr(medians.a, medians.b).statistic)
                     if len(medians) >= 10 else np.nan)
        core_counts = a[a.region == "core"].groupby("islet_key").size()
        median_core = float(core_counts.median()) if len(core_counts) else np.nan

        # The verdict names WHICH thing is unstable; it does not pass or fail the donor.
        if np.isfinite(rho_profile) and rho_profile < RHO_REPRODUCIBLE:
            # Distinguish "a different draw gave a different axis" from "the SAME cells gave a
            # different axis". The second is not a sampling property at all.
            verdict = ("axis_not_reproducible_across_seeds_on_the_same_cells"
                       if np.isfinite(overlap) and overlap >= 0.9
                       else "axis_not_reproducible_across_subsamples")
        elif np.isfinite(rho_islet) and rho_islet < RHO_REPRODUCIBLE:
            verdict = "axis_reproducible_islet_median_thin"
        else:
            verdict = "reproducible"

        status = str(a["donor_status"].iloc[0])
        report[donor] = {
            "donor_status": status,
            "shared_cells": int(len(shared)),
            "subsample_overlap_fraction": None if not np.isfinite(overlap) else round(overlap, 4),
            "rho_shared_cells": None if not np.isfinite(rho_cells) else round(rho_cells, 4),
            "rho_marker_profile": round(rho_profile, 4),
            "islets_compared": int(len(medians)),
            "rho_islet_core_median": None if not np.isfinite(rho_islet) else round(rho_islet, 4),
            "median_sampled_core_cells_per_islet": median_core,
            "verdict": verdict,
        }
        rows.append(verdict)
        print(f"  {donor:>6} {status:>5} {len(shared):>8,} {overlap:>6.3f} {rho_cells:>10.3f} "
              f"{rho_profile:>12.3f} {len(medians):>7,} {rho_islet:>10.3f} {median_core:>11.1f}  "
              f"{verdict}")

    tally = pd.Series(rows).value_counts().to_dict()
    print(f"\n  {tally}")
    target = TRAJ / "reproducibility.json"
    target.write_text(json.dumps({
        "what_this_is": ("Whether the within-donor ordering reproduces itself on an independently "
                         "drawn subsample of the same donor's cells. Reported per donor. It gates "
                         "nothing and excludes no donor."),
        "threshold_for_the_label": RHO_REPRODUCIBLE,
        "primary": str(primary_path.relative_to(S.REPO)),
        "replicate": str(replicate_path.relative_to(S.REPO)),
        "tally": tally,
        "per_donor": report,
    }, indent=2))
    print(f"-> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
