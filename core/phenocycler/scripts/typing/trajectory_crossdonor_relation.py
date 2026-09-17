#!/usr/bin/env python3
"""Relate the twenty within-donor orderings to each other -- as a separate, labelled step.

``trajectory_singlecell.R`` fits ONE diffusion map and ONE slingshot curve per donor. No donor's
cells enter another donor's kernel, so nothing in that file is a cross-donor claim and the donor
eta-squared it prints is zero by construction: every donor is rank-uniform on its own axis.

This is where the donors are compared, and it is kept separate for exactly that reason. The
question it answers is **are the twenty orderings the same ordering** -- do they order cells by the
same marker gradient -- and then, if they differ, **does the difference track acquisition**.

    profile[donor, marker] = Spearman( within-donor pseudotime , marker ) over that donor's
                             ordered cells

That profile is the donor's axis described in a vocabulary every donor shares, which the arc length
is not. Two donors whose orderings agree have similar profiles whatever their curves' units were.
The relation is then read off the profile matrix: its principal components, and the pairwise
correlations behind them.

**THE NULL.** Every cross-donor number here is reported with a permutation null, because the cohort
is small and aliased: twenty donors over seven acquisition runs and nine staining groups, with
every ND donor in the first acquisition batch and every second-batch donor Aab+ or T1D. An
eta-squared over seven groups and twenty donors is large *on random labels*, and quoting the raw
value without the null would read as structure. The unit permuted is the DONOR, because the donor
is the unit that carries a batch label; permuting cells would manufacture significance out of
per-donor cell counts.

**NOTHING HERE GATES ANYTHING.** If the relation is a batch axis, that is written down and the
within-donor orderings are used exactly as they were. A cross-donor ordering is not produced, not
used to sort the donors, and not offered to the app.

    .venv/bin/python core/phenocycler/scripts/typing/trajectory_crossdonor_relation.py
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import rankdata

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

TRAJ = S.OUT / "trajectory"
N_PERMUTATIONS = 10_000
SEED = 20260917

#: The coarse acquisition batch, derived rather than restated: `cohort.BATCH2_REEXPORTED` IS the
#: second batch, and it is the level at which batch is aliased with disease status.
BATCH2 = frozenset(COHORT.BATCH2_REEXPORTED)


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    """SS_between / SS_total. Descriptive; never a gate."""
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(values)
    values, groups = values[ok], np.asarray(groups)[ok]
    if values.size < 2:
        return float("nan")
    grand = values.mean()
    total = float(((values - grand) ** 2).sum())
    if total <= 0:
        return float("nan")
    between = sum(
        (groups == level).sum() * (values[groups == level].mean() - grand) ** 2
        for level in np.unique(groups)
    )
    return float(between / total)


def permutation_null(values: np.ndarray, groups: np.ndarray, *, rng: np.random.Generator,
                     n: int = N_PERMUTATIONS) -> dict:
    """Observed eta-squared against eta-squared on permuted DONOR labels.

    p is the usual add-one estimate, so it is never reported as exactly zero: with 10,000
    permutations the smallest attainable value is 1/10,001, and that is what it says.
    """
    observed = eta_squared(values, groups)
    if not np.isfinite(observed):
        return {"eta_squared": None, "state": "unavailable",
                "reason": "eta_squared_not_estimable(constant_or_single_group)"}
    groups = np.asarray(groups)
    null = np.array([eta_squared(values, rng.permutation(groups)) for _ in range(n)])
    null = null[np.isfinite(null)]
    p = (1 + int((null >= observed).sum())) / (1 + null.size)
    return {
        "eta_squared": round(float(observed), 4),
        "null_median": round(float(np.median(null)), 4),
        "null_p95": round(float(np.quantile(null, 0.95)), 4),
        "p_permutation": round(float(p), 5),
        "n_permutations": int(null.size),
        "n_donors": int(np.isfinite(np.asarray(values, dtype=float)).sum()),
        "n_groups": int(len(np.unique(groups))),
        "state": "measured", "reason": "",
    }


def donor_profiles(cells: pd.DataFrame, markers: list[str]) -> pd.DataFrame:
    """Spearman(pseudotime, marker) per donor -- the axis in a shared vocabulary."""
    rows = {}
    for donor, group in cells.groupby("donor_id", observed=True):
        pseudo = rankdata(group["pseudotime"].to_numpy())
        out = {}
        for marker in markers:
            values = group[marker].to_numpy(dtype=float)
            ok = np.isfinite(values)
            # A marker unavailable on this donor's cells is unavailable, not zero and not a
            # correlation of zero: the profile entry stays absent and the PCA is told so.
            if ok.sum() < 3 or np.nanstd(values[ok]) == 0:
                out[marker] = np.nan
                continue
            out[marker] = float(np.corrcoef(rankdata(values[ok]), pseudo[ok])[0, 1])
        rows[donor] = out
    return pd.DataFrame(rows).T.sort_index()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    args = parser.parse_args()
    rng = np.random.default_rng(SEED)

    meta = json.loads((TRAJ / "singlecell_trajectory.json").read_text())
    if meta.get("fit_scope") != "within_donor":
        raise SystemExit(
            "refusing: singlecell_trajectory.json does not record fit_scope 'within_donor'. "
            "This script relates per-donor orderings and has nothing to say about a pooled one.")
    markers = list(meta["markers"])

    cells = pd.read_parquet(TRAJ / "singlecell_pseudotime.parquet",
                            columns=["donor_id", "object_id", "islet_key", "region",
                                     "donor_status", "stain_group", "pseudotime"])
    values = pq.read_table(TRAJ / "singlecell_features.parquet",
                           columns=["donor_id", "object_id"] + markers).to_pandas()
    cells = cells.merge(values, on=["donor_id", "object_id"], how="left", validate="one_to_one")
    print(f"{len(cells):,} ordered cells over {cells['donor_id'].nunique()} donors")

    runs = pd.read_csv(S.REPO / "data/donor_acquisition_runs.csv").rename(
        columns={"donor_id": "case_id"}).drop_duplicates("case_id")
    runs["case_id"] = runs["case_id"].astype(str)
    acq = runs.set_index("case_id")["acq_run"]

    profiles = donor_profiles(cells, markers)
    donors = list(profiles.index)
    factors = pd.DataFrame({
        "acq_run": [str(acq.get(d, "unknown")) for d in donors],
        "stain_group": [str(cells.loc[cells.donor_id == d, "stain_group"].iloc[0]) for d in donors],
        "coarse_batch": ["batch2" if d in BATCH2 else "batch1" for d in donors],
        "donor_status": [str(cells.loc[cells.donor_id == d, "donor_status"].iloc[0])
                         for d in donors],
    }, index=donors)

    usable = [m for m in markers if profiles[m].notna().all()]
    print(f"profile matrix {len(donors)} donors x {len(usable)} markers "
          f"({len(markers) - len(usable)} markers absent on at least one donor and left out)")
    matrix = profiles[usable].to_numpy(float)
    centred = matrix - matrix.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    scores = u * s
    variance = (s ** 2) / (s ** 2).sum()

    pairwise = np.corrcoef(matrix)
    off = pairwise[np.triu_indices_from(pairwise, k=1)]
    print(f"\npairwise profile correlation between donors: median {np.median(off):+.3f}  "
          f"min {off.min():+.3f}  max {off.max():+.3f}")
    print(f"profile PC variance explained: PC1 {variance[0]:.3f}  PC2 {variance[1]:.3f}  "
          f"PC3 {variance[2]:.3f}")

    relation: dict[str, dict] = {}
    print(f"\ncross-donor relation, permutation null over DONORS ({args.permutations:,} draws)")
    print(f"  {'quantity':<28} {'factor':<14} {'eta2':>7} {'null med':>9} {'null p95':>9} {'p':>8}")
    quantities = {
        "profile_PC1": scores[:, 0],
        "profile_PC2": scores[:, 1],
    }
    # The donor-level summary of the quantity the card leads on: the median islet's peri-minus-core
    # on that donor's own axis. It is a within-islet difference, so the donors can be compared on
    # it without any pooled embedding -- but the comparison is still a cross-donor quantity and
    # still gets the null.
    summary_path = TRAJ / "singlecell_islet_summary.parquet"
    if summary_path.is_file():
        summary = pd.read_parquet(summary_path)
        summary["donor_id"] = summary["case_id"].astype(str)
        per_donor = summary.groupby("donor_id")["sc_pseudotime_peri_minus_core"].median()
        quantities["median_peri_minus_core"] = np.array(
            [per_donor.get(d, np.nan) for d in donors], dtype=float)

    for name, vector in quantities.items():
        relation[name] = {}
        for factor in ("acq_run", "stain_group", "coarse_batch", "donor_status"):
            result = permutation_null(vector, factors[factor].to_numpy(),
                                      rng=rng, n=args.permutations)
            relation[name][factor] = result
            if result["state"] == "measured":
                print(f"  {name:<28} {factor:<14} {result['eta_squared']:>7.4f} "
                      f"{result['null_median']:>9.4f} {result['null_p95']:>9.4f} "
                      f"{result['p_permutation']:>8.5f}")
            else:
                print(f"  {name:<28} {factor:<14} {'Unavailable':>7}  {result['reason']}")

    profiles.assign(**factors).to_parquet(TRAJ / "crossdonor_axis_profiles.parquet")
    (TRAJ / "crossdonor_relation.json").write_text(json.dumps({
        "what_this_is": (
            "The relation BETWEEN the twenty within-donor orderings. It is not an ordering, it "
            "does not sort the donors, and nothing downstream reads it."),
        "unit_permuted": "donor",
        "n_permutations": int(args.permutations),
        "seed": SEED,
        "profile": ("Spearman(within-donor pseudotime, marker) per donor; markers absent on any "
                    "donor are left out of the PCA rather than filled"),
        "n_donors": len(donors), "donors": donors,
        "n_markers_in_profile": len(usable),
        "markers_left_out": sorted(set(markers) - set(usable)),
        "pairwise_profile_correlation": {
            "median": round(float(np.median(off)), 4),
            "min": round(float(off.min()), 4), "max": round(float(off.max()), 4)},
        "pc_variance_explained": [round(float(v), 4) for v in variance[:5]],
        "relation": relation,
    }, indent=2))
    print(f"\n-> {TRAJ / 'crossdonor_relation.json'}")
    print(f"-> {TRAJ / 'crossdonor_axis_profiles.parquet'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
