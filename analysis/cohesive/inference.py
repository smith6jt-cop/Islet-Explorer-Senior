"""Donor-level inference for the fixed six-hypothesis family, under the frozen protocol as amended.

The donor is the unit: H1-H5 compare equal-weight donor endpoints, T1D minus ND, by two-sided label
permutation with a stratified donor bootstrap; AAb+ is shown, never tested. H6 is the paired-donor
partial rank association of pLN architecture with pancreatic H3, conditioned on disease and both mean
T states, by within-disease Freedman-Lane permutation. BH over the six. No minimum donor count: a
contrast is UNAVAILABLE only when it cannot be computed, and its bookkeeping p of 1 is never shown.
"""
from __future__ import annotations

from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .common import Context, read_json, require, write_json

PREDICTED = {"H1": "positive", "H2": "positive", "H3": "positive", "H4": "negative", "H5": "positive", "H6": "positive"}
COMPONENTS = {"H1": ("vascular_excess", "neural_excess"), "H2": ("fragmentation_z", "interface_z")}


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "infer", "Unknown inference operation")
    return []


def endpoints(context: Context, hypothesis: str) -> pd.DataFrame:
    table = pd.read_csv(context.task(hypothesis.lower()) / "donor_endpoints.csv", dtype={"donor_id": str})
    require(not table.duplicated(["donor_id", "configuration"]).any(), f"Duplicate donor endpoints: {hypothesis}")
    return table


def difference(t1d, nd):
    return float(np.mean(t1d) - np.mean(nd))


def contrast(table: pd.DataFrame, permutations: int, bootstraps: int, seed: int, column="donor_effect",
             strata: pd.Series | None = None) -> dict:
    """T1D minus ND over available donors; exact enumeration when it is no larger than the Monte Carlo."""
    usable = table.loc[table.status.eq("available") & np.isfinite(table[column].astype(float))]
    t1d = usable.loc[usable.disease_status.eq("T1D"), column].to_numpy(float)
    nd = usable.loc[usable.disease_status.eq("ND"), column].to_numpy(float)
    out = {"n_T1D": len(t1d), "n_ND": len(nd), "effect": np.nan, "p": np.nan, "CI95_low": np.nan,
           "CI95_high": np.nan, "status": "unavailable", "reason": ""}
    if not len(t1d) or not len(nd):
        computable = ", ".join(f"{g} {int((usable.disease_status == g).sum())} of {int((table.disease_status == g).sum())}"
                               for g in ("ND", "T1D"))
        return dict(out, reason=f"no computable {'T1D' if not len(t1d) else 'ND'} donor endpoint ({computable} computable)")
    observed = difference(t1d, nd)
    values = np.r_[t1d, nd]
    n = len(values)
    rng = np.random.default_rng(seed)
    if strata is None and comb(n, len(t1d)) <= permutations:
        null = np.array([values[list(c)].mean() - np.delete(values, list(c)).mean()
                         for c in combinations(range(n), len(t1d))])
        p = float(np.mean(np.abs(null) >= abs(observed) - 1e-12))
        method = "exact enumeration of disease labels"
    else:
        order = np.r_[usable.loc[usable.disease_status.eq("T1D")].index, usable.loc[usable.disease_status.eq("ND")].index]
        labels = np.tile(np.r_[np.ones(len(t1d), bool), np.zeros(len(nd), bool)], (permutations, 1))
        blocks = ([np.flatnonzero(strata.loc[order].to_numpy() == s) for s in np.unique(strata.loc[order])]
                  if strata is not None else [np.arange(n)])
        for block in blocks:
            labels[:, block] = rng.permuted(labels[:, block], axis=1)
        k = labels.sum(axis=1)
        null = (labels @ values) / k - ((~labels) @ values) / (n - k)
        p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-12)) / (permutations + 1))
        method = "Monte Carlo label permutation" + (" within acquisition strata" if strata is not None else "")
    boot = rng.choice(t1d, (bootstraps, len(t1d))).mean(axis=1) - rng.choice(nd, (bootstraps, len(nd))).mean(axis=1)
    low, high = np.quantile(boot, [0.025, 0.975])
    return dict(out, effect=observed, p=p, CI95_low=float(low), CI95_high=float(high), status="available",
                method=method, n_null=len(null), mean_T1D=float(t1d.mean()), mean_ND=float(nd.mean()))


def descriptive(table: pd.DataFrame) -> dict:
    usable = table.loc[table.status.eq("available")]
    out = {}
    for group in ("ND", "AAb+", "T1D"):
        values = usable.loc[usable.disease_status.eq(group), "donor_effect"].astype(float)
        out[f"n_available_{group}"] = len(values)
        out[f"n_total_{group}"] = int(table.disease_status.eq(group).sum())
        out[f"mean_{group}"] = float(values.mean()) if len(values) else np.nan
    return out


def h6_design(z: pd.DataFrame):
    y = rankdata(z.H3.to_numpy(float))
    x = rankdata(z.architecture_excess.to_numpy(float))
    columns = [np.ones(len(z))]
    columns.extend(pd.get_dummies(z.disease_status, drop_first=True, dtype=float).to_numpy().T)
    columns.append(rankdata(z.mean_memory_state.to_numpy(float)))
    columns.append(rankdata(z.mean_cytotoxic_state.to_numpy(float)))
    X = np.column_stack(columns)
    return x, y, X


def partial(x, y, X):
    rank = int(np.linalg.matrix_rank(X))
    if len(x) <= rank + 1:
        return None, rank
    rx = x - X @ np.linalg.lstsq(X, x, rcond=None)[0]
    ry = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    if min(rx.std(), ry.std()) < 1e-10:
        return None, rank
    return (float(np.corrcoef(rx, ry)[0, 1]), rx, ry), rank


def h6_association(paired: pd.DataFrame, permutations: int, bootstraps: int, seed: int, exclude=()) -> tuple[dict, pd.DataFrame | None]:
    z = paired.loc[paired.paired_status.eq("available") & paired.disease_status.isin(["ND", "T1D"])
                   & ~paired.donor_id.isin(exclude)].reset_index(drop=True)
    out = {"n_pairs": len(z), "n_ND": int(z.disease_status.eq("ND").sum()), "n_T1D": int(z.disease_status.eq("T1D").sum()),
           "effect": np.nan, "p": np.nan, "CI95_low": np.nan, "CI95_high": np.nan, "status": "unavailable", "reason": ""}
    if len(z) < 3:
        return dict(out, reason=f"{len(z)} paired ND/T1D donors: a partial correlation needs more pairs than its design rank"), None
    x, y, X = h6_design(z)
    result, rank = partial(x, y, X)
    out["conditioning_design_rank"] = rank
    if result is None:
        return dict(out, reason="the conditioning design leaves no residual variation in one endpoint"), None
    r, rx, ry = result
    rng = np.random.default_rng(seed)
    group = z.disease_status.to_numpy()
    slots = [np.flatnonzero(group == g) for g in np.unique(group)]
    null = np.empty(permutations)
    for b in range(permutations):
        shuffled = ry.copy()
        for index in slots:
            shuffled[index] = ry[rng.permutation(index)]
        shuffled = shuffled - X @ np.linalg.lstsq(X, shuffled, rcond=None)[0]
        null[b] = np.corrcoef(rx, shuffled)[0, 1] if shuffled.std() > 1e-12 else 0.0
    p = float((1 + np.sum(np.abs(null) >= abs(r) - 1e-12)) / (permutations + 1))
    boots = []
    for _ in range(bootstraps):
        index = np.concatenate([rng.choice(s, len(s)) for s in slots])
        resampled, _ = partial(*h6_design(z.iloc[index]))
        if resampled is not None:
            boots.append(resampled[0])
    low, high = np.quantile(boots, [0.025, 0.975]) if boots else (np.nan, np.nan)
    out.update(effect=r, p=p, CI95_low=float(low), CI95_high=float(high), status="available", n_null=permutations,
               n_valid_bootstrap=len(boots), method="within-disease Freedman-Lane residual permutation; stratified donor bootstrap")
    return out, z.assign(pLN_rank_residual=rx, H3_rank_residual=ry)


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    q[order] = np.minimum(1, np.minimum.accumulate((p[order] * m / np.arange(1, m + 1))[::-1])[::-1])
    return q


def interpret(row) -> str:
    if row["status"] != "available":
        return "UNAVAILABLE: " + row["reason"]
    if row["q_BH6"] >= 0.05:
        return "Inconclusive for the fixed prediction; nonsignificance is not refutation"
    agrees = row["effect"] < 0 if PREDICTED[row["hypothesis"]] == "negative" else row["effect"] > 0
    if not agrees:
        return "Contradicts the specified directional prediction"
    text = "Supports the specified conditional spatial prediction"
    if row["hypothesis"] == "H2" and not (row.get("fragmentation_z_effect", np.nan) > 0 and row.get("interface_z_effect", np.nan) > 0):
        return "Joint score differs in the predicted direction, but the two components do not both agree"
    if row["hypothesis"] == "H1":
        text += "; relative vascular-versus-neural alignment only"
    if row["hypothesis"] == "H4":
        text += "; the higher-order term alone does not establish net checkpoint slope weakening"
    return text


def infer(context: Context, stage: Path) -> dict:
    protocol = read_json(context.task("protocol") / "frozen_protocol.json")
    statistics = protocol["statistics"]
    permutations, bootstraps, seed = (statistics["donor_permutations"], statistics["bootstrap_replicates"], statistics["seed"])
    names = {h["id"]: h["name"] for h in protocol["hypotheses"]}
    tests, sensitivities, donors = [], [], []
    for h in ("H1", "H2", "H3", "H4", "H5"):
        table = endpoints(context, h)
        primary = table.loc[table.configuration.eq("primary")].reset_index(drop=True)
        donors.append(primary.assign(hypothesis=h))
        row = {"hypothesis": h, "name": names[h], "contrast": "T1D minus ND", "predicted_direction": PREDICTED[h],
               **contrast(primary, permutations, bootstraps, seed), **descriptive(primary)}
        for component in COMPONENTS.get(h, ()):
            row[f"{component}_effect"] = contrast(primary, permutations, bootstraps, seed, column=component)["effect"]
        tests.append(row)
        for configuration, z in table.groupby("configuration"):
            if configuration != "primary":
                sensitivities.append({"hypothesis": h, "sensitivity": configuration,
                                      **contrast(z.reset_index(drop=True), permutations, bootstraps, seed)})
        if "stain_group" in primary:
            strata = primary.stain_group.astype(str)
            sensitivities.append({"hypothesis": h, "sensitivity": "permutation_within_stain_group",
                                  **contrast(primary, permutations, bootstraps, seed, strata=strata)})
    paired = pd.read_csv(context.task("h6") / "paired_donors.csv", dtype={"donor_id": str})
    architecture = pd.read_csv(context.task("h6") / "pln_donor_architecture.csv", dtype={"donor_id": str})
    result, pairs = h6_association(paired, permutations, bootstraps, seed)
    tests.append({"hypothesis": "H6", "name": names["H6"], "contrast": "paired-donor partial rank association",
                  "predicted_direction": "positive", **result})
    if pairs is not None:
        pairs.to_csv(stage / "H6_primary_pairs.csv", index=False)
    h3 = endpoints(context, "H3").loc[lambda t: t.configuration.eq("primary"), ["donor_id", "donor_effect", "status"]]
    for configuration, z in architecture.groupby("configuration"):
        if configuration == "primary":
            continue
        pairs_c = z.merge(h3.rename(columns={"donor_effect": "H3", "status": "H3_status"}), on="donor_id", how="inner")
        pairs_c["paired_status"] = np.where(pairs_c.status.eq("available") & pairs_c.H3_status.eq("available"), "available", "unavailable")
        sensitivities.append({"hypothesis": "H6", "sensitivity": configuration,
                              **h6_association(pairs_c, permutations, bootstraps, seed)[0]})
    tests = pd.DataFrame(tests)
    bookkeeping = tests.p.where(tests.status.eq("available"), statistics["unavailable_p"]).astype(float)
    tests["q_BH6"] = bh(bookkeeping.to_numpy())
    tests.loc[~tests.status.eq("available"), ["p", "q_BH6"]] = np.nan
    tests["interpretation"] = [interpret(row) for row in tests.to_dict("records")]
    tests.to_csv(stage / "primary_tests.csv", index=False)
    sensitivities = pd.DataFrame(sensitivities).assign(
        multiplicity="declared sensitivity; outside the six-test family; unadjusted")
    sensitivities.to_csv(stage / "sensitivity_tests.csv", index=False)
    pd.concat(donors, ignore_index=True).to_csv(stage / "donor_endpoints_primary.csv", index=False)
    acquisition = pd.concat(donors, ignore_index=True)
    if "stain_group" in acquisition:
        (acquisition.loc[acquisition.hypothesis.eq("H1")].groupby(["stain_group", "disease_status"]).size()
         .unstack(fill_value=0).to_csv(stage / "acquisition_by_disease.csv"))
    write_json(stage / "inference_plan.json", {
        "unit": "independent donor, equal weight", "primary_contrast": "T1D minus ND (H1-H5)",
        "AAb+": "descriptive only", "permutations": permutations, "bootstrap_replicates": bootstraps, "seed": seed,
        "multiplicity": "BH over the six primaries; an unavailable primary enters as p=1 and is displayed as unavailable",
        "minimum_support": "none (ruled 2026-10-17); only computability",
        "literature_overlap": "No donor is removed (ruled 2026-10-18). 6551, imaged in Peters et al., anchors the "
                              "cross-study accuracy report (crossstudy stage) instead.",
        "acquisition": "stain_group-stratified permutation as a sensitivity; disease and acquisition alias at the coarse level, so no causal disease claim"})
    return {"tests": {row["hypothesis"]: {"status": row["status"], "effect": row["effect"], "p": row["p"],
                                          "q_BH6": row["q_BH6"], "interpretation": row["interpretation"]}
                      for row in tests.replace({np.nan: None}).to_dict("records")}}
