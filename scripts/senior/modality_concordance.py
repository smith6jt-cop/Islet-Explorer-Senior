#!/usr/bin/env python
"""Protein (PhenoCycler) vs RNA (Xenium) lineage-composition concordance.

WHAT THIS MEASURES
------------------
Per donor, per lineage, the share of that lineage within the *comparable* lineage
set, computed independently in each modality, and the difference between them.

These are DIFFERENT SECTIONS OF DIFFERENT BLOCKS.  They are not two views of one
cell, and no cell identifier crosses: a PhenoCycler ``object_id`` is a QuPath UUID
and a Xenium ``cell_id`` is 10x segmentation, unique only within a bundle.  The
join is on ``donor_id`` and nothing else.  Composition against composition is the
only comparison these data license.

Nothing here is tuned to make the modalities agree.  Agreement is informative and
disagreement is informative; the deliverable is which lineages may carry a
cross-modality claim and which may not.

THE DISAGREEMENT THRESHOLD IS FIXED IN THIS FILE BEFORE ANY NUMBER IS COMPUTED.
See ``THRESHOLD_*`` below and the rationale beside them.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: E402

CROSSWALK = REPO / "core/phenocycler/data/integration/vocab_crosswalk.csv"
ACQ_RUNS = REPO / "data/donor_acquisition_runs.csv"
DONOR_META = REPO / "data/development/typing/app/donor_metadata.csv"

PLN_PHENO = REPO / "data/development/typing/pln/pheno"
PLN_XEN = REPO / "data/development/xenium_typing/release_pln_immune_v2"
PANC_PHENO = REPO / "data/development/typing/cells"
PANC_XEN_CANDIDATES = [
    REPO / "data/development/xenium_typing/release_panc_v1",   # fill-panc-rna track
    REPO / "data/development/xenium_typing/release",           # what exists today
]

OUT = REPO / "data/development/integration/concordance"

# ---------------------------------------------------------------------------
# THRESHOLD, FIXED BEFORE THE NUMBERS
# ---------------------------------------------------------------------------
# A lineage share is a part of a composition over K comparable lineages, so the
# mean share is 100/K percentage points.  For the pLN that is three parts and a
# mean share of 33.3 pp.
#
# These sections are different blocks of the same donor, so *some* disagreement
# is sectioning and not method: which part of a node a section cuts -- cortex,
# paracortex, medulla, capsule -- moves a B-vs-T share, and no modality can be
# blamed for that.  What sectioning cannot do is reorder the composition.
#
#   |median over donors of (protein share - RNA share)| <= 10 pp  -> CONCORDANT
#       10 pp is under a third of the mean share and inside the range regional
#       sampling of one node plausibly produces.
#   10 pp < |median| <= 20 pp                                     -> GREY
#       usable only with the offset stated alongside every number.
#   |median| > 20 pp                                              -> STRUCTURAL
#       20 pp is 60% of the mean share; a shift that large changes which lineage
#       dominates, which sectioning does not do.
#
# Independently of the median, a DISPERSION test: if the inter-quartile range of
# the per-donor difference exceeds 20 pp the disagreement is donor-dependent, and
# a donor-dependent offset cannot be carried by any pooled cross-modality claim
# even when it happens to cancel in the median.
#
# And a RANK test: if the two modalities disagree about the rank order of the
# comparable lineages in more than half the donors, the composition is structurally
# different whatever the per-lineage medians say.
THRESHOLD_CONCORDANT_PP = 10.0
THRESHOLD_STRUCTURAL_PP = 20.0
THRESHOLD_DISPERSION_IQR_PP = 20.0
THRESHOLD_RANK_DISAGREE_FRAC = 0.50

N_PERM = 10000
RNG_SEED = 20260917

# PhenoCycler pLN immune_class -> crosswalk common lineage.
PLN_PHENO_CLASS_TO_COMMON = {
    "T_cell": "T_NK",
    "NK_ILC": "T_NK",
    "B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Granulocyte": "Granulocyte",
}

# PhenoCycler pancreas cell_class -> common lineage (non-immune).
PANC_PHENO_CLASS_TO_COMMON = {
    "Beta": "Endocrine",
    "Alpha": "Endocrine",
    "Delta": "Endocrine",
    "Exocrine": "Exocrine",
    "Blood vessel": "Vascular",
    "Nerve": "Neural",
}
# PhenoCycler pancreas immune_definition -> common lineage.  The crosswalk records
# T_NK / B_Plasma / Myeloid as pheno_resolution == "subclass": recovered from the
# per-marker gates, not from the broad label, which is exactly this column.
PANC_PHENO_IMMDEF_TO_COMMON = {
    "Imm-T-cell": "T_NK",
    "Imm-NK-cell": "T_NK",
    "Imm-B-cell": "B_Plasma",
    "Imm-plasma-cell": "B_Plasma",
    "Imm-macrophage": "Myeloid",
    "Imm-dendritic-cell": "Myeloid",
    "Imm-neutrophil": "Granulocyte",
}


def log(msg: str) -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# crosswalk
# ---------------------------------------------------------------------------
def load_crosswalk() -> pd.DataFrame:
    cw = pd.read_csv(CROSSWALK, dtype=str).fillna("")
    return cw[cw["kind"] == "lineage"].copy()


def comparable_lineages(cw: pd.DataFrame, *, immune_only: bool) -> tuple[list[str], dict[str, str]]:
    """The lineages a protein-vs-RNA comparison may use, and why each other one is out.

    A lineage is comparable only when the crosswalk gives it a name on BOTH sides
    and neither side's resolution disqualifies it.  Granulocyte is excluded by
    ``xenium_resolution == fine_only`` -- it has no broad class on the RNA side, so
    a Granulocyte share is not the same quantity in the two modalities.
    """
    keep: list[str] = []
    dropped: dict[str, str] = {}
    for _, r in cw.iterrows():
        name = r["common"]
        if immune_only and r["is_immune"] != "yes":
            dropped[name] = f"not an immune lineage (is_immune={r['is_immune']!r})"
            continue
        if not r["phenocycler"]:
            dropped[name] = f"no PhenoCycler label (pheno_resolution={r['pheno_resolution']!r})"
            continue
        if not r["xenium"]:
            dropped[name] = f"no Xenium label (xenium_resolution={r['xenium_resolution']!r})"
            continue
        if r["xenium_resolution"] == "fine_only":
            dropped[name] = "xenium_resolution=fine_only; no broad RNA class, so not the same quantity"
            continue
        if r["pheno_resolution"] in {"no", "never"}:
            dropped[name] = f"pheno_resolution={r['pheno_resolution']!r}"
            continue
        keep.append(name)
    return keep, dropped


def asymmetry_notes(cw: pd.DataFrame, lineages: list[str]) -> dict[str, str]:
    notes = {}
    for _, r in cw.iterrows():
        if r["common"] not in lineages:
            continue
        bits = []
        if r["pheno_resolution"] != "yes":
            bits.append(f"protein side is {r['pheno_resolution']} ({r['pheno_note']})")
        if r["xenium_resolution"] != "yes":
            bits.append(f"RNA side is {r['xenium_resolution']} ({r['xenium_note']})")
        if r["collapsed"]:
            bits.append(f"folds {r['collapsed']}")
        if bits:
            notes[r["common"]] = "; ".join(bits)
    return notes


# ---------------------------------------------------------------------------
# readers (retry: a concurrent workflow owns these trees and may be mid-write)
# ---------------------------------------------------------------------------
def read_parquet_retry(path: Path, columns=None, attempts: int = 5):
    import time
    last = None
    for i in range(attempts):
        try:
            return pd.read_parquet(path, columns=columns)
        except Exception as exc:  # partial write -> retry, never record a partial count
            last = exc
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"could not read {path} after {attempts} attempts: {last}")


def pln_protein_counts(donors_excluded: set[str]) -> tuple[pd.DataFrame, dict]:
    status = pd.read_csv(PLN_PHENO / "per_donor_status.csv", dtype={"donor_id": str})
    rows, states = [], {}
    for _, s in status.iterrows():
        d = str(s["donor_id"])
        if d in donors_excluded:
            states[d] = "donor centrally excluded"
            continue
        if s["state"] != "typed":
            states[d] = f"{s['state']}: {s['reason']}"
            continue
        p = PLN_PHENO / f"donor_id={d}" / "roi=pln_1" / "data.parquet"
        if not p.exists():
            states[d] = "typed in status table but no partition present"
            continue
        df = read_parquet_retry(p, columns=["object_id", "is_immune", "immune_class"])
        ic = df["immune_class"].astype(str)
        common = ic.map(PLN_PHENO_CLASS_TO_COMMON)
        vc = common.value_counts(dropna=True)
        rec = {"donor_id": d, "n_cells": len(df), "n_immune": int(df["is_immune"].sum())}
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)
        rec["n_immune_unnamed"] = int(df["is_immune"].sum() - vc.sum())
        rows.append(rec)
        states[d] = "typed"
        log(f"  protein pLN {d}: {len(df):,} cells")
    return pd.DataFrame(rows), states


def pln_rna_counts(donors_excluded: set[str]) -> tuple[pd.DataFrame, dict]:
    rows, states = [], {}
    for dd in sorted((PLN_XEN / "calls").glob("donor_id=*")):
        d = dd.name.split("=", 1)[1]
        if d in donors_excluded:
            states[d] = "donor centrally excluded"
            continue
        p = dd / "roi=pln_1" / "data.parquet"
        if not p.exists():
            states[d] = "no roi=pln_1 partition"
            continue
        df = read_parquet_retry(
            p, columns=["cell_id", "common_lineage", "lineage", "nucleus_count",
                        "nucleus_measured", "depth", "y_probe_detected"]
        )
        n_raw = len(df)
        # A cell with a MEASURED absent nucleus is excluded from every calculation.
        # A cell whose nucleus was never measured is a missing measurement -- not a
        # measured zero, and not evidence that the cell has no nucleus.  Some sections
        # carry no nucleus table at all (nucleus_count is the -1 sentinel for every
        # cell); dropping those would discard whole sections on the strength of an
        # absent measurement.  Those cells are kept and the section is flagged, and
        # the flag is carried into every downstream table.
        measured = df["nucleus_measured"].fillna(False)
        no_nuc = int((measured & (df["nucleus_count"] <= 0)).sum())
        unmeasured = int((~measured).sum())
        df = df[~(measured & (df["nucleus_count"] <= 0))]
        vc = df["common_lineage"].replace("", np.nan).value_counts(dropna=True)
        rec = {
            "donor_id": d, "n_cells": n_raw, "n_retained": len(df),
            "n_no_nucleus_excluded": no_nuc, "n_nucleus_unmeasured": unmeasured,
            "nucleus_measured_frac": float(measured.mean()) if n_raw else np.nan,
            "median_depth": float(df["depth"].median()) if len(df) else np.nan,
            "y_probe_rate": float(df["y_probe_detected"].mean()) if len(df) else np.nan,
        }
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)
        rows.append(rec)
        states[d] = "called"
        log(f"  RNA pLN {d}: {n_raw:,} cells, {len(df):,} retained, nucleus measured on {100*float(measured.mean()):.0f}%")
    return pd.DataFrame(rows), states


def panc_protein_counts(donors: set[str], donors_excluded: set[str]) -> pd.DataFrame:
    rows = []
    for d in sorted(donors):
        if d in donors_excluded:
            continue
        p = PANC_PHENO / f"donor_id={d}" / "data.parquet"
        if not p.exists():
            continue
        df = read_parquet_retry(
            p, columns=["object_id", "cell_class", "immune_definition", "qc_analysis_eligible"]
        )
        df = df[df["qc_analysis_eligible"].fillna(False)]
        cc = df["cell_class"].astype(str)
        common = cc.map(PANC_PHENO_CLASS_TO_COMMON)
        imm = cc.eq("Immune")
        common = common.where(~imm, df["immune_definition"].astype(str).map(PANC_PHENO_IMMDEF_TO_COMMON))
        vc = common.value_counts(dropna=True)
        rec = {"donor_id": d, "n_cells": len(df),
               "n_immune_unnamed": int((imm & common.isna()).sum())}
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)
        rows.append(rec)
        log(f"  protein panc {d}: {len(df):,} eligible cells")
    return pd.DataFrame(rows)


def panc_rna_counts(donors_excluded: set[str]) -> tuple[pd.DataFrame, Path | None, list[str]]:
    notes = []
    root = None
    for cand in PANC_XEN_CANDIDATES:
        if (cand / "calls").exists():
            root = cand
            break
    if root is None:
        return pd.DataFrame(), None, ["no pancreas Xenium calls tree present at all"]
    if root != PANC_XEN_CANDIDATES[0]:
        notes.append(
            f"release_panc_v1 (the fill-panc-rna track) was not present when this ran; "
            f"measured instead on {root.relative_to(REPO)}"
        )
    rel = json.loads((root / "release.json").read_text())
    notes.append(f"pancreas RNA release_id = {rel.get('release_id')}, lane = {rel.get('producer', {}).get('lane')}")
    rows = []
    for dd in sorted((root / "calls").glob("donor_id=*")):
        d = dd.name.split("=", 1)[1]
        if d in donors_excluded:
            continue
        parts = list(dd.glob("roi=*/data.parquet"))
        if not parts:
            continue
        df = pd.concat([read_parquet_retry(p) for p in parts], ignore_index=True)
        cols = set(df.columns)
        if "nucleus_count" in cols and "nucleus_measured" in cols:
            df = df[df["nucleus_measured"] & (df["nucleus_count"] > 0)]
            nuc = "nucleus filter applied"
        else:
            nuc = "NO nucleus column in this release; nucleated-only filter could NOT be applied"
        vc = df["common_lineage"].replace("", np.nan).value_counts(dropna=True)
        rec = {"donor_id": d, "n_cells": len(df), "nucleus_handling": nuc}
        for k, v in vc.items():
            rec[f"n_{k}"] = int(v)
        rows.append(rec)
        log(f"  RNA panc {d}: {len(df):,} cells ({nuc})")
    return pd.DataFrame(rows), root, notes


# ---------------------------------------------------------------------------
# composition
# ---------------------------------------------------------------------------
def shares(counts: pd.DataFrame, lineages: list[str]) -> pd.DataFrame:
    if counts.empty:
        return pd.DataFrame()
    cols = [f"n_{l}" for l in lineages]
    present = [c for c in cols if c in counts.columns]
    sub = counts.reindex(columns=["donor_id"] + present).copy()
    for c in cols:
        if c not in sub.columns:
            sub[c] = np.nan  # the label does not exist in this modality: unavailable, not zero
    sub["n_comparable"] = sub[present].sum(axis=1)
    for l in lineages:
        sub[f"share_{l}"] = 100.0 * sub[f"n_{l}"] / sub["n_comparable"]
    return sub


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    ok = ~np.isnan(values)
    v, g = values[ok], groups[ok]
    if len(v) < 3 or len(np.unique(g)) < 2:
        return np.nan
    gm = v.mean()
    ss_tot = ((v - gm) ** 2).sum()
    if ss_tot <= 0:
        return np.nan
    ss_b = sum(len(v[g == k]) * (v[g == k].mean() - gm) ** 2 for k in np.unique(g))
    return float(ss_b / ss_tot)


def eta_with_permutation(values: np.ndarray, groups: np.ndarray, rng) -> tuple[float, float, float]:
    """eta-squared, its permutation-null median, and p.

    The null is the same values with the batch labels shuffled over donors, which
    is the only null that answers "is this more batch structure than a random
    partition of these donors into these group sizes would show".
    """
    obs = eta_squared(values, groups)
    if np.isnan(obs):
        return np.nan, np.nan, np.nan
    ok = ~np.isnan(values)
    v, g = values[ok], groups[ok].copy()
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        rng.shuffle(g)
        null[i] = eta_squared(v, g)
    p = (1 + int((null >= obs).sum())) / (1 + N_PERM)
    return obs, float(np.nanmedian(null)), float(p)


MIN_DONORS_FOR_VERDICT = 4


def verdict(med: float, iqr: float, rank_disagree: float, n_donors: int) -> tuple[str, str]:
    if np.isnan(med):
        return "UNAVAILABLE", "one modality does not produce this label"
    if n_donors < MIN_DONORS_FOR_VERDICT:
        return "INSUFFICIENT DONORS FOR A VERDICT", (
            f"{n_donors} donor(s); a median and an inter-quartile range over fewer than "
            f"{MIN_DONORS_FOR_VERDICT} donors is not a measurement of a cohort, and the "
            "numbers below are reported as observations only"
        )
    reasons = []
    if abs(med) > THRESHOLD_STRUCTURAL_PP:
        reasons.append(f"median offset {med:+.1f} pp exceeds the {THRESHOLD_STRUCTURAL_PP:.0f} pp structural threshold")
    if iqr > THRESHOLD_DISPERSION_IQR_PP:
        reasons.append(f"IQR {iqr:.1f} pp exceeds the {THRESHOLD_DISPERSION_IQR_PP:.0f} pp dispersion threshold; the offset is donor-dependent")
    if rank_disagree > THRESHOLD_RANK_DISAGREE_FRAC:
        reasons.append(f"the modalities disagree on composition rank order in {100*rank_disagree:.0f}% of donors")
    if reasons:
        return "MAY NOT CARRY A CROSS-MODALITY CLAIM", "; ".join(reasons)
    if abs(med) <= THRESHOLD_CONCORDANT_PP:
        return "MAY CARRY A CROSS-MODALITY CLAIM", f"median offset {med:+.1f} pp, IQR {iqr:.1f} pp, both inside the thresholds fixed in advance"
    return "MAY CARRY A CROSS-MODALITY CLAIM ONLY WITH THE OFFSET STATED", f"median offset {med:+.1f} pp is in the grey band ({THRESHOLD_CONCORDANT_PP:.0f}-{THRESHOLD_STRUCTURAL_PP:.0f} pp)"


def compare(tissue: str, prot: pd.DataFrame, rna: pd.DataFrame, lineages: list[str],
            batch: pd.DataFrame, rng) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if prot.empty or rna.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    sp = shares(prot, lineages).set_index("donor_id")
    sr = shares(rna, lineages).set_index("donor_id")
    donors = sorted(set(sp.index) & set(sr.index), key=lambda x: (len(x), x))
    long = []
    for d in donors:
        for l in lineages:
            a, b = sp.loc[d, f"share_{l}"], sr.loc[d, f"share_{l}"]
            long.append({
                "tissue": tissue, "donor_id": d, "lineage": l,
                "n_protein": sp.loc[d, f"n_{l}"], "n_rna": sr.loc[d, f"n_{l}"],
                "protein_comparable_total": sp.loc[d, "n_comparable"],
                "rna_comparable_total": sr.loc[d, "n_comparable"],
                "share_protein_pct": a, "share_rna_pct": b,
                "diff_pp": a - b,
            })
    long = pd.DataFrame(long)

    # Rank-order disagreement.  Two readings, both reported:
    #   whole_rank_frac  -- fraction of donors where the FULL ordering of the
    #                       comparable set differs between modalities.  This is the
    #                       criterion as first written; it is one number for the
    #                       whole composition, so it cannot separate one lineage from
    #                       another and it is reported as context, not as a verdict.
    #   per-lineage      -- fraction of donors in which THIS lineage sits at a
    #                       different position in its own modality's ordering.  This
    #                       is the lineage-specific reading of the same criterion and
    #                       is what the per-lineage verdict uses.  The 50% threshold
    #                       is unchanged.
    whole, per_lineage = {}, {l: [] for l in lineages}
    for d in donors:
        pa = np.array([sp.loc[d, f"share_{l}"] for l in lineages], dtype=float)
        rb = np.array([sr.loc[d, f"share_{l}"] for l in lineages], dtype=float)
        if np.isnan(pa).any() or np.isnan(rb).any():
            continue
        op, orr = list(np.argsort(pa)), list(np.argsort(rb))
        whole[d] = int(op != orr)
        rp, rr = np.argsort(np.argsort(pa)), np.argsort(np.argsort(rb))
        for i, l in enumerate(lineages):
            per_lineage[l].append(int(rp[i] != rr[i]))
    whole_rank_frac = float(np.mean(list(whole.values()))) if whole else np.nan

    bat = batch.set_index("donor_id")
    summary, eta_rows = [], []
    for l in lineages:
        sub = long[long.lineage == l]
        dif = sub["diff_pp"].to_numpy(dtype=float)
        med = float(np.nanmedian(dif)) if len(dif) else np.nan
        iqr = float(np.nanpercentile(dif, 75) - np.nanpercentile(dif, 25)) if np.isfinite(dif).sum() > 1 else np.nan
        rank_frac = float(np.mean(per_lineage[l])) if per_lineage[l] else np.nan
        v, why = verdict(med, iqr, rank_frac, int(np.isfinite(dif).sum()))
        summary.append({
            "tissue": tissue, "lineage": l, "n_donors": int(np.isfinite(dif).sum()),
            "median_share_protein_pct": float(np.nanmedian(sub["share_protein_pct"])),
            "median_share_rna_pct": float(np.nanmedian(sub["share_rna_pct"])),
            "median_diff_pp": med, "iqr_diff_pp": iqr,
            "min_diff_pp": float(np.nanmin(dif)) if np.isfinite(dif).any() else np.nan,
            "max_diff_pp": float(np.nanmax(dif)) if np.isfinite(dif).any() else np.nan,
            "rank_disagreement_frac": rank_frac,
            "whole_composition_rank_disagreement_frac": whole_rank_frac,
            "verdict": v, "verdict_reason": why,
        })
        # batch eta-squared with a permutation null, for every cross-donor quantity
        for qname, series in (("share_protein_pct", sub.set_index("donor_id")["share_protein_pct"]),
                              ("share_rna_pct", sub.set_index("donor_id")["share_rna_pct"]),
                              ("diff_pp", sub.set_index("donor_id")["diff_pp"])):
            for bvar in ("acq_run", "stain_group"):
                idx = [d for d in series.index if d in bat.index]
                if len(idx) < 3:
                    continue
                vals = series.loc[idx].to_numpy(dtype=float)
                grp = bat.loc[idx, bvar].astype(str).to_numpy()
                obs, nullmed, p = eta_with_permutation(vals, grp, rng)
                eta_rows.append({
                    "tissue": tissue, "lineage": l, "quantity": qname, "batch_var": bvar,
                    "n_donors": len(idx), "n_groups": int(len(set(grp))),
                    "eta_squared": obs, "perm_null_median": nullmed, "p_perm": p,
                    "n_perm": N_PERM,
                })
    return long, pd.DataFrame(summary), pd.DataFrame(eta_rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    excluded = set(DONOR_EXCLUSIONS)
    log(f"centrally excluded donors (from cohort.DONOR_EXCLUSIONS): {sorted(excluded) or 'none'}")

    cw = load_crosswalk()
    pln_lin, pln_dropped = comparable_lineages(cw, immune_only=True)
    panc_lin, panc_dropped = comparable_lineages(cw, immune_only=False)
    log(f"comparable pLN lineages: {pln_lin}")
    log(f"comparable pancreas lineages: {panc_lin}")

    batch = pd.read_csv(ACQ_RUNS, dtype={"donor_id": str})[["donor_id", "acq_run", "stain_group"]]
    batch = batch.drop_duplicates("donor_id")
    meta = pd.read_csv(DONOR_META, dtype={"case_id": str}).rename(columns={"case_id": "donor_id"})
    meta = meta[["donor_id", "donor_status"]].drop_duplicates("donor_id")

    log("reading pLN protein ...")
    pln_p, pln_p_state = pln_protein_counts(excluded)
    log("reading pLN RNA ...")
    pln_r, pln_r_state = pln_rna_counts(excluded)

    pln_long, pln_sum, pln_eta = compare("pln", pln_p, pln_r, pln_lin, batch, rng)

    # Controls, reported beside the comparison rather than used to gate it.
    # The Y-probe rate is a control only on donors with no Y template, so it is
    # reported against sex.  Transcript depth is the flag that matters for a 5K
    # panel call; a low-depth section is flagged on depth, not on a control rate.
    controls = pln_r.merge(meta, on="donor_id", how="left")
    sex = pd.read_csv(DONOR_META, dtype={"case_id": str}).rename(columns={"case_id": "donor_id"})
    controls = controls.merge(sex[["donor_id", "sex"]].drop_duplicates("donor_id"), on="donor_id", how="left")
    controls = controls[["donor_id", "sex", "donor_status", "n_cells", "n_retained",
                         "n_no_nucleus_excluded", "nucleus_measured_frac",
                         "median_depth", "y_probe_rate"]]
    controls.to_csv(out / "pln_rna_controls.csv", index=False)
    log("\n=== pLN RNA controls (Y-probe rate is a control only where there is no Y template) ===")
    log(controls.sort_values(["sex", "y_probe_rate"]).to_string(index=False))

    log("reading pancreas RNA ...")
    panc_r, panc_root, panc_notes = panc_rna_counts(excluded)
    panc_donors = set(panc_r["donor_id"]) if not panc_r.empty else set()
    log(f"pancreas RNA donors available: {sorted(panc_donors) or 'none'}")
    log("reading pancreas protein ...")
    panc_p = panc_protein_counts(panc_donors, excluded)

    # a lineage with no count column on one side is unavailable there, not zero
    panc_lin_eff = [l for l in panc_lin
                    if (f"n_{l}" in panc_p.columns) and (f"n_{l}" in panc_r.columns)]
    panc_unavail = {l: "no count produced by one modality on this tissue"
                    for l in panc_lin if l not in panc_lin_eff}
    panc_long, panc_sum, panc_eta = compare("pancreas", panc_p, panc_r, panc_lin_eff, batch, rng)

    for name, df in (("pln_per_donor_counts_protein", pln_p),
                     ("pln_per_donor_counts_rna", pln_r),
                     ("pln_per_donor_per_lineage", pln_long),
                     ("pln_lineage_summary", pln_sum),
                     ("pln_batch_eta_squared", pln_eta),
                     ("pancreas_per_donor_counts_protein", panc_p),
                     ("pancreas_per_donor_counts_rna", panc_r),
                     ("pancreas_per_donor_per_lineage", panc_long),
                     ("pancreas_lineage_summary", panc_sum),
                     ("pancreas_batch_eta_squared", panc_eta)):
        if isinstance(df, pd.DataFrame) and not df.empty:
            df.to_csv(out / f"{name}.csv", index=False)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "script": "scripts/senior/modality_concordance.py",
        "latch": "development; no production typing is declared and nothing here is a release",
        "join": "donor_id ONLY. No cell identifier crosses modalities: PhenoCycler object_id is a "
                "QuPath UUID, Xenium cell_id is 10x segmentation unique only within a bundle.",
        "claim_scope": "COMPOSITIONAL ONLY. Different sections of different blocks. No distance, "
                       "adjacency, neighbourhood or shared scatter is computed or licensed here.",
        "thresholds_fixed_before_measurement": {
            "concordant_pp": THRESHOLD_CONCORDANT_PP,
            "structural_pp": THRESHOLD_STRUCTURAL_PP,
            "dispersion_iqr_pp": THRESHOLD_DISPERSION_IQR_PP,
            "rank_disagreement_frac": THRESHOLD_RANK_DISAGREE_FRAC,
        },
        "donor_exclusions": {d: DONOR_EXCLUSIONS[d] for d in sorted(excluded)},
        "crosswalk": str(CROSSWALK.relative_to(REPO)),
        "crosswalk_is_generated": "run_vocab in phenocycler/integration/vocab.py writes it; never hand-edited here",
        "pln": {
            "comparable_lineages": pln_lin,
            "excluded_lineages": pln_dropped,
            "asymmetries": asymmetry_notes(cw, pln_lin),
            "protein_source": str(PLN_PHENO.relative_to(REPO)),
            "rna_source": str(PLN_XEN.relative_to(REPO)),
            "protein_donor_state": pln_p_state,
            "rna_donor_state": pln_r_state,
            "donors_compared": sorted(set(pln_long.donor_id)) if not pln_long.empty else [],
        },
        "pancreas": {
            "comparable_lineages": panc_lin_eff,
            "excluded_lineages": {**panc_dropped, **panc_unavail},
            "asymmetries": asymmetry_notes(cw, panc_lin_eff),
            "protein_source": str(PANC_PHENO.relative_to(REPO)),
            "rna_source": str(panc_root.relative_to(REPO)) if panc_root else None,
            "rna_notes": panc_notes,
            "donors_compared": sorted(set(panc_long.donor_id)) if not panc_long.empty else [],
        },
        "detection": "Xenium detection is lognorm >= log1p(1) as applied upstream by the calling "
                     "lane; no raw-count floor is introduced here.",
        "nucleus_rule": "A cell whose nucleus was MEASURED absent is excluded from every "
                        "calculation.  A cell whose nucleus was never measured is a missing "
                        "measurement, not a measured zero; it is retained and its section is "
                        "flagged by nucleus_measured_frac, because discarding it would record an "
                        "unavailable measurement as an absent nucleus.",
        "min_donors_for_verdict": MIN_DONORS_FOR_VERDICT,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")

    for title, df in (("pLN lineage summary", pln_sum), ("pancreas lineage summary", panc_sum)):
        log(f"\n=== {title} ===")
        if df.empty:
            log("  (nothing to compare)")
        else:
            log(df.to_string(index=False))
    for title, df in (("pLN batch eta-squared", pln_eta), ("pancreas batch eta-squared", panc_eta)):
        log(f"\n=== {title} ===")
        if df.empty:
            log("  (none)")
        else:
            log(df.to_string(index=False))
    log(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
