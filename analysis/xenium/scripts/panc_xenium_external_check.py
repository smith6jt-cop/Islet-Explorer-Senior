#!/usr/bin/env python3
"""External check: the new pancreas calls against the existing 2-donor scANVI release.

It is EXTERNAL because nothing in this lane was fitted to it, tuned against it or gated on it.
It is not a verification target and it is not evidence that either lane is right: the two lanes
type different cells by different means (scANVI on a per-donor model there, within-section
cluster naming here) from different h5ads of the same sections, and their label spaces differ.
The agreement is reported as a number, with the confusion that produced it, and nothing is
adjusted on the strength of it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OLD = REPO / "data/development/xenium_typing/release"

#: The old space collapses B and plasma into Immune_B and both stromal cores into Stromal, and
#: carries no granulocyte and no NK label at all.  Mapping the NEW space onto the OLD one is the
#: only direction that loses no distinction the old release ever made.
TO_OLD = {
    "Endocrine": "Endocrine", "Exocrine": "Exocrine",
    "Immune_T": "Immune_T", "Immune_NK_ILC": "Immune_T",
    "Immune_B_Plasma": "Immune_B", "Myeloid": "Myeloid",
    "Granulocyte": "Myeloid",
    "Vascular": "Vascular",
    "Stromal_Pericyte_SMC": "Stromal", "Stromal_Fibroblast_Stellate": "Stromal",
    "Neural": "Neural",
}
#: why each collapse: the old lane had no such label, so a cell of that lineage could only have
#: landed in the listed one
WHY = {
    "Immune_NK_ILC->Immune_T": "the old space has no NK label; T_NK is the shared parent",
    "Granulocyte->Myeloid": "the old space has no granulocyte label; both are myeloid-lineage",
    "Stromal_Pericyte_SMC->Stromal": "the old space collapses both stromal cores",
    "Stromal_Fibroblast_Stellate->Stromal": "the old space collapses both stromal cores",
    "Immune_B_Plasma->Immune_B": "the old label is Immune_B",
}


def check(new_release: Path, donors: list[str]) -> dict:
    out = {"what_this_is": __doc__.strip().splitlines()[2], "mapping": TO_OLD,
           "why_each_collapse": WHY, "donors": {}}
    for donor in donors:
        a = pd.read_parquet(new_release / f"calls/donor_id={donor}/roi=panc/data.parquet",
                            columns=["cell_id", "lineage"])
        b = pd.read_parquet(OLD / f"calls/donor_id={donor}/roi=panc/data.parquet",
                            columns=["cell_id", "lineage_native"])
        m = a.merge(b, on="cell_id", how="inner")
        m["new_in_old_space"] = m["lineage"].map(TO_OLD)
        both = m[m["new_in_old_space"].notna()
                 & ~m["lineage_native"].isin(["Unavailable"])].copy()
        agree = float((both["new_in_old_space"] == both["lineage_native"]).mean())
        conf = (both.groupby(["new_in_old_space", "lineage_native"]).size()
                .rename("n").reset_index().sort_values("n", ascending=False))
        per_lab = (both.assign(ok=both["new_in_old_space"] == both["lineage_native"])
                   .groupby("new_in_old_space").agg(n=("ok", "size"), agree=("ok", "mean")))
        out["donors"][donor] = {
            "n_cells_new": len(a), "n_cells_old": len(b),
            "n_cell_ids_in_both": len(m),
            "n_comparable": len(both),
            "n_new_uncalled_or_no_nucleus": int(m["new_in_old_space"].isna().sum()),
            "n_old_unavailable": int((m["lineage_native"] == "Unavailable").sum()),
            "agreement_on_comparable_cells": round(agree, 4),
            "agreement_by_label": {k: {"n": int(v.n), "agree": round(float(v.agree), 4)}
                                   for k, v in per_lab.iterrows()},
            "top_confusions": [
                {"new": r.new_in_old_space, "old": r.lineage_native, "n": int(r.n)}
                for r in conf[conf.new_in_old_space != conf.lineage_native].head(12)
                .itertuples()],
            "new_label_counts_native": {k: int(v) for k, v in
                                        a.lineage.value_counts().items()},
            "old_label_counts": {k: int(v) for k, v in b.lineage_native.value_counts().items()},
        }
    return out


#: raw count >= 1 on these is direct evidence a cell is a T cell; every one is on the 5K panel
T_EVIDENCE = ["CD3E", "CD3G", "CD2", "CD8A", "CD247", "LCK"]


def t_evidence(new_release: Path, donors: list[str]) -> dict:
    """How much direct T evidence each lane's T calls carry, on the same cells.

    This does not adjudicate between the lanes and nothing is adjusted on it. The two label sets
    disagree most about T cells, and this says what is measured in the cells they disagree about.
    """
    import anndata as ad
    import scipy.sparse as sp

    out = {}
    for donor in donors:
        h5 = Path("/data/xenium/h5ad_out") / f"{donor}_PancXenium.h5ad"
        adata = ad.read_h5ad(h5)
        ids = adata.obs_names.to_numpy().astype(str)
        var = list(adata.var_names)
        X = adata.X.tocsc() if sp.issparse(adata.X) else sp.csc_matrix(adata.X)
        det = pd.DataFrame(
            {g: (np.asarray(X[:, var.index(g)].todense()).ravel() >= 1)
             for g in T_EVIDENCE + ["PTPRC"] if g in var}, index=ids)
        del X, adata
        det["any_T_gene"] = det[[g for g in T_EVIDENCE if g in det]].any(axis=1)
        new = pd.read_parquet(new_release / f"calls/donor_id={donor}/roi=panc/data.parquet",
                              columns=["cell_id", "lineage"])
        old = pd.read_parquet(OLD / f"calls/donor_id={donor}/roi=panc/data.parquet",
                              columns=["cell_id", "lineage_native"])
        m = new.merge(old, on="cell_id").set_index("cell_id").join(det, how="inner")

        def block(sel, name):
            s = m[sel]
            return {"set": name, "n": int(len(s)),
                    **{g: round(float(s[g].mean()), 4) for g in det.columns}}

        out[donor] = [
            block(m.lineage == "Immune_T", "called Immune_T by this lane"),
            block(m.lineage_native == "Immune_T", "called Immune_T by the scANVI lane"),
            block((m.lineage == "Exocrine") & (m.lineage_native == "Immune_T"),
                  "Exocrine here and Immune_T there"),
            block((m.lineage == "Immune_T") & (m.lineage_native == "Immune_T"),
                  "both lanes call Immune_T"),
            block(m.index.notna(), "the whole section"),
        ]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", default="/data/xenium_typing/panc_2026-09-17/release")
    ap.add_argument("--donors", nargs="*", default=["6436", "6442"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--evidence", action="store_true",
                    help="also measure T evidence per lane; reads each donor's h5ad")
    args = ap.parse_args()
    res = check(Path(args.release), args.donors)
    if args.evidence:
        res["t_evidence_in_the_cells_the_lanes_disagree_about"] = t_evidence(
            Path(args.release), args.donors)
    text = json.dumps(res, indent=1)
    print(text)
    dest = Path(args.out) if args.out else Path(args.release) / "external_check_vs_scanvi.json"
    dest.write_text(text)


if __name__ == "__main__":
    main()
