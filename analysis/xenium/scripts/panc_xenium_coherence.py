#!/usr/bin/env python3
"""Where the naming evidence of a lineage ended up, per section.

The release already records whether a lineage received ZERO cells while its own naming genes
were detected.  That test is too coarse on its own: a section that gives Immune_B_Plasma 24 cells
while MS4A1 is detected in ~590 of its cells has the same defect as one that gives it none.

This asks the finer question, still entirely inside one section: for each naming gene, how many of
the cells that DETECT it sit in a subcluster called that lineage, and where the rest sit.  Every
quantity is that section against itself.  There is no expected fraction anywhere in this file, no
target, and nothing here gates or tunes a call -- it says where the evidence went so that a
scoring or resolution failure is visible instead of being reported as biology.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

#: gene -> the lineage whose naming set it belongs to
EVIDENCE_OF = {
    "CD3E": "Immune_T", "CD3G": "Immune_T", "CD8A": "Immune_T", "CD2": "Immune_T",
    "CD247": "Immune_T", "LCK": "Immune_T", "ZAP70": "Immune_T",
    "MS4A1": "Immune_B_Plasma", "CD79A": "Immune_B_Plasma", "CD19": "Immune_B_Plasma",
    "NKG7": "Immune_NK_ILC", "KLRD1": "Immune_NK_ILC", "KLRF1": "Immune_NK_ILC",
    "CD68": "Myeloid", "CSF1R": "Myeloid", "ITGAM": "Myeloid",
    "FCGR3B": "Granulocyte", "CSF3R": "Granulocyte", "ELANE": "Granulocyte",
}


def per_section(work: Path, section: str) -> list[dict]:
    gs = pd.read_parquet(work / f"{section}.gene_stats.parquet",
                         columns=["cluster_id", "gene", "det_frac"])
    sub = pd.read_parquet(work / f"{section}.subclusters.parquet",
                          columns=["subcluster_id", "n_cells", "lineage", "median_depth"])
    sub = sub.rename(columns={"subcluster_id": "cluster_id"})
    gs = gs[gs.gene.isin(EVIDENCE_OF)].merge(sub, on="cluster_id", how="inner")
    gs["n_detecting"] = gs.det_frac * gs.n_cells
    out = []
    for gene, blk in gs.groupby("gene"):
        home = EVIDENCE_OF[gene]
        tot = float(blk.n_detecting.sum())
        by = blk.groupby("lineage").n_detecting.sum().sort_values(ascending=False)
        out.append({
            "section_key": section, "gene": gene, "evidence_for": home,
            "n_cells_detecting": int(round(tot)),
            "n_in_that_lineage": int(round(float(by.get(home, 0.0)))),
            "frac_in_that_lineage": round(float(by.get(home, 0.0)) / tot, 4) if tot else None,
            "n_cells_called_that_lineage": int(sub.loc[sub.lineage == home, "n_cells"].sum()),
            "where_the_rest_went": {k: int(round(v)) for k, v in by.head(4).items()},
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="/data/xenium_typing/panc_2026-09-17/_work")
    ap.add_argument("--out", default="/data/xenium_typing/panc_2026-09-17/release/"
                                     "coherence_where_the_evidence_went.csv")
    args = ap.parse_args()
    work = Path(args.work)
    rows = []
    for rp in sorted(work.glob("*.report.json")):
        rows += per_section(work, rp.name.replace(".report.json", ""))
    frame = pd.DataFrame(rows)
    frame.to_csv(args.out, index=False)
    pd.set_option("display.width", 220)
    key = frame[frame.gene.isin(["CD3E", "MS4A1", "CD68", "NKG7", "FCGR3B"])]
    piv = key.pivot(index="section_key", columns="gene", values="frac_in_that_lineage")
    print("fraction of the cells detecting each naming gene that sit in that gene's lineage")
    print(piv.to_string())
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
