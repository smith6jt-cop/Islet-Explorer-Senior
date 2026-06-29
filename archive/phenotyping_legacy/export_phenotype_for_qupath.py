#!/usr/bin/env python3
"""
Export the Step-1 broad-lineage phenotype as per-image CSVs for visual inspection in QuPath,
in the format consumed by `scripts/groovy/import clusters and cell neighborhoods.groovy`
(CytoMAP-style): column 0 = Centroid X µm, column 1 = Centroid Y µm, column 2 = Image name,
columns 3+ = NUMERIC measurements matched to each detection by centroid.

The importer matches cells by `round(Centroid µm) == X/Y` — and data/cells centroids are already
stored at QuPath's round() precision (<1000 → 2 dp, ≥1000 → 1 dp), so they match directly.

Per cell we write: broad_lineage (integer code 1–6, see legend), the six lineage scores
(RESTORE _norm), assign_margin, and epi_default (1 = sub-threshold acinar default). All numeric
so the Groovy `as double` import works; the phenotype string→code legend is written alongside.

  python scripts/senior/export_phenotype_for_qupath.py
  -> data/phenotype/qupath_import/pheno_broad_<donor>.csv  (one per image)
  -> data/phenotype/qupath_import/broad_lineage_legend.csv
In QuPath: run the import groovy, select the qupath_import folder. Color/classify detections by
the `broad_lineage` measurement (1=Epithelial … 6=Muscle) to inspect against the image.
"""
from __future__ import annotations
import glob
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2] / "data"
OUT = BASE / "phenotype" / "qupath_import"
LNAMES = ["Epithelial", "Fibroblast", "Immune", "Endocrine", "Endothelial", "Muscle"]
CODE = {L: i + 1 for i, L in enumerate(LNAMES)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"code": list(CODE.values()), "broad_lineage": list(CODE.keys())}).to_csv(
        OUT / "broad_lineage_legend.csv", index=False)
    score_cols = [f"score_{L}" for L in LNAMES]
    total = 0
    for bf in sorted(glob.glob(str(BASE / "phenotype" / "broad" / "donor_id=*" / "*.parquet"))):
        donor = bf.split("donor_id=")[1].split("/")[0]
        b = pd.read_parquet(bf, columns=["object_id", "broad_lineage", "assign_margin",
                                         "epi_default"] + score_cols)
        c = pd.read_parquet(sorted(glob.glob(str(BASE / "cells" / f"donor_id={donor}" / "*.parquet")))[0],
                            columns=["object_id", "X_centroid", "Y_centroid", "image"])
        m = c.merge(b, on="object_id", how="inner")
        out = pd.DataFrame({
            "Centroid X µm": m["X_centroid"].to_numpy(),
            "Centroid Y µm": m["Y_centroid"].to_numpy(),
            "Image": m["image"].astype(str).to_numpy(),
            "broad_lineage": m["broad_lineage"].map(CODE).astype(int).to_numpy(),
            "assign_margin": m["assign_margin"].astype(float).round(4).to_numpy(),
            "epi_default": m["epi_default"].astype(int).to_numpy(),
        })
        for L in LNAMES:
            out[f"score_{L}"] = m[f"score_{L}"].astype(float).round(4).to_numpy()
        csv = OUT / f"pheno_broad_{donor}.csv"
        out.to_csv(csv, index=False)
        total += len(out)
        vc = m["broad_lineage"].value_counts()
        print(f"[{donor}] {len(out):,} cells -> {csv.name} | " +
              " ".join(f"{L[:3]}={int(vc.get(L,0)):,}" for L in LNAMES), flush=True)
    print(f"\n[done] {total:,} cells across {len(glob.glob(str(OUT/'pheno_broad_*.csv')))} images -> {OUT}")
    print("legend:", {v: k for k, v in CODE.items()})


if __name__ == "__main__":
    raise SystemExit(main())
