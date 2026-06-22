#!/usr/bin/env python3
"""
Export Step-1 broad lineage as tiny per-image CSVs keyed by QuPath's detection UUID, for the
companion `scripts/groovy/import_broad_lineage.groovy` (UUID match → set PathClass directly).

This is more robust than the CytoMAP centroid-matching importer: QuPath's PathObject.getID()
is exactly our `object_id`, so matching is exact (no centroid rounding) and the class is set
in one step (no separate measurement→class conversion).

  python scripts/senior/export_broad_class_for_qupath.py
  -> data/phenotype/qupath_class/pheno_class_<donor>.csv   (object_id, broad_lineage, image)
"""
from __future__ import annotations
import glob
from pathlib import Path
import pandas as pd

BASE = Path(__file__).resolve().parents[2] / "data"
OUT = BASE / "phenotype" / "qupath_class"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    total = 0
    for bf in sorted(glob.glob(str(BASE / "phenotype" / "broad" / "donor_id=*" / "*.parquet"))):
        donor = bf.split("donor_id=")[1].split("/")[0]
        b = pd.read_parquet(bf, columns=["object_id", "broad_lineage"])
        img = pd.read_parquet(sorted(glob.glob(str(BASE / "cells" / f"donor_id={donor}" / "*.parquet")))[0],
                              columns=["object_id", "image"])
        out = b.merge(img.astype({"object_id": str}), on="object_id", how="inner")
        out = out[["object_id", "broad_lineage", "image"]]
        csv = OUT / f"pheno_class_{donor}.csv"
        out.to_csv(csv, index=False)
        total += len(out)
        print(f"[{donor}] {len(out):,} cells -> {csv.name}", flush=True)
    print(f"\n[done] {total:,} cells across {len(glob.glob(str(OUT/'pheno_class_*.csv')))} images -> {OUT}")


if __name__ == "__main__":
    raise SystemExit(main())
