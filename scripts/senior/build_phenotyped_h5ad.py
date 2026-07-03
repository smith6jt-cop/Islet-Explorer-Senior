#!/usr/bin/env python3
"""
Merge broad-lineage cell types into the single-cell AnnData (2026-07 gap-fill).

No current script produced `data/singlecell_protein_phenotyped.h5ad`, yet every Phase-2
app-build script reads it (the old one came from the now-archived scimap typers). This
joins `data/phenotype/broad/` (assign_broad_lineage.py output) into
`data/singlecell_protein.h5ad` by `object_id` and writes `obs['phenotype'] =
broad_lineage` (the 6 broad classes: Epithelial / Fibroblast / Immune / Endocrine /
Endothelial / Muscle). X (58-common raw MFI) is unchanged.

  IN : data/singlecell_protein.h5ad , data/phenotype/broad/donor_id=*/data_0.parquet
  OUT: data/singlecell_protein_phenotyped.h5ad

⚠️ POSITIVITY CONVENTION: `compute_islet_markers.py` and `build_composition_and_tissue.py`
currently threshold `X > 0.5` (the old scimap-rescaled 0-1 X). Here X is RAW MFI, so `>0.5`
is NOT a valid positivity call and those scripts must be reconciled to the RESTORE `_pos`
convention (join `data/restore_gated_redsea/`) — RESTORE covers 10 markers; a per-marker
threshold is needed for the rest. Tracked as a Phase-2 follow-up (see the plan / scripts/CLAUDE.md).

Run in `scvi-env` (or `scimap`).
"""
from __future__ import annotations
import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

SD = Path(__file__).resolve().parents[2] / "data"
IN = SD / "singlecell_protein.h5ad"
BROAD = SD / "phenotype" / "broad"
OUT = SD / "singlecell_protein_phenotyped.h5ad"


def load_broad(broad_dir: Path) -> pd.DataFrame:
    files = sorted(glob.glob(str(broad_dir / "donor_id=*" / "data_0.parquet")))
    if not files:
        raise FileNotFoundError(f"no broad-lineage parquet under {broad_dir} — run assign_broad_lineage.py")
    cols = ["object_id", "broad_lineage", "epi_default", "assign_margin"]
    df = pd.concat([pd.read_parquet(f, columns=cols) for f in files], ignore_index=True)
    df["object_id"] = df["object_id"].astype(str)
    return df.drop_duplicates("object_id").set_index("object_id")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", type=Path, default=IN)
    ap.add_argument("--broad", type=Path, default=BROAD)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)

    adata = ad.read_h5ad(a.inp)
    oid = (adata.obs["Object ID"].astype(str).to_numpy() if "Object ID" in adata.obs
           else adata.obs_names.astype(str).to_numpy())
    broad = load_broad(a.broad)

    idx = broad.index.get_indexer(pd.Index(oid))
    found = idx >= 0
    n_missing = int((~found).sum())
    if n_missing:
        # assign_broad_lineage types every cell (0 "Unassigned"); a gap means the broad
        # set is stale vs this h5ad (e.g. rebuilt cells but not re-typed).
        print(f"[warn] {n_missing}/{len(oid)} cells have no broad label -> 'Unassigned' "
              "(broad output is stale vs singlecell_protein.h5ad; re-run assign_broad_lineage.py)")
    safe = idx.clip(min=0)
    lineage = np.where(found, broad["broad_lineage"].to_numpy()[safe], "Unassigned")
    adata.obs["phenotype"] = pd.Categorical(lineage)
    adata.obs["epi_default"] = np.where(found, broad["epi_default"].to_numpy()[safe], False)
    adata.obs["assign_margin"] = np.where(found, broad["assign_margin"].to_numpy()[safe].astype(float), np.nan)
    adata.uns["phenotype_scheme"] = "broad_lineage_step1"

    a.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(a.out)
    print("[done]", a.out, adata.shape)
    print("phenotype:", adata.obs["phenotype"].value_counts().to_dict())
    print("[ok] 0 Unassigned — every cell typed" if n_missing == 0
          else f"[check] {n_missing} Unassigned")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
