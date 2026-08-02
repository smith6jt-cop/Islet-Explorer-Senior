#!/usr/bin/env python3
"""
Merge broad-lineage cell types into the single-cell AnnData (2026-07 gap-fill).

No current script produced `data/singlecell_protein_phenotyped.h5ad`, yet every Phase-2
app-build script reads it (the old one came from the now-archived scimap typers). This
joins `data/phenotype/broad/` (assign_broad_lineage.py output) into
`data/singlecell_protein.h5ad` by `object_id` and writes `obs['phenotype'] =
broad_lineage` (the 7 broad classes: Epithelial / Fibroblast / Immune / Endocrine /
Endothelial / Muscle / Neural — Neutrophil/MPO folded into Immune). It ALSO bakes the REDSEA-corrected 58-marker MFI into `X` + `.raw`
(joined by object_id from data/cells_redsea) and drops the ~27.7k zero-area NaN edge cells,
making this file the single canonical, row-aligned input for scVI / aggregation / drill-down.

  IN : data/singlecell_protein.h5ad , data/phenotype/broad/donor_id=*/data_0.parquet ,
       data/cells_redsea/donor_id=*/data_0.parquet (corrected MFI)
  OUT: data/singlecell_protein_phenotyped.h5ad  (X = raw = REDSEA-corrected MFI)

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

    # --- bake REDSEA-corrected 58-marker MFI into X + .raw; drop zero-area NaN edge cells ---
    # This is the single canonical row set for the whole trajectory chain: the phenotyped h5ad,
    # the scVI input (raw/X), the per-cell latent, and aggregate_islets' row-order attach all key
    # off exactly these rows. train_scvi_senior.py / extract_islet_cells.py read raw/X + raw/var
    # (neither build set .raw before) and aggregate_islets reads a.raw.X — all now = corrected MFI.
    markers = list(adata.var_names)                       # 58 common; identical in cells_redsea
    red_files = sorted(glob.glob(str(SD / "cells_redsea" / "donor_id=*" / "data_0.parquet")))
    if not red_files:
        raise FileNotFoundError(f"no cells_redsea parquet under {SD/'cells_redsea'} — run redsea_full.py")
    cor = pd.concat([pd.read_parquet(f, columns=["object_id", *markers]) for f in red_files],
                    ignore_index=True)
    cor["object_id"] = cor["object_id"].astype(str)
    cor = cor.drop_duplicates("object_id").set_index("object_id")
    ridx = cor.index.get_indexer(pd.Index(oid))
    if (ridx < 0).any():                                  # cells ⊆ cells_redsea → expect no gap
        raise KeyError(f"{int((ridx < 0).sum())} cells absent from cells_redsea (object_id join gap)")
    Xc = cor.iloc[ridx][markers].to_numpy(np.float32)     # reindexed to var_names order
    finite = ~np.isnan(Xc).any(axis=1)                    # zero-area REDSEA edge cells -> NaN
    n_drop = int((~finite).sum())
    print(f"[redsea] baking corrected MFI into X + .raw; dropping {n_drop} NaN edge cells "
          f"({100 * n_drop / len(finite):.3f}%)")
    adata = adata[finite].copy()
    adata.X = Xc[finite]
    adata.raw = adata                                     # -> raw/X + raw/var (58 corrected markers)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(a.out)
    print("[done]", a.out, adata.shape)
    print("phenotype:", adata.obs["phenotype"].value_counts().to_dict())
    print("[ok] 0 Unassigned — every cell typed" if n_missing == 0
          else f"[check] {n_missing} Unassigned")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
