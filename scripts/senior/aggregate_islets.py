#!/usr/bin/env python3
"""
3B (Islet Explorer Senior): aggregate phenotyped single cells -> islet-level
datasets, reusing Junior's create_islet_dataset_fixed (copied alongside this file).

Run in the `scimap` conda env (anndata 0.11). Input is the phenotyped single-cell
h5ad; expression aggregated is the RAW MFI (from .raw), grouped by imageid+Islet_N
(core) / _exp20um (peri), require_paired=True (same convention as Junior's 5,214).

Outputs:
  ~/IO60panc2nd/senior_data/islets_core_senior.h5ad
  ~/IO60panc2nd/senior_data/islets_peri_senior.h5ad
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fixed_islet_aggregation import create_islet_dataset_fixed

BASE = Path.home() / "IO60panc2nd" / "senior_data"
PHENO = BASE / "singlecell_protein_phenotyped.h5ad"

# map Senior's morphology column names -> the names create_islet_dataset_fixed expects
RENAME = {
    "nucleus_area": "Nucleus Area", "nucleus_length": "Nucleus Length",
    "nucleus_circularity": "Nucleus Circularity", "nucleus_solidity": "Nucleus Solidity",
    "nucleus_maxdiam": "Nucleus: Max diameter µm", "nucleus_mindiam": "Nucleus: Min diameter µm",
    "cell_area": "Cell Area", "cell_length": "Cell Length",
    "cell_circularity": "Cell Circularity", "cell_solidity": "Cell Solidity",
    "cell_maxdiam": "Cell: Max diameter µm", "cell_mindiam": "Cell: Min diameter µm",
    "dist_tissue": "Dist to Closest Tissue",
}


def main():
    a = ad.read_h5ad(PHENO)
    raw = a.raw.X if a.raw is not None else a.X
    if not isinstance(raw, np.ndarray):
        raw = raw.toarray()
    sub = ad.AnnData(X=np.asarray(raw, dtype=np.float32), obs=a.obs.copy(), var=a.var.copy())
    sub.obs = sub.obs.rename(columns=RENAME)
    if "Gender" not in sub.obs.columns and "Sex" in sub.obs.columns:
        sub.obs["Gender"] = sub.obs["Sex"].values
    print(f"[3B] aggregating {sub.n_obs:,} cells x {sub.n_vars} markers (raw MFI means)")

    # Attach the per-cell scVI latent (trained by train_scvi_senior.py, row-aligned to
    # the phenotyped h5ad) so create_islet_dataset_fixed averages it into X_scVI_mean —
    # the load-bearing trajectory input (rebuild_trajectory.ipynb).
    latent = BASE / "X_scVI_singlecell.npy"
    if latent.exists():
        z = np.load(latent)
        assert z.shape[0] == sub.n_obs, f"latent rows {z.shape[0]} != cells {sub.n_obs} (re-train)"
        sub.obsm["X_scVI"] = z.astype(np.float32)
        print(f"[3B] attached X_scVI {z.shape} -> X_scVI_mean will be computed per islet")
    else:
        print(f"[3B] WARNING: {latent} not found — no X_scVI_mean (trajectory cannot be fixed)")

    core, peri = create_islet_dataset_fixed(sub, region="separate",
                                            min_cells=0, require_paired=True)
    core.write_h5ad(BASE / "islets_core_senior.h5ad")
    peri.write_h5ad(BASE / "islets_peri_senior.h5ad")
    print(f"\n[3B] core islets: {core.n_obs} | peri islets: {peri.n_obs}")
    print("[3B] phenotypes:", core.uns.get("phenotype_names"))
    if "prop_Beta cell" in core.obs.columns:
        print("[3B] mean Beta-cell proportion (core) by donor status:")
        print(core.obs.groupby("donor_status", observed=True)["prop_Beta cell"].mean().round(3).to_string())
    if "prop_Alpha cell" in core.obs.columns:
        print("[3B] mean Alpha-cell proportion (core) by donor status:")
        print(core.obs.groupby("donor_status", observed=True)["prop_Alpha cell"].mean().round(3).to_string())
    print("[done]", BASE / "islets_core_senior.h5ad")


if __name__ == "__main__":
    main()
