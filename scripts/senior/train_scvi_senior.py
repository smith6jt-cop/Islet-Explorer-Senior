#!/usr/bin/env python3
"""
Train scVI on the Senior 59-plex single cells, reproducing the Junior recipe in
single_cell_analysis/SCVI_CODEX_v2.ipynb EXACTLY, to produce the per-cell latent
adata.obsm['X_scVI'] that the islet trajectory needs (its per-islet mean,
X_scVI_mean, is the load-bearing input to rebuild_trajectory.ipynb's DPT).

Junior recipe (cells 8-25):
  - input = raw MFI as the scVI "counts" layer
  - setup_anndata(layer="counts", categorical_covariate_keys=["Age","Gender"])  # NO batch_key
  - SCVI(n_latent=10, n_layers=2, n_hidden=128, dropout_rate=0.1,
         dispersion='gene-batch', gene_likelihood='nb')
  - train(max_epochs=400, batch_size=2048, early_stopping, patience=45, lr=1e-3)
  - obsm['X_scVI'] = get_latent_representation()

Deviation (justified): Senior cells were already QC'd upstream (build_cells_parquet
--min-cell-area 5, NaN dropped), so the notebook's extra per-cell QC filtering is
SKIPPED to keep the latent row-aligned 1:1 with the 23.28M-cell phenotyped h5ad
(so it can be attached by row order during aggregation). CLR is viz-only in the
notebook and unused by scVI, so it is omitted.

Reads via h5py (the 23 GB file is anndata-0.11 format; scvi-env's anndata is older).
Writes the latent to data/X_scVI_singlecell.npy (+ saves the model dir).

Run in the `scvi-env` conda env (scvi-tools 1.4, torch+cuda). GPU: RTX 4000 Ada.
    python scripts/senior/train_scvi_senior.py --smoke   # 150k cells, 2 epochs (API/GPU check)
    python scripts/senior/train_scvi_senior.py           # full 23.28M-cell, 400-epoch run
"""
import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

SD = Path(__file__).resolve().parents[2] / "data"
PHENO = SD / "singlecell_protein_phenotyped.h5ad"
LATENT_OUT = SD / "X_scVI_singlecell.npy"
MODEL_DIR = SD / "scvi_model_senior"


def cat_col(grp):
    cats = grp["categories"][()]
    cats = np.array([c.decode() if isinstance(c, (bytes, bytearray)) else c for c in cats], dtype=object)
    codes = grp["codes"][()]
    out = np.empty(codes.shape, dtype=object)
    out[codes >= 0] = cats[codes[codes >= 0]]
    out[codes < 0] = ""
    return out


def load_anndata(smoke=False):
    f = h5py.File(PHENO, "r")
    vg = f["raw/var"]
    vidx = vg.attrs.get("_index", "_index"); vidx = vidx.decode() if isinstance(vidx, bytes) else vidx
    markers = [v.decode() if isinstance(v, bytes) else v for v in vg[vidx][()]]
    obs = f["obs"]
    age = np.asarray(obs["Age"][()], dtype=float)
    sex = cat_col(obs["Sex"]) if isinstance(obs["Sex"], h5py.Group) else np.asarray(obs["Sex"][()]).astype(str)
    imageid = cat_col(obs["imageid"]) if isinstance(obs["imageid"], h5py.Group) else np.asarray(obs["imageid"][()]).astype(str)
    n = f["raw/X"].shape[0]
    if smoke:
        idx = np.linspace(0, n - 1, 150_000, dtype=int)
        X = f["raw/X"][:][idx]                      # read all then subsample (simple, RAM ok)
        age, sex, imageid = age[idx], sex[idx], imageid[idx]
    else:
        X = f["raw/X"][:]                           # 23.28M x 59 raw MFI (~5.5 GB)
    f.close()
    obs_df = pd.DataFrame({
        "Age": age,                                  # categorical covariate (Junior passed Age as categorical)
        "Gender": pd.Categorical(sex),               # Sex -> Gender (Junior covariate name)
        "imageid": pd.Categorical(imageid),          # kept for reference; NOT a batch_key (matches notebook)
    })
    a = ad.AnnData(X=np.asarray(X, dtype=np.float32), obs=obs_df)
    a.var_names = markers
    a.layers["counts"] = a.X                          # scVI input = raw MFI as "counts"
    print(f"[load] AnnData {a.shape} | raw MFI range [{a.X.min():.1f}, {a.X.max():.1f}] "
          f"median {np.median(a.X):.1f} | donors {obs_df['imageid'].nunique()} | ages {obs_df['Age'].nunique()}")
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="150k cells, 2 epochs (validate API/GPU)")
    ap.add_argument("--max-epochs", type=int, default=400)
    args = ap.parse_args()

    scvi.settings.seed = 0
    print(f"[scvi] version {scvi.__version__} | smoke={args.smoke}")
    a = load_anndata(smoke=args.smoke)

    scvi.model.SCVI.setup_anndata(
        a, layer="counts",
        categorical_covariate_keys=["Age", "Gender"],   # NO batch_key (matches SCVI_CODEX_v2.ipynb)
        continuous_covariate_keys=None,
    )
    model = scvi.model.SCVI(
        a, n_latent=10, n_layers=2, n_hidden=128, dropout_rate=0.1,
        dispersion="gene-batch", gene_likelihood="nb",
    )
    print(model)
    model.train(
        max_epochs=2 if args.smoke else args.max_epochs,
        batch_size=2048,
        train_size=0.9,
        early_stopping=True, early_stopping_patience=45, check_val_every_n_epoch=1,
        accelerator="gpu", devices=1,
        plan_kwargs={"lr": 1e-3},
        datasplitter_kwargs={"num_workers": 0 if args.smoke else 8, "pin_memory": True},
    )
    latent = model.get_latent_representation()
    print(f"[scvi] latent {latent.shape}")

    if args.smoke:
        print("[smoke] OK — API + GPU + covariates + NB-on-MFI all work. Not saving.")
        return
    np.save(LATENT_OUT, latent.astype(np.float32))
    model.save(str(MODEL_DIR), overwrite=True)
    print(f"[done] saved latent -> {LATENT_OUT} | model -> {MODEL_DIR}")


if __name__ == "__main__":
    main()
