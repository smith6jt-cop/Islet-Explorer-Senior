#!/usr/bin/env python3
"""LEGACY LANE, NOT RUNNABLE HERE -- recorded 2026-09-18.

``scripts/retire/dispositions.py`` records this file as a verified orphan (code tier 1):
nothing in the tree imports it, and the inputs named below are absent from this checkout.
It is kept as the record of how that artifact was built, not as a description of anything
current.

The current PhenoCycler typing is the subpanel classification under
``core/phenocycler/scripts/typing/`` -- context from QuPath objects, then one question per
context against reference levels drawn from the ordinary tissue, fitted PER DONOR.
``data/development/typing/cells/`` holds 25,230,109 cells over 20 donors (2026-09-15 build,
eight classes plus ``Excluded``). Nothing below describes that method.

RESTORE WAS RETIRED 2026-08-14 and is on the never-revive list, together with the flat
``cells_redsea*`` / ``restore*`` roots. Every RESTORE ``_pos`` / positivity statement below is
part of that retired lane. The current typing never asks "is this cell <marker>-positive" as a
standalone question.

Train scVI on the Senior 59-plex single cells to produce the per-cell latent
adata.obsm['X_scVI'] that the islet trajectory needs (its per-islet mean,
X_scVI_mean, is the load-bearing input to compute_trajectory_senior.py's DPT).

BATCH-CORRECTION FIX (2026-07): the trajectory pseudotime was found to be a
donor/slide axis, not disease progression (batch eta^2 = 0.988 on the deployed
data). Root cause: the Junior recipe (SCVI_CODEX_v2.ipynb) set only
categorical_covariate_keys=["Age","Gender"] and NO batch_key, so scVI's batch
decoder never activated and per-image staining/gain drift leaked into the latent.
We now:
  - set batch_key="imageid" (activates the batch decoder — the actual fix),
  - drop Age/Gender covariates (each donor has one Age/Sex → collinear with
    imageid, redundant and prone to degenerate fits once batch is modelled),
  - keep dispersion='gene-batch' (now meaningful with >1 batch).

INPUT MUST BE REDSEA+RESTORE-CORRECTED MFI. scVI's NB likelihood + batch_key
handle per-image intensity drift, but pixel-level lateral spillover must already
be removed upstream (REDSEA) so the batch key is not doing double duty. Supply
the corrected raw-scale means either as PHENO's raw/X, or overlay them from the
REDSEA parquet with --corrected-mfi-dir data/cells_redsea (joined by object_id).

Recipe (unchanged otherwise): input = (corrected) MFI as the scVI "counts" layer;
SCVI(n_latent=10, n_layers=2, n_hidden=128, dropout_rate=0.1, gene_likelihood='nb');
train(max_epochs=400, batch_size=2048, early_stopping, patience=45, lr=1e-3);
obsm['X_scVI'] = get_latent_representation().

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


def _overlay_corrected_mfi(a, obj_ids, corrected_dir, markers):
    """Replace a.X (and the 'counts' layer) with REDSEA-corrected MFI joined by
    object_id from data/cells_redsea/donor_id=*/*.parquet. Off unless requested;
    raises loudly on any misalignment so a bad join can never silently ship."""
    corrected_dir = Path(corrected_dir)
    files = sorted(corrected_dir.glob("donor_id=*/*.parquet"))
    if not files:
        raise FileNotFoundError(f"--corrected-mfi-dir {corrected_dir} has no donor_id=*/*.parquet")
    if obj_ids is None:
        raise KeyError("PHENO obs lacks 'object_id' — cannot join corrected MFI by UUID")
    frames = [pd.read_parquet(fp) for fp in files]
    cor = pd.concat(frames, ignore_index=True)
    missing = [m for m in markers if m not in cor.columns]
    if missing:
        raise KeyError(f"corrected parquet missing markers: {missing[:5]}... ({len(missing)} total)")
    cor = cor.drop_duplicates("object_id").set_index("object_id")
    idx = cor.index.get_indexer(pd.Index(obj_ids))
    if (idx < 0).any():
        raise KeyError(f"{int((idx < 0).sum())} cells have no corrected MFI (object_id join gap)")
    Xc = cor.iloc[idx][markers].to_numpy(dtype=np.float32)
    if np.isnan(Xc).any():
        # cells_redsea carries a few NaN edge-cells per donor (documented gotcha);
        # scVI cannot ingest NaN — fail rather than impute silently.
        raise ValueError(f"{int(np.isnan(Xc).any(axis=1).sum())} corrected cells contain NaN "
                         "(REDSEA edge cells) — filter upstream before training")
    a.X = Xc
    a.layers["counts"] = a.X
    print(f"[load] overlaid REDSEA-corrected MFI from {corrected_dir} ({len(files)} donors)")
    return a


def load_anndata(smoke=False, corrected_mfi_dir=None):
    f = h5py.File(PHENO, "r")
    vg = f["raw/var"]
    vidx = vg.attrs.get("_index", "_index"); vidx = vidx.decode() if isinstance(vidx, bytes) else vidx
    markers = [v.decode() if isinstance(v, bytes) else v for v in vg[vidx][()]]
    obs = f["obs"]
    age = np.asarray(obs["Age"][()], dtype=float)
    sex = cat_col(obs["Sex"]) if isinstance(obs["Sex"], h5py.Group) else np.asarray(obs["Sex"][()]).astype(str)
    imageid = cat_col(obs["imageid"]) if isinstance(obs["imageid"], h5py.Group) else np.asarray(obs["imageid"][()]).astype(str)
    obj_ids = None
    if "object_id" in obs:
        obj_ids = (cat_col(obs["object_id"]) if isinstance(obs["object_id"], h5py.Group)
                   else np.asarray(obs["object_id"][()]).astype(str))
    n = f["raw/X"].shape[0]
    if smoke:
        idx = np.linspace(0, n - 1, 150_000, dtype=int)
        X = f["raw/X"][:][idx]                      # read all then subsample (simple, RAM ok)
        age, sex, imageid = age[idx], sex[idx], imageid[idx]
        if obj_ids is not None:
            obj_ids = obj_ids[idx]
    else:
        X = f["raw/X"][:]                           # 23.28M x 59 (corrected) MFI (~5.5 GB)
    f.close()
    obs_df = pd.DataFrame({
        "Age": age,                                  # kept for reference (NOT a covariate — collinear with imageid)
        "Gender": pd.Categorical(sex),               # kept for reference (NOT a covariate — collinear with imageid)
        "imageid": pd.Categorical(imageid),          # batch_key: activates scVI's batch decoder (the fix)
    })
    a = ad.AnnData(X=np.asarray(X, dtype=np.float32), obs=obs_df)
    a.var_names = markers
    a.layers["counts"] = a.X                          # scVI input = (corrected) MFI as "counts"
    if corrected_mfi_dir is not None:
        a = _overlay_corrected_mfi(a, obj_ids, corrected_mfi_dir, markers)
    print(f"[load] AnnData {a.shape} | MFI range [{a.X.min():.1f}, {a.X.max():.1f}] "
          f"median {np.median(a.X):.1f} | donors {obs_df['imageid'].nunique()} | ages {obs_df['Age'].nunique()}")
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="150k cells, 2 epochs (validate API/GPU)")
    ap.add_argument("--max-epochs", type=int, default=400)
    ap.add_argument("--corrected-mfi-dir", default=None,
                    help="overlay REDSEA-corrected MFI from this dir (e.g. data/cells_redsea), "
                         "joined by object_id; omit if PHENO's raw/X is already corrected")
    args = ap.parse_args()

    scvi.settings.seed = 0
    print(f"[scvi] version {scvi.__version__} | smoke={args.smoke}")
    a = load_anndata(smoke=args.smoke, corrected_mfi_dir=args.corrected_mfi_dir)

    scvi.model.SCVI.setup_anndata(
        a, layer="counts",
        batch_key="imageid",          # activate the batch decoder — the batch-correction fix
        categorical_covariate_keys=None,   # Age/Gender dropped: collinear with imageid, now redundant
        continuous_covariate_keys=None,
    )
    model = scvi.model.SCVI(
        a, n_latent=10, n_layers=2, n_hidden=128, dropout_rate=0.1,
        dispersion="gene-batch", gene_likelihood="nb",   # gene-batch now meaningful (>1 batch)
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
