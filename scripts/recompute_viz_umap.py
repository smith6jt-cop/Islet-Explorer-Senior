#!/usr/bin/env python
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

Recompute visualization UMAP from raw marker expression.

The trajectory UMAP (from X_scVI_mean) shows poor disease-stage separation
because scVI's Age+Gender covariates correct out donor-level variation.
This script computes a separate UMAP from raw marker PCA, which preserves
INS/GCG/immune differences that distinguish ND vs Aab+ vs T1D.

Reads:  data/adata_ins_root.h5ad  (trajectory H5AD with X_scVI_mean UMAP)
Writes: data/adata_ins_root.h5ad  (overwrites X_umap only; DPT untouched)
        data/islet_explorer.h5ad  (rebuilt via build_h5ad_for_app.py)
"""
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent
traj_path = ROOT / "data" / "adata_ins_root.h5ad"

print(f"Loading {traj_path}")
adata = sc.read_h5ad(traj_path)
print(f"  Shape: {adata.shape}")
print(f"  Old UMAP range: X [{adata.obsm['X_umap'][:,0].min():.2f}, {adata.obsm['X_umap'][:,0].max():.2f}]"
      f"  Y [{adata.obsm['X_umap'][:,1].min():.2f}, {adata.obsm['X_umap'][:,1].max():.2f}]")

# PCA of raw marker means (31 features -> 15 PCs)
print("\nComputing PCA of raw marker expression...")
sc.pp.pca(adata, n_comps=15)

# Separate neighbor graph for visualization (does NOT touch DPT neighbors)
print("Computing visualization neighbors (PCA, euclidean, k=30)...")
sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca', metric='euclidean',
                key_added='neighbors_viz')

# UMAP with parameters that reveal global disease-stage structure
print("Computing visualization UMAP (min_dist=0.5, spread=2.0)...")
sc.tl.umap(adata, neighbors_key='neighbors_viz', min_dist=0.5, spread=2.0)
print(f"  New UMAP range: X [{adata.obsm['X_umap'][:,0].min():.2f}, {adata.obsm['X_umap'][:,0].max():.2f}]"
      f"  Y [{adata.obsm['X_umap'][:,1].min():.2f}, {adata.obsm['X_umap'][:,1].max():.2f}]")

# Verify DPT pseudotime is untouched
pt = adata.obs['dpt_pseudotime'].values
print(f"\nDPT pseudotime range: [{np.nanmin(pt):.3f}, {np.nanmax(pt):.3f}] (unchanged)")

# Quick separation check: mean UMAP-1 by donor status
for status in ['ND', 'Aab+', 'T1D']:
    mask = adata.obs['donor_status'] == status
    u1 = adata.obsm['X_umap'][mask, 0].mean()
    u2 = adata.obsm['X_umap'][mask, 1].mean()
    print(f"  {status:5s}: UMAP centroid = ({u1:.2f}, {u2:.2f})  n={mask.sum()}")

# Clean up temporary keys before saving
for key in ['neighbors_viz']:
    if key in adata.uns:
        del adata.uns[key]
for key in ['neighbors_viz_distances', 'neighbors_viz_connectivities']:
    if key in adata.obsp:
        del adata.obsp[key]
if 'X_pca' in adata.obsm:
    del adata.obsm['X_pca']
if 'pca' in adata.uns:
    del adata.uns['pca']

print(f"\nSaving {traj_path}")
adata.write_h5ad(traj_path)
print("Done. Now run: python scripts/build_h5ad_for_app.py")
