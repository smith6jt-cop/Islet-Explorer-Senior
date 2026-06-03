#!/usr/bin/env python3
"""
3B(+) Senior trajectory + Leiden on the islet-level dataset (9,130 islets).

No scVI for the new panel, so embed on PCA of scaled log-mean expression
(mirrors reaggregate_islets' intent). DPT rooted at the highest-INS ND islet;
Leiden at 4 resolutions (best-effort). Writes X_umap / X_diffmap / dpt_pseudotime
/ DC1,DC2 / leiden_* into islets_core_senior.h5ad. Run in the `scimap` env.
"""
from pathlib import Path
import numpy as np
import scanpy as sc
import anndata as ad
from scipy.stats import spearmanr

SD = Path.home() / "IO60panc2nd" / "senior_data"
F = SD / "islets_core_senior.h5ad"


def main():
    a = ad.read_h5ad(F)
    w = a.copy()
    sc.pp.log1p(w)
    sc.pp.scale(w, max_value=10)
    sc.pp.pca(w, n_comps=min(30, w.n_vars - 1))
    sc.pp.neighbors(w, n_neighbors=15, use_rep="X_pca", metric="cosine")
    sc.tl.umap(w, min_dist=0.1)

    ins = list(a.var_names).index("INS")
    nd = (a.obs["donor_status"].astype(str) == "ND").to_numpy()
    Xins = np.asarray(a.X[:, ins]).ravel()
    w.uns["iroot"] = int(np.where(nd)[0][np.argmax(Xins[nd])])
    sc.tl.diffmap(w, n_comps=10)
    sc.tl.dpt(w)

    leiden_ok = False
    try:
        for res in [0.3, 0.5, 0.8, 1.0]:
            sc.tl.leiden(w, resolution=res, key_added=f"leiden_{res}")
        leiden_ok = True
    except Exception as e:
        print("[leiden skipped]:", type(e).__name__, e)

    a.obsm["X_umap"] = w.obsm["X_umap"]
    a.obsm["X_diffmap"] = w.obsm["X_diffmap"]
    a.obs["dpt_pseudotime"] = w.obs["dpt_pseudotime"].to_numpy()
    a.obs["DC1"] = w.obsm["X_diffmap"][:, 1]
    a.obs["DC2"] = w.obsm["X_diffmap"][:, 2]
    if leiden_ok:
        for res in [0.3, 0.5, 0.8, 1.0]:
            a.obs[f"leiden_{res}"] = w.obs[f"leiden_{res}"].to_numpy()
    a.write_h5ad(F)

    print("[traj] pseudotime by donor status:")
    print(a.obs.groupby("donor_status", observed=True)["dpt_pseudotime"].mean().round(3).to_string())
    if "posfrac_INS" in a.obs.columns:
        r, _ = spearmanr(a.obs["dpt_pseudotime"], a.obs["posfrac_INS"])
        print(f"[traj] spearman(pseudotime, INS+ frac) = {r:.3f}  (expect negative)")
    if leiden_ok:
        print("[traj] leiden_0.5 clusters:", a.obs["leiden_0.5"].nunique())


if __name__ == "__main__":
    main()
