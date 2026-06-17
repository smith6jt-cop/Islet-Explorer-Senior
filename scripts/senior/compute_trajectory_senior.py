#!/usr/bin/env python3
"""
3D Senior trajectory + islet Leiden — faithful reproduction of the PROVEN Junior
notebook recipe (notebooks/rebuild_trajectory.ipynb), now that aggregation emits
obsm['X_scVI_mean'] (the per-islet mean of the single-cell scVI latent).

Trajectory (rebuild_trajectory.ipynb cells 10-13), NO Harmony:
  neighbors(X_scVI_mean, n=15, cosine) -> PAGA(donor_status) -> UMAP(init='paga',
  min_dist=0.01, spread=2.5) -> root = ND islet argmax INS -> diffmap(10) -> DPT.
  DC1/DC2 = X_diffmap[:,0]/[:,1].

Islet Leiden (Spatial tab) — reaggregate_islets.py recipe: Harmony(theta=2,
batch=imageid) on X_scVI_mean -> separate neighbors graph -> leiden at 4 resolutions.

Validation (the notebook's own gates): spearman(pt,INS) < -0.3, GCG > 0.2,
ND < Aab+ < T1D, no donor dominates the pseudotime quintiles. Run in `scimap` env.
"""
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

SD = Path(__file__).resolve().parents[2] / "data"
F = SD / "islets_core_senior.h5ad"


def eta2(y, groups):
    s = pd.Series(np.asarray(y, dtype=float)); grand = s.mean(); sst = ((s - grand) ** 2).sum()
    ssb = sum(len(v) * (v.mean() - grand) ** 2 for _, v in s.groupby(pd.Series(np.asarray(groups))))
    return float(ssb / sst) if sst > 0 else float("nan")


def main():
    a = ad.read_h5ad(F)
    assert "X_scVI_mean" in a.obsm, (
        "X_scVI_mean missing from islets_core_senior.h5ad — re-run aggregate_islets.py "
        "with the scVI latent attached (train_scvi_senior.py must run first)."
    )
    ds = a.obs["donor_status"].astype(str).to_numpy()
    don = a.obs["imageid"].astype(str).to_numpy()
    ins_idx = list(a.var_names).index("INS")

    # ---- Trajectory: neighbors on X_scVI_mean (batch-corrected representation) ----
    print("[traj] neighbors on X_scVI_mean (n=15, cosine)")
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_scVI_mean", metric="cosine")
    sc.tl.paga(a, groups="donor_status")
    sc.pl.paga(a, plot=False)
    sc.tl.umap(a, init_pos="paga", min_dist=0.01, spread=2.5)

    # root = ND islet with highest INS (argmax — matches rebuild_trajectory.ipynb)
    nd_indices = np.where(ds == "ND")[0]
    ins_vals = np.asarray(a.X[nd_indices, ins_idx]).ravel()
    root_idx = int(nd_indices[np.argmax(ins_vals)])
    a.uns["iroot"] = root_idx
    print(f"[traj] root islet idx={root_idx} donor={don[root_idx]} INS={a.X[root_idx, ins_idx]:.2f}")

    sc.tl.diffmap(a, n_comps=10)
    sc.tl.dpt(a)
    a.obs["DC1"] = a.obsm["X_diffmap"][:, 0]
    a.obs["DC2"] = a.obsm["X_diffmap"][:, 1]

    # ---- Islet Leiden: Harmony(theta=2) on X_scVI_mean, separate graph ----
    leiden_ok = False
    try:
        import harmonypy as hm
        meta = a.obs.copy(); meta["_imageid_str"] = pd.Categorical(don)
        ho = hm.run_harmony(np.asarray(a.obsm["X_scVI_mean"], dtype=np.float64),
                            meta, "_imageid_str", theta=2.0, max_iter_harmony=20)
        a.obsm["X_scVI_harmony"] = np.asarray(ho.Z_corr).T          # harmonypy returns (d, N)
        sc.pp.neighbors(a, use_rep="X_scVI_harmony", n_neighbors=15, metric="cosine",
                        key_added="leiden_nn")
        for res in [0.3, 0.5, 0.8, 1.0]:
            sc.tl.leiden(a, resolution=res, key_added=f"leiden_{res}", neighbors_key="leiden_nn")
        leiden_ok = True
        print("[leiden] Harmony(theta=2) on X_scVI_mean -> leiden_{0.3,0.5,0.8,1.0}")
    except Exception as e:
        print("[leiden skipped]:", type(e).__name__, e)

    a.write_h5ad(F)

    # ---- Validation (rebuild_trajectory.ipynb gates) ----
    pt = a.obs["dpt_pseudotime"].to_numpy()
    ins_expr = np.asarray(a.X[:, ins_idx]).ravel()
    r_ins = spearmanr(pt, ins_expr, nan_policy="omit")[0]
    print("\n" + "=" * 56 + "\nTRAJECTORY VALIDATION")
    print(f"1. spearman(pt, INS) = {r_ins:+.3f}  (expect < -0.3)  -> {'PASS' if r_ins < -0.3 else 'FAIL'}")
    if "GCG" in a.var_names:
        r_gcg = spearmanr(pt, np.asarray(a.X[:, list(a.var_names).index('GCG')]).ravel(), nan_policy="omit")[0]
        print(f"2. spearman(pt, GCG) = {r_gcg:+.3f}  (expect > 0.2)   -> {'PASS' if r_gcg > 0.2 else 'FAIL'}")
    means = {s: float(pt[ds == s].mean()) for s in ["ND", "Aab+", "T1D"] if (ds == s).any()}
    print("3. pseudotime by status:", {k: round(v, 3) for k, v in means.items()})
    ordered = list(means.values())
    print(f"   ND < Aab+ < T1D -> {'PASS' if ordered == sorted(ordered) else 'FAIL'}")
    print(f"4. eta2 disease={eta2(pt, ds):.3f}  eta2 BATCH={eta2(pt, don):.3f}  "
          f"(was disease=0.103 batch=0.988; batch must drop)")
    if leiden_ok:
        print(f"5. leiden_0.5 clusters={a.obs['leiden_0.5'].nunique()} | "
              f"eta2(leiden_0.5 donor-purity proxy) — check Spatial tab")
    print("[done]", F)


if __name__ == "__main__":
    main()
