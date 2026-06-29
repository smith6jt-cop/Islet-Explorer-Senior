#!/usr/bin/env python3
"""
Marker-reliability audit for the Senior 59-plex panel (run in `scvi-env`).

Empirically, several markers sabotage gating AND Astir on this panel (SOX2/EpCAM/Ker8_18/Vimentin/
SMA over-gated ~90% positive everywhere; Pan_Cytokeratin under-gated ~1-5%). This script scores every
marker so the over/under-gated and batch-unstable channels can be dropped or down-weighted — the single
fix that improves gating, Astir, and clustering feature selection at once.

Two independent signals per marker (so we can tell 'bad antibody' from 'legitimately common'):
  GATE-BASED  : per-(marker,donor) 2-comp GMM gate -> global gate-positive fraction, GMM separation,
                and the cross-donor CV of the positive fraction (batch stability).
  GATE-FREE   : cluster the cells (arcsinh -> Harmony -> Leiden); a real marker is concentrated in a
                FEW clusters (high CV of cluster means, low entropy); a high-background marker is
                uniformly elevated (low CV). This needs no gate, so it breaks the tie that gating can't:
                EpCAM (common but cluster-structured) vs SOX2 (uniformly high = artifact).

Verdict heuristic per marker:
  high_background : pos_frac>0.5 AND cluster_cv low  (uniformly positive, not cluster-specific) -> DROP
  under_gated     : pos_frac<0.05 AND a cluster mean is clearly high  (signal exists, gate too high)
  batch_unstable  : cross-donor pos_frac CV high
  informative     : otherwise (cluster-specific, well-separated)         -> KEEP

Outputs: <out_prefix>_marker_audit.csv (ranked table) + <out_prefix>_distributions.pdf (per-marker
per-donor log1p densities with gate lines) + <out_prefix>_cluster_heatmap.png.

Usage:
  python scripts/senior/marker_reliability_audit.py --per-donor-cells 50000 \
         --out-prefix data/marker_audit
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.mixture import GaussianMixture

BASE = Path(__file__).resolve().parents[2] / "data"
INP = BASE / "singlecell_protein.h5ad"


def harmony_embed(pca, donor_labels, donor_col, theta=2.0):
    """Batch-correct a PCA embedding with harmonypy, robust to the GPU/PyTorch build's
    Z_corr orientation (newer harmonypy returns cells×dims; scanpy's wrapper wrongly assumes
    dims×cells and transposes -> shape error). We orient by matching the cell count."""
    import harmonypy
    meta = pd.DataFrame({donor_col: np.asarray(donor_labels).astype(str)})
    ho = harmonypy.run_harmony(np.asarray(pca, dtype=np.float64), meta, [donor_col], theta=theta)
    Z = ho.Z_corr
    if hasattr(Z, "cpu"):
        Z = Z.cpu().numpy()
    Z = np.asarray(Z)
    if Z.shape[0] != pca.shape[0]:
        Z = Z.T
    assert Z.shape[0] == pca.shape[0], f"harmony shape {Z.shape} vs pca {pca.shape}"
    return Z.astype(np.float32)


def gmm_gate_sep(x, min_sep=2.0):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if x.size < 50:
        return float("nan"), 0.0
    lx = np.log1p(np.clip(x, 0, None))
    sub = lx[np.linspace(0, lx.size - 1, min(40000, lx.size)).astype(int)].reshape(-1, 1)
    try:
        g = GaussianMixture(2, n_init=1, random_state=0, reg_covar=1e-3).fit(sub)
    except Exception:
        return float(np.quantile(lx, 0.99)), 0.0
    mu = g.means_.ravel(); sd = np.sqrt(g.covariances_.ravel())
    hi = int(np.argmax(mu)); lo = 1 - hi
    pooled = max(0.5 * (sd[hi] + sd[lo]), 1e-6)
    sep = float((mu[hi] - mu[lo]) / pooled)
    if (mu[hi] - mu[lo]) < min_sep * pooled:
        return float(np.quantile(lx, 0.99)), sep
    grid = np.linspace(mu[lo], mu[hi], 256).reshape(-1, 1)
    p = g.predict_proba(grid)[:, hi]
    cr = np.where(p >= 0.5)[0]
    return float(grid[cr[0], 0] if cr.size else mu[hi]), sep


def load_strat(path, per_donor, donor_col):
    print(f"[load] reading {path} ...")
    a = ad.read_h5ad(path)
    dvals = a.obs[donor_col].astype(str).to_numpy()
    keep = []
    for d in pd.unique(dvals):
        di = np.where(dvals == d)[0]
        if per_donor and di.size > per_donor:
            di = di[np.linspace(0, di.size - 1, per_donor).astype(int)]
        keep.append(di)
    idx = np.sort(np.concatenate(keep))
    a = a[idx].copy()
    print(f"[load] stratified subsample {a.shape} | {a.obs[donor_col].nunique()} donors")
    return a


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=INP)
    ap.add_argument("--donor-col", default="imageid")
    ap.add_argument("--per-donor-cells", type=int, default=50000)
    ap.add_argument("--cofactor", type=float, default=5.0)
    ap.add_argument("--resolution", type=float, default=2.0)
    ap.add_argument("--out-prefix", type=Path, required=True)
    a = ap.parse_args(argv)

    adata = load_strat(a.input, a.per_donor_cells, a.donor_col)
    markers = list(adata.var_names)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    X = np.asarray(X, np.float64)
    didx = adata.obs[a.donor_col].astype(str).to_numpy()
    donors = list(pd.unique(didx))

    # --- gate-based stats per (marker, donor)
    print("[gates] per-(marker,donor) GMM ...")
    pos_frac = pd.DataFrame(index=markers, columns=donors, dtype=float)
    seps = pd.DataFrame(index=markers, columns=donors, dtype=float)
    for d in donors:
        rows = np.where(didx == d)[0]
        lxd = np.log1p(np.clip(X[rows], 0, None))
        for j, m in enumerate(markers):
            g, sep = gmm_gate_sep(X[rows, j])
            pos_frac.loc[m, d] = float((lxd[:, j] >= g).mean())
            seps.loc[m, d] = sep

    # --- gate-free cluster-informativeness
    print("[cluster] arcsinh -> Harmony -> Leiden ...")
    arc = np.arcsinh(X / a.cofactor)
    adata.layers["arcsinh"] = arc.astype(np.float32)
    mu = arc.mean(0); sd = arc.std(0); sd[sd == 0] = 1
    adata.X = ((arc - mu) / sd).astype(np.float32)
    sc.pp.pca(adata, n_comps=30, random_state=0)
    rep = "X_pca"
    if adata.obs[a.donor_col].nunique() > 1:
        adata.obsm["X_pca_harmony"] = harmony_embed(adata.obsm["X_pca"], adata.obs[a.donor_col].values, a.donor_col)
        rep = "X_pca_harmony"
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata, resolution=a.resolution, random_state=0, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False)
    leiden = adata.obs["leiden"].astype(str).to_numpy()
    clusters = sorted(pd.unique(leiden), key=lambda s: int(s))
    cmean = np.vstack([arc[leiden == cl].mean(0) for cl in clusters])      # (C, M)
    cmean_df = pd.DataFrame(cmean, index=clusters, columns=markers)
    # informativeness: CV of cluster means + how concentrated (top cluster share of total signal)
    cm_mean = cmean.mean(0); cm_std = cmean.std(0)
    cluster_cv = cm_std / np.where(cm_mean == 0, 1, cm_mean)
    z = (cmean - cm_mean) / np.where(cm_std == 0, 1, cm_std)
    n_clusters_high = (z > 1.0).sum(0)
    top_cluster_z = z.max(0)

    # --- assemble + verdict
    pf_global = pos_frac.mean(1).values
    pf_cv = (pos_frac.std(1) / pos_frac.mean(1).replace(0, np.nan)).fillna(0).values
    sep_med = seps.median(1).values
    df = pd.DataFrame({
        "marker": markers,
        "pos_frac_global": pf_global.round(3),
        "pos_frac_cv_donor": pf_cv.round(3),
        "gmm_sep_median": np.round(sep_med, 2),
        "cluster_cv": np.round(cluster_cv, 3),
        "n_clusters_high": n_clusters_high,
        "top_cluster_z": np.round(top_cluster_z, 2),
    }).set_index("marker")

    def verdict(r):
        if r.pos_frac_global > 0.5 and r.cluster_cv < 0.12:
            return "high_background_DROP"
        if r.pos_frac_global < 0.05 and r.top_cluster_z > 2.0:
            return "under_gated"
        if r.pos_frac_cv_donor > 1.0:
            return "batch_unstable"
        if r.cluster_cv >= 0.12 and r.n_clusters_high >= 1:
            return "informative_KEEP"
        return "weak_review"
    df["verdict"] = df.apply(verdict, axis=1)
    df = df.sort_values(["verdict", "cluster_cv"], ascending=[True, False])

    pd.set_option("display.width", 200)
    print("\n=== MARKER RELIABILITY AUDIT ===")
    print(df.to_string())
    print("\n=== verdict counts ===")
    print(df["verdict"].value_counts().to_string())

    a.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(f"{a.out_prefix}_marker_audit.csv")
    cmean_df.to_csv(f"{a.out_prefix}_cluster_means.csv")

    # per-marker per-donor density PDF
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        ncol = 5; nrow = int(np.ceil(len(donors) / ncol))
        with PdfPages(f"{a.out_prefix}_distributions.pdf") as pdf:
            for j, m in enumerate(markers):
                fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 2 * nrow), squeeze=False)
                for k, d in enumerate(donors):
                    ax = axes[k // ncol][k % ncol]
                    rows = np.where(didx == d)[0]
                    ax.hist(np.log1p(np.clip(X[rows, j], 0, None)), bins=70, color="#4477AA", alpha=.8)
                    ax.set_title(f"{d} pf={pos_frac.loc[m,d]:.2f}", fontsize=7); ax.set_yticks([])
                for k in range(len(donors), nrow * ncol):
                    axes[k // ncol][k % ncol].axis("off")
                fig.suptitle(f"{m} | {df.loc[m,'verdict']} | clusterCV={df.loc[m,'cluster_cv']:.2f} "
                             f"posfrac={df.loc[m,'pos_frac_global']:.2f}", fontsize=10)
                fig.tight_layout(rect=[0, 0, 1, 0.96]); pdf.savefig(fig); plt.close(fig)
        print(f"[out] wrote {a.out_prefix}_distributions.pdf")
    except Exception as e:
        print("[warn] pdf failed:", e)
    print(f"[done] wrote {a.out_prefix}_marker_audit.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
