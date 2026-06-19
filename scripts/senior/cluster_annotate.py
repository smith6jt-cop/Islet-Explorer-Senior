#!/usr/bin/env python3
"""
Phase 1 (cross-check A) — Leiden clustering + automated hierarchical cluster annotation
(run in the `scvi-env` conda env: scanpy 1.11, leidenalg, harmonypy).

The MULTIVARIATE fix for the univariate-gating failures (e.g. SOX2 over-called positive in
~90% of cells → 72% of cells mislabelled Ductal in Phase 0). Here cells are clustered on the
JOINT marker profile (batch-corrected across donors), and each cluster is annotated by walking
the SAME rules hierarchy as scimap — but scoring cluster-mean z-scores and descending into the
best-matching child at each level. A SOX2-high-but-otherwise-acinar cell clusters with acinar,
so that cluster's MEAN SOX2 stays low and the cluster is correctly called Acinar.

Pipeline:
  raw MFI -> arcsinh(x/cofactor) -> per-marker z-score -> PCA
          -> [Harmony(imageid) if >1 donor] -> kNN -> Leiden(resolution)
          -> cluster-mean z-scores -> hierarchical rules annotation -> obs['phenotype_cluster']

Outputs (per --out-prefix): <prefix>.h5ad? no — writes obs columns to <prefix>_cells.parquet
(cell_id, leiden, phenotype_cluster, annot_margin), a cluster->phenotype table CSV, a
cluster x marker mean-z heatmap PNG, and obsm embeddings as .npy for reuse/propagation.

Smoke (single donor, no batch correction needed):
  python scripts/senior/cluster_annotate.py --donors 6374 --max-cells-per-donor 200000 \
         --resolution 1.5 --out-prefix /tmp/clust_6374
Subsample-cluster for full-cohort (stratified, then propagate separately):
  python scripts/senior/cluster_annotate.py --max-cells-per-donor 150000 --harmony \
         --resolution 1.5 --out-prefix data/cluster_subsample
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

BASE = Path(__file__).resolve().parents[2] / "data"
INP = BASE / "singlecell_protein.h5ad"
RULES = Path(__file__).resolve().parent / "phenotype_rules_senior.csv"

POS_KW = {"pos", "anypos", "allpos"}
NEG_KW = {"neg", "allneg", "anyneg"}


def harmony_embed(pca, donor_labels, donor_col, theta=2.0):
    """Batch-correct a PCA embedding with harmonypy, robust to the GPU/PyTorch build's Z_corr
    orientation (scanpy's harmony_integrate wrapper mis-transposes it -> shape error)."""
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


# ------------------------------------------------------------------ rules tree
def parse_rules(rules_path):
    """Return (children, sig): children[parent] -> [phenotype...]; sig[pheno] -> (pos,neg)
    of that phenotype's OWN markers (anypos/pos/allpos -> pos; neg/allneg -> neg). Hierarchy
    via the parent_class column ('all' is the root)."""
    df = pd.read_csv(rules_path)
    markers = list(df.columns[2:])
    children = {}
    sig = {}
    for _, r in df.iterrows():
        parent = str(r["parent_class"]).strip()
        name = str(r["phenotype_name"]).strip()
        children.setdefault(parent, []).append(name)
        pos, neg = [], []
        for m in markers:
            v = r[m]
            if pd.isna(v):
                continue
            v = str(v).strip().lower()
            if v in POS_KW:
                pos.append(m)
            elif v in NEG_KW:
                neg.append(m)
        sig[name] = (pos, neg)
    return children, sig


# Markers the reliability audit flags as background/non-specific (high positive-fraction, low
# cluster specificity) or housekeeping/proliferation -- excluded from annotation SIGNATURES so a
# bad marker (e.g. SOX2 positive in ~90% of cells) cannot drive the call.
AUDIT_BACKGROUND = {"SOX2", "EpCAM", "E_cadherin", "DAPI", "Beta_actin", "Bcl_2", "PCNA", "Ki67"}


def annotate_cluster(zrow, children, sig, lineage_floor=0.3, lineage_margin=0.05,
                     drop=AUDIT_BACKGROUND, root="all"):
    """Annotate a cluster's mean z-profile against the rules hierarchy. Gives a best-guess for
    EVERY cluster (a starting point for heatmap QC) rather than rejecting most to 'Unassigned'.

    1. Dominant lineage = argmax of mean(z[lineage own pos]); pos-only, since negative penalties
       backfire on high-background clusters. Audit-flagged background markers are dropped.
    2. If the winning lineage's strength < lineage_floor (nothing actually expressed) ->
       'Unassigned (low signal)' (else low-signal clusters default to whatever wins by a hair).
    3. Otherwise ALWAYS descend to a leaf by per-level argmax of own-marker score
       mean(z[own_pos]) - mean(z[own_neg]) -- OWN markers, so Beta/Alpha/Delta are distinguished
       by INS/GCG/SST (inheriting the lineage union would make them identical).

    Returns (label, path, lin_margin, n_strong, sub_alts) where sub_alts is the ranked top-3
    (subtype, score) at the final branch -- the alternatives a human weighs on the heatmap."""
    def opos(n):
        return [m for m in sig.get(n, ([], []))[0] if m not in drop]

    def oneg(n):
        return [m for m in sig.get(n, ([], []))[1] if m not in drop]

    def own_score(n):
        p, q = opos(n), oneg(n)
        sp = np.mean([zrow[m] for m in p]) if p else 0.0
        sn = np.mean([zrow[m] for m in q]) if q else 0.0
        return float(sp - sn)

    def lin_strength(L):
        p = opos(L)
        return float(np.mean([zrow[m] for m in p])) if p else -np.inf

    lineages = children.get(root, [])
    ranked_L = sorted(((L, lin_strength(L)) for L in lineages), key=lambda kv: kv[1], reverse=True)
    n_strong = sum(1 for L in lineages if max((zrow[m] for m in opos(L)), default=-9) > 1.0)
    topL, topS = ranked_L[0]
    # margin only vs lineages with >=2 own-pos markers: a degenerate single-marker
    # lineage (e.g. Neural=B3TUBB) must not spuriously win the near-tie test and
    # dump a real cluster to "ambiguous".
    others = [s for L, s in ranked_L if L != topL and len(opos(L)) >= 2]
    lin_margin = float(topS - (max(others) if others else topS - 1))
    if topS < lineage_floor:
        return "Unassigned (low signal)", ["Unassigned"], lin_margin, n_strong, []
    if lin_margin < lineage_margin:    # near-tie across lineages -> genuinely ambiguous/non-specific
        return "Unassigned (ambiguous)", ["Unassigned"], lin_margin, n_strong, []

    node = topL
    path = [node]
    sub_alts = []
    while node in children:
        ranked = sorted(((c, own_score(c)) for c in children[node]), key=lambda kv: kv[1], reverse=True)
        sub_alts = [(c, round(s, 2)) for c, s in ranked[:3]]
        path.append(ranked[0][0]); node = ranked[0][0]
    return path[-1], path, lin_margin, n_strong, sub_alts


# ------------------------------------------------------------------ IO + normalize
def load_subset(path, donors, max_cells_per_donor, donor_col):
    print(f"[load] reading {path} into memory ...")
    a = ad.read_h5ad(path)
    print(f"[load] full {a.shape}")
    if donors is None and max_cells_per_donor is None:
        return a
    dvals = a.obs[donor_col].astype(str).to_numpy()
    sel = np.ones(dvals.size, dtype=bool) if donors is None else np.isin(dvals, [str(x) for x in donors])
    idx = np.where(sel)[0]
    if max_cells_per_donor is not None:
        keep = []
        for d in pd.unique(dvals[idx]):
            di = idx[dvals[idx] == d]
            if di.size > max_cells_per_donor:
                di = di[np.linspace(0, di.size - 1, max_cells_per_donor).astype(int)]
            keep.append(di)
        idx = np.sort(np.concatenate(keep))
    a = a[idx].copy()
    print(f"[load] subset {a.shape} | donors {pd.unique(a.obs[donor_col]).tolist()}")
    return a


def normalize_embed(adata, donor_col, cofactor, n_pcs, use_harmony, resolution, seed=0,
                    exclude_markers=None, winsor_pct=99.5):
    X = np.asarray(adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray(), dtype=np.float64)
    # Tail-pulling: this panel's raw MFI is heavy-tailed, so a few very bright cells
    # inflate the per-marker mean/SD and wash the dim majority toward 0 (binary lineage
    # markers become uninterpretable). Winsorize each marker's bright tail, then arcsinh
    # with a cofactor matched to the MFI scale (~150, not 5 — medians are 10s-100s).
    if winsor_pct is not None and winsor_pct < 100:
        hi = np.percentile(X, winsor_pct, axis=0)
        X = np.minimum(X, hi)
    X = np.arcsinh(X / cofactor)
    adata.layers["arcsinh"] = X.astype(np.float32)       # keep ALL markers for annotation scoring
    # per-marker z-score for the embedding; high-background markers (per the reliability audit) may be
    # dropped from the clustering FEATURES while staying in the arcsinh layer for rule scoring.
    feat = [i for i, m in enumerate(adata.var_names) if not (exclude_markers and m in exclude_markers)]
    if exclude_markers:
        print(f"[features] excluding {sorted(set(exclude_markers) & set(adata.var_names))} from clustering")
    Xf = X[:, feat]
    mu = Xf.mean(0); sd = Xf.std(0); sd[sd == 0] = 1.0
    adata.obsm["X_feat"] = ((Xf - mu) / sd).astype(np.float32)
    from sklearn.decomposition import PCA
    n_pcs = min(n_pcs, len(feat) - 1)
    adata.obsm["X_pca"] = PCA(n_components=n_pcs, random_state=seed).fit_transform(adata.obsm["X_feat"])
    rep = "X_pca"
    n_donors = adata.obs[donor_col].nunique()
    if use_harmony and n_donors > 1:
        print(f"[harmony] integrating {n_donors} donors on {donor_col} ...")
        adata.obsm["X_pca_harmony"] = harmony_embed(adata.obsm["X_pca"], adata.obs[donor_col].values, donor_col)
        rep = "X_pca_harmony"
    print(f"[neighbors+leiden] rep={rep} res={resolution}")
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=15, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=seed, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    return adata, rep


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=INP)
    ap.add_argument("--rules", type=Path, default=RULES)
    ap.add_argument("--donor-col", default="imageid")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--max-cells-per-donor", type=int, default=None)
    ap.add_argument("--cofactor", type=float, default=150.0,
                    help="arcsinh cofactor matched to MFI scale (raw medians are 10s-100s, so ~150, not 5)")
    ap.add_argument("--winsor-pct", type=float, default=99.5,
                    help="per-marker upper-tail winsorize percentile (pulls heavy tails in; <100 to enable)")
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--harmony", action="store_true")
    ap.add_argument("--resolution", type=float, default=2.0)
    ap.add_argument("--exclude-markers", nargs="*", default=None,
                    help="markers to drop from clustering FEATURES (e.g. audit high_background), still scored for annotation")
    ap.add_argument("--out-prefix", type=Path, required=True)
    a = ap.parse_args(argv)

    adata = load_subset(a.input, a.donors, a.max_cells_per_donor, a.donor_col)
    children, sig = parse_rules(a.rules)
    rule_markers = sorted({m for pn in sig for m in (sig[pn][0] + sig[pn][1])})
    missing = [m for m in rule_markers if m not in set(adata.var_names)]
    if missing:
        raise SystemExit(f"[err] rule markers missing: {missing}")

    adata, rep = normalize_embed(adata, a.donor_col, a.cofactor, a.n_pcs, a.harmony, a.resolution,
                                 exclude_markers=a.exclude_markers, winsor_pct=a.winsor_pct)

    # cluster-mean arcsinh, z-scored across clusters
    arc = adata.layers["arcsinh"]
    mk = list(adata.var_names)
    leiden = adata.obs["leiden"].astype(str).values
    clusters = sorted(pd.unique(leiden), key=lambda s: int(s))
    cmean = pd.DataFrame(index=clusters, columns=mk, dtype=float)
    for cl in clusters:
        cmean.loc[cl] = arc[leiden == cl].mean(0)
    z = (cmean - cmean.mean(0)) / cmean.std(0).replace(0, 1.0)

    # annotate each cluster: best-guess + subtype alternatives + flag (for heatmap QC)
    rows = []
    for cl in clusters:
        leaf, path, margin, nstrong, alts = annotate_cluster(z.loc[cl], children, sig)
        flag = bool(str(leaf).startswith("Unassigned") or nstrong >= 4 or margin < 0.15)
        rows.append({"leiden": cl, "phenotype_cluster": leaf, "lineage": path[0] if path else "",
                     "subtype_alts": "; ".join(f"{n}({v})" for n, v in alts),
                     "lin_margin": round(margin, 3), "n_strong_lineages": nstrong,
                     "review_flag": flag, "n_cells": int((leiden == cl).sum())})
    ann = pd.DataFrame(rows).set_index("leiden")
    print("\n[annotation] cluster -> phenotype (review_flag clusters need heatmap QC):")
    print(ann.sort_values("n_cells", ascending=False).to_string())

    # map to cells
    cell_pheno = pd.Series(leiden, index=adata.obs_names).map(ann["phenotype_cluster"])
    cell_margin = pd.Series(leiden, index=adata.obs_names).map(ann["lin_margin"])
    print("\n[annotation] phenotype_cluster cell counts:")
    print(cell_pheno.value_counts().to_string())

    # outputs
    a.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame({
        "cell_id": adata.obs_names.astype(str),
        a.donor_col: adata.obs[a.donor_col].astype(str).values,
        "leiden": leiden,
        "phenotype_cluster": cell_pheno.values,
        "annot_margin": cell_margin.values,
    })
    out.to_parquet(f"{a.out_prefix}_cells.parquet")
    ann.to_csv(f"{a.out_prefix}_cluster_annotation.csv")
    cmean.to_csv(f"{a.out_prefix}_cluster_mean_arcsinh.csv")
    np.save(f"{a.out_prefix}_{rep}.npy", adata.obsm[rep])

    # heatmap of cluster z-scores over rule markers
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        zr = z[rule_markers]
        fig, ax = plt.subplots(figsize=(0.32 * len(rule_markers) + 3, 0.32 * len(clusters) + 2))
        im = ax.imshow(zr.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set_xticks(range(len(rule_markers))); ax.set_xticklabels(rule_markers, rotation=90, fontsize=7)
        ax.set_yticks(range(len(clusters)))
        ax.set_yticklabels([f"{cl}:{ann.loc[cl,'phenotype_cluster']}" for cl in clusters], fontsize=7)
        fig.colorbar(im, ax=ax, fraction=0.02); fig.tight_layout()
        fig.savefig(f"{a.out_prefix}_heatmap.png", dpi=150); plt.close(fig)
        print(f"[out] wrote {a.out_prefix}_heatmap.png")
    except Exception as e:
        print("[warn] heatmap failed:", e)
    print(f"[done] wrote {a.out_prefix}_cells.parquet / _cluster_annotation.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
