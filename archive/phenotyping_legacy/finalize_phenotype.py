#!/usr/bin/env python3
"""
Phase 1 (finalize) — propagate cluster labels to all 23 M cells + multi-method consensus
(run in `scvi-env`).

Clustering+annotation is done on a stratified per-donor subsample (cluster_annotate.py). This script
scales those labels to the full cohort and fuses the independent typers into a confidence-scored
per-cell phenotype:

  1. Train a classifier on the SUBSAMPLE (arcsinh marker features -> curated cluster phenotype). The
     subsample spans all 15 donors, so the classifier learns a donor-robust marker->phenotype map.
     'Unassigned*' clusters are excluded from training (don't learn from noise).
  2. Predict the phenotype for all 23 M cells in chunks (fast; LogReg/HGB), with the max class
     probability as a per-cell confidence.
  3. Fuse with the other typers by cell_id when available (Phase-0 gate label, Astir label): where the
     propagated label agrees with Astir -> high confidence; disagreements / low probability -> flagged
     for the cluster/QC review. Output obs['phenotype'] (+ _confidence, _flag, and each method's column).

Inputs:
  --full-h5ad            the 23 M-cell raw h5ad (singlecell_protein.h5ad)
  --subsample-parquet    cluster_annotate <prefix>_cells.parquet (cell_id, imageid, leiden, phenotype_cluster)
  --cluster-labels       optional CSV (leiden,phenotype) — your curated relabelling of the heatmap
                         (defaults to the auto phenotype_cluster)
  --gate-h5ad            optional Phase-0 phenotyped h5ad (obs['phenotype']) joined as phenotype_gate
  --astir-parquet        optional run_astir <prefix>_cells.parquet (phenotype_astir, astir_prob)
  --exclude-markers      markers excluded from the clustering features (match cluster_annotate)

Usage:
  python scripts/senior/finalize_phenotype.py \
     --full-h5ad data/singlecell_protein.h5ad \
     --subsample-parquet data/cluster_subsample_cells.parquet \
     --astir-parquet data/astir_full_cells.parquet \
     --out data/singlecell_protein_phenotyped.h5ad
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--full-h5ad", type=Path, required=True)
    ap.add_argument("--subsample-parquet", type=Path, required=True)
    ap.add_argument("--cluster-labels", type=Path, default=None)
    ap.add_argument("--gate-h5ad", type=Path, default=None)
    ap.add_argument("--astir-parquet", type=Path, default=None)
    ap.add_argument("--cofactor", type=float, default=5.0)
    ap.add_argument("--exclude-markers", nargs="*", default=None)
    ap.add_argument("--classifier", choices=["logreg", "hgb"], default="hgb")
    ap.add_argument("--no-endocrine-split", action="store_true",
                    help="disable splitting Endocrine/Beta cells into Beta/Alpha/Delta by dominant hormone")
    ap.add_argument("--min-prob", type=float, default=0.5, help="below this max-proba -> flagged")
    ap.add_argument("--chunk", type=int, default=2_000_000)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    sub = pd.read_parquet(a.subsample_parquet)
    if a.cluster_labels:
        cl = pd.read_csv(a.cluster_labels)
        cmap = dict(zip(cl["leiden"].astype(str), cl["phenotype"].astype(str)))
        sub["phenotype_cluster"] = sub["leiden"].astype(str).map(cmap).fillna(sub["phenotype_cluster"])
    sub["cell_id"] = sub["cell_id"].astype(str)
    train = sub[~sub["phenotype_cluster"].astype(str).str.startswith("Unassigned")].copy()
    print(f"[train] {len(train):,} labelled subsample cells | {train['phenotype_cluster'].nunique()} classes")

    print(f"[load] {a.full_h5ad}")
    adata = ad.read_h5ad(a.full_h5ad)
    markers = list(adata.var_names)
    feat = [i for i, m in enumerate(markers) if not (a.exclude_markers and m in a.exclude_markers)]
    feat_names = [markers[i] for i in feat]
    Xfull = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    Xfull = np.arcsinh(np.asarray(Xfull, np.float64)[:, feat] / a.cofactor)
    mu = Xfull.mean(0); sd = Xfull.std(0); sd[sd == 0] = 1.0
    Xfull = (Xfull - mu) / sd                                   # z on full-cohort stats
    ids = adata.obs_names.astype(str).to_numpy()
    id_pos = {c: i for i, c in enumerate(ids)}

    tr_rows = np.array([id_pos[c] for c in train["cell_id"] if c in id_pos])
    ytr = train.set_index("cell_id").loc[ids[tr_rows], "phenotype_cluster"].to_numpy()
    Xtr = Xfull[tr_rows]
    print(f"[train] matched {len(tr_rows):,} cells to full h5ad")

    if a.classifier == "logreg":
        from sklearn.linear_model import LogisticRegression
        clf = LogisticRegression(max_iter=400, n_jobs=-1, class_weight="balanced")
    else:
        from sklearn.ensemble import HistGradientBoostingClassifier
        clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, max_depth=None,
                                             l2_regularization=1.0, early_stopping=True)
    clf.fit(Xtr, ytr)
    classes = clf.classes_

    print(f"[predict] {adata.n_obs:,} cells in chunks of {a.chunk:,}")
    pred = np.empty(adata.n_obs, dtype=object)
    conf = np.empty(adata.n_obs, dtype=np.float32)
    for s in range(0, adata.n_obs, a.chunk):
        e = min(s + a.chunk, adata.n_obs)
        P = clf.predict_proba(Xfull[s:e])
        j = P.argmax(1)
        pred[s:e] = classes[j]
        conf[s:e] = P[np.arange(e - s), j]
    adata.obs["phenotype_cluster"] = pd.Categorical(pred)
    adata.obs["phenotype_cluster_prob"] = conf

    # fuse other methods by cell_id
    methods = ["phenotype_cluster"]
    if a.gate_h5ad and Path(a.gate_h5ad).exists():
        g = ad.read_h5ad(a.gate_h5ad, backed="r")
        gp = pd.Series(g.obs["phenotype"].astype(str).values, index=g.obs_names.astype(str))
        adata.obs["phenotype_gate"] = gp.reindex(ids).values
        methods.append("phenotype_gate")
    if a.astir_parquet and Path(a.astir_parquet).exists():
        asr = pd.read_parquet(a.astir_parquet)
        asr["cell_id"] = asr["cell_id"].astype(str)
        am = asr.set_index("cell_id")["phenotype_astir"].reindex(ids).values
        adata.obs["phenotype_astir"] = am
        methods.append("phenotype_astir")

    # consensus: cluster-propagation is primary (the robust grouper); confidence boosted by agreement
    final = adata.obs["phenotype_cluster"].astype(str).to_numpy()
    agree = np.zeros(adata.n_obs, dtype=int)
    for c in ["phenotype_astir", "phenotype_gate"]:
        if c in adata.obs:
            agree += (adata.obs[c].astype(str).to_numpy() == final).astype(int)
    n_other = sum(c in adata.obs for c in ["phenotype_astir", "phenotype_gate"])
    confidence = conf.copy()
    if n_other:
        confidence = 0.6 * conf + 0.4 * (agree / n_other)
    flag = (conf < a.min_prob) | (agree < n_other) if n_other else (conf < a.min_prob)

    # Endocrine subtype split: clusters lump beta/alpha/delta together, but each islet endocrine
    # cell IS defined by its dominant hormone — INS=Beta, GCG=Alpha, SST=Delta. These three markers
    # are the most reliable in the panel (audit top_cluster_z 3.4/2.9/2.5), so split per-cell by the
    # argmax of the standardized hormone (Xfull is already per-marker z-scored arcsinh).
    if not a.no_endocrine_split:
        endo_labels = {"Endocrine", "Beta cell", "Alpha cell", "Delta cell"}
        emask = pd.Series(final).isin(endo_labels).to_numpy()
        if emask.any():
            feat_pos = {gi: j for j, gi in enumerate(feat)}
            hidx = [feat_pos[markers.index(h)] for h in ("INS", "GCG", "SST") if markers.index(h) in feat_pos]
            if len(hidx) == 3:
                sub = np.array(["Beta cell", "Alpha cell", "Delta cell"])[Xfull[:, hidx].argmax(1)]
                final[emask] = sub[emask]
                print(f"[endocrine] split {int(emask.sum()):,} endocrine cells -> "
                      f"{pd.Series(final[emask]).value_counts().to_dict()}")
    adata.obs["phenotype"] = pd.Categorical(final)
    adata.obs["phenotype_confidence"] = confidence.astype(np.float32)
    adata.obs["phenotype_flag"] = flag

    print("\n[final] phenotype distribution:")
    print(adata.obs["phenotype"].value_counts().to_string())
    print(f"\n[final] mean confidence={float(confidence.mean()):.3f} | flagged={float(flag.mean()):.1%} "
          f"| methods fused: {methods}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(a.out)
    print(f"[done] wrote {a.out} {adata.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
