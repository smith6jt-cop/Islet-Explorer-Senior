#!/usr/bin/env python3
"""
Phase 2 — validate a phenotyping WITHOUT ground-truth labels (run in `scvi-env`).

No expert-labelled cells exist, so accuracy is assessed by biological coherence instead:
  (A) Marker coherence  — each phenotype's mean arcsinh profile should be enriched (z>0 across
      phenotypes) for its rule-defining positive markers (Beta=INS, T cell=CD3e, ...).
  (B) Spatial coherence — homotypic k-NN fraction from X/Y_centroid (cells of a type cluster
      physically; endocrine in islet cores). Computed per donor on a subsample.
  (C) Per-donor composition — phenotype frequencies per donor / Donor Status, with the expected
      directional trend (beta-cell loss + immune gain ND -> Aab+ -> T1D) as a SANITY check.
  (D) Cross-method concordance — if several label columns are present (phenotype_gate /
      phenotype_cluster / phenotype_astir), their agreement matrix; ~60-80% is realistic.

Inputs: a phenotyped h5ad with obs[--label] (default 'phenotype'), obs imageid / Donor Status /
X_centroid / Y_centroid, and X = raw MFI (for marker coherence). Optionally other label columns.

Usage:
  python scripts/senior/validate_phenotypes.py --input <phenotyped.h5ad> --label phenotype \
         --out-prefix data/validation
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

RULES = Path(__file__).resolve().parent / "phenotype_rules_senior.csv"
POS_KW = {"pos", "anypos", "allpos"}
ENDO = {"Beta cell", "Alpha cell", "Delta cell"}


def own_pos_markers(rules_path):
    df = pd.read_csv(rules_path); mk = list(df.columns[2:])
    out = {}
    for _, r in df.iterrows():
        nm = str(r["phenotype_name"]).strip()
        out[nm] = [m for m in mk if not pd.isna(r[m]) and str(r[m]).strip().lower() in POS_KW]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--label", default="phenotype")
    ap.add_argument("--rules", type=Path, default=RULES)
    ap.add_argument("--donor-col", default="imageid")
    ap.add_argument("--status-col", default="Donor Status")
    ap.add_argument("--cofactor", type=float, default=5.0)
    ap.add_argument("--spatial-cells", type=int, default=20000, help="per-donor subsample for k-NN")
    ap.add_argument("--out-prefix", type=Path, required=True)
    a = ap.parse_args(argv)

    adata = ad.read_h5ad(a.input)
    obs = adata.obs
    lab = obs[a.label].astype(str)
    mk = list(adata.var_names)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    arc = np.arcsinh(np.asarray(X, float) / a.cofactor)
    pos_markers = own_pos_markers(a.rules)
    a.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    report = []

    # (A) marker coherence
    phenos = [p for p in lab.value_counts().index if p in pos_markers and pos_markers[p]]
    pmean = pd.DataFrame({p: arc[(lab == p).values].mean(0) for p in lab.unique() if (lab == p).sum() > 0},
                         index=mk).T
    z = (pmean - pmean.mean(0)) / pmean.std(0).replace(0, 1)
    report.append("=== (A) MARKER COHERENCE (z of phenotype mean across phenotypes; want >0 for own pos markers) ===")
    npass = 0
    for p in phenos:
        zs = [z.loc[p, m] for m in pos_markers[p] if m in z.columns]
        score = float(np.mean(zs)) if zs else float("nan")
        ok = score > 0
        npass += ok
        report.append(f"  {p:24s} pos={','.join(pos_markers[p]):28s} mean_z={score:+.2f}  {'PASS' if ok else 'WARN'}")
    report.append(f"  -> {npass}/{len(phenos)} phenotypes coherent")

    # (B) spatial coherence (per-donor subsample k-NN homotypic fraction)
    report.append("\n=== (B) SPATIAL COHERENCE (homotypic 10-NN fraction; want > chance) ===")
    if {"X_centroid", "Y_centroid"}.issubset(obs.columns):
        from sklearn.neighbors import NearestNeighbors
        fracs = []
        for d in pd.unique(obs[a.donor_col].astype(str)):
            m = (obs[a.donor_col].astype(str) == d).values
            idx = np.where(m)[0]
            if idx.size > a.spatial_cells:
                idx = idx[np.linspace(0, idx.size - 1, a.spatial_cells).astype(int)]
            xy = obs.loc[obs.index[idx], ["X_centroid", "Y_centroid"]].to_numpy(float)
            ld = lab.values[idx]
            nn = NearestNeighbors(n_neighbors=11).fit(xy)
            _, ind = nn.kneighbors(xy)
            same = np.mean([(ld[ind[i, 1:]] == ld[i]).mean() for i in range(len(idx))])
            fracs.append(same)
        chance = float((lab.value_counts(normalize=True) ** 2).sum())
        report.append(f"  mean homotypic 10-NN fraction = {np.mean(fracs):.3f} (chance≈{chance:.3f})  "
                      f"{'PASS' if np.mean(fracs) > 2 * chance else 'WARN'}")
    else:
        report.append("  (no centroids in obs — skipped)")

    # (C) per-donor composition + status trend
    report.append("\n=== (C) PER-DONOR COMPOSITION + status trend (sanity) ===")
    comp = pd.crosstab(obs[a.donor_col], lab, normalize="index")
    comp.to_csv(f"{a.out_prefix}_composition_per_donor.csv")
    if a.status_col in obs.columns:
        st = obs.groupby(a.donor_col)[a.status_col].first()
        order = {"ND": 0, "Aab+": 1, "T1D": 2}
        for target in ["Beta cell"]:
            if target in comp.columns:
                m = comp[target].groupby(st).mean()
                report.append(f"  mean {target} fraction by status: " +
                              ", ".join(f"{s}={m.get(s, float('nan')):.3f}" for s in ["ND", "Aab+", "T1D"] if s in m.index))
        # immune gain
        imm = [c for c in comp.columns if c in {"CD8 T cell", "CD4 T cell", "Regulatory T cell",
                                                "Macrophage", "M2 macrophage", "B cell", "Granulocyte"}]
        if imm:
            im = comp[imm].sum(1).groupby(st).mean()
            report.append(f"  mean immune fraction by status: " +
                          ", ".join(f"{s}={im.get(s, float('nan')):.3f}" for s in ["ND", "Aab+", "T1D"] if s in im.index))

    # (D) cross-method concordance
    others = [c for c in ["phenotype_gate", "phenotype_cluster", "phenotype_astir"] if c in obs.columns and c != a.label]
    if others:
        report.append("\n=== (D) CROSS-METHOD CONCORDANCE ===")
        for c in others:
            agree = float((obs[c].astype(str).values == lab.values).mean())
            report.append(f"  {a.label} vs {c}: {agree:.1%} agree")

    txt = "\n".join(report)
    print(txt)
    Path(f"{a.out_prefix}_report.txt").write_text(txt)

    # marker coherence heatmap
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        rule_mk = sorted({m for p in pos_markers for m in pos_markers[p]})
        zr = z.reindex(columns=rule_mk)
        fig, ax = plt.subplots(figsize=(0.32 * len(rule_mk) + 3, 0.32 * len(z) + 2))
        im = ax.imshow(zr.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set_xticks(range(len(rule_mk))); ax.set_xticklabels(rule_mk, rotation=90, fontsize=7)
        ax.set_yticks(range(len(zr.index))); ax.set_yticklabels(zr.index, fontsize=7)
        fig.colorbar(im, ax=ax, fraction=0.02); fig.tight_layout()
        fig.savefig(f"{a.out_prefix}_marker_coherence.png", dpi=150); plt.close(fig)
        print(f"[out] wrote {a.out_prefix}_marker_coherence.png")
    except Exception as e:
        print("[warn] heatmap failed:", e)
    print(f"[done] wrote {a.out_prefix}_report.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
