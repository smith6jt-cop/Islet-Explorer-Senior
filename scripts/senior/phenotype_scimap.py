#!/usr/bin/env python3
"""
scimap phenotyping for Islet Explorer Senior (run in the `scimap` conda env).

Mirrors Junior's CODEX_Panc_scimap_Analysis.ipynb:
  raw MFI  ->  per-marker gates  ->  sm.pp.rescale (0..1, 0.5=gate)  ->  sm.tl.phenotype_cells

Two modes around the human-in-the-loop gate refinement:
  seed   Auto-seed per-marker gates (2-component GMM on log1p, crossover -> raw MFI)
         into adata.uns['gates']; set adata.raw + sm.pp.log1p; save the *seeded* h5ad.
         => open it in sm.pl.napariGater to refine the gates, then re-save.
  apply  sm.pp.rescale(adata, gate=adata.uns['gates']) + sm.tl.phenotype_cells(rules)
         => phenotyped h5ad with obs['phenotype'] (+ 'lineage' parent via the rules).

The seed is only a starting point; napariGater refinement by a human is expected
(per-cell auto-thresholding alone was shown to be unreliable on this panel).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scimap as sm
from sklearn.mixture import GaussianMixture

BASE = Path(__file__).resolve().parents[2] / "data"
INP = BASE / "singlecell_protein.h5ad"
SEEDED = BASE / "singlecell_protein_seeded.h5ad"
PHENO = BASE / "singlecell_protein_phenotyped.h5ad"
RULES = Path(__file__).resolve().parent / "phenotype_rules_senior.csv"


def gmm_gate_log(x, min_sep=2.0):
    """Seed gate in the LOG1P scale — the scale sm.pp.rescale compares against
    (rescale log1p-transforms adata.raw.X internally when log=True). 2-comp GMM
    on log1p, threshold at the posterior crossover. Fallback to the 99th
    percentile (log1p) when the marker is effectively unimodal (napariGater needs
    a finite starting value)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 50:
        return float("nan")
    lx = np.log1p(np.clip(x, 0, None))
    sub = lx[np.linspace(0, lx.size - 1, min(50000, lx.size)).astype(int)].reshape(-1, 1)
    try:
        g = GaussianMixture(2, n_init=1, random_state=0, reg_covar=1e-3).fit(sub)
    except Exception:
        return float(np.quantile(lx, 0.99))
    mu = g.means_.ravel(); sd = np.sqrt(g.covariances_.ravel())
    hi = int(np.argmax(mu)); lo = 1 - hi
    if (mu[hi] - mu[lo]) < min_sep * max(0.5 * (sd[hi] + sd[lo]), 1e-6):
        return float(np.quantile(lx, 0.99))
    grid = np.linspace(mu[lo], mu[hi], 256).reshape(-1, 1)
    p = g.predict_proba(grid)[:, hi]
    cr = np.where(p >= 0.5)[0]
    return float(grid[cr[0], 0] if cr.size else mu[hi])


def run_seed(a):
    adata = ad.read_h5ad(a.input or INP)
    markers = list(adata.var_names)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    print(f"[seed] {adata.n_obs:,} cells x {len(markers)} markers (GMM seed gates)")
    vals = []
    for i, m in enumerate(markers):
        gv = gmm_gate_log(X[:, i], a.min_sep)
        vals.append(gv)
        print(f"  {m:18s} seed gate(log1p) = {gv:.3f}")
    gates_df = pd.DataFrame({"gates": vals}, index=pd.Index(markers))
    adata.raw = adata                      # keep raw for napariGater layer='raw'
    adata = sm.pp.log1p(adata)             # working X -> log1p (display/rescale)
    adata.uns["gates"] = gates_df
    adata.write_h5ad(a.out or SEEDED)
    print("[seed] wrote", a.out or SEEDED, "— refine in napariGater, then run `apply`")


def run_apply(a):
    adata = ad.read_h5ad(a.input or SEEDED)
    if "gates" not in adata.uns:
        raise SystemExit("adata.uns['gates'] missing — run `seed` (+ napariGater refine) first.")
    adata = sm.pp.rescale(adata, gate=adata.uns["gates"])
    phenotype = pd.read_csv(RULES)
    adata = sm.tl.phenotype_cells(adata, phenotype=phenotype, label="phenotype")
    adata.write_h5ad(a.out or PHENO)
    print("[apply] phenotype counts:")
    print(adata.obs["phenotype"].value_counts().to_string())
    print("[apply] wrote", a.out or PHENO)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["seed", "apply"])
    ap.add_argument("--input", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--min-sep", type=float, default=2.0)
    a = ap.parse_args(argv)
    (run_seed if a.mode == "seed" else run_apply)(a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
