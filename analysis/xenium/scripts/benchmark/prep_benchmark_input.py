#!/usr/bin/env python3
"""Step 1 of the spatial-method benchmark: build ONE canonical input per
(sample, variant) consumed IDENTICALLY by the Python (BANKSY/UTAG) and R
(KODAMA/MISTy) runners.

Runs in the `xenium_analysis` env (stable scanpy 1.12 / anndata 0.12). It only
READS the phenotyped h5ad (+ the spatial-analysis sidecar for the niche
baseline, if present) and WRITES new files under data/processed/benchmark/{s}/.
Nothing existing is modified.

Per (sample, variant) it writes:
  {s}_bench_{v}_input.h5ad     genes = var['embedding_panel']; X = lognorm;
                               layers counts+lognorm; obsm['spatial']; obs labels
  {s}_bench_{v}_expr.mtx       sparse lognorm, ROWS=cells x COLS=genes
  {s}_bench_{v}_genes.csv      gene names (column order of the MTX)
  {s}_bench_{v}_coords.csv     x,y microns (row order of the MTX)
  {s}_bench_{v}_meta.csv       cell_id, celltype, celltype_broad, celltype_lineage,
                               (spatial_niche if available)
  {s}_bench_{v}_barcodes.csv   cell_id (row order) — fairness/parity guard
  {s}_bench_{v}_misty_markers.txt   tissue MARKER_DICT union ∩ panel
  {s}_bench_{v}_cropinfo.json  n, mode, anchor, extent, celltype composition

Subsampling (seed-pinned, decision: preserve spatial structure):
  crop<N>k : the N cells nearest the robust spatial centre (median x,y) — a
             compact contiguous region that preserves the micron neighbourhood
             graph the radius/kNN methods depend on. Composition is REPORTED
             (no proportion gate — Tier-C).
  thin<N>k : a uniform-random N cells (sampling-robustness control).
  full     : all cells.

Usage:
  python scripts/benchmark/prep_benchmark_input.py --sample 0041323 --variants crop50k thin50k
  python scripts/benchmark/prep_benchmark_input.py --sample 0041435 --variants crop50k
"""
from __future__ import annotations

import argparse
import importlib
import json
import sys
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[2]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[4]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.benchmark import bench_common as bc  # noqa: E402


def _load_backed(sample: str):
    import scanpy as sc
    p = bc.DATA_ROOT / f"processed/{sample}/{sample}_phenotyped.h5ad"
    if not p.exists():
        raise FileNotFoundError(p)
    return sc.read_h5ad(p, backed="r")


def _niche_baseline(sample: str, index: pd.Index) -> pd.Series | None:
    """Pull the stage-03 spatial_niche, aligned to `index`, if the sidecar exists."""
    import anndata as ad
    p = bc.DATA_ROOT / f"processed/{sample}/{sample}_spatial_analysis.h5ad"
    if not p.exists():
        return None
    try:
        a = ad.read_h5ad(p, backed="r")
        if "spatial_niche" not in a.obs.columns:
            return None
        s = a.obs["spatial_niche"].astype(str)
        return s.reindex(index)
    except Exception as e:  # noqa: BLE001
        print(f"  [warn] could not read spatial_niche: {e}")
        return None


def _choose_idx(xy: np.ndarray, variant: str) -> np.ndarray:
    n = xy.shape[0]
    target = bc.variant_to_n(variant)
    mode = bc.variant_mode(variant)
    if mode == "full" or target is None or target >= n:
        return np.arange(n)
    if mode == "thin":
        rng = np.random.default_rng(bc.SEED)
        idx = rng.choice(n, size=target, replace=False)
        return np.sort(idx)
    # crop: N nearest the robust spatial centre (deterministic, contiguous disc)
    anchor = np.median(xy, axis=0)
    d2 = ((xy - anchor[None, :]) ** 2).sum(axis=1)
    idx = np.argsort(d2, kind="stable")[:target]
    return np.sort(idx)


def _misty_markers(tissue: str, genes_present: set[str]) -> list[str]:
    """Union of the tissue MARKER_DICT panels ∩ exported genes (keeps MISTy's
    per-target-marker RF cost tractable on the 5k panel)."""
    mod = importlib.import_module(f"scripts.tissue_markers.{tissue}")
    union: list[str] = []
    seen = set()
    for genes in mod.MARKER_DICT.values():
        for g in genes:
            if g in genes_present and g not in seen:
                seen.add(g)
                union.append(g)
    return union


def prep_one(sample: str, tissue: str, variant: str, gene_panel: str) -> None:
    print(f"\n=== prep {sample} ({tissue}) :: {variant}  [panel={gene_panel}] ===")
    paths = bc.input_paths(sample, variant)
    paths["h5ad"].parent.mkdir(parents=True, exist_ok=True)

    adata = _load_backed(sample)
    xy_all = np.asarray(adata.obsm["spatial"])
    idx = _choose_idx(xy_all, variant)
    print(f"  cells: {len(idx):,} / {adata.n_obs:,}  ({bc.variant_mode(variant)})")

    # Materialise the subset, then restrict to the chosen curated gene panel
    # (default embedding_panel = the data-driven lineage-discriminating set,
    # decision #2 — NOT HVG-variance ranking). All methods see identical genes.
    sub = adata[idx].to_memory()
    if gene_panel not in sub.var.columns:
        raise KeyError(f"{sample}: var['{gene_panel}'] missing")
    panel = sub.var[gene_panel].to_numpy().astype(bool)
    print(f"  genes: {int(panel.sum()):,} (var['{gene_panel}'])")
    sub = sub[:, panel].copy()
    # X := lognorm (the normalised space the gene-based methods consume).
    if "lognorm" not in sub.layers:
        raise KeyError(f"{sample}: layers['lognorm'] missing")
    sub.X = sub.layers["lognorm"].copy()
    # Keep only what the runners need; drop heavy obsm/uns to shrink the file.
    keep_obsm = {k: sub.obsm[k] for k in ("spatial", "X_scvi") if k in sub.obsm}
    keep_obs_cols = [c for c in ("cell_id", "celltype", "celltype_broad",
                                 "celltype_lineage", "total_counts", "cell_area",
                                 "donor_sex") if c in sub.obs.columns]
    import anndata as ad
    out = ad.AnnData(
        X=sub.X,
        obs=sub.obs[keep_obs_cols].copy(),
        var=sub.var[[c for c in ("embedding_panel", "celltype_panel",
                                  "panel_for_scoring") if c in sub.var.columns]].copy(),
        layers={"counts": sub.layers["counts"].copy(),
                "lognorm": sub.layers["lognorm"].copy()},
        obsm=keep_obsm,
    )
    out.obs_names = sub.obs_names
    out.var_names = sub.var_names

    # Niche baseline (optional reference; the head-to-head recomputes it anyway).
    niche = _niche_baseline(sample, out.obs_names)
    if niche is not None:
        out.obs["spatial_niche_ref"] = pd.Categorical(niche.values)
        print(f"  joined spatial_niche reference ({niche.notna().sum():,} non-null)")

    # ---- write h5ad (Python methods) ----
    out.write_h5ad(paths["h5ad"])

    # ---- write R bridge (identical cell/gene order, by construction) ----
    Xln = out.layers["lognorm"]
    Xln = sparse.csr_matrix(Xln) if not sparse.issparse(Xln) else Xln.tocsr()
    mmwrite(str(paths["expr_mtx"]), Xln)  # ROWS=cells, COLS=genes
    pd.Series(out.var_names, name="gene").to_csv(paths["genes"], index=False)
    cid = (out.obs["cell_id"].astype(str).to_numpy()
           if "cell_id" in out.obs.columns else out.obs_names.to_numpy())
    pd.DataFrame({"cell_id": cid}).to_csv(paths["barcodes"], index=False)
    xy = np.asarray(out.obsm["spatial"])
    pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}).to_csv(paths["coords"], index=False)
    meta_cols = {"cell_id": cid}
    for c in ("celltype", "celltype_broad", "celltype_lineage", "spatial_niche_ref"):
        if c in out.obs.columns:
            meta_cols[c] = out.obs[c].astype(str).to_numpy()
    pd.DataFrame(meta_cols).to_csv(paths["meta"], index=False)

    # ---- MISTy marker subset ----
    markers = _misty_markers(tissue, set(map(str, out.var_names)))
    paths["misty_markers"].write_text("\n".join(markers) + "\n")
    print(f"  MISTy markers: {len(markers)}")

    # ---- cropinfo (composition REPORTED, never gated) ----
    comp = (out.obs["celltype"].astype(str).value_counts()
            if "celltype" in out.obs.columns else pd.Series(dtype=int))
    info = {
        "sample": sample, "tissue": tissue, "variant": variant,
        "n_cells": int(out.n_obs), "n_genes": int(out.n_vars),
        "gene_panel": gene_panel,
        "mode": bc.variant_mode(variant), "seed": bc.SEED,
        "anchor_median_xy": [float(v) for v in np.median(xy_all, axis=0)],
        "extent_um": {"x": [float(xy[:, 0].min()), float(xy[:, 0].max())],
                      "y": [float(xy[:, 1].min()), float(xy[:, 1].max())]},
        "celltype_composition": {k: int(v) for k, v in comp.items()},
        "n_misty_markers": len(markers),
    }
    paths["cropinfo"].write_text(json.dumps(info, indent=2))

    # ---- parity assert (the fairness guard) ----
    assert len(cid) == out.n_obs == xy.shape[0], "row-count mismatch"
    assert Xln.shape == (out.n_obs, out.n_vars), "MTX shape mismatch"
    print(f"  wrote -> {paths['h5ad'].parent}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--tissue", choices=("pancreas", "thymus"))
    ap.add_argument("--variants", nargs="+", default=[bc.HEADLINE_VARIANT])
    ap.add_argument("--gene-panel", default="embedding_panel",
                    choices=("embedding_panel", "panel_for_embedding",
                             "celltype_panel", "panel_for_scoring"),
                    help="var mask selecting the gene set all methods consume")
    args = ap.parse_args()
    tissue = args.tissue or bc.SAMPLE_TISSUE.get(args.sample)
    if tissue is None:
        ap.error(f"unknown tissue for sample {args.sample}; pass --tissue")
    for v in args.variants:
        prep_one(args.sample, tissue, v, args.gene_panel)


if __name__ == "__main__":
    main()
