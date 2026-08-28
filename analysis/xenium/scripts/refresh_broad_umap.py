#!/usr/bin/env python3
"""Refresh the DIAGNOSTIC broad-tier UMAP from the already-stored X_scvi_broad.

The broad UMAP (uns['umap_broad'] = {index, coords}, a 200k subsample) is
diagnostic-only — downstream (stage 03 / insulitis / zonation) drops the broad
embedding; only the `celltype_broad` LABEL is load-bearing. When the broad UMAP
reuses the fine/flat UMAP params (min_dist=0.05, spread=2.5, k=10) it over-
compresses the ~8-lineage latent into a blob + detached rings. This helper
recomputes that UMAP from the existing `obsm['X_scvi_broad']` with the looser
BROAD_* params, persists a small `.npz` sidecar next to the h5ad (NEVER rewrites
the 60 GB h5ad, so labels are untouched), and re-renders the broad-UMAP figures.

It does NOT retrain scANVI. The subsample index is reproduced with the same RNG
+ seed as hierarchical_phenotype.py so it aligns with obs in every renderer.

Run from a rapids-aware shell (GPU neighbors — else silent CPU fallback):
    conda activate xenium_analysis && source scripts/env_rapids.sh
    python scripts/refresh_broad_umap.py                       # all 4, locked params
    python scripts/refresh_broad_umap.py --samples 0041323     # one sample
    python scripts/refresh_broad_umap.py --sweep 0041323       # param grid, render-only
"""
from __future__ import annotations

import argparse
import sys
import time
from importlib import import_module
import os
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import scanpy as sc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

H = import_module("scripts.hierarchical_phenotype")          # scvi/torch imports are lazy → cheap
RS2 = import_module("scripts.render_stage02_phenotyping")

DEFAULT_TISSUE = {"0041323": "pancreas", "0041326": "pancreas",
                  "0041435": "thymus", "0041534": "thymus"}
DEFAULT_SAMPLES = tuple(DEFAULT_TISSUE)
SWEEP_GRID = [(0.3, 30), (0.3, 50), (0.5, 30), (0.5, 50)]   # (min_dist, n_neighbors); spread=1.0
SWEEP_SUBSAMPLE = 50_000   # sweep is a render-only param PICK → smaller subsample is fast
                           # and shows the same blob-vs-separated / spots effect; the full
                           # refresh below uses the real 200k UMAP_FIG_SUBSAMPLE index.


def h5ad_path(sid):
    return DATA_ROOT / f"processed/{sid}/{sid}_phenotyped.h5ad"


def sidecar_path(sid):
    return DATA_ROOT / f"processed/{sid}/{sid}_umap_broad.npz"


def figdir(sid):
    d = DATA_ROOT / f"figures/{sid}/02_phenotyping"
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_latent_subsample(sid, n_sub=None):
    """Read X_scvi_broad (backed → obsm/obs load to memory, X stays on disk) and
    reproduce the exact diagnostic subsample index (same RNG + seed as
    hierarchical_phenotype.py:713-714). n_sub defaults to UMAP_FIG_SUBSAMPLE (200k);
    the sweep passes a smaller value for speed."""
    a = ad.read_h5ad(h5ad_path(sid), backed="r")
    n = int(a.n_obs)
    n_sub = min(H.UMAP_FIG_SUBSAMPLE if n_sub is None else n_sub, n)
    sidx = np.sort(np.random.default_rng(H.RANDOM_STATE).choice(n, n_sub, replace=False))
    lat = np.asarray(a.obsm["X_scvi_broad"][sidx], np.float32)
    obs = a.obs.copy()
    try:
        a.file.close()
    except Exception:
        pass
    return n, sidx, lat, obs


def recompute_umap(lat_sub, n_neighbors, min_dist, spread, rsc):
    usub = ad.AnnData(X=lat_sub)
    usub.obsm["X_scvi_broad"] = lat_sub
    H.build_neighbors(usub, "X_scvi_broad", "broad", rsc, n_neighbors=n_neighbors)
    t0 = time.time()
    sc.tl.umap(usub, neighbors_key="broad", min_dist=min_dist, spread=spread,
               random_state=H.RANDOM_STATE)
    print(f"      UMAP[broad] CPU on {lat_sub.shape[0]:,}-subsample: "
          f"{time.time() - t0:.1f}s", flush=True)
    return np.asarray(usub.obsm["X_umap"], np.float32)


def refresh_one(sid, n_neighbors, min_dist, spread, rsc):
    tissue = DEFAULT_TISSUE[sid]
    mod = RS2.load_tissue_module(tissue)
    print(f"[{sid}] ({tissue}) refresh broad UMAP "
          f"min_dist={min_dist} spread={spread} k={n_neighbors}", flush=True)
    n, sidx, lat, obs = load_latent_subsample(sid)
    coords = recompute_umap(lat, n_neighbors, min_dist, spread, rsc)

    # 1) sidecar — tiny, atomic, carries its provenance; h5ad is never rewritten.
    out = sidecar_path(sid)
    np.savez(out, index=sidx.astype(np.int64), coords=coords,
             min_dist=np.float64(min_dist), spread=np.float64(spread),
             n_neighbors=np.int64(n_neighbors), random_state=np.int64(H.RANDOM_STATE),
             n_obs=np.int64(n))
    print(f"  wrote sidecar {out}", flush=True)

    # 2) hierarchical_umap.png — feed the new coords through the canonical renderer
    #    via a minimal obs-only AnnData (render_figure reads obs.iloc[idx] + uns).
    padata = ad.AnnData(obs=obs)
    padata.uns["umap_broad"] = {"index": sidx.astype(np.int64), "coords": coords}
    H.render_figure(padata, sid, figdir(sid), mod)

    # 3) latent + lineage figures — Sample() picks up the sidecar we just wrote.
    S = RS2.Sample(sid, tissue)
    RS2.fig_latent_diagnostic(S, figdir(sid))
    RS2.fig_lineage_panels(S, figdir(sid))
    print(f"[{sid}] done", flush=True)


def sweep_one(sid, rsc):
    """Render-only param grid on ONE sample (no sidecar). Compare visually,
    then lock the winner into the BROAD_* constants in hierarchical_phenotype.py."""
    tissue = DEFAULT_TISSUE[sid]
    mod = RS2.load_tissue_module(tissue)
    n, sidx, lat, obs = load_latent_subsample(sid, SWEEP_SUBSAMPLE)
    cb = obs["celltype_broad"].astype(str).to_numpy()[sidx]
    cats = sorted(set(cb))
    pal = RS2.palette_for(mod, cats, broad=True)
    rng = np.random.default_rng(H.RANDOM_STATE)
    order = rng.permutation(len(cb))   # de-bias draw order

    fig, axes = plt.subplots(2, 2, figsize=(22, 20))
    for ax, (md, nn) in zip(axes.ravel(), SWEEP_GRID):
        xy = recompute_umap(lat, nn, md, 1.0, rsc)
        for c in cats:
            m = (cb == c)[order]
            ax.scatter(xy[order][m, 0], xy[order][m, 1], s=1.0, c=[pal[c]],
                       label=f"{c} ({m.sum():,})", rasterized=True, linewidths=0)
        ax.set_title(f"{sid} — min_dist={md}, k={nn}, spread=1.0", fontsize=16)
        ax.set_xticks([]); ax.set_yticks([])
    axes[0, 0].legend(markerscale=10, fontsize=11, loc="center left",
                      bbox_to_anchor=(-0.55, 0.5), frameon=False)
    plt.tight_layout()
    out = figdir(sid) / "broad_umap_param_sweep.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[{sid}] wrote sweep grid {out}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", nargs="+", default=list(DEFAULT_SAMPLES))
    ap.add_argument("--sweep", metavar="SAMPLE", default=None,
                    help="render-only param grid on one sample (no sidecar)")
    ap.add_argument("--min-dist", type=float, default=H.BROAD_UMAP_MIN_DIST)
    ap.add_argument("--spread", type=float, default=H.BROAD_UMAP_SPREAD)
    ap.add_argument("--n-neighbors", type=int, default=H.BROAD_NEIGHBORS_K)
    ap.add_argument("--cpu", action="store_true", help="force CPU neighbors")
    args = ap.parse_args()

    rsc = None if args.cpu else H.try_import_rsc()
    if args.sweep:
        sweep_one(args.sweep, rsc)
        return
    for sid in args.samples:
        refresh_one(sid, args.n_neighbors, args.min_dist, args.spread, rsc)


if __name__ == "__main__":
    main()
