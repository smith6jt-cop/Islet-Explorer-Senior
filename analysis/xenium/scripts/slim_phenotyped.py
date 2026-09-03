"""Write slim versions of phenotyped h5ads to disk so 05 integration can load
both into the 240 GB cgroup without OOM during concat.

Keeps only what's needed for 05:
  - X + layers['counts']     (scVI needs counts; X preserved for downstream)
  - obs (celltype, sample_id, QC cols, niche if present)
  - var
  - obsm['spatial']           (for spatial ops later)
  - obsm['X_umap']            (per-sample UMAP for QC)

Drops: all obsp, scVI/PCA obsm, log_normalized/scaled layers, heavy uns.

Usage: python scripts/slim_phenotyped.py 0041323 0041326
"""
import gc
import sys
import os
from pathlib import Path

import scanpy as sc

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
PROCESSED_ROOT = DATA_ROOT / "processed"

DROP_OBSM  = {"X_pca", "X_scvi", "X_umap",
              "X_scvi_broad", "X_umap_broad"}  # 05 rebuilds UMAP; drop all embeddings
DROP_OBSP  = {"connectivities", "distances", "scvi_connectivities",
              "scvi_distances", "spatial_connectivities", "spatial_distances"}
DROP_LAYER = {"log_normalized", "scaled"}
DROP_UNS   = {"neighbors", "scvi", "pca", "umap",
              "_scvi_manager_uuid", "_scvi_uuid",
              "dendrogram_leiden_1.5", "leiden_1.5_nhood_enrichment",
              "leiden_1.5_centrality_scores",
              "rank_genes_leiden", "hvg", "log1p", "spatial_neighbors",
              "celltype_nhood_enrichment", "co_occurrence", "moranI",
              "ripley_celltype", "celltype_ligrec"}

def slim(sample: str) -> Path:
    in_path = PROCESSED_ROOT / sample / f"{sample}_phenotyped.h5ad"
    out_path = PROCESSED_ROOT / sample / f"{sample}_phenotyped_slim.h5ad"
    if not in_path.exists():
        raise FileNotFoundError(in_path)
    print(f"[{sample}] loading {in_path}")
    ad = sc.read_h5ad(in_path)
    print(f"[{sample}] shape={ad.shape}; slimming ...")
    for k in list(ad.obsm.keys()):
        if k in DROP_OBSM: del ad.obsm[k]
    for k in list(ad.obsp.keys()):
        if k in DROP_OBSP: del ad.obsp[k]
    for k in list(ad.layers.keys()):
        if k in DROP_LAYER: del ad.layers[k]
    for k in list(ad.uns.keys()):
        if k in DROP_UNS: del ad.uns[k]
    # Replace X with counts layer (smaller: raw counts sparse integer instead of scaled float)
    # and drop the scaled 'X' to halve per-sample memory at load time.
    if 'counts' in ad.layers:
        ad.X = ad.layers['counts']  # view, no copy
        del ad.layers['counts']     # avoid double storage
    # Sanitize any remaining uns dicts with tuple keys
    for k in list(ad.uns.keys()):
        v = ad.uns[k]
        if isinstance(v, dict) and any(not isinstance(ki, str) for ki in v.keys()):
            del ad.uns[k]
    ad.obs['sample_id'] = sample
    gc.collect()
    print(f"[{sample}] writing {out_path}")
    ad.write_h5ad(out_path)
    del ad
    gc.collect()
    print(f"[{sample}] slim size: {out_path.stat().st_size / 1024**3:.1f} GB")
    return out_path


if __name__ == "__main__":
    samples = sys.argv[1:] if len(sys.argv) > 1 else ["0041323", "0041326"]
    for s in samples:
        slim(s)
    print("All slims done.")
