"""Cortical-medullary zonation + AIRE+ mTEC neighborhood for the thymus
pipeline.

2026-05-21 — algorithm replaced. The previous DBSCAN-on-sparse-marker-seeds
approach (`thymic_zonation_analysis.py.bak.bullseye`) produced small
bullseye centroids and concentric distance-bin zones — wrong for thymic
anatomy where medulla is a broad contiguous central region inside each
lobule. New algorithm uses the post-5-fix cell-type labels directly:

Inputs:
    data/processed/{sample}/{sample}_phenotyped.h5ad  (post-refinement)

Algorithm (per sample):
    1. Source-class assignment from obs['celltype']:
         cortex_source  = {cTEC, DP_thymocyte}
         medulla_source = {mTEC, AIRE_mTEC, Hassalls,
                           SP4_thymocyte, SP8_thymocyte, Treg, DC, pDC}
       Everything else is neutral (counts toward the kNN denominator but
       not toward either numerator).
    2. cuML brute-force kNN (k=50, fallback sklearn) on obsm['spatial'].
    3. Pass-1 local fractions: per cell, cortex_frac_1 = fraction of kNN
       in cortex_source; medulla_frac_1 same for medulla_source.
    4. Pass-2 smoothing on the same kNN graph: average the pass-1
       fractions over each cell's neighbors. Smooths ~150 μm radius
       without smearing lobule boundaries.
    5. Per-cell zone:
         medulla_frac > MEDULLA_THRESHOLD                → medulla
         cortex_frac  > CORTEX_THRESHOLD                 → cortex
         both > CMJ_BOTH_THRESHOLD                       → cortico_medullary_junction
         medulla_frac > CMJ_BOTH_THRESHOLD               → peri_medullary
         otherwise                                        → cortex (default)
    6. Connected components: BFS over the kNN graph restricted to
       `zone == 'medulla'`. Components < MIN_MEDULLA_CELLS are demoted
       to 'peri_medullary'.
    7. Composition + thymocyte-gradient + AIRE+ neighborhood as before.

Outputs:
    data/processed/medulla_zones.csv                 (renamed; rows = lobule-scale CCs)
    data/processed/thymic_zone_composition.csv       (schema unchanged)
    data/processed/thymocyte_maturation_gradient.csv (schema unchanged)
    data/processed/AIRE_neighborhood.csv             (schema unchanged)
    figures/{sample}/07_thymic_zonation/
        medulla_spatial_overview.png    (broad zones, not bullseyes)
        medulla_lobules_overview.png    NEW — CCs colored individually
        zone_composition_heatmap.png
        thymocyte_gradient_per_zone.png
        aire_neighborhood_composition.png
"""
import argparse
import sys
import time
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree

ROOT = _REPO_ROOT
DEFAULT_SAMPLES = ("0041435", "0041534")

# Source classes — driven by post-5-fix celltype labels, not markers.
CORTEX_SOURCE_TYPES = ("cTEC", "DP_thymocyte")
MEDULLA_SOURCE_TYPES = (
    "mTEC", "AIRE_mTEC", "Hassalls",
    "SP4_thymocyte", "SP8_thymocyte", "Treg",
    "DC", "pDC",
)

# kNN-density defaults. The fraction is the proportion of cell-type-
# labeled neighbors in the source class; both are computed against the
# full kNN denominator (so neutral cells dilute both equally).
#
# We use a SIGNED dominance score:
#   score = medulla_frac - cortex_frac
# Thresholding on this score (rather than absolute fractions per class)
# is symmetric and robust to sample-level differences in overall medulla
# vs cortex proportions (e.g., 0041534 has more cortex source cells
# than medulla, so absolute-threshold approaches under-call medulla).
KNN_K                  = 50
SCORE_MEDULLA          = 0.05   # score >  this  → medulla
SCORE_CORTEX           = -0.10  # score <  this  → cortex
PRESENCE_FLOOR         = 0.15   # require max(med, cor) above this to be on-tissue
MIN_MEDULLA_CELLS      = 200
SMOOTH_PASSES          = 2      # 2 = local fraction + averaging step (recommended)

ZONE_ORDER = ("medulla", "cortico_medullary_junction", "peri_medullary", "cortex")
ZONE_PALETTE = {
    "medulla":                    "#d62728",   # red
    "cortico_medullary_junction": "#ff7f0e",   # orange
    "peri_medullary":             "#fdd0a2",   # peach
    "cortex":                     "#1f77b4",   # blue
}

AIRE_NEIGHBOR_RADIUS_UM = 50.0

THYMOCYTE_SUBTYPES = ("DN_thymocyte", "DP_thymocyte",
                       "SP4_thymocyte", "SP8_thymocyte", "Treg")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _has_cuml():
    try:
        import cuml.neighbors  # noqa: F401
        return True
    except Exception:
        return False


def source_class_vectors(adata):
    """Return (cortex_source, medulla_source) bool arrays of length n_obs."""
    label_col = "celltype_lineage" if "celltype_lineage" in adata.obs.columns else "celltype"
    labels = adata.obs[label_col].astype(str).values
    cortex_src  = np.isin(labels, CORTEX_SOURCE_TYPES)
    medulla_src = np.isin(labels, MEDULLA_SOURCE_TYPES)
    print(f"    source classes: cortex={int(cortex_src.sum()):,}, "
          f"medulla={int(medulla_src.sum()):,}, "
          f"neutral={int((~cortex_src & ~medulla_src).sum()):,}")
    return cortex_src, medulla_src, label_col


def knn_indices(coords, k):
    """k+1 nearest neighbors (index 0 is self). cuML brute, sklearn fallback."""
    n = coords.shape[0]
    if _has_cuml():
        try:
            import cuml.neighbors
            print(f"    cuML kNN (k={k+1}, brute) on {n:,} cells")
            nn = cuml.neighbors.NearestNeighbors(n_neighbors=k+1,
                                                  algorithm="brute",
                                                  metric="euclidean")
            nn.fit(coords.astype(np.float32))
            _, idx = nn.kneighbors(coords.astype(np.float32))
            return np.asarray(idx, dtype=np.int64)
        except Exception as e:
            print(f"    cuML kNN failed: {e!r}; falling back to sklearn")
    from sklearn.neighbors import NearestNeighbors
    print(f"    sklearn kNN (k={k+1}) on {n:,} cells")
    nn = NearestNeighbors(n_neighbors=k+1, algorithm="auto",
                          metric="euclidean", n_jobs=-1).fit(coords)
    _, idx = nn.kneighbors(coords)
    return np.asarray(idx, dtype=np.int64)


def knn_density(nn_idx, cortex_src, medulla_src, passes):
    """Pass 1: local fraction of source-class neighbors per cell.
    Pass 2 (optional): smooth those fractions over the same kNN graph.
    Returns (cortex_frac, medulla_frac)."""
    nbrs = nn_idx[:, 1:]  # exclude self
    k = nbrs.shape[1]
    c = cortex_src[nbrs].sum(axis=1) / k
    m = medulla_src[nbrs].sum(axis=1) / k
    if passes >= 2:
        c = c[nbrs].mean(axis=1)
        m = m[nbrs].mean(axis=1)
    return c, m


def assign_zone(cortex_frac, medulla_frac,
                score_med, score_cor, presence_floor):
    """Per-cell zone label vector based on signed dominance score
    `score = medulla_frac - cortex_frac`:
        score > score_med                                  → medulla
        score_cor ≤ score ≤ score_med  AND on-tissue       → cortico_medullary_junction
        score < score_cor                                  → cortex
        max(med, cor) < presence_floor                     → cortex (default; low-signal)
    """
    score = medulla_frac - cortex_frac
    on_tissue = np.maximum(medulla_frac, cortex_frac) >= presence_floor
    zone = np.full(cortex_frac.shape, "cortex", dtype=object)
    zone[(score < score_cor) & on_tissue]                    = "cortex"
    cmj_mask = (score >= score_cor) & (score <= score_med) & on_tissue
    zone[cmj_mask]                                            = "cortico_medullary_junction"
    # peri-medullary: scoring favors medulla but below the medulla gate
    peri_mask = (score > 0) & (score <= score_med) & on_tissue & ~cmj_mask
    zone[peri_mask]                                           = "peri_medullary"
    zone[(score > score_med) & on_tissue]                    = "medulla"
    return zone


def connected_components_on_knn(nn_idx, member_mask, min_size):
    """scipy.sparse.csgraph.connected_components on the kNN graph
    restricted to member_mask.

    Returns a label vector of length n_cells:
        -1 for non-member cells (or members demoted by min_size filter)
        component_id (>= 0) for members in a CC of size >= min_size
    Also returns (n_kept_components, raw_component_sizes).
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    n = nn_idx.shape[0]
    out = np.full(n, -1, dtype=np.int64)
    if not member_mask.any():
        return out, 0, np.array([], dtype=np.int64)

    # Build edges (i → j) for j in i's kNN (excluding self), but only
    # where BOTH endpoints are members. Use symmetric closure via
    # `directed=False` in connected_components.
    nbrs = nn_idx[:, 1:]  # drop self column
    k = nbrs.shape[1]
    row = np.repeat(np.arange(n, dtype=np.int64), k)
    col = nbrs.ravel().astype(np.int64)
    edge_keep = member_mask[row] & member_mask[col]
    row = row[edge_keep]; col = col[edge_keep]
    data = np.ones(row.size, dtype=np.uint8)
    g = csr_matrix((data, (row, col)), shape=(n, n))

    n_cc, raw_labels = connected_components(g, directed=False, return_labels=True)
    # raw_labels has a value for every cell — non-members get a singleton
    # component each (since they have no edges). Restrict to members.
    member_labels = raw_labels[member_mask]
    sizes = np.bincount(member_labels)
    # Reverse-map: cc_id → new sequential id if size >= min_size
    keep = sizes >= min_size
    remap = np.full(sizes.size, -1, dtype=np.int64)
    remap[keep] = np.arange(int(keep.sum()))
    out[member_mask] = remap[member_labels]
    final_n_cc = int(keep.sum())
    # Return only the sizes of components that contain at least one
    # member (filter out non-member singletons to keep diagnostics clean)
    member_cc_ids = np.flatnonzero(sizes > 0)
    return out, final_n_cc, sizes[member_cc_ids]


# -----------------------------------------------------------------------------
# AIRE neighborhood (unchanged)
# -----------------------------------------------------------------------------

def aire_neighborhood(a, coords, label_col):
    """For each AIRE_mTEC cell, count non-self cells within RADIUS_UM
    by lineage. Returns a DataFrame (one row per AIRE_mTEC cell)."""
    aire_mask = a.obs[label_col].astype(str).values == "AIRE_mTEC"
    if aire_mask.sum() == 0:
        print("    no AIRE_mTEC cells; AIRE neighborhood skipped")
        return pd.DataFrame()
    aire_idx = np.where(aire_mask)[0]
    print(f"    AIRE+ mTEC neighborhood: {len(aire_idx):,} centers, "
          f"radius={AIRE_NEIGHBOR_RADIUS_UM} μm")
    tree = cKDTree(coords)
    neighbor_ids = tree.query_ball_point(coords[aire_idx], r=AIRE_NEIGHBOR_RADIUS_UM)
    label = a.obs[label_col].astype(str).values
    rows = []
    for ci, nbrs in zip(aire_idx, neighbor_ids):
        nbrs = [n for n in nbrs if n != ci]
        comp = pd.Series([label[n] for n in nbrs]).value_counts().to_dict()
        rows.append({"cell_id": a.obs.index[ci],
                     "n_within_50um": len(nbrs),
                     **comp})
    return pd.DataFrame(rows).fillna(0)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def fig_dir_for(sample):
    d = ROOT / f"figures/{sample}/07_thymic_zonation"
    d.mkdir(parents=True, exist_ok=True)
    return d


def plot_medulla_spatial(coords, zone, sample, out_path):
    fig, ax = plt.subplots(figsize=(8, 8))
    for z in ZONE_ORDER:
        m = zone == z
        if m.sum() == 0:
            continue
        ax.scatter(coords[m, 0], coords[m, 1], s=0.02,
                   c=ZONE_PALETTE.get(z, "grey"),
                   label=f"{z} ({int(m.sum()):,})",
                   alpha=0.55, rasterized=True)
    ax.set_title(f"{sample} — cortical / medullary zones (kNN density)")
    ax.set_aspect("equal")
    ax.set_xlabel("x (μm)"); ax.set_ylabel("y (μm)")
    ax.legend(markerscale=20, fontsize=8, loc="upper right", framealpha=0.7)
    fig.tight_layout(); fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_medulla_lobules(coords, medulla_id, sample, out_path):
    """One color per medullary connected component; grey background."""
    fig, ax = plt.subplots(figsize=(8, 8))
    # background: non-medulla cells
    bg = medulla_id < 0
    ax.scatter(coords[bg, 0], coords[bg, 1], s=0.02,
               c="#dddddd", alpha=0.35, rasterized=True, label=None)
    # foreground: each lobule a separate color (tab20 cycle)
    cc_ids = np.unique(medulla_id[~bg])
    cmap = cm.get_cmap("tab20", max(20, len(cc_ids)))
    for i, cid in enumerate(sorted(cc_ids.tolist())):
        m = medulla_id == cid
        ax.scatter(coords[m, 0], coords[m, 1], s=0.04,
                   c=[cmap(i % cmap.N)], alpha=0.85, rasterized=True,
                   label=f"lobule {int(cid)} ({int(m.sum()):,})")
    ax.set_title(f"{sample} — medullary lobules ({len(cc_ids)} CCs ≥ {MIN_MEDULLA_CELLS} cells)")
    ax.set_aspect("equal")
    ax.set_xlabel("x (μm)"); ax.set_ylabel("y (μm)")
    if len(cc_ids) <= 20:
        ax.legend(markerscale=20, fontsize=7, loc="upper right",
                  framealpha=0.7, ncol=1)
    fig.tight_layout(); fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_zone_composition(comp_df, sample, out_path):
    pivot = comp_df.pivot_table(index="zone", columns="celltype",
                                 values="frac_zone", aggfunc="mean", fill_value=0.0)
    zone_order = [z for z in ZONE_ORDER if z in pivot.index]
    pivot = pivot.reindex(zone_order)
    pivot = pivot.loc[:, pivot.sum(axis=0).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(min(20, 1.5 + 0.6*pivot.shape[1]), 4.5))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis", vmin=0, vmax=pivot.values.max())
    ax.set_xticks(range(pivot.shape[1])); ax.set_xticklabels(pivot.columns, rotation=60, ha="right")
    ax.set_yticks(range(pivot.shape[0])); ax.set_yticklabels(pivot.index)
    ax.set_title(f"{sample} — composition fraction by zone")
    plt.colorbar(im, ax=ax, label="fraction of cells in zone")
    fig.tight_layout(); fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_thymocyte_gradient(thy_df, sample, out_path):
    pivot = thy_df.pivot_table(index="zone", columns="thymocyte_subtype",
                               values="frac", aggfunc="mean", fill_value=0.0)
    zone_order = [z for z in ZONE_ORDER if z in pivot.index]
    pivot = pivot.reindex(zone_order)
    cols = [c for c in THYMOCYTE_SUBTYPES if c in pivot.columns]
    pivot = pivot[cols]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bottom = np.zeros(pivot.shape[0])
    colors = plt.cm.tab10(np.arange(len(cols)))
    for ci, col in enumerate(cols):
        ax.bar(pivot.index, pivot[col].values, bottom=bottom,
               label=col, color=colors[ci])
        bottom += pivot[col].values
    ax.set_ylabel("fraction of thymocytes in zone")
    ax.set_title(f"{sample} — thymocyte maturation gradient")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    fig.tight_layout(); fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_aire_neighborhood(neigh_df, sample, out_path):
    if neigh_df.empty:
        return
    drop_cols = {"sample", "cell_id", "n_within_50um"}
    comp_cols = [c for c in neigh_df.columns if c not in drop_cols]
    if not comp_cols:
        return
    totals = neigh_df[comp_cols].sum(axis=0)
    totals = pd.to_numeric(totals, errors="coerce").dropna()
    if totals.empty:
        return
    totals = totals.sort_values(ascending=False).head(15)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(range(len(totals)), totals.values, color="#1f77b4")
    ax.set_xticks(range(len(totals))); ax.set_xticklabels(totals.index, rotation=60, ha="right")
    ax.set_ylabel("total cells in any AIRE+ 50 μm neighborhood")
    ax.set_title(f"{sample} — AIRE+ mTEC neighborhood composition "
                 f"({len(neigh_df):,} centers)")
    fig.tight_layout(); fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main per-sample pipeline
# -----------------------------------------------------------------------------

def process_sample(sample, k, score_med, score_cor, presence_floor,
                    min_cells, passes):
    print(f"\n{'='*72}\n=== {sample}\n{'='*72}")
    p = ROOT / f"data/processed/{sample}/{sample}_phenotyped.h5ad"
    print(f"  loading {p.name}...")
    t = time.time()
    a = sc.read_h5ad(p)
    print(f"    loaded {a.shape} in {time.time()-t:.1f}s")

    coords = np.asarray(a.obsm["spatial"], dtype=np.float32)

    # 1. Source classes
    cortex_src, medulla_src, label_col = source_class_vectors(a)
    if medulla_src.sum() == 0:
        print("  no medulla-source cells; aborting zonation")
        return None

    # 2. kNN
    t = time.time()
    nn_idx = knn_indices(coords, k)
    print(f"    kNN done in {time.time()-t:.1f}s, shape={nn_idx.shape}")

    # 3+4. kNN density + zone assignment
    cortex_frac, medulla_frac = knn_density(nn_idx, cortex_src, medulla_src, passes)
    score = medulla_frac - cortex_frac
    print(f"    density (pass={passes}): "
          f"medulla_frac med={np.median(medulla_frac):.3f}  "
          f"p90={np.quantile(medulla_frac, 0.90):.3f}  max={medulla_frac.max():.3f} | "
          f"cortex_frac med={np.median(cortex_frac):.3f}  "
          f"p90={np.quantile(cortex_frac, 0.90):.3f}  max={cortex_frac.max():.3f}")
    print(f"    score (med-cor): "
          f"min={score.min():.3f}  p10={np.quantile(score,0.10):.3f}  "
          f"med={np.median(score):.3f}  p90={np.quantile(score,0.90):.3f}  max={score.max():.3f}")
    zone = assign_zone(cortex_frac, medulla_frac,
                       score_med, score_cor, presence_floor)
    n_med0 = int((zone == "medulla").sum())
    print(f"    pre-CC: medulla={n_med0:,}  CMJ={int((zone == 'cortico_medullary_junction').sum()):,}  "
          f"peri={int((zone == 'peri_medullary').sum()):,}  cortex={int((zone == 'cortex').sum()):,}")

    # 5. Connected components on the medulla mask
    t = time.time()
    med_mask = zone == "medulla"
    cc_labels, n_cc, sizes = connected_components_on_knn(nn_idx, med_mask, min_cells)
    print(f"    CC done in {time.time()-t:.1f}s: {n_cc} medullary lobule(s) ≥ {min_cells} cells "
          f"(raw CCs={len(sizes)})")
    # Demote small-CC members (in med_mask but assigned -1 by the filter)
    demoted = med_mask & (cc_labels < 0)
    if demoted.any():
        zone[demoted] = "peri_medullary"
        print(f"    demoted {int(demoted.sum()):,} small-CC members → peri_medullary")
    n_med_final = int((cc_labels >= 0).sum())
    print(f"    post-CC medulla cells: {n_med_final:,}")

    # Persist on obs for any downstream consumer
    a.obs["medulla_zone"] = pd.Categorical(zone, categories=list(ZONE_ORDER))
    a.obs["medulla_id"] = cc_labels  # -1 = non-medulla

    # Zone metadata (one row per CC)
    zones_meta = []
    for cid in range(n_cc):
        m = cc_labels == cid
        cx, cy = coords[m].mean(axis=0)
        x0, x1 = coords[m, 0].min(), coords[m, 0].max()
        y0, y1 = coords[m, 1].min(), coords[m, 1].max()
        zones_meta.append({"sample": sample,
                           "medulla_id": int(cid),
                           "centroid_x": float(cx),
                           "centroid_y": float(cy),
                           "n_cells": int(m.sum()),
                           "x_extent_um": float(x1 - x0),
                           "y_extent_um": float(y1 - y0)})
    zones_meta_df = pd.DataFrame(zones_meta)

    # 6. Composition (zone × celltype)
    label = a.obs[label_col].astype(str).values
    comp = (pd.DataFrame({"zone": zone, "celltype": label})
            .groupby(["zone", "celltype"], observed=True)
            .size().rename("n_cells").reset_index())
    zone_totals = comp.groupby("zone", observed=True)["n_cells"].transform("sum")
    comp["frac_zone"] = comp["n_cells"] / zone_totals
    comp.insert(0, "sample", sample)
    print(f"    composition rows: {len(comp)}")

    # 7. Thymocyte maturation gradient
    thy_mask = pd.Series(label).isin(THYMOCYTE_SUBTYPES).values
    if thy_mask.any():
        thy_df = (pd.DataFrame({
            "zone": zone[thy_mask],
            "thymocyte_subtype": label[thy_mask],
        })
        .groupby(["zone", "thymocyte_subtype"], observed=True)
        .size().rename("n_cells").reset_index())
        thy_totals = thy_df.groupby("zone", observed=True)["n_cells"].transform("sum")
        thy_df["frac"] = thy_df["n_cells"] / thy_totals
        thy_df.insert(0, "sample", sample)
    else:
        thy_df = pd.DataFrame(columns=["sample", "zone", "thymocyte_subtype",
                                       "n_cells", "frac"])
    print(f"    thymocyte-gradient rows: {len(thy_df)}")

    # 8. AIRE+ neighborhood
    aire_df = aire_neighborhood(a, coords, label_col)
    if not aire_df.empty:
        aire_df.insert(0, "sample", sample)
    print(f"    AIRE neighborhood rows: {len(aire_df)}")

    # Figures
    fd = fig_dir_for(sample)
    print("  rendering figures ...")
    plot_medulla_spatial(coords, zone, sample, fd / "medulla_spatial_overview.png")
    plot_medulla_lobules(coords, cc_labels, sample, fd / "medulla_lobules_overview.png")
    plot_zone_composition(comp, sample, fd / "zone_composition_heatmap.png")
    if not thy_df.empty:
        plot_thymocyte_gradient(thy_df, sample, fd / "thymocyte_gradient_per_zone.png")
    if not aire_df.empty:
        plot_aire_neighborhood(aire_df, sample, fd / "aire_neighborhood_composition.png")
    print(f"  figures saved → {fd}")

    return {"zones_meta": zones_meta_df,
            "composition": comp,
            "thymocyte_gradient": thy_df,
            "aire_neighborhood": aire_df}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main(samples, k, score_med, score_cor, presence_floor, min_cells, passes):
    parts = []
    for s in samples:
        r = process_sample(s, k, score_med, score_cor, presence_floor,
                           min_cells, passes)
        if r is None:
            continue
        parts.append(r)
    if not parts:
        print("no sample produced zones; nothing to write")
        return

    def _concat(key):
        dfs = [p[key] for p in parts if not p[key].empty]
        if not dfs:
            return pd.DataFrame()
        return pd.concat(dfs, ignore_index=True)

    out_dir = ROOT / "data/processed"
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {
        "medulla_zones.csv":                 _concat("zones_meta"),
        "thymic_zone_composition.csv":       _concat("composition"),
        "thymocyte_maturation_gradient.csv": _concat("thymocyte_gradient"),
        "AIRE_neighborhood.csv":             _concat("aire_neighborhood"),
    }
    for name, df in files.items():
        path = out_dir / name
        df.to_csv(path, index=False)
        print(f"  wrote {path}  ({len(df):,} rows)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", default=list(DEFAULT_SAMPLES))
    ap.add_argument("--knn", type=int, default=KNN_K,
                    help=f"k for kNN density (default {KNN_K})")
    ap.add_argument("--score-medulla", type=float, default=SCORE_MEDULLA,
                    help=f"score (med-cor) above this → medulla (default {SCORE_MEDULLA})")
    ap.add_argument("--score-cortex", type=float, default=SCORE_CORTEX,
                    help=f"score below this → cortex (default {SCORE_CORTEX})")
    ap.add_argument("--presence-floor", type=float, default=PRESENCE_FLOOR,
                    help=f"require max(med, cor) ≥ this to score (default {PRESENCE_FLOOR})")
    ap.add_argument("--min-medulla-cells", type=int, default=MIN_MEDULLA_CELLS,
                    help=f"minimum CC size to keep as medulla (default {MIN_MEDULLA_CELLS})")
    ap.add_argument("--smooth-passes", type=int, default=SMOOTH_PASSES, choices=[1, 2],
                    help=f"1 = local fraction only, 2 = also smooth (default {SMOOTH_PASSES})")
    args = ap.parse_args()
    print(f"thymic_zonation_analysis: samples={args.samples}  "
          f"k={args.knn}  passes={args.smooth_passes}  "
          f"score_med>{args.score_medulla}  score_cor<{args.score_cortex}  "
          f"presence>={args.presence_floor}  min_cells={args.min_medulla_cells}")
    main(args.samples, args.knn,
         args.score_medulla, args.score_cortex, args.presence_floor,
         args.min_medulla_cells, args.smooth_passes)
    print("\nDONE")
