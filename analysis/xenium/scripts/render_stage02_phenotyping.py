#!/usr/bin/env python3
"""Consolidated stage-2 (hierarchical scANVI phenotyping) figure renderer.

Reads the clean ``{sample}_phenotyped.h5ad`` (per CLAUDE.md decision #16 schema:
``obsm['X_scvi_broad']`` + ``obsm['spatial']``, broad UMAP in
``uns['umap_broad']`` = {index, coords}, fine call in ``obs['celltype']``, broad
in ``obs['celltype_broad']``, plus per-type ``{ct}_score`` and the confidence
columns) and renders the figure set that tells the
**preprocessed -> scANVI -> cell types** story:

    marker_dotplot.png            tight markers x celltype (mean expr + % cells)
    marker_heatmap.png            z-scored celltype x marker means
    scvi_latent_diagnostic.png    X_scvi_broad dim usage + UMAP by depth/area/lineage
    broad_to_fine_composition.png fine-within-broad stacked bars + abundance
    lineage_umap_panels.png       per-broad-lineage UMAP small-multiples (Tier-2)
    spatial_celltypes_crop.png    a CROPPED, readable spatial cell-type map
    confidence_summary.png        broad/subtype confidence violins per celltype
    marker_panel_coverage.png     #on-panel markers + detection per celltype (fine)
    subpanel_coverage.png         broad-call x subpanel detected-fraction matrix + bars
    program_offpanel_coverage.png off-subpanel genes annotated by pathway/program (ANNOT)

Design rules (deliberate, to avoid the usual small-font / overlapping-legend
problems): base font 14pt, titles 18pt; every categorical legend is anchored
OUTSIDE the axes (``bbox_to_anchor``) and saved with ``bbox_inches='tight'`` so
nothing clips; scatter collections are rasterized; the spatial map is a crop
(the full 1.2M-cell tissue render is unreadable).

Usage
-----
    python scripts/render_stage02_phenotyping.py                 # all 4 samples
    python scripts/render_stage02_phenotyping.py 0041323         # one sample
    python scripts/render_stage02_phenotyping.py --tissue thymus 0041435
    python scripts/render_stage02_phenotyping.py --only dotplot heatmap
"""
import argparse
import importlib
import sys
import os
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

DEFAULT_TISSUE = {"0041323": "pancreas", "0041326": "pancreas",
                  "0041435": "thymus", "0041534": "thymus"}
RS = 0
STRAT_CAP = 8000          # per-celltype cap for the expression subsample (keeps
                          # rare types whole, bounds memory/plot cost)
LOG_PRESENT = float(np.log1p(1.0))

# ---- global style: generous fonts, no clipping ----------------------------
plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "font.size": 14,
    "axes.titlesize": 18,
    "axes.titleweight": "bold",
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "legend.title_fontsize": 13,
    "figure.facecolor": "white",
})

GREY = "#d9d9d9"

# subpanel-coverage (broad) — per-cell own-subpanel detected-fraction segments:
#   {==0}, (0, 5%], (5%, 15%], (>15%).  Only the interior thresholds are listed; the
#   ==0 and >last segments are implicit. The subpanel pools are 168-991 genes, so a
#   COUNT-based bin (as in fig_panel_coverage) saturates instantly — fraction bins keep
#   signal. Tunable.
SUBPANEL_FRAC_EDGES = [0.05, 0.15]
SUBPANEL_BIN_LABELS = ["0", "(0,5%]", "(5,15%]", ">15%"]
SUBPANEL_BIN_COLORS = ["#e0e0e0", "#bdd7e7", "#6baed6", "#08519c"]   # grey->Blues


# ===========================================================================
# helpers
# ===========================================================================
def load_tissue_module(tissue):
    return importlib.import_module(f"scripts.tissue_markers.{tissue}")


def palette_for(mod, cats, broad=False):
    """Return {category: hex}, falling back to tab20 for anything the tissue
    palette omits so scanpy/our scatter never KeyError on a category."""
    base = getattr(mod, "BROAD_PALETTE" if broad else "PALETTE", {}) or {}
    cats = [str(c) for c in cats]
    out, extra = {}, []
    for c in cats:
        if c in base:
            out[c] = base[c]
        else:
            extra.append(c)
    if extra:
        cmap = plt.get_cmap("tab20")
        for i, c in enumerate(sorted(extra)):
            out[c] = matplotlib.colors.to_hex(cmap(i % 20))
    return out


def fine_order(mod, present):
    """Fine celltypes ordered by broad lineage (so related types are adjacent),
    restricted to those actually present."""
    present = list(present)
    b2f = getattr(mod, "BROAD_TO_FINE", {}) or {}
    ordered, seen = [], set()
    for broad in b2f:
        for ct in sorted(b2f[broad]):
            if ct in present and ct not in seen:
                ordered.append(ct); seen.add(ct)
    for ct in sorted(present):              # rare/post-hoc types (Gamma, Epsilon, Hassalls…)
        if ct not in seen:
            ordered.append(ct); seen.add(ct)
    return ordered


def marker_table(mod, var_names):
    """Ordered unique on-panel marker list + per-gene owning lineage block, from
    MARKER_DICT, ordered by broad lineage."""
    md = mod.MARKER_DICT
    present = set(var_names)
    genes, owner, seen = [], [], set()
    for ct in fine_order(mod, md.keys()):
        for g in md.get(ct, []):
            if g in present and g not in seen:
                genes.append(g); owner.append(ct); seen.add(g)
    return genes, owner


def outside_legend(ax, handles, labels, title=None, ncol=1, markerscale=2.0):
    leg = ax.legend(handles, labels, title=title, loc="upper left",
                    bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0,
                    frameon=False, ncol=ncol, markerscale=markerscale,
                    handletextpad=0.4, labelspacing=0.4)
    return leg


def cat_handles(cats, pal, marker="o", size=9):
    return [Line2D([0], [0], marker=marker, color="none",
                   markerfacecolor=pal[str(c)], markeredgecolor="none",
                   markersize=size, label=str(c)) for c in cats]


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    print(f"  [fig] wrote {path.name}")


# ===========================================================================
# data container
# ===========================================================================
class Sample:
    def __init__(self, sid, tissue):
        self.sid, self.tissue = sid, tissue
        self.mod = load_tissue_module(tissue)
        path = DATA_ROOT / f"processed/{sid}/{sid}_phenotyped.h5ad"
        a = ad.read_h5ad(path, backed="r")
        self.obs = a.obs.copy()
        self.var_names = list(a.var_names)
        self.var = a.var.copy()
        self.spatial = np.asarray(a.obsm["spatial"])
        self.xscvi = np.asarray(a.obsm["X_scvi_broad"], np.float32)
        ub = a.uns.get("umap_broad")
        self.umap_idx = np.asarray(ub["index"]) if ub is not None else None
        self.umap_xy = np.asarray(ub["coords"], np.float32) if ub is not None else None
        # Sidecar override: scripts/refresh_broad_umap.py recomputes the diagnostic
        # broad UMAP (looser params) into {sid}_umap_broad.npz WITHOUT rewriting the
        # h5ad. Prefer it when present; fall back to the in-h5ad uns['umap_broad'].
        side = DATA_ROOT / f"processed/{sid}/{sid}_umap_broad.npz"
        if side.exists():
            z = np.load(side)
            self.umap_idx = np.asarray(z["index"])
            self.umap_xy = np.asarray(z["coords"], np.float32)
        self.n = a.n_obs

        # stratified expression subsample (lognorm), capped per celltype
        ct = self.obs["celltype"].astype(str).to_numpy()
        rng = np.random.default_rng(RS)
        keep = []
        for c in np.unique(ct):
            ix = np.where(ct == c)[0]
            if len(ix) > STRAT_CAP:
                ix = rng.choice(ix, STRAT_CAP, replace=False)
            keep.append(ix)
        strat = np.sort(np.concatenate(keep))
        sub = a[strat].to_memory()
        a.file.close()
        self.expr = sub                        # AnnData, lognorm in layers
        self.expr_ct = ct[strat]
        self.fine_present = fine_order(self.mod, np.unique(ct))
        self.pal_fine = palette_for(self.mod, self.fine_present)
        self.broad_present = list(pd.unique(self.obs["celltype_broad"].astype(str)))
        self.pal_broad = palette_for(self.mod, self.broad_present, broad=True)

    def lognorm(self, genes):
        idx = [self.expr.var_names.get_loc(g) for g in genes]
        X = self.expr.layers["lognorm"][:, idx]
        return np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)


# ===========================================================================
# figures
# ===========================================================================
def fig_dotplot(S, figdir):
    genes, owner = marker_table(S.mod, S.var_names)
    cts = S.fine_present
    X = S.lognorm(genes)
    ct = S.expr_ct
    frac = np.zeros((len(cts), len(genes)))
    mean = np.zeros((len(cts), len(genes)))
    for i, c in enumerate(cts):
        m = ct == c
        if not m.any():
            continue
        block = X[m]
        frac[i] = (block > 0).mean(0)
        mean[i] = block.mean(0)
    # scale mean per gene to [0,1] for color comparability
    mn, mx = mean.min(0), mean.max(0)
    scaled = np.where(mx > mn, (mean - mn) / (mx - mn + 1e-9), 0.0)

    cell_w, cell_h = 0.30, 0.42
    fig, ax = plt.subplots(figsize=(max(10, len(genes) * cell_w + 4),
                                    max(6, len(cts) * cell_h + 2.5)))
    gx, gy = np.meshgrid(np.arange(len(genes)), np.arange(len(cts)))
    sizes = (frac * 220) + 4
    sc = ax.scatter(gx.ravel(), gy.ravel(), s=sizes.ravel(),
                    c=scaled.ravel(), cmap="Reds", vmin=0, vmax=1,
                    edgecolors="#bbbbbb", linewidths=0.25)
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=10)
    ax.set_yticks(np.arange(len(cts)))
    ax.set_yticklabels(cts, fontsize=12)
    ax.set_xlim(-0.5, len(genes) - 0.5)
    ax.set_ylim(-0.5, len(cts) - 0.5)
    ax.invert_yaxis()
    ax.set_title(f"{S.sid} — marker expression by cell type "
                 f"(color = scaled mean lognorm, size = % cells)")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    # colorbar
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("scaled mean expression", fontsize=12)
    # size legend
    sl = [Line2D([0], [0], marker="o", color="none", markerfacecolor="#999999",
                 markersize=np.sqrt((f * 220) + 4), label=f"{int(f*100)}%")
          for f in (0.1, 0.25, 0.5, 1.0)]
    ax.legend(handles=sl, title="% cells\nexpressing", loc="upper left",
              bbox_to_anchor=(1.06, 1.0), frameon=False, labelspacing=1.3,
              borderaxespad=0.0)
    save(fig, figdir / "marker_dotplot.png")


def fig_heatmap(S, figdir):
    genes, owner = marker_table(S.mod, S.var_names)
    cts = S.fine_present
    X = S.lognorm(genes)
    ct = S.expr_ct
    mat = np.zeros((len(cts), len(genes)))
    for i, c in enumerate(cts):
        m = ct == c
        if m.any():
            mat[i] = X[m].mean(0)
    # z-score per gene (column)
    z = (mat - mat.mean(0)) / (mat.std(0) + 1e-9)
    z = np.clip(z, -2.5, 2.5)

    fig, ax = plt.subplots(figsize=(max(10, len(genes) * 0.30 + 4),
                                    max(6, len(cts) * 0.42 + 2.5)))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, fontsize=10)
    ax.set_yticks(np.arange(len(cts)))
    ax.set_yticklabels(cts, fontsize=12)
    # lineage block separators (where owner changes)
    for j in range(1, len(owner)):
        if owner[j] != owner[j - 1]:
            ax.axvline(j - 0.5, color="#333333", lw=0.6, alpha=0.5)
    ax.set_title(f"{S.sid} — z-scored mean marker expression by cell type")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("z-score (per gene)", fontsize=12)
    save(fig, figdir / "marker_heatmap.png")


def _umap_scatter(ax, xy, values, pal=None, continuous=False, size=None,
                  cmap="viridis", title=""):
    n = xy.shape[0]
    size = size if size is not None else max(1.5, 5e5 / max(n, 1))
    if continuous:
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=values, s=size, cmap=cmap,
                        linewidths=0, rasterized=True)
        return sc
    for c in pd.unique(values):
        m = values == c
        ax.scatter(xy[m, 0], xy[m, 1], s=size, c=pal[str(c)], linewidths=0,
                   rasterized=True, label=str(c))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return None


def fig_latent_diagnostic(S, figdir):
    if S.umap_idx is None:
        print("  [fig] latent_diagnostic skipped (no umap_broad)"); return
    idx = S.umap_idx
    xy = S.umap_xy
    obs = S.obs.iloc[idx]
    fig, axes = plt.subplots(2, 2, figsize=(20, 17))

    # A: per-dim variance of X_scvi_broad (how much each latent dim is used)
    var = S.xscvi.var(0)
    axA = axes[0, 0]
    axA.bar(np.arange(1, len(var) + 1), var, color="#4c72b0")
    axA.set_xlabel("scANVI latent dimension")
    axA.set_ylabel("variance across cells")
    axA.set_title("X_scvi_broad — per-dimension usage")
    axA.set_xticks(np.arange(1, len(var) + 1))

    # B: UMAP by total_counts (depth — descriptive, NOT a gate)
    depth = obs["total_counts"].to_numpy().astype(float)
    scB = _umap_scatter(axes[0, 1], xy, np.log1p(depth), continuous=True,
                        cmap="magma", title="broad UMAP — log1p(total_counts) [depth]")
    fig.colorbar(scB, ax=axes[0, 1], fraction=0.04, pad=0.02, label="log1p depth")

    # C: UMAP by cell_area (segmentation size — descriptive)
    axC = axes[1, 0]
    if "cell_area" in obs:
        area = obs["cell_area"].to_numpy().astype(float)
        vlo, vhi = np.percentile(area, [2, 98])
        scC = _umap_scatter(axC, xy, np.clip(area, vlo, vhi), continuous=True,
                            cmap="cividis", title="broad UMAP — cell_area")
        fig.colorbar(scC, ax=axC, fraction=0.04, pad=0.02, label="cell_area (2–98%ile)")
    else:
        axC.axis("off")

    # D: UMAP by celltype_broad (the structure the latent encodes)
    axD = axes[1, 1]
    cb = obs["celltype_broad"].astype(str).to_numpy()
    _umap_scatter(axD, xy, cb, pal=S.pal_broad,
                  title="broad UMAP — celltype_broad (lineage structure)")
    bcats = [c for c in S.broad_present if c in set(cb)]
    outside_legend(axD, cat_handles(bcats, S.pal_broad), bcats,
                   title="broad lineage")
    fig.suptitle(f"{S.sid} — scANVI latent diagnostic "
                 f"(depth/area are descriptive, not gated)", fontsize=20, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    save(fig, figdir / "scvi_latent_diagnostic.png")


def fig_composition(S, figdir):
    obs = S.obs
    broad = obs["celltype_broad"].astype(str)
    fine = obs["celltype"].astype(str)
    b2f = getattr(S.mod, "BROAD_TO_FINE", {}) or {}
    broad_order = [b for b in b2f if b in set(broad)]
    for b in pd.unique(broad):                    # any broad not in b2f
        if b not in broad_order:
            broad_order.append(b)

    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(20, max(7, 0.5 * len(S.fine_present) + 3)),
        gridspec_kw={"width_ratios": [1.15, 1.0]})

    # Left: fine-within-broad stacked horizontal bars
    ct = pd.crosstab(broad, fine)
    ct = ct.reindex(index=broad_order)
    fracs = ct.div(ct.sum(1), axis=0).fillna(0.0)
    ypos = np.arange(len(broad_order))
    for b_i, b in enumerate(broad_order):
        left = 0.0
        # order segments by fine_order for stable coloring
        for ft in [f for f in S.fine_present if f in fracs.columns]:
            w = fracs.loc[b, ft]
            if w <= 0:
                continue
            axL.barh(b_i, w, left=left, color=S.pal_fine[str(ft)],
                     edgecolor="white", linewidth=0.5)
            if w >= 0.06:                          # label only wide-enough segments
                axL.text(left + w / 2, b_i, ft, ha="center", va="center",
                         fontsize=10, color="white", fontweight="bold")
            left += w
    axL.set_yticks(ypos); axL.set_yticklabels(broad_order, fontsize=13)
    axL.set_xlim(0, 1); axL.set_xlabel("fraction of broad lineage")
    axL.set_title("fine-subtype composition within each broad lineage")
    axL.invert_yaxis()
    for sp in ("top", "right"):
        axL.spines[sp].set_visible(False)

    # Right: total abundance per fine celltype (the color key)
    counts = fine.value_counts()
    order = [c for c in S.fine_present if c in counts.index]
    vals = [counts[c] for c in order]
    ypos2 = np.arange(len(order))
    axR.barh(ypos2, vals, color=[S.pal_fine[str(c)] for c in order],
             edgecolor="#444444", linewidth=0.4)
    axR.set_yticks(ypos2); axR.set_yticklabels(order, fontsize=12)
    axR.invert_yaxis()
    axR.set_xscale("log")
    axR.set_xlabel("cells (log scale)")
    axR.set_title("cell-type abundance")
    for i, v in enumerate(vals):
        axR.text(v * 1.05, i, f"{v:,}", va="center", fontsize=10)
    for sp in ("top", "right"):
        axR.spines[sp].set_visible(False)
    fig.suptitle(f"{S.sid} — two-tier structure: broad lineage → fine subtype",
                 fontsize=20, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    save(fig, figdir / "broad_to_fine_composition.png")


def fig_lineage_panels(S, figdir):
    if S.umap_idx is None:
        print("  [fig] lineage_panels skipped (no umap_broad)"); return
    idx = S.umap_idx
    xy = S.umap_xy
    obs = S.obs.iloc[idx]
    broad = obs["celltype_broad"].astype(str).to_numpy()
    fine = obs["celltype"].astype(str).to_numpy()
    b2f = getattr(S.mod, "BROAD_TO_FINE", {}) or {}
    lineages = [b for b in b2f if b in set(broad)]
    for b in pd.unique(broad):
        if b not in lineages:
            lineages.append(b)
    nb = len(lineages)
    ncol = 4 if nb > 6 else 3
    nrow = int(np.ceil(nb / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 5.0 * nrow))
    axes = np.atleast_1d(axes).ravel()
    size = max(2.0, 4e5 / xy.shape[0])
    for k, b in enumerate(lineages):
        ax = axes[k]
        ax.scatter(xy[:, 0], xy[:, 1], s=size * 0.6, c=GREY, linewidths=0,
                   rasterized=True)
        m = broad == b
        sub_fine = fine[m]
        fcats = [c for c in S.fine_present if c in set(sub_fine)]
        for c in fcats:
            mm = m & (fine == c)
            ax.scatter(xy[mm, 0], xy[mm, 1], s=size * 1.6, c=S.pal_fine[str(c)],
                       linewidths=0, rasterized=True, label=c)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{b}  (n={int(m.sum()):,})")
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.legend(handles=cat_handles(fcats, S.pal_fine), labels=fcats,
                  loc="lower left", bbox_to_anchor=(0.0, -0.02), frameon=False,
                  fontsize=10, markerscale=1.4, handletextpad=0.3,
                  labelspacing=0.3, ncol=1)
    for k in range(nb, len(axes)):
        axes[k].axis("off")
    fig.suptitle(f"{S.sid} — Tier-2 fine subtypes within each broad lineage "
                 f"(broad UMAP, grey = other lineages)", fontsize=20, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, figdir / "lineage_umap_panels.png")


def _pick_crop(xy, ct, target=22000):
    """Choose a square crop that sits INSIDE densely-packed tissue (not an empty
    edge). The window side L is sized from the dense-region density (80th
    percentile of occupied-bin density) rather than the global density, so for
    sparse/irregular tissues the window is small enough to fit a packed area.
    Among windows whose cell-type diversity is within 2 of the best, the DENSEST
    (max cell count) is returned. Returns (x0,x1,y0,y1)."""
    x, y = xy[:, 0], xy[:, 1]
    xspan, yspan = x.max() - x.min(), y.max() - y.min()
    H, _, _ = np.histogram2d(x, y, bins=60)
    bin_area = (xspan / 60) * (yspan / 60)
    occ = H[H > 0] / (bin_area + 1e-9)
    dens_hi = float(np.percentile(occ, 80))            # dense-region density
    L = float(np.sqrt(target / dens_hi))
    L = float(np.clip(L, 0.05 * max(xspan, yspan), 0.35 * max(xspan, yspan)))
    steps = 22
    cx = np.linspace(x.min() + L / 2, x.max() - L / 2, steps)
    cy = np.linspace(y.min() + L / 2, y.max() - L / 2, steps)
    cands = []                                         # (ndiv, n, window)
    for ux in cx:
        for uy in cy:
            m = (x >= ux - L / 2) & (x < ux + L / 2) & \
                (y >= uy - L / 2) & (y < uy + L / 2)
            n = int(m.sum())
            if n < target * 0.6:                       # reject mostly-empty
                continue
            ndiv = len(np.unique(ct[m]))
            cands.append((ndiv, n, (ux - L / 2, ux + L / 2,
                                    uy - L / 2, uy + L / 2)))
    if not cands:                                      # fallback: median-centered
        ux, uy = float(np.median(x)), float(np.median(y))
        return (ux - L / 2, ux + L / 2, uy - L / 2, uy + L / 2)
    maxdiv = max(c[0] for c in cands)
    near = [c for c in cands if c[0] >= maxdiv - 2]    # within 2 types of best
    best = max(near, key=lambda c: c[1])               # densest among them
    return best[2]


def fig_spatial_crop(S, figdir, crop=None):
    xy = S.spatial
    fine = S.obs["celltype"].astype(str).to_numpy()
    broad = S.obs["celltype_broad"].astype(str).to_numpy()
    if crop is None:
        x0, x1, y0, y1 = _pick_crop(xy, fine)
    else:
        x0, y0, sz = crop; x1, y1 = x0 + sz, y0 + sz
    m = (xy[:, 0] >= x0) & (xy[:, 0] < x1) & (xy[:, 1] >= y0) & (xy[:, 1] < y1)
    n = int(m.sum())
    cx, cy = xy[m, 0], xy[m, 1]
    size = max(4.0, 9e4 / max(n, 1))

    fig, (axB, axF) = plt.subplots(1, 2, figsize=(22, 11))
    # broad
    bcats = [c for c in S.broad_present if c in set(broad[m])]
    for c in bcats:
        mm = m & (broad == c)
        axB.scatter(xy[mm, 0], xy[mm, 1], s=size, c=S.pal_broad[str(c)],
                    linewidths=0, rasterized=True)
    axB.set_title("celltype_broad")
    outside_legend(axB, cat_handles(bcats, S.pal_broad), bcats,
                   title="broad lineage")
    # fine
    fcats = [c for c in S.fine_present if c in set(fine[m])]
    for c in fcats:
        mm = m & (fine == c)
        axF.scatter(xy[mm, 0], xy[mm, 1], s=size, c=S.pal_fine[str(c)],
                    linewidths=0, rasterized=True)
    axF.set_title("celltype (fine)")
    outside_legend(axF, cat_handles(fcats, S.pal_fine), fcats, title="cell type",
                   ncol=1 if len(fcats) <= 22 else 2)
    for ax in (axB, axF):
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.invert_yaxis()
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.suptitle(f"{S.sid} — spatial cell-type map "
                 f"(crop {x1-x0:.0f}×{y1-y0:.0f} µm, n={n:,} cells)",
                 fontsize=20, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, figdir / "spatial_celltypes_crop.png")


def fig_confidence(S, figdir):
    obs = S.obs
    cts = [c for c in S.fine_present if c in set(obs["celltype"].astype(str))]
    metrics = [("broad_confidence_prob", "scANVI posterior (0–1)"),
               ("subtype_confidence", "Tier-2 score margin")]
    metrics = [(c, lab) for c, lab in metrics if c in obs.columns]
    if not metrics:
        print("  [fig] confidence skipped (no confidence cols)"); return
    fig, axes = plt.subplots(len(metrics), 1, sharex=True,
                             figsize=(max(12, 0.62 * len(cts) + 3),
                                      4.8 * len(metrics)))
    axes = np.atleast_1d(axes)
    ctarr = obs["celltype"].astype(str).to_numpy()
    for ax, (col, lab) in zip(axes, metrics):
        vals = obs[col].to_numpy().astype(float)
        data = [vals[ctarr == c] for c in cts]
        data = [d[np.isfinite(d)] for d in data]
        parts = ax.violinplot([d if len(d) else [0.0] for d in data],
                              positions=np.arange(len(cts)), showmedians=True,
                              widths=0.85)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(S.pal_fine[str(cts[i])]); pc.set_alpha(0.85)
            pc.set_edgecolor("#444444")
        for key in ("cmedians", "cbars", "cmins", "cmaxes"):
            if key in parts:
                parts[key].set_color("#222222"); parts[key].set_linewidth(1.0)
        ax.set_ylabel(lab)
        ax.grid(axis="y", alpha=0.25)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[-1].set_xticks(np.arange(len(cts)))
    axes[-1].set_xticklabels(cts, rotation=90, fontsize=12)
    fig.suptitle(f"{S.sid} — assignment confidence per cell type",
                 fontsize=20, y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    save(fig, figdir / "confidence_summary.png")


def fig_panel_coverage(S, figdir):
    md = S.mod.MARKER_DICT
    cts = [c for c in S.fine_present if c in md]
    present = set(S.var_names)
    ctarr = S.expr_ct
    n_on, dist = [], []                # dist[i] = [%0, %1, %2, %>=3] markers detected
    for c in cts:
        genes = [g for g in md[c] if g in present]
        n_on.append(len(genes))
        m = ctarr == c
        if not genes or not m.any():
            dist.append([1.0, 0.0, 0.0, 0.0]); continue
        d = (S.lognorm(genes)[m] >= LOG_PRESENT).sum(1)
        n = d.size
        dist.append([(d == 0).sum() / n, (d == 1).sum() / n,
                     (d == 2).sum() / n, (d >= 3).sum() / n])
    dist = np.asarray(dist)

    ypos = np.arange(len(cts))
    fig, (axA, axB) = plt.subplots(1, 2, sharey=True,
                                   figsize=(18, max(6, 0.5 * len(cts) + 3.2)))
    # --- A: panel size ---
    axA.barh(ypos, n_on, color=[S.pal_fine[str(c)] for c in cts],
             edgecolor="#444444", linewidth=0.4)
    axA.set_yticks(ypos); axA.set_yticklabels(cts, fontsize=12)
    axA.invert_yaxis()
    axA.set_xlabel("# markers on panel")
    axA.set_title("marker panel size")
    for i, v in enumerate(n_on):
        axA.text(v + 0.1, i, str(v), va="center", fontsize=10)

    # --- B: per-cell detection distribution (stacked) ---
    bin_labels = ["0", "1", "2", "≥3"]
    bin_colors = ["#e0e0e0", "#bdd7e7", "#6baed6", "#08519c"]   # grey→Blues
    left = np.zeros(len(cts))
    for j in range(4):
        axB.barh(ypos, dist[:, j], left=left, color=bin_colors[j],
                 edgecolor="white", linewidth=0.5)
        left += dist[:, j]
    # annotate the % of cells detecting ZERO markers (the cautionary number)
    for i in range(len(cts)):
        z = dist[i, 0]
        if z > 0.005:
            axB.text(min(z / 2, 0.5), i, f"{z*100:.0f}%", va="center",
                     ha="center", fontsize=9,
                     color="#444444" if z < 0.5 else "white")
    axB.set_xlim(0, 1)
    axB.set_xlabel("fraction of the type's cells")
    axB.set_title("per-cell marker detection (lognorm ≥ log1p 1)")
    axB.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=bin_colors[j],
                                      ec="white") for j in range(4)],
               labels=bin_labels, title="# own markers\ndetected / cell",
               loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    for ax in (axA, axB):
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.suptitle(f"{S.sid} — marker evidence per cell type", fontsize=20, y=1.0)
    # honesty caption: detection sparsity is NOT assignment confidence
    fig.text(0.5, -0.02,
             "Cells are assigned by the scANVI broad latent + per-lineage "
             "cluster-mean score argmax — NOT by per-cell marker thresholds. "
             "This panel shows Xenium per-cell detection sparsity (dropout), "
             "not label correctness: a high 0-marker fraction (e.g. small "
             "low-RNA thymocytes, or panel-marginal NK markers) is expected "
             "and does not by itself imply a miscall.",
             ha="center", va="top", fontsize=11, wrap=True, color="#333333")
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    save(fig, figdir / "marker_panel_coverage.png")


def fig_subpanel_coverage(S, figdir):
    """Broad-call coverage of the LARGE curated subpanel pools (var['subpanel_<L>']).

    Sibling to marker_panel_coverage.png, but grouped by the scANVI broad call
    (obs['celltype_broad']) against the 168-991-gene subpanel masks instead of the tight
    fine MARKER_DICT. Left: broad-call x subpanel matrix of mean per-cell detected
    fraction (boxed diagonal = own subpanel; off-diagonal = cross-talk). Right: each
    call's per-cell detection of its OWN subpanel, binned by fraction. Descriptive only;
    enumerates subpanel_<L> columns dynamically so it is tissue-agnostic.
    """
    subpanel_cols = sorted(c for c in S.var.columns if c.startswith("subpanel_"))
    if not subpanel_cols:
        print("  [fig] subpanel skipped (no subpanel_<L> masks in var)"); return
    vnames = np.asarray(S.var_names)
    expr_genes = set(S.expr.var_names)
    genes_by_lin = {}
    for col in subpanel_cols:
        lin = col[len("subpanel_"):]
        genes_by_lin[lin] = [g for g in vnames[S.var[col].to_numpy()] if g in expr_genes]

    broad_sub = S.expr.obs["celltype_broad"].astype(str).to_numpy()
    present_broad = set(np.unique(broad_sub))
    lin_set = set(genes_by_lin)
    b2f = getattr(S.mod, "BROAD_TO_FINE", {}) or {}
    rows = [b for b in b2f if b in present_broad and b in lin_set]
    rows += [b for b in sorted(present_broad & lin_set) if b not in rows]
    if not rows:
        print("  [fig] subpanel skipped (no celltype_broad ∩ subpanel_<L>)"); return
    # order columns so each row's own subpanel lands on the visual diagonal
    cols = list(rows) + [lin for lin in genes_by_lin if lin not in set(rows)]

    nrow, ncol = len(rows), len(cols)
    row_mask = {b: (broad_sub == b) for b in rows}
    M = np.full((nrow, ncol), np.nan)
    own_fc = {}                                  # row lineage -> per-cell own-subpanel frac
    for j, lin in enumerate(cols):
        genes = genes_by_lin[lin]
        if not genes:
            continue
        fcell = (S.lognorm(genes) >= LOG_PRESENT).mean(axis=1)   # one dense block at a time
        for i, b in enumerate(rows):
            m = row_mask[b]
            if m.any():
                M[i, j] = float(fcell[m].mean())
        if lin in row_mask:
            own_fc[lin] = fcell

    e0, e1 = SUBPANEL_FRAC_EDGES
    dist = np.zeros((nrow, 4))
    for i, b in enumerate(rows):
        fc = own_fc.get(b); m = row_mask[b]
        if fc is None or not m.any():
            dist[i] = [1.0, 0.0, 0.0, 0.0]; continue
        v = fc[m]; n = v.size
        dist[i] = [(v == 0).sum() / n, ((v > 0) & (v <= e0)).sum() / n,
                   ((v > e0) & (v <= e1)).sum() / n, (v > e1).sum() / n]

    vmax = float(np.nanmax(M)) if np.isfinite(M).any() else 1.0
    fig, (axM, axB) = plt.subplots(
        1, 2, sharey=True,
        figsize=(max(18, ncol * 1.05 + 7), max(6, 0.62 * nrow + 3.2)),
        gridspec_kw={"width_ratios": [1.3, 1.0]})

    # --- left: broad-call x subpanel mean detected fraction ---
    cmap = plt.get_cmap("Blues").copy(); cmap.set_bad("#f0f0f0")
    im = axM.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
    axM.set_xticks(np.arange(ncol))
    axM.set_xticklabels([f"{lin}\n(n={len(genes_by_lin[lin])})" for lin in cols],
                        rotation=45, ha="right", fontsize=10)
    axM.set_yticks(np.arange(nrow)); axM.set_yticklabels(rows, fontsize=12)
    for i in range(nrow):
        for j in range(ncol):
            if np.isfinite(M[i, j]):
                axM.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=8,
                         color="white" if M[i, j] > 0.55 * vmax else "#333333")
        axM.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False,
                                    edgecolor="#000000", lw=2.0))   # own subpanel (diagonal)
    axM.set_title("broad call × subpanel\n(mean per-cell detected fraction)")
    cb = fig.colorbar(im, ax=axM, fraction=0.030, pad=0.02)
    cb.set_label("mean fraction of subpanel detected", fontsize=11)

    # --- right: own-subpanel per-cell detection (stacked fraction bars) ---
    ypos = np.arange(nrow); left = np.zeros(nrow)
    for k in range(4):
        axB.barh(ypos, dist[:, k], left=left, color=SUBPANEL_BIN_COLORS[k],
                 edgecolor="white", linewidth=0.5)
        left += dist[:, k]
    for i in range(nrow):
        z = dist[i, 0]
        if z > 0.005:
            axB.text(min(z / 2, 0.5), i, f"{z*100:.0f}%", va="center", ha="center",
                     fontsize=9, color="#444444" if z < 0.5 else "white")
    axB.set_xlim(0, 1); axB.set_xlabel("fraction of the lineage's cells")
    axB.set_title("own-subpanel per-cell detection")
    axB.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=SUBPANEL_BIN_COLORS[k], ec="white")
                        for k in range(4)], labels=SUBPANEL_BIN_LABELS,
               title="fraction of own\nsubpanel detected / cell", loc="upper left",
               bbox_to_anchor=(1.02, 1.0), frameon=False)
    axM.invert_yaxis()                           # sharey -> both flip (broad order top->bottom)
    for sp in ("top", "right"):
        axB.spines[sp].set_visible(False)

    fig.suptitle(f"{S.sid} — subpanel coverage by broad lineage "
                 f"(large curated pools, descriptive)", fontsize=20, y=1.0)
    fig.text(0.5, -0.02,
             "Rows are cells grouped by their scANVI broad call; columns are the curated "
             "subpanel gene pools (var['subpanel_<L>'], up to ~1k genes). The boxed "
             "diagonal is each call's OWN subpanel. These pools are large, overlapping and "
             "housekeeping-laden, so high OFF-diagonal fractions are EXPECTED and do not "
             "imply miscalls — this describes per-cell detection, not label correctness. "
             "Rows are reweighted by the per-fine-type 8,000-cell subsample cap.",
             ha="center", va="top", fontsize=11, wrap=True, color="#333333")
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    save(fig, figdir / "subpanel_coverage.png")


def fig_program_coverage(S, figdir):
    """Annotate the OFF-subpanel genes (panel_for_embedding & ~celltype_panel) by the
    vendored pathway/program panels. ANNOTATION ONLY — these panels are never embedding
    or phenotyping input (CLAUDE.md #2). Left: # off-subpanel genes per program (sorted)
    + a coverage summary. Right: program x broad-lineage mean per-cell detection of each
    program's off-subpanel genes. The program loader is imported locally so a missing
    programs/ snapshot never breaks the other figures.
    """
    try:
        from scripts.program_panels import PROGRAMS_DIR, load_program_panels  # noqa: F401
        programs = load_program_panels()
    except FileNotFoundError as exc:
        print(f"  [fig] programs skipped ({exc})"); return
    except Exception as exc:                                  # noqa: BLE001
        print(f"  [fig] programs skipped (load failed: {exc!r})"); return

    if not {"panel_for_embedding", "celltype_panel"} <= set(S.var.columns):
        print("  [fig] programs skipped (no panel_for_embedding/celltype_panel)"); return
    vnames = np.asarray(S.var_names)
    off = set(vnames[S.var["panel_for_embedding"].to_numpy()
                     & ~S.var["celltype_panel"].to_numpy()])
    N = len(off)
    if N == 0:
        print("  [fig] programs skipped (no off-subpanel genes)"); return

    expr_genes = set(S.expr.var_names)
    prog_off = {}
    for p, genes in programs.items():
        gs = [g for g in (genes & off) if g in expr_genes]
        if gs:
            prog_off[p] = gs
    if not prog_off:
        print("  [fig] programs skipped (no off-subpanel gene covered by any program)"); return
    order = sorted(prog_off, key=lambda p: len(prog_off[p]), reverse=True)
    covered = set().union(*(set(prog_off[p]) for p in order))
    M_cov, U = len(covered), N - len(covered)

    broad_sub = S.expr.obs["celltype_broad"].astype(str).to_numpy()
    present_broad = set(np.unique(broad_sub))
    b2f = getattr(S.mod, "BROAD_TO_FINE", {}) or {}
    rows_broad = [b for b in b2f if b in present_broad]
    rows_broad += [b for b in sorted(present_broad) if b not in rows_broad]
    row_mask = {b: (broad_sub == b) for b in rows_broad}

    P, B = len(order), len(rows_broad)
    H = np.full((P, B), np.nan)
    for pi, p in enumerate(order):
        fcell = (S.lognorm(prog_off[p]) >= LOG_PRESENT).mean(axis=1)
        for bj, b in enumerate(rows_broad):
            m = row_mask[b]
            if m.any():
                H[pi, bj] = float(fcell[m].mean())

    fig, (axBar, axH) = plt.subplots(
        1, 2, sharey=True,
        figsize=(max(16, B * 1.1 + 9), max(7, 0.42 * P + 3.0)),
        gridspec_kw={"width_ratios": [1.0, 1.0]})
    ypos = np.arange(P)
    counts = [len(prog_off[p]) for p in order]
    axBar.barh(ypos, counts, color="#6baed6", edgecolor="#444444", linewidth=0.4)
    axBar.set_yticks(ypos); axBar.set_yticklabels(order, fontsize=9)
    axBar.invert_yaxis()
    axBar.set_xlabel("# off-subpanel genes in program")
    axBar.set_title("off-subpanel genes per program")
    for i, v in enumerate(counts):
        axBar.text(v + max(counts) * 0.01, i, str(v), va="center", fontsize=8)
    axBar.text(0.97, 0.02,
               f"{M_cov} of {N} off-subpanel genes\ncovered by ≥1 program; {U} by none",
               transform=axBar.transAxes, ha="right", va="bottom", fontsize=10,
               bbox=dict(boxstyle="round", fc="white", ec="#999999"))
    for sp in ("top", "right"):
        axBar.spines[sp].set_visible(False)

    vmax = float(np.nanmax(H)) if np.isfinite(H).any() else 1.0
    cmap = plt.get_cmap("magma").copy(); cmap.set_bad("#f0f0f0")
    im = axH.imshow(np.ma.masked_invalid(H), aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
    axH.set_xticks(np.arange(B))
    axH.set_xticklabels(rows_broad, rotation=45, ha="right", fontsize=11)
    axH.set_title("mean per-cell detection of a program's\noff-subpanel genes, by broad lineage")
    cb = fig.colorbar(im, ax=axH, fraction=0.030, pad=0.02)
    cb.set_label("mean fraction detected", fontsize=11)

    fig.suptitle(f"{S.sid} — off-subpanel genes annotated by pathway/program "
                 f"(NOT embedding input)", fontsize=20, y=1.0)
    fig.text(0.5, -0.02,
             (f"The {N} genes in panel_for_embedding but in NO cell-type subpanel "
              f"(celltype_panel) are classified here by the {len(programs)} vendored program panels "
              "(data/reference/subpanels/programs/). These panels are ANNOTATION ONLY: "
              "never in panel_for_embedding / celltype_panel / any subpanel_<L> mask, they "
              "drive no obs score and cross cell-type boundaries by design — which is "
              "exactly why they are excluded from the embedding (CLAUDE.md #2). Programs "
              "overlap, so per-program counts are not a partition."),
             ha="center", va="top", fontsize=11, wrap=True, color="#333333")
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    save(fig, figdir / "program_offpanel_coverage.png")


FIGS = {
    "dotplot": fig_dotplot,
    "heatmap": fig_heatmap,
    "latent": fig_latent_diagnostic,
    "composition": fig_composition,
    "lineage": fig_lineage_panels,
    "spatial": fig_spatial_crop,
    "confidence": fig_confidence,
    "coverage": fig_panel_coverage,
    "subpanel": fig_subpanel_coverage,
    "programs": fig_program_coverage,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("samples", nargs="*", default=list(DEFAULT_TISSUE))
    ap.add_argument("--tissue", default=None,
                    help="override tissue (else inferred per sample)")
    ap.add_argument("--only", nargs="+", choices=list(FIGS),
                    help="render only these figures")
    args = ap.parse_args()
    which = args.only or list(FIGS)
    for sid in args.samples:
        tissue = args.tissue or DEFAULT_TISSUE.get(sid)
        if tissue is None:
            print(f"!! unknown tissue for {sid}; pass --tissue"); continue
        print(f"\n=== {sid} ({tissue}) ===")
        figdir = DATA_ROOT / f"figures/{sid}/02_phenotyping"
        S = Sample(sid, tissue)
        for key in which:
            try:
                FIGS[key](S, figdir)
            except Exception as exc:
                import traceback
                print(f"  [fig] {key} FAILED: {exc!r}")
                traceback.print_exc()


if __name__ == "__main__":
    main()
