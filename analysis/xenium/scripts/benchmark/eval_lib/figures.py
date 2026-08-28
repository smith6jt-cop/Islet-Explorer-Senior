"""Figure rendering for the benchmark. Agg backend, defensive (skips a figure
if its inputs are missing). Reuses the repo's shared palette/sizing helpers."""
from __future__ import annotations

import sys
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[3]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[5]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
try:
    from utils.analysis_utils import umap_dot_size
except Exception:  # noqa: BLE001
    def umap_dot_size(n, k=9e6, lo=2.0, hi=12.0):
        return float(np.clip(k / max(n, 1), lo, hi))

DPI = 200


def _save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"    fig -> {path.name}")


def spatial_domain_maps(xy, label_dict: dict, save: Path, title=""):
    methods = list(label_dict)
    k = len(methods)
    ncol = min(k, 3)
    nrow = int(np.ceil(k / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 5 * nrow), squeeze=False)
    axes = axes.ravel()
    s = max(0.5, min(4.0, umap_dot_size(len(xy))))
    for ax, m in zip(axes, methods):
        cats = pd.Categorical(np.asarray(label_dict[m]).astype(str))
        ax.scatter(xy[:, 0], xy[:, 1], c=cats.codes, cmap="tab20", s=s,
                   linewidths=0, rasterized=True)
        ax.set_title(f"{m}  (n={cats.categories.size})", fontsize=11)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    for ax in axes[k:]:
        ax.axis("off")
    fig.suptitle(title or "Spatial domains by method", fontsize=13)
    _save(fig, save)


def domain_quality_bars(metrics_df: pd.DataFrame, cols, save: Path):
    cols = [c for c in cols if c in metrics_df.columns]
    if not cols:
        return
    n = len(cols)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.6), squeeze=False)
    axes = axes[0]
    for ax, c in zip(axes, cols):
        d = metrics_df[c].dropna()
        ax.bar(range(len(d)), d.values, color="#4C72B0")
        ax.set_xticks(range(len(d)))
        ax.set_xticklabels(d.index, rotation=45, ha="right", fontsize=8)
        ax.set_title(c, fontsize=10)
    fig.suptitle("Domain-quality metrics (read with n_domains; see bias notes)", fontsize=11)
    _save(fig, save)


def agreement_heatmap(mat: pd.DataFrame, save: Path, title="cross-method ARI"):
    if mat is None or mat.empty:
        return
    fig, ax = plt.subplots(figsize=(1.1 * len(mat) + 2, 1.0 * len(mat) + 1.5))
    im = ax.imshow(mat.values.astype(float), cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(mat))); ax.set_xticklabels(mat.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index)
    for i in range(len(mat)):
        for j in range(len(mat)):
            v = mat.values[i, j]
            if pd.notna(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if v < 0.6 else "black", fontsize=8)
    ax.set_title(title); fig.colorbar(im, ax=ax, fraction=0.046)
    _save(fig, save)


def stability_bar(stab: dict, save: Path):
    items = [(m, v) for m, v in stab.items() if v.get("n_seeds", 0) >= 2]
    if not items:
        return
    fig, ax = plt.subplots(figsize=(1.0 * len(items) + 2, 3.5))
    names = [m for m, _ in items]
    means = [v["mean_ari"] for _, v in items]
    sds = [v["sd_ari"] for _, v in items]
    ax.bar(range(len(names)), means, yerr=sds, capsize=4, color="#55A868")
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_ylim(0, 1.02); ax.set_ylabel("mean pairwise ARI across seeds")
    ax.set_title("Seed stability (higher = more reproducible)")
    _save(fig, save)


def enrichment_heatmap(pearson: pd.DataFrame, save: Path, method=""):
    if pearson is None or pearson.empty:
        return
    fig, ax = plt.subplots(figsize=(0.45 * pearson.shape[1] + 3,
                                    0.45 * pearson.shape[0] + 2))
    vmax = float(np.nanpercentile(np.abs(pearson.values), 98)) or 1.0
    im = ax.imshow(pearson.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(pearson.shape[1]))
    ax.set_xticklabels(pearson.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(pearson.shape[0]))
    ax.set_yticklabels(pearson.index, fontsize=8)
    ax.set_title(f"{method}: domain x celltype enrichment (Pearson residual)", fontsize=10)
    ax.set_xlabel("celltype"); ax.set_ylabel("domain")
    fig.colorbar(im, ax=ax, fraction=0.046)
    _save(fig, save)


def scalability_curves(metrics_long: pd.DataFrame, save: Path):
    """metrics_long: columns [method, n_cells, wall_s, peak_rss_mb, status]."""
    if metrics_long is None or metrics_long.empty:
        return
    d = metrics_long.dropna(subset=["n_cells", "wall_s"])
    if d.empty:
        return
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    for m, g in d.groupby("method"):
        g = g.sort_values("n_cells")
        a1.plot(g["n_cells"], g["wall_s"], "o-", label=m)
        if "peak_rss_mb" in g:
            a2.plot(g["n_cells"], g["peak_rss_mb"], "o-", label=m)
    for a, ylab, title in ((a1, "wall (s)", "runtime vs cells"),
                           (a2, "peak RSS (MB)", "memory vs cells")):
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlabel("cells"); a.set_ylabel(ylab); a.set_title(title)
        a.legend(fontsize=8); a.grid(True, which="both", alpha=0.3)
    fig.suptitle("Scalability (capped runs excluded)", fontsize=12)
    _save(fig, save)


# ---- MISTy figures ----------------------------------------------------------
def misty_improvement(gain: pd.Series, save: Path):
    if gain is None or len(gain) == 0:
        return
    g = gain.sort_values()
    fig, ax = plt.subplots(figsize=(6, max(3, 0.16 * len(g))))
    ax.barh(range(len(g)), g.values, color="#C44E52")
    ax.set_yticks(range(len(g))); ax.set_yticklabels(g.index, fontsize=6)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("gain.R2 (multiview - intraview)")
    ax.set_title("MISTy per-target multiview gain")
    _save(fig, save)


def misty_contributions(by_view: dict, save: Path):
    if not by_view:
        return
    fig, ax = plt.subplots(figsize=(4, 4))
    views = list(by_view); vals = [by_view[v] for v in views]
    ax.bar(range(len(views)), vals, color=["#4C72B0", "#DD8452", "#55A868"][:len(views)])
    ax.set_xticks(range(len(views))); ax.set_xticklabels(views, rotation=30, ha="right")
    ax.set_ylabel("mean contribution fraction")
    ax.set_title("MISTy view contributions")
    _save(fig, save)


def misty_importance_heatmap(mat: pd.DataFrame, save: Path, view=""):
    if mat is None or mat.empty:
        return
    m = mat.fillna(0.0)
    fig, ax = plt.subplots(figsize=(0.25 * m.shape[1] + 3, 0.25 * m.shape[0] + 3))
    im = ax.imshow(m.values, cmap="magma", aspect="auto")
    ax.set_title(f"MISTy importances ({view})", fontsize=10)
    ax.set_xlabel("target"); ax.set_ylabel("predictor")
    ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046)
    _save(fig, save)


def misty_bridge_scatter(table: pd.DataFrame, rho: float, save: Path):
    if table is None or table.empty:
        return
    fig, ax = plt.subplots(figsize=(4.5, 4.2))
    ax.scatter(table["moran_I"], table["gain_R2"], s=14, color="#8172B3")
    ax.set_xlabel("per-marker Moran's I (spatial autocorrelation)")
    ax.set_ylabel("per-marker MISTy gain.R2")
    ax.set_title(f"Bridge: gain.R2 vs Moran's I  (Spearman rho={rho:.2f})")
    _save(fig, save)
