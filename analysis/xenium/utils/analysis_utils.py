"""
Utility functions for Xenium spatial transcriptomics analysis
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def load_xenium_data(
    data_path: Path,
    sample_name: str,
    file_type: str = "h5ad"
) -> sc.AnnData:
    """
    Load Xenium data from various file formats.

    Parameters
    ----------
    data_path : Path
        Path to data directory
    sample_name : str
        Sample identifier
    file_type : str
        File format ('h5ad', 'h5', 'csv')

    Returns
    -------
    AnnData object
    """
    file_path = data_path / f"{sample_name}.{file_type}"

    if file_type == "h5ad":
        adata = sc.read_h5ad(file_path)
    elif file_type == "h5":
        adata = sc.read_10x_h5(file_path)
    elif file_type == "csv":
        adata = sc.read_csv(file_path)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    return adata


def add_spatial_coordinates(
    adata: sc.AnnData,
    coords_file: Path | None = None,
    x_col: str = "x",
    y_col: str = "y"
) -> sc.AnnData:
    """
    Add spatial coordinates to AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    coords_file : Path, optional
        Path to CSV file with coordinates
    x_col : str
        Column name for x coordinates
    y_col : str
        Column name for y coordinates

    Returns
    -------
    AnnData with spatial coordinates in obsm['spatial']
    """
    if coords_file is not None:
        coords = pd.read_csv(coords_file, index_col=0)
        spatial_coords = coords[[x_col, y_col]].values
    else:
        # Try to get from obs
        if x_col in adata.obs and y_col in adata.obs:
            spatial_coords = adata.obs[[x_col, y_col]].values
        else:
            raise ValueError("Spatial coordinates not found")

    adata.obsm['spatial'] = spatial_coords
    return adata


def calculate_qc_metrics(
    adata: sc.AnnData,
    mt_prefix: str = "MT-",
    ribo_prefix: tuple[str, str] = ("RPS", "RPL")
) -> sc.AnnData:
    """
    Calculate comprehensive QC metrics.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    mt_prefix : str
        Prefix for mitochondrial genes
    ribo_prefix : tuple
        Prefixes for ribosomal genes

    Returns
    -------
    AnnData with QC metrics in obs
    """
    # Identify gene types
    adata.var['mt'] = adata.var_names.str.startswith(mt_prefix)
    adata.var['ribo'] = adata.var_names.str.startswith(ribo_prefix)

    # Calculate metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    return adata


def filter_cells_and_genes(
    adata: sc.AnnData,
    min_genes: int = 10,
    min_cells: int = 3,
    max_mt_pct: float = 20.0,
    min_counts: int = 50
) -> sc.AnnData:
    """
    Filter low-quality cells and genes.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    min_genes : int
        Minimum genes per cell
    min_cells : int
        Minimum cells per gene
    max_mt_pct : float
        Maximum mitochondrial percentage
    min_counts : int
        Minimum UMI counts per cell

    Returns
    -------
    Filtered AnnData
    """
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filter
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if 'pct_counts_mt' in adata.obs:
        adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :].copy()

    print(f"Filtered: {n_cells_before - adata.n_obs} cells, "
          f"{n_genes_before - adata.n_vars} genes")

    return adata


def normalize_and_hvg(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    flavor: str = 'seurat_v3'
) -> sc.AnnData:
    """
    Normalize data and identify highly variable genes.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    target_sum : float
        Target sum for normalization
    n_top_genes : int
        Number of highly variable genes
    flavor : str
        Method for HVG selection

    Returns
    -------
    Normalized AnnData with HVG annotation
    """
    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers['log1p_norm'] = adata.X.copy()

    # Identify HVGs
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        subset=False
    )

    return adata


def plot_highly_variable_genes(
    var: pd.DataFrame,
    save_path: Path,
    max_labels: int = 30,
    rel_floor: float = 1.5,
    dpi: int = 200,
    figsize: tuple[float, float] = (16, 7),
    show: bool = False,
) -> Path:
    """
    Two-panel highly-variable-gene plot with adaptive, non-overlapping labels.

    Renders mean-counts vs variance (left) and mean-counts vs normalized
    variance (right), both log-log, with seurat_v3 "highly variable" genes in
    red over a grey background of the rest. The gene names that get labelled
    are chosen *relative to the data* rather than by a fixed cutoff:

    * a gene is a label *candidate* only if its normalized variance is at least
      ``rel_floor`` x the sample median normalized variance — a data-relative
      gate that keeps bulk genes unlabelled when a sample has a flat tail;
    * at most ``max_labels`` candidates are labelled, taken in descending
      normalized variance.

    The cap (``max_labels``) is what guarantees legibility; the median-multiple
    floor is the relative outlier gate. Labels are de-overlapped with
    ``adjustText`` and joined to their point by a thin leader line, so names
    never collide no matter how tightly the top genes cluster.

    Parameters
    ----------
    var : pd.DataFrame
        An ``adata.var``-style frame with columns ``highly_variable``,
        ``means``, ``variances`` and (optionally) ``variances_norm``; the index
        holds gene names. Pass ``adata.var`` directly.
    save_path : Path
        Destination PNG path (parent dirs are created).
    max_labels : int
        Maximum gene names annotated per panel (the legibility cap).
    rel_floor : float
        Label candidates need normalized variance >= ``rel_floor`` *
        median(normalized variance). Relative, sample-adaptive outlier gate.
    dpi : int
        Figure resolution.
    figsize : tuple
        Figure size in inches (the right panel is dropped, width halved, if no
        ``variances_norm`` column is present).
    show : bool
        If True call ``plt.show()`` (notebook use); otherwise close the figure
        after saving (batch use).

    Returns
    -------
    Path
        ``save_path``.
    """
    import matplotlib.pyplot as plt
    from adjustText import adjust_text

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    is_hv = var["highly_variable"].astype(bool).to_numpy()
    means = var["means"].to_numpy()
    var_raw = var["variances"].to_numpy()
    var_norm = (
        var["variances_norm"].to_numpy()
        if "variances_norm" in var.columns else None
    )
    names = var.index.to_numpy()

    # --- adaptive, data-relative label selection (shared by both panels) ---
    # Rank by normalized variance (the meaningful HVG metric) when available,
    # else fall back to raw variance for the single-panel case.
    rank_metric = var_norm if var_norm is not None else var_raw
    floor = rel_floor * float(np.median(rank_metric))
    candidates = np.flatnonzero(rank_metric >= floor)
    label_idx = candidates[np.argsort(-rank_metric[candidates])][:max_labels]

    def _scatter(ax, y, ylabel, title):
        ax.scatter(means[~is_hv], y[~is_hv], s=8, c="#cfcfcf", alpha=0.55,
                   linewidths=0, label=f"Other ({(~is_hv).sum():,})")
        ax.scatter(means[is_hv], y[is_hv], s=14, c="#d62728", alpha=0.85,
                   linewidths=0, label=f"Highly variable ({is_hv.sum():,})")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("mean counts")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc="best", framealpha=0.9)
        ax.grid(True, alpha=0.3, linestyle="--")

    def _label(ax, y, long_lines=False):
        texts = [
            ax.text(means[i], y[i], str(names[i]), fontsize=7, color="#222")
            for i in label_idx
        ]
        if not texts:
            return
        # The selected genes pile into a tight corner of the diagonal on the
        # raw-variance ("Not normalized") panel. There we push labels much
        # harder off their anchors and let them travel further (force_static /
        # force_explode / max_move), so the leader lines are LONGER and the
        # labels fan out into the empty regions. The normalized panel already
        # has vertical room, so it keeps shorter lines. min_arrow_len=0 draws a
        # thin leader line for every label, even ones that barely move.
        if long_lines:
            kw = dict(expand=(1.5, 2.0), force_static=(1.2, 1.6),
                      force_explode=(2.0, 3.0), max_move=(300, 300))
        else:
            kw = dict(expand=(1.3, 1.6), force_static=(0.3, 0.5),
                      force_explode=(0.5, 1.0), max_move=(60, 60))
        adjust_text(
            texts, ax=ax,
            force_text=(0.5, 0.8),
            min_arrow_len=0,
            arrowprops=dict(arrowstyle="-", color="#888", lw=0.5),
            **kw,
        )

    # Build all panels first, lay out, THEN place labels — adjustText must run
    # against the final axes geometry (tight_layout would otherwise shift it).
    # long_lines=True for the raw-variance panel (cramped diagonal) only.
    if var_norm is None:
        fig, ax = plt.subplots(figsize=(figsize[0] / 2, figsize[1]))
        _scatter(ax, var_raw, "variance", "Highly variable genes")
        panels = [(ax, var_raw, True)]
    else:
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        _scatter(axes[0], var_raw, "variance", "Not normalized")
        _scatter(axes[1], var_norm, "normalized variance", "Normalized")
        panels = [(axes[0], var_raw, True), (axes[1], var_norm, False)]
    fig.tight_layout()
    for ax, y, long_lines in panels:
        _label(ax, y, long_lines=long_lines)

    # Honest framing: HVG is a QC/reference diagnostic, NOT a feature filter.
    # PCA and scVI run on the curated `panel_for_embedding` gene set (the panel
    # minus audited tissue-restricted + sex-chromosome genes), so the red "highly
    # variable" points are informational only — they are not the genes that drive
    # the embedding or the cell-type calls. (Mirrors the caption on the PCA scree
    # plot; keep the two consistent.)
    fig.text(
        0.5, -0.02,
        "Highly variable genes (seurat_v3) are computed for QC/reference only and "
        "do NOT subset features.\nPCA & scVI use the curated panel_for_embedding "
        "gene set, not this highly-variable subset.",
        transform=fig.transFigure, ha="center", va="top", fontsize=8.5,
        color="#555555",
        bbox=dict(boxstyle="round", fc="white", ec="#cccccc", alpha=0.9),
    )

    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return save_path


def plot_scvi_latent_diagnostic(
    adata,
    save_path: Path,
    scvi_key: str = "X_scvi",
    umap_key: str | None = None,
    counts_col: str = "total_counts",
    max_cells: int = 200_000,
    leak_thresh: float = 0.3,
    dpi: int = 200,
    figsize: tuple[float, float] = (14, 6),
    show: bool = False,
    random_state: int = 0,
) -> Path:
    """
    Diagnostic for the scVI latent — the representation phenotyping ACTUALLY uses.

    This replaces the old gene-space PCA scree (``plot_pca_variance_ratio``),
    which was misleading on two counts: (1) phenotyping never touches ``X_pca`` —
    cell typing runs on the scVI latent (``X_scvi`` → cosine kNN → ``leiden_scvi``
    + ``X_umap``) plus ``score_genes`` on the lognorm layer; and (2) the
    gene-space scree ran on z-scored data over ~5k genes, so the total-variance
    denominator was dominated by the unit-variance noise floor of the (mostly
    undetected) panel genes. In low-detection samples (e.g. thymus 0041534:
    median 90 counts/cell, 93% of genes detected in <5% of cells) that produced a
    near-flat spectrum (PC1 ~0.55%, cum-at-30 ~5%) that looked alarming but was a
    scaling artifact, not low intrinsic dimensionality.

    A PCA on the 30-dim scVI latent has no such noise floor: the cumulative
    explained variance reaches ~1.0 within the 30 components, and the
    per-component decay reflects the *real* effective dimensionality of the
    embedding that drives clustering. Two panels:

    * **Panel A** — scree of PCA on ``X_scvi``: per-component explained-variance
      ratio (left, log) + cumulative (right, linear, →1.0), with the
      effective-dimensionality marker (components to reach 90%).
    * **Panel B** — per-latent-dimension std (how each of the 30 scVI dims is
      used), plus the library-size leak gate: Spearman(UMAP axis, total_counts)
      should have ``|rho| < leak_thresh`` (the same gate computed inline in
      stage-02 cell 17).

    Parameters
    ----------
    adata : AnnData
        Phenotyped object; must carry ``obsm[scvi_key]``. ``obsm[umap_key]`` and
        ``obs[counts_col]`` enable Panel B's leak gate (skipped if absent).
    save_path : Path
        Destination PNG path (parent dirs are created).
    scvi_key : str
        ``obsm`` key for the scVI latent (default ``"X_scvi"``).
    umap_key : str | None
        ``obsm`` key for the UMAP used in the leak gate. ``None`` auto-selects
        ``"X_umap_scvi"`` if present, else ``"X_umap"``.
    counts_col : str
        ``obs`` column with per-cell total counts for the leak gate.
    max_cells : int
        Subsample cap for the PCA fit + Spearman (reproducible; matches the
        pipeline's 200k-subsample convention).
    leak_thresh : float
        ``|rho|`` threshold above which library size is flagged as driving the
        embedding.
    dpi, figsize, show : see ``plot_highly_variable_genes``.
    random_state : int
        Seed for the subsample RNG.

    Returns
    -------
    Path
        ``save_path``.
    """
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    if scvi_key not in adata.obsm:
        raise KeyError(
            f"obsm[{scvi_key!r}] not found; cannot render scVI-latent diagnostic"
        )

    if umap_key is None:
        umap_key = "X_umap_scvi" if "X_umap_scvi" in adata.obsm else "X_umap"

    n = adata.n_obs
    if n > max_cells:
        sel = np.random.default_rng(random_state).choice(n, max_cells, replace=False)
        sel.sort()
    else:
        sel = slice(None)

    Z = np.asarray(adata.obsm[scvi_key])[sel]
    n_dims = Z.shape[1]

    # Panel A: PCA on the latent (center only — scVI dims are already comparable).
    pca = PCA(n_components=n_dims, random_state=random_state)
    pca.fit(Z)
    vr = pca.explained_variance_ratio_
    cum = np.cumsum(vr)
    x = np.arange(1, n_dims + 1)
    eff_dim = int(np.searchsorted(cum, 0.90)) + 1  # components to reach 90%

    # Per-dimension std (Panel B): on the (subsampled) latent.
    dim_std = Z.std(axis=0)

    # Library-size leak gate (Panel B annotation): reuse the stage-02 cell-17 logic.
    leak_txt, leak_ok = None, None
    if umap_key in adata.obsm and counts_col in adata.obs:
        from scipy.stats import spearmanr

        um = np.asarray(adata.obsm[umap_key])[sel]
        counts = np.asarray(adata.obs[counts_col])[sel]
        r1 = spearmanr(um[:, 0], counts).correlation
        r2 = spearmanr(um[:, 1], counts).correlation
        rho = max(abs(r1), abs(r2))
        leak_ok = rho < leak_thresh
        leak_txt = (
            f"library-size leak gate (Spearman UMAP vs {counts_col}):\n"
            f"  axis-0 rho = {r1:+.3f}    axis-1 rho = {r2:+.3f}\n"
            f"  max|rho| = {rho:.3f}  "
            + (f"PASS (< {leak_thresh:.1f})" if leak_ok
               else f"WARN (>= {leak_thresh:.1f})")
        )

    fig, (axA, axB) = plt.subplots(1, 2, figsize=figsize)

    # --- Panel A: latent scree -------------------------------------------------
    axA.plot(x, vr, "o-", color="#1f77b4", lw=1.5, markersize=4,
             label="variance ratio (per latent component)")
    axA.set_xlabel("scVI latent component (PCA-ordered)")
    axA.set_ylabel("explained variance ratio (log)", color="#1f77b4")
    axA.set_yscale("log")
    axA.tick_params(axis="y", labelcolor="#1f77b4")
    axA.grid(True, alpha=0.3, linestyle="--")

    axA2 = axA.twinx()
    axA2.plot(x, cum, "-", color="#d62728", lw=2, label="cumulative variance")
    axA2.set_ylabel("cumulative variance (linear, →1.0)", color="#d62728")
    axA2.tick_params(axis="y", labelcolor="#d62728")
    axA2.set_ylim(0, 1.02)
    axA2.axhline(0.90, color="#666", linestyle=":", lw=0.8, alpha=0.6)

    axA.axvline(eff_dim, color="#2ca02c", linestyle="-", lw=2.2, alpha=0.9,
                label=f"90% variance @ {eff_dim} of {n_dims} dims")

    _lA, _llA = axA.get_legend_handles_labels()
    _lA2, _llA2 = axA2.get_legend_handles_labels()
    axA.legend(_lA + _lA2, _llA + _llA2, loc="center right",
               framealpha=0.9, fontsize=8)
    axA.set_title(f"scVI latent scree — {n_dims} dims "
                  f"(cumulative reaches {cum[-1] * 100:.0f}%)")

    # --- Panel B: per-dimension usage + leak gate ------------------------------
    order = np.argsort(dim_std)[::-1]
    axB.bar(np.arange(1, n_dims + 1), dim_std[order], color="#7f7f7f")
    axB.set_xlabel("scVI latent dimension (rank-ordered by std)")
    axB.set_ylabel("per-dimension std")
    axB.set_title("scVI latent dimension usage")
    axB.grid(True, axis="y", alpha=0.3, linestyle="--")

    if leak_txt is not None:
        axB.text(0.97, 0.95, leak_txt, transform=axB.transAxes,
                 fontsize=8, va="top", ha="right", family="monospace",
                 bbox=dict(boxstyle="round",
                           fc=("#e6ffe6" if leak_ok else "#ffe6e6"),
                           ec=("#2ca02c" if leak_ok else "#d62728"), alpha=0.95))
    else:
        axB.text(0.97, 0.95, "no UMAP / counts available\n(leak gate skipped)",
                 transform=axB.transAxes, fontsize=8, va="top", ha="right",
                 color="#888", family="monospace")

    fig.suptitle("scVI-latent diagnostic — the representation phenotyping uses "
                 "(cell types come from X_scvi, not X_pca)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return save_path


# -----------------------------------------------------------------------------
# Shared UMAP rendering helpers (used by BOTH notebooks 01/02 and the
# scripts/regen_*/rerender_* figure-refresh scripts, so a fresh notebook run and
# a post-hoc re-render produce identical figures — no notebook-vs-script drift).
# -----------------------------------------------------------------------------

def umap_dot_size(n_cells: int, k: float = 9_000_000.0,
                  lo: float = 10.0, hi: float = 20.0) -> float:
    """Adaptive scatter point size so UMAP panels with different cell counts
    look equally dense.

    Point ink ~ n_cells * area, so for constant visual density area ~ 1/n.
    We use ``size = k / n_cells`` (matplotlib ``s`` is ~area in points²) clamped
    to ``[lo, hi]``: full-population panels (e.g. 1.2M-cell pancreas) land near
    ``lo`` so dense cores don't saturate, while lineage-subset / rare-type
    panels (far fewer cells) land near ``hi`` so they stay visible. Tune
    ``k``/``lo``/``hi`` to taste — they are the single knobs for every UMAP in
    the pipeline.
    """
    n = max(int(n_cells), 1)
    return float(min(hi, max(lo, k / n)))


# 12-color maximally-distinct qualitative palette for celltype-subset UMAPs.
_HIGH_CONTRAST = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231",
    "#911eb4", "#000000", "#42d4f4", "#f032e6",
    "#bfef45", "#fabed4", "#469990", "#9a6324",
]


def high_contrast_palette(types) -> dict:
    """Return ``{type: hex}`` cycling the 12-color high-contrast palette."""
    return {t: _HIGH_CONTRAST[i % len(_HIGH_CONTRAST)]
            for i, t in enumerate(types)}


def lineage_subsets_for(tissue):
    """Resolve a tissue's ``LINEAGE_SUBSETS`` into
    ``[{'title': str, 'types': [celltype, ...]}, ...]``.

    ``tissue`` may be the already-imported ``scripts.tissue_markers.{x}`` module
    or a tissue-name string. Returns ``None`` if the module defines no
    ``LINEAGE_SUBSETS`` (callers then fall back to a single-panel layout).
    """
    import importlib
    mod = (tissue if hasattr(tissue, "LINEAGE_MAP")
           else importlib.import_module(f"scripts.tissue_markers.{tissue}"))
    if not hasattr(mod, "LINEAGE_SUBSETS"):
        return None
    lineage_map = dict(mod.LINEAGE_MAP)
    out = []
    for sub in mod.LINEAGE_SUBSETS:
        wanted = set(sub["lineages"])
        types = sorted(t for t, lin in lineage_map.items() if lin in wanted)
        out.append({"title": sub["title"], "types": types})
    return out


def render_celltype_subset_panels(adata, axes, subsets, color: str = "celltype",
                                  dot_size=None, legend_fontsize: int = 9):
    """Fill ``axes`` with one lineage-subset UMAP each: the subset's cell types
    in a high-contrast palette over a light-grey background of all other cells
    (scanpy ``groups=``). Splitting a dense ~15-25-category celltype UMAP into a
    few lineage panels keeps individual types distinguishable.

    Point size defaults to ``umap_dot_size(adata.n_obs)`` so these panels match
    the figure's adaptive sizing (all cells are drawn, so size is by the full N,
    not the subset — otherwise the grey background would blob).
    """
    import scanpy as sc
    if dot_size is None:
        dot_size = umap_dot_size(adata.n_obs)
    col = adata.obs[color]
    categories = (list(col.cat.categories) if hasattr(col, "cat")
                  else sorted(set(col.astype(str).unique())))
    present = set(categories)
    for ax, sub in zip(axes, subsets):
        types = [t for t in sub["types"] if t in present]
        if not types:
            ax.set_axis_off()
            continue
        sub_pal = high_contrast_palette(types)
        # scanpy needs a palette entry for every category, not just `groups=`.
        full_pal = {c: sub_pal.get(c, "#cccccc") for c in categories}
        n_in = int(col.astype(str).isin(types).sum())
        sc.pl.umap(
            adata, color=color, groups=types, palette=full_pal,
            ax=ax, show=False, size=dot_size, frameon=False,
            legend_loc="right margin", legend_fontsize=legend_fontsize,
            title=f"{sub['title']} (n={n_in:,})",
        )


def render_major_lineage_panel(adata, ax, lineage_map, palette,
                               color: str = "celltype", dot_size=None,
                               legend_fontsize: int = 10, title=None):
    """Color one UMAP panel by MAJOR lineage instead of fine cell type.

    The fine ``adata.obs[color]`` labels (~15-25 types) are collapsed through
    ``lineage_map`` (e.g. ``tissue_markers.LINEAGE_MAP`` — Endocrine / Exocrine /
    Immune / Stromal / Vascular / Neural / Cycling for pancreas) into ~7 groups,
    drawn with the lineage-level entries of ``palette``
    (``tissue_markers.PALETTE``). This de-crowds the dense fine-type
    ``predicted_celltypes`` / ``scvi_vs_markers`` UMAPs so the major regions
    (acinar bulk vs islets vs immune/stroma) read as distinct, rather than a
    confetti of 20 colors that hides the embedding structure.

    Adds a transient ``_major_lineage`` obs column, so pass a *non-view*
    AnnData (a ``.copy()`` or a freshly read object — assigning obs columns to a
    view does not persist). Point size defaults to ``umap_dot_size(adata.n_obs)``
    to match the rest of the pipeline's adaptive sizing.

    Returns the ordered list of lineage categories actually drawn.
    """
    import pandas as pd
    import scanpy as sc
    if dot_size is None:
        dot_size = umap_dot_size(adata.n_obs)
    lin = adata.obs[color].astype(str).map(dict(lineage_map)).fillna("Other")
    # Order categories by frequency for a stable, readable legend.
    order = list(pd.Series(lin).value_counts().index)
    adata.obs["_major_lineage"] = pd.Categorical(lin, categories=order)
    pal = {c: dict(palette).get(c, "#cccccc") for c in order}
    sc.pl.umap(
        adata, color="_major_lineage", ax=ax, show=False, size=dot_size,
        palette=pal, frameon=False, legend_loc="right margin",
        legend_fontsize=legend_fontsize,
        title=title or f"Major lineage ({len(order)} groups)",
    )
    return order


def run_standard_workflow(
    adata: sc.AnnData,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    resolution: float = 0.5
) -> sc.AnnData:
    """
    Run standard scanpy workflow: scaling, PCA, neighbors, UMAP, clustering.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object
    n_pcs : int
        Number of principal components
    n_neighbors : int
        Number of neighbors for graph
    resolution : float
        Resolution for Leiden clustering

    Returns
    -------
    AnnData with embeddings and clusters
    """
    # Scale
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)

    # Neighbors and UMAP
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(30, n_pcs))
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata, resolution=resolution)

    return adata


def export_to_csv(
    adata: sc.AnnData,
    output_dir: Path,
    sample_name: str,
    include_spatial: bool = True
) -> None:
    """
    Export AnnData to CSV files.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    output_dir : Path
        Output directory
    sample_name : str
        Sample identifier
    include_spatial : bool
        Whether to include spatial coordinates
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Export metadata
    adata.obs.to_csv(output_dir / f"{sample_name}_metadata.csv")

    # Export expression matrix
    expr_df = pd.DataFrame(
        adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    expr_df.to_csv(output_dir / f"{sample_name}_expression.csv")

    # Export spatial coordinates if available
    if include_spatial and 'spatial' in adata.obsm:
        spatial_df = pd.DataFrame(
            adata.obsm['spatial'],
            index=adata.obs_names,
            columns=['x', 'y']
        )
        spatial_df.to_csv(output_dir / f"{sample_name}_spatial_coords.csv")

    print(f"Exported data to {output_dir}")


def create_summary_report(
    adata: sc.AnnData,
    sample_name: str,
    output_file: Path
) -> pd.DataFrame:
    """
    Create a summary report of the dataset.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    sample_name : str
        Sample identifier
    output_file : Path
        Output file path

    Returns
    -------
    DataFrame with summary statistics
    """
    summary = {
        'Sample': sample_name,
        'Total cells': adata.n_obs,
        'Total genes': adata.n_vars,
        'Mean counts per cell': adata.obs['total_counts'].mean() if 'total_counts' in adata.obs else np.nan,
        'Median counts per cell': adata.obs['total_counts'].median() if 'total_counts' in adata.obs else np.nan,
        'Mean genes per cell': adata.obs['n_genes_by_counts'].mean() if 'n_genes_by_counts' in adata.obs else np.nan,
        'Has spatial coords': 'spatial' in adata.obsm,
        'Has cell types': 'celltype' in adata.obs,
    }

    summary_df = pd.DataFrame([summary]).T
    summary_df.columns = ['Value']

    if output_file:
        summary_df.to_csv(output_file)

    return summary_df
