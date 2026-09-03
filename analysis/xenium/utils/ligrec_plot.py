"""Shared, tissue-agnostic rendering for squidpy ligand-receptor (ligrec)
dotplots at the LINEAGE level.

Both the stage-03 notebook (`notebooks/03_spatial_analysis.ipynb`, cell 14) and
the ad-hoc re-render script (`scripts/regen_ligrec_lineage.py`) import from here,
so the two colormap/legend fixes below run identically in every code path and
for every tissue.

The two fixes (originally only in regen_ligrec_lineage.py, commit 59a0f15):

1. Colormap scale — squidpy's CustomDotplot bakes vmin = nanmin(dot_color_df)
   into both scanpy's color norm and its own colorbar ticks. Because squidpy
   fills filtered-out cells with 0 (line 319 of squidpy/pl/_ligrec.py), the
   colorbar always starts at 0 even when the smallest plotted log2(mean+1) is
   ~0.38. That leaves the displayed dots clustered in the upper half of viridis,
   with a few outliers stretching the top. `_RescaledDotplot` injects vmin/vmax
   (upper end capped at the 95th percentile of plotted means) into both the norm
   and the colorbar tick range, spreading the colormap over the range where the
   data actually lives.
2. Legend spacing — squidpy spreads the three right-side legend axes (size dots,
   significance true/false, colorbar) vertically across the full figure height,
   so they overlap. `compact_ligrec_legends` pulls them into a tight stack and
   centers each axis's dot row.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import squidpy as sq
import squidpy.pl._ligrec as _sq_ligrec_mod

_OrigCustomDotplot = _sq_ligrec_mod.CustomDotplot


class _RescaledDotplot(_OrigCustomDotplot):
    """Forward injected vmin/vmax into scanpy DotPlot, and use them for the
    colorbar tick range too (squidpy's parent uses nanmin/nanmax of the
    dot_color_df, which is contaminated by 0-filled empty cells)."""

    _injected_vmin: float | None = None
    _injected_vmax: float | None = None

    def __init__(self, *args, **kwargs):
        if _RescaledDotplot._injected_vmin is not None:
            kwargs["vmin"] = _RescaledDotplot._injected_vmin
        if _RescaledDotplot._injected_vmax is not None:
            kwargs["vmax"] = _RescaledDotplot._injected_vmax
        super().__init__(*args, **kwargs)

    def _plot_colorbar(self, color_legend_ax, normalize) -> None:
        from matplotlib.colorbar import ColorbarBase

        cmap = plt.colormaps[self.cmap]
        if _RescaledDotplot._injected_vmin is not None:
            vmin = _RescaledDotplot._injected_vmin
        else:
            vmin = float(np.nanmin(self.dot_color_df.values))
        if _RescaledDotplot._injected_vmax is not None:
            vmax = _RescaledDotplot._injected_vmax
        else:
            vmax = float(np.nanmax(self.dot_color_df.values))
        ColorbarBase(
            color_legend_ax,
            orientation="horizontal",
            cmap=cmap,
            norm=normalize,
            ticks=np.linspace(vmin, vmax, self.DEFAULT_NUM_COLORBAR_TICKS),
            format="%.2f",
        )
        color_legend_ax.set_title(self.color_legend_title, fontsize="small")
        color_legend_ax.xaxis.set_tick_params(labelsize="small")


# Install the monkeypatch on import so any sq.pl.ligrec call below uses it.
_sq_ligrec_mod.CustomDotplot = _RescaledDotplot


def compact_ligrec_legends(fig) -> None:
    """squidpy ligrec spreads the three right-side legend axes (size dots,
    significance true/false, colorbar) vertically across the full figure
    height. Pull them into a tight stack just above the colorbar so the
    legend reads as one block.

    Why this is non-trivial: squidpy places the size-legend dots at data
    y=1 with ylim ≈ (-1.2, 4), while the significance dots sit at data
    y=0 with auto-fit ylim ≈ (-0.05, 0.05). The two axes scale their
    data coords completely differently, so naively applying the same
    ylim to both leaves size dots near the top and significance dots
    near the bottom of their boxes (the user-reported overlap). Instead,
    detect each axis's actual dot y-position and choose a per-axis ylim
    that centers the dots, so titles and labels line up across both
    sub-legends."""
    cands = []
    for ax in fig.axes:
        bb = ax.get_position()
        # Right-margin small axes only (excludes main heatmap + cluster
        # bracket). Width ~0.10 fig-units, height < 0.30.
        if bb.x0 > 0.75 and bb.width < 0.20 and bb.height < 0.30:
            cands.append((ax, bb))
    if len(cands) < 2:
        return
    short = [c for c in cands if c[1].height < 0.10]   # colorbar
    tall = [c for c in cands if c[1].height >= 0.10]   # size + significance
    if not (short and tall):
        return
    sb = short[0][1]
    new_h = 0.05
    target_y = sb.y0 + sb.height + 0.14
    for ax, bb in tall:
        dot_ys = []
        for coll in ax.collections:
            offs = coll.get_offsets()
            if offs is not None and len(offs) > 0:
                dot_ys.extend(np.asarray(offs)[:, 1].tolist())
        ymid = float(np.median(dot_ys)) if dot_ys else 0.5
        ax.set_position([bb.x0, target_y, bb.width, new_h])
        # Symmetric ylim around the dot row → dots always at axes-y=0.5.
        # Half-width 1.0 leaves ~half the box above for the title and
        # ~half below for the xtick labels.
        ax.set_ylim(ymid - 1.0, ymid + 1.0)
        if ax.get_title():
            # axes-y=1.40 sits the title ~40% above the axes top — high
            # enough that a 2-line title (e.g. "significant\np=0.0001")
            # doesn't collide with the dot row, matched across both
            # sub-legends so the visible titles align horizontally.
            ax.title.set_position((0.5, 1.40))
            ax.title.set_size("small")


def render_ligrec_per_source(
    sample_id: str,
    lr: dict,
    lineages: list[str],
    fig_dir,
    *,
    top_n: int = 30,
    min_mean: float = 0.3,
    pval_thresh: float = 0.001,
    alpha: float = 1e-4,
    p_floor: float = 1e-12,
    min_pairs: int = 2,
    dpi: int = 130,
) -> list:
    """Render one dotplot per source lineage from a squidpy ligrec result dict
    (``{'means', 'pvalues', 'metadata'}``). Applies the colormap-scale and
    legend-spacing fixes. Returns the list of written PNG paths.

    Tissue-agnostic: `lineages` is whatever set of source categories the caller
    passes (7 pancreas / 8 thymus lineages, etc.)."""
    from pathlib import Path

    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    means = lr["means"]
    pvals = lr["pvalues"]
    metadata = lr.get("metadata")

    pvals_clean = np.where(
        np.isfinite(pvals.to_numpy(dtype=np.float64)),
        np.maximum(pvals.to_numpy(dtype=np.float64), p_floor),
        1.0,
    )
    means_arr = np.where(
        np.isfinite(means.to_numpy(dtype=np.float64)),
        means.to_numpy(dtype=np.float64),
        0.0,
    )

    written = []
    for src in lineages:
        col_mask = np.array([c1 == src for c1, c2 in means.columns])
        if not col_mask.any():
            continue
        m_src = means.loc[:, col_mask]
        p_src = pvals.loc[:, col_mask]
        means_src = means_arr[:, col_mask]
        pvals_src = pvals_clean[:, col_mask]

        sig_any = ((pvals_src < pval_thresh) & (means_src >= min_mean)).any(axis=1)
        # A dotplot with a single significant pair (n_lr == 1) makes squidpy's
        # legend collections degenerate — matplotlib then raises "need at least
        # one array to concatenate" when computing extents at draw time. Such a
        # plot carries no comparative information anyway, so skip sources below
        # min_pairs. (Pancreas never hit this; several thymus lineages do.)
        if sig_any.sum() < min_pairs:
            print(f"[{sample_id}] {src}: {int(sig_any.sum())} significant "
                  f"L-R pair(s) < min_pairs={min_pairs} — skipping")
            continue
        m_src = m_src.loc[sig_any]
        p_src = p_src.loc[sig_any]
        score = (-np.log10(pvals_src[sig_any])).max(axis=1)
        order = np.argsort(-score)[:top_n]
        m_top = m_src.iloc[order].copy()
        p_top = p_src.iloc[order].copy()

        # Color-scale targeting: in raw form a handful of outliers stretch
        # the colorbar all the way to the data max while the bulk of the
        # plotted dots sit in a narrow band, so most cells end up dark.
        # Cap the upper end at the 95th percentile of cells that pass both
        # the mean and p-value filters; everything above clips to the same
        # bright color, but the mid-range now gets the full viridis spread.
        m_top_vals = m_top.to_numpy(dtype=np.float64)
        p_top_vals = p_top.to_numpy(dtype=np.float64)
        plotted_mask = np.isfinite(m_top_vals) & (m_top_vals >= min_mean) & \
            np.isfinite(p_top_vals) & (p_top_vals < pval_thresh)
        plotted_vals = m_top_vals[plotted_mask]
        if plotted_vals.size:
            upper_cap = float(np.quantile(plotted_vals, 0.95))
            # Squidpy will log2(m+1)-transform inside, so the colorbar shows
            # log2(min_mean+1) on the left and log2(upper_cap+1) on the right.
            inj_vmin = float(np.log2(min_mean + 1.0))
            inj_vmax = float(np.log2(upper_cap + 1.0))
            # Clip outliers so they don't pull the colormap norm.
            m_top_plot = m_top.clip(upper=upper_cap)
        else:
            inj_vmin = None
            inj_vmax = None
            m_top_plot = m_top

        res = {"means": m_top_plot, "pvalues": p_top}
        if metadata is not None:
            res["metadata"] = metadata.reindex(m_top.index).copy()

        # With swap_axes=True the data is transposed for display:
        # x-axis has the L-R pairs (n_lr columns -> width), y-axis has the
        # cluster-pairs (n_cp rows -> height). Width scales with n_lr,
        # height scales with n_cp + small fixed pad for labels.
        # No subplots_adjust — squidpy manages the gridspec internally; an
        # external adjust pushed the suptitle off-canvas.
        n_lr = m_top.shape[0]
        n_cp = m_top.shape[1]
        fig_w = max(10.0, 0.26 * n_lr + 4.0)
        fig_h = max(4.5, 0.40 * n_cp + 3.0)

        _RescaledDotplot._injected_vmin = inj_vmin
        _RescaledDotplot._injected_vmax = inj_vmax
        try:
            sq.pl.ligrec(
                res,
                means_range=(min_mean, np.inf),
                pvalue_threshold=pval_thresh,
                alpha=alpha,
                remove_empty_interactions=True,
                swap_axes=True,
                figsize=(fig_w, fig_h),
                dpi=dpi,
            )
        finally:
            _RescaledDotplot._injected_vmin = None
            _RescaledDotplot._injected_vmax = None
        # squidpy's `title` arg renders inconsistently across versions —
        # set the suptitle explicitly after the layout is built. Then
        # compact the right-side legends so they sit close together.
        fig = plt.gcf()
        fig.suptitle(
            f"{sample_id} — {src} → all targets  (top {n_lr} L-R interactions)",
            fontsize=12, y=1.0,
        )
        compact_ligrec_legends(fig)
        out = fig_dir / f"ligand_receptor_dotplot_{src}.png"
        # Some sources end up with a degenerate dotplot after
        # remove_empty_interactions strips zero rows (e.g. n_cp shrinks to 1).
        # matplotlib then raises "need at least one array to concatenate" in
        # PathCollection.get_extents at draw time. Skip those rather than
        # killing the whole sample.
        try:
            plt.savefig(out, dpi=dpi, bbox_inches="tight")
        except ValueError as exc:
            print(f"[{sample_id}] {src}: SKIP (degenerate dotplot — "
                  f"{type(exc).__name__}: {exc})")
            plt.close("all")
            continue
        plt.close("all")
        size_mb = out.stat().st_size / 1e6
        print(f"[{sample_id}] {src}: wrote {out.name}  "
              f"({n_lr} pairs × {n_cp} targets, {size_mb:.2f} MB)")
        written.append(out)

    # Remove the legacy monolithic file if present so users don't open the
    # superseded subtype-level plot. Keep the per-source PNGs as the new index.
    legacy = fig_dir / "ligand_receptor_dotplot.png"
    if legacy.exists():
        legacy.unlink()
        print(f"[{sample_id}] removed legacy {legacy.name}")
    return written
