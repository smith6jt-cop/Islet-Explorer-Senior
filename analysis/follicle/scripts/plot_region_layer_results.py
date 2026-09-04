#!/usr/bin/env python
"""Figures for the cell-derived region layer: basic results first, then method checks.

    python scripts/plot_region_layer_results.py

Writes PNGs to analysis/figures/. Every chart ships direct value labels and has
a CSV counterpart in analysis/tables/ — the palette's aqua slot sits below 3:1
on the light surface, so the relief rule applies.
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from shapely.ops import unary_union

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT / "analysis"))
import data_utils as du  # noqa: E402

# Reference palette, light mode. Slots 1-3 are the documented all-pairs-safe
# subset; RedPulp takes a neutral rather than a fourth hue, which is also what
# it is semantically — the tissue remainder, not a peer category.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8983"
NEUTRAL, GRID = "#d5d4cd", "#e8e7e1"

GENO_COLOR = {"C/C": BLUE, "C/T": ORANGE, "T/T": AQUA}
REGION_COLOR = {"Follicle": BLUE, "PALS": ORANGE,
                "Perifollicular": AQUA, "RedPulp": NEUTRAL}
# Platform is a secondary encoding — shape, not hue — so it never competes with
# genotype for the reader's colour channel.
PLATFORM_MARKER = {"Phenocycler": "o", "CODEX": "s"}
PLATFORM_SHORT = {"Phenocycler": "Pheno", "CODEX": "CODEX"}
REGIONS = ["Follicle", "PALS", "Perifollicular", "RedPulp"]


def style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": NEUTRAL, "axes.linewidth": 0.8,
        "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.titlesize": 11, "axes.labelsize": 10,
        "grid.color": GRID, "grid.linewidth": 0.8,
        "legend.frameon": False, "legend.fontsize": 9,
        "font.size": 10, "figure.dpi": 150,
    })


def despine(ax, keep=("left", "bottom")):
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def box_strip(ax, data, value, title, ylabel):
    """Box + per-donor points, coloured by genotype, shaped by platform."""
    groups = [g for g in du.GENO_ORDER if (data["Genotype"] == g).any()]
    for i, g in enumerate(groups):
        sub = data[data["Genotype"] == g]
        bp = ax.boxplot([sub[value].values], positions=[i], widths=0.55,
                        patch_artist=True, showfliers=False, zorder=1)
        for patch in bp["boxes"]:
            patch.set(facecolor=GENO_COLOR[g], alpha=0.16,
                      edgecolor=GENO_COLOR[g], linewidth=1.6)
        for part in ("whiskers", "caps"):
            for ln in bp[part]:
                ln.set(color=GENO_COLOR[g], linewidth=1.2)
        for ln in bp["medians"]:
            ln.set(color=GENO_COLOR[g], linewidth=2.2)
        rng = np.random.default_rng(0)
        jitter = rng.uniform(-0.13, 0.13, len(sub))
        for (_, r), jx in zip(sub.iterrows(), jitter):
            ax.plot(i + jx, r[value], PLATFORM_MARKER.get(r["Platform"], "o"),
                    markersize=8, markerfacecolor=GENO_COLOR[g],
                    markeredgecolor=SURFACE, markeredgewidth=1.6,
                    alpha=0.95, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{g}\nn={(data['Genotype'] == g).sum()}" for g in groups])
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=8)
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)
    despine(ax)


def main():
    style()
    d = pd.read_csv(du.TABLES_DIR / "P1_vessel_density_cell_regions.csv")
    stats = pd.read_csv(du.TABLES_DIR / "P1_vessel_density_stats.csv")
    qc = pd.read_csv(du.TABLES_DIR / "region_layer_qc.csv")
    d["Genotype"] = pd.Categorical(d["Genotype"], du.GENO_ORDER, ordered=True)

    def p_of(region, test):
        r = stats[(stats.Region == region) & (stats.Test.str.contains(test))]
        return float(r["p"].iloc[0]) if len(r) else np.nan

    # ---- F1 basic: the primary endpoint, one bar per donor ---------------
    fol = d[d.Region == "Follicle"].dropna(subset=["Density_per_mm2"])
    fol = fol.sort_values("Density_per_mm2")
    fig, ax = plt.subplots(figsize=(8.5, 6))
    y = np.arange(len(fol))
    ax.barh(y, fol["Density_per_mm2"], height=0.68,
            color=[GENO_COLOR[g] for g in fol["Genotype"]], zorder=2)
    for yi, (_, r) in zip(y, fol.iterrows()):
        ax.text(r["Density_per_mm2"] + 8, yi, f"{r['Density_per_mm2']:.0f}",
                va="center", ha="left", fontsize=8.5, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r.Sample}  ({PLATFORM_SHORT[r.Platform]})" for _, r in fol.iterrows()],
        fontsize=9)
    ax.set_xlabel("SmallVessel density in Follicle (per mm²)")
    ax.set_title("Per-donor follicle vessel density, cell-derived regions",
                 loc="left", pad=10)
    ax.set_xlim(0, fol["Density_per_mm2"].max() * 1.14)
    ax.grid(axis="x", zorder=0)
    ax.set_axisbelow(True)
    despine(ax)
    ax.legend(handles=[Patch(facecolor=GENO_COLOR[g], label=g) for g in du.GENO_ORDER],
              loc="lower right", title="rs3184504", title_fontsize=9)
    du.save_figure(fig, "P1_F1_follicle_density_per_donor")

    # ---- F2 basic: the primary test --------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5))
    box_strip(ax, fol, "Density_per_mm2",
              "Follicle vessel density by rs3184504", "vessels per mm²")
    kw, rho = p_of("Follicle", "Kruskal"), p_of("Follicle", "Spearman")
    ax.text(0.02, 0.97, f"Kruskal-Wallis p = {kw:.3f}\nSpearman dosage p = {rho:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=9, color=INK2)
    ax.legend(handles=[Line2D([], [], marker=m, linestyle="none", color=MUTED,
                              markersize=8, label=p)
                       for p, m in PLATFORM_MARKER.items()],
              loc="upper right", title="platform", title_fontsize=9)
    du.save_figure(fig, "P1_F2_follicle_density_by_genotype")

    # ---- F3 detail: every region -----------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, region in zip(axes.ravel(), REGIONS):
        sub = d[d.Region == region].dropna(subset=["Density_per_mm2"])
        box_strip(ax, sub, "Density_per_mm2", region, "vessels per mm²")
        ax.text(0.02, 0.97,
                f"KW p = {p_of(region, 'Kruskal'):.3f}   n = {len(sub)}",
                transform=ax.transAxes, va="top", ha="left",
                fontsize=9, color=INK2)
    fig.suptitle("Vessel density by genotype, all cell-derived regions",
                 x=0.02, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    du.save_figure(fig, "P1_F3_density_by_region_genotype", tight=False)

    # ---- F4 detail: dosage trend -----------------------------------------
    dose = {"C/C": 0, "C/T": 1, "T/T": 2}
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
    for ax, region in zip(axes.ravel(), REGIONS):
        sub = d[d.Region == region].dropna(subset=["Density_per_mm2"])
        x = sub["Genotype"].map(dose).astype(float)
        rng = np.random.default_rng(1)
        for (_, r), jx in zip(sub.iterrows(), rng.uniform(-0.11, 0.11, len(sub))):
            ax.plot(dose[r["Genotype"]] + jx, r["Density_per_mm2"],
                    PLATFORM_MARKER.get(r["Platform"], "o"), markersize=8,
                    markerfacecolor=GENO_COLOR[r["Genotype"]],
                    markeredgecolor=SURFACE, markeredgewidth=1.6, zorder=3)
        if len(sub) > 2:
            b, a = np.polyfit(x, sub["Density_per_mm2"], 1)
            xs = np.linspace(-0.25, 2.25, 20)
            ax.plot(xs, a + b * xs, "-", color=MUTED, linewidth=2, zorder=2)
        row = stats[(stats.Region == region) & (stats.Test.str.contains("Spearman"))]
        eff = row["Effect_Size"].iloc[0] if len(row) else ""
        ax.text(0.02, 0.97, f"{eff}   p = {p_of(region, 'Spearman'):.3f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9, color=INK2)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(du.GENO_ORDER)
        ax.set_xlabel("T-allele dosage")
        ax.set_ylabel("vessels per mm²")
        ax.set_title(region, loc="left", pad=8)
        ax.grid(axis="y", zorder=0)
        ax.set_axisbelow(True)
        despine(ax)
    fig.suptitle("Gene-dosage trend, rs3184504", x=0.02, ha="left",
                 fontsize=12, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    du.save_figure(fig, "P1_F4_dosage_trend", tight=False)

    # ---- F5 detail: the platform confound --------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    ax = axes[0]
    for i, plat in enumerate(["Phenocycler", "CODEX"]):
        sub = fol[fol.Platform == plat]
        bp = ax.boxplot([sub["Density_per_mm2"].values], positions=[i], widths=0.5,
                        patch_artist=True, showfliers=False, zorder=1)
        for patch in bp["boxes"]:
            patch.set(facecolor=MUTED, alpha=0.12, edgecolor=MUTED, linewidth=1.4)
        for part in ("whiskers", "caps"):
            for ln in bp[part]:
                ln.set(color=MUTED, linewidth=1.2)
        for ln in bp["medians"]:
            ln.set(color=INK2, linewidth=2.2)
        rng = np.random.default_rng(2)
        for (_, r), jx in zip(sub.iterrows(), rng.uniform(-0.12, 0.12, len(sub))):
            ax.plot(i + jx, r["Density_per_mm2"], PLATFORM_MARKER[plat],
                    markersize=8, markerfacecolor=GENO_COLOR[r["Genotype"]],
                    markeredgecolor=SURFACE, markeredgewidth=1.6, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Phenocycler\nn={(fol.Platform=='Phenocycler').sum()}",
                        f"CODEX\nn={(fol.Platform=='CODEX').sum()}"])
    ax.set_ylabel("vessels per mm²")
    ax.set_title("Follicle density differs by platform", loc="left", pad=8)
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)
    despine(ax)

    ax = axes[1]
    donors = qc.drop_duplicates("Sample")
    ct = pd.crosstab(donors["Platform"], donors["Genotype"]).reindex(
        columns=du.GENO_ORDER, fill_value=0)
    left = np.zeros(len(ct))
    for g in du.GENO_ORDER:
        ax.barh(np.arange(len(ct)), ct[g], left=left, height=0.55,
                color=GENO_COLOR[g], label=g, zorder=2)
        for i, (v, l) in enumerate(zip(ct[g], left)):
            if v:
                ax.text(l + v / 2, i, str(int(v)), ha="center", va="center",
                        fontsize=10, color=SURFACE, fontweight="bold")
        left = left + ct[g].values
    ax.set_yticks(np.arange(len(ct)))
    ax.set_yticklabels(ct.index)
    ax.set_xlabel("donors")
    ax.set_title("Genotype now spans both platforms", loc="left", pad=8)
    ax.legend(title="rs3184504", loc="lower right", title_fontsize=9)
    despine(ax, keep=("left",))
    ax.set_xticks([])
    fig.tight_layout()
    du.save_figure(fig, "P1_F5_platform_check", tight=False)

    # ---- F6 verification: every donor at ONE scale ------------------------
    layer = du.load_region_layer()
    by_sample = {du.extract_sample_id(k): v for k, v in layer.items()}
    order = qc.sort_values(["Platform", "Genotype", "Sample"])["Sample"].tolist()

    # The main view is a fixed window the size of the smallest whole section, so
    # a µm renders to the same number of pixels on every panel and structures
    # are directly comparable. Scaling each panel to its own section instead
    # drew a 226 mm2 donor's follicles 19x smaller than a 12 mm2 donor's. The
    # whole section moves to a small corner locator showing where the window is.
    spans = []
    for s_ in order:
        allp = [q for ps in by_sample.get(s_, {}).values() for q in ps]
        if allp:
            x0, y0, x1, y1 = unary_union(allp).bounds
            spans.append(max(x1 - x0, y1 - y0))
    WIN = min(spans) if spans else 4000.0

    def draw(ax, polys, lw=0.5):
        for region in ["RedPulp", "PALS", "Perifollicular", "Follicle"]:
            for poly in polys.get(region, []):
                xs, ys = poly.exterior.xy
                ax.fill(xs, ys, facecolor=REGION_COLOR[region], edgecolor=SURFACE,
                        linewidth=lw, alpha=0.95 if region != "RedPulp" else 0.5)

    ncol = 4
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.6 * nrow))
    for ax, sample in zip(axes.ravel(), order):
        polys = by_sample.get(sample, {})
        allp = [q for ps in polys.values() for q in ps]
        if not allp:
            ax.set_visible(False)
            continue
        x0, y0, x1, y1 = unary_union(allp).bounds

        # centre the window on the follicle-densest part
        fols = polys.get("Follicle", [])
        if fols:
            cs = np.array([[p.centroid.x, p.centroid.y] for p in fols])
            counts = [((np.abs(cs[:, 0] - c[0]) < WIN / 2) &
                       (np.abs(cs[:, 1] - c[1]) < WIN / 2)).sum() for c in cs]
            cx, cy = cs[int(np.argmax(counts))]
        else:
            cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        cx = np.clip(cx, x0 + WIN / 2, max(x1 - WIN / 2, x0 + WIN / 2))
        cy = np.clip(cy, y0 + WIN / 2, max(y1 - WIN / 2, y0 + WIN / 2))

        draw(ax, polys, lw=0.6)
        ax.set_xlim(cx - WIN / 2, cx + WIN / 2)
        ax.set_ylim(cy + WIN / 2, cy - WIN / 2)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); despine(ax, keep=())

        row = qc[qc.Sample == sample].iloc[0]
        ax.set_title(f"{sample}  ·  {row.Genotype}  ·  {PLATFORM_SHORT[row.Platform]}",
                     loc="left", fontsize=9.5, pad=4)
        nb = row.get("Perifollicular_follicles_with_band", np.nan)
        band = f" · {int(nb)}/{int(row.Follicle_n)} with band" if nb == nb else ""
        ax.text(0.02, 0.02,
                f"{int(row.Follicle_n)} follicles · {row.tissue_mm2:.0f} mm² total{band}",
                transform=ax.transAxes, fontsize=8, color=INK2, va="bottom")
        ax.plot([cx + WIN / 2 - 1150, cx + WIN / 2 - 150],
                [cy + WIN / 2 - 150] * 2, "-", color=INK2, linewidth=2.5,
                solid_capstyle="butt")
        ax.text(cx + WIN / 2 - 650, cy + WIN / 2 - 200, "1 mm", ha="center",
                va="bottom", fontsize=7.5, color=INK2)

        # corner locator: the whole section, with the window marked
        if max(x1 - x0, y1 - y0) > WIN * 1.05:
            axl = ax.inset_axes([0.72, 0.72, 0.28, 0.28])
            draw(axl, polys, lw=0.0)
            axl.add_patch(plt.Rectangle((cx - WIN / 2, cy - WIN / 2), WIN, WIN,
                                        fill=False, edgecolor=INK2, linewidth=1.2))
            axl.set_xlim(x0, x1); axl.set_ylim(y1, y0)
            axl.set_aspect("equal"); axl.set_xticks([]); axl.set_yticks([])
            axl.set_facecolor(SURFACE)
            for sp in axl.spines.values():
                sp.set(edgecolor=MUTED, linewidth=0.8)

    for ax in axes.ravel()[len(order):]:
        ax.set_visible(False)
    fig.legend(handles=[Patch(facecolor=REGION_COLOR[r], label=r) for r in REGIONS],
               loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(f"Regions per donor — every panel is the same {WIN/1000:.1f} mm "
                 "window; corner locator shows where it sits in the section",
                 x=0.02, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    du.save_figure(fig, "P1_F6_region_maps", tight=False)

    # ---- F7 verification: cell-derived vs QuPath -------------------------
    af = pd.read_csv(du.DATA_ROOT / "Measurements" / "AnnotationsFinal.csv",
                     usecols=["Image", "Classification", "Area µm^2"])
    f3 = pd.read_csv(du.DATA_ROOT / "Measurements" / "ForSH2B3.csv",
                     usecols=["Image", "Classification", "Area µm^2"])
    qp = pd.concat([af, f3])
    qp["Sample"] = qp["Image"].apply(du.extract_sample_id)
    qpf = (qp[qp.Classification == "Follicle"].groupby("Sample")["Area µm^2"].sum() / 1e6)
    cmp = qc[["Sample", "Platform", "Genotype", "Follicle_mm2"]].copy()
    cmp["QuPath_mm2"] = cmp["Sample"].map(qpf)
    cmp = cmp.dropna()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    ax = axes[0]
    lim = [min(cmp.QuPath_mm2.min(), cmp.Follicle_mm2.min()) * 0.6,
           max(cmp.QuPath_mm2.max(), cmp.Follicle_mm2.max()) * 1.6]
    ax.plot(lim, lim, "-", color=NEUTRAL, linewidth=2, zorder=1)
    ax.text(lim[1] * 0.62, lim[1] * 0.72, "equal", fontsize=8.5,
            color=MUTED, rotation=45)
    for _, r in cmp.iterrows():
        ax.plot(r.QuPath_mm2, r.Follicle_mm2, PLATFORM_MARKER[r.Platform],
                markersize=9, markerfacecolor=GENO_COLOR[r.Genotype],
                markeredgecolor=SURFACE, markeredgewidth=1.6, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("QuPath follicle area (mm²)")
    ax.set_ylabel("cell-derived follicle area (mm²)")
    rho = cmp[["Follicle_mm2", "QuPath_mm2"]].corr(method="spearman").iloc[0, 1]
    ax.set_title(f"Same ranking, larger scale  (Spearman ρ = {rho:.3f})",
                 loc="left", pad=8)
    ax.grid(zorder=0); ax.set_axisbelow(True); despine(ax)

    ax = axes[1]
    cmp2 = cmp.assign(ratio=cmp.Follicle_mm2 / cmp.QuPath_mm2).sort_values("ratio")
    yy = np.arange(len(cmp2))
    ax.barh(yy, cmp2.ratio, height=0.66,
            color=[GENO_COLOR[g] for g in cmp2.Genotype], zorder=2)
    ax.axvline(1, color=MUTED, linewidth=1.5, zorder=3)
    for yi, v in zip(yy, cmp2.ratio):
        ax.text(v + 0.08, yi, f"{v:.2f}×", va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(yy)
    ax.set_yticklabels(cmp2.Sample, fontsize=8.5)
    ax.set_xlabel("cell-derived ÷ QuPath follicle area")
    ax.set_title("Inflation varies by donor, not by platform", loc="left", pad=8)
    ax.set_xlim(0, cmp2.ratio.max() * 1.2)
    ax.grid(axis="x", zorder=0); ax.set_axisbelow(True); despine(ax)
    fig.tight_layout()
    du.save_figure(fig, "P1_F7_cell_vs_qupath_follicle", tight=False)

    # ---- F8 verification: region composition -----------------------------
    comp = qc.set_index("Sample")[[f"{r}_mm2" for r in REGIONS]]
    comp.columns = REGIONS
    frac = comp.div(comp.sum(axis=1), axis=0).loc[order]
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    left = np.zeros(len(frac))
    for region in REGIONS:
        ax.barh(np.arange(len(frac)), frac[region], left=left, height=0.62,
                color=REGION_COLOR[region], label=region,
                edgecolor=SURFACE, linewidth=1.4, zorder=2)
        left = left + frac[region].values
    for i, s in enumerate(frac.index):
        ax.text(frac.loc[s, "Follicle"] / 2, i, f"{100*frac.loc[s,'Follicle']:.0f}%",
                ha="center", va="center", fontsize=8, color=SURFACE,
                fontweight="bold")
    ax.set_yticks(np.arange(len(frac)))
    ax.set_yticklabels(frac.index, fontsize=9)
    ax.set_xlabel("fraction of tissue")
    ax.set_xlim(0, 1)
    ax.set_title("Region composition per donor (Follicle % labelled)",
                 loc="left", pad=10)
    ax.legend(loc="lower right", ncol=2)
    despine(ax)
    fig.tight_layout()
    du.save_figure(fig, "P1_F8_region_composition", tight=False)

    # ---- F9 detail: is the inflation size-dependent? ---------------------
    from scipy.stats import spearmanr, mannwhitneyu
    cmp2 = cmp.assign(ratio=cmp.Follicle_mm2 / cmp.QuPath_mm2)
    rs, ps = spearmanr(cmp2.QuPath_mm2, cmp2.ratio)
    fig, ax = plt.subplots(figsize=(6.5, 5))
    for _, r in cmp2.iterrows():
        ax.plot(r.QuPath_mm2, r.ratio, PLATFORM_MARKER[r.Platform], markersize=9,
                markerfacecolor=GENO_COLOR[r.Genotype], markeredgecolor=SURFACE,
                markeredgewidth=1.6, zorder=3)
        ax.annotate(r.Sample, (r.QuPath_mm2, r.ratio), textcoords="offset points",
                    xytext=(8, 3), fontsize=7.5, color=MUTED)
    ax.axhline(1, color=NEUTRAL, linewidth=2, zorder=1)
    ax.set_xscale("log")
    ax.set_xlabel("QuPath follicle area (mm²)")
    ax.set_ylabel("cell-derived ÷ QuPath")
    ax.set_title(f"Inflation shrinks as follicles get larger\n"
                 f"Spearman ρ = {rs:.2f}, p = {ps:.3f}", loc="left", pad=8)
    ax.grid(zorder=0); ax.set_axisbelow(True); despine(ax)
    fig.tight_layout()
    du.save_figure(fig, "P1_F9_inflation_vs_size", tight=False)

    # ---- F10 detail: is the perifollicular band platform-dependent? ------
    q = qc.copy()
    q["PFZ_frac"] = q["Perifollicular_mm2"] / q["tissue_mm2"]
    q["mean_follicle_mm2"] = q["Follicle_mm2"] / q["Follicle_n"].replace(0, np.nan)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, (col, lab) in zip(axes, [("PFZ_frac", "Perifollicular ÷ tissue"),
                                     ("mean_follicle_mm2", "mean follicle area (mm²)")]):
        for i, plat in enumerate(["Phenocycler", "CODEX"]):
            sub = q[q.Platform == plat]
            bp = ax.boxplot([sub[col].dropna().values], positions=[i], widths=0.5,
                            patch_artist=True, showfliers=False, zorder=1)
            for patch in bp["boxes"]:
                patch.set(facecolor=MUTED, alpha=0.12, edgecolor=MUTED, linewidth=1.4)
            for part in ("whiskers", "caps"):
                for ln in bp[part]:
                    ln.set(color=MUTED, linewidth=1.2)
            for ln in bp["medians"]:
                ln.set(color=INK2, linewidth=2.2)
            rng = np.random.default_rng(3)
            for (_, r), jx in zip(sub.iterrows(), rng.uniform(-0.12, 0.12, len(sub))):
                ax.plot(i + jx, r[col], PLATFORM_MARKER[plat], markersize=8,
                        markerfacecolor=GENO_COLOR[r.Genotype],
                        markeredgecolor=SURFACE, markeredgewidth=1.6, zorder=3)
        a = q[q.Platform == "Phenocycler"][col].dropna()
        b = q[q.Platform == "CODEX"][col].dropna()
        p = mannwhitneyu(a, b)[1]
        ax.text(0.02, 0.97, f"Mann-Whitney p = {p:.3f}", transform=ax.transAxes,
                va="top", fontsize=9, color=INK2)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Phenocycler\nn={len(a)}", f"CODEX\nn={len(b)}"])
        ax.set_ylabel(lab)
        ax.set_title(lab, loc="left", pad=8)
        ax.grid(axis="y", zorder=0); ax.set_axisbelow(True); despine(ax)
    fig.suptitle("Does the perifollicular band really differ by platform?",
                 x=0.02, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    du.save_figure(fig, "P1_F10_perifollicular_platform_check", tight=False)
    q[["Sample", "Platform", "Genotype", "PFZ_frac", "mean_follicle_mm2"]].to_csv(
        du.TABLES_DIR / "P1_perifollicular_platform_check.csv", index=False)

    # ---- F11: radial marker composition around follicles -----------------
    prof_path = du.TABLES_DIR / "radial_profiles.csv"
    if prof_path.exists():
        prof = pd.read_csv(prof_path)
        mcols = [c for c in prof.columns
                 if c not in ("Sample", "Platform", "Genotype", "distance_um", "n_cells")]
        med = prof.groupby("distance_um")[mcols].median()
        med = med.loc[(med.index >= -240) & (med.index <= 450)]
        show = ["CD20", "CD21", "SMA", "HLA-DR", "CD45RO", "Ki67"]
        show = [m for m in show if m in med.columns]
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

        ax = axes[0]
        for i, m in enumerate(show):
            ax.plot(med.index, med[m], "-", linewidth=2,
                    color=[BLUE, AQUA, ORANGE, MUTED, NEUTRAL, INK2][i % 6], label=m)
            ax.annotate(m, (med.index[0], med[m].iloc[0]), fontsize=8.5,
                        color=INK2, xytext=(3, 0), textcoords="offset points")
        ax.axvline(0, color=INK2, linewidth=1.4)
        ax.axhline(1, color=NEUTRAL, linewidth=1.4)
        ax.set_yscale("log")
        ax.set_xlabel("distance from follicle boundary (µm) — negative = inside")
        ax.set_ylabel("fold over red-pulp baseline")
        ax.set_title("Marker composition in 30 µm radial bands", loc="left", pad=8)
        ax.grid(zorder=0); ax.set_axisbelow(True); despine(ax)

        ax = axes[1]
        sm = med["SMA"] if "SMA" in med.columns else None
        if sm is not None:
            o = sm.loc[sm.index > -60]
            ax.plot(o.index, o.values, "-", linewidth=2.5, color=ORANGE)
            ax.axvline(0, color=INK2, linewidth=1.4)
            ax.axhline(1, color=NEUTRAL, linewidth=1.4)
            wpath = du.TABLES_DIR / "radial_pfz_width.csv"
            if wpath.exists():
                w = pd.read_csv(wpath)
                wmed = w.loc[w.PFZ_width_um <= 200, "PFZ_width_um"].median()
                ax.axvspan(0, wmed, color=AQUA, alpha=0.18, zorder=0)
                ax.text(wmed + 20, ax.get_ylim()[1] * 0.96,
                        f"perifollicular\n0–{wmed:.0f} µm", ha="left", va="top",
                        fontsize=9, color=INK2)
            ax.set_xlabel("distance from follicle boundary (µm)")
            ax.set_ylabel("SMA, fold over red pulp")
            ax.set_title("SMA peaks just outside the follicle — the band it defines",
                         loc="left", pad=8)
            ax.grid(zorder=0); ax.set_axisbelow(True); despine(ax)
        fig.tight_layout()
        du.save_figure(fig, "P2_F11_radial_profiles", tight=False)

    print("\nfigures written to", du.FIGURES_DIR)


if __name__ == "__main__":
    main()
