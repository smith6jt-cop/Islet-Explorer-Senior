"""Figures, source tables and the figure-first deck for the main cohort package.

Every figure is drawn only from tables the earlier stages published, and the exact values it plots are
written beside it. Unavailable donors are shown as unavailable, never as zero. The deck (PPTX, then a
PDF rendered from it) is assembled in an isolated presentation environment.
"""
from __future__ import annotations

import json
import os
import subprocess
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from .common import Context, binding, read_json, require, write_json

PRESENTATION_PYTHON = Path("/data/io60-analysis/env/presentation/bin/python")
COLORS = {"ND": "#3a6e94", "AAb+": "#bd8a29", "T1D": "#a74151"}
GROUPS = ("ND", "AAb+", "T1D")
SUBTYPE_COLORS = {"Beta": "#3f7e92", "Alpha": "#d39b2a", "Delta": "#7a5a96"}
LABEL_COLORS = {"Beta": "#3f7e92", "Alpha": "#d39b2a", "Delta": "#7a5a96", "Immune": "#c0392b",
                "Endocrine (subtype unresolved)": "#9aa5ad", "": "#d9d9d9"}
SLIDE = (13.33, 7.5)
plt.rcParams.update({"font.size": 15, "axes.titlesize": 18, "axes.labelsize": 16, "xtick.labelsize": 14,
                     "ytick.labelsize": 14, "legend.fontsize": 13, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none", "pdf.fonttype": 42})


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "render", "Unknown figures operation")
    here = Path(__file__).parent
    return [binding(here / "deck.py"), binding(PRESENTATION_PYTHON.resolve()),
            *(binding(context.task(t) / "completion.json") for t in ("h1", "h2", "h3", "h4", "h5", "h6",
                                                                       "crossstudy", "cd8"))]


class Package:
    def __init__(self, stage: Path):
        self.stage = stage
        for folder in ("figures", "source_data"):
            (stage / folder).mkdir()
        self.slides = []

    def save(self, fig, name: str, title: str, caption: str, notes: str, sources: dict[str, pd.DataFrame]):
        for extension in ("svg", "pdf", "png"):
            fig.savefig(self.stage / "figures" / f"{name}.{extension}", dpi=200 if extension == "png" else None)
        plt.close(fig)
        written = []
        for label, frame in sources.items():
            path = self.stage / "source_data" / f"{name}__{label}.csv"
            frame.to_csv(path, index=False)
            written.append(str(path.relative_to(self.stage)))
        self.slides.append({"name": name, "title": title, "image": f"figures/{name}.png", "caption": caption,
                            "notes": notes, "source_data": written})

    def add_rendered(self, name: str, folder: Path, stem: str, title: str, caption: str, notes: str, tables: list[Path]):
        """A figure another stage already drew: copied as-is, with that stage's tables beside it."""
        import shutil
        for extension in ("svg", "pdf", "png"):
            shutil.copyfile(folder / f"{stem}.{extension}", self.stage / "figures" / f"{name}.{extension}")
        written = []
        for table in tables:
            target = self.stage / "source_data" / f"{name}__{table.name}"
            shutil.copyfile(table, target)
            written.append(str(target.relative_to(self.stage)))
        self.slides.append({"name": name, "title": title, "image": f"figures/{name}.png", "caption": caption,
                            "notes": notes, "source_data": written})


def strip(ax, table: pd.DataFrame, column: str, ylabel: str, label_donors=True):
    """Donor dots by disease group; unavailable donors are listed under their group, not plotted."""
    for x, group in enumerate(GROUPS):
        z = table.loc[table.disease_status.eq(group)].sort_values("donor_id")
        shown = z.loc[np.isfinite(z[column].astype(float))]
        offsets = np.linspace(-0.22, 0.22, len(shown)) if len(shown) > 1 else np.zeros(len(shown))
        ax.scatter(x + offsets, shown[column], s=70, color=COLORS[group], edgecolor="white", linewidth=0.8, zorder=3)
        if len(shown):
            ax.hlines(shown[column].mean(), x - 0.3, x + 0.3, color=COLORS[group], linewidth=2.5, zorder=2)
        if label_donors:
            for dx, (_, row) in zip(offsets, shown.iterrows()):
                ax.annotate(row.donor_id, (x + dx, row[column]), xytext=(4, 2), textcoords="offset points", fontsize=9,
                            color="#555555")
        ax.text(x, -0.13, f"{len(shown)} of {len(z)}", transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=12, color="#555555")
    ax.set_xticks(range(len(GROUPS)), GROUPS)
    ax.set_xlim(-0.6, len(GROUPS) - 0.4)
    ax.set_ylabel(ylabel)


def unavailable_text(table: pd.DataFrame, column="donor_effect") -> str:
    missing = table.loc[~np.isfinite(table[column].astype(float))]
    if missing.empty:
        return ""
    reasons = missing.groupby("reason").donor_id.apply(lambda s: ", ".join(sorted(s)))
    return "\n".join(f"Unavailable ({ids}): {reason[:110]}" for reason, ids in reasons.items())


def outline(ax, geom, color="#222222", width=1.6):
    for part in (geom.geoms if hasattr(geom, "geoms") else [geom]):
        if part.geom_type == "Polygon":
            ax.plot(*part.exterior.xy, color=color, linewidth=width)
            for ring in part.interiors:
                ax.plot(*ring.xy, color=color, linewidth=width)
        elif hasattr(part, "geoms"):
            outline(ax, part, color, width)


def effect_panel(ax, test: dict, units: str):
    if test["status"] != "available":
        ax.text(0.5, 0.5, f"{test['contrast']}\n\nUNAVAILABLE\n" + textwrap.fill(str(test["reason"]), 45),
                transform=ax.transAxes, ha="center", va="center", fontsize=14)
        ax.axis("off")
        return
    ax.axvline(0, color="#888888", linestyle="--", linewidth=1)
    if test["status"] == "available":
        ax.errorbar(test["effect"], 0, xerr=[[test["effect"] - test["CI95_low"]], [test["CI95_high"] - test["effect"]]],
                    fmt="o", color="#222222", markersize=10, capsize=8, linewidth=2)
        ax.set_title(f"{test['contrast']}\np = {test['p']:.3g}   BH q = {test['q_BH6']:.3g}", fontsize=15)
        width = max(abs(test["CI95_low"]), abs(test["CI95_high"]), abs(test["effect"])) * 1.3 or 1
        ax.set_xlim(-width, width)
    else:
        ax.text(0.5, 0.5, "UNAVAILABLE\n" + str(test["reason"]), transform=ax.transAxes, ha="center", va="center",
                fontsize=14, wrap=True)
        ax.set_title(test["contrast"], fontsize=15)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.set_xlabel(units)
    ax.spines["left"].set_visible(False)


def pick_example(endpoints: pd.DataFrame, examples: Path, column="n_islets_available"):
    """The T1D (else ND) donor with an example whose support is the group median; never chosen by effect."""
    for group in ("T1D", "ND"):
        z = endpoints.loc[endpoints.disease_status.eq(group) & endpoints.status.eq("available")]
        z = z.loc[[(examples / f"donor_id={d}.json").is_file() for d in z.donor_id]]
        if len(z):
            median = z[column].median()
            row = z.assign(gap=(z[column] - median).abs()).sort_values(["gap", "donor_id"]).iloc[0]
            return row, read_json(examples / f"donor_id={row.donor_id}.json")
    return None, None


# ---------------------------------------------------------------- atlas

def atlas_figures(context: Context, package: Package):
    atlas = context.task("atlas")
    coverage = pd.read_csv(atlas / "cohort_coverage.csv", dtype={"donor_id": str})
    coverage["group_order"] = coverage.disease_status.map({g: i for i, g in enumerate(GROUPS)})
    coverage = coverage.sort_values(["group_order", "donor_id"]).reset_index(drop=True)
    columns = [("protein_pancreas_cells", "Protein\ncells (k)", 1e-3), ("protein_true_islets", "True\nIslets", 1),
               ("protein_immune_cells", "Protein immune\ncells (k)", 1e-3), ("rna_pancreas_cells", "RNA pancreas\ncells (k)", 1e-3),
               ("rna_pLN_sections", "RNA pLN\nsections", 1), ("postIF_islet_masks", "Post-IF\nislet masks", 1)]
    fig, ax = plt.subplots(figsize=SLIDE)
    fig.subplots_adjust(left=0.02, right=0.99, top=0.9, bottom=0.08)
    for i, row in coverage.iterrows():
        ax.text(-0.6, i, row.donor_id, ha="right", va="center", fontsize=12, color=COLORS[row.disease_status], weight="bold")
        for j, (field, _, factor) in enumerate(columns):
            value = row[field]
            ok = bool(value) and value > 0
            text = f"{value * factor:,.0f}" if ok else ("none" if field != "postIF_islet_masks" else "no mask")
            ax.add_patch(plt.Rectangle((j - 0.45, i - 0.42), 0.9, 0.84, color="#e8f0f4" if ok else "#f3e3e3"))
            ax.text(j, i, text, ha="center", va="center", fontsize=11)
    ax.set_xlim(-1.4, len(columns) - 0.5)
    ax.set_ylim(len(coverage) - 0.5, -1.2)
    for j, (_, label, _) in enumerate(columns):
        ax.text(j, -0.75, label, ha="center", va="bottom", fontsize=13, weight="bold")
    ax.axis("off")
    handles = [Line2D([], [], marker="s", linestyle="", color=COLORS[g], markersize=12, label=g) for g in GROUPS]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, -0.06), ncol=3, frameon=False)
    package.save(fig, "01_cohort_coverage", "Cohort and assay coverage",
                 f"{len(coverage)} donors (ND {int(coverage.disease_status.eq('ND').sum())}, AAb+ "
                 f"{int(coverage.disease_status.eq('AAb+').sum())}, T1D {int(coverage.disease_status.eq('T1D').sum())}). "
                 "Each donor is one biological unit; sections and modalities are views of it, not extra donors.",
                 "Protein: admitted, nucleated cells of the 2026-10-18 typing generation, where epithelium is excluded "
                 "by E-cadherin, EpCAM and pan-cytokeratin on every donor (Ker8_18 is not on the four same-slide panels, "
                 "so it is a veto nowhere). True Islets: original QuPath Islet outlines. RNA: approved release QC "
                 "universe. Post-IF masks: exact-serial registered islet masks used for RNA localization only.",
                 {"coverage": coverage.drop(columns="group_order")})

    summary = pd.read_csv(atlas / "donor_true_islet_summary.csv", dtype={"donor_id": str})
    summary["group_order"] = summary.disease_status.map({g: i for i, g in enumerate(GROUPS)})
    summary = summary.sort_values(["group_order", "donor_id"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=SLIDE)
    bottom = np.zeros(len(summary))
    for subtype, color in SUBTYPE_COLORS.items():
        values = summary[f"{subtype.lower()}_fraction_pooled"].to_numpy(float)
        ax.bar(range(len(summary)), values, bottom=bottom, color=color, label=subtype, width=0.8)
        bottom += np.nan_to_num(values)
    ax.set_xticks(range(len(summary)), summary.donor_id, rotation=60, fontsize=12)
    for tick, group in zip(ax.get_xticklabels(), summary.disease_status):
        tick.set_color(COLORS[group])
    for i, n in enumerate(summary.true_islets):
        ax.text(i, 1.01, str(n), ha="center", va="bottom", fontsize=9, color="#555555", rotation=90)
    ax.set_ylabel("Fraction of subtyped endocrine cells\nin true Islets (pooled per donor)")
    ax.set_ylim(0, 1.12)
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False)
    package.save(fig, "02_islet_endocrine_subtypes", "Endocrine subtypes inside true Islets",
                 "Beta, Alpha and Delta fractions of refined, subtyped endocrine cells inside each donor's original true "
                 "Islets (numbers above bars: true Islets). Donor labels are coloured by disease group.",
                 "Identity is hierarchical: the accepted broad endocrine label kept by refinement, then the subtype. "
                 "Pooled within donor over all true-Islet members; per-islet means are in the source table.",
                 {"donor_summary": summary.drop(columns="group_order")})

    fig, axes = plt.subplots(1, 3, figsize=SLIDE)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, bottom=0.2, wspace=0.42)
    strip(axes[0], summary, "immune_per_mm2_mean_of_islets", "Immune cells per mm² of islet\n(mean of per-islet densities)")
    strip(axes[1], summary, "immune_per_mm2_pooled", "Immune cells per mm² of islet\n(pooled: total / total area)")
    strip(axes[2], summary, "fraction_of_islets_with_immune", "Fraction of true Islets with\n≥1 refined immune cell")
    package.save(fig, "03_islet_immune_density", "Immune cells inside true Islets",
                 "Refined immune cells within original true-Islet outlines, per donor. The mean of per-islet densities and "
                 "the pooled density are different estimands; both are shown. Bars mark group means; n = available of total donors.",
                 "Density divides by the exact outline area of each original true Islet (no expansion). Descriptive only; "
                 "the disease tests are the six spatial hypotheses.",
                 {"donor_summary": summary.drop(columns="group_order")})

    immune_columns = [c for c in summary.columns if c.startswith("Imm-") and c.endswith("_per_mm2_pooled")]
    fig, ax = plt.subplots(figsize=SLIDE)
    bottom = np.zeros(len(summary))
    palette = plt.get_cmap("tab10")
    for k, column in enumerate(immune_columns):
        values = summary[column].to_numpy(float)
        ax.bar(range(len(summary)), np.nan_to_num(values), bottom=bottom, color=palette(k), width=0.8,
               label=column.replace("_per_mm2_pooled", "").replace("Imm-", ""))
        bottom += np.nan_to_num(values)
    ax.set_xticks(range(len(summary)), summary.donor_id, rotation=60, fontsize=12)
    for tick, group in zip(ax.get_xticklabels(), summary.disease_status):
        tick.set_color(COLORS[group])
    ax.set_ylabel("Immune cells per mm² of true Islet\n(pooled, by immune subclass)")
    ax.legend(ncol=4, loc="upper left", frameon=False, fontsize=12)
    package.save(fig, "04_islet_immune_subclasses", "Immune subclasses inside true Islets",
                 "Pooled within-Islet density of each refined immune subclass per donor.",
                 "Subclass is the immune_class of a cell whose broad Immune label is kept by refinement.",
                 {"donor_summary": summary[["donor_id", "disease_status", *immune_columns]]})

    rna = pd.read_csv(atlas / "rna_section_composition.csv", dtype={"donor_id": str})
    pancreas = rna.loc[rna.lane.eq("pancreas_RNA") & rna.named_lineage_eligible]
    table = pancreas.pivot_table(index=["donor_id", "disease_status"], columns="annotation", values="cells",
                                 aggfunc="sum", fill_value=0)
    fractions = table.div(table.sum(axis=1), axis=0).reset_index()
    fractions["group_order"] = fractions.disease_status.map({g: i for i, g in enumerate(GROUPS)})
    fractions = fractions.sort_values(["group_order", "donor_id"]).reset_index(drop=True)
    lineages = [c for c in table.columns if c != "Exocrine"]
    fig, ax = plt.subplots(figsize=SLIDE)
    bottom = np.zeros(len(fractions))
    cmap = plt.get_cmap("tab20")
    for k, lineage in enumerate(lineages):
        ax.bar(range(len(fractions)), fractions[lineage], bottom=bottom, color=cmap(k), width=0.8, label=lineage)
        bottom += fractions[lineage].to_numpy()
    ax.set_xticks(range(len(fractions)), fractions.donor_id, rotation=60, fontsize=12)
    for tick, group in zip(ax.get_xticklabels(), fractions.disease_status):
        tick.set_color(COLORS[group])
    ax.set_ylabel("Fraction of named-lineage RNA cells\n(Exocrine omitted from the stack)")
    ax.legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=12)
    fig.subplots_adjust(left=0.09, right=0.74, bottom=0.18, top=0.97)
    package.save(fig, "05_rna_pancreas_lineages", "RNA lineages in each pancreas section",
                 "Approved Xenium pancreas lineages per section as a fraction of all named-lineage cells; the remainder "
                 "to 1 is Exocrine (one label; no acinar/ductal claim).",
                 "Named-lineage eligible cells of the current approved release; unresolved and uncalled units are retained "
                 "in the source table but are not named here.",
                 {"fractions": fractions.drop(columns="group_order")})


# ---------------------------------------------------------------- hypotheses

def h1_figure(context, package, test):
    stage = context.task("h1")
    endpoints = pd.read_csv(stage / "donor_endpoints.csv", dtype={"donor_id": str})
    fig, axes = plt.subplots(2, 2, figsize=SLIDE, constrained_layout=True)
    row, example = pick_example(endpoints, stage / "examples")
    ax = axes[0, 0]
    if example:
        for wkt in example["vascular_wkt"]:
            outline(ax, shapely.from_wkt(wkt), "#4b8fbc", 0.8)
        for wkt in example["neural_wkt"]:
            outline(ax, shapely.from_wkt(wkt), "#c29a43", 1.0)
        xy, label = np.asarray(example["cell_xy"]), np.asarray(example["cell_label"])
        for name in ("", "Beta", "Alpha", "Delta", "Immune"):
            k = label == name
            ax.scatter(xy[k, 0], xy[k, 1], s=6 if name != "Immune" else 22, color=LABEL_COLORS[name],
                       label=name or "other", zorder=3 if name == "Immune" else 2)
        outline(ax, shapely.from_wkt(example["outline_wkt"]))
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(f"{row.donor_id} ({row.disease_status}): vessel (blue) and nerve (gold) objects", fontsize=13)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(fontsize=10, ncol=1, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, markerscale=1.5)
        ax = axes[0, 1]
        arc = np.asarray(example["arc_fraction"])
        ax.plot(arc, example["vascular_profile"], color="#4b8fbc", label="vascular proximity")
        ax.plot(arc, example["neural_profile"], color="#c29a43", label="neural proximity")
        twin = ax.twinx()
        twin.plot(arc, example["residual"], color="#c0392b", label="immune − all-cell exposure")
        twin.set_ylabel("immune exposure residual", fontsize=13)
        ax.set_xlabel("fraction of the islet's outer boundary")
        ax.set_ylabel("exp(−distance / cell diameter)", fontsize=13)
        ax.set_title("Null: every circular shift of the structural profile\n(red: immune − all-cell exposure)", fontsize=13)
        ax.legend(fontsize=10, loc="upper left", bbox_to_anchor=(0.0, -0.22), ncol=2, frameon=False)
    strip(axes[1, 0], endpoints, "donor_effect", "vascular − neural alignment\nexcess (donor mean)")
    effect_panel(axes[1, 1], test, "difference in donor means")
    package.save(fig, "07_H1_neurovascular_boundary", "H1 · Neurovascular boundary access",
                 "Do islet-boundary immune cells align with vascular more than neural structure, and more in T1D? "
                 "Per islet, immune-minus-all-cell boundary exposure is scored against the vascular-minus-neural proximity "
                 "profile and compared with all circular shifts of that profile.",
                 "Example chosen by support (median islets available) and median islet area, never by effect. Kernel: "
                 "exp(−distance to the outer boundary / donor mean cell diameter). Multipart outlines have no single "
                 "arclength and are unavailable. " + unavailable_text(endpoints),
                 {"donor_endpoints": endpoints})


def h2_figure(context, package, test):
    stage = context.task("h2")
    endpoints = pd.read_csv(stage / "donor_endpoints.csv", dtype={"donor_id": str})
    fig, axes = plt.subplots(2, 2, figsize=SLIDE, constrained_layout=True)
    row, example = pick_example(endpoints, stage / "examples")
    if example:
        ax = axes[0, 0]
        xy, label = np.asarray(example["xy"]), np.asarray(example["labels"])
        for i, j in np.asarray(example["core_edges"], int):
            if {label[i], label[j]} == {"Beta", "Immune"}:
                ax.plot(xy[[i, j], 0], xy[[i, j], 1], color="#c0392b", linewidth=1.2, alpha=0.7)
        for name in ("", "Beta", "Alpha", "Delta", "Immune"):
            k = label == name
            ax.scatter(xy[k, 0], xy[k, 1], s=26, color=LABEL_COLORS[name], label=name or "other", zorder=3)
        outline(ax, shapely.from_wkt(example["outline_wkt"]))
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{row.donor_id} ({row.disease_status}): Beta–immune edges in red", fontsize=13)
        ax.legend(fontsize=10, ncol=1, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        ax = axes[0, 1]
        ax.hist(example["null_joint"], bins=25, color="#b7c2c8")
        ax.axvline(example["observed"], color="#c0392b", linewidth=2.5)
        ax.set_xlabel("joint standardized topology score")
        ax.set_ylabel("shuffled arrangements")
        ax.set_title("Null: subtype labels shuffled within radial strata", fontsize=13)
    strip(axes[1, 0], endpoints, "donor_effect", "joint topology excess\n(donor mean of islets)")
    effect_panel(axes[1, 1], test, "difference in donor means")
    components = (f"Component effects (T1D − ND): fragmentation {test.get('fragmentation_z_effect', np.nan):.3g}, "
                  f"immune interface {test.get('interface_z_effect', np.nan):.3g}.")
    package.save(fig, "08_H2_beta_topology", "H2 · Beta topology beyond depletion",
                 "Is the Beta-cell graph more fragmented and more immune-facing in T1D than its own composition and radial "
                 "placement predict? " + components,
                 "Graph: k=6 nearest centroids inside the outline, edges weighted by exp(−length / cell diameter). The joint "
                 "score needs both components; islets without immune cells or without Beta edges cannot vary the interface "
                 "and are unavailable for the joint. " + unavailable_text(endpoints),
                 {"donor_endpoints": endpoints})


def h3_figure(context, package, test):
    stage = context.task("h3")
    endpoints = pd.read_csv(stage / "donor_endpoints.csv", dtype={"donor_id": str})
    primary = endpoints.loc[endpoints.configuration.eq("primary")]
    fig, axes = plt.subplots(2, 2, figsize=SLIDE, constrained_layout=True)
    row, example = pick_example(primary.rename(columns={"n_islets": "n_islets_available"}), stage / "examples")
    if example:
        ax = axes[0, 0]
        xy, residual = np.asarray(example["xy"]), np.asarray(example["hla_residual"], float)
        limit = np.nanmax(np.abs(residual)) or 1
        points = ax.scatter(xy[:, 0], xy[:, 1], c=residual, cmap="coolwarm", vmin=-limit, vmax=limit, s=36)
        outline(ax, shapely.from_wkt(example["outline_wkt"]))
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(points, ax=ax, shrink=0.8, label="HLA-A residual")
        ax.set_title(f"{row.donor_id} ({row.disease_status}): residual HLA-A after nuisance fit", fontsize=13)
    ax = axes[0, 1]
    sens = endpoints.pivot_table(index="donor_id", columns="configuration", values="donor_effect")
    for configuration in [c for c in sens.columns if c != "primary"]:
        ax.scatter(sens["primary"], sens[configuration], s=30, label=configuration)
    values = sens.to_numpy()[np.isfinite(sens.to_numpy())]
    low, high = (values.min(), values.max()) if values.size else (-1, 1)
    pad = 0.05 * (high - low or 1)
    ax.plot([low - pad, high + pad], [low - pad, high + pad], color="#888888", linewidth=1)
    ax.set_xlabel("primary (k=4, 50 µm)")
    ax.set_ylabel("declared sensitivity")
    ax.set_title("Graph k and neighbourhood sensitivities", fontsize=13)
    ax.legend(fontsize=10, frameon=False)
    strip(axes[1, 0], primary, "donor_effect", "residual HLA-A autocorrelation\nexcess (donor mean)")
    effect_panel(axes[1, 1], test, "difference in donor means")
    package.save(fig, "09_H3_HLA_A_domains", "H3 · Endocrine HLA-A domains beyond contact",
                 "Is endocrine HLA-A spatially clustered among same-subtype neighbours beyond what local immune exposure, "
                 "subtype, islet position and structure explain, and more so in T1D?",
                 "Freedman–Lane: nuisance residuals shuffled within islet and subtype, then refitted; homotypic k-NN graph "
                 "inside the outline with an edge cap of two mean cell diameters. " + unavailable_text(primary),
                 {"donor_endpoints": endpoints})


def h4_figure(context, package, test):
    stage = context.task("h4")
    endpoints = pd.read_csv(stage / "donor_endpoints.csv", dtype={"donor_id": str})
    primary = endpoints.loc[endpoints.configuration.eq("primary")]
    fig, axes = plt.subplots(2, 2, figsize=SLIDE, constrained_layout=True)
    slopes = ["lowHLA_lowCyto", "highHLA_lowCyto", "lowHLA_highCyto", "highHLA_highCyto"]
    ax = axes[0, 0]
    available = primary.loc[primary.status.eq("available")]
    for group in ("ND", "T1D"):
        z = available.loc[available.disease_status.eq(group)]
        if len(z):
            means = [z[f"PD1_simple_slope_{s}"].mean() for s in slopes]
            ax.plot(range(4), means, marker="o", color=COLORS[group], label=f"{group} (n={len(z)})", linewidth=2)
    ax.set_xticks(range(4), ["low HLA\nlow cyto", "high HLA\nlow cyto", "low HLA\nhigh cyto", "high HLA\nhigh cyto"], fontsize=11)
    ax.axhline(0, color="#888888", linestyle="--", linewidth=1)
    ax.set_ylabel("PD-1-neighbour slope\non endocrine PD-L1")
    ax.set_title("Simple slopes at ±1 SD (all lower-order terms)", fontsize=13)
    ax.legend(frameon=False, fontsize=11)
    ax = axes[0, 1]
    sens = endpoints.pivot_table(index="donor_id", columns="configuration", values="donor_effect")
    for configuration in [c for c in sens.columns if c != "primary"]:
        ax.scatter(sens["primary"], sens[configuration], s=30, label=configuration)
    values = sens.to_numpy()[np.isfinite(sens.to_numpy())]
    low, high = (values.min(), values.max()) if values.size else (-1, 1)
    pad = 0.05 * (high - low or 1)
    ax.plot([low - pad, high + pad], [low - pad, high + pad], color="#888888", linewidth=1)
    ax.set_xlabel("primary (50 µm)")
    ax.set_ylabel("declared sensitivity")
    ax.set_title("Neighbourhood radius sensitivities", fontsize=13)
    ax.legend(fontsize=10, frameon=False)
    strip(axes[1, 0], primary, "donor_effect", "triple-interaction coefficient\n(null-centred)")
    effect_panel(axes[1, 1], test, "difference in donor means")
    package.save(fig, "10_H4_checkpoint_interaction", "H4 · Non-additive checkpoint alignment",
                 "Does the PD-1-neighbour × endocrine HLA-A × local cytotoxic-state interaction predict endocrine PD-L1 "
                 "beyond all lower-order terms, and less so in T1D (predicted negative)?",
                 "Islet fixed effects; null retains the pairwise-only fit and flips residual signs by whole islet. The "
                 "interaction alone does not establish net checkpoint slope weakening. " + unavailable_text(primary),
                 {"donor_endpoints": endpoints})


def h5_figure(context, package, test):
    stage = context.task("h5")
    endpoints = pd.read_csv(stage / "donor_endpoints.csv", dtype={"donor_id": str})
    primary = endpoints.loc[endpoints.configuration.eq("primary")]
    fig, axes = plt.subplots(2, 2, figsize=SLIDE, constrained_layout=True)
    ax = axes[0, 0]
    reasons = primary.reason.fillna("").replace("", "available").value_counts()
    ax.barh(range(len(reasons)), reasons.to_numpy(), color="#477b8b")
    ax.set_yticks(range(len(reasons)), [r[:55] for r in reasons.index], fontsize=10)
    ax.set_xlabel("pancreas RNA sections")
    ax.set_title("Computability across the 25 sections", fontsize=13)
    ax = axes[0, 1]
    nodes = sorted((stage / "organizer_nodes").glob("*.parquet"))
    if nodes:
        frame = pd.read_parquet(nodes[0])
        ax.scatter(frame.interaction_predictor, frame.T_neighbour_fraction, s=12, color="#7f8c8d")
        ax.set_xlabel("CCL19/21 state × B-cell LTB exposure (z)")
        ax.set_ylabel("T-neighbour fraction")
        ax.set_title(f"Organizers of donor {frame.donor_id.iloc[0]}", fontsize=13)
    else:
        ax.text(0.5, 0.5, "No computable organizer motif", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
    if np.isfinite(primary.donor_effect.astype(float)).any():
        strip(axes[1, 0], primary, "donor_effect", "organizer × B-LTB interaction\non T adjacency")
    else:
        sensitivity = endpoints.loc[endpoints.status.eq("available")].groupby("configuration").donor_id.apply(
            lambda d: ", ".join(sorted(d)))
        extra = "; ".join(f"{c}: {d}" for c, d in sensitivity.items()) or "none"
        axes[1, 0].text(0.5, 0.5, f"0 of {len(primary)} donors computable at the primary settings\n\n"
                        f"Computable under a declared sensitivity:\n{extra}", ha="center", va="center",
                        transform=axes[1, 0].transAxes, fontsize=15, wrap=True)
        axes[1, 0].axis("off")
    effect_panel(axes[1, 1], test, "difference in donor means")
    package.save(fig, "11_H5_organizer_motifs", "H5 · Higher-order organizer motifs",
                 "Do CCL19/21 organizer cells inside islets with LTB-expressing B-cell neighbours have more T-cell neighbours "
                 "than their pairwise effects predict, and more so in T1D?",
                 "Organizers inside exact-serial post-IF islet masks (localization only); released identities, no marker "
                 "gate; whole-islet Rademacher null with the pairwise baseline. " + unavailable_text(primary),
                 {"donor_endpoints": endpoints, "computability": reasons.rename_axis("reason").reset_index(name="sections")})


def h6_figure(context, package, test):
    stage = context.task("h6")
    inference = context.task("inference")
    sections = pd.read_csv(stage / "pln_section_architecture.csv", dtype={"donor_id": str})
    paired = pd.read_csv(stage / "paired_donors.csv", dtype={"donor_id": str})
    fig, axes = plt.subplots(2, 2, figsize=SLIDE, constrained_layout=True)
    states = sorted((stage / "pln_states").glob("*.parquet"))
    ax = axes[0, 0]
    if states:
        primary = sections.loc[sections.configuration.eq("primary") & sections.status.eq("available")]
        median = primary.n_T_cells.median()
        chosen = primary.assign(gap=(primary.n_T_cells - median).abs()).sort_values(["gap", "section_id"]).iloc[0]
        frame = pd.read_parquet(stage / "pln_states" / f"{chosen.section_id}.parquet")
        limit = np.nanquantile(np.abs(frame.residual_contrast), 0.98)
        points = ax.scatter(frame.x_um, frame.y_um, c=frame.residual_contrast, cmap="PuOr", vmin=-limit, vmax=limit, s=1)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(points, ax=ax, shrink=0.8, label="memory − cytotoxic\n(residual)")
        ax.set_title(f"pLN section {chosen.section_id} ({chosen.disease_status}): T-state map", fontsize=13)
    architecture = paired.rename(columns={"architecture_excess": "donor_effect"})
    strip(axes[0, 1], architecture.assign(donor_effect=architecture.donor_effect.where(architecture.status.eq("available"))),
          "donor_effect", "pLN T-state assortativity\nexcess (donor)")
    ax = axes[1, 0]
    pairs_path = inference / "H6_primary_pairs.csv"
    if pairs_path.is_file():
        pairs = pd.read_csv(pairs_path, dtype={"donor_id": str})
        for group, z in pairs.groupby("disease_status"):
            ax.scatter(z.pLN_rank_residual, z.H3_rank_residual, s=70, color=COLORS[group], label=group)
            for _, r in z.iterrows():
                ax.annotate(r.donor_id, (r.pLN_rank_residual, r.H3_rank_residual), xytext=(4, 2), textcoords="offset points", fontsize=9)
        ax.set_xlabel("pLN architecture (rank residual)")
        ax.set_ylabel("pancreatic H3 (rank residual)")
        ax.legend(frameon=False)
        ax.set_title("Paired ND/T1D donors after conditioning", fontsize=13)
    else:
        ax.text(0.5, 0.5, "No computable paired association", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
    effect_panel(axes[1, 1], test, "partial rank correlation")
    package.save(fig, "12_H6_linked_pLN_islet", "H6 · Linked pLN and islet architecture",
                 "Across paired donors, does pLN T-cell state architecture track pancreatic HLA-A domain excess (H3) beyond "
                 "disease and mean T states?",
                 "pLN: TCF7/LEF1 minus TOX/GZMB/PRF1/NKG7 residual contrast; assortativity over 8 nearest T neighbours "
                 "(≤50 µm) against a QC-stratified label null; sections averaged, coordinates never merged. Association: "
                 "within-disease Freedman–Lane permutation. Cross-sectional; no trafficking or clonality claim.",
                 {"pLN_sections": sections, "paired_donors": paired})


def direction(row) -> str:
    if row["status"] != "available":
        return "—"
    predicted = row["effect"] < 0 if row["predicted_direction"] == "negative" else row["effect"] > 0
    return f"as predicted ({row['predicted_direction']})" if predicted else f"opposite (predicted {row['predicted_direction']})"


def summary_figure(context, package, tests: pd.DataFrame):
    tests = tests.assign(direction=[direction(r) for r in tests.to_dict("records")])
    rows = []
    for row in tests.to_dict("records"):
        donors = (f"{row.get('n_pairs', 0):.0f} pairs" if row["hypothesis"] == "H6"
                  else f"{row.get('n_ND', 0):.0f} ND, {row.get('n_T1D', 0):.0f} T1D")
        available = row["status"] == "available"
        rows.append([textwrap.fill(f"{row['hypothesis']}  {row['name']}", 30), donors,
                     f"{row['effect']:.3g}\n[{row['CI95_low']:.3g}, {row['CI95_high']:.3g}]" if available else "unavailable",
                     f"{row['p']:.3g}" if available else "—", f"{row['q_BH6']:.3g}" if available else "—",
                     row["direction"].replace(" (", "\n("),
                     ("inconclusive" if row["q_BH6"] >= 0.05 else textwrap.fill(row["interpretation"].split(";")[0], 26))
                     if available else f"unavailable\n(see the {row['hypothesis']} slide)"])
    fig, ax = plt.subplots(figsize=SLIDE)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.98, bottom=0.12)
    ax.axis("off")
    table = ax.table(cellText=rows, colLabels=["Hypothesis", "Donors", "T1D − ND\n[95% CI]", "p", "BH q",
                                               "Direction", "Reading"],
                     colWidths=[0.25, 0.11, 0.16, 0.06, 0.06, 0.17, 0.19], loc="upper center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(13)
    table.scale(1, 4.6)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if r == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#eef2f5")
        elif tests.iloc[r - 1].status != "available":
            cell.get_text().set_color("#8a3b3b")
    fig.text(0.01, 0.02, "BH over the fixed six. H6 effect is a partial rank correlation; the others are differences in donor "
             "means,\nin each hypothesis's own units. Inconclusive is not refutation.", fontsize=12, color="#555555")
    package.save(fig, "06_six_hypotheses_summary", "Six spatial hypotheses: donor-level results",
                 "One fixed family of six; BH across the six. An unavailable hypothesis is shown as unavailable, never as a "
                 "negative result. Nonsignificance is inconclusive, not refutation.",
                 "Retrospective/exploratory: these data were inspected before the protocol was frozen. Effects are in each "
                 "hypothesis's own units (see its slide). Acquisition and disease alias at the coarse level; stain-group "
                 "stratified permutations are in sensitivity_tests.csv.",
                 {"primary_tests": tests.drop(columns="direction")})


# ---------------------------------------------------------------- deeper phenotypes, 6551, CD8

STATE_COLORS = {"TCF1+TOX-": "#4b8fbc", "TCF1+TOX+": "#8e5ea2", "TCF1-TOX+": "#c0392b", "TCF1-TOX-": "#c7c7c7"}


def ordered(table: pd.DataFrame) -> pd.DataFrame:
    return table.assign(order=table.disease_status.map({g: i for i, g in enumerate(GROUPS)})).sort_values(
        ["order", "donor_id"]).drop(columns="order").reset_index(drop=True)


def deeper_figures(context, package):
    atlas = context.task("atlas")
    states = ordered(pd.read_csv(atlas / "protein_t_cell_states.csv", dtype={"donor_id": str}))
    fig, axes = plt.subplots(2, 1, figsize=SLIDE, sharex=True, constrained_layout=True)
    for ax, compartment in zip(axes, ("inside true Islets", "outside true Islets")):
        z = states.loc[states.compartment.eq(compartment)].reset_index(drop=True)
        bottom = np.zeros(len(z))
        for name, color in STATE_COLORS.items():
            values = np.nan_to_num(z[f"share_{name}"].to_numpy(float))
            ax.bar(range(len(z)), values, bottom=bottom, color=color, width=0.8, label=name)
            bottom += values
        for i, n in enumerate(z.T_cells_TCF1_TOX_measured):
            ax.text(i, 1.01, f"{n:,}", ha="center", va="bottom", fontsize=8, rotation=90, color="#555555")
        ax.set_ylim(0, 1.18)
        ax.set_ylabel(f"T cells {compartment}\n(share of TCF1/TOX-measured)", fontsize=13)
        ax.set_xticks(range(len(z)), z.donor_id, rotation=60, fontsize=11)
        for tick, group in zip(ax.get_xticklabels(), z.disease_status):
            tick.set_color(COLORS[group])
    axes[0].legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.28), frameon=False)
    package.save(fig, "13_protein_T_cell_states", "Protein T-cell states inside and outside true Islets",
                 "TCF1/TOX states of refined T cells per donor (numbers: T cells with both markers measured). Positive "
                 "means at least 0.10 on the cohort scale, the level the immune definitions use.",
                 "T cells: broad Immune kept by refinement, immune_class Imm-T-cell. No distance band is drawn: a cell "
                 "is inside an original true-Islet outline or it is not. Descriptive; the frozen primaries are unchanged.",
                 {"t_cell_states": states})

    subtypes = pd.read_csv(atlas / "protein_immune_subtypes.csv", dtype={"donor_id": str})
    inside = subtypes.loc[subtypes.compartment.eq("inside true Islets")].copy()
    # Subtype rules exist only under T cells and macrophages; every other cell is shown by its class.
    no_subtype = inside.subtype.eq("no subtype rule fired")
    inside["group"] = np.where(no_subtype, inside.immune_class.str.replace("Imm-", "", regex=False) + " (class)",
                               inside.subtype)
    table = inside.pivot_table(index=["donor_id", "disease_status"], columns="group", values="cells",
                               aggfunc="sum", fill_value=0)
    fractions = ordered(table.div(table.sum(axis=1), axis=0).reset_index())
    names = [c for c in table.columns]
    fig, ax = plt.subplots(figsize=SLIDE)
    bottom = np.zeros(len(fractions))
    cmap = plt.get_cmap("tab20")
    for k, name in enumerate(names):
        ax.bar(range(len(fractions)), fractions[name], bottom=bottom, color=cmap(k), width=0.8, label=name)
        bottom += fractions[name].to_numpy()
    ax.set_xticks(range(len(fractions)), fractions.donor_id, rotation=60, fontsize=12)
    for tick, group in zip(ax.get_xticklabels(), fractions.disease_status):
        tick.set_color(COLORS[group])
    ax.set_ylabel("Share of refined immune cells inside true Islets")
    ax.legend(ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=11)
    fig.subplots_adjust(left=0.08, right=0.74, bottom=0.18, top=0.97)
    package.save(fig, "14_protein_immune_subtypes_in_islets", "Immune subtypes inside true Islets",
                 "Second-pass immune subtypes of refined immune cells inside original true Islets, per donor.",
                 "Subtypes come from the immune pass (immune_subtype.py) over the Immune class, which has subtype rules "
                 "only under T cells and macrophages; '(class)' marks cells shown by their immune class because no subtype "
                 "rule named them. CD8-defined subtypes are provisional until the CD8 ruling.",
                 {"immune_subtypes": subtypes})

    rna = pd.read_csv(atlas / "rna_t_cell_states.csv", dtype={"donor_id": str})
    fig, axes = plt.subplots(1, 4, figsize=SLIDE, constrained_layout=True)
    for column, (ax, lane, field, label) in enumerate(zip(
            axes, ("pancreas_RNA", "pancreas_RNA", "pLN_RNA", "pLN_RNA"),
            ("mean_memory_program", "mean_cytotoxic_program") * 2,
            ("pancreas T: memory\n(TCF7/LEF1)", "pancreas T: cytotoxic\n(TOX/GZMB/PRF1/NKG7)",
             "pLN T: memory\n(TCF7/LEF1)", "pLN T: cytotoxic\n(TOX/GZMB/PRF1/NKG7)"))):
        z = rna.loc[rna.lane.eq(lane)].groupby(["donor_id", "disease_status"], as_index=False)[field].mean()
        strip(ax, z, field, label if column == 0 else label, label_donors=False)
        ax.set_ylabel(label, fontsize=12)
    package.save(fig, "15_rna_T_cell_states", "RNA T-cell programs per donor",
                 "Mean log-normalized memory and cytotoxic programs of released RNA T cells; sections of one donor are "
                 "averaged. Pancreas and pLN are separate sections and are never pooled.",
                 "Programs: mean of log1p(1e4 × count / library) over the genes. Identity is the approved release, with "
                 "no marker gate. Descriptive.", {"rna_t_states": rna})


def crossstudy_figure(context, package):
    stage = context.task("crossstudy")
    shares = pd.read_csv(stage / "tcf1_tox_shares_by_tissue.csv")
    distance = pd.read_csv(stage / "tcf1_tox_shares_by_islet_distance.csv")
    agreement = pd.read_csv(stage / "agreement.csv")
    studies = list(dict.fromkeys(shares.study))
    fig, axes = plt.subplots(1, 3, figsize=SLIDE, constrained_layout=True)
    for ax, tissue in zip(axes[:2], ("Pancreas", "pLN")):
        z = shares.loc[shares.tissue.eq(tissue)].set_index("study")
        width = 0.8 / len(z)
        for k, (study, row) in enumerate(z.iterrows()):
            ax.bar(np.arange(3) + k * width - 0.4 + width / 2, [row[s] for s in ("TCF1+TOX-", "TCF1+TOX+", "TCF1-TOX+")],
                   width=width, label=f"{study} (n={int(row.n):,})",
                   color=["#222222", "#3a6e94", "#9cc1dc"][studies.index(study) % 3])
        ax.set_xticks(range(3), ["TCF1+\nTOX−", "TCF1+\nTOX+", "TCF1−\nTOX+"])
        ax.set_ylabel("% of TCF1- and/or TOX-positive T cells", fontsize=13)
        ax.set_ylim(0, 80)
        ax.set_title(f"6551 {tissue}", fontsize=15)
        ax.legend(fontsize=9, frameon=False, loc="upper left")
    ax = axes[2]
    for study, z in distance.groupby("study", sort=False):
        ax.plot(range(len(z)), z["TCF1+TOX+"], marker="o", label=study, linewidth=2,
                color="#222222" if study.startswith("Peters") else "#3a6e94")
    ax.set_xticks(range(len(z)), z.distance_bin, rotation=30, fontsize=11)
    ax.set_xlabel("distance to islet, µm (< 0 inside)", fontsize=13)
    ax.set_ylabel("% TCF1+TOX+ of TCF1/TOX-positive T cells", fontsize=13)
    ax.set_title("6551 pancreas by islet distance", fontsize=15)
    ax.legend(fontsize=9, frameon=False)
    text = "; ".join(f"{r.comparison.split(' vs ')[0].replace('this pipeline, ', '')} ({r.comparison.rsplit(', ', 1)[1]}): "
                     f"TVD {r.total_variation_distance_percent_points:.1f} pts" for r in agreement.itertuples())
    package.save(fig, "16_crossstudy_6551", "Cross-study accuracy on donor 6551",
                 "The same donor analysed by Peters et al. (Figure 7) and by this pipeline: T-cell TCF1/TOX state shares "
                 "per tissue and by distance to islets. " + text,
                 "Audit context only: agreement is reported, never used to tune or gate. Different segmentation, gates and "
                 "islet definitions are part of what is compared. TVD = total variation distance between the two state "
                 "distributions, in percentage points.", {"shares": shares, "distance": distance, "agreement": agreement})


def cd8_figure(context, package):
    stage = context.task("cd8")
    package.add_rendered("17_cd8_roc_decision", stage, "cd8_roc_per_donor", "CD8 per donor: for the shift-or-not ruling",
                         "Each donor's ROC for unshifted CD8 (CD3e+CD4− T cells against CD4+ T cells and against "
                         "hormone+ CD3e− cells), with the cohort level marked under UniFORM's shift (red) and with no shift (square).",
                         "Decision support for the maintainer; nothing is decided by this slide. CD3e+CD4− includes "
                         "double-negative T cells, so each curve is a lower bound on what CD8 can do.",
                         [stage / "cd8_roc_summary.csv", stage / "cd8_roc_definition.json"])


def render(context: Context, stage: Path) -> dict:
    package = Package(stage)
    atlas_figures(context, package)
    tests = pd.read_csv(context.task("inference") / "primary_tests.csv")
    summary_figure(context, package, tests)
    by_id = {row["hypothesis"]: row for row in tests.to_dict("records")}
    for h, function in (("H1", h1_figure), ("H2", h2_figure), ("H3", h3_figure), ("H4", h4_figure),
                        ("H5", h5_figure), ("H6", h6_figure)):
        function(context, package, by_id[h])
    deeper_figures(context, package)
    crossstudy_figure(context, package)
    cd8_figure(context, package)
    package.slides.sort(key=lambda s: s["name"])
    write_json(stage / "slides.json", {"title": "IO60: atlas and six spatial hypotheses",
                                       "subtitle": "Retrospective/exploratory; donor is the unit of inference",
                                       "slides": package.slides})
    environment = {k: v for k, v in os.environ.items() if not k.startswith("PYTHON")}
    subprocess.run([str(PRESENTATION_PYTHON), "-I", str(Path(__file__).with_name("deck.py")), str(stage)],
                   check=True, env=environment)
    return {"figures": len(package.slides), "deck": "deck.pptx", "pdf": "deck.pdf"}
