#!/usr/bin/env python3
"""Plot donor-level protein literature endpoints from the protein endpoint builder."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PETERS_MARKERS = ("CD8", "TCF1", "TOX", "Granzyme_B", "LAG3", "PD_1", "CD39", "IFNG", "ICOS")
INFLAMMATORY_MARKERS = ("IFNG", "IDO1", "HLA_DR", "HLA_A", "VISTA", "PD_L1")
COLORS = {"pancreas": "#337f8c", "pLN": "#cf7b4b"}


def points_and_summary(ax, frame: pd.DataFrame, marker_order: tuple[str, ...], *,
                       y_column: str, summary: pd.DataFrame, family: str,
                       tissue: str | None = None, reference: float | None = None,
                       ylabel: str, title: str, y_limits: tuple[float, float] | None = None):
    for marker_i, marker in enumerate(marker_order):
        if tissue is None:
            part = frame[frame.marker.eq(marker)]
            tissues = sorted(part.tissue.unique())
        else:
            part = frame[frame.marker.eq(marker) & frame.tissue.eq(tissue)]
            tissues = [tissue]
        for tissue_i, tissue_name in enumerate(tissues):
            vals = part.loc[part.tissue.eq(tissue_name), y_column].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            offset = (-0.14 if tissue_name == "pancreas" else 0.14) if tissue is None else 0.0
            x = marker_i + offset
            if len(vals):
                jitter = np.linspace(-0.045, 0.045, len(vals)) if len(vals) > 1 else np.array([0.0])
                ax.scatter(np.full(len(vals), x) + jitter, vals, color=COLORS.get(tissue_name, "#777777"),
                           s=14, alpha=0.62, edgecolors="none", zorder=2)
            if family:
                selected = summary[(summary.endpoint_family == family)
                                   & summary.marker.eq(marker)
                                   & summary.tissue.eq(tissue_name)]
                if len(selected):
                    s = selected.iloc[0]
                    mean = float(s.donor_balanced_mean)
                    lo, hi = float(s.bootstrap_95ci_low), float(s.bootstrap_95ci_high)
                    if np.isfinite(mean):
                        ax.errorbar(x, mean, yerr=[[max(0, mean-lo)], [max(0, hi-mean)]],
                                    fmt="D", color="#222222", markersize=4, capsize=2,
                                    linewidth=1, zorder=3)
    if reference is not None:
        ax.axhline(reference, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_xticks(np.arange(len(marker_order)), labels=[m.replace("_", " ") for m in marker_order],
                  rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    if y_limits:
        ax.set_ylim(*y_limits)
    ax.grid(axis="y", color="#e2e2e2", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--analysis", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True, help="figure basename; PNG and SVG are written")
    args = ap.parse_args()
    donor = pd.read_csv(args.analysis / "protein_literature_donor_marker_endpoints.csv")
    inflammation = pd.read_csv(args.analysis / "protein_literature_donor_inflammation_endpoints.csv")
    summary = pd.read_csv(args.analysis / "protein_literature_donor_balanced_summary.csv")
    manifest = json.loads((args.analysis / "protein_literature_analysis_manifest.json").read_text())
    fig, axes = plt.subplots(3, 1, figsize=(13.0, 11.0), constrained_layout=False,
                             gridspec_kw={"height_ratios": [1.0, 1.0, 0.85], "hspace": 0.42})
    points_and_summary(
        axes[0], donor, PETERS_MARKERS, y_column="marker_auc_strict_cd8_vs_known_noncd8_t",
        summary=summary, family="Peters_T_cell_marker_AUC", reference=0.5,
        ylabel="Within-donor marker AUC", title="Peters: marker signal in native strict CD8 vs known non-CD8 T-cell calls",
        y_limits=(0, 1))
    axes[0].legend(handles=[plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[t],
                                       label=t, markersize=6) for t in ("pancreas", "pLN")]
                   + [plt.Line2D([0], [0], marker="D", color="#222222", label="Donor-balanced mean and 95% CI",
                                 markersize=4)], frameon=False, ncols=3, loc="upper left",
                   bbox_to_anchor=(0, 1.02))
    points_and_summary(
        axes[1], donor[donor.tissue.eq("pancreas")], PETERS_MARKERS,
        y_column="distance_spearman_marker_vs_signed_nearest_islet_um",
        summary=summary, family="Peters_continuous_distance_association", tissue="pancreas",
        ylabel="Within-donor Spearman ρ", title="Peters: continuous marker association with signed distance to nearest islet",
        y_limits=(-1, 1))
    axes[1].axhline(0, color="#555555", linestyle="--", linewidth=0.8)
    points_and_summary(
        axes[2], inflammation, INFLAMMATORY_MARKERS,
        y_column="marker_auc_immune_vs_other_native_typed", summary=summary,
        family="Medina_pancreatic_inflammatory_marker_distribution", tissue="pancreas",
        reference=0.5, ylabel="Within-donor marker AUC",
        title="Medina-Serpas: pancreatic inflammatory marker signal in Immune vs other native typed cells",
        y_limits=(0, 1))
    status = str(manifest.get("status", "unspecified"))
    fig.suptitle(f"Protein literature endpoints ({status})", x=0.07, y=0.995, ha="left",
                 fontsize=14, fontweight="bold")
    fig.text(0.07, 0.96,
             "Points are donors; diamonds are donor-balanced means with donor-bootstrap 95% CIs.\n"
             "Cross-tissue panels use within-tissue rank AUC, not absolute intensity. Native calls are not independent truth labels.\n"
             "Distance panels use continuous signed micrometres; no distance bins.",
             ha="left", va="top", fontsize=8.5, color="#555555")
    fig.subplots_adjust(left=0.09, right=0.98, top=0.90, bottom=0.08)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(args.out.with_suffix(".svg"), facecolor="white")
    plt.close(fig)
    print(f"wrote donor-level protein endpoint figure -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
