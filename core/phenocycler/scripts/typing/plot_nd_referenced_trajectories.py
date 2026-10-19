#!/usr/bin/env python3
"""Plot donor-specific paths in the frozen ND MFA representation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm, colors


OBJECT_CLASS_MARKERS = {
    "Islet": "o",
    "EndocrineCluster": "s",
    "EndocrineSingle": "^",
}

GROUP_LABELS = {
    "endocrine": "endocrine",
    "vascular": "vascular",
    "lymphatic": "lymphatic",
    "neural": "neural",
    "immune": "immune",
    "stromal": "stromal",
}


def component_axis_label(manifest: dict, component: int) -> str:
    """Name each frozen MFA axis by its two largest measured group contributions."""
    contributions = manifest.get("mfa_group_contributions_by_component", {})
    ranked = []
    for group, values in contributions.items():
        if component < len(values):
            value = float(values[component])
            if np.isfinite(value):
                ranked.append((value, GROUP_LABELS.get(group, group)))
    ranked.sort(reverse=True)
    if not ranked:
        return f"MFA {component + 1}"
    groups = " + ".join(label for _, label in ranked[:3])
    return f"MFA {component + 1} · {groups}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectories", required=True, type=Path,
                        help="donor_endocrine_object_trajectories.parquet")
    parser.add_argument("--out", required=True, type=Path,
                        help="figure basename; PNG and SVG are written")
    parser.add_argument("--view", choices=("zoom", "donor"), default="zoom",
                        help="both views use identical shared pooled 0.5th-99.5th percentile limits; donor includes fitted donors only")
    args = parser.parse_args()
    frame = pd.read_parquet(args.trajectories)
    manifest_path = args.trajectories.parent / "trajectory_manifest.json"
    trajectory_manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    representation_rank = trajectory_manifest.get("representation_rank")
    nd_center_label = (f"frozen {int(representation_rank)}D ND center"
                       if representation_rank is not None else "frozen ND reference center")
    x_axis_label = component_axis_label(trajectory_manifest, 0)
    y_axis_label = component_axis_label(trajectory_manifest, 1)
    x_col, y_col = "mfa_coordinate_1", "mfa_coordinate_2"
    required = {"donor_id", "donor_status", "component_id", "original_object_class",
                "endocrine_object_id", "component_root_object_id", "dpt_pseudotime",
                x_col, y_col, "within_donor_component_normalized_time", "trajectory_status"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"trajectory table is missing plot columns: {missing}")
    donors = sorted(frame.donor_id.astype(str).unique())
    if args.view == "donor":
        has_path = frame.trajectory_status.astype(str).eq("ordered")
        has_xy = (np.isfinite(pd.to_numeric(frame[x_col], errors="coerce"))
                  & np.isfinite(pd.to_numeric(frame[y_col], errors="coerce")))
        donors = sorted(frame.loc[has_path & has_xy, "donor_id"].astype(str).unique())
        ncols = min(7, max(1, len(donors)))
    else:
        ncols = 5
    nrows = int(np.ceil(len(donors) / ncols))
    if args.view == "donor":
        fig = plt.figure(figsize=(16, 2.5 * nrows))
        grid = fig.add_gridspec(nrows, ncols, left=0.055, right=0.925, top=0.82,
                                bottom=0.12, wspace=0.30, hspace=0.24)
    else:
        fig = plt.figure(figsize=(16, 2.5 * nrows))
        grid = fig.add_gridspec(nrows, ncols, left=0.065, right=0.925, top=0.90,
                                bottom=0.12, wspace=0.22, hspace=0.28)
    axes = []
    for row in range(nrows):
        axis_row = []
        for col in range(ncols):
            axis_row.append(fig.add_subplot(grid[row, col]))
        axes.append(axis_row)
    axes = np.asarray(axes, dtype=object)
    cmap = plt.get_cmap("viridis")
    norm = colors.Normalize(vmin=0, vmax=1)
    coordinate_mask = (
        np.isfinite(pd.to_numeric(frame[x_col], errors="coerce"))
        & np.isfinite(pd.to_numeric(frame[y_col], errors="coerce"))
        & frame.trajectory_status.astype(str).eq("ordered"))
    # For display only, shift the shared axes to a robust ND center: first take
    # each ND donor's median object coordinate, then the component-wise median
    # across ND donors. This keeps donor offsets comparable while centering the
    # reference cloud; stored MFA coordinates and DPT scores are untouched.
    nd_reference = frame.loc[
        coordinate_mask & frame.donor_status.astype(str).eq("ND")
    ].groupby("donor_id")[[x_col, y_col]].median()
    if nd_reference.empty:
        raise ValueError("no ordered ND coordinates are available to center the plot")
    display_center = nd_reference.median(axis=0).to_numpy(dtype=float)
    raw_coordinates = frame.loc[coordinate_mask, [x_col, y_col]].to_numpy(dtype=float)
    if not len(raw_coordinates):
        raise ValueError("no ordered trajectory coordinates are available to plot")
    # Use one pooled robust scale per component and the same transformation for
    # every donor. This display-only standardization lets both axes fill the
    # panel while preserving within-axis positions and cross-donor comparability.
    display_scale = np.quantile(np.abs(raw_coordinates - display_center), 0.975, axis=0)
    display_scale = np.maximum(display_scale, 1e-6)
    frame["_plot_mfa_1"] = (pd.to_numeric(frame[x_col], errors="coerce") - display_center[0]) / display_scale[0]
    frame["_plot_mfa_2"] = (pd.to_numeric(frame[y_col], errors="coerce") - display_center[1]) / display_scale[1]
    coordinates = frame.loc[coordinate_mask, ["_plot_mfa_1", "_plot_mfa_2"]].to_numpy(dtype=float)
    # The shared robust scale sets the 97.5th percentile absolute deviation of
    # each displayed coordinate to one; more distant coordinates stay in tables.
    plot_radius = np.ones(2, dtype=float)
    omitted = (np.abs(coordinates) > plot_radius).any(axis=1)
    omitted_count = int(omitted.sum())
    x_radius, y_radius = map(float, plot_radius)
    missing_donors = sorted(set(frame.donor_id.astype(str)) - set(donors))
    for ax, donor in zip(axes.flat, donors):
        block = frame[frame.donor_id.astype(str).eq(donor)]
        finite = (np.isfinite(pd.to_numeric(block[x_col], errors="coerce"))
                  & np.isfinite(pd.to_numeric(block[y_col], errors="coerce"))
                  & np.isfinite(pd.to_numeric(block.within_donor_component_normalized_time,
                                               errors="coerce"))
                  & block.trajectory_status.astype(str).eq("ordered"))
        plot = block.loc[finite]
        for _, component in plot.groupby("component_id", sort=True):
            for object_class, part in component.groupby("original_object_class", sort=True):
                rank = pd.to_numeric(
                    part.within_donor_component_normalized_time).to_numpy(float)
                rgba = cmap(norm(rank))
                rgba[:, 3] = 0.72
                marker = OBJECT_CLASS_MARKERS.get(str(object_class), "o")
                ax.scatter(part["_plot_mfa_1"], part["_plot_mfa_2"], c=rgba, s=5, marker=marker,
                           linewidths=0, rasterized=True)
            # Mark the algorithmic root and furthest DPT object. The graph is
            # higher-dimensional than this display; do not imply a straight 2D path.
            root_id = str(component.component_root_object_id.iloc[0])
            roots = component[component.endocrine_object_id.astype(str).eq(root_id)]
            if not roots.empty:
                ax.scatter(roots["_plot_mfa_1"], roots["_plot_mfa_2"], marker="*", s=95,
                           facecolor="white", edgecolor="#111111", linewidth=0.75,
                           zorder=5)
            terminal = component.loc[pd.to_numeric(component.dpt_pseudotime,
                                                   errors="coerce").idxmax()]
            ax.scatter([terminal["_plot_mfa_1"]], [terminal["_plot_mfa_2"]], marker="D", s=35,
                       facecolor="none", edgecolor="#111111", linewidth=0.8,
                       zorder=5)
        donor_status = str(block.donor_status.iloc[0]) if len(block) else "unknown"
        donor_status = donor_status if donor_status in {"ND", "AAb+", "T1D"} else "NA"
        object_count = len(block)
        component_count = int(block.loc[block.component_id.notna(), "component_id"].nunique()) if len(block) else 0
        component_note = str(component_count) if component_count else "NA"
        ax.set_title(f"{donor} · {donor_status}  |  N={object_count:,} · k={component_note}",
                     fontsize=8 if args.view != "donor" else 7.4, loc="left")
        if plot.empty:
            ax.text(0.5, 0.5, "No ordering fitted\nnot enough shared measured features\nobjects retained in data",
                    transform=ax.transAxes, ha="center", va="center", fontsize=7,
                    color="#555555")
        ax.axhline(0, color="#dedede", linewidth=0.45, zorder=0)
        ax.axvline(0, color="#dedede", linewidth=0.45, zorder=0)
        ax.set_xlim(-x_radius, x_radius)
        ax.set_ylim(-y_radius, y_radius)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax in list(axes.flat)[len(donors):]:
        ax.set_visible(False)
    view_label = {"zoom": "ND-centered robust shared-scale view",
                  "donor": "fitted-donor view (same scale)"}[args.view]
    fig.suptitle(f"Endocrine-object composition maps by donor ({view_label})",
                 x=0.06, y=0.995, ha="left", fontsize=14, fontweight="bold")
    note = (
        "Each dot is an endocrine object; shape shows Islet, Cluster, or Single. Color is within-donor/component order, not elapsed time or a shared clock.\n"
        f"★ marks the graph root nearest the {nd_center_label}; ◇ marks the highest-DPT object. Panel titles give total objects (N) and connected components (k).\n"
        "Axis labels name the top contributing measured groups, not validated biological programs. All panels share one display transform; axes are rescaled separately, so plotted distances are not original MFA distances.\n"
        f"Origin is the donor-balanced ND median; each axis uses its pooled 97.5th-percentile deviation. {omitted_count:,} distant objects are clipped but retained in tables."
    )
    fig.text(0.06, 0.972, note, fontsize=8.0, color="#444444", va="top")
    fig.text(0.49, 0.083, x_axis_label, ha="center", va="center", fontsize=9,
             color="#333333")
    fig.text(0.016, 0.49, y_axis_label, ha="center", va="center", rotation=90,
             fontsize=9, color="#333333")
    scalar = cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    colorbar_axis = fig.add_axes([0.947, 0.22, 0.015, 0.55])
    fig.colorbar(scalar, cax=colorbar_axis,
                 label="Relative pseudotime order within donor/component")
    from matplotlib.lines import Line2D
    legend = [
        Line2D([0], [0], marker=OBJECT_CLASS_MARKERS["Islet"], color="none",
               markerfacecolor="#777777", markersize=5, label="Islet"),
        Line2D([0], [0], marker=OBJECT_CLASS_MARKERS["EndocrineCluster"], color="none",
               markerfacecolor="#777777", markersize=5, label="Cluster"),
        Line2D([0], [0], marker=OBJECT_CLASS_MARKERS["EndocrineSingle"], color="none",
               markerfacecolor="#777777", markersize=5, label="Single"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor="white",
               markeredgecolor="#111111", markersize=9, label="ND-centered root"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="none",
               markeredgecolor="#111111", markersize=6, label="highest DPT object"),
    ]
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.49, 0.025),
               frameon=False, ncol=5, fontsize=8)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(args.out.with_suffix(".svg"), facecolor="white")
    plt.close(fig)
    print(f"wrote donor trajectory figure for {len(donors)} donors -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
