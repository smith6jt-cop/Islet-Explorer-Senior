"""Display each donor's fitted islet graph and its measured inflammation."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import PowerNorm
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

DENSITY = "combined_immune_per_mm2"
GROUP_ORDER = ["ND", "AAb+", "T1D"]


def load_trajectories(data_path):
    path = Path(data_path)
    points = pd.read_parquet(path / "donor_islet_trajectory_points.parquet")
    nodes = pd.read_parquet(path / "donor_islet_trajectory_nodes.parquet")
    edges = pd.read_parquet(path / "donor_islet_trajectory_edges.parquet")
    measured = pd.read_parquet(path / "islet_inflammation_ranks.parquet")
    measured.index.name = "islet_id"
    for table in (points, nodes, edges, measured):
        table["donor"] = table.donor.astype(str)
    points = points.merge(measured.reset_index()[["donor", "islet_id", "core_immune_per_mm2",
                                                  "peri_immune_per_mm2"]], on=["donor", "islet_id"],
                          validate="one_to_one")
    if len(points) != len(measured) or points.islet_id.duplicated().any():
        raise ValueError("The trajectory and measured islet rosters differ")
    density_max = measured[DENSITY].max()
    return points, nodes, edges, PowerNorm(gamma=.5, vmin=0, vmax=density_max)


def _graph(ax, inputs, donor, compact=False):
    points, nodes, edges, norm = inputs
    data = points.loc[points.donor.eq(donor)]
    fitted = data.dropna(subset=["graph_x", "graph_y", "graph_dpt"]).sort_values(DENSITY)
    if fitted.empty:
        ax.text(.5, .5, "No fitted graph", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return data, fitted
    centers = nodes.loc[nodes.donor.eq(donor)].copy()
    centers.index = centers.graph_cluster.astype(str)
    lines = edges.loc[edges.donor.eq(donor) & edges.in_tree]
    if len(lines):
        xy = centers[["graph_x","graph_y"]]
        segments = np.stack([xy.loc[lines.source.astype(str)].to_numpy(),
                             xy.loc[lines.target.astype(str)].to_numpy()], axis=1)
        ax.add_collection(LineCollection(segments, linewidths=1.1 if compact else 2.4,
                                         color="#617582", alpha=.8, zorder=2))
    ax.scatter(fitted.graph_x, fitted.graph_y, c=fitted[DENSITY], cmap="viridis", norm=norm,
               s=6 if compact else 23, edgecolors="none", alpha=.87, zorder=1, rasterized=True)
    root = fitted.loc[fitted.graph_root]
    if len(root) != 1:
        raise ValueError(f"Expected one root islet for donor {donor}")
    ax.scatter(root.graph_x, root.graph_y, marker="*", s=72 if compact else 155, facecolor="#fff5df",
               edgecolor="#251b19", lw=.9, zorder=3)
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(.13)
    ax.set_axis_off()
    return data, fitted


def plot_donor_trajectory(inputs, donor_id):
    """Plot graph topology and immune density along this donor's diffusion pseudotime."""
    donor = str(donor_id)
    points, _, _, norm = inputs
    if donor not in set(points.donor):
        raise KeyError(f"No islet measurements for donor {donor}")
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9), gridspec_kw={"width_ratios":[1.05,1.25,1.25]},
                             layout="constrained")
    data, fitted = _graph(axes[0], inputs, donor)
    axes[0].set_title("Islet composition graph")
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=axes[0],
                       orientation="horizontal", shrink=.8, pad=.05, aspect=38)
    ticks = [t for t in (0,1000,5000,10000) if t <= norm.vmax]
    bar.set_ticks(ticks, labels=[f"{t:,}" for t in ticks])
    bar.set_label("Core + peri-islet immune cells/mm²", fontsize=9)
    fields = [("core_immune_per_mm2", "Core", "#25627a"), ("peri_immune_per_mm2","20 µm rim","#ab553d")]
    ymax=points[["core_immune_per_mm2","peri_immune_per_mm2"]].max().max()
    for ax, ordering, xfield, xlabel, title in [
            (axes[1], data, DENSITY, "Combined immune cells/mm²", "Measured inflammation trajectory"),
            (axes[2], fitted, "graph_dpt", "Diffusion pseudotime within donor", "Inflammation along graph pseudotime")]:
        for field, name, ink in fields:
            observed=ordering.dropna(subset=[field,xfield]).sort_values(xfield)
            ax.scatter(observed[xfield],observed[field],s=8,alpha=.19,color=ink,lw=0,rasterized=True)
            # Ten bins on the common horizontal scale; gaps show unsupported intervals.
            if xfield==DENSITY:
                bounds=np.linspace(0,np.sqrt(norm.vmax),11)**2
                bounds[-1]=norm.vmax
            else:
                bounds=np.linspace(0,1,11)
            bins=pd.cut(observed[xfield],bounds,include_lowest=True)
            by_bin=observed.groupby(bins,observed=False)[[xfield,field]].median()
            ax.plot(by_bin[xfield],by_bin[field],marker="o",ms=3.4,lw=1.6,color=ink,label=name)
        transforms=(lambda x:np.sign(x)*np.sqrt(np.abs(x)),lambda x:np.sign(x)*np.square(x))
        ax.set_yscale("function", functions=transforms)
        ax.set_ylim(-ymax*.0025,ymax*1.04)
        ax.set_yticks([t for t in (0,1000,5000,10000,25000,50000) if t<=ymax])
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_:f"{v:,.0f}"))
        if xfield==DENSITY:
            ax.set_xscale("function",functions=transforms)
            ax.set_xlim(-norm.vmax*.001,norm.vmax*1.04)
            ax.set_xticks(ticks,labels=[f"{t:,}" for t in ticks])
        else:
            ax.set_xlim(-.025,1.025)
        ax.set(xlabel=xlabel,ylabel="Core or rim immune cells/mm²",title=title)
        ax.grid(axis="y",alpha=.19)
    axes[1].legend(frameon=False)
    fig.suptitle(f"Donor {donor}  ·  {data.disease.iloc[0]}  ·  {len(fitted):,} of {len(data):,} islets fitted")
    fig.supxlabel("Star = low-inflammation root     ·     Graph lines = PAGA tree     ·     Curves = median per interval",
                  fontsize=9, color="#374852")
    return fig


def plot_trajectory_atlas(inputs):
    points, _, _, norm = inputs
    roster = points[["disease","donor"]].drop_duplicates().assign(
        order=lambda d:d.disease.map({g:i for i,g in enumerate(GROUP_ORDER)}),
        donor_number=lambda d:d.donor.astype(int)).sort_values(["order","donor_number"])
    fig, axes = plt.subplots(int(np.ceil(len(roster)/5)),5,figsize=(13.8,12.8),squeeze=False,layout="constrained")
    for ax, row in zip(axes.flat,roster.itertuples()):
        data,fitted=_graph(ax,inputs,row.donor,compact=True)
        ax.set_title(f"{row.donor} · {row.disease}   ({len(fitted):,}/{len(data):,})",fontsize=9)
    for ax in list(axes.flat)[len(roster):]:
        ax.set_visible(False)
    cbar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap="viridis"),ax=axes,orientation="horizontal",
                      shrink=.46,pad=.02,aspect=52)
    ticks=[t for t in (0,1000,5000,10000) if t<=norm.vmax]
    cbar.set_ticks(ticks,labels=[f"{v:,}" for v in ticks])
    cbar.set_label("Combined immune density (cells/mm²)")
    fig.suptitle("Independently fitted donor islet trajectories\nStar = low-inflammation root  ·  Lines = PAGA tree  ·  Fitted / native islets",fontsize=14)
    return fig
