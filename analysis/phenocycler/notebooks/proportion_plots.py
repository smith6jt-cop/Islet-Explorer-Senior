"""Plot full-section percentages in compact facets with explicitly labelled percent axes."""
import math
import textwrap

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np


def percent_windows(values):
    """Return one continuous range or two ranges separated by an endpoint-free gap."""
    raw = np.asarray(values, dtype=float).ravel()
    if np.isinf(raw).any():
        raise ValueError("A cell proportion cannot be infinite")
    points = np.unique(raw[np.isfinite(raw)])
    if points.size and (points.min() < 0 or points.max() > 100 + 1e-9):
        raise ValueError("Percentages must lie between 0 and 100")
    if not points.size or points.max() == 0:
        return ((0, 1),)
    maximum = min(100., points.max())
    top = min(100., max(maximum + 0.04 * maximum, maximum + .001))
    if maximum >= 40:
        top = float(min(100, math.ceil(top / 5) * 5))
    if points.size < 2:
        return ((0, top),)
    # Keep the base of the scale and select the widest truly empty interval.
    # Restrict breaks to a substantial interval separating a dominant population.
    candidates = [(b - a, a, b) for a, b in zip(points[:-1], points[1:])
                  if a > 0 and b >= 40 and b-a >= max(20, .33 * maximum)]
    if not candidates:
        return ((0, top),)
    gap, a, b = max(candidates)
    step = 5 if gap >= 35 else 2
    low = math.ceil((a + max(.15, .02 * gap)) / step) * step
    high = math.floor((b - max(.15, .02 * gap)) / step) * step
    if high - low < .22 * maximum:
        return ((0, top),)
    return ((0, float(low)), (float(high), top))


def _percent_ticks(maximum):
    """Choose a short linear percentage scale, including for sub-percent values."""
    locator = MaxNLocator(nbins=4, min_n_ticks=2, steps=[1, 2, 2.5, 5, 10])
    ticks = locator.tick_values(0, min(100, maximum * 1.06))
    return ticks[(ticks >= 0) & (ticks <= 100)]


def _pack_facets(bands):
    """Keep a scale local to the types that use it; balance large bands across rows."""
    rows = []
    for band, labels in bands:
        pages = math.ceil(len(labels) / 11)
        for page in np.array_split(np.arange(len(labels)), pages):
            block = [labels[i] for i in page]
            width = max(2.5, len(block) + 1)
            if not rows or sum(item[2] for item in rows[-1]) + width > 15:
                rows.append([])
            rows[-1].append((band, block, width))
    return rows


def _bar_axes(fig, spec, broken):
    if not broken:
        return [fig.add_subplot(spec)]
    nested = spec.subgridspec(2, 1, height_ratios=[3, 1], hspace=.05)
    upper = fig.add_subplot(nested[0])
    lower = fig.add_subplot(nested[1], sharex=upper)
    upper.spines["bottom"].set_visible(False)
    upper.tick_params(axis="x", bottom=False, labelbottom=False)
    for axis, y in [(upper, 0), (lower, 1)]:
        axis.plot(0, y, transform=axis.transAxes, marker=[(-1,-.6),(1,.6)], markersize=7,
                  linestyle="none", color="#343b40", markeredgewidth=.95, clip_on=False)
    return [upper, lower]


def _status_lines(summary, labels, groups, format_label, char_limit):
    """List unplotted zero/unavailable calls with their explicit series support."""
    by_state = {}
    for label in labels:
        zero = [str(g) for g in groups if summary.at[label, g] == 0]
        if len(zero) == len(groups):
            state = "Measured zero"
        elif not zero:
            state = "Unavailable"
        else:
            missing = [str(g) for g in groups if not np.isfinite(summary.at[label, g])]
            state = f"Zero ({', '.join(zero)}), unavailable ({', '.join(missing)})"
        by_state.setdefault(state, []).append(format_label(label))
    lines = []
    for state, names in by_state.items():
        line = state + "  ·  " + ", ".join(names)
        lines.extend(textwrap.wrap(line, width=char_limit, break_long_words=False))
    return lines


def plot_percentages(summary, *, colors, windows, title, format_label=str,
                     groups=None, observations=None, group_colors=None):
    """Plot absolute percentages, grouped by observed frequency, in one figure.

    Summary rows are cell types and columns are series. Donor observations
    use a disease index level and one column per type. Abundance panels are
    based on the highest displayed percentage for each type.
    """
    if not summary.index.is_unique or not summary.columns.is_unique:
        raise ValueError("Duplicate percentage-series labels")
    groups = list(summary.columns if groups is None else groups)
    labels = list(summary.index)
    if not groups:
        raise ValueError("A proportion figure needs at least one series")
    series = summary.loc[:,groups].to_numpy(dtype=float).T
    arrays = [series]
    if observations is not None:
        displayed = observations.loc[observations.index.get_level_values("disease").isin(groups), labels]
        arrays.append(displayed.to_numpy(dtype=float))
    data = np.vstack(arrays)
    finite = np.isfinite(data)
    if np.isinf(data).any() or (data[finite] < 0).any() or (data[finite] > 100+1e-9).any():
        raise ValueError("A percentage must be between 0 and 100 or missing")
    covered = np.zeros(data.shape, bool)
    for low, high in windows:
        covered |= (data >= low-1e-9) & (data <= high+1e-9)
    if (finite & ~covered).any():
        raise ValueError("A proportion falls in the supplied axis break or outside the plot")
    maxima = np.max(np.where(finite, data, 0), axis=0)
    bands = {}
    inactive = []
    for label, peak in zip(labels, maxima):
        if peak == 0:
            inactive.append(label)
            continue
        # Keep the discontinuity local to bars that actually reach its upper segment.
        band = math.inf if len(windows)==2 and peak >= windows[1][0] else math.ceil(math.log10(peak))
        bands.setdefault(band, []).append(label)
    rows = _pack_facets(sorted(bands.items(), reverse=True))
    widths = [sum(item[2] for item in row) for row in rows]
    width = min(13.2, max(4.1, .71*max(widths, default=3) + .6))
    title_lines = textwrap.wrap(title, max(36, int(width*10.7)), break_long_words=False)
    notes = _status_lines(summary, inactive, groups, format_label, max(45,int(width*15)))
    active = np.asarray([label not in inactive for label in labels])
    active_data = series[:,active]
    show_zero = np.any(active_data == 0)
    show_missing = np.isnan(active_data).any()
    legend = [Line2D([],[],linestyle="none",marker="s",color=group_colors[g],label=g)
              for g in groups] if group_colors else []
    if show_zero:
        legend.append(Line2D([],[],linestyle="none",marker="_",color="#1f292b",label="Measured zero"))
    if show_missing:
        legend.append(Line2D([],[],linestyle="none",marker="x",color="#6c7274",label="Unavailable"))
    legend_cols = min(len(legend), max(2, int(width*1.25))) if legend else 1
    legend_rows = math.ceil(len(legend)/legend_cols)
    header_height = .22*len(title_lines) + .20*legend_rows + .06
    # Reserve the header in inches so tall multi-panel figures don't acquire a
    # larger title gap. The layout engine handles ticks and rotated labels.
    data_height = 1.50
    notes_height = .18*len(notes) if notes else 0
    ratios = [data_height]*len(rows) + ([notes_height] if notes else [])
    if not ratios:
        raise ValueError("A proportion figure needs at least one cell type")
    figure_height = max(1.0, header_height + (data_height+1.05)*len(rows) + notes_height + .20)
    fig = plt.figure(figsize=(width,figure_height), layout="constrained")
    fig.get_layout_engine().set(w_pad=.06, h_pad=.04, wspace=.03, hspace=.025,
                                rect=(0, 0, 1, 1-header_height/figure_height))
    grid = fig.add_gridspec(len(ratios),1,height_ratios=ratios)
    header_x = .06 / width
    fig.text(header_x,1-.035/figure_height,"\n".join(title_lines),va="top",weight="bold",fontsize=11)
    if legend:
        fig.legend(handles=legend,loc="upper left",
                   bbox_to_anchor=(header_x,1-(.04+.22*len(title_lines))/figure_height),borderaxespad=0,
                   frameon=False,ncol=legend_cols,fontsize=8,handletextpad=.25,columnspacing=.9)
    maxima_by_label = dict(zip(labels, maxima))
    for row_index, row in enumerate(rows):
        slots = grid[row_index].subgridspec(1,len(row),width_ratios=[item[2] for item in row],wspace=.07)
        for col, (band, block, _) in enumerate(row):
            broken = math.isinf(band)
            axes = _bar_axes(fig, slots[col], broken)
            peak = max(maxima_by_label[label] for label in bands[band])
            if broken:
                top_lo, top_hi = windows[1]
                axes[0].set_ylim(top_lo,top_hi)
                axes[0].set_yticks([top_lo, (top_lo+top_hi)/2, top_hi])
                axes[1].set_ylim(*windows[0]); axes[1].set_yticks(windows[0])
            else:
                ticks = _percent_ticks(peak)
                axes[0].set_ylim(0, max(peak,ticks[-1])); axes[0].set_yticks(ticks)
            positions = np.arange(len(block))
            bar_width = .82/len(groups)
            for axis in axes:
                axis.spines[["top","right"]].set_visible(False)
                axis.set_axisbelow(True);axis.grid(False)
                axis.grid(axis="y",linewidth=.45,color="#dfe2e3")
                axis.tick_params(axis="both",labelsize=8,pad=2)
                axis.yaxis.set_major_formatter(FuncFormatter(lambda v,_: f"{v:g}"))
                axis.set_xlim(-.55, len(block)-.45)
            for i,g in enumerate(groups):
                mids = positions + (i-(len(groups)-1)/2)*bar_width
                heights = summary.loc[block,g].to_numpy(dtype=float)
                bar_colors = [group_colors[g] if group_colors else colors[label] for label in block]
                valid = np.isfinite(heights)
                for axis in axes:
                    axis.bar(mids[valid],heights[valid],width=.91*bar_width,
                             color=np.asarray(bar_colors)[valid],edgecolor="none",alpha=.75 if group_colors else 1)
                for mask,marker in [(valid & (heights==0),"_"),(~valid,"x")]:
                    if mask.any():
                        axes[-1].scatter(mids[mask],np.zeros(mask.sum()),marker=marker,s=18,c="#44494c",
                                         linewidths=.9,zorder=6,clip_on=False)
                if observations is not None:
                    donor_values = displayed.loc[displayed.index.get_level_values("disease")==g,block].sort_index()
                    jitter = (np.linspace(-.29,.29,len(donor_values)) if len(donor_values)>1 else np.zeros(len(donor_values))) * bar_width
                    for j,label in enumerate(block):
                        values = donor_values[label].to_numpy(dtype=float)
                        for axis in axes:
                            lo, hi = axis.get_ylim()
                            measured = np.isfinite(values) & (values>=lo) & (values<=hi)
                            axis.scatter(mids[j]+jitter[measured],values[measured],s=10,linewidth=.28,
                                         edgecolor="white",color=bar_colors[j],zorder=4,clip_on=False)
            rotation = 40 if len(block)>1 else 0
            axes[-1].set_xticks(positions,[textwrap.fill(format_label(label), width=22, break_long_words=False)
                                           for label in block],rotation=rotation,ha="right" if rotation else "center",
                                rotation_mode="anchor")
    if rows:
        fig.supylabel("Admitted cells (%)",fontsize=9)
    if notes:
        ax=fig.add_subplot(grid[-1]);ax.axis("off")
        ax.text(0,1,"\n".join(notes),transform=ax.transAxes,va="top",fontsize=8,linespacing=1.3)
    return fig
