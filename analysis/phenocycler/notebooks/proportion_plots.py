"""Plot absolute cell-type percentages with a shared break in an unoccupied interval."""
import math
import textwrap

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
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


def _axes(n_categories, windows, height=4.2):
    broken = len(windows) == 2
    fig, axes = plt.subplots(len(windows), 1, figsize=(max(7, .57*n_categories + 2.4), height),
                             sharex=True, squeeze=False,
                             gridspec_kw={"height_ratios":[1, 3.2] if broken else [1], "hspace":.1})
    axes = list(axes[:,0])
    for ax, limits in zip(axes, reversed(windows)):
        ax.set_ylim(*limits)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_axisbelow(True)
        ax.grid(axis='y', linewidth=.45, color='#dfe2e3')
        locator = MaxNLocator(nbins=2 if broken and ax is axes[0] else 5, min_n_ticks=2, steps=[1, 2, 2.5, 5, 10])
        ticks = locator.tick_values(*limits)
        ticks = ticks[(ticks >= limits[0]) & (ticks <= limits[1])]
        if broken:
            endpoint = limits[0] if ax is axes[0] else limits[1]
            ticks = np.unique(np.append(ticks[np.abs(ticks - endpoint) > (.24 if ax is axes[0] else .10)*np.ptp(limits)], endpoint))
        ax.set_yticks(ticks)
        ax.tick_params(labelsize=8)
    if broken:
        axes[0].spines['bottom'].set_visible(False)
        axes[0].tick_params(axis='x', bottom=False)
        for ax, y in [(axes[0], 0), (axes[1], 1)]:
            ax.plot(0, y, transform=ax.transAxes, marker=[(-1, -.6), (1, .6)], markersize=8,
                    linestyle='none', color='#343b40', markeredgewidth=1.0, clip_on=False)
    return fig, axes


def plot_percentages(summary, *, colors, windows, title, format_label=str,
                     groups=None, observations=None, group_colors=None):
    """Plot one bar per series/category. Optional individual donors use the same absolute axis.

    summary has cell types as its index, and tissue/donor or disease groups as columns.
    observations has one row per donor, an index level 'disease', and the same cell-type columns.
    """
    if not summary.index.is_unique or not summary.columns.is_unique:
        raise ValueError("Duplicate percentage-series labels")
    groups = list(summary.columns if groups is None else groups)
    labels = list(summary.index)
    values = [summary[groups].to_numpy(dtype=float).ravel()]
    if observations is not None:
        values.append(observations.loc[:, labels].to_numpy(dtype=float).ravel())
    finite = np.concatenate(values)
    finite = finite[np.isfinite(finite)]
    covered = np.zeros(len(finite), bool)
    for low, high in windows:
        covered |= (finite >= low-1e-9) & (finite <= high+1e-9)
    if not covered.all():
        raise ValueError("A proportion would fall in the axis break or outside the plot")
    fig, axes = _axes(len(labels), windows, height=5.4)
    xs = np.arange(len(labels))
    width = .78/len(groups)
    missing_entries = False
    zero_entries = False
    for i, group in enumerate(groups):
        mid = xs + (i - (len(groups)-1)/2) * width
        vals = summary[group].to_numpy(dtype=float)
        valid = np.isfinite(vals)
        fill = [group_colors[group] if group_colors else colors[label] for label in labels]
        for ax in axes:
            ax.bar(mid[valid], vals[valid], width=.90*width,
                   color=np.asarray(fill)[valid], edgecolor='none', alpha=.72 if group_colors else 1)
        zero = valid & (vals == 0)
        if zero.any():
            axes[-1].scatter(mid[zero], np.zeros(zero.sum()), marker='_', s=16, c='#1f292b',
                             zorder=5, clip_on=False)
            zero_entries = True
        if (~valid).any():
            axes[-1].scatter(mid[~valid], np.zeros((~valid).sum()), marker='x', s=16,
                             c='#6c7274', linewidth=.9, zorder=5, clip_on=False)
            missing_entries = True
        if observations is not None:
            cohort = observations.xs(group, level='disease').sort_index()
            offsets = np.linspace(-.27, .27, len(cohort)) * width if len(cohort)>1 else np.zeros(len(cohort))
            for j,label in enumerate(labels):
                measured = cohort[label].notna()
                for ax in axes:
                    ax.scatter(mid[j]+offsets[measured], cohort.loc[measured, label], s=10, linewidth=.25,
                               edgecolor='white', color=group_colors[group], zorder=4, clip_on=True)
    for ax in axes:
        ax.set_xlim(-.62, len(labels)-.38)
    axes[-1].set_xticks(xs, [textwrap.fill(format_label(label), width=26, break_long_words=False)
                             for label in labels], rotation=48, ha='right', rotation_mode='anchor')
    axes[0].set_title(title, loc='left', fontsize=11, fontweight='bold', pad=8)
    fig.supylabel('Admitted cells (%)', x=.012, fontsize=10)
    legend = []
    if group_colors:
        legend.extend(Line2D([],[],linestyle='none',marker='s',color=group_colors[g],label=g) for g in groups)
    if zero_entries:
        legend.append(Line2D([],[],linestyle='none',marker='_',color='#1f292b',label='Measured zero'))
    if missing_entries:
        legend.append(Line2D([],[],linestyle='none',marker='x',color='#6c7274',label='Unavailable'))
    if legend:
        fig.legend(handles=legend, loc='upper right', bbox_to_anchor=(.99,.998), ncol=len(legend),
                   frameon=False, fontsize=8, handletextpad=.35, columnspacing=1.)
    fig.subplots_adjust(left=.075, right=.99, top=.86, bottom=.31)
    return fig
