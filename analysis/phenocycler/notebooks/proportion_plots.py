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


def _band(peak):
    """Return a reproducible frequency range and file suffix for a positive percentage."""
    if peak >= 10:
        return 1, "above_10pct", "≥10%"
    if peak >= 1:
        return 0, "common", "1% to <10%"
    power = math.floor(math.log10(peak)) + 1
    upper = 10.0 ** power
    slug = np.format_float_positional(upper, trim="-").replace(".", "p")
    return power - 1, f"below_{slug}pct", f"{upper/10:g}% to <{upper:g}%"


def _nice_axis(peak):
    """Round just above a measured maximum, even when it is below 0.001%."""
    locator = MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10])
    ticks = locator.tick_values(0, min(100., peak * 1.035))
    ticks = ticks[(ticks >= 0) & (ticks <= 100)]
    if ticks[-1] < peak:
        ticks = np.append(ticks, min(100., peak))
    # A full extra interval can dominate a small panel if the maximum is
    # already a round value. Keep the whitespace to at most 10%.
    top = min(ticks[-1], peak * 1.10, 100.)
    return top, ticks[ticks <= top + np.spacing(top)]


def _frequency_pages(summary, groups, observations, windows):
    if not summary.index.is_unique or not summary.columns.is_unique:
        raise ValueError("Duplicate percentage-series labels")
    if not groups or summary.empty:
        raise ValueError("A proportion figure needs at least one cell type and series")
    summary = summary.loc[:, groups]
    values = summary.to_numpy(dtype=float).T
    if observations is not None:
        mask = observations.index.get_level_values("disease").isin(groups)
        observations = observations.loc[mask, summary.index]
        values = np.vstack((values, observations.to_numpy(dtype=float)))
    finite = np.isfinite(values)
    if np.isinf(values).any() or (values[finite] < 0).any() or (values[finite] > 100 + 1e-9).any():
        raise ValueError("A percentage must be between 0 and 100 or missing")
    covered = np.zeros(values.shape, dtype=bool)
    for low, high in windows:
        covered |= (values >= low - 1e-9) & (values <= high + 1e-9)
    if (finite & ~covered).any():
        raise ValueError("A proportion falls in the supplied axis break or outside the plot")

    peaks = np.max(np.where(finite, values, 0), axis=0)
    medians = summary.max(axis=1)
    bins = {}
    for label, peak in zip(summary.index, peaks):
        typical = medians[label]
        if peak == 0:
            key = (2, 0, "zero_or_unavailable", "Measured zero or unavailable")
        elif np.isfinite(typical) and typical > 0:
            order, slug, size = _band(typical)
            if len(windows) == 2 and typical >= windows[1][0]:
                order, slug, size = 2, "dominant", f"≥{windows[1][0]:g}%"
            description = ("Largest disease median " if observations is not None else "Frequency ") + size
            key = (0, -order, slug, description)
        else:
            # Zero disease medians can coexist with a substantial individual donor.
            # Separate these by donor frequency so each positive observation is visible.
            order, slug, size = _band(peak)
            status = "Median 0%" if np.isfinite(typical) else "Median unavailable"
            key = (1, -order, ("zero_median_" if np.isfinite(typical) else "no_median_") + slug,
                   f"{status} · highest individual donor {size}")
        bins.setdefault(key, []).append(label)
    pages = {}
    for (_, _, slug, description), types in sorted(bins.items()):
        n = 1 if slug == "zero_or_unavailable" else math.ceil(len(types) / 10)
        for i, indices in enumerate(np.array_split(np.arange(len(types)), n)):
            part = [types[j] for j in indices]
            suffix = f"_{i+1}" if n > 1 else ""
            if slug.startswith(("zero_median_", "no_median_")):
                status = "Median 0%" if slug.startswith("zero_median_") else "Median unavailable"
                description = f"{status} · individual donors up to {observations[part].max().max():.4g}%"
            subtitle = f"{description} · {i+1}/{n}" if n > 1 else description
            pages[slug + suffix] = (part, subtitle)
    return summary, observations, pages


def _axis_style(axis):
    axis.spines[["top", "right"]].set_visible(False)
    axis.set_axisbelow(True)
    axis.grid(False)
    axis.grid(axis="y", linewidth=.45, color="#dfe2e3")
    axis.tick_params(axis="both", labelsize=8, pad=2)
    axis.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def _draw_page(fig, *, summary, observations, groups, colors, group_colors, format_label,
               title, description, page, windows, width):
    labels = list(summary.index)
    title_lines = textwrap.wrap(title, width=int(width * 11.3), break_long_words=False)
    subtitle = textwrap.wrap(description, width=int(width * 10.6), break_long_words=False)
    fig.suptitle("\n".join(title_lines + subtitle), x=.015, ha="left", fontsize=10,
                 fontweight="bold")
    if page == "zero_or_unavailable":
        lines = _status_lines(summary, labels, groups, format_label, int(width * 14))
        axis = fig.add_subplot()
        axis.set_axis_off()
        axis.text(.008, .96, "\n".join(lines), transform=axis.transAxes, va="top", fontsize=8.6)
        return

    if len(labels) == 1 and len(groups) == 1 and observations is None and page != "dominant":
        label = labels[0]
        value = float(summary.iloc[0, 0])
        limit, ticks = _nice_axis(value)
        axis = fig.add_subplot()
        axis.barh(0, value, height=.44, color=group_colors[groups[0]] if group_colors else colors[label])
        axis.set(xlim=(0, limit), xticks=ticks, ylim=(-.5, .5), yticks=[])
        axis.spines[["left", "right", "top"]].set_visible(False)
        axis.grid(False)
        axis.tick_params(axis="x", labelsize=8, pad=1)
        axis.set_xlabel("Admitted cells (%)", fontsize=8, labelpad=1)
        axis.set_title(f"{format_label(label)}  ·  {value:.4g}%", fontsize=9, loc="left", pad=1)
        return

    raw = summary.to_numpy(dtype=float)
    donor_peak = 0 if observations is None else observations.max().max()
    positive_maximum = max(np.nanmax(np.nan_to_num(raw)), np.nan_to_num(donor_peak))
    peak = np.nanmax(np.nan_to_num(raw))
    bounded = len(windows) == 2 and peak >= windows[1][0]
    # A zero median may be supported by many zero observations plus a small,
    # tightly grouped set of positive donors. Show only the observed ranges.
    positive_donors = observations.to_numpy(dtype=float) if peak == 0 and observations is not None else np.array([])
    positive_donors = positive_donors[np.isfinite(positive_donors) & (positive_donors > 0)]
    zero_gap = positive_donors.size > 0 and positive_donors.min() >= .65 * positive_maximum
    if zero_gap:
        spread = positive_maximum - positive_donors.min()
        pad = max(.025 * positive_maximum, .2 * spread)
        locator = MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10])
        upper_low = max(0., positive_donors.min() - pad)
        upper_high = min(100., positive_maximum + pad)
        ticks = locator.tick_values(upper_low, upper_high)
        focus_top = min(.04 * positive_maximum, upper_low / 2)
        focus_ticks = [0]
        upper_label = "Positive individual donors"
        bounded = True
    elif bounded:
        focus_top, focus_ticks = windows[0][1], list(windows[0])
        # Give the common populations below the exocrine break readable ticks.
        if focus_top > 0:
            focus_ticks = np.linspace(0, focus_top, 3)
        upper_low, upper_high = windows[1]
        upper_label = None
    else:
        focus_top, focus_ticks = _nice_axis(peak or positive_maximum)
        upper_low, upper_high = focus_top, _nice_axis(positive_maximum)[0]
        upper_label = "Higher individual donors"
    has_upper = bounded or (peak > 0 and donor_peak > focus_top)
    handles = [Line2D([], [], linestyle="none", marker="s", color=group_colors[g], label=g)
               for g in groups] if group_colors else []
    if np.any(raw == 0):
        handles.append(Line2D([], [], linestyle="none", marker="_", color="#273038", label="Measured zero"))
    if np.isnan(raw).any():
        handles.append(Line2D([], [], linestyle="none", marker="x", color="#555b5d", label="Unavailable"))
    cols = min(len(handles), max(2, int(width * 1.1))) if handles else 1
    heights = ([.17 * math.ceil(len(handles)/cols)] if handles else []) + ([1.1 if bounded else .56] if has_upper else []) + [.38 if bounded else 1.18]
    grid = fig.add_gridspec(len(heights), 1, height_ratios=heights, hspace=.025)
    if handles:
        header = fig.add_subplot(grid[0])
        header.set_axis_off()
        header.legend(handles=handles, loc="center left", borderaxespad=0, ncol=cols,
                      fontsize=8, frameon=False, handletextpad=.25, columnspacing=.9)
    base = fig.add_subplot(grid[-1])
    base.set_ylim(0, focus_top)
    base.set_yticks(focus_ticks)
    if bounded:
        fig.supylabel("Admitted cells (%)", fontsize=8.5)
    else:
        base.set_ylabel("Admitted cells (%)", fontsize=8.5)
    axes = [base]
    if has_upper:
        upper = fig.add_subplot(grid[-2], sharex=base)
        upper.set_ylim(upper_low, upper_high)
        # Keep the small upper panel legible. Numeric ticks are independent
        # of the lower axis and include the highest individual data point.
        if zero_gap:
            upper.set_yticks(ticks[(ticks >= upper_low) & (ticks <= upper_high)])
        elif bounded:
            upper.set_yticks([upper_low, (upper_low + upper_high) / 2, upper_high])
        else:
            upper.yaxis.set_major_locator(MaxNLocator(nbins=2, min_n_ticks=2))
        upper.spines["bottom"].set_visible(False)
        upper.tick_params(axis="x", bottom=False, labelbottom=False)
        if upper_label:
            upper.set_title(upper_label, fontsize=8, loc="left", pad=1)
        for ax, edge in ([(upper, 0), (base, 1)] if bounded else []):
            ax.plot(0, edge, transform=ax.transAxes, marker=[(-1,-.6),(1,.6)], markersize=6,
                    linestyle="none", color="#343b40", markeredgewidth=.8, clip_on=False)
        axes.append(upper)
    positions = np.arange(len(labels))
    bar_width = .82 / len(groups)
    for ax in axes:
        _axis_style(ax)
        ax.set_xlim(-.54, len(labels) - .46)
    for i, g in enumerate(groups):
        center = positions + (i-(len(groups)-1)/2) * bar_width
        heights = summary[g].to_numpy(dtype=float)
        inks = np.asarray([group_colors[g] if group_colors else colors[label] for label in labels])
        valid = np.isfinite(heights)
        for ax in (axes if bounded else [base]):
            ax.bar(center[valid], heights[valid], width=.9*bar_width, color=inks[valid],
                   alpha=.72 if group_colors else 1, edgecolor="none")
        for mask, marker in [(valid & (heights == 0), "_"), (~valid, "x")]:
            if mask.any():
                base.scatter(center[mask], np.zeros(mask.sum()), marker=marker, s=18,
                             c="#343b40", linewidths=.9, zorder=6, clip_on=False)
        if observations is None:
            continue
        donors = observations.loc[observations.index.get_level_values("disease") == g].sort_index()
        jitter = (np.linspace(-.32, .32, len(donors)) if len(donors) > 1 else np.zeros(len(donors))) * bar_width
        for j,label in enumerate(labels):
            y = donors[label].to_numpy(dtype=float)
            for k, ax in enumerate(axes):
                lo, hi = ax.get_ylim()
                inside = np.isfinite(y) & (y >= lo) & (y <= hi)
                if k and not bounded:  # The two zoom ranges meet at one value.
                    inside &= y > focus_top
                ax.scatter(center[j] + jitter[inside], y[inside], s=13, linewidth=.32,
                           edgecolor="white", color=inks[j], zorder=4, clip_on=False)
    wrap = min(28, max(11, int((width - .55) * 12 / len(labels))))
    xticks = [textwrap.fill(format_label(label), width=wrap, break_long_words=False,
                            break_on_hyphens=True) for label in labels]
    base.set_xticks(positions, xticks, ha="center")


def _page_size(title, page, labels, description, summary, groups, observations, format_label):
    if page == "zero_or_unavailable":
        width = min(8., max(4.3, 1.1+max(len(title), len(description)) / 12.))
        notes = _status_lines(summary, labels, groups, format_label, int(width * 14))
        return width, .49 + .17 * (len(textwrap.wrap(title, int(width * 11.3))) + len(notes))
    if len(labels) == 1 and len(groups) == 1 and observations is None and page != "dominant":
        width = min(7.0, max(4.3, len(title) / 11.7))
        title_rows = len(textwrap.wrap(title, int(width * 11.3)))
        return width, 1.39 + .18 * title_rows
    width = min(13.3, max(4.1, .95 + len(labels) * (1.13 if observations is not None else .99)))
    title_rows = (len(textwrap.wrap(title, int(width * 11.3)))
                  + len(textwrap.wrap(description, int(width * 10.6))) - 1)
    wrap = min(28, max(11, int((width - .55) * 12 / len(labels))))
    label_rows = max(len(textwrap.wrap(format_label(label), width=wrap, break_long_words=False)) for label in labels)
    peak = np.nanmax(np.nan_to_num(summary.to_numpy(dtype=float)))
    donor_peak = observations.max().max() if observations is not None else 0
    has_upper = peak > 0 and donor_peak > _nice_axis(peak)[0]
    # Exocrine sections also require space for the upper measured segment.
    if peak >= 40:
        has_upper = True
    return width, 1.85 + .17 * (title_rows + label_rows) + .77 * has_upper + .14 * (observations is not None)


def plot_percentage_pages(summary, *, colors, windows, title, format_label=str,
                          groups=None, observations=None, group_colors=None):
    """Return separate labelled figures for observed frequencies and disease median ranges.

    Per-donor percentages use all admitted cells. For group summaries, the largest
    measured disease median determines membership in a figure. Donors exceeding
    the median plot's scale appear on the aligned upper axis at their actual value.
    """
    groups = list(summary.columns if groups is None else groups)
    summary, observations, pages = _frequency_pages(summary, groups, observations, windows)
    figures = {}
    for slug, (labels, description) in pages.items():
        selected = summary.loc[labels]
        donors = observations.loc[:, labels] if observations is not None else None
        size = _page_size(title, slug, labels, description, selected, groups, donors, format_label)
        fig = plt.figure(figsize=size, layout="constrained")
        fig.get_layout_engine().set(h_pad=.035, w_pad=.045, hspace=.025)
        _draw_page(fig, summary=selected, observations=donors, groups=groups, colors=colors,
                   group_colors=group_colors, format_label=format_label, title=title,
                   description=description, page=slug, windows=windows, width=size[0])
        figures[slug] = fig
    return figures


def plot_percentages(summary, *, colors, windows, title, format_label=str, groups=None,
                     observations=None, group_colors=None):
    """Return one figure containing all ranges for callers using the original API."""
    groups = list(summary.columns if groups is None else groups)
    summary, observations, pages = _frequency_pages(summary, groups, observations, windows)
    sizes = [_page_size(title, slug, labels, desc, summary.loc[labels], groups,
                        observations.loc[:, labels] if observations is not None else None, format_label)
             for slug, (labels, desc) in pages.items()]
    width = max(w for w, _ in sizes)
    heights = [h for _, h in sizes]
    fig = plt.figure(figsize=(width, sum(heights)), layout="constrained")
    sections = fig.subfigures(len(pages), 1, height_ratios=heights, squeeze=False)
    for section, (slug, (labels, description)) in zip(sections.flat, pages.items()):
        _draw_page(section, summary=summary.loc[labels], observations=(
            observations.loc[:, labels] if observations is not None else None),
            groups=groups, colors=colors, group_colors=group_colors, format_label=format_label,
            title=title, description=description, page=slug, windows=windows, width=width)
    return fig
