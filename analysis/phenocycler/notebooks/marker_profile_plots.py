"""Quantitative protein profiles for marker comparisons with limited heatmap contrast."""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np


def _plot_points(ax, values, positions, *, color, shape, upper):
    """Label each calibrated measurement and distinguish missing values from zero."""
    finite = np.isfinite(values)
    for x, y in zip(values[finite], positions[finite]):
        ax.scatter(x, y, marker=shape, s=27, color=color if x != 0 else "white",
                   edgecolor=color, linewidth=1.1, clip_on=False, zorder=3)
        near_right = x >= .86 * upper
        ax.annotate(f"{x:.3g}", (x, y), xytext=(-6 if near_right else 6, 0),
                    textcoords="offset points", va="center", ha="right" if near_right else "left",
                    color=color, fontsize=8)
    for y in positions[~finite]:
        ax.text(.015, y, "Unavailable", transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=8, color=color)


def plot_calibrated_profiles(mean, *, title, scope_note, reference_max=None):
    """Draw ordered, untransformed calibrated means from the same table as the heatmap.

    One marker is plotted against all row labels. Two rows can be compared on a
    single numeric axis across markers. A larger matrix retains its row and column
    labels with point area encoding the calibrated mean. An optional reference
    maximum sets the same numeric range across companion identity and QC views.
    """
    values = mean.to_numpy(dtype=float)
    finite = np.isfinite(values)
    if not mean.shape[0] or not mean.shape[1]:
        raise ValueError("A profile needs at least one row and marker")
    if finite.any() and ((values[finite] < 0).any() or (values[finite] > 1 + 1e-7).any()):
        raise ValueError("Profile means must be release-calibrated values on the 0–1 scale")

    single_marker = len(mean.columns) == 1
    paired_rows = len(mean.index) == 2
    if single_marker or paired_rows:
        labels = list(mean.index if single_marker else mean.columns)
        width = max(6.2, min(10, 4.35 + .057 * max(len(str(label)) for label in labels)))
        height = 1.75 + (.23 if single_marker else .52) * len(labels)
        if single_marker:
            fig, ax = plt.subplots(figsize=(width, max(3.2, height)), layout="constrained")
        else:
            fig, (ax, key) = plt.subplots(2, 1, figsize=(max(width, 7), max(3.4, height+.56)),
                                          gridspec_kw={"height_ratios":[max(1, .52*len(labels)), .45]},
                                          layout="constrained")
            key.set_axis_off()
        positions = np.arange(len(labels))
        peak = np.max(values[finite]) if finite.any() else 0.
        if reference_max is not None and np.isfinite(reference_max):
            if reference_max < 0 or reference_max > 1 + 1e-7:
                raise ValueError("The reference maximum must be a calibrated value on the 0–1 scale")
            peak = max(peak, reference_max)
        upper = min(1., 1.12 * peak) if peak > 0 else 1.
        ax.set_xlim(0, max(upper, peak))
        ax.set_ylim(len(labels)-.5, -.5)
        ax.set_yticks(positions, labels)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=2))
        ax.set_xlabel("Mean calibrated protein signal", fontsize=9)
        ax.tick_params(axis="both", labelsize=8)
        ax.set_axisbelow(True)
        ax.grid(False)
        ax.grid(axis="x", color="#dce2e5", linewidth=.55)
        ax.spines[["top", "right"]].set_visible(False)
        if single_marker:
            color = "#235d72"
            _plot_points(ax, values[:, 0], positions, color=color, shape="o", upper=upper)
            ax.set_title(f"{title}\n{mean.columns[0]} · {scope_note}", fontsize=10)
        else:
            handles = []
            for i, (color, shape, offset) in enumerate([("#245673", "o", -.19), ("#9d4824", "^", .19)]):
                _plot_points(ax, values[i], positions + offset, color=color, shape=shape, upper=upper)
                handles.append(Line2D([], [], marker=shape, color=color, linestyle="none",
                                      markersize=5, label=str(mean.index[i])))
            key.legend(handles=handles, loc="center left", frameon=False, fontsize=8)
            ax.set_title(f"{title}\n{scope_note}", fontsize=10)
    else:
        nrows, ncols = values.shape
        fig, ax = plt.subplots(figsize=(max(6.8,.34*ncols+3.6), max(3.2,.28*nrows+1.9)),
                               layout="constrained")
        yy, xx = np.indices(values.shape)
        positive = finite & (values > 0)
        measured_zero = finite & (values == 0)
        ax.scatter(xx[positive], yy[positive], s=125*values[positive], c="#245e74", linewidths=0)
        ax.scatter(xx[measured_zero], yy[measured_zero], s=11, facecolors="none", edgecolors="#245e74", linewidths=.8)
        ax.scatter(xx[~finite], yy[~finite], s=12, c="#858b91", marker="x", linewidths=.7)
        ax.set_xticks(np.arange(ncols), mean.columns, rotation=90)
        ax.set_yticks(np.arange(nrows), mean.index)
        ax.set_xlim(-.5,ncols-.5)
        ax.set_ylim(nrows-.5,-.5)
        ax.tick_params(labelsize=8)
        ax.grid(False)
        handles = [ax.scatter([],[], s=125*x, color="#245e74", label=f"{x:g}") for x in [.1,.5,1.]]
        if measured_zero.any():
            handles.append(Line2D([],[],marker="o",mec="#245e74",mfc="none",ms=4,ls="none",label="Zero"))
        if (~finite).any():
            handles.append(Line2D([],[],marker="x",color="#858b91",ms=4,ls="none",label="Unavailable"))
        ax.legend(handles=handles, title="Mean calibrated signal", loc="center left",
                  bbox_to_anchor=(1, .5), fontsize=8, title_fontsize=8, frameon=False)
        ax.set_title(f"{title}\n{scope_note}",fontsize=10)
    plt.show()
    return fig
