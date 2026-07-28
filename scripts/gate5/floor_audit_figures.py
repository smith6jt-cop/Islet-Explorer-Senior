#!/usr/bin/env python
"""Review figures for the raw input floor audit — one page per marker.

The decision surface for `audit_marker_floors.py`. Each page answers one question about ONE marker:
does its production floor sit at a real shoulder in every donor's own distribution, or does it land in
background on some of them?

Three panels per page:

  LEFT    every donor's own log-intensity histogram, stacked as a ridgeline, with that donor's
          PRODUCTION floor drawn on it. This is the working rule's "validate against that donor's own
          histogram" made routine. A floor in background is visible as a line sitting on the left-hand
          noise peak rather than between two modes.
  MIDDLE  positive fraction by donor, one line per candidate method, production method emphasised.
          The evidence a floor has failed is the marker's OWN donor-to-donor instability -- a marker
          cannot be 3.6% positive on one donor and 90.3% on another.
  RIGHT   separation (Otsu's criterion at the chosen threshold, between-class / total variance) by
          donor and method. Frequency-free: it answers "is there a real split here at all", without
          any assumption about how abundant the marker should be.

No percentile black points and no gamma lift anywhere -- histograms are drawn on the marker's own
measured values, per `review-figures-must-show-background-honestly`. Every page writes its source rows
to `figures/<marker>.csv` so no claim here rests on reading a colour off a line.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")                     # headless; the default tkagg blocks on this machine
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages   # noqa: E402

# Paul Tol bright — colourblind-safe, consistent with the app's DONOR_COLORS convention.
METHOD_COLORS = {
    "linear_otsu": "#EE6677", "linear_triangle": "#4477AA", "linear_yen": "#228833",
    "log_otsu": "#CCBB44", "log_triangle": "#66CCEE", "log_yen": "#AA3377",
    "log_li": "#BBBBBB", "log_gmm_crossover": "#000000",
}
PRODUCTION_LW, OTHER_LW = 2.6, 1.2


def _ridgeline(ax, hist: pd.DataFrame, marker: str, prod: pd.DataFrame, sep_median: float) -> None:
    """Each donor's own distribution with its own production floor drawn on it."""
    donors = sorted(hist.donor.unique())
    if not donors:
        ax.text(0.5, 0.5, "no histogram data", ha="center", va="center", transform=ax.transAxes)
        return
    for i, donor in enumerate(donors):
        h = hist[hist.donor == donor]
        y = h["density"].to_numpy(float)
        y = y / y.max() * 0.9 if y.max() > 0 else y
        ax.fill_between(h["log_center"], i, i + y, color="#AAAAAA", alpha=0.55, linewidth=0)
        ax.plot(h["log_center"], i + y, color="#555555", lw=0.6)
        f = prod[prod.donor == donor]
        if len(f) and np.isfinite(f["threshold"].iloc[0]):
            x = np.log10(max(f["threshold"].iloc[0], 1e-3))
            # Red where THIS donor's floor finds a much weaker split than the same marker's floor finds
            # on its typical donor. A marker compared against itself -- not against a frequency, and not
            # against other markers. An absolute cutoff cannot do this job: it flags EpCAM, which is
            # stable on all 20 donors, and misses 6566 INS at 23% positive.
            weak = f["separation"].iloc[0] < 0.5 * sep_median
            ax.plot([x, x], [i, i + 1.0], color="#EE6677" if weak else "#228833",
                    lw=2.0, solid_capstyle="butt")
    ax.set_yticks(np.arange(len(donors)) + 0.4)
    ax.set_yticklabels(donors, fontsize=8)
    ax.set_xlabel("log10 raw whole-cell mean (pre-REDSEA)", fontsize=10)
    ax.set_title(f"{marker} — each donor's own distribution with its PRODUCTION floor\n"
                 f"(red = this donor's split is under half the strength of this marker's typical donor)",
                 fontsize=11)
    ax.set_ylim(-0.2, len(donors) + 0.2)


def _by_donor(ax, d: pd.DataFrame, column: str, prod_method: str, title: str, ylabel: str,
              logy: bool = False) -> None:
    donors = sorted(d.donor.unique())
    x = np.arange(len(donors))
    for method, g in d.groupby("method"):
        g = g.set_index("donor").reindex(donors)
        is_prod = method == prod_method
        ax.plot(x, g[column].to_numpy(float), marker="o", ms=3.5,
                lw=PRODUCTION_LW if is_prod else OTHER_LW,
                color=METHOD_COLORS.get(method, "#888888"),
                alpha=1.0 if is_prod else 0.75, zorder=5 if is_prod else 2,
                label=f"{method}  ← PRODUCTION" if is_prod else method)
    ax.set_xticks(x)
    ax.set_xticklabels(donors, rotation=90, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11)
    if logy:
        ax.set_yscale("log")
    ax.grid(alpha=0.25, linewidth=0.5)


def marker_page(marker: str, per_donor: pd.DataFrame, hist: pd.DataFrame, summary: pd.DataFrame):
    d = per_donor[per_donor.marker == marker]
    prod_method = d["production_method"].iloc[0] or "(not a gate marker)"
    prod = d[d.method == prod_method] if prod_method in set(d.method) else d.iloc[:0]

    sep_median = float(prod["separation"].median()) if len(prod) else float("nan")
    fig, axes = plt.subplots(1, 3, figsize=(22, 9),
                             gridspec_kw={"width_ratios": [1.15, 1.0, 1.0]})
    _ridgeline(axes[0], hist[hist.marker == marker], marker, prod, sep_median)
    _by_donor(axes[1], d, "positive_fraction", prod_method,
              "Positive fraction by donor and method", "fraction of cells above floor")
    axes[1].axhline(0.5, color="#EE6677", ls=":", lw=1.2)
    _by_donor(axes[2], d, "separation", prod_method,
              "Separation — between-class / total variance\n(frequency-free: is there a real split?)",
              "separation")
    axes[2].set_ylim(0, 1)
    if np.isfinite(sep_median):
        # The line the red/green call in the left panel is made against.
        axes[2].axhline(0.5 * sep_median, color="#EE6677", ls="--", lw=1.2,
                        label="half this marker's median")
    axes[2].legend(fontsize=8, loc="upper right", framealpha=0.9)

    s = summary[(summary.marker == marker) & (summary.method == prod_method)]
    if len(s):
        r = s.iloc[0]
        head = (f"{marker}  —  production floor: {prod_method}  |  positive fraction median "
                f"{r.positive_fraction_median:.3f} (range {r.positive_fraction_min:.3f}–"
                f"{r.positive_fraction_max:.3f})  |  floor spread {r.floor_spread:.1f}×  |  "
                f"separation {r.separation_median:.3f} (min {r.separation_min:.3f})  |  "
                f"CONSISTENCY {r.separation_consistency:.2f}, collapsed on "
                f"{int(r.n_donors_separation_collapsed)}/{int(r.n_donors)} donors")
    else:
        head = f"{marker}  —  not a RESTORE gate marker (no production floor assigned)"
    fig.suptitle(head, fontsize=13, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    return fig


def build(per_donor: pd.DataFrame, summary: pd.DataFrame, outdir: Path,
          hist: pd.DataFrame | None = None) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    if hist is None:
        hp = outdir / "floor_audit_histograms.csv"
        hist = pd.read_csv(hp) if hp.exists() else pd.DataFrame(
            columns=["donor", "marker", "log_center", "density"])

    # Gate markers first, then the rest; within each, worst instability first — the reviewer should
    # meet the failures before the clean markers.
    # Worst consistency first within the gate markers: the reviewer should meet the failures before
    # the markers whose floors are sound.
    order = (summary[summary.method == summary.production_method]
             .sort_values(["is_gate_marker", "separation_consistency"],
                          ascending=[False, True])["marker"].tolist())
    order += [m for m in sorted(per_donor.marker.unique()) if m not in set(order)]

    pdf_path = outdir / "floor_audit.pdf"
    with PdfPages(pdf_path) as pdf:
        for marker in order:
            fig = marker_page(marker, per_donor, hist, summary)
            fig.savefig(figdir / f"{marker}.png", dpi=110)
            pdf.savefig(fig)
            plt.close(fig)
            per_donor[per_donor.marker == marker].to_csv(figdir / f"{marker}.csv", index=False)
    print(f"[figures] {len(order)} marker pages -> {pdf_path}")
