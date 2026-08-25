"""Figures for Step 1: one sheet per marker, plus a cohort verdict sheet.

Drawing rules this project holds to, and why they matter here:

* **Window on the workflow's own measured floor.** These values are clamped at 0 upstream
  (unsigned at the scanner, then again at ``redsea.py``), so 0 is where the axis starts. A
  percentile black point would hide the zero-inflation, which is 3–22% of cells depending on
  marker and donor and is itself a finding.
* **No gamma < 1.** The x-axis is ``arcsinh``, which is stated on the axis, monotone, and
  handles the zeros; nothing else is applied.
* **Draw the per-donor spread, never only the pooled number.** The whole reason this lane
  exists is that donors differ; a pooled curve would erase the thing being measured.
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")                     # headless: no display on this box
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
import pyarrow.parquet as pq              # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                              # noqa: E402
import usability as U                     # noqa: E402

GROUP_COLOR = {"2.3.1/CD38": "#1f77b4", "2.4.0/CD38": "#4fa3d1",
               "2.4.1/b_Catenin1": "#d62728", "2.4.1/neither": "#ff9896"}

#: The cohort sheet draws MEASUREMENTS, not a verdict.  An earlier version drew a
#: usable/degraded/unusable grid; every marker passes on every donor, so the grid carried no
#: information, and the version before it mislabelled markers by folding "this pair does not
#: mark the same cells" into "this marker is bad".  What replaces it is the continuous
#: quantity those labels were derived from, so a reader can see the margin rather than a
#: thresholded shadow of it.


def _asinh_axis(ax, hi: float):
    ax.set_xscale("asinh", linear_width=max(hi / 1e3, 1.0))
    ax.set_xlim(0, hi)                    # the measured floor, not a percentile


def load_samples(markers: list, donors: list, n_cells: int = 200_000,
                 seed: int = 0) -> dict:
    """A capped per-donor sample of every marker, read once per donor rather than per figure.

    The sheets only ever draw a density, so a 200k-cell sample per donor is visually
    identical to the full 25M and turns 1,120 column reads into 20 file reads.
    """
    rng = np.random.default_rng(seed)
    out = {}
    for d in donors:
        t = pq.read_table(panel.EXPRESSION / f"donor_id={d}" / "data.parquet",
                          columns=["qc_analysis_eligible"] + markers).to_pandas()
        t = t[t["qc_analysis_eligible"]]
        if len(t) > n_cells:
            t = t.iloc[rng.choice(len(t), n_cells, replace=False)]
        out[d] = {m: t[m].to_numpy("float64") for m in markers}
        print(f"  loaded {d}: {len(t):,} sampled cells", flush=True)
    return out


def marker_sheet(marker: str, stats: pd.DataFrame, out_dir: pathlib.Path,
                 samples: dict) -> pathlib.Path:
    """Three panels for one marker: per-donor distribution, tails, spatial structure."""
    s = stats[stats["marker"] == marker].sort_values("donor_id")
    coh = panel.cohort().set_index("donor_id")

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6), dpi=140)
    hi = float(np.nanmax(s["p99.9"])) * 1.5

    # (1) per-donor distribution, on the measured floor
    ax = axes[0]
    for _, r in s.iterrows():
        d = r["donor_id"]
        v = samples[d][marker]
        v = v[np.isfinite(v)]
        grid = np.arcsinh(np.linspace(0, hi, 400))
        cnt, edges = np.histogram(np.arcsinh(v), bins=grid)
        ax.plot(np.sinh(edges[:-1]), cnt / max(cnt.max(), 1),
                color=GROUP_COLOR.get(coh.loc[d, "acq_group"], "#777"), lw=1.2, alpha=0.85)
    _asinh_axis(ax, hi)
    ax.set_xlabel(f"{marker}  (arcsinh, from the measured floor 0)")
    ax.set_ylabel("density, scaled per donor")
    ax.set_title(f"{marker} — per-donor distribution", fontweight="bold")
    for g, c in GROUP_COLOR.items():
        ax.plot([], [], color=c, lw=2, label=g)
    ax.legend(fontsize=7, frameon=False)

    # (2) tails: declared partners against derived non-partners
    ax = axes[1]
    for i, (_, r) in enumerate(s.iterrows()):
        for name, val in json.loads(r["partners"]).items():
            ax.scatter(i, max(val, 1e-2), s=26, c="#2ca02c", marker="o", zorder=3)
        for name, val in json.loads(r["controls"]).items():
            ax.scatter(i, max(val, 1e-2), s=30, facecolors="none", edgecolors="#d62728",
                       marker="s", zorder=3)
    ax.axhline(1.0, color="k", lw=1, ls="--")     # chance; the only line worth drawing here
    ax.set_yscale("log")
    ax.set_xticks(range(len(s)))
    ax.set_xticklabels(s["donor_id"], rotation=90, fontsize=7)
    ax.set_ylabel(f"tail overlap at top {U.BARS['tail_q']:.1%}, fold over chance")
    ax.set_title("● declared partner   □ non-partner control", fontweight="bold")

    # (3) spatial structure against each donor's own permutation null
    ax = axes[2]
    x = np.arange(len(s))
    ax.plot(x, s["morans_rank"], "o-", c="#1f77b4", ms=4, label="Moran's I, ranks")
    ax.plot(x, s["morans_top1"], "s-", c="#9467bd", ms=4, label="Moran's I, top 1%")
    ax.plot(x, s["morans_null"], "-", c="k", lw=1, label="permutation null (q0.99)")
    ax.fill_between(x, 0, s["morans_null"], color="k", alpha=0.12)
    ax.set_xticks(x)
    ax.set_xticklabels(s["donor_id"], rotation=90, fontsize=7)
    ax.set_ylabel("Moran's I")
    ax.set_title("spatial coherence vs its own null", fontweight="bold")
    ax.legend(fontsize=7, frameon=False)

    pb, cw = s["partner_best"].median(), s["control_worst"].median()
    bits = [f"acquisition-group ratio {s['acq_ratio'].iloc[0]:.2f}"]
    if np.isfinite(pb):
        bits.append(f"partner tail {pb:.1f}x vs worst control {cw:.1f}x (donor medians)")
    else:
        bits.append("no declared same-cell partner in this panel")
    fig.suptitle(f"{marker}   —   " + "   |   ".join(bits), fontweight="bold", y=1.02)
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / f"{marker}.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def cohort_sheet(stats: pd.DataFrame, out_dir: pathlib.Path) -> pathlib.Path:
    """How far each marker's declared partner tail beats its non-partner controls, and
    whether its level steps with the acquisition group."""
    coh = panel.cohort().sort_values(["acq_group", "donor_id"])
    donors = coh["donor_id"].tolist()
    markers = sorted(stats["marker"].unique())

    P = stats.pivot(index="marker", columns="donor_id", values="partner_best")
    C = stats.pivot(index="marker", columns="donor_id", values="control_worst")
    M = np.log10((P.reindex(index=markers, columns=donors).to_numpy() + 0.1)
                 / (C.reindex(index=markers, columns=donors).to_numpy() + 0.1))

    fig, axes = plt.subplots(1, 2, figsize=(19, 13), dpi=140,
                             gridspec_kw={"width_ratios": [3, 1]})
    ax = axes[0]
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                   interpolation="nearest")
    ax.set_xticks(range(len(donors)))
    ax.set_xticklabels([f"{d}\n{coh.set_index('donor_id').loc[d, 'acq_group']}"
                        for d in donors], rotation=90, fontsize=6)
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels(markers, fontsize=7)
    for i, g in enumerate(coh["acq_group"]):
        if i and g != coh["acq_group"].iloc[i - 1]:
            ax.axvline(i - 0.5, color="k", lw=1.5)
    ax.set_title("log10( declared partner tail / worst non-partner control )  —  RED means "
                 "the pair marks the same cells and the controls do not\npale means it does "
                 "not resolve; blank means no same-cell partner exists in this panel. "
                 "Vertical rules separate acquisition groups", fontweight="bold", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="log10 ratio")

    ax = axes[1]
    r = stats.groupby("marker")["acq_ratio"].first().reindex(markers)
    ax.barh(range(len(markers)), r.to_numpy(),
            color=["#d62728" if v >= U.BARS["acq_ratio_flag"] else "#7f7f7f" for v in r])
    ax.axvline(U.BARS["acq_ratio_flag"], color="k", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_yticks(range(len(markers)))
    ax.set_yticklabels([])
    ax.set_ylim(len(markers) - 0.5, -0.5)
    ax.set_xlabel("acquisition-group ratio  (between / worst within)")
    ax.set_title("does the level step with the group?", fontweight="bold")

    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / "_cohort_overview.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> int:
    stats = pd.read_parquet(panel.OUT / "marker_donor_usability.parquet")
    out = panel.OUT / "figures_step1"
    markers = sorted(stats["marker"].unique())
    donors = panel.cohort()["donor_id"].tolist()
    print(cohort_sheet(stats, out))
    samples = load_samples(markers, donors)
    for m in markers:
        print(marker_sheet(m, stats, out, samples), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
