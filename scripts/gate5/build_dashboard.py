#!/usr/bin/env python3
"""
Gate-5 evidence package — the 20-donor donor x marker dashboard, the regional anatomy control, and the
per-donor-marker divisor context. Writes to ``data/gate5/``.

Gate 5's exit criterion is *explicit per-donor acceptance*, and nothing in the repo currently supports
that: ``run_lineage`` emits one cohort PNG, and ``apply_manifest.csv`` is 480 rows x 85 columns. This
builds the reviewable surface.

Everything here READS artifacts already on disk (``apply_manifest.csv``, ``data/phenotype/broad/``,
``data/cells_redsea/``). Nothing is re-evaluated, no divisor is recomputed, no production output is
touched.

WORKING RULE (root CLAUDE.md): a donor is never flagged for differing from the cohort. Donor biology is
genuinely variable, so every flag drawn here is DONOR-LOCAL -- the pair's own state, its own low-SBR
diagnostic, its own reference sizing, its own seed reproducibility -- and every call is shown against
that donor's OWN REDSEA distribution. The cohort spread appears only as context a reader may ignore,
never as a threshold.

Colour policy follows ``scripts/anatomy_check/figures.py``: compartments use ``COMPARTMENT_COLORS``
unchanged so the figure agrees with ``lineage._composition_figure`` and the QuPath overlays the
maintainer inspects. That palette does not pass the CVD validator cleanly (chroma floor on #4477AA /
#B5838D, lightness on #66CCEE), so identity is carried redundantly: fixed compartment order, a legend
always present, and direct labels on every segment >= 3%.

    python scripts/gate5/build_dashboard.py [--config config.ini] [--outdir data/gate5]

Outputs
    donor_marker_dashboard.csv / .png   20 donors x 24 targets: call rate, state, donor-local flags
    composition_by_region.csv / .png    donor x {islet core, peri-islet, tissue} x compartment
    divisor_context.csv / .pdf          per donor-marker: the donor's own REDSEA distribution + divisor
    gate5_signoff.csv                   the maintainer's one editable surface (blank disposition/notes)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                                # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages                           # noqa: E402
from matplotlib.colors import PowerNorm                                        # noqa: E402
from matplotlib.lines import Line2D                                            # noqa: E402
from matplotlib.patches import Patch, Rectangle                                # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
from phenocycler.config import (PipelineConfig, load_config, COMPARTMENT_COLORS,   # noqa: E402
                                COMPARTMENT_GATES, COMPARTMENT_ORDER, OTHER_LABEL,
                                UNRESOLVED_LABEL, STATUS_ORDER)
from phenocycler.lineage import status_map                                      # noqa: E402

CLASSES = COMPARTMENT_ORDER + [OTHER_LABEL, UNRESOLVED_LABEL]
REGIONS = ["core", "peri", "tissue"]
REGION_LABEL = {"core": "islet core", "peri": "peri-islet", "tissue": "tissue (non-islet)"}
STATUS_LABEL = {"ND": "ND", "AAB": "Aab+", "T1D": "T1D"}
STATUS_COLORS = {"ND": "#4477AA", "AAB": "#CC6633", "T1D": "#228833"}       # Paul Tol bright, as the app

# Gate markers grouped by the compartment they gate, in COMPARTMENT_ORDER. The dashboard's column order
# IS the ordered-residual priority, so a reader sees which gate had first claim on a cell.
MARKER_GROUPS = [(c, COMPARTMENT_GATES[c]) for c in COMPARTMENT_ORDER]
MARKERS = [m for _, group in MARKER_GROUPS for m in group]

N_HIST_BINS = 60

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "legend.frameon": False,
})


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def load_manifest(cfg: PipelineConfig) -> pd.DataFrame:
    """The Gate-4 cohort QC disposition, one row per (donor, accepted pair), + a derived call rate."""
    path = cfg.restore_apply_manifest_csv
    if not path.is_file():
        raise SystemExit(f"[err] no apply manifest at {path} — run phenocycler.restore_apply first")
    m = pd.read_csv(path, dtype={"donor": str}, low_memory=False)
    called = m[["called_positive_n", "called_negative_n", "called_unavailable_n"]].sum(axis=1)
    m["n_called"] = called
    m["call_rate"] = 100.0 * m["called_positive_n"] / called.replace(0, np.nan)
    m["unavailable_rate"] = 100.0 * m["called_unavailable_n"] / called.replace(0, np.nan)
    missing = sorted(set(MARKERS) - set(m["target"]))
    if missing:
        raise SystemExit(
            f"[err] manifest is missing gate markers {missing}; config.COMPARTMENT_GATES and the "
            "accepted pair set have drifted apart"
        )
    return m


def load_composition(cfg: PipelineConfig, donors: list[str]) -> pd.DataFrame:
    """Per (donor, cell_region, compartment) cell counts, streamed one donor at a time."""
    rows = []
    for donor in donors:
        t = ds.dataset(cfg.broad_dir / f"donor_id={donor}", format="parquet").to_table(
            columns=["cell_region", "compartment"]).to_pandas()
        counts = t.value_counts(["cell_region", "compartment"], dropna=False).reset_index(name="n")
        counts.insert(0, "donor", donor)
        rows.append(counts)
        print(f"[composition] {donor}: {len(t):,} cells", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out["cell_region"] = out["cell_region"].fillna("tissue")
    return out


def load_marker_histograms(cfg: PipelineConfig, donors: list[str], manifest: pd.DataFrame) -> pd.DataFrame:
    """Per (donor, marker) histogram of that donor's OWN REDSEA target values, on shared log bins.

    Bins are shared PER MARKER so one page is readable; the counts, the divisor and the call rate are
    each entirely the donor's own. Zeros are counted separately -- REDSEA clamps negative results to
    zero and does so unevenly across markers (E-cadherin 8.0% of tissue cells against EpCAM 0.5%), so
    folding them into the first bin would hide a marker-specific floor of forced negatives."""
    edges = _shared_log_edges(cfg, donors)
    rows = []
    for donor in donors:
        v = ds.dataset(cfg.cells_redsea_dir / f"donor_id={donor}", format="parquet").to_table(
            columns=MARKERS).to_pandas()
        for marker in MARKERS:
            x = v[marker].to_numpy(dtype=float)
            finite = np.isfinite(x)
            zero = finite & (x <= 0)
            counts, _ = np.histogram(x[finite & ~zero], bins=edges[marker])
            rows.append({
                "donor": donor, "marker": marker,
                "n_finite": int(finite.sum()), "n_zero": int(zero.sum()),
                **{f"bin_{i}": int(c) for i, c in enumerate(counts)},
            })
        print(f"[histograms] {donor}: {len(v):,} cells x {len(MARKERS)} markers", flush=True)
    hist = pd.DataFrame(rows)
    edge_table = pd.DataFrame({m: edges[m] for m in MARKERS})
    return hist, edge_table


def _shared_log_edges(cfg: PipelineConfig, donors: list[str]) -> dict[str, np.ndarray]:
    """Per-marker log-spaced bin edges spanning the cohort's own positive REDSEA range.

    This is an AXIS, not a threshold: it decides only what a page's x-limits are. Each donor's counts,
    divisor and call rate are computed from its own values alone."""
    lo = {m: np.inf for m in MARKERS}
    hi = {m: 0.0 for m in MARKERS}
    for donor in donors:
        v = ds.dataset(cfg.cells_redsea_dir / f"donor_id={donor}", format="parquet").to_table(
            columns=MARKERS).to_pandas()
        for marker in MARKERS:
            x = v[marker].to_numpy(dtype=float)
            x = x[np.isfinite(x) & (x > 0)]
            if x.size:
                lo[marker] = min(lo[marker], float(np.quantile(x, 0.001)))
                hi[marker] = max(hi[marker], float(np.quantile(x, 0.9999)))
        print(f"[bins] scanned {donor}", flush=True)
    return {
        m: np.geomspace(max(lo[m], 1e-3), max(hi[m], max(lo[m], 1e-3) * 10), N_HIST_BINS + 1)
        for m in MARKERS
    }


# --------------------------------------------------------------------------- #
# 1 — donor x marker dashboard
# --------------------------------------------------------------------------- #
def dashboard_table(manifest: pd.DataFrame, smap: dict) -> pd.DataFrame:
    """One row per (donor, gate marker) with the call and every donor-local QC flag beside it."""
    keep = ["donor", "target", "reference", "state", "state_reason", "expression_review",
            "call_rate", "unavailable_rate", "called_positive_n", "n_called", "canonical_divisor",
            "target_fold", "reference_fold", "target_low_sbr", "reference_undersized",
            "reference_to_target_ratio", "divisor_reproducibility", "anchor_recovery",
            "feature_compartments", "control_definition", "divisor_statistic"]
    d = manifest[keep].copy()
    d["status"] = [smap.get(x, "?") for x in d["donor"]]
    d["compartment_gated"] = d["target"].map(
        {m: c for c, group in MARKER_GROUPS for m in group})
    # Donor-local flags only. None of these compares the donor to another donor.
    d["flag_not_valid"] = d["state"] != "VALID"
    d["flag_low_sbr"] = d["target_low_sbr"].astype(bool)
    d["flag_undersized_reference"] = d["reference_undersized"].astype(bool)
    d["flag_unstable_divisor"] = d["divisor_reproducibility"] < 0.95
    # ADVISORY, and mine -- not a pipeline state. The method's own floor is 2x (below it the pair is
    # `target_low_sbr`), so this marks pairs sitting just above that floor: CD11b's cohort-median
    # target_fold is 2.89 and EpCAM's 3.15, against 20-27x for the other keratins. It is reported in the
    # CSV and named in the drill-down, and it is deliberately NOT drawn on the heatmap, so it is never
    # mistaken for a verdict the frozen method reached.
    d["flag_weak_arm_separation"] = d["target_fold"] < 3.0     # the pair's own arms, not a frequency
    d["n_flags"] = d[[c for c in d.columns if c.startswith("flag_")]].sum(axis=1)
    d["marker_order"] = d["target"].map({m: i for i, m in enumerate(MARKERS)})
    d["status_order"] = d["status"].map({s: i for i, s in enumerate(STATUS_ORDER)}).fillna(99)
    return d.sort_values(["status_order", "donor", "marker_order"], ignore_index=True)


def figure_dashboard(d: pd.DataFrame, out_png: Path) -> None:
    donors = list(dict.fromkeys(d.sort_values(["status_order", "donor"])["donor"]))
    grid = d.pivot(index="donor", columns="target", values="call_rate").reindex(
        index=donors, columns=MARKERS)

    fig, ax = plt.subplots(figsize=(1 + 0.62 * len(MARKERS), 1.6 + 0.42 * len(donors)))
    # Sequential, ONE hue, light->dark (never a rainbow). gamma < 1 so a 0.4% marker and a 60% marker
    # are both legible on one absolute scale -- the numbers are printed in-cell regardless.
    im = ax.imshow(grid.to_numpy(dtype=float), cmap="Blues", aspect="auto",
                   norm=PowerNorm(gamma=0.45, vmin=0, vmax=float(np.nanmax(grid.to_numpy()))))
    ax.set_xticks(range(len(MARKERS)))
    ax.set_xticklabels(MARKERS, rotation=90, fontsize=8)
    ax.set_yticks(range(len(donors)))
    ax.set_yticklabels([f"{x}  {STATUS_LABEL.get(d.loc[d.donor == x, 'status'].iloc[0], '?')}"
                        for x in donors], fontsize=8)
    ax.grid(False)

    flags = d.set_index(["donor", "target"])
    for yi, donor in enumerate(donors):
        for xi, marker in enumerate(MARKERS):
            row = flags.loc[(donor, marker)]
            value = row["call_rate"]
            if np.isfinite(value):
                ax.text(xi, yi, f"{value:.1f}", ha="center", va="center", fontsize=6.4,
                        color="white" if value > 22 else "#111111", zorder=3)
            if row["flag_not_valid"]:                      # state is the loudest signal -> a full ring
                ax.add_patch(Rectangle((xi - 0.5, yi - 0.5), 1, 1, fill=False,
                                       edgecolor="#CC3311", linewidth=2.0, zorder=4))
            if row["flag_low_sbr"]:                        # weak arms -> a corner wedge, top-left
                ax.plot([xi - 0.5, xi - 0.16], [yi - 0.16, yi - 0.5], color="#EE7733", lw=1.6, zorder=4)
            if row["flag_unstable_divisor"]:               # seed spread -> corner wedge, top-right
                ax.plot([xi + 0.16, xi + 0.5], [yi - 0.5, yi - 0.16], color="#AA3377", lw=1.6, zorder=4)
            if row["flag_undersized_reference"]:           # reference smaller than the target population
                ax.plot([xi], [yi + 0.34], marker=".", color="#555555", ms=3.5, zorder=4)

    for boundary in np.cumsum([len(g) for _, g in MARKER_GROUPS])[:-1]:
        ax.axvline(boundary - 0.5, color="#333333", lw=1.2, zorder=5)
    xpos = 0
    for compartment, group in MARKER_GROUPS:
        ax.text(xpos + (len(group) - 1) / 2, -1.15, compartment, ha="center", va="bottom", fontsize=9,
                color=COMPARTMENT_COLORS[compartment], fontweight="bold")
        xpos += len(group)
    seen = 0
    for status in STATUS_ORDER:
        n = sum(1 for x in donors if d.loc[d.donor == x, "status"].iloc[0] == status)
        if n and seen:
            ax.axhline(seen - 0.5, color="#333333", lw=1.2, zorder=5)
        seen += n

    fig.colorbar(im, ax=ax, fraction=0.018, pad=0.012).set_label("% of the donor's cells called +", fontsize=9)
    ax.set_title("Gate 5 — RESTORE calls per donor x gate marker, with each pair's own QC flags\n"
                 "columns in ordered-residual priority; every flag is donor-local (no cohort comparison)",
                 loc="left", fontsize=12, pad=34)
    fig.legend(handles=[
        Line2D([], [], color="#CC3311", lw=2.0, label="state != VALID"),
        Line2D([], [], color="#EE7733", lw=1.6, label="target_low_sbr (weak arm separation)"),
        Line2D([], [], color="#AA3377", lw=1.6, label="divisor reproducibility < 0.95"),
        Line2D([], [], color="#555555", marker=".", ls="", label="reference smaller than target population"),
    ], loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02), fontsize=9)
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[saved] {out_png}")


# --------------------------------------------------------------------------- #
# 2 — composition by anatomical region (the external control)
# --------------------------------------------------------------------------- #
def composition_table(comp: pd.DataFrame, smap: dict) -> pd.DataFrame:
    total = comp.groupby(["donor", "cell_region"])["n"].transform("sum")
    out = comp.assign(pct=100.0 * comp["n"] / total, region_n=total)
    out["status"] = [smap.get(x, "?") for x in out["donor"]]
    out["status_order"] = out["status"].map({s: i for i, s in enumerate(STATUS_ORDER)}).fillna(99)
    return out.sort_values(["status_order", "donor", "cell_region"], ignore_index=True)


def figure_composition(c: pd.DataFrame, out_png: Path, *, min_region_n: int = 5000) -> None:
    """Stacked composition per donor x region.

    A region with few cells is shown but marked: donor 6539 has 124 peri-islet cells and donor 6380 has
    1,294 islet-core cells, so those bars are percentages of almost nothing and must not be read as
    composition. ``min_region_n`` only controls that MARK -- no donor is excluded and no value is
    altered."""
    donors = list(dict.fromkeys(c.sort_values(["status_order", "donor"])["donor"]))
    fig, axes = plt.subplots(len(REGIONS), 1, figsize=(2.2 + 0.62 * len(donors), 3.6 * len(REGIONS)),
                             sharex=True)
    for ax, region in zip(np.atleast_1d(axes), REGIONS):
        piv = (c[c.cell_region == region]
               .pivot_table(index="donor", columns="compartment", values="pct", aggfunc="sum")
               .reindex(index=donors, columns=CLASSES).fillna(0.0))
        n_region = (c[c.cell_region == region].groupby("donor")["region_n"].first()
                    .reindex(donors).fillna(0).astype(int))
        xs = np.arange(len(donors))
        bottom = np.zeros(len(donors))
        for compartment in CLASSES:
            v = piv[compartment].to_numpy()
            ax.bar(xs, v, bottom=bottom, width=0.74, color=COMPARTMENT_COLORS[compartment],
                   edgecolor="white", linewidth=1.4, zorder=3)
            for xi, (b, h) in enumerate(zip(bottom, v)):
                if h >= 3.0:                               # direct labels: identity is never colour alone
                    ax.text(xs[xi], b + h / 2, f"{h:.0f}", ha="center", va="center", fontsize=7,
                            color="#111111" if compartment in ("Endothelial", OTHER_LABEL) else "white",
                            zorder=4)
            bottom += v
        ax.set_ylim(0, 128)                                # headroom so the count row never meets the title
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.set_xlim(-0.6, len(donors) - 0.4)
        ax.set_ylabel(f"{REGION_LABEL[region]}\n% of cells")
        for xi, n in enumerate(n_region):                  # a 124-cell peri population is not evidence
            thin = n < min_region_n
            ax.text(xs[xi], 103, f"{n:,}" + (" ⚠" if thin else ""), ha="center", va="bottom",
                    fontsize=6.5, color="#CC3311" if thin else "#555555", rotation=90)
            if thin:                                       # hatch the bar itself, so the mark is unmissable
                ax.add_patch(Rectangle((xs[xi] - 0.37, 0), 0.74, 100, fill=False, hatch="///",
                                       edgecolor="#CC3311", linewidth=0.0, zorder=5))
        ax.axhline(100, color="#999999", lw=0.8, zorder=6)
        ax.set_xticks(xs)
    axes[-1].set_xticklabels(
        [f"{x}  {STATUS_LABEL.get(c.loc[c.donor == x, 'status'].iloc[0], '?')}" for x in donors],
        rotation=90, fontsize=8)
    axes[0].set_title(
        "Gate 5 — compartment composition per donor, by anatomical region\n"
        "Islet-core cells are endocrine by QuPath annotation and hormone-validated, so islet core "
        "SHOULD be epithelial-dominant and endothelium SHOULD be enriched there. This is the one check "
        "independent of every marker threshold.\n"
        f"Grey numbers above each bar = cells in that region for that donor; red + hatched = fewer than "
        f"{min_region_n:,}, i.e. a percentage of too little to read as composition.",
        loc="left", fontsize=11, pad=16)
    fig.legend(handles=[Patch(facecolor=COMPARTMENT_COLORS[x], label=x) for x in CLASSES],
               loc="lower center", ncol=7, bbox_to_anchor=(0.5, -0.015), fontsize=9)
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[saved] {out_png}")


# --------------------------------------------------------------------------- #
# 3 — divisor against each donor's own distribution
# --------------------------------------------------------------------------- #
def divisor_context_pdf(hist: pd.DataFrame, edges: pd.DataFrame, d: pd.DataFrame, out_pdf: Path) -> None:
    """One page per gate marker; 20 small panels, each a donor's OWN REDSEA histogram + its OWN divisor."""
    bin_cols = [c for c in hist.columns if c.startswith("bin_")]
    donors = list(dict.fromkeys(d.sort_values(["status_order", "donor"])["donor"]))
    flags = d.set_index(["donor", "target"])
    with PdfPages(out_pdf) as pdf:
        for marker in MARKERS:
            e = edges[marker].to_numpy()
            centres = np.sqrt(e[:-1] * e[1:])
            ncol = 5
            nrow = int(np.ceil(len(donors) / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.35 * nrow))
            axes = np.atleast_1d(axes).ravel()
            for ax, donor in zip(axes, donors):
                row = hist[(hist.donor == donor) & (hist.marker == marker)]
                info = flags.loc[(donor, marker)]
                counts = row[bin_cols].to_numpy(dtype=float).ravel()
                ax.fill_between(centres, counts, step="mid", color="#BBBBBB", zorder=2)
                ax.set_xscale("log")
                ax.set_yscale("log")
                divisor = float(info["canonical_divisor"])
                if np.isfinite(divisor):
                    ax.axvline(divisor, color="#CC3311", lw=1.6, zorder=4)
                zero_pct = 100.0 * int(row["n_zero"].iloc[0]) / max(int(row["n_finite"].iloc[0]), 1)
                subtitle = (f"{info['call_rate']:.1f}%+" if np.isfinite(info["call_rate"]) else "no call")
                ax.set_title(f"{donor} {STATUS_LABEL.get(info['status'], '?')} · {subtitle}"
                             f" · {zero_pct:.1f}% zero", fontsize=8.5,
                             color="#CC3311" if info["flag_not_valid"] else "#111111")
                ax.tick_params(labelsize=7)
                marks = [k.replace("flag_", "") for k in
                         ("flag_low_sbr", "flag_unstable_divisor", "flag_undersized_reference")
                         if bool(info[k])]
                if marks:
                    ax.text(0.02, 0.95, "\n".join(marks), transform=ax.transAxes, fontsize=6.5,
                            va="top", color="#EE7733")
            for ax in axes[len(donors):]:
                ax.set_visible(False)
            fig.suptitle(
                f"{marker} ← {flags.loc[(donors[0], marker), 'reference']} — each donor's own "
                f"REDSEA distribution with its own divisor (red)\n"
                "x = REDSEA whole-cell mean (log) · y = cells (log) · zeros excluded from the curve and "
                "reported per panel (REDSEA clamps negatives to 0)",
                fontsize=12, y=1.0)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig, bbox_inches="tight")     # never clip the suptitle or the outer tick labels
            plt.close(fig)
            print(f"[divisor-context] {marker}", flush=True)
    print(f"[saved] {out_pdf}")


# --------------------------------------------------------------------------- #
# 4 — the sign-off surface
# --------------------------------------------------------------------------- #
def signoff_table(d: pd.DataFrame, c: pd.DataFrame, batch: pd.DataFrame | None) -> pd.DataFrame:
    """One row per donor, with the facts a sign-off needs and two blank columns for the maintainer."""
    per_donor = d.groupby("donor").agg(
        status=("status", "first"),
        n_pairs=("target", "size"),
        n_not_valid=("flag_not_valid", "sum"),
        n_low_sbr=("flag_low_sbr", "sum"),
        n_unstable_divisor=("flag_unstable_divisor", "sum"),
        n_undersized_reference=("flag_undersized_reference", "sum"),
        n_weak_arm_separation=("flag_weak_arm_separation", "sum"),
    )
    flagged = (d[d["n_flags"] > 0].groupby("donor")["target"]
               .apply(lambda s: " ".join(sorted(s))).rename("flagged_markers"))
    per_donor["flagged_markers"] = flagged.reindex(per_donor.index).fillna("")

    region_n = c.pivot_table(index="donor", columns="cell_region", values="region_n", aggfunc="first")
    out = per_donor.join(region_n.rename(columns={r: f"n_{r}" for r in REGIONS}))
    for label in ("tissue", "core"):
        frame = (c[c.cell_region == label]
                 .pivot_table(index="donor", columns="compartment", values="pct", aggfunc="sum")
                 .reindex(index=out.index, columns=CLASSES))
        for compartment in CLASSES:
            out[f"{label}_{compartment}"] = frame[compartment].round(2)
    if batch is not None:
        out = out.join(batch.set_index("donor"), how="left")
    out["disposition"] = ""            # ACCEPT | REJECT | HOLD  <- the maintainer fills these two
    out["notes"] = ""
    return out.reset_index()


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=REPO / "config.ini")
    ap.add_argument("--outdir", type=Path, default=REPO / "data" / "gate5")
    ap.add_argument("--skip-histograms", action="store_true",
                    help="reuse a cached marker_histograms.csv instead of re-reading cells_redsea")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    donors = cfg.discover_donors(cfg.broad_dir)
    if not donors:
        raise SystemExit(f"[err] no broad output under {cfg.broad_dir} — run phenocycler.lineage first")
    print(f"[gate5] {len(donors)} donors, {len(MARKERS)} gate markers -> {out}")

    manifest = load_manifest(cfg)
    smap = status_map(cfg)

    d = dashboard_table(manifest, smap)
    d.to_csv(out / "donor_marker_dashboard.csv", index=False)
    figure_dashboard(d, out / "donor_marker_dashboard.png")

    comp = load_composition(cfg, donors)
    c = composition_table(comp, smap)
    c.to_csv(out / "composition_by_region.csv", index=False)
    figure_composition(c, out / "composition_by_region.png")

    hist_csv, edge_csv = out / "marker_histograms.csv", out / "marker_histogram_edges.csv"
    if a.skip_histograms and hist_csv.is_file() and edge_csv.is_file():
        hist = pd.read_csv(hist_csv, dtype={"donor": str})
        edges = pd.read_csv(edge_csv)
    else:
        hist, edges = load_marker_histograms(cfg, donors, manifest)
        hist.to_csv(hist_csv, index=False)
        edges.to_csv(edge_csv, index=False)
    divisor_context_pdf(hist, edges, d, out / "divisor_context.pdf")

    # Acquisition batch is a Gate-5 covariate: ND is aliased with Fusion 2.4.1 (Fisher p = 0.0147), so a
    # per-donor disposition has to be readable against the software/panel variant, not only the label.
    batch_path = cfg.data_dir / "donor_acquisition_batch.csv"
    batch = None
    if batch_path.is_file():
        batch = (pd.read_csv(batch_path, dtype={"donor_id": str})
                 .rename(columns={"donor_id": "donor"})
                 [["donor", "acquisition_software", "n_channels", "panel_variant"]])
    signoff = signoff_table(d, c, batch)
    signoff.to_csv(out / "gate5_signoff.csv", index=False)
    print(f"[saved] {out / 'gate5_signoff.csv'} — the maintainer's editable surface")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
