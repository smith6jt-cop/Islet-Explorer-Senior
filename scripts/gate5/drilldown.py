#!/usr/bin/env python3
"""
Gate-5 per-donor drill-down — one page per eligible donor, for the per-donor sign-off.

``build_dashboard.py`` organises the same evidence BY MARKER, which is how a marker-level problem shows
up (SST weak on 8 donors, CD11b's arms barely separated). This organises it BY DONOR, which is how a
donor is signed off. Both read the tables `build_dashboard.py` already wrote; neither recomputes a
divisor or touches production output.

    python scripts/gate5/drilldown.py [--config config.ini] [--outdir data/gate5]
    -> data/gate5/donor_drilldown.pdf   (20 pages)

WORKING RULE (root CLAUDE.md): the cohort range is drawn as CONTEXT and is explicitly labelled as such.
A donor sitting outside it is not thereby wrong -- donor biology is genuinely variable, and the only
question a page asks is whether that donor's own calls are defensible against that donor's own data.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                                # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages                           # noqa: E402
from matplotlib.patches import Patch, Rectangle                                # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenocycler.config import load_config, COMPARTMENT_COLORS                  # noqa: E402
from build_dashboard import (CLASSES, MARKERS, MARKER_GROUPS, REGIONS, REGION_LABEL,   # noqa: E402
                             STATUS_LABEL, plt as _plt)

MARKER_COMPARTMENT = {m: c for c, group in MARKER_GROUPS for m in group}
FLAG_LABEL = {
    "flag_not_valid": "state != VALID",
    "flag_low_sbr": "target_low_sbr",
    "flag_unstable_divisor": "divisor reproducibility < 0.95",
    "flag_undersized_reference": "reference < target population",
    "flag_weak_arm_separation": "target_fold < 3x (advisory, not a pipeline state)",
}
# Short codes for the table, so a row fits the panel width; expanded in the panel's own footer.
FLAG_CODE = {
    "flag_not_valid": "NOT_VALID",
    "flag_low_sbr": "LOW_SBR",
    "flag_unstable_divisor": "REPRO",
    "flag_undersized_reference": "REF_SMALL",
    "flag_weak_arm_separation": "FOLD<3",
}
THIN_REGION_N = 5000


def _page(donor: str, d: pd.DataFrame, c: pd.DataFrame, batch: pd.DataFrame | None) -> plt.Figure:
    dd = d[d.donor == donor].sort_values("marker_order")
    cc = c[c.donor == donor]
    status = dd["status"].iloc[0]
    n_cells = int(cc[cc.cell_region == "tissue"]["region_n"].iloc[0]) if len(cc) else 0
    total = int(cc.groupby("cell_region")["region_n"].first().sum())

    fig = plt.figure(figsize=(16.5, 10.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.25], height_ratios=[1.0, 1.05],
                          hspace=0.42, wspace=0.20)
    ax_comp, ax_calls = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[:, 1])
    ax_flags = fig.add_subplot(gs[1, 0])

    # --- A: this donor's composition by region -------------------------------------------------
    piv = (cc.pivot_table(index="cell_region", columns="compartment", values="pct", aggfunc="sum")
           .reindex(index=REGIONS, columns=CLASSES).fillna(0.0))
    region_n = cc.groupby("cell_region")["region_n"].first().reindex(REGIONS).fillna(0).astype(int)
    xs = np.arange(len(REGIONS))
    bottom = np.zeros(len(REGIONS))
    for compartment in CLASSES:
        v = piv[compartment].to_numpy()
        ax_comp.bar(xs, v, bottom=bottom, width=0.68, color=COMPARTMENT_COLORS[compartment],
                    edgecolor="white", linewidth=1.4, zorder=3)
        for xi, (b, h) in enumerate(zip(bottom, v)):
            if h >= 3.0:
                ax_comp.text(xs[xi], b + h / 2, f"{h:.0f}", ha="center", va="center", fontsize=8.5,
                             color="#111111" if compartment in ("Endothelial", "Other") else "white",
                             zorder=4)
        bottom += v
    ax_comp.set_xticks(xs)
    ax_comp.set_xticklabels([f"{REGION_LABEL[r]}\n{region_n[r]:,} cells"
                             + ("\n⚠ too few to read" if region_n[r] < THIN_REGION_N else "")
                             for r in REGIONS], fontsize=8.5)
    for xi, r in enumerate(REGIONS):
        if region_n[r] < THIN_REGION_N:
            ax_comp.add_patch(Rectangle((xs[xi] - 0.34, 0), 0.68, 100, fill=False, hatch="///",
                                        edgecolor="#CC3311", linewidth=0.0, zorder=5))
    ax_comp.set_ylim(0, 100)
    ax_comp.set_ylabel("% of cells")
    ax_comp.set_title("A — compartment composition by anatomical region\n"
                      "islet core is endocrine by annotation → expect Epithelial-dominant, "
                      "Endothelial enriched", fontsize=10, loc="left")
    ax_comp.legend(handles=[Patch(facecolor=COMPARTMENT_COLORS[x], label=x) for x in CLASSES],
                   loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=4, fontsize=8)

    # --- B: this donor's calls, with the cohort range as context only ---------------------------
    ys = np.arange(len(MARKERS))[::-1]
    span = d.groupby("target")["call_rate"].agg(["min", "max"]).reindex(MARKERS)
    for yi, marker in zip(ys, MARKERS):
        lo, hi = span.loc[marker, "min"], span.loc[marker, "max"]
        if np.isfinite(lo) and np.isfinite(hi):
            ax_calls.plot([lo, hi], [yi, yi], color="#CCCCCC", lw=6, solid_capstyle="butt", zorder=1)
    values = dd.set_index("target")["call_rate"].reindex(MARKERS)
    by_target = dd.set_index("target")
    ax_calls.barh(ys, values.fillna(0.0).to_numpy(), height=0.62, zorder=3,
                  color=[COMPARTMENT_COLORS[MARKER_COMPARTMENT[m]] for m in MARKERS])
    for yi, marker in zip(ys, MARKERS):
        v = values[marker]
        text = f"{v:.1f}%" if np.isfinite(v) else "no call"
        ax_calls.text(max(v if np.isfinite(v) else 0, 0) + 0.7, yi, text, va="center", fontsize=8)
    ax_calls.set_yticks(ys)
    ax_calls.set_yticklabels(MARKERS, fontsize=9)
    # A flagged pair is marked by recolouring its tick label, not by a glyph beside the axis — the glyph
    # column collided with the labels and said less than the colour does.
    flagged_targets = set(dd.loc[dd["n_flags"] > 0, "target"])
    for label, marker in zip(ax_calls.get_yticklabels(), MARKERS):
        if marker in flagged_targets:
            label.set_color("#CC3311")
            label.set_fontweight("bold")
    ax_calls.set_xlabel("% of this donor's cells called positive")
    ax_calls.set_xlim(0, max(float(np.nanmax(d["call_rate"])) * 1.08, 5))
    ax_calls.set_title("B — RESTORE calls for this donor (bars, coloured by the gate they feed)\n"
                       "grey band = the range across all 20 donors, CONTEXT ONLY — a donor outside it "
                       "is not thereby wrong.\nred bold marker name = that pair carries QC flags (see C)",
                       fontsize=10, loc="left")

    # --- C: every flagged pair for this donor, with the numbers behind the flag ------------------
    ax_flags.axis("off")
    flagged = dd[dd["n_flags"] > 0]
    header = f"C — flagged pairs for donor {donor}  ({len(flagged)} of {len(dd)})"
    footer = ("NOT_VALID = pair state is not VALID · LOW_SBR = target_low_sbr (weak arm separation)\n"
              "REPRO = divisor reproducibility < 0.95 · REF_SMALL = reference smaller than the target\n"
              "population · FOLD<3 = target_fold below 3x — ADVISORY, not a state the method reached")
    if flagged.empty:
        ax_flags.text(0, 1.0, header + "\n\nnone — every pair VALID with no donor-local flag",
                      va="top", fontsize=9, family="monospace")
    else:
        lines = [f"{'marker':<16}{'ref':<16}{'state':<20}{'tgt_f':>7}{'ref_f':>7}{'repro':>7}  flags"]
        for row in flagged.itertuples():
            marks = " ".join(FLAG_CODE[k] for k in FLAG_CODE if bool(getattr(row, k)))
            repro = row.divisor_reproducibility
            lines.append(
                f"{row.target:<16}{row.reference:<16}{row.state[:19]:<20}"
                f"{row.target_fold:>7.2f}{row.reference_fold:>7.2f}"
                f"{(f'{repro:.2f}' if np.isfinite(repro) else '   -'):>7}  {marks}")
        ax_flags.text(0, 1.0, header + "\n\n" + "\n".join(lines), va="top", fontsize=6.6,
                      family="monospace")
    ax_flags.text(0, 0.0, footer, va="bottom", fontsize=6.4, family="monospace", color="#555555")

    software = panel = "—"
    if batch is not None and donor in set(batch["donor"]):
        b = batch[batch.donor == donor].iloc[0]
        software, panel = b["acquisition_software"], b["panel_variant"]
    fig.suptitle(
        f"Gate 5 drill-down — donor {donor}  ·  {STATUS_LABEL.get(status, status)}  ·  {total:,} cells "
        f"({n_cells:,} non-islet)  ·  {software}  ·  panel {panel}",
        fontsize=15, y=0.985)
    return fig


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=REPO / "config.ini")
    ap.add_argument("--outdir", type=Path, default=REPO / "data" / "gate5")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    out = Path(a.outdir)
    dash, comp = out / "donor_marker_dashboard.csv", out / "composition_by_region.csv"
    for p in (dash, comp):
        if not p.is_file():
            raise SystemExit(f"[err] {p} not found — run scripts/gate5/build_dashboard.py first")
    d = pd.read_csv(dash, dtype={"donor": str})
    c = pd.read_csv(comp, dtype={"donor": str})
    batch_path = cfg.data_dir / "donor_acquisition_batch.csv"
    batch = None
    if batch_path.is_file():
        batch = (pd.read_csv(batch_path, dtype={"donor_id": str})
                 .rename(columns={"donor_id": "donor"}))

    donors = list(dict.fromkeys(d.sort_values(["status_order", "donor"])["donor"]))
    out_pdf = out / "donor_drilldown.pdf"
    with PdfPages(out_pdf) as pdf:
        for donor in donors:
            fig = _page(donor, d, c, batch)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            print(f"[drilldown] {donor}", flush=True)
    print(f"[saved] {out_pdf} ({len(donors)} pages)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
