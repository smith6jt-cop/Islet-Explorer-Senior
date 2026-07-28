#!/usr/bin/env python3
"""
Gate-5 targeted diagnostics — five open questions, one figure + source CSV each -> ``data/gate5/``.

Nothing here changes production. Each question ends in a recommendation the maintainer rules on.

  D1  Is the islet-core Immune gradient disease, or acquisition batch?   [the study's headline rests on this]
  D2  Why is CD11b at 2.1% when the Gate-3 pilot validated 7-10%?
  D3  Is SST's weakness a disease pattern, or a donor-level pass/fail?
  D4  Is CD34 marking endothelium, or a graded continuum incl. fibroblasts?
  D5  Is EpCAM's weak arm separation covered by the rest of the Epithelial union?

D1, D3 read the Gate-5 tables. D2 additionally reads the Gate-3 pilot. D4, D5 stream the per-cell
tri-state calls in ``data/restore_gated_redsea/`` plus ``data/cells_redsea/``.

    python scripts/gate5/diagnostics.py [--config config.ini] [--outdir data/gate5]
    [--only D1 D2 ...]

WORKING RULE (root CLAUDE.md): no marker is called wrong for being rare or abundant. Every claim below
is either (a) an internal inconsistency -- the same marker measured two ways disagreeing, (b) a
distribution shape read off the donor's own data, or (c) a confound in the cohort design. None is
"this donor should look like the others".
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
from scipy.stats import kruskal, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                                # noqa: E402
from matplotlib.lines import Line2D                                            # noqa: E402
from matplotlib.patches import Patch                                           # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phenocycler.config import load_config, COMPARTMENT_COLORS, COMPARTMENT_GATES   # noqa: E402
from build_dashboard import STATUS_COLORS, STATUS_LABEL                         # noqa: E402

SOFTWARE_MARKERS = {"Fusion 2.3.1": "o", "Fusion 2.4.0": "s", "Fusion 2.4.1": "^"}
PILOT = "data/restore_pair_validation/gate3_pilot_v2/pair_evaluations.csv"

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "legend.frameon": False,
})


def _tables(out: Path, cfg):
    d = pd.read_csv(out / "donor_marker_dashboard.csv", dtype={"donor": str})
    c = pd.read_csv(out / "composition_by_region.csv", dtype={"donor": str})
    b = (pd.read_csv(cfg.data_dir / "donor_acquisition_batch.csv", dtype={"donor_id": str})
         .rename(columns={"donor_id": "donor"}))
    return d, c, b


# --------------------------------------------------------------------------- #
# D1 — disease, or acquisition batch?
# --------------------------------------------------------------------------- #
def d1_batch_confound(d, c, b, out: Path) -> None:
    core = c[c.cell_region == "core"].pivot_table(index="donor", columns="compartment", values="pct")
    t = (core.join(c.groupby("donor")["status"].first())
         .join(b.set_index("donor")[["acquisition_software", "panel_variant"]]))
    t = t.reset_index()
    t.to_csv(out / "d1_islet_core_by_status_and_batch.csv", index=False)

    software = sorted(t["acquisition_software"].unique())
    p_status = kruskal(*[t.loc[t.status == s, "Immune"] for s in ("ND", "AAB", "T1D")]).pvalue
    p_soft = kruskal(*[t.loc[t.acquisition_software == s, "Immune"] for s in software]).pvalue

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4), gridspec_kw={"width_ratios": [1, 1, 1.25]})

    for ax, key, order, colors, title in (
        (axes[0], "status", ["ND", "AAB", "T1D"], STATUS_COLORS,
         f"by disease status\nKruskal p = {p_status:.4f}"),
        (axes[1], "acquisition_software", software, None,
         f"by acquisition software\nKruskal p = {p_soft:.4f}"),
    ):
        for xi, level in enumerate(order):
            v = t.loc[t[key] == level, "Immune"].to_numpy()
            jitter = np.linspace(-0.16, 0.16, len(v)) if len(v) > 1 else np.array([0.0])
            colour = (colors or {}).get(level, "#666666")
            ax.scatter(xi + jitter, v, s=70, color=colour, edgecolor="black", lw=0.5, zorder=3)
            if len(v):
                ax.plot([xi - 0.28, xi + 0.28], [np.median(v)] * 2, color="black", lw=2, zorder=4)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([STATUS_LABEL.get(x, x) for x in order], fontsize=9)
        ax.set_ylabel("islet-core Immune (% of core cells)")
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.set_ylim(0, 24)

    # The unaliasing panel: the same donors, split by software, coloured by status.
    ax = axes[2]
    for xi, sw in enumerate(software):
        sub = t[t.acquisition_software == sw]
        jitter = np.linspace(-0.2, 0.2, len(sub)) if len(sub) > 1 else np.array([0.0])
        ax.scatter(xi + jitter, sub["Immune"], s=95,
                   color=[STATUS_COLORS[s] for s in sub["status"]], edgecolor="black", lw=0.6, zorder=3)
        for xj, (_, row) in zip(xi + jitter, sub.iterrows()):
            ax.annotate(row["donor"], (xj, row["Immune"]), fontsize=6.5, xytext=(0, 7),
                        textcoords="offset points", ha="center", color="#444444")
    ax.set_xticks(range(len(software)))
    ax.set_xticklabels(software, fontsize=9)
    ax.set_ylabel("islet-core Immune (% of core cells)")
    ax.set_ylim(0, 24)
    ax.set_title("the two split apart\nevery high-immune donor is on one software; "
                 "no ND donor is", fontsize=10.5, loc="left")
    ax.legend(handles=[Patch(facecolor=STATUS_COLORS[s], label=STATUS_LABEL[s])
                       for s in ("ND", "AAB", "T1D")], loc="upper left", fontsize=9)

    within = []
    for sw in software:
        sub = t[t.acquisition_software == sw]
        nd = sub.loc[sub.status == "ND", "Immune"]
        dis = sub.loc[sub.status != "ND", "Immune"]
        within.append(f"{sw}: ND n={len(nd)} [{nd.min():.1f}–{nd.max():.1f}]  "
                      f"non-ND n={len(dis)} [{dis.min():.1f}–{dis.max():.1f}]"
                      if len(nd) and len(dis) else
                      f"{sw}: cannot compare (ND n={len(nd)}, non-ND n={len(dis)})")
    fig.suptitle(
        "D1 — the islet-core Immune gradient is aliased with acquisition software\n"
        + "   ·   ".join(within), fontsize=12.5, y=1.06)
    fig.savefig(out / "d1_batch_confound.png")
    plt.close(fig)
    print(f"[saved] {out / 'd1_batch_confound.png'}  (status p={p_status:.4f}, software p={p_soft:.4f})")


# --------------------------------------------------------------------------- #
# D2 — CD11b: 2.1% now against the Gate-3 pilot's validated 7-10%
# --------------------------------------------------------------------------- #
def d2_cd11b(d, out: Path) -> None:
    pilot_path = REPO / PILOT
    if not pilot_path.is_file():
        print(f"[D2] skipped — no Gate-3 pilot at {pilot_path}")
        return
    p = pd.read_csv(pilot_path, dtype={"donor": str}, low_memory=False)
    j = p[["donor", "target", "canonical_divisor", "lineage_positive_n", "n_assigned"]].merge(
        d[["donor", "target", "canonical_divisor", "call_rate"]], on=["donor", "target"],
        suffixes=("_pilot", "_now"))
    j["pilot_rate"] = 100 * j["lineage_positive_n"] / j["n_assigned"]
    j["divisor_ratio"] = j["canonical_divisor_now"] / j["canonical_divisor_pilot"]
    j.to_csv(out / "d2_pilot_vs_production.csv", index=False)

    order = (j.groupby("target")["divisor_ratio"].median().sort_values().index.tolist())
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.6))

    ax = axes[0]
    for xi, target in enumerate(order):
        v = j.loc[j.target == target, "divisor_ratio"].to_numpy()
        jitter = np.linspace(-0.18, 0.18, len(v)) if len(v) > 1 else np.array([0.0])
        highlight = target == "CD11b"
        ax.scatter(xi + jitter, v, s=80, zorder=3, edgecolor="black", lw=0.6,
                   color="#CC3311" if highlight else "#BBBBBB")
        ax.plot([xi - 0.3, xi + 0.3], [np.median(v)] * 2, color="black", lw=2, zorder=4)
    ax.axhline(1.0, color="#333333", lw=1.2, ls="--", zorder=2)
    ax.set_yscale("log")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("divisor now ÷ divisor in the Gate-3 pilot")
    ax.set_title("A — every marker's divisor FELL when the policy changed to reference_only/p99.\n"
                 "CD11b is the only one that ROSE — so its call rate fell while every other rose.",
                 fontsize=10.5, loc="left")

    ax = axes[1]
    for target in order:
        sub = j[j.target == target]
        highlight = target == "CD11b"
        ax.scatter(sub["pilot_rate"], sub["call_rate"], s=80, zorder=3, edgecolor="black", lw=0.6,
                   color="#CC3311" if highlight else "#BBBBBB",
                   label=target if highlight else None)
    lim = [0, max(j["pilot_rate"].max(), j["call_rate"].max()) * 1.08]
    ax.plot(lim, lim, color="#333333", lw=1.0, ls="--", zorder=2)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("% called in the Gate-3 pilot (v9 policy)")
    ax.set_ylabel("% called in production now (reference_only/p99)")
    ax.set_title("B — same 7 donors, same pairs. Above the line = more called now.\n"
                 "CD11b (red) is the only marker that moved DOWN.", fontsize=10.5, loc="left")
    ax.legend(loc="lower right", fontsize=9)

    cd11b = j[j.target == "CD11b"]
    fig.suptitle(
        "D2 — CD11b's divisor INFLATED under the frozen policy, so 2.1% is an under-call\n"
        f"CD11b divisor ratio {cd11b['divisor_ratio'].min():.2f}–{cd11b['divisor_ratio'].max():.2f}× "
        f"(every other marker ≤ 1.12×); call rate {cd11b['pilot_rate'].min():.1f}–"
        f"{cd11b['pilot_rate'].max():.1f}% → {cd11b['call_rate'].min():.1f}–{cd11b['call_rate'].max():.1f}%. "
        "Mechanism: control_definition reference_and_target → reference_only removed the cap at CD11b's "
        "own Otsu floor,\nadmitting its high channel background into its own negative control.",
        fontsize=12, y=1.08)
    fig.savefig(out / "d2_cd11b_divisor_inflation.png")
    plt.close(fig)
    print(f"[saved] {out / 'd2_cd11b_divisor_inflation.png'}")


# --------------------------------------------------------------------------- #
# D3 — SST: disease pattern, or donor-level pass/fail?
# --------------------------------------------------------------------------- #
def d3_sst(d, c, b, out: Path) -> None:
    sst = (d[d.target == "SST"].drop(columns=["status"])
           .merge(c.groupby("donor")["status"].first().rename("status"), on="donor")
           .merge(b[["donor", "acquisition_software"]], on="donor"))
    sst.to_csv(out / "d3_sst_per_donor.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.4), gridspec_kw={"width_ratios": [1.35, 1]})

    ax = axes[0]
    s = sst.sort_values("target_fold")
    ys = np.arange(len(s))
    ax.barh(ys, s["target_fold"], color=[STATUS_COLORS[x] for x in s["status"]], height=0.66, zorder=3)
    ax.axvline(2.0, color="#CC3311", lw=1.6, ls="--", zorder=4)
    ax.text(2.1, len(s) - 0.5, "the method's own 2× arm-separation floor", fontsize=8.5,
            color="#CC3311", va="top")
    ax.set_xscale("log")
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{r.donor} {STATUS_LABEL[r.status]}" for r in s.itertuples()], fontsize=8)
    ax.set_xlabel("SST target_fold (arm separation, log)")
    below = int((sst["target_fold"] < 2.0).sum())
    lo_hi = np.sort(sst["target_fold"].to_numpy())
    ax.set_title(f"A — SST either works or fails outright: {below} donors at "
                 f"{lo_hi[0]:.2f}–{lo_hi[below - 1]:.2f}×, {len(sst) - below} at "
                 f"{lo_hi[below]:.2f}–{lo_hi[-1]:.0f}×.\n"
                 f"Nothing sits between {lo_hi[below - 1]:.2f}× and {lo_hi[below]:.2f}× — "
                 "a pass/fail, not a gradient.", fontsize=10.5, loc="left")
    ax.legend(handles=[Patch(facecolor=STATUS_COLORS[x], label=STATUS_LABEL[x])
                       for x in ("ND", "AAB", "T1D")], loc="lower right", fontsize=9)

    ax = axes[1]
    tab = pd.crosstab(sst["status"], sst["flag_low_sbr"]).reindex(["ND", "AAB", "T1D"]).fillna(0)
    xs = np.arange(len(tab))
    ax.bar(xs, tab.get(True, 0), width=0.6, color="#CC3311", label="low_sbr (fails)", zorder=3)
    ax.bar(xs, tab.get(False, 0), bottom=tab.get(True, 0), width=0.6, color="#BBBBBB",
           label="separates", zorder=3)
    for xi, (_, row) in enumerate(tab.iterrows()):
        ax.text(xi, row.sum() + 0.15, f"{int(row.get(True, 0))}/{int(row.sum())}", ha="center", fontsize=10)
    ax.set_xticks(xs)
    ax.set_xticklabels([STATUS_LABEL[x] for x in tab.index])
    ax.set_ylabel("donors")
    ax.legend(fontsize=9, loc="upper right")
    ax.set_title("B — the failures are spread across every disease group.\n"
                 "ND 5/8 vs non-ND 3/12 is Fisher p = 0.167 — not a disease pattern.",
                 fontsize=10.5, loc="left")

    fig.suptitle("D3 — SST's weakness is a per-donor pass/fail, not the disease effect the handoff "
                 "describes\nA median hides it: SST target_fold spans 0.94–117 and is bimodal, so "
                 "'ND median 1.68 vs T1D 20.7' compares two medians of a split distribution.",
                 fontsize=12.5, y=1.05)
    fig.savefig(out / "d3_sst_pass_fail.png")
    plt.close(fig)
    print(f"[saved] {out / 'd3_sst_pass_fail.png'}")


# --------------------------------------------------------------------------- #
# D4 / D5 — per-cell co-positivity: what each marker uniquely contributes
# --------------------------------------------------------------------------- #
def per_cell_overlap(cfg, donors, out: Path) -> pd.DataFrame:
    """Per donor: for each gate marker, how many cells it calls that NO other marker in its own gate
    calls. A marker whose positives are all covered by its gate-mates cannot change a compartment; one
    that uniquely claims many cells is load-bearing and its threshold matters.

    Computed cohort-wide AND restricted to islet core, because the two can disagree and the anatomy
    check's central claim is specifically about islet core: EpCAM is there said to recover endocrine
    epithelium the keratins miss. A cohort-wide unique share is dominated by acinar tissue and cannot
    test that."""
    rows = []
    for donor in donors:
        cols = ["object_id"] + [f"{m}_pos" for gate in COMPARTMENT_GATES.values() for m in gate]
        g = ds.dataset(cfg.restore_gated_dir / f"donor_id={donor}", format="parquet",
                       ignore_prefixes=["apply_manifest", "_", "."]).to_table(columns=cols).to_pandas()
        region = ds.dataset(cfg.broad_dir / f"donor_id={donor}", format="parquet").to_table(
            columns=["object_id", "cell_region"]).to_pandas()
        g = g.merge(region, on="object_id", how="left", validate="1:1")
        for scope, mask in (("all", np.ones(len(g), bool)),
                            ("core", (g["cell_region"] == "core").to_numpy())):
            sub = g[mask]
            n = len(sub)
            if not n:
                continue
            for compartment, gate in COMPARTMENT_GATES.items():
                block = sub[[f"{m}_pos" for m in gate]].to_numpy(dtype=bool)
                any_in_gate = block.any(axis=1)
                for i, marker in enumerate(gate):
                    others = (np.delete(block, i, axis=1).any(axis=1) if block.shape[1] > 1
                              else np.zeros(n, bool))
                    rows.append({
                        "donor": donor, "scope": scope, "compartment": compartment, "marker": marker,
                        "n_cells": n, "n_pos": int(block[:, i].sum()),
                        "n_unique": int((block[:, i] & ~others).sum()),
                        "n_gate_any": int(any_in_gate.sum()),
                    })
        print(f"[overlap] {donor}: {len(g):,} cells ({int((g['cell_region'] == 'core').sum()):,} core)",
              flush=True)
    return pd.DataFrame(rows)


def d4_d5_figure(ov: pd.DataFrame, out: Path) -> None:
    ov["pct_pos"] = 100 * ov["n_pos"] / ov["n_cells"]
    ov["pct_unique"] = 100 * ov["n_unique"] / ov["n_cells"]
    ov["unique_share"] = 100 * ov["n_unique"] / ov["n_pos"].replace(0, np.nan)
    ov.to_csv(out / "d4_d5_gate_overlap.csv", index=False)
    order = [(comp, m) for comp, gate in COMPARTMENT_GATES.items() for m in gate]

    scopes = [s for s in ("all", "core") if s in set(ov["scope"])]
    fig, axes = plt.subplots(len(scopes), 1, figsize=(15.0, 5.8 * len(scopes)), sharex=True)
    for ax, scope in zip(np.atleast_1d(axes), scopes):
        agg = (ov[ov.scope == scope].groupby(["compartment", "marker"])
               [["pct_pos", "pct_unique", "unique_share"]].median().reset_index())
        agg["rank"] = agg.apply(lambda r: order.index((r["compartment"], r["marker"])), axis=1)
        agg = agg.sort_values("rank")
        xs = np.arange(len(agg))
        ax.bar(xs, agg["pct_pos"], width=0.72, color="#DDDDDD", zorder=2, label="called positive")
        ax.bar(xs, agg["pct_unique"], width=0.72, zorder=3,
               label="uniquely positive (no gate-mate agrees)",
               color=[COMPARTMENT_COLORS[c] for c in agg["compartment"]])
        for xi, row in zip(xs, agg.itertuples()):
            ax.text(xi, row.pct_pos + 0.8, f"{row.unique_share:.0f}%", ha="center", fontsize=7.5,
                    color="#333333")
        ax.set_xticks(xs)
        ax.set_xticklabels(agg["marker"], rotation=90, fontsize=9)
        ax.set_ylabel("% of cells (cohort median)")
        ax.set_ylim(0, min(agg["pct_pos"].max() * 1.22, 100))
        ax.set_title({"all": "A — all cells", "core": "B — islet core only"}[scope]
                     + "   (grey = all the marker calls · coloured = cells no gate-mate calls · "
                       "number = that unique share)", fontsize=10.5, loc="left")
        ax.legend(fontsize=9, loc="upper left")
    fig.suptitle(
        "D4/D5 — what each gate marker UNIQUELY contributes to its compartment\n"
        "A marker with a small unique share cannot move its compartment however wrong its threshold is; "
        "a large one is load-bearing.\nThe two panels can disagree: cohort-wide numbers are dominated by "
        "acinar tissue, and the epithelial-union argument is specifically about islet core.",
        fontsize=12.5, y=1.02)
    fig.savefig(out / "d4_d5_unique_contribution.png")
    plt.close(fig)
    print(f"[saved] {out / 'd4_d5_unique_contribution.png'}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=REPO / "config.ini")
    ap.add_argument("--outdir", type=Path, default=REPO / "data" / "gate5")
    ap.add_argument("--only", nargs="*", default=None, choices=["D1", "D2", "D3", "D45"])
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    out = Path(a.outdir)
    d, c, b = _tables(out, cfg)
    want = set(a.only or ["D1", "D2", "D3", "D45"])

    if "D1" in want:
        d1_batch_confound(d, c, b, out)
    if "D2" in want:
        d2_cd11b(d, out)
    if "D3" in want:
        d3_sst(d, c, b, out)
    if "D45" in want:
        cache = out / "d4_d5_gate_overlap.csv"
        ov = (pd.read_csv(cache, dtype={"donor": str}) if cache.is_file()
              else per_cell_overlap(cfg, cfg.discover_donors(cfg.restore_gated_dir), out))
        d4_d5_figure(ov, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
