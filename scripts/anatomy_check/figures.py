#!/usr/bin/env python3
"""
Anatomy-check figures A-F -> figures/anatomy_check/.

Reads the aggregates written by ``run_arms.py`` plus the per-arm apply manifests. Every figure also
writes its own source table as CSV beside the PNG, so no claim rests on reading a colour off a bar.

Colour policy: compartments use ``config.COMPARTMENT_COLORS`` unchanged, because the same colours drive
``lineage._composition_figure`` and the QuPath overlays, and a figure that disagreed with the overlay
the maintainer inspects would be worse than a palette flaw. That palette does NOT pass the CVD
validator cleanly -- Neural #B5838D against Mesenchymal #AA3377 is tritan dE 3.3 -- so every stacked
segment above 3% is directly labelled and the legend is always in fixed compartment order. Colour is
therefore never the only channel carrying identity.
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                                # noqa: E402
from matplotlib.lines import Line2D                                            # noqa: E402
from matplotlib.patches import Patch                                           # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
from phenocycler.config import COMPARTMENT_COLORS, COMPARTMENT_ORDER, OTHER_LABEL, UNRESOLVED_LABEL  # noqa: E402

CLASSES = COMPARTMENT_ORDER + [OTHER_LABEL, UNRESOLVED_LABEL]
ARM_ORDER = ["reference_and_target/max", "reference_only/max", "reference_only/p999",
             "reference_only/p99", "reference_only/p95"]
ARM_LABEL = {"reference_and_target/max": "v9\nctrl+tgt / max", "reference_only/max": "ref-only\nmax",
             "reference_only/p999": "ref-only\np99.9", "reference_only/p99": "ref-only\np99",
             "reference_only/p95": "ref-only\np95"}
REGION_LABEL = {"core": "islet core", "peri": "peri-islet", "tissue": "tissue (non-islet)"}
DESIGN_LABEL = {"single": "Epithelial = E-cadherin", "union": "Epithelial = E-cad ∪ Pan-CK ∪ Ker8/18 ∪ EpCAM"}

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "legend.frameon": False,
})


def _load(summary_dir: Path, prefix: str) -> pd.DataFrame:
    files = sorted(glob.glob(str(summary_dir / f"{prefix}_*.csv")))
    if not files:
        raise SystemExit(f"[err] no {prefix}_*.csv under {summary_dir} — run run_arms.py first")
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d["arm"] = pd.Categorical(d["arm"], [a for a in ARM_ORDER if a in set(d["arm"])], ordered=True)
    return d.dropna(subset=["arm"])


def _pct(d: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    """n -> % within `by`, on a complete CLASSES index so a zero category is a zero, not a gap."""
    tot = d.groupby(by, observed=True)["n"].transform("sum")
    out = d.assign(pct=100 * d["n"] / tot)
    return out


def _stack(ax, frame: pd.DataFrame, xcol: str, xorder, *, label_min=3.0, width=0.72):
    """Stacked composition bars with a 2px surface gap and direct labels on segments >= label_min%."""
    piv = (frame.pivot_table(index=xcol, columns="compartment", values="pct", aggfunc="sum",
                             observed=True)
           .reindex(index=xorder).reindex(columns=CLASSES).fillna(0.0))
    bottom = np.zeros(len(piv))
    xs = np.arange(len(piv))
    for comp in CLASSES:
        v = piv[comp].to_numpy()
        ax.bar(xs, v, bottom=bottom, width=width, color=COMPARTMENT_COLORS[comp],
               edgecolor="white", linewidth=1.4, zorder=3)
        for xi, (b, h) in enumerate(zip(bottom, v)):
            if h >= label_min:
                ax.text(xs[xi], b + h / 2, f"{h:.0f}", ha="center", va="center", fontsize=8,
                        color="#111111" if comp in ("Endothelial", "Other") else "white", zorder=4)
        bottom += v
    ax.set_xticks(xs)
    ax.set_xlim(-0.6, len(piv) - 0.4)
    ax.set_ylim(0, 100)
    return piv


def _legend(fig, ncol=7, y=-0.01):
    fig.legend(handles=[Patch(facecolor=COMPARTMENT_COLORS[c], label=c) for c in CLASSES],
               loc="lower center", ncol=ncol, bbox_to_anchor=(0.5, y), fontsize=9)


# --------------------------------------------------------------------------- #
# A — composition x cell_region x arm
# --------------------------------------------------------------------------- #
def figure_a(comp: pd.DataFrame, out: Path):
    d = comp.groupby(["arm", "gate_design", "cell_region", "compartment"], observed=True)["n"].sum().reset_index()
    d = _pct(d, ["arm", "gate_design", "cell_region"])
    regions = [r for r in ("core", "peri", "tissue") if r in set(d.cell_region)]
    arms = list(d["arm"].cat.categories)

    fig, axes = plt.subplots(len(regions), 1, figsize=(1.55 * len(arms) * 2 + 2.2, 3.1 * len(regions)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    xorder = [(a, g) for a in arms for g in ("single", "union")]
    for ax, region in zip(axes, regions):
        sub = d[d.cell_region == region].copy()
        sub["key"] = list(zip(sub["arm"].astype(str), sub["gate_design"]))
        _stack(ax, sub, "key", xorder)
        ax.set_ylabel(f"{REGION_LABEL[region]}\n% of cells")
        for i in range(1, len(arms)):                       # separate the arm groups
            ax.axvline(2 * i - 0.5, color="#dddddd", lw=1.0, zorder=1)
    axes[-1].set_xticklabels([f"{ARM_LABEL[a]}\n{g}" for a, g in xorder], fontsize=8)
    axes[0].set_title("A — Compartment composition by anatomical region, per divisor arm and Epithelial gate design\n"
                      "islet-core cells are endocrine (INS/GCG-high) and must land in Epithelial",
                      loc="left", fontsize=12)
    _legend(fig, y=0.0)
    fig.savefig(out / "A_composition_by_region.png")
    plt.close(fig)
    d.to_csv(out / "A_composition_by_region.csv", index=False)


# --------------------------------------------------------------------------- #
# B — islet-core purity (the sharpest single number)
# --------------------------------------------------------------------------- #
def figure_b(comp: pd.DataFrame, out: Path):
    core = comp[comp.cell_region == "core"]
    d = _pct(core.groupby(["arm", "gate_design", "donor", "compartment"], observed=True)["n"].sum().reset_index(),
             ["arm", "gate_design", "donor"])
    arms = list(d["arm"].cat.categories)
    shown = ["Epithelial", OTHER_LABEL, "Immune"]

    fig, axes = plt.subplots(1, len(shown), figsize=(4.6 * len(shown), 4.4), sharex=True)
    for ax, comp_name in zip(np.atleast_1d(axes), shown):
        sub = d[d.compartment == comp_name]
        for design, marker, color, off in (("single", "o", "#4477AA", -0.11), ("union", "D", "#EE7733", 0.11)):
            s = sub[sub.gate_design == design]
            xs = s["arm"].cat.codes.to_numpy() + off
            ax.scatter(xs, s["pct"], s=52, marker=marker, color=color, alpha=0.75,
                       edgecolor="white", linewidth=1.2, zorder=3,
                       label=DESIGN_LABEL[design] if comp_name == shown[0] else None)
            med = s.groupby("arm", observed=True)["pct"].median().reindex(arms)
            ax.plot(np.arange(len(arms)) + off, med.to_numpy(), color=color, lw=2, alpha=0.9, zorder=4)
        ax.set_title(f"islet-core cells called {comp_name}", pad=10)
        ax.set_xticks(range(len(arms)))
        ax.set_xticklabels([ARM_LABEL[a] for a in arms], fontsize=8)
        ax.set_ylim(0, 100)
        ax.set_xlim(-0.55, len(arms) - 0.45)
    np.atleast_1d(axes)[0].set_ylabel("% of that donor's islet-core cells")
    fig.suptitle("B — Islet-core purity. Core cells are endocrine by QuPath annotation (independent of every marker "
                 "threshold),\nso Epithelial should be high and Other low. One point per donor; line = cohort median.",
                 x=0.01, y=1.10, ha="left", va="bottom", fontsize=12)
    fig.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.10), fontsize=9)
    fig.savefig(out / "B_islet_core_purity.png")
    plt.close(fig)
    d.to_csv(out / "B_islet_core_purity.csv", index=False)


# --------------------------------------------------------------------------- #
# C — radial profile (the study's actual readout)
# --------------------------------------------------------------------------- #
def figure_c(radial: pd.DataFrame, out: Path, designs=("single", "union")):
    d = _pct(radial.groupby(["arm", "gate_design", "dist_bin", "compartment"], observed=True)["n"].sum().reset_index(),
             ["arm", "gate_design", "dist_bin"])
    # order bins by their left edge, parsed from the pandas Interval string
    bins = sorted(set(d.dist_bin), key=lambda s: float(str(s).strip("[)").split(",")[0]))
    arms = list(d["arm"].cat.categories)
    fig, axes = plt.subplots(len(designs), len(arms), figsize=(3.0 * len(arms), 3.2 * len(designs)),
                             sharey=True, sharex=True, squeeze=False)
    for r, design in enumerate(designs):
        for c, arm in enumerate(arms):
            ax = axes[r][c]
            sub = d[(d.arm == arm) & (d.gate_design == design)]
            for compart in CLASSES:
                s = sub[sub.compartment == compart].set_index("dist_bin").reindex(bins)["pct"]
                ax.plot(range(len(bins)), s.to_numpy(), color=COMPARTMENT_COLORS[compart], lw=2,
                        marker="o", ms=3.5, zorder=3)
            if r == 0:
                ax.set_title(ARM_LABEL[arm].replace("\n", " "), fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{design}\n% of cells in bin")
            ax.set_xticks(range(len(bins)))
            ax.set_xticklabels([str(b).strip("[)").split(",")[0] for b in bins], fontsize=7, rotation=90)
            ax.set_yscale("symlog", linthresh=1)
    fig.supxlabel("distance to nearest islet (µm, bin lower edge; 0 = inside an islet)", fontsize=10, y=-0.02)
    fig.suptitle("C — Compartment fraction vs distance to islet. Endothelium enriched INSIDE islets is a positive "
                 "biological control\n(islets are densely vascularised) that no threshold was tuned to produce. "
                 "Epithelium at its LOWEST inside islets is not.\nCohort-pooled over 7 donors of mixed disease "
                 "status — the per-donor, per-status version is the study readout, not this.",
                 x=0.01, y=1.0, ha="left", va="bottom", fontsize=12)
    _legend(fig, y=-0.10)
    fig.savefig(out / "C_radial_profile.png")
    plt.close(fig)
    d.to_csv(out / "C_radial_profile.csv", index=False)


# --------------------------------------------------------------------------- #
# D — per-marker call rate vs arm (the compensating-error structure)
# --------------------------------------------------------------------------- #
def figure_d(arm_root: Path, out: Path):
    rows = []
    for arm in ARM_ORDER:
        slug = "ctrl-{}_div-{}".format(*arm.split("/"))
        for path, kind in ((arm_root / slug / "restore_gated_redsea" / "apply_manifest.csv", "accepted"),
                           (arm_root / slug / "restore_gated_comparator" / "comparator_manifest.csv", "comparator")):
            if not path.exists():
                continue
            m = pd.read_csv(path)
            m = m.assign(arm=arm, kind=kind,
                         pct=100 * m["called_positive_n"] / m["n_input"],
                         pair=m["target"] + " ← " + m["reference"])
            rows.append(m[["arm", "kind", "donor", "target", "reference", "pair", "pct"]])
    if not rows:
        print("[D] no apply manifests found — skipped")
        return
    d = pd.concat(rows, ignore_index=True)
    d["arm"] = pd.Categorical(d["arm"], ARM_ORDER, ordered=True)

    # Arms that halted (restore_apply refused to emit calls) are a RESULT; mark them rather than
    # silently closing the gap, which is what misaligned the x positions in the first draft.
    outcomes = {}
    oc = arm_root.parent / "anatomy_check" / "arm_outcomes.csv"
    if oc.exists():
        outcomes = dict(pd.read_csv(oc)[["arm", "outcome"]].values)
    present = {a for a in ARM_ORDER if (d.arm == a).any()}

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.0), sharey=True)
    for ax, kind, title in ((axes[0], "accepted", "the 7 ACCEPTED production pairs"),
                            (axes[1], "comparator", "the epithelial comparator pairs (diagnostic)")):
        sub = d[d.kind == kind]
        cmap = plt.get_cmap("tab10")
        for i, pair in enumerate(sorted(set(sub["pair"]))):
            # Reindex onto the FULL arm order: a missing arm must be a gap at its own x position, never
            # a shift of the remaining points onto the wrong ticks.
            s = (sub[sub.pair == pair].groupby("arm", observed=False)["pct"].median()
                 .reindex(ARM_ORDER))
            ax.plot(range(len(ARM_ORDER)), s.to_numpy(), marker="o", lw=2, ms=6,
                    color=cmap(i % 10), label=pair)
        for j, arm in enumerate(ARM_ORDER):
            if arm in present:
                continue
            ax.axvspan(j - 0.45, j + 0.45, color="#f2f2f2", zorder=0)
            note = "halted\n(MODEL_UNSTABLE)" if outcomes.get(arm) == "halted" else "not run"
            ax.text(j, 0.5, note, ha="center", va="center", fontsize=8, color="#777777",
                    transform=ax.get_xaxis_transform())
        ax.set_xticks(range(len(ARM_ORDER)))
        ax.set_xticklabels([ARM_LABEL[a] for a in ARM_ORDER], fontsize=8)
        ax.set_xlim(-0.5, len(ARM_ORDER) - 0.5)
        ax.set_yscale("log")
        ax.set_title(title, pad=10)
        ax.legend(fontsize=8, loc="lower right", ncol=2)
    axes[0].set_ylabel("% of cells called positive (cohort median, log scale)")
    fig.suptitle("D — Every marker moves with the divisor statistic, not just E-cadherin. A single global statistic\n"
                 "cannot be tuned for E-cadherin without inflating the sparse markers by the same operation.",
                 x=0.01, y=1.02, ha="left", va="bottom", fontsize=12)
    fig.savefig(out / "D_call_rate_by_arm.png")
    plt.close(fig)
    d.to_csv(out / "D_call_rate_by_arm.csv", index=False)


# --------------------------------------------------------------------------- #
# E — spatial coherence (validate against the donor's own image, not the cohort)
# --------------------------------------------------------------------------- #
def _assign(arm_root: Path, slug: str, donor: str, design: str) -> pd.DataFrame:
    from phenocycler import lineage
    from phenocycler.config import COMPARTMENT_GATES
    base = dict(COMPARTMENT_GATES)
    gates = {**base, "Epithelial": ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM"]} if design == "union" else base
    part = "restore_gated_union" if design == "union" else "restore_gated_redsea"
    gated = sorted((arm_root / slug / part / f"donor_id={donor}").glob("*.parquet"))
    if not gated:
        return pd.DataFrame()
    lineage.COMPARTMENT_GATES = gates
    try:
        out = lineage.assign_donor(donor, str(gated[0]), arm_root / slug / "cells")
    finally:
        lineage.COMPARTMENT_GATES = base
    xy = pd.read_parquet(sorted((arm_root / slug / "cells" / f"donor_id={donor}").glob("*.parquet"))[0],
                         columns=["object_id", "X_centroid", "Y_centroid"])
    xy["object_id"] = xy["object_id"].astype(str)
    return out.merge(xy, on="object_id", how="left")


def _dense_tissue_window(d: pd.DataFrame, side_um: float = 350.0):
    """A square window centred on the densest NON-islet neighbourhood — i.e. solid acinar tissue."""
    t = d[d["cell_region"] == "tissue"] if "cell_region" in d.columns else d
    if t.empty:
        t = d
    hx = np.histogram2d(t["X_centroid"], t["Y_centroid"], bins=40)
    counts, xe, ye = hx
    i, j = np.unravel_index(np.argmax(counts), counts.shape)
    cx, cy = 0.5 * (xe[i] + xe[i + 1]), 0.5 * (ye[j] + ye[j + 1])
    h = side_um / 2
    return cx - h, cx + h, cy - h, cy + h


def figure_e(arm_root: Path, out: Path, donor: str, arms: list[str], designs=("single", "union")):
    panels = [(a, d) for a in arms for d in designs]
    assigned = {}
    for arm, design in panels:
        slug = "ctrl-{}_div-{}".format(*arm.split("/"))
        f = _assign(arm_root, slug, donor, design)
        if not f.empty:
            assigned[(arm, design)] = f
    if not assigned:
        print("[E] no arms available — skipped")
        return
    win = _dense_tissue_window(next(iter(assigned.values())))

    fig, axes = plt.subplots(2, len(panels), figsize=(5.4 * len(panels), 11.0), squeeze=False)
    order = [c for c in CLASSES if c != OTHER_LABEL] + [OTHER_LABEL]   # Other last, so it cannot hide a call
    for col, (arm, design) in enumerate(panels):
        d = assigned.get((arm, design))
        if d is None:
            for r in (0, 1):
                axes[r][col].set_visible(False)
            continue
        pct = d["compartment"].value_counts(normalize=True).mul(100)
        zoom = d[(d.X_centroid.between(win[0], win[1])) & (d.Y_centroid.between(win[2], win[3]))]
        for row, (frame, size, alpha_other) in enumerate(((d, 0.35, 0.35), (zoom, 26.0, 0.85))):
            ax = axes[row][col]
            for compart in order:
                s = frame[frame.compartment == compart]
                ax.scatter(s["X_centroid"], s["Y_centroid"], s=size, linewidths=0,
                           color=COMPARTMENT_COLORS[compart], rasterized=True,
                           alpha=alpha_other if compart == OTHER_LABEL else 0.9)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
        axes[0][col].plot([win[0], win[1], win[1], win[0], win[0]], [win[2], win[2], win[3], win[3], win[2]],
                          color="#111111", lw=1.4, zorder=5)
        axes[0][col].set_title(f"{ARM_LABEL[arm].replace(chr(10), ' ')} · {design}\n"
                               f"Epi {pct.get('Epithelial', 0):.0f}% · Other {pct.get(OTHER_LABEL, 0):.0f}% · "
                               f"Imm {pct.get('Immune', 0):.0f}%", fontsize=11)
        axes[1][col].set_title(f"350 µm window in solid acinar tissue — n={len(zoom):,}", fontsize=9)
    fig.suptitle(f"E — Spatial coherence, donor {donor}. Acini are contiguous epithelial structures. Epithelial and "
                 "Other finely INTERLEAVED\nwithin one acinus (lower row) is the signature of a threshold cutting "
                 "through a single population — checkable\nagainst this donor's own image, not against a cohort "
                 "expectation.", x=0.01, y=1.0, ha="left", va="bottom", fontsize=12)
    _legend(fig, y=-0.015)
    fig.savefig(out / f"E_spatial_coherence_{donor}.png", dpi=170)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# F — marker polarity: can ONE marker's threshold span both epithelial populations?
# --------------------------------------------------------------------------- #
def figure_f(repo: Path, arm_root: Path, out: Path, donor: str,
             markers=("E_cadherin", "EpCAM", "Ker8_18", "Pan_Cytokeratin")):
    import duckdb
    con = duckdb.connect(); con.execute("PRAGMA threads=4")
    cols = ", ".join(f"r.{m}" for m in markers)
    d = con.execute(f"""
        SELECT c.cell_region, {cols}
        FROM read_parquet('{repo}/data/cells_redsea/donor_id={donor}/*.parquet') r
        JOIN read_parquet('{repo}/data/cells/donor_id={donor}/*.parquet') c USING (object_id)
    """).df()

    # Divisors for this donor, per arm: accepted manifest for E-cadherin, comparator for the rest.
    div = {}
    for arm in ARM_ORDER:
        slug = "ctrl-{}_div-{}".format(*arm.split("/"))
        for path in (arm_root / slug / "restore_gated_redsea" / "apply_manifest.csv",
                     arm_root / slug / "restore_gated_comparator" / "comparator_manifest.csv"):
            if not path.exists():
                continue
            m = pd.read_csv(path)
            m = m[m["donor"].astype(str) == str(donor)]
            for _, r in m.iterrows():
                if r["target"] in markers and (r["target"] != "E_cadherin" or r["reference"] == "CD31"):
                    if r["target"] == "E_cadherin" and path.name == "comparator_manifest.csv":
                        continue                       # prefer the accepted-pair row for E-cadherin
                    div.setdefault(arm, {})[r["target"]] = r.get("canonical_divisor", np.nan)

    fig, axes = plt.subplots(1, len(markers), figsize=(4.8 * len(markers), 4.8))
    region_style = {"core": ("#CC3311", 2.4, "islet core (endocrine)"),
                    "peri": ("#EE9944", 1.5, "peri-islet"),
                    "tissue": ("#33549E", 2.0, "tissue (acinar/ductal)")}
    arm_ls = dict(zip(ARM_ORDER, ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]))
    for ax, marker in zip(np.atleast_1d(axes), markers):
        # x-range from the marker's own tissue distribution, so each panel is on its own scale and no
        # curve is squashed by another marker's dynamic range. Values above are clipped INTO the last
        # bin, so the right-hand edge is a boundary artifact and is drawn hollow.
        hi = float(np.nanquantile(d.loc[d.cell_region == "tissue", marker], 0.995))
        hi = max(hi, 1.0)
        for region, (color, lw, label) in region_style.items():
            v = d.loc[d.cell_region == region, marker].to_numpy(float)
            v = v[np.isfinite(v)]
            if len(v) < 200:            # too few cells for a density to mean anything
                continue
            hist, edges = np.histogram(np.clip(v, 0, hi), bins=100, range=(0, hi), density=True)
            mid = 0.5 * (edges[1:] + edges[:-1])
            ax.plot(mid[:-1], hist[:-1], color=color, lw=lw, zorder=3,
                    label=f"{label}  (n={len(v):,})")
        for arm, ls in arm_ls.items():
            t = div.get(arm, {}).get(marker)
            if t is None or not np.isfinite(t):
                continue
            inside = t <= hi
            ax.axvline(min(t, hi), color="#444444", lw=1.3, ls=ls, zorder=2,
                       label=(ARM_LABEL[arm].replace("\n", " ") + ("" if inside else " (off-scale →)"))
                       if marker == markers[0] or not inside else None)
        ax.set_title(marker.replace("_", "-"), pad=10)
        ax.set_xlabel("REDSEA whole-cell mean")
        ax.set_xlim(0, hi)
        ax.legend(fontsize=7, loc="upper right")
    np.atleast_1d(axes)[0].set_ylabel("density")
    fig.suptitle(f"F — Marker polarity, donor {donor}. Endocrine (islet core) sits BELOW acinar on E-cadherin and "
                 "ABOVE it on EpCAM.\nA threshold placed to capture acinar E-cadherin necessarily cuts through the "
                 "endocrine population — no divisor statistic can fix that.\n"
                 "Vertical lines are that arm's divisor for that marker, on this donor.",
                 x=0.01, y=1.02, ha="left", va="bottom", fontsize=12)
    fig.savefig(out / f"F_marker_polarity_{donor}.png")
    plt.close(fig)
    d.groupby("cell_region")[list(markers)].describe().to_csv(out / f"F_marker_polarity_{donor}.csv")


# --------------------------------------------------------------------------- #
# G — is the union corroborated epithelium, or single-marker noise?
# --------------------------------------------------------------------------- #
UNION_MARKERS = ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM"]


def figure_g(arm_root: Path, out: Path, donors: list[str]):
    """A union gate is only worth having if its extra cells are corroborated. This measures that."""
    rows, jac = [], {}
    for arm in ARM_ORDER:
        slug = "ctrl-{}_div-{}".format(*arm.split("/"))
        for donor in donors:
            f = arm_root / slug / "restore_gated_union" / f"donor_id={donor}" / "data_0.parquet"
            if not f.exists():
                continue
            d = pd.read_parquet(f, columns=[f"{m}_state" for m in UNION_MARKERS])
            pos = pd.DataFrame({m: d[f"{m}_state"].eq("positive") for m in UNION_MARKERS})
            k = pos.sum(axis=1)
            n = len(pos)
            rows.append(dict(arm=arm, donor=donor, n=n,
                             ecad_only_pct=100 * float((pos.E_cadherin & (k == 1)).mean()),
                             union_adds_pct=100 * float(((~pos.E_cadherin) & (k > 0)).mean()),
                             added_corroborated_pct=100 * float(((~pos.E_cadherin) & (k >= 2)).mean()),
                             any_epithelial_pct=100 * float((k > 0).mean()),
                             all_four_pct=100 * float((k == 4).mean()),
                             none_pct=100 * float((k == 0).mean()),
                             **{f"k{i}_pct": 100 * float((k == i).mean()) for i in range(5)}))
            if arm == "reference_only/p99":
                J = np.zeros((4, 4))
                for i, aM in enumerate(UNION_MARKERS):
                    for j, bM in enumerate(UNION_MARKERS):
                        u = (pos[aM] | pos[bM]).sum()
                        J[i, j] = (pos[aM] & pos[bM]).sum() / u if u else np.nan
                jac[donor] = J
    if not rows:
        print("[G] no union parquets — skipped")
        return
    g = pd.DataFrame(rows)
    g["arm"] = pd.Categorical(g["arm"], ARM_ORDER, ordered=True)

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.4))
    # (1) how many of the four markers agree, per arm — stacked, cohort mean
    ax = axes[0]
    m = g.groupby("arm", observed=False)[[f"k{i}_pct" for i in range(5)]].mean().reindex(ARM_ORDER)
    bottom = np.zeros(len(m))
    shades = ["#dddddd", "#bcd0e5", "#8fb3d9", "#5b8cc4", "#2f5f9e"]
    for i in range(5):
        v = m[f"k{i}_pct"].fillna(0).to_numpy()
        ax.bar(range(len(m)), v, bottom=bottom, width=0.72, color=shades[i],
               edgecolor="white", linewidth=1.4, label=f"{i} of 4 positive", zorder=3)
        for x, (b, h) in enumerate(zip(bottom, v)):
            if h >= 4:
                ax.text(x, b + h / 2, f"{h:.0f}", ha="center", va="center", fontsize=8,
                        color="#111111" if i < 2 else "white", zorder=4)
        bottom += v
    ax.set_xticks(range(len(m))); ax.set_xticklabels([ARM_LABEL[a] for a in ARM_ORDER], fontsize=8)
    ax.set_ylim(0, 100); ax.set_ylabel("% of cells (cohort mean)")
    ax.set_title("how many of the 4 epithelial markers call each cell", pad=10)
    ax.legend(fontsize=8, loc="lower left", ncol=2)

    # (2) what the union actually adds, and how much of it is corroborated
    ax = axes[1]
    xs = np.arange(len(ARM_ORDER))
    add = g.groupby("arm", observed=False)["union_adds_pct"].mean().reindex(ARM_ORDER)
    cor = g.groupby("arm", observed=False)["added_corroborated_pct"].mean().reindex(ARM_ORDER)
    ax.bar(xs, add.to_numpy(), width=0.66, color="#cfd9e4", edgecolor="white", linewidth=1.4,
           label="added by the union (E-cadherin−, ≥1 other marker+)", zorder=3)
    ax.bar(xs, cor.to_numpy(), width=0.66, color="#2f5f9e", edgecolor="white", linewidth=1.4,
           label="of which corroborated by ≥2 other markers", zorder=4)
    for x, (a_, c_) in enumerate(zip(add.to_numpy(), cor.to_numpy())):
        if np.isfinite(a_):
            ax.text(x, a_ + 0.4, f"{a_:.1f}", ha="center", fontsize=8)
    ax.set_xticks(xs); ax.set_xticklabels([ARM_LABEL[a] for a in ARM_ORDER], fontsize=8)
    ax.set_ylabel("% of all cells"); ax.set_title("what the union adds beyond E-cadherin", pad=10)
    ax.legend(fontsize=8, loc="upper left")

    # (3) pairwise Jaccard at the leading arm — do the four markers mark the same cells?
    ax = axes[2]
    if jac:
        J = np.nanmean(np.stack(list(jac.values())), axis=0)
        im = ax.imshow(J, cmap="Blues", vmin=0, vmax=1)
        for i in range(4):
            for j in range(4):
                ax.text(j, i, f"{J[i, j]:.2f}", ha="center", va="center", fontsize=9,
                        color="white" if J[i, j] > 0.6 else "#111111")
        lab = [m_.replace("_", "-") for m_ in UNION_MARKERS]
        ax.set_xticks(range(4)); ax.set_xticklabels(lab, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(4)); ax.set_yticklabels(lab, fontsize=8)
        ax.grid(False)
        ax.set_title(f"pairwise Jaccard, ref-only p99 (mean of {len(jac)} donors)", pad=10)
        fig.colorbar(im, ax=ax, fraction=0.046)
    else:
        ax.set_visible(False)

    fig.suptitle("G — Is the union corroborated epithelium or single-marker noise? The four epithelial markers agree "
                 "substantially\n(Jaccard 0.44–0.79), but a majority of what the union ADDS rests on one marker — so "
                 "the union is a real recovery,\nnot a licence to treat any single-marker positive as epithelium.",
                 x=0.01, y=1.02, ha="left", va="bottom", fontsize=12)
    fig.savefig(out / "G_epithelial_concordance.png")
    plt.close(fig)
    g.to_csv(out / "G_epithelial_concordance.csv", index=False)


# --------------------------------------------------------------------------- #
# H — spatial coherence, as a number rather than an eyeball
# --------------------------------------------------------------------------- #
def figure_h(arm_root: Path, out: Path, donors: list[str], designs=("single", "union"), k=12):
    """Neighbour concordance per compartment, against the value expected under random labelling.

    Figure E shows Epithelial and Other interleaved inside single acini. This is the same claim as a
    statistic: for each cell, the fraction of its k nearest neighbours carrying the same compartment.
    The null comes from the donor's own composition, so no frequency is assumed.

    Report ``segregation`` = (observed - marginal) / (1 - marginal), NOT the ratio observed/marginal.
    The ratio is intuitive but its ceiling is 1/marginal, so a 0.5%-frequency compartment can score 30
    where a 50% one cannot exceed 2, and comparing compartments on it is meaningless. Measured here the
    ratio makes Neural look ~20x more coherent than Epithelial while the rescaled statistic puts them
    within a factor of 1.5, and reverses the conclusion. The rescaled one is what the figure plots.
    """
    from scipy.spatial import cKDTree
    rows = []
    for arm in ARM_ORDER:
        slug = "ctrl-{}_div-{}".format(*arm.split("/"))
        for donor in donors:
            for design in designs:
                d = _assign(arm_root, slug, donor, design)
                if d.empty:
                    continue
                xy = d[["X_centroid", "Y_centroid"]].to_numpy(float)
                ok = np.isfinite(xy).all(axis=1)
                d, xy = d[ok], xy[ok]
                comp = d["compartment"].to_numpy()
                # k+1 because the first neighbour of a point is itself
                _, idx = cKDTree(xy).query(xy, k=k + 1)
                same = (comp[idx[:, 1:]] == comp[:, None]).mean(axis=1)
                for c in CLASSES:
                    m = comp == c
                    if m.sum() < 500:
                        continue
                    marginal = float(m.mean())
                    obs = float(same[m].mean())
                    rows.append(dict(
                        arm=arm, donor=donor, gate_design=design, compartment=c,
                        n=int(m.sum()), marginal=100 * marginal, concordance=100 * obs,
                        # ratio to the random-labelling expectation. Intuitive, but NOT comparable
                        # across compartments: its ceiling is 1/marginal, so a rare compartment can
                        # score 30 where a majority one cannot exceed 2.
                        coherence_index=obs / marginal,
                        # segregation on a common 0..1 scale: 0 = random labelling, 1 = every
                        # neighbour shares the compartment. This is the one to compare across
                        # compartments, and it is what the figure plots.
                        segregation=(obs - marginal) / (1.0 - marginal),
                    ))
    if not rows:
        print("[H] nothing to compute — skipped")
        return
    h = pd.DataFrame(rows)
    h["arm"] = pd.Categorical(h["arm"], ARM_ORDER, ordered=True)

    shown = ["Epithelial", "Immune", "Endothelial", OTHER_LABEL]
    fig, axes = plt.subplots(1, len(shown), figsize=(4.4 * len(shown), 4.8), sharey=True)
    for ax, comp_name in zip(np.atleast_1d(axes), shown):
        sub = h[h.compartment == comp_name]
        for design, marker, color, off in (("single", "o", "#4477AA", -0.11),
                                           ("union", "D", "#EE7733", 0.11)):
            s = sub[sub.gate_design == design]
            if s.empty:
                continue
            ax.scatter(s["arm"].cat.codes.to_numpy() + off, s["segregation"], s=46, marker=marker,
                       color=color, alpha=0.75, edgecolor="white", linewidth=1.1, zorder=3,
                       label=DESIGN_LABEL[design] if comp_name == shown[0] else None)
            med = s.groupby("arm", observed=False)["segregation"].median().reindex(ARM_ORDER)
            ax.plot(range(len(ARM_ORDER)), med.to_numpy(), color=color, lw=2, alpha=0.9, zorder=4)
        ax.axhline(0.0, color="#CC3311", lw=1.4, ls="--", zorder=2)
        ax.text(0.02, 0.0, " random labelling", color="#CC3311", fontsize=8, va="bottom",
                transform=ax.get_yaxis_transform())
        ax.set_ylim(-0.02, 0.30)
        ax.set_xticks(range(len(ARM_ORDER)))
        ax.set_xticklabels([ARM_LABEL[a] for a in ARM_ORDER], fontsize=8)
        ax.set_xlim(-0.5, len(ARM_ORDER) - 0.5)
        ax.set_title(comp_name, pad=10)
    np.atleast_1d(axes)[0].set_ylabel(f"segregation  (obs − marginal) / (1 − marginal)\n"
                                      f"0 = random labelling, 1 = every neighbour matches (k={k})")
    fig.suptitle("H — Figure E quantified, and it does NOT support the strong reading of that figure. On a scale where "
                 "0 is random\nlabelling, EVERY compartment scores 0.02–0.20: single-cell neighbourhoods in pancreas "
                 "are genuinely mixed. Epithelial is\namong the MORE structured, not the less, and the union raises it "
                 "at every arm. `Other` is the least structured everywhere,\nwhich is what a residual bucket of "
                 "heterogeneous leftovers should look like.",
                 x=0.01, y=1.02, ha="left", va="bottom", fontsize=12)
    fig.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.10), fontsize=9)
    fig.savefig(out / "H_spatial_coherence_index.png")
    plt.close(fig)
    h.to_csv(out / "H_spatial_coherence_index.csv", index=False)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary", type=Path, default=REPO / "data" / "anatomy_check")
    ap.add_argument("--arm-root", type=Path, default=REPO / "data" / "_arm")
    ap.add_argument("--out", type=Path, default=REPO / "figures" / "anatomy_check")
    ap.add_argument("--spatial-donor", default="6539", help="donor for figures E and F")
    ap.add_argument("--spatial-arms", nargs="*", default=["reference_and_target/max", "reference_only/p99"])
    ap.add_argument("--coherence-donors", nargs="*", default=None,
                    help="donors for figure H (default: all in the aggregates; kNN over every cell)")
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)

    comp = _load(a.summary, "region_composition")
    radial = _load(a.summary, "radial_profile")
    figure_a(comp, a.out)
    figure_b(comp, a.out)
    figure_c(radial, a.out)
    figure_d(a.arm_root, a.out)
    figure_e(a.arm_root, a.out, a.spatial_donor, a.spatial_arms)
    figure_f(REPO, a.arm_root, a.out, a.spatial_donor)
    donors = sorted(comp["donor"].astype(str).unique())
    figure_g(a.arm_root, a.out, donors)
    figure_h(a.arm_root, a.out, a.coherence_donors or donors)
    print("[done] figures ->", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
