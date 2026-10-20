"""Feasibility figure: what the concordance gate actually sees, on every executed arm.

The 2026-08-21 draft of this figure drew two of the four arms that had been run, six
hand-picked marker->class rows out of the 28 the registry supplies, and one label space.
It therefore could not show the two things that decide this lane:

  B  the argmax throws away what the probability vector knows -- E-cadherin's TOP decile
     is called Epithelial less often than its BOTTOM decile under the 2026-06-15
     checkpoint, an inversion no summary of the calls alone can reveal;
  C  a quarter of that checkpoint's probability mass sits on classes a pancreas cannot
     contain, which needs no marker anchors at all and is the cleanest single
     discriminator between the two checkpoints.

Every panel is drawn against the permutation null (no information) and, where the scale
allows, the 0.920 plane-identity ceiling (the same quantity measured twice), because a
correlation with neither end marked is unreadable.

    python scripts/dct/figure_feasibility.py            # from core/multiplex5000/
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")                       # headless
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from concordance_gate import CEILING, GateInputs, load_anchors, run_gate  # noqa: E402

ROOT = HERE.parents[3]
D = ROOT / "data" / "development" / "deepcell_types_review" / "6539"
CLASSES = tuple(json.loads((HERE / "broad_class_crosswalk.json").read_text())["entries"])
ANCHORS = load_anchors()

ARMS = [("2026-06-15",             "feasibility_tile_no-dapi.parquet",             "#c1442e"),
        ("2026-06-15 +DAPI",       "feasibility_tile_with-dapi.parquet",           "#e08a76"),
        ("2026-06-15 median-tile", "feasibility_tile_no-dapi_median-tile.parquet", "#8c6d1f"),
        ("2026-06-23-ptft",        "feasibility_tile_no-dapi_ptft.parquet",        "#2e6fc1")]


def score(fname):
    d = pd.read_parquet(D / fname)
    P = d[list(CLASSES)].to_numpy(float)
    P = P / P.sum(1, keepdims=True)          # guard float32 round-off from the parquet
    meas = {c: d[c].to_numpy(float) for c in d.columns
            if c not in CLASSES and pd.api.types.is_numeric_dtype(d[c])}
    inp = GateInputs(P, CLASSES, meas)
    return run_gate(inp, restricted=True), run_gate(inp, restricted=False), len(d)


def main():
    scored = {n: score(f) for n, f, _ in ARMS}
    fig, ax = plt.subplots(2, 2, figsize=(16.5, 12.5))
    fig.subplots_adjust(hspace=0.45, wspace=0.3)

    # --- A: every anchor, restricted label space -----------------------------
    a = ax[0, 0]
    order = sorted(ANCHORS, key=lambda t: (t[1], t[0]))
    y = np.arange(len(order))
    for k, (name, _, col) in enumerate(ARMS):
        r = scored[name][0]
        vals = [r.rho.get(m, {}).get(c, np.nan) for m, c in order]
        a.barh(y + (k - 1.5) * 0.2, vals, height=0.19, color=col, label=name)
    null = scored[ARMS[0][0]][0].null["rho_abs"]
    a.axvspan(-null, null, color="0.85", zorder=0)
    a.axvline(0, color="0.2", lw=1)
    a.set_yticks(y)
    a.set_yticklabels([f"{m}  → {c}" for m, c in order], fontsize=7)
    a.invert_yaxis()
    a.set_xlabel("Spearman rho: measured marker vs its compartment's probability")
    a.set_title("A  All 28 registry-derived anchors, restricted label space (26 of 51 classes)\n"
                f"grey band = permutation null (|rho| <= {null:.3f});  "
                f"same-quantity ceiling = {CEILING}", fontsize=10, loc="left")
    a.legend(fontsize=7, loc="lower right", frameon=False)

    # --- B: the argmax throws away what the probability knows ----------------
    b = ax[0, 1]
    show = [m for m, _ in order
            if any(np.isfinite(scored[n][0].coherence.get(m, (np.nan,))[0]) for n, _, _ in ARMS)][:14]
    yy = np.arange(len(show))
    for k, (name, _, col) in enumerate(ARMS):
        r = scored[name][0]
        top = np.array([100 * r.coherence.get(m, (np.nan, np.nan))[0] for m in show])
        bot = np.array([100 * r.coherence.get(m, (np.nan, np.nan))[1] for m in show])
        b.barh(yy + (k - 1.5) * 0.2, top - bot, height=0.19, color=col)
    b.axvline(0, color="0.2", lw=1)
    b.set_yticks(yy); b.set_yticklabels(show, fontsize=7); b.invert_yaxis()
    b.set_xlabel("percentage points: cells called this compartment,\n"
                 "marker's TOP decile minus its BOTTOM decile")
    b.set_title("B  Does the argmax carry what the probability knows?\n"
                "negative = the call is inverted against the marker that defines it",
                fontsize=10, loc="left")

    # --- C: probability on classes a pancreas cannot contain ------------------
    c = ax[1, 0]
    names = [n for n, _, _ in ARMS]
    mass = [scored[n][0].impossible_mass_mean for n in names]
    outw = [100 * scored[n][0].impossible_outweighed for n in names]
    xx = np.arange(len(names))
    # One axis, one unit.  A twin axis here put two different scales under adjacent bars,
    # which is exactly how a reader mistakes a 34 % for a 0.25.
    c.bar(xx - 0.19, 100 * np.array(mass), width=0.36, color=[col for _, _, col in ARMS],
          label="mean probability mass on excluded classes")
    c.bar(xx + 0.19, outw, width=0.36, color="none", edgecolor="0.25", hatch="////",
          label="cells where removed mass beats the survivor")
    c.axhline(10, color="0.3", ls="--", lw=1)
    c.text(len(names) - 0.45, 10.6, "gate bar, both statistics", fontsize=7, ha="right", color="0.3")
    c.set_xticks(xx); c.set_xticklabels(names, fontsize=8, rotation=12, ha="right")
    c.set_ylabel("percent")
    c.legend(fontsize=7, loc="upper right", frameon=False)
    c.set_title("C  Probability spent on classes this tissue cannot contain\n"
                "needs no marker anchors -- the cleanest checkpoint discriminator",
                fontsize=10, loc="left")

    # --- D: does confidence predict agreement? --------------------------------
    d = ax[1, 1]
    for name, _, col in ARMS:
        r = scored[name][0]
        pts = [(lo, hi) for lo, hi in r.confidence_auroc.values()
               if np.isfinite(lo) and np.isfinite(hi)]
        if pts:
            d.scatter(*zip(*pts), s=26, color=col, alpha=0.8, label=name)
    d.plot((0, 1), (0, 1), color="0.3", lw=1, ls="--")
    d.set_xlim(0, 1); d.set_ylim(0, 1)
    d.set_xlabel("anchor AUROC in the LOW max-softmax tercile")
    d.set_ylabel("anchor AUROC in the HIGH tercile")
    d.set_title("D  Does the model's confidence predict its agreement?\n"
                "below the line = high confidence is WORSE; the 2026-08-21 protocol's\n"
                "stratum design assumed everything sits above it", fontsize=10, loc="left")
    d.legend(fontsize=7, loc="lower right", frameon=False)

    for a_ in ax.ravel():
        a_.spines[["top", "right"]].set_visible(False)
    n = scored[ARMS[0][0]][2]
    fig.suptitle("DeepCellTypes feasibility -- donor 6539, 2048-px tile at (16384, 18432), "
                 f"{n:,} cells, stock predict(); all four executed arms", fontsize=12, y=0.98)
    out = D / "feasibility_gate_panels.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print("wrote", out)

    for name, _, _ in ARMS:
        r, u, _ = scored[name]
        print(f"\n{name}")
        for k in r.criteria:
            print(f"   restricted [{'PASS' if r.criteria[k][0] else 'FAIL'}]  "
                  f"full [{'PASS' if u.criteria[k][0] else 'FAIL'}]  {k}")
        print(f"   VERDICT restricted={'PASS' if r.passed else 'FAIL'}  "
              f"full={'PASS' if u.passed else 'FAIL'}")


if __name__ == "__main__":
    main()
