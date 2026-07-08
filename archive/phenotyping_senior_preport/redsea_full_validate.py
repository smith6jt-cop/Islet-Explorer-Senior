#!/usr/bin/env python3
"""
Validate pixel-level REDSEA with ONE figure (donor smoke-test).

Spillover makes cells falsely positive for markers of OTHER lineages, creating impossible
cross-compartment co-expression. A correct validation must use MUTUALLY-EXCLUSIVE broad
compartments -- NOT the 9 individual lineage markers, because several markers are not
lineage-exclusive (Vimentin is pan-mesenchymal: endothelium/muscle/immune/fibroblasts are
all Vimentin+; CD31+SMA+Vimentin = a real blood-vessel wall, not an artifact). Using the 9
markers over-counts "pan-positive" cells with real vascular/stromal biology.

Broad mutually-exclusive compartments (a cell should be ~exactly one):
    Epithelial      = Pan_Cytokeratin
    Immune          = CD3e | CD20 | CD163
    Endocrine       = INS  | GCG
    Stromal/vascular= CD31 | SMA  | Vimentin
Genuine impossible artifacts:
    cross-compartment >= 2/3      (e.g. epithelial + immune)
    multi-immune (>= 2 of T/B/Mac)
    Beta+Alpha in one cell        (mostly SEGMENTATION doublets in dense islets, not spillover)

Gating is held FIXED (raw RESTORE SSC thresholds applied to raw AND corrected) so we measure
REDSEA's effect on the VALUES, not RESTORE threshold drift. Reads saved intermediates +
mask, so it sweeps the subtraction strength alpha cheaply for the tradeoff curve.

  python scripts/senior/redsea_full_validate.py --donor 6539 --alpha 1.0
"""
from __future__ import annotations
import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import redsea_full as R  # noqa: E402

BASE = R.BASE
MARK = ["Pan_Cytokeratin", "CD3e", "CD20", "CD163", "INS", "GCG", "CD31", "SMA", "Vimentin"]
COMP = {"Epithelial": ["Pan_Cytokeratin"], "Immune": ["CD3e", "CD20", "CD163"],
        "Endocrine": ["INS", "GCG"], "Stromal/vascular": ["CD31", "SMA", "Vimentin"]}


def compartment_counts(pos):
    cmat = np.vstack([np.any([pos[m] for m in ms], axis=0) for ms in COMP.values()]).T
    nC = cmat.sum(1)
    multi_imm = (pos["CD3e"].astype(int) + pos["CD20"].astype(int) + pos["CD163"].astype(int)) >= 2
    beta_alpha = pos["INS"] & pos["GCG"]
    return nC, multi_imm, beta_alpha


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", default="6539")
    ap.add_argument("--thresholds", type=Path, default=BASE / "restore_raw6539_thresholds.csv")
    ap.add_argument("--model", default="SSC")
    ap.add_argument("--comp-mode", type=int, default=0)
    ap.add_argument("--gap-bridge", type=int, default=1)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--alphas", type=float, nargs="*", default=[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5])
    ap.add_argument("--out", type=Path, default=BASE / "redsea_full_validation.png")
    a = ap.parse_args(argv)

    oids, cols, data, edge, sizes, M = R.load_intermediates(a.donor)
    cols = list(cols); ci = {c: k for k, c in enumerate(cols)}
    if a.gap_bridge > 0 and (R.MASK_DIR / f"{a.donor}.tif").exists():
        mask = tifffile.imread(R.MASK_DIR / f"{a.donor}.tif")
        M = R.contact_matrix_gapbridge(mask, len(oids), a.gap_bridge); del mask

    td = pd.read_csv(a.thresholds); td = td[td.model == a.model]
    thr = {r.marker: float(r.threshold) for r in td.itertuples()}
    rawp = pd.read_parquet(sorted(glob.glob(str(BASE / f"cells/donor_id={a.donor}" / "*.parquet")))[0],
                           columns=["object_id"] + MARK)
    rv = rawp.set_index("object_id").reindex(pd.Index(oids))
    posr = {m: (rv[m].to_numpy() >= thr[m]) for m in MARK}
    nCr, mir, bar = compartment_counts(posr)

    def posc_at(alpha):
        cor, _, _ = R.compensate(data, edge, sizes, M, a.comp_mode, alpha)
        return {m: (np.nan_to_num(cor[:, ci[m]]) >= thr[m]) for m in MARK}

    cur = {"a": [], "x3": [], "mi": [], "one": [], "x2": []}
    for al in a.alphas:
        nC, mi, ba = compartment_counts(posc_at(al))
        cur["a"].append(al); cur["x3"].append(100 * (nC >= 3).mean())
        cur["mi"].append(100 * mi.mean()); cur["one"].append(100 * (nC == 1).mean())
        cur["x2"].append(100 * (nC >= 2).mean())

    posc = posc_at(a.alpha)
    nCc, mic, bac = compartment_counts(posc)

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(21, 6))

    # A: genuine impossible artifacts, raw vs corrected
    labels = ["cross-comp\n≥2", "cross-comp\n≥3", "multi-immune\n(≥2 T/B/Mac)", "Beta+Alpha\n(doublet)"]
    rawv = [100 * (nCr >= 2).mean(), 100 * (nCr >= 3).mean(), 100 * mir.mean(), 100 * bar.mean()]
    corv = [100 * (nCc >= 2).mean(), 100 * (nCc >= 3).mean(), 100 * mic.mean(), 100 * bac.mean()]
    x = np.arange(len(labels)); bw = 0.4
    ax[0].bar(x - bw / 2, rawv, bw, label="raw", color="#999999")
    ax[0].bar(x + bw / 2, corv, bw, label=f"REDSEA α={a.alpha}", color="#3366CC")
    for xi, (r, c) in enumerate(zip(rawv, corv)):
        if r > 0.01:
            ax[0].text(xi, max(r, c) + max(rawv) * 0.02, f"−{100*(1-c/r):.0f}%", ha="center", fontsize=8, color="#3366CC")
    ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=9)
    ax[0].set_ylabel("% of cells"); ax[0].legend()
    ax[0].set_title("Genuine impossible artifacts collapse\n(broad mutually-exclusive compartments)")

    # B: compartment-count distribution
    xk = np.arange(5)
    ax[1].bar(xk - bw / 2, [100 * (nCr == k).mean() for k in xk], bw, label="raw", color="#999999")
    ax[1].bar(xk + bw / 2, [100 * (nCc == k).mean() for k in xk], bw, label=f"REDSEA α={a.alpha}", color="#3366CC")
    ax[1].axvspan(1.5, 4.5, color="red", alpha=0.06)
    ax[1].set_xticks(xk); ax[1].set_xlabel("# mutually-exclusive compartments positive")
    ax[1].set_ylabel("% of cells"); ax[1].legend()
    ax[1].set_title(f"Compartment-count: {100*(nCc<=1).mean():.1f}% of cells now clean (0–1);\n"
                    f"single-compartment preserved {100*(nCr==1).mean():.1f}%→{100*(nCc==1).mean():.1f}%")

    # C: genuine-artifact vs alpha (note residual ≥2 is mostly real Vimentin biology)
    ax[2].plot(cur["a"], cur["x3"], "o-", color="#CC3311", label="cross-compartment ≥3 (impossible)")
    ax[2].plot(cur["a"], cur["mi"], "s-", color="#EE7733", label="multi-immune (impossible)")
    ax[2].plot(cur["a"], cur["one"], "^-", color="#009988", label="single-compartment (preserve)")
    ax[2].axvline(a.alpha, color="#3366CC", ls="--", lw=1.5, label=f"operating point α={a.alpha}")
    ax[2].set_xlabel("subtraction strength α"); ax[2].set_ylabel("% of cells")
    ax[2].legend(fontsize=9); ax[2].set_title("Tradeoff: remove impossible artifacts vs preserve true cells")

    fig.suptitle(f"Pixel-level REDSEA validation — donor {a.donor} "
                 f"(subtract-only, 1-px band, gap-bridged; gate fixed at raw RESTORE thresholds)\n"
                 f"cross-compartment≥3 {100*(nCr>=3).mean():.3f}%→{100*(nCc>=3).mean():.3f}%, "
                 f"multi-immune {100*mir.mean():.3f}%→{100*mic.mean():.3f}%  "
                 f"(residual ≥2 is mostly real Vimentin⁺ vascular/stromal biology, not spillover)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150, bbox_inches="tight")
    print(f"[saved] {a.out}")

    print(f"\n=== operating point α={a.alpha} (subtract-only, 1px, gap-bridged) ===")
    print(f"  cross-compartment >=2 : {100*(nCr>=2).mean():.3f}% -> {100*(nCc>=2).mean():.3f}%  (residual ~ Vimentin biology)")
    print(f"  cross-compartment >=3 : {100*(nCr>=3).mean():.3f}% -> {100*(nCc>=3).mean():.3f}%")
    print(f"  multi-immune          : {100*mir.mean():.3f}% -> {100*mic.mean():.3f}%")
    print(f"  Beta+Alpha (doublet)  : {100*bar.mean():.3f}% -> {100*bac.mean():.3f}%  (segmentation, not spillover)")
    print(f"  single-compartment    : {100*(nCr==1).mean():.3f}% -> {100*(nCc==1).mean():.3f}%  (preserve)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
