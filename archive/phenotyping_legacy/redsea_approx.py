#!/usr/bin/env python3
"""
Table-level REDSEA approximation (lateral spillover correction) + validation.

REDSEA (Bai et al. 2021) corrects signal that bleeds across shared cell boundaries
from neighbours. The pixel-accurate version (labsyspharm/redseapy) needs segmentation
masks + the multichannel image. This is a cheap *table-level approximation* that uses
only the per-cell centroids we already have, to TEST whether spillover is the lever
behind the "pan-positive bright" cells before committing to full redseapy.

Model, per image (donor): build cell adjacency from centroids via Delaunay; weight each
edge by 1/distance (proxy for shared-boundary length); row-normalise -> W. Spillover-
corrected signal per marker:  corrected = max(0, raw - alpha * (W @ raw)).
(W @ raw = distance-weighted neighbour mean; alpha = total spillover fraction.)

VALIDATION (the point): mutually-exclusive marker pairs (e.g. Pan_Cytokeratin vs CD3e)
should essentially never be double-positive in clean data; lateral spillover inflates
the double-positive fraction. If correction collapses these co-positives while preserving
true single-positives, spillover is confirmed -> escalate to full redseapy.

Run in `scvi-env`. Read-only on data/cells/ and data/restore/; writes only a figure.
"""
from __future__ import annotations
import glob
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix, diags
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2] / "data"

# the 9 RESTORE markers (per-image thresholds available) = the mutually-exclusive set
MARKERS = ["CD3e", "CD20", "CD163", "Pan_Cytokeratin", "Vimentin", "INS", "GCG", "CD31", "SMA"]
# impossible co-occurrences (lineage A vs lineage B); should be ~0% in clean data
PAIRS = [
    ("Pan_Cytokeratin", "CD3e"), ("Pan_Cytokeratin", "CD163"), ("Pan_Cytokeratin", "Vimentin"),
    ("Pan_Cytokeratin", "CD31"), ("Pan_Cytokeratin", "INS"), ("Pan_Cytokeratin", "CD20"),
    ("CD3e", "CD163"), ("CD3e", "INS"), ("CD3e", "SMA"), ("CD31", "Vimentin"),
]
ALPHAS = [0.1, 0.25, 0.5]


def row_norm_adjacency(xy: np.ndarray):
    """Delaunay neighbour graph, edge weight 1/dist, row-normalised (sparse)."""
    tri = Delaunay(xy)
    s = tri.simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    i, j = e[:, 0], e[:, 1]
    d = np.linalg.norm(xy[i] - xy[j], axis=1)
    w = 1.0 / np.maximum(d, 1.0)
    n = len(xy)
    I = np.concatenate([i, j]); J = np.concatenate([j, i]); W = np.concatenate([w, w])
    A = csr_matrix((W, (I, J)), shape=(n, n))
    rs = np.asarray(A.sum(1)).ravel(); rs[rs == 0] = 1.0
    return diags(1.0 / rs) @ A


def main():
    # per-image applied RESTORE thresholds (image x marker -> threshold)
    pf = pd.read_csv(BASE / "restore" / "positive_fractions.csv")
    pf["image"] = pf["image"].astype(str)
    THR = {(r.image, r.marker): float(r.threshold) for r in pf.itertuples()}

    pair_rows, single_rows = [], []
    files = sorted(glob.glob(str(BASE / "cells" / "donor_id=*" / "*.parquet")))
    for f in files:
        donor = f.split("donor_id=")[1].split("/")[0]
        df = pd.read_parquet(f, columns=["X_centroid", "Y_centroid"] + MARKERS)
        xy = df[["X_centroid", "Y_centroid"]].to_numpy(np.float64)
        raw = df[MARKERS].to_numpy(np.float64)
        thr = np.array([THR.get((donor, m), np.inf) for m in MARKERS])
        Wn = row_norm_adjacency(xy)
        nbr = Wn @ raw                              # weighted neighbour mean per marker
        pos_raw = raw >= thr
        n = len(raw)
        print(f"[{donor}] n={n:,}", flush=True)
        for a in ALPHAS:
            # REDSEA subtract + REINFORCE (table form): sharpen cell vs neighbours.
            # spillover recipients (dim-among-bright) drop; true positives preserved/boosted;
            # homogeneous same-type regions (raw ~= nbr) unchanged -> no over-erosion.
            pos_c = np.maximum(0.0, (1.0 + a) * raw - a * nbr) >= thr
            for m1, m2 in PAIRS:
                k1, k2 = MARKERS.index(m1), MARKERS.index(m2)
                pair_rows.append(dict(donor=donor, pair=f"{m1}×{m2}", alpha=a, n=n,
                                      dp_raw=int((pos_raw[:, k1] & pos_raw[:, k2]).sum()),
                                      dp_corr=int((pos_c[:, k1] & pos_c[:, k2]).sum())))
            for k, m in enumerate(MARKERS):
                single_rows.append(dict(donor=donor, marker=m, alpha=a, n=n,
                                        sp_raw=int(pos_raw[:, k].sum()),
                                        sp_corr=int(pos_c[:, k].sum())))

    P = pd.DataFrame(pair_rows); S = pd.DataFrame(single_rows)
    # aggregate across all donors -> cohort fractions
    pa = P.groupby(["pair", "alpha"]).agg(n=("n", "sum"), dp_raw=("dp_raw", "sum"),
                                          dp_corr=("dp_corr", "sum")).reset_index()
    pa["pct_raw"] = 100 * pa.dp_raw / pa.n
    pa["pct_corr"] = 100 * pa.dp_corr / pa.n
    sa = S.groupby(["marker", "alpha"]).agg(n=("n", "sum"), sp_raw=("sp_raw", "sum"),
                                            sp_corr=("sp_corr", "sum")).reset_index()
    sa["pct_raw"] = 100 * sa.sp_raw / sa.n
    sa["pct_corr"] = 100 * sa.sp_corr / sa.n

    # ---- figure: (1) impossible double-positives before/after; (2) single-pos preservation ----
    pairs = list(dict.fromkeys(pa["pair"]))
    fig, ax = plt.subplots(1, 2, figsize=(17, 6))
    x = np.arange(len(pairs)); bw = 0.2
    raw_pct = [pa[(pa.pair == p) & (pa.alpha == ALPHAS[0])].pct_raw.iloc[0] for p in pairs]
    ax[0].bar(x - 1.5 * bw, raw_pct, bw, label="raw", color="#444")
    for n_a, a in enumerate(ALPHAS):
        v = [pa[(pa.pair == p) & (pa.alpha == a)].pct_corr.iloc[0] for p in pairs]
        ax[0].bar(x + (n_a - 0.5) * bw, v, bw, label=f"corr α={a}")
    ax[0].set_xticks(x); ax[0].set_xticklabels(pairs, rotation=45, ha="right", fontsize=8)
    ax[0].set_ylabel("% cells double-positive (cohort)")
    ax[0].set_title("Impossible (mutually-exclusive) co-positives — should collapse")
    ax[0].legend(fontsize=8)
    mk = MARKERS
    xm = np.arange(len(mk))
    rawm = [sa[(sa.marker == m) & (sa.alpha == ALPHAS[0])].pct_raw.iloc[0] for m in mk]
    ax[1].bar(xm - 1.5 * bw, rawm, bw, label="raw", color="#444")
    for n_a, a in enumerate(ALPHAS):
        v = [sa[(sa.marker == m) & (sa.alpha == a)].pct_corr.iloc[0] for m in mk]
        ax[1].bar(xm + (n_a - 0.5) * bw, v, bw, label=f"corr α={a}")
    ax[1].set_xticks(xm); ax[1].set_xticklabels(mk, rotation=45, ha="right", fontsize=8)
    ax[1].set_ylabel("% cells positive (cohort)")
    ax[1].set_title("Single-positive fraction — should be preserved (not over-corrected)")
    ax[1].legend(fontsize=8)
    fig.suptitle("Table-level REDSEA (subtract + reinforce) — spillover correction validation (15 donors)", fontsize=13)
    fig.tight_layout()
    out = BASE / "redsea_approx_validation.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("[saved]", out)

    # concise numeric summary to stdout (mean over pairs)
    print("\n=== mean % double-positive over the 10 impossible pairs ===")
    print("  raw          :", round(pa[pa.alpha == ALPHAS[0]].pct_raw.mean(), 3))
    for a in ALPHAS:
        print(f"  corrected α={a}:", round(pa[pa.alpha == a].pct_corr.mean(), 3))
    print("\n=== mean single-positive % over 9 markers (preservation) ===")
    print("  raw          :", round(sa[sa.alpha == ALPHAS[0]].pct_raw.mean(), 3))
    for a in ALPHAS:
        print(f"  corrected α={a}:", round(sa[sa.alpha == a].pct_corr.mean(), 3))


if __name__ == "__main__":
    main()
