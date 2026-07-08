#!/usr/bin/env python3
"""
Cohort validation of pixel-level REDSEA across all 15 donors (one figure).

Two questions:
 (1) Does REDSEA remove the genuine impossible artifacts cohort-wide (not just on 6539)?
 (2) Does the corrected data STILL recover the T1D biology (beta-cell loss + T-cell
     infiltration across ND -> Aab+ -> T1D)? — the ultimate sanity check that the
     correction sharpened signal rather than destroying it.

Inputs:
  data/restore/positive_fractions.csv            raw RESTORE per-image fractions+thresholds (robust)
  data/restore_redsea/positive_fractions.csv     re-RESTORE on corrected cells (run first)
  data/cells/ , data/cells_redsea/               per-cell raw / corrected MFI
  ~/IO60panc2nd/donor_metadata_panc.xlsx         disease.status per donor

  # 1. re-RESTORE the corrected cohort (robust, all 15):
  python scripts/senior/restore_normalize.py --cells-dir data/cells_redsea \
      --out-dir data/restore_redsea --thresh-csv data/restore_redsea_thresholds.csv \
      --gated-dir data/restore_gated_redsea
  # 2. this figure:
  python scripts/senior/redsea_cohort_validate.py
"""
from __future__ import annotations
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2] / "data"
META = Path("/home/smith6jt/IO60panc2nd/donor_metadata_panc.xlsx")
STATUS_ORDER = ["ND", "AAB", "T1D"]
STATUS_COL = {"ND": "#4477AA", "AAB": "#CC6633", "T1D": "#228833"}
MARK = ["Pan_Cytokeratin", "CD3e", "CD20", "CD163", "INS", "GCG", "CD31", "SMA", "Vimentin"]


def status_map():
    m = pd.read_excel(META)
    return dict(zip(m["donor_id"].astype(str), m["disease.status"].astype(str)))


def frac_table(pf_csv):
    pf = pd.read_csv(pf_csv); pf["image"] = pf["image"].astype(str)
    return pf


def artifacts_per_donor(thr_lut):
    """Gate raw & corrected cells (fixed raw thresholds) -> cross-comp>=3 + multi-immune per donor."""
    rows = []
    for f in sorted(glob.glob(str(BASE / "cells" / "donor_id=*" / "*.parquet"))):
        d = f.split("donor_id=")[1].split("/")[0]
        cf = BASE / "cells_redsea" / f"donor_id={d}" / "data_0.parquet"
        if not cf.exists():
            continue
        raw = pd.read_parquet(f, columns=["object_id"] + MARK)
        cor = pd.read_parquet(cf, columns=["object_id"] + MARK).set_index("object_id").reindex(raw["object_id"])
        thr = {m: thr_lut.get((d, m), np.inf) for m in MARK}
        def metrics(df):
            P = {m: (df[m].to_numpy() >= thr[m]) for m in MARK}
            nC = (P["Pan_Cytokeratin"].astype(int) + (P["CD3e"] | P["CD20"] | P["CD163"]).astype(int)
                  + (P["INS"] | P["GCG"]).astype(int) + (P["CD31"] | P["SMA"] | P["Vimentin"]).astype(int))
            mi = (P["CD3e"].astype(int) + P["CD20"].astype(int) + P["CD163"].astype(int)) >= 2
            return 100 * (nC >= 3).mean(), 100 * mi.mean()
        r3, rmi = metrics(raw)
        c3, cmi = metrics(cor.reset_index())
        rows.append(dict(donor=d, raw_x3=r3, cor_x3=c3, raw_mi=rmi, cor_mi=cmi))
    return pd.DataFrame(rows)


def main():
    smap = status_map()
    raw_pf = frac_table(BASE / "restore" / "positive_fractions.csv")
    cor_pf = frac_table(BASE / "restore_redsea" / "positive_fractions.csv")
    # robust applied thresholds for the fixed-gate artifact metric (from raw RESTORE)
    thr_lut = {(r.image, r.marker): float(r.threshold) for r in raw_pf.itertuples()}

    def by_status(pf, marker):
        sub = pf[pf.marker == marker].copy()
        sub["status"] = sub.image.map(smap)
        return sub
    art = artifacts_per_donor(thr_lut)
    art["status"] = art.donor.map(smap)

    fig, ax = plt.subplots(1, 3, figsize=(21, 6))
    # Panel A & B: biology — INS (beta) and CD3e (T) fraction by status, raw vs corrected
    for axi, (marker, title, expect) in [(0, ("INS", "β-cell (INS⁺) fraction", "loss ND→T1D")),
                                          (1, ("CD3e", "T-cell (CD3e⁺) fraction", "gain ND→T1D"))]:
        rb = by_status(raw_pf, marker); cb = by_status(cor_pf, marker)
        xs = np.arange(len(STATUS_ORDER)); bw = 0.35
        rm = [100 * rb[rb.status == s].frac_pos.mean() for s in STATUS_ORDER]
        cm = [100 * cb[cb.status == s].frac_pos.mean() for s in STATUS_ORDER]
        ax[axi].bar(xs - bw / 2, rm, bw, label="raw RESTORE", color="#bbbbbb")
        ax[axi].bar(xs + bw / 2, cm, bw, label="REDSEA+RESTORE", color="#3366CC")
        # per-donor dots
        for s_i, s in enumerate(STATUS_ORDER):
            cv = cb[cb.status == s].frac_pos * 100
            ax[axi].scatter(np.full(len(cv), s_i + bw / 2), cv, color=STATUS_COL[s], s=18, zorder=3, edgecolor="k", lw=0.3)
        ax[axi].set_xticks(xs); ax[axi].set_xticklabels(STATUS_ORDER)
        ax[axi].set_ylabel("% cells positive"); ax[axi].set_title(f"{title}\n(expect {expect})"); ax[axi].legend()

    # Panel C: genuine artifacts removed cohort-wide
    art2 = art.sort_values("raw_x3", ascending=False)
    xs = np.arange(len(art2)); bw = 0.4
    ax[2].bar(xs - bw / 2, art2.raw_x3, bw, label="raw cross-comp≥3", color="#cc9999")
    ax[2].bar(xs + bw / 2, art2.cor_x3, bw, label="REDSEA cross-comp≥3", color="#3366CC")
    ax[2].set_xticks(xs); ax[2].set_xticklabels(art2.donor, rotation=90, fontsize=8)
    ax[2].set_ylabel("% cells (impossible cross-comp≥3)")
    ax[2].set_title(f"Genuine artifacts removed (all donors)\n"
                    f"cohort mean {art.raw_x3.mean():.3f}%→{art.cor_x3.mean():.3f}%, "
                    f"multi-immune {art.raw_mi.mean():.3f}%→{art.cor_mi.mean():.3f}%")
    ax[2].legend()
    fig.suptitle("REDSEA cohort validation (15 donors): biology preserved + artifacts removed", fontsize=14)
    fig.tight_layout(); fig.savefig(BASE / "redsea_cohort_validation.png", dpi=150, bbox_inches="tight")
    print(f"[saved] {BASE/'redsea_cohort_validation.png'}")

    # numeric summary
    for marker, name in [("INS", "β (INS)"), ("CD3e", "T (CD3e)")]:
        rb = by_status(raw_pf, marker); cb = by_status(cor_pf, marker)
        print(f"\n{name} % positive by status (raw / corrected):")
        for s in STATUS_ORDER:
            print(f"  {s:4s} raw {100*rb[rb.status==s].frac_pos.mean():.2f}  corrected {100*cb[cb.status==s].frac_pos.mean():.2f}")
    print(f"\ncross-comp≥3 cohort: raw {art.raw_x3.mean():.3f}% -> corrected {art.cor_x3.mean():.3f}%")
    print(f"multi-immune cohort: raw {art.raw_mi.mean():.3f}% -> corrected {art.cor_mi.mean():.3f}%")
    print(art[["donor", "status", "raw_x3", "cor_x3", "raw_mi", "cor_mi"]].to_string(index=False))


if __name__ == "__main__":
    raise SystemExit(main())
