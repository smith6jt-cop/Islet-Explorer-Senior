"""Prove the normalisation did exactly what it claims, on the written parquets.

Two claims carry the whole design and both are checkable rather than arguable:

1. the shift is a **pure per-(donor, marker) multiplicative rescaling** -- if it were
   anything else, distribution shape and within-donor rank would not be preserved and the
   "standardize the unit, never the distribution" contract would be broken;
2. therefore **every rank statistic from Step 1 is invariant under it** -- which is what lets
   Step 1's tail numbers check Step 2 without circularity.

    python scripts/celltype/verify_step2.py        # from core/phenocycler/
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                             # noqa: E402
import usability as U                    # noqa: E402
import run_step2 as R2                   # noqa: E402

CHECK_PAIRS = [("CD20", "CD79a"), ("CD3e", "CD8"), ("CD3e", "FOXP3"), ("MPO", "CD66"),
               ("Ker8_18", "Pan_Cytokeratin"), ("CD31", "CD34"), ("B3TUBB", "CD56"),
               ("CD20", "SMA"), ("CD20", "INS"), ("FOXP3", "SMA"), ("TOX", "E_cadherin")]


def main() -> int:
    meta = json.loads((panel.OUT / "shifts.json").read_text())
    shifts, donors, markers = meta["shifts"], meta["donors"], meta["markers"]
    fails = []

    print("1. the shift is a pure per-(donor, marker) multiplicative rescaling")
    for d in donors:
        raw = pq.read_table(panel.EXPRESSION / f"donor_id={d}" / "data.parquet",
                            columns=["object_id"] + markers).to_pandas()
        nrm = pq.read_table(R2.OUT / f"donor_id={d}" / "data.parquet",
                            columns=["object_id"] + markers).to_pandas()
        if len(raw) != len(nrm) or not (raw["object_id"] == nrm["object_id"]).all():
            fails.append(f"{d}: object_id order or count changed")
            continue
        worst = 0.0
        for m in markers:
            a, b = raw[m].to_numpy("float64"), nrm[m].to_numpy("float64")
            ok = np.isfinite(a) & np.isfinite(b) & (a > 0)
            if not ok.any():
                continue
            r = b[ok] / a[ok]
            worst = max(worst, float(np.nanmax(np.abs(r - shifts[m]["scale"][d]))))
            if not np.allclose(r, shifts[m]["scale"][d], rtol=1e-9, atol=0):
                fails.append(f"{d}/{m}: ratio is not the recorded constant")
        print(f"   {d}: {len(raw):,} rows, object_id identical, "
              f"max |ratio - recorded scale| = {worst:.2e}")

    print("\n2. every rank statistic is therefore invariant "
          f"(tail overlap at top {U.BARS['tail_q']:.1%}, donor 115 + 6534)")
    print(f"   {'pair':<28}{'raw':>10}{'normalised':>13}")
    for d in ("115", "6534"):
        raw = pq.read_table(panel.EXPRESSION / f"donor_id={d}" / "data.parquet",
                            columns=["qc_analysis_eligible"] + markers).to_pandas()
        nrm = pq.read_table(R2.OUT / f"donor_id={d}" / "data.parquet",
                            columns=["qc_analysis_eligible"] + markers).to_pandas()
        raw, nrm = raw[raw.qc_analysis_eligible], nrm[nrm.qc_analysis_eligible]
        for a, b in CHECK_PAIRS:
            x = U.tail_overlap(raw[a].to_numpy(), raw[b].to_numpy(), U.BARS["tail_q"])
            y = U.tail_overlap(nrm[a].to_numpy(), nrm[b].to_numpy(), U.BARS["tail_q"])
            same = (np.isnan(x) and np.isnan(y)) or x == y
            print(f"   {d} {a} & {b:<18}{x:10.1f}{y:13.1f}   {'' if same else '<-- CHANGED'}")
            if not same:
                fails.append(f"{d}/{a}&{b}: tail overlap changed under a monotone shift")

    print("\n3. what the run recorded")
    prov = pd.Series([v["provenance"] for v in shifts.values()]).value_counts()
    print(prov.to_string())
    print(f"   background-mode SD, median {np.median([v['mode_sd_after'] for v in shifts.values()]):.4f}")

    print("\n" + ("FAILED: " + "; ".join(fails[:5]) if fails else "all checks passed"))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
