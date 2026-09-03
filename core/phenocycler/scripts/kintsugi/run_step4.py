#!/usr/bin/env python
"""Step 4 driver: measure this lane against the cosine lane on the cells both call.

    python scripts/kintsugi/run_step4.py --calls 6551_armA_n25
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
sys.path.insert(0, str(HERE))

import concordance                                  # noqa: E402
import panel                                        # noqa: E402

OUT = panel.REPO / "data" / "development" / "kintsugi" / "concordance"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="stem under calls/, e.g. 6551_armA_n25")
    args = ap.parse_args()

    mine = pq.read_table(
        panel.REPO / "data" / "development" / "kintsugi" / "calls"
        / f"{args.calls}_cells.parquet").to_pandas()[["donor_id", "object_id", "call"]]
    donors = sorted(mine["donor_id"].unique())

    theirs = pq.read_table(concordance.COSINE_ISLET).to_pandas()
    theirs = theirs[theirs["donor_id"].isin(donors)][["donor_id", "object_id", "call"]]

    matrix, summary = concordance.compare(mine, theirs)
    summary["refusal_overlap"] = concordance.refusal_overlap(mine, theirs)
    summary["donors"] = donors
    summary["population"] = ("islet-associated cells only (core + peri), from "
                             "celltype/islet_residents_all.parquet")

    OUT.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(OUT / f"{args.calls}_matrix.csv")
    (OUT / f"{args.calls}_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    print(f"  joined {summary['n_joined']:,} islet-associated cells")
    print(f"  both called          {summary['n_both_called']:,}")
    print(f"  label agreement      {summary['label_agreement']:.4f}")
    print(f"  branch agreement     {summary['branch_agreement']:.4f}")
    ro = summary["refusal_overlap"]
    print(f"  refused by kintsugi  {ro['kintsugi_refused']:,}")
    print(f"  refused by cosine    {ro['cosine_refused']:,}")
    print(f"  refused by both      {ro['refused_by_both']:,}  (Jaccard {ro['jaccard']:.4f})")
    print(f"\nwrote {OUT}/{args.calls}_*.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
