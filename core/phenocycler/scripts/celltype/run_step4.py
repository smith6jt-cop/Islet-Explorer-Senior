"""Step 4: UTAG tissue domains, per donor, beside the per-cell calls.

    python scripts/celltype/run_step4.py --donors 6450        # one donor
    python scripts/celltype/run_step4.py                      # whole cohort

Domains are discovered per donor.  Pooling would let one donor's tissue define another's
neighbourhoods, and this cohort's donors differ in exactly the way that would matter --
peri-islet content spans 0.03% to 8.35%.  Comparing domains across donors is then a matter of
what each one contains, not of a shared cluster id.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402
import utag_adapter as UA                        # noqa: E402

#: Cells per donor handed to UTAG.  The neighbour graph and the Leiden step, not the message
#: passing, set the cost, and it does not scale linearly: 40k cells took 100 s for two
#: resolutions, while 120k took over 12 minutes for two of three.  A cap also keeps every
#: donor comparable in graph resolution rather than letting the largest slide get a finer one,
#: and the realised count is recorded per run.
MAX_CELLS = 100_000
RESOLUTIONS = (0.1, 0.3, 1.0)


def run_donor(donor: str, markers: list, max_cells: int, resolutions, seed: int = 0):
    utag = UA._import_utag()
    a = UA.donor_adata(donor, markers, max_cells=max_cells, seed=seed)
    t = time.time()
    res = utag(a, slide_key=None, max_dist=UA.MAX_DIST_UM, normalization_mode="l1_norm",
               apply_clustering=True, clustering_method=["leiden"],
               resolutions=list(resolutions), parallel=False, return_copy=True)
    keys = [c for c in res.obs.columns if c.startswith("UTAG")]
    out = res.obs[["object_id", "cell_region", "islet_num", "donor_id"] + keys].copy()
    out = out.rename(columns={k: k.replace("UTAG Label_leiden_", "domain_r") for k in keys})
    return out, time.time() - t


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--max-cells", type=int, default=MAX_CELLS)
    ap.add_argument("--resolutions", nargs="*", type=float, default=list(RESOLUTIONS))
    args = ap.parse_args()

    donors = args.donors or panel.cohort()["donor_id"].tolist()
    markers = [m for m in panel.markers() if m != "DAPI"]
    print(f"Step 4 — UTAG @ {UA.UTAG_COMMIT[:12]}, {len(donors)} donors, {len(markers)} markers, "
          f"max_dist {UA.MAX_DIST_UM} um, up to {args.max_cells:,} cells/donor")

    UA.OUT.mkdir(parents=True, exist_ok=True)
    conc = []
    for d in donors:
        out, secs = run_donor(d, markers, args.max_cells, args.resolutions)
        p = UA.OUT / f"donor_id={d}"
        p.mkdir(parents=True, exist_ok=True)
        out.to_parquet(p / "data.parquet", index=False)
        cols = [c for c in out.columns if c.startswith("domain_r")]
        print(f"  {d}: {len(out):,} cells, "
              + ", ".join(f"{c} -> {out[c].nunique()} domains" for c in cols)
              + f"   {secs:.0f}s", flush=True)
        for c in cols:
            t = UA.domain_islet_concordance(out, c)
            t.insert(0, "resolution", c.replace("domain_r", ""))
            t.insert(0, "donor_id", d)
            conc.append(t)
    C = pd.concat(conc, ignore_index=True)
    C.to_parquet(panel.OUT / "domain_islet_concordance.parquet", index=False)

    print("\nDoes any domain recover the islets QuPath drew?  `cell_region` was never an input.")
    print("Best domain per donor and resolution, by islet-core enrichment over that donor's base rate:\n")
    best = C.sort_values("enrichment", ascending=False).groupby(
        ["donor_id", "resolution"], as_index=False).first()
    print(best[["donor_id", "resolution", "domain", "n", "core_frac", "base_rate",
                "enrichment"]].to_string(index=False, float_format=lambda v: f"{v:8.2f}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
