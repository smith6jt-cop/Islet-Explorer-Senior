#!/usr/bin/env python3
"""
Build the donor acquisition-batch covariate table from the qptiff headers.

Acquisition batch is confounded with disease group in this cohort: no ND donor was acquired on
Fusion 2.4.1. It cannot be corrected away, only disclosed — so it has to exist as data before Gate 5,
and it currently exists nowhere: no schema, no model, and the composition figure runs Kruskal-Wallis
and Mann-Whitney with no batch term.

Reads ONLY each qptiff's first-page description (no pixel data), so it is safe to run against /data.

    python scripts/build_acquisition_batch_table.py \
        --images /path/to/Images \
        --donor-metadata /path/to/donor_metadata_panc.xlsx
    -> data/donor_acquisition_batch.csv
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
from phenocycler.cohort import DONOR_EXCLUSIONS                                # noqa: E402

# The 2026-07 cohort expansion: 15 donors acquired first, then these 7. A separate staining/imaging run
# with a physical panel change, i.e. a real acquisition batch. Cross-checked against panel_variant below.
BATCH2_DONORS = frozenset({"6450", "6505", "6523", "6534", "6566", "6591", "6623"})


def _first_page_description(path: str) -> str:
    import tifffile
    with tifffile.TiffFile(path) as t:
        return t.pages[0].description or ""


def _channel_names(path: str) -> list[str]:
    """Channel biomarker names, one per page, in acquisition order."""
    import tifffile
    names = []
    with tifffile.TiffFile(path) as t:
        for page in t.series[0].pages:
            d = page.description or ""
            m = re.search(r"<Biomarker>(.*?)</Biomarker>", d)
            if m:
                names.append(m.group(1))
    return names


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--images", type=Path, required=True,
                    help="Directory of source .qptiff files (only first-page headers are read).")
    ap.add_argument("--donor-metadata", type=Path, required=True,
                    help="Excel file with donor_id and disease.status columns.")
    ap.add_argument("--out", type=Path, default=REPO / "data" / "donor_acquisition_batch.csv")
    a = ap.parse_args(argv)

    rows = []
    for f in sorted(glob.glob(str(a.images / "*.qptiff"))):
        base = os.path.basename(f)
        m = re.match(r"(?:[A-Za-z]+)?(\d+)", base)          # HDL115pancLN_... -> 115
        if not m:
            # Fail closed: a null donor_id would read as eligible and batch1 below, silently
            # admitting an unknown donor and bypassing the central DONOR_EXCLUSIONS.
            raise SystemExit(f"[err] cannot parse donor id from qptiff filename: {base!r}")
        desc = _first_page_description(f)
        soft = re.search(r"<AcquisitionSoftware>(.*?)</AcquisitionSoftware>", desc)
        chans = _channel_names(f)
        rows.append({
            "donor_id": m.group(1),
            "image_file": base,
            "acquisition_software": soft.group(1) if soft else None,
            "n_channels": len(chans),
            "has_CD38": "CD38" in chans,
            "has_b_Catenin1": any(c.lower().startswith("b-catenin") or c.lower().startswith("b_catenin")
                                  for c in chans),
            "channels": ";".join(chans),
        })
    t = pd.DataFrame(rows)

    meta = pd.read_excel(a.donor_metadata)[["donor_id", "disease.status"]]
    meta["donor_id"] = meta["donor_id"].astype(str)
    t = t.merge(meta, on="donor_id", how="left")
    t["eligible"] = ~t["donor_id"].isin(DONOR_EXCLUSIONS)
    # Panel variant: the batch-2 CD38 <-> b-Catenin1 swap, recorded rather than inferred downstream.
    t["panel_variant"] = t.apply(
        lambda r: "CD38" if r.has_CD38 and not r.has_b_Catenin1
        else ("b_Catenin1" if r.has_b_Catenin1 and not r.has_CD38
              else ("both" if r.has_CD38 and r.has_b_Catenin1 else "neither")), axis=1)
    # ACQUISITION BATCH — the covariate to analyse against. The two batches are separate staining and
    # imaging runs months apart with a physical panel change (batch 1 carries CD38, batch 2 carries
    # b-Catenin1 or neither), which is what makes them a batch rather than a label. `panel_variant`
    # and `acquisition_software` are each perfect proxies for it on this cohort — every batch-2 donor
    # is Fusion 2.4.1 and no batch-1 donor is — but the batch is the thing, not the software version.
    t["batch"] = t["donor_id"].map(lambda d: "batch2" if d in BATCH2_DONORS else "batch1")
    for probe, name in ((t["panel_variant"] == "CD38", "panel_variant"),):
        drift = t.loc[t.eligible & (probe != (t["batch"] == "batch1"))]
        if len(drift):                        # fail loudly rather than let the batch definition rot
            raise SystemExit(
                f"[err] {name} no longer agrees with BATCH2_DONORS for {drift['donor_id'].tolist()}; "
                "the batch definition is hardcoded and must be re-derived if the cohort changes")

    t = t.sort_values("donor_id").reset_index(drop=True)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    t.drop(columns=["channels"]).to_csv(a.out, index=False)
    t.to_csv(a.out.with_name(a.out.stem + "_with_channels.csv"), index=False)

    e = t[t.eligible]
    ct = pd.crosstab(e["acquisition_software"], e["disease.status"])
    print(ct.to_string(), "\n")
    print(pd.crosstab(e["acquisition_software"], e["panel_variant"]).to_string(), "\n")

    # The specific confound: ND vs the rest, on Fusion 2.4.1 vs the rest.
    from scipy.stats import fisher_exact
    is241 = e["acquisition_software"].eq("Fusion 2.4.1")
    isnd = e["disease.status"].eq("ND")
    table = [[int((is241 & isnd).sum()), int((is241 & ~isnd).sum())],
             [int((~is241 & isnd).sum()), int((~is241 & ~isnd).sum())]]
    odds, p = fisher_exact(table)
    print(f"ND x Fusion-2.4.1 contingency (eligible donors only): {table}")
    print(f"Fisher exact p = {p:.4f}  -> batch is {'CONFOUNDED with' if p < 0.05 else 'independent of'} "
          "disease group; disclose it, do not correct it away")
    print(f"\n[done] {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
