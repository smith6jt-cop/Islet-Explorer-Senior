#!/usr/bin/env python3
"""Per-islet composition, in a form the app writers can each format from.

Deliberately one intermediate rather than several near-duplicates: the app wants the same counts
as proportions in a CSV, as hormone co-expression in one workbook sheet, and as structure density
in another, and deriving all three from one table is what keeps them consistent.

Grain: one row per (islet, region, class).  Regions are the four from ``build_islet_geometry``;
classes are the display labels.  ``n_cells`` is a count and ``prop`` is that count over the
region's total, so a reader can always recover the denominator -- which matters here because the
denominator changed when cells without a nucleus were excluded.

    python core/phenocycler/scripts/typing/aggregate_islets.py [--donors ...]
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

REGIONS = ("core", "expansion", "non_endocrine", "union")
HORMONE = {"Beta_cell": "Ins", "Alpha_cell": "Glu", "Delta_cell": "Stt"}


def _region_membership(cells: pd.DataFrame) -> dict[str, pd.Series]:
    core = cells["cell_region_20um"] == "core"
    expansion = cells["cell_region_20um"] == "expansion"
    union = core | expansion
    return {"core": core, "expansion": expansion, "union": union,
            "non_endocrine": union & ~cells["is_endocrine"]}


def build_donor(donor: str) -> pd.DataFrame:
    cells = pq.read_table(
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        columns=["islet_key", "combined_islet_id", "cell_region_20um", "phenotype",
                 "is_endocrine", "object_context", "cell_class"]).to_pandas()
    # The signed distance and QuPath's own parent column disagree at the islet boundary for a
    # small number of cells: inside by distance, assigned to no islet by QuPath.  They cannot be
    # aggregated per islet, so they are dropped -- counted here rather than vanishing.
    orphan = cells["islet_key"].isna() & cells["cell_region_20um"].isin(["core", "expansion"])
    cells = cells[cells["islet_key"].notna()]

    # `Excluded` is the exocrine-inside-an-islet-border set the maintainer ruled out of analysis.
    # Out of the DENOMINATOR, not carried as a class: an islet is an endocrine structure, so
    # "exocrine inside one" is a contradiction to be excluded, and leaving it in the composition
    # both invents a phenotype and understates every real one.  Counted, never silently dropped.
    excluded = int((cells["cell_class"] == "Excluded").sum())
    cells = cells[cells["cell_class"] != "Excluded"]
    membership = _region_membership(cells)

    frames = []
    for region, mask in membership.items():
        subset = cells[mask]
        if subset.empty:
            continue
        counts = (subset.groupby(["combined_islet_id", "islet_key", "phenotype"])
                  .size().rename("n_cells").reset_index())
        totals = counts.groupby("combined_islet_id")["n_cells"].transform("sum")
        counts["region"] = region
        counts["region_total"] = totals
        counts["prop"] = counts["n_cells"] / totals
        frames.append(counts)

    out = pd.concat(frames, ignore_index=True)
    out["donor_id"] = donor
    out.attrs["orphan_cells"] = int(orphan.sum())
    out.attrs["excluded_cells"] = excluded
    return out


def hormone_composition(donor: str) -> pd.DataFrame:
    """Per islet core, the hormone breakdown the workbook's composition sheet carries.

    ``Triple_Neg`` means a hormone-negative ENDOCRINE cell.  This method produces none by
    construction -- a cell in an islet that clears no hormone is not called endocrine at all, it
    continues down the ordinary chain -- so the column is present and zero rather than absent, and
    the non-endocrine intruders are reported as themselves instead of as failed beta cells.
    """
    cells = pq.read_table(
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        columns=["combined_islet_id", "islet_key", "cell_region_20um", "phenotype",
                 "is_endocrine", "cell_class"]).to_pandas()
    core = cells[(cells["cell_region_20um"] == "core") & cells["islet_key"].notna()
                 & (cells["cell_class"] != "Excluded")]
    rows = []
    for (cid, key), group in core.groupby(["combined_islet_id", "islet_key"]):
        row = {"donor_id": donor, "combined_islet_id": cid, "islet_key": key,
               "cells_total": int(len(group)),
               "n_endocrine": int(group["is_endocrine"].sum()),
               "Multi_Pos": 0, "Triple_Neg": 0}
        for phenotype, short in HORMONE.items():
            n = int((group["phenotype"] == phenotype).sum())
            row[f"{short}_single"] = n
            row[f"{short}_any"] = n
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    composition, hormones = [], []
    for donor in donors:
        composition.append(build_donor(donor))
        hormones.append(hormone_composition(donor))
        n = composition[-1].combined_islet_id.nunique()
        orphan = composition[-1].attrs["orphan_cells"]
        excluded = composition[-1].attrs["excluded_cells"]
        print(f"  {donor}: {n:>5,} islets   "
              f"({excluded:,} exocrine-inside-border excluded, "
              f"{orphan:,} boundary cells in no islet)", flush=True)

    target = S.OUT / "app"
    target.mkdir(parents=True, exist_ok=True)
    comp = pd.concat(composition, ignore_index=True)
    horm = pd.concat(hormones, ignore_index=True)
    comp.to_parquet(target / "islet_composition_long.parquet", index=False)
    horm.to_parquet(target / "islet_hormone_composition.parquet", index=False)

    core = comp[comp.region == "core"]
    print(f"\n{comp.combined_islet_id.nunique():,} islets, {len(donors)} donors -> {target}")
    # Pivot with fill_value before averaging: in the long form an islet lacking a phenotype has no
    # row for it, so a plain groupby().mean() averages each class over a different set of islets
    # and the "proportions" sum to well over 1.
    wide = core.pivot_table(index="combined_islet_id", columns="phenotype",
                            values="prop", fill_value=0.0)
    print("\ncohort mean per-islet proportion, islet core (absent counts as zero):")
    print(wide.mean().sort_values(ascending=False).to_string())
    print(f"  sum = {wide.mean().sum():.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
