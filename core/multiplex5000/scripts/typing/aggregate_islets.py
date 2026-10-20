#!/usr/bin/env python3
"""Per-islet composition, in a form the app writers can each format from.

Deliberately one intermediate rather than several near-duplicates: the app wants the same counts
as proportions in a CSV, as hormone co-expression in one workbook sheet, and as structure density
in another, and deriving all three from one table is what keeps them consistent.

Grain: one row per (islet, region, class).  Regions are the four from ``build_islet_geometry``;
classes are the display labels.  ``n_cells`` is a count and ``prop`` is that count over the
region's total, so a reader can always recover the denominator -- which matters here because the
denominator changed when cells without a nucleus were excluded.

    python core/multiplex5000/scripts/typing/aggregate_islets.py [--donors ...]
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import build_cells as CELLS
import cohort as COHORT
import subpanel_io as S

REGIONS = ("core", "expansion", "non_endocrine", "union")
HORMONE = {"Beta_cell": "Ins", "Alpha_cell": "Glu", "Delta_cell": "Stt"}

#: Every class that can reach ``phenotype``, in the app's spelling, read out of the ONE crosswalk
#: ``build_cells.DISPLAY`` rather than restated here -- so ``Unresolved (no lineage evidence)``,
#: which ``refine.py`` would produce and which no cell on disk carries yet, is already in this
#: list and does not have to be remembered on the day the pass is wired.  A class with no cells
#: is a ZERO, and a zero has to be printed: an omitted line and a line that was never asked about
#: read the same, and the reader cannot tell which one they are looking at.
#:
#: ``Excluded`` is deliberately NOT here.  ``build_donor`` removes those cells from the
#: DENOMINATOR by name below, so they are not zero, they are gone; printing ``Excluded 0.0``
#: would state a proportion about cells that were never counted.
COUNTED_CLASSES = tuple(sorted(set(CELLS.DISPLAY.values()) - {CELLS.DISPLAY["Excluded"]}))

#: What a cell whose ``phenotype`` is NULL is called when it has to become a column.  It is not a
#: class and it is not a zero -- it is one cell whose class is UNAVAILABLE, kept visible under a
#: name that says so.  Folding it into a named class would invent a call; dropping it would take
#: the cell out of a denominator it is already counted in.
UNMAPPED = "(class unavailable)"


def _full_class_columns(observed) -> list[str]:
    """The declared classes first, then anything observed that the crosswalk does not name.

    Nothing observed is ever dropped -- a phenotype this file has never heard of (including the
    null one a broken join would leave) keeps its column and shows up, rather than being tidied
    away by a reindex.
    """
    return list(COUNTED_CLASSES) + [c for c in observed if c not in COUNTED_CLASSES]


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
        # dropna=False: a cell whose phenotype did not map is counted in a null group of its
        # own instead of being dropped.  It is the same group that supplies ``region_total``
        # three lines down, so without this a cell leaves the NUMERATOR AND THE DENOMINATOR in
        # one step and every remaining proportion silently rises.  The null group is NOT folded
        # into a named class: "this cell's class is unavailable" is not "zero cells of class X".
        # No zero ROW is written here -- the long form carries a row only where cells exist, by
        # design (see the pivot in main()), and a row per (islet, class) would be an invention
        # rather than a count.
        counts = (subset.groupby(["combined_islet_id", "islet_key", "phenotype"], dropna=False)
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


def cohort_core_summary(core: pd.DataFrame) -> pd.Series:
    """Mean per-islet proportion per class over the islet cores, every declared class printed.

    Pivot with fill_value before averaging: in the long form an islet lacking a phenotype has no
    row for it, so a plain groupby().mean() averages each class over a different set of islets and
    the "proportions" sum to well over 1.

    The pivot is then reindexed over every declared class.  A class with no cells anywhere has no
    column to pivot and would leave this summary with no line at all -- which a reader cannot tell
    from a class that was never asked about, and which is how a class leaves a denominator without
    a word.  Every islet here WAS measured for every class (the chain gives each eligible cell
    exactly one), so "none of these cells is that class" is a real 0.0.  An islet with no core
    cells is a different thing: it has no row at all and is absent, not zero.

    ``dropna=False`` is deliberately NOT passed to ``pivot_table``.  On THIS index -- one level --
    it changes nothing (measured on 115 + 6442: 1,787 x 8 either way), but on ``build_workbook``'s
    three-level index the same argument makes pandas reindex the INDEX over the cartesian product
    of its levels and invent islets (1,787 -> 4,646,200 rows).  So the null phenotype is kept by
    NAMING it :data:`UNMAPPED` before the pivot, which works on both shapes and reads as what it
    is, rather than by an argument whose effect depends on how many levels the index happens to
    have.
    """
    core = core.assign(phenotype=core["phenotype"].fillna(UNMAPPED))
    wide = core.pivot_table(index="combined_islet_id", columns="phenotype",
                            values="prop", fill_value=0.0)
    wide = wide.reindex(columns=_full_class_columns(wide.columns), fill_value=0.0)
    return wide.mean().sort_values(ascending=False)


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
    summary = cohort_core_summary(core)
    print("\ncohort mean per-islet proportion, islet core (absent counts as zero):")
    print(summary.to_string())
    print(f"  sum = {summary.sum():.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
