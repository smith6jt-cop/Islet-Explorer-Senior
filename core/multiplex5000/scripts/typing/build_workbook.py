#!/usr/bin/env python3
"""The master workbook, the per-islet proportions and the spatial lookup.

Three artifacts because the app reads the same per-islet facts three ways, and one source keeps
them consistent.  Conventions taken from the code rather than from the legacy files:

  * ``islet_key`` is never stored -- ``add_unit_key()`` re-derives it by regex from ``region``
    (Islet_Targets), ``name`` (Islet_Composition) or ``region`` (Islet_Markers), so each of those
    columns must carry a parseable ``Islet_<n>``.  The Targets sheet's band rows therefore read
    ``Islet_<n>_band``, which the same regex still resolves.
  * ``Islet_Targets.type`` is ``islet_core`` / ``islet_band`` / ``islet_union``; the app
    synthesises nothing, so the union rows are written rather than inferred.
  * Microenvironment is AREA-valued: ``area_um2 / region_um2``.  An islet with no measured
    structure area gets NA, never 0 -- 0 asserts an absence that was not measured.
  * ``case_id`` in the lookup cannot carry zero padding: ``read.csv`` type-converts per column, so
    ``0115`` reads back as the integer 115 whatever the file says.  Written unpadded.

    python core/multiplex5000/scripts/typing/build_workbook.py
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import aggregate_islets as AGG
import build_cells as CELLS
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"
STRUCTURES = ("Capillary", "Nerve", "Lymphatic")

#: Every class that can reach ``phenotype``, in the app's spelling, taken from the ONE crosswalk
#: ``build_cells.DISPLAY`` instead of being restated -- ``Unresolved (no lineage evidence)``
#: included, although ``refine.py`` is not wired and no cell on disk carries it.  The columns of
#: ``islet_celltype_proportions.csv`` ARE the app's vocabulary: a class with no cells gets no
#: pivot column at all, and a column that is not there is not a zero, it is a control the app
#: cannot offer and a reader who never learns the class was asked about.
#:
#: ``Excluded`` is not here: ``aggregate_islets`` takes those cells out of the denominator by
#: name upstream, so they are absent from this table rather than zero in it.
COUNTED_CLASSES = tuple(sorted(set(CELLS.DISPLAY.values()) - {CELLS.DISPLAY["Excluded"]}))


def _full_class_columns(observed) -> list[str]:
    """The declared classes first, then any observed class the crosswalk does not name."""
    return list(COUNTED_CLASSES) + [c for c in observed if c not in COUNTED_CLASSES]


#: Islet_Targets `type` per geometry region.  non_endocrine is this build's own addition.
TYPE_OF = {"core": "islet_core", "expansion": "islet_band",
           "union": "islet_union", "non_endocrine": "islet_non_endocrine"}
#: The region label add_unit_key() must still resolve an Islet_<n> out of.
REGION_SUFFIX = {"core": "", "expansion": "_band", "union": "_union",
                 "non_endocrine": "_nonendo"}


def _donor_meta() -> pd.DataFrame:
    meta = pd.read_csv(APP / "donor_metadata.csv")
    meta = meta.rename(columns={"case_id": "Case ID", "donor_status": "Donor Status",
                                "age": "Age", "sex": "Gender", "n_aab": "Aabs"})
    for assay in ("AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_mIAA"):
        if assay not in meta:
            meta[assay] = pd.NA
    meta["AAb_IAA"] = meta["AAb_mIAA"]      # the app checks both spellings
    keep = ["Case ID", "Donor Status", "Age", "Gender", "Aabs",
            "AAb_GADA", "AAb_IA2A", "AAb_ZnT8A", "AAb_IAA", "AAb_mIAA"]
    return meta[keep]


def build_targets(geometry: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, islet in geometry.iterrows():
        region = islet["region"]
        label = f"{islet['islet_key']}{REGION_SUFFIX[region]}"
        for structure in STRUCTURES:
            area = islet.get(f"area_um2_{structure}", np.nan)
            rows.append({
                "Case ID": int(islet["donor_id"]), "region": label, "type": TYPE_OF[region],
                "class": structure, "area_um2": area, "region_um2": islet["region_um2"],
                "area_density": (area / islet["region_um2"]
                                 if pd.notna(area) and islet["region_um2"] > 0 else np.nan),
                "count": np.nan,   # these targets are area-valued; a count would measure
                                   # segmentation morphology, not biology
            })
    frame = pd.DataFrame(rows)
    return frame.merge(meta, on="Case ID", how="left", validate="many_to_one")


def build_composition(hormones: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    frame = hormones.rename(columns={"islet_key": "name"}).copy()
    frame["Case ID"] = frame["donor_id"].astype(int)
    frame["islet_id"] = frame["name"].str.extract(r"Islet_(\d+)").astype(int)
    keep = ["Case ID", "name", "islet_id", "cells_total", "Ins_single", "Glu_single",
            "Stt_single", "Multi_Pos", "Triple_Neg", "Ins_any", "Glu_any", "Stt_any"]
    return frame[keep].merge(meta, on="Case ID", how="left", validate="many_to_one")


def build_proportions(composition: pd.DataFrame) -> pd.DataFrame:
    """One row per islet, counts and proportions wide, core and peri side by side."""
    out = []
    for region, prefix in (("core", ""), ("expansion", "peri_")):
        here = composition[composition["region"] == region]
        if here.empty:
            continue
        # A null phenotype is NAMED before the pivot rather than dropped by it.  ``region_total``
        # -- and therefore ``total_cells_core`` / ``total_cells_peri`` below -- already counts
        # those cells, so letting pivot_table drop them would leave the per-class counts no
        # longer summing to the total they are printed beside.  ``dropna=False`` is NOT the fix
        # here: on this three-level index pandas would also reindex the INDEX over the cartesian
        # product of its levels and invent islets.  Measured on donors 115 and 6442: the core
        # pivot alone goes 1,787 -> 4,646,200 rows and this function returns 4,678,752 islets
        # instead of 1,794.  The named column says "class unavailable", which is not a class and
        # not a zero.
        here = here.assign(phenotype=here["phenotype"].fillna(AGG.UNMAPPED))
        counts = here.pivot_table(index=["combined_islet_id", "donor_id", "islet_key"],
                                  columns="phenotype", values="n_cells", fill_value=0)
        props = here.pivot_table(index=["combined_islet_id", "donor_id", "islet_key"],
                                 columns="phenotype", values="prop", fill_value=0.0)
        # Reindexed over every declared class, so a class with no cells anywhere still gets its
        # column.  These column NAMES are the app's vocabulary -- ``mod_plot_server`` greps
        # ``^prop_`` and matches the stem against the configured phenotypes -- so a class with no
        # column is a control the app cannot offer, which is indistinguishable from a class that
        # was never asked about.  For an islet that HAS cells in this region the fill is a real
        # measurement: the chain gave every eligible cell exactly one class, so "none of them is
        # this class" is 0 cells and 0.0 of the region.  For an islet with NO cells in this
        # region there is no row here at all -- it arrives NaN through the join below, and the
        # block after this loop keeps the two apart: counts become 0, proportions STAY NaN,
        # because 0 % of no cells is not a measurement.
        counts = counts.reindex(columns=_full_class_columns(counts.columns), fill_value=0)
        props = props.reindex(columns=_full_class_columns(props.columns), fill_value=0.0)
        totals = here.groupby(["combined_islet_id", "donor_id", "islet_key"])["region_total"].first()
        counts.columns = [f"{prefix}count_{c}" for c in counts.columns]
        props.columns = [f"{prefix}prop_{c}" for c in props.columns]
        frame = pd.concat([totals.rename(f"total_cells_{'core' if not prefix else 'peri'}"),
                           counts, props], axis=1)
        out.append(frame)
    wide = out[0].join(out[1:], how="outer") if len(out) > 1 else out[0]
    wide = wide.reset_index().rename(columns={"donor_id": "case_id"})
    # An islet present in one region and absent from the other has genuine zeros in the counts:
    # the region exists, it simply contains no cells of that class.  A *proportion* over an empty
    # region does not -- three islets have no peri cell at all, and writing 0.0 there would state
    # "0 % of this islet's peri cells are Alpha" about an islet with no peri cells to be 0 % of,
    # pulling every Card A scatter and Card B median toward zero.  Those are unavailable.
    for column in wide.columns:
        if column.startswith(("count_", "peri_count_", "total_cells_")):
            wide[column] = wide[column].fillna(0)
    for prefix, total in (("prop_", "total_cells_core"), ("peri_prop_", "total_cells_peri")):
        columns = [c for c in wide.columns if c.startswith(prefix)]
        if not columns or total not in wide.columns:
            continue
        empty = wide[total].fillna(0) == 0
        wide.loc[~empty, columns] = wide.loc[~empty, columns].fillna(0.0)
        wide.loc[empty, columns] = float("nan")
    return wide


def build_lookup(geometry: pd.DataFrame) -> pd.DataFrame:
    core = geometry[geometry["region"] == "core"].copy()
    return pd.DataFrame({
        "case_id": core["donor_id"].astype(int),
        "imageid": core["donor_id"].astype(int),
        "islet_key": core["islet_key"],
        "combined_islet_id": core["combined_islet_id"],
        "centroid_x_um": core["centroid_x_um"].round(1),
        "centroid_y_um": core["centroid_y_um"].round(1),
        "area_um2": core["region_um2"].round(1),
    })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    geometry = pd.read_parquet(APP / "islet_geometry.parquet")
    geometry = geometry[geometry["donor_id"].isin(donors)]
    composition = pd.read_parquet(APP / "islet_composition_long.parquet")
    composition = composition[composition["donor_id"].isin(donors)]
    hormones = pd.read_parquet(APP / "islet_hormone_composition.parquet")
    hormones = hormones[hormones["donor_id"].isin(donors)]
    meta = _donor_meta()

    targets = build_targets(geometry, meta)
    comp_sheet = build_composition(hormones, meta)
    proportions = build_proportions(composition)
    lookup = build_lookup(geometry)

    workbook = APP / "master_results_run10.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        comp_sheet.to_excel(writer, sheet_name="Islet_Composition", index=False)
        targets.to_excel(writer, sheet_name="Islet_Targets", index=False)
    proportions.to_csv(APP / "islet_celltype_proportions.csv", index=False)
    lookup.to_csv(APP / "islet_spatial_lookup.csv", index=False)

    print(f"Islet_Composition {len(comp_sheet):>7,} rows")
    print(f"Islet_Targets     {len(targets):>7,} rows  ({targets['type'].nunique()} region types "
          f"x {targets['class'].nunique()} classes)")
    print(f"proportions       {len(proportions):>7,} islets x {proportions.shape[1]} cols")
    print(f"spatial lookup    {len(lookup):>7,} islets")
    print(f"\n-> {APP}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
