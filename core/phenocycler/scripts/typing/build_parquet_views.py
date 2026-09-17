#!/usr/bin/env python3
"""The four Parquet files the app's DuckDB views read.

`data_loading.R::load_master_duckdb()` is tried before the H5AD and before Excel, so these are the
fast path.  They are not a different dataset -- every value comes from the same tables the workbook
is written from, which is the point: two readers of one source rather than two sources.

What each must carry, taken from the code that reads it rather than from the legacy file:

  islets.parquet        the `.obs` mirror.  `imageid` is read as
                        `as.integer(as.character(imageid))` and `base_islet_id` through
                        `gsub("^Islet_Islet_", "Islet_", ...)`, so they are written plain: the
                        donor id and `Islet_<n>`.  Its `prop_*` columns become the phenotype
                        merge, its `age`/`gender` the donor demographics, and the
                        `peri_prop_*` / `peri_count_*` / `enrich_z_*` / `min_dist_*` /
                        `total_cells_*` group the neighbourhood merge that Spatial Cards A-C read.
  uns_targets           the Islet_Targets sheet, unchanged.
  uns_composition       the Islet_Composition sheet, unchanged.
  uns_markers           NOT written.  This dataset has no per-marker positivity, and the view is
                        registered only if the file exists -- so its absence is the honest signal
                        rather than an empty table that reads as "measured, all zero".

    python core/phenocycler/scripts/typing/build_parquet_views.py
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"


def build_islets(donors: list[str]) -> pd.DataFrame:
    proportions = pd.read_csv(APP / "islet_celltype_proportions.csv")
    proportions = proportions[proportions["case_id"].astype(str).isin(donors)]

    geometry = pd.read_parquet(APP / "islet_geometry.parquet")
    core = geometry[(geometry["region"] == "core") & geometry["donor_id"].isin(donors)]
    core = core[["donor_id", "islet_key", "region_um2", "islet_diam_um",
                 "centroid_x_um", "centroid_y_um"]]
    core = core.rename(columns={"region_um2": "core_region_um2"})
    expansion = geometry[(geometry["region"] == "expansion") & geometry["donor_id"].isin(donors)]
    expansion = expansion[["donor_id", "islet_key", "region_um2"]].rename(
        columns={"region_um2": "peri_region_um2"})

    donor_meta = pd.read_csv(APP / "donor_metadata.csv")
    donor_meta = donor_meta[["case_id", "age", "sex", "donor_status"]].rename(
        columns={"sex": "gender"})

    # The neighbourhood metrics: Card A's enrichment, Card B's z bars and heatmap, Card C's
    # left plot.  It recomputes peri_prop_/peri_count_/total_cells_* from the same per-cell
    # table the proportions come from, so the overlap is an independent second computation of
    # the same quantity -- asserted equal here, then dropped, so the merge adds only what is new.
    neighborhood = pd.read_csv(APP / "neighborhood_metrics.csv")
    neighborhood = neighborhood[neighborhood["case_id"].astype(str).isin(donors)]

    frame = proportions.copy()
    frame["donor_id"] = frame["case_id"].astype(str)
    frame = frame.merge(core, on=["donor_id", "islet_key"], how="left", validate="one_to_one")
    frame = frame.merge(expansion, on=["donor_id", "islet_key"], how="left", validate="one_to_one")
    frame = frame.merge(donor_meta, on="case_id", how="left", validate="many_to_one")

    shared = [c for c in neighborhood.columns
              if c in frame.columns and c not in ("case_id", "islet_key")]
    check = frame[["case_id", "islet_key"] + shared].merge(
        neighborhood[["case_id", "islet_key"] + shared],
        on=["case_id", "islet_key"], suffixes=("_prop", "_nbr"), validate="one_to_one")
    for column in shared:
        left, right = check[f"{column}_prop"], check[f"{column}_nbr"]
        if not left.astype(float).round(9).equals(right.astype(float).round(9)):
            raise ValueError(f"{column}: the proportions and neighbourhood tables disagree")
    frame = frame.merge(neighborhood.drop(columns=shared), on=["case_id", "islet_key"],
                        how="left", validate="one_to_one")

    # Islet-level Leiden + its UMAP, from trajectory_islet.py. The Spatial tab's "Islet-Level
    # UMAPs" card greps `^leiden_` off the islet table, which is why these ride here rather than
    # only in the trajectory H5AD.
    leiden_path = APP / "islet_leiden.parquet"
    if leiden_path.exists():
        leiden = pd.read_parquet(leiden_path)
        leiden = leiden[leiden["case_id"].astype(str).isin(donors)].copy()
        # The proportions table types case_id as int; these writers type it from a string donor id.
        # Merging int64 against str raises rather than silently matching nothing, which is the
        # better failure, but both sides still have to agree.
        leiden["case_id"] = leiden["case_id"].astype(int)
        frame = frame.merge(leiden, on=["case_id", "islet_key"], how="left",
                            validate="one_to_one")

    # The galectin layer (correction J), as per-islet proportions so it reaches the Plot tab's
    # measure list rather than only a sheet. Class names are distinct from the eight phenotype
    # names, which is what lets the app group them separately -- these are QuPath object-classifier
    # verdicts, not subpanel classes, and the two must not read as one list.
    galectin_path = APP / "islet_galectin.parquet"
    if galectin_path.exists():
        galectin = pd.read_parquet(galectin_path)
        galectin = galectin[galectin["case_id"].astype(str).isin(donors)].copy()
        galectin["case_id"] = galectin["case_id"].astype(int)
        classes = [c[len("prop_"):] for c in galectin.columns if c.startswith("prop_")]
        for region, prefix in (("core", "prop_"), ("peri", "peri_prop_")):
            here = galectin[galectin["region"] == region]
            here = here.set_index(["case_id", "islet_key"])[[f"prop_{c}" for c in classes]]
            here.columns = [f"{prefix}{c}" for c in classes]
            frame = frame.merge(here.reset_index(), on=["case_id", "islet_key"],
                                how="left", validate="one_to_one")

    # The two identity columns the loader reads, spelled the way it expects.
    frame["imageid"] = frame["case_id"].astype(int)
    frame["base_islet_id"] = frame["islet_key"]
    return frame.drop(columns=["donor_id"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    target = APP / "parquet"
    target.mkdir(parents=True, exist_ok=True)

    islets = build_islets(donors)
    islets.to_parquet(target / "islets.parquet", index=False)

    workbook = pd.ExcelFile(APP / "master_results_run10.xlsx")
    sheets = {"targets": "Islet_Targets", "composition": "Islet_Composition"}
    written = {}
    for view, sheet in sheets.items():
        frame = workbook.parse(sheet)
        frame.to_parquet(target / f"uns_{view}.parquet", index=False)
        written[view] = len(frame)

    # If a stale uns_markers is lying around from another build it would register as a view and
    # read as "measured, all zero"; absence is the honest signal.
    stale = target / "uns_markers.parquet"
    if stale.exists():
        stale.unlink()
        print("  removed a stale uns_markers.parquet -- this dataset has no per-marker positivity")

    families = {name: len([c for c in islets.columns if c.startswith(name)])
                for name in ("prop_", "peri_prop_", "enrich_z_", "min_dist_")}
    print(f"islets.parquet        {len(islets):>7,} islets x {islets.shape[1]} cols ("
          + ", ".join(f"{n} {name}*" for name, n in families.items()) + ")")
    for view, n in written.items():
        print(f"uns_{view:<16s} {n:>7,} rows")
    print(f"uns_markers           not written (no per-marker positivity in this dataset)")
    print(f"\n-> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
