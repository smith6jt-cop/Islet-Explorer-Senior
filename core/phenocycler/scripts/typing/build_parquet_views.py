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

islets.parquet is a MIRROR, not a projection: nothing here selects a subset of what its sources
carry, so a column present in a source and missing from the file is a defect, not a design.

  Fixed 2026-09-18.  Two islet-level rollups -- `crossed_*` from `summarise_crossed_fraction.py`
  and `sc_*` from `summarise_singlecell_trajectory.py` -- merged themselves into the
  islets.parquet this script had already written.  Rewriting the file from scratch therefore
  removed them, and nothing logged or raised: measured on the file of 2026-09-17 18:04, 104
  columns on disk against 83 a rerun produced, 21 gone (14 `crossed_*`, 7 `sc_*`).  A consumer
  could not tell a column that was never produced from one this script had wiped, and
  `scripts/senior/build_layers.py` had to encode a node ordering to work around it.  Both
  rollups are now read from their own summary files here (:data:`ROLLUPS`), the same merge on the
  same keys those scripts do, so the two orders agree.  Anything still missing at write time is
  named -- refused by default, recorded in `islets_columns.json` beside the file either way.

    python core/phenocycler/scripts/typing/build_parquet_views.py
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import pathlib
import sys

import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"

#: Islet-level rollups written by their own scripts and merged in HERE, not only by those scripts
#: afterwards.  Each script still drops its own stale columns and remerges, so running one after
#: this is a no-op on the values; reading the summary file is what makes the order not matter.
#: Keyed `(family, label, source, producer)`.
ROLLUPS = (
    ("crossed_", "crossed fractions",
     S.OUT / "crossed/crossed_islet_summary.parquet",
     "summarise_crossed_fraction.py"),
    ("sc_", "single-cell pseudotime",
     S.OUT / "trajectory/singlecell_islet_summary.parquet",
     "summarise_singlecell_trajectory.py"),
)

#: Written beside islets.parquet.  A reader of the Parquet sees only the columns that are there;
#: this says which rollups went in, and names the ones that did not and why.
COLUMN_MANIFEST = "islets_columns.json"


def merge_rollups(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    """Merge every rollup in :data:`ROLLUPS`, and report each one present or absent BY NAME.

    Absent is a first-class outcome: the family is recorded with the path that does not exist and
    the script that would write it, rather than the columns simply not appearing.  No value is
    computed here -- the summary files carry their own thinning (an islet under the producer's
    minimum cell count is already NaN in them), and this reads them as written.
    """
    notes: list[dict] = []
    for family, label, path, producer in ROLLUPS:
        if not path.exists():
            notes.append({"family": f"{family}*", "label": label, "present": False,
                          "source": str(path), "producer": producer,
                          "reason": f"{path.name} does not exist -- run {producer}"})
            continue
        rollup = pd.read_parquet(path)
        # The proportions table types case_id as int; these writers type it from a string donor id.
        rollup["case_id"] = rollup["case_id"].astype(int)
        before = set(frame.columns)
        frame = frame.merge(rollup, on=["case_id", "islet_key"], how="left", validate="one_to_one")
        notes.append({"family": f"{family}*", "label": label, "present": True,
                      "source": str(path), "producer": producer,
                      "columns": [c for c in frame.columns if c not in before]})
    return frame, notes


def refuse_silent_drop(islets: pd.DataFrame, target: pathlib.Path,
                       allow_drop: bool) -> list[str]:
    """Name every column the islets.parquet being replaced carries and this build does not.

    Three scripts write this one path.  Overwriting it with fewer columns than it had is how the
    21 `crossed_*`/`sc_*` columns used to vanish, and the only thing wrong with that was the
    silence: a consumer reading the view could not tell a wiped column from one never produced.
    So the overwrite is refused by default and the names are printed; `--allow-drop` proceeds and
    records them in the manifest.  Whether a dropped column SHOULD come back is not decided here.
    """
    previous = target / "islets.parquet"
    if not previous.exists():
        return []
    have = set(islets.columns)
    dropped = [c for c in pq.read_schema(previous).names if c not in have]
    if not dropped:
        return []
    named = "\n    ".join(dropped)
    if not allow_drop:
        raise SystemExit(
            f"REFUSING to overwrite {previous}: it carries {len(dropped)} column(s) this build "
            f"does not produce, and replacing it would remove them with no record.\n"
            f"    {named}\n"
            "Run the producer that writes them (see ROLLUPS in this file) and rerun, or pass "
            "--allow-drop to overwrite anyway -- the names are then recorded in "
            f"{COLUMN_MANIFEST}.")
    print(f"  DROPPING {len(dropped)} column(s) the islets.parquet being replaced carries: "
          + ", ".join(dropped))
    return dropped


def build_islets(donors: list[str]) -> tuple[pd.DataFrame, list[dict]]:
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
    frame = frame.drop(columns=["donor_id"])

    # The two rollups this script used to write out from under.  Last, so their columns sit where
    # the merging scripts leave them and the file reads the same whichever ran last.
    return merge_rollups(frame)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--allow-drop", action="store_true",
                        help="overwrite islets.parquet even though it carries columns this build "
                             "does not produce; they are named on stdout and in "
                             f"{COLUMN_MANIFEST} rather than disappearing")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    target = APP / "parquet"
    target.mkdir(parents=True, exist_ok=True)

    islets, rollups = build_islets(donors)
    dropped = refuse_silent_drop(islets, target, args.allow_drop)
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

    # What is in the file, and what is not and why.  stdout scrolls away; this sits beside the
    # Parquet, so a consumer that finds no `sc_pseudotime_core` can tell "never produced" from
    # "produced and dropped here".
    (target / COLUMN_MANIFEST).write_text(json.dumps({
        "written": dt.datetime.now().isoformat(timespec="seconds"),
        "donors": donors,
        "rows": int(len(islets)),
        "columns": list(islets.columns),
        "rollups": rollups,
        "dropped_from_previous": dropped,
    }, indent=2) + "\n")

    families = {name: len([c for c in islets.columns if c.startswith(name)])
                for name in ("prop_", "peri_prop_", "enrich_z_", "min_dist_",
                             "crossed_", "sc_")}
    print(f"islets.parquet        {len(islets):>7,} islets x {islets.shape[1]} cols ("
          + ", ".join(f"{n} {name}*" for name, n in families.items()) + ")")
    for note in rollups:
        if note["present"]:
            print(f"  {note['family']:<10s} {len(note['columns']):>2} cols from "
                  f"{pathlib.Path(note['source']).name}")
        else:
            print(f"  {note['family']:<10s} NOT IN THE VIEW -- {note['reason']}")
    for view, n in written.items():
        print(f"uns_{view:<16s} {n:>7,} rows")
    print(f"uns_markers           not written (no per-marker positivity in this dataset)")
    print(f"\n-> {target}   (columns and absences: {COLUMN_MANIFEST})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
