#!/usr/bin/env python3
"""Assemble the per-cell table every app artifact is aggregated from.

Joins the subpanel calls to the run's own coordinates and QC by ``(donor_id, object_id)`` -- never
row order, never a rounded centroid, never a filename -- and applies the corrections that the
calls on disk do not carry:

  G  no nucleus, no cell.  The calls already cover only ``qc_analysis_eligible`` cells, and a
     detection without a nucleus is never eligible, so this is satisfied by construction rather
     than by a filter.  The dropped count is recorded per donor so the denominator stays auditable.
  F  every peri-islet expansion is 20 um.  Derived from the signed distance rather than read from
     the annotation, because one donor's annotation is drawn at 23.1 um and is named ``_exp20um``
     like every other, so nothing downstream can see the difference.
  C  which object a cell sits in, recorded.  The classification already uses object membership --
     including inside islets -- but ``context`` overwrites it with ``islet`` for every islet cell,
     so the evidence the call was made on is not in the output.  Kept as its own column.

Writes one partition per donor, so a donor can be rebuilt without touching the others.

    python core/phenocycler/scripts/typing/build_cells.py [--donors 6442 ...] [--include-marginal]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import context as CX
import subpanel_io as S

DISTANCE = S.REPO / "data/development/kintsugi/distance"
OUT = S.OUT / "cells"

#: Display label per call, and the region a cell is reported in.  ``Excluded`` keeps its reason.
DISPLAY = {
    "Alpha": "Alpha_cell", "Beta": "Beta_cell", "Delta": "Delta_cell",
    "Blood vessel": "Blood_Vessel", "Lymphatic": "Lymphatic", "Nerve": "Neural",
    "Immune": "Immune", "Exocrine": "Exocrine", "Excluded": "Excluded",
}
ENDOCRINE = ("Alpha", "Beta", "Delta")


def _manifest() -> tuple[dict, dict]:
    raw = json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    pixel = {r["donor_id"]: r["pixel_size_um_x"] for r in raw["images"]}
    stem = {r["donor_id"]: r["image_id"].split(" - ")[0].replace(".qptiff", "")
            for r in raw["images"]}
    return pixel, stem


def build_donor(donor: str, pixel: dict, stem: dict, runs: pd.DataFrame) -> pd.DataFrame:
    calls = pq.read_table(
        S.OUT / "subpanel/subpanel_calls.parquet",
        columns=["donor_id", "object_id", "cell_class", "context", "excluded_reason"],
        filters=[("donor_id", "=", donor)]).to_pandas()

    expression = pq.read_table(
        S.RUN / f"expression/donor_id={donor}/data.parquet",
        columns=["object_id", "X_centroid", "Y_centroid", "islet_num",
                 "cell_region", "qc_analysis_eligible"]).to_pandas()

    geometry = pq.read_table(
        S.RUN / f"geometry_qc/donor_id={donor}/data.parquet",
        columns=["object_id", "analysis_reason"]).to_pandas()

    distance = pq.read_table(
        DISTANCE / f"donor_id={donor}/data_0.parquet",
        columns=["object_id", "signed_dist_islet_um"]).to_pandas()

    # --- correction G, as an accounting not a filter -------------------------------------
    no_nucleus = int((geometry["analysis_reason"] == "missing_nucleus_geometry").sum())
    ineligible = int((~expression["qc_analysis_eligible"]).sum())
    if not expression.loc[expression["object_id"].isin(calls["object_id"]),
                          "qc_analysis_eligible"].all():
        raise AssertionError(f"{donor}: a typed cell is not analysis-eligible -- correction G "
                             "cannot be satisfied by construction, so it must become a filter")

    frame = calls.merge(expression, on="object_id", how="left", validate="one_to_one")
    frame = frame.merge(distance, on="object_id", how="left", validate="one_to_one")
    if frame["X_centroid"].isna().any():
        raise AssertionError(f"{donor}: {int(frame['X_centroid'].isna().sum())} typed cells have "
                             "no coordinate -- the join on object_id did not resolve")
    unmatched = int(frame["signed_dist_islet_um"].isna().sum())
    if unmatched:
        raise AssertionError(f"{donor}: {unmatched} typed cells have no signed distance")

    # --- correction F: the expansion is 20 um by construction ----------------------------
    d = frame["signed_dist_islet_um"].to_numpy()
    frame["cell_region_20um"] = np.where(d <= 0, "core", np.where(d <= 20, "expansion", "tissue"))

    # --- correction C: keep the object the call was actually made on ----------------------
    frame["object_context"] = CX.assign(
        frame["X_centroid"].to_numpy(), frame["Y_centroid"].to_numpy(),
        np.zeros(len(frame), dtype=bool),   # ask the object question of every cell, islet included
        stem[donor], pixel[donor])

    frame["phenotype"] = frame["cell_class"].map(DISPLAY)
    if frame["phenotype"].isna().any():
        missing = sorted(frame.loc[frame["phenotype"].isna(), "cell_class"].unique())
        raise KeyError(f"{donor}: no display label for {missing}")
    frame["is_endocrine"] = frame["cell_class"].isin(ENDOCRINE)
    # islet_num is a string in the run, empty for a cell in no islet -- not a nullable int.
    islet_num = frame["islet_num"].fillna("").astype(str).str.strip()
    frame["islet_num"] = islet_num.replace("", None)
    frame["islet_key"] = np.where(islet_num != "", "Islet_" + islet_num, None)
    frame["combined_islet_id"] = np.where(islet_num != "", f"{donor}_Islet_" + islet_num, None)

    run = runs.loc[runs["donor_id"] == donor]
    if run.empty:
        raise KeyError(f"{donor}: absent from donor_acquisition_runs.csv")
    frame["acq_run"] = run["acq_run"].iloc[0]
    frame["stain_group"] = run["stain_group"].iloc[0]

    frame.attrs["dropped_no_nucleus"] = no_nucleus
    frame.attrs["dropped_ineligible"] = ineligible
    return frame


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    pixel, stem = _manifest()
    runs = pd.read_csv(S.REPO / "data/donor_acquisition_runs.csv", dtype={"donor_id": str})

    OUT.mkdir(parents=True, exist_ok=True)
    ledger = []
    for donor in donors:
        frame = build_donor(donor, pixel, stem, runs)
        target = OUT / f"donor_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pandas(frame, preserve_index=False),
                       target / "data.parquet", compression="zstd")
        counts = frame["cell_region_20um"].value_counts()
        ledger.append({
            "donor_id": donor, "typed_cells": len(frame),
            "dropped_no_nucleus": frame.attrs["dropped_no_nucleus"],
            "dropped_ineligible": frame.attrs["dropped_ineligible"],
            "core": int(counts.get("core", 0)), "expansion": int(counts.get("expansion", 0)),
            "tissue": int(counts.get("tissue", 0)),
            "in_object": int((frame["object_context"] != "tissue").sum()),
            "endocrine": int(frame["is_endocrine"].sum()),
        })
        print(f"  {donor}: {len(frame):>9,} cells "
              f"(core {ledger[-1]['core']:,} / expansion {ledger[-1]['expansion']:,})", flush=True)

    book = pd.DataFrame(ledger)
    book.to_csv(S.OUT / "cells/build_ledger.csv", index=False)
    print(f"\n{len(donors)} donors, {book.typed_cells.sum():,} cells -> {OUT}")
    print(f"excluded and recorded: {book.dropped_no_nucleus.sum():,} without a nucleus, "
          f"{book.dropped_ineligible.sum():,} ineligible overall")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
