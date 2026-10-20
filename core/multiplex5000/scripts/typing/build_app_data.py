#!/usr/bin/env python3
"""Write the app's data artifacts from the per-cell table.

Each writer is built against a contract derived from the app's own R source rather than from the
legacy files, because an exemplar shows one valid instance while the code shows what is actually
required.  Three conventions that only the code reveals, and that a copy of the exemplar would
have got wrong:

  * The app's region vocabulary is core / **peri** / tissue.  This lane calls the middle one
    ``expansion``; the app artifacts use the app's word.
  * ``islet_key`` is never stored in the workbook -- ``add_unit_key()`` re-derives it from
    ``region`` / ``name`` by regex, so those columns must carry a parseable ``Islet_<n>``.
  * A drill-down CSV with a header and no rows makes ``load_islet_cells()`` return NULL and the
    panel falls back to boundaries silently, so an empty islet gets no file at all.

``object_id`` IS THE FIRST COLUMN OF ``cells/`` AND ``donors/``.  It was absent from both until
2026-09-18, which meant a drill-down cell traced back to nothing: the app could draw it and never
say which QuPath detection it was.  The repo rule is that the UUID is preserved throughout and is
never stood in for by centroid matching, row order or a filename, so it is written, and it is
asserted unique and non-null within the donor before either writer emits a row.

``parquet/tissue/`` and ``parquet/cell_distance/`` still carry no ``object_id``: they are the
DuckDB aggregate views, the app selects their columns by name, and nothing has asked to trace a
bin back to a cell.  Stated here so its absence reads as a decision rather than an oversight.

    python core/multiplex5000/scripts/typing/build_app_data.py [--only tissue donors cells outlines]
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

APP = S.OUT / "app"
#: The app's own word for the 20 um zone; this lane says "expansion" everywhere else.
REGION_TO_APP = {"core": "core", "expansion": "peri", "tissue": "tissue"}


def _cells(donor: str, columns: list[str]) -> pd.DataFrame:
    return pq.read_table(S.OUT / f"cells/donor_id={donor}/data.parquet",
                         columns=columns).to_pandas()


def check_object_ids(donor: str, ids: pd.Series, artifact: str) -> None:
    """Refuse to write a per-cell artifact whose UUIDs cannot identify a cell.

    Two separate failures, reported separately because they have different causes.  A null
    ``object_id`` means the upstream table lost the UUID for that row; a duplicated one means
    two rows claim the same QuPath detection, and a reader joining on it would silently
    multiply cells.  Either makes the column useless as the join key it exists to be, so the
    writer stops rather than emitting a file that looks traceable and is not.
    """
    missing = int(ids.isna().sum()) + int((ids.astype("string").fillna("") == "").sum())
    if missing:
        raise ValueError(f"{donor}: {missing:,} of {len(ids):,} {artifact} rows carry no "
                         "object_id; the UUID is the only handle back to the QuPath detection")
    duplicated = int(ids.duplicated().sum())
    if duplicated:
        example = ids[ids.duplicated(keep=False)].iloc[0]
        raise ValueError(f"{donor}: {duplicated:,} duplicated object_id in {artifact} "
                         f"(e.g. {example!r}); one UUID must name one cell")


def _app_frame(donor: str, *, with_object_id: bool = False) -> pd.DataFrame:
    """The columns every per-cell app artifact shares, in the app's vocabulary.

    ``with_object_id`` puts the QuPath UUID first.  It is checked here, on the frame that is
    actually written, rather than on the source table: the ``Excluded`` filter runs first, and
    a duplicate that both rows survive is a different fact from one where only one row does.
    """
    columns = ["X_centroid", "Y_centroid", "phenotype", "cell_region_20um",
               "islet_key", "cell_class"]
    frame = _cells(donor, (["object_id"] if with_object_id else []) + columns)
    frame = frame[frame["cell_class"] != "Excluded"].copy()
    frame["cell_region"] = frame["cell_region_20um"].map(REGION_TO_APP)
    frame["islet_name"] = frame["islet_key"]
    keep = ["X_centroid", "Y_centroid", "phenotype", "cell_region", "islet_name"]
    if with_object_id:
        check_object_ids(donor, frame["object_id"], "donors")
        keep = ["object_id"] + keep
    return frame[keep]


def write_tissue_parquet(donors: list[str]) -> None:
    """parquet/tissue/case_id=<donor>/ -- the DuckDB `tissue` view, one row per cell."""
    root = APP / "parquet/tissue"
    for donor in donors:
        frame = _app_frame(donor)
        target = root / f"case_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pandas(frame, preserve_index=False),
                       target / "part-0.parquet", compression="zstd")
        print(f"  tissue  {donor}: {len(frame):>9,} cells", flush=True)


def write_drilldown_cells(donors: list[str]) -> None:
    """cells/<case_id>_<islet_key>.csv -- core + peri cells of one islet, identity by filename."""
    root = APP / "cells"
    root.mkdir(parents=True, exist_ok=True)
    markers = [m for m in S.load_scale().index if m not in ("DAPI", "IAPP")]
    written = skipped = 0
    for donor in donors:
        frame = _cells(donor, ["object_id", "X_centroid", "Y_centroid", "phenotype",
                               "cell_region_20um", "islet_key", "cell_class"])
        frame = frame[(frame["cell_class"] != "Excluded") & frame["islet_key"].notna()
                      & frame["cell_region_20um"].isin(["core", "expansion"])]
        if frame.empty:
            continue
        # Before the merge, not after: `validate="one_to_one"` would raise on a duplicated
        # UUID first, and its message says nothing about which donor or how many.
        check_object_ids(donor, frame["object_id"], "cells")
        have = [m for m in markers
                if m in pq.ParquetFile(S.ALIGNED / f"donor_id={donor}/data.parquet").schema_arrow.names]
        values = pq.read_table(S.ALIGNED / f"donor_id={donor}/data.parquet",
                               columns=["object_id"] + have).to_pandas()
        frame = frame.merge(values, on="object_id", how="left", validate="one_to_one")
        frame["cell_region"] = frame["cell_region_20um"].map(REGION_TO_APP)
        keep = ["object_id", "X_centroid", "Y_centroid", "phenotype", "cell_region"] + have
        for islet, group in frame.groupby("islet_key"):
            if group.empty:            # a header-only file reads back as NULL, so write none
                skipped += 1
                continue
            group[keep].to_csv(root / f"{donor}_{islet}.csv", index=False)
            written += 1
        print(f"  cells   {donor}: {frame['islet_key'].nunique():>5,} islets", flush=True)
    print(f"  -> {written:,} drill-down files ({skipped} empty islets skipped)")


def write_outlines(donors: list[str]) -> None:
    """json/<case_id>.geojson and gson/<case_id>.geojson.gz -- islet and expansion outlines.

    The app matches `classification.name` case-insensitively and expects the expansion object to
    be named `<islet_key>_exp20um`, so both are carried through from the QuPath export unchanged.
    """
    raw = json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    stem = {r["donor_id"]: r["image_id"].split(" - ")[0].replace(".qptiff", "")
            for r in raw["images"]}
    (APP / "json").mkdir(parents=True, exist_ok=True)
    (APP / "gson").mkdir(parents=True, exist_ok=True)
    for donor in donors:
        matches = [p for p in S.ISLETS.glob("*.geojson") if stem[donor] in p.name]
        if len(matches) != 1:
            raise FileNotFoundError(f"{donor}: {len(matches)} islet outline files")
        target = APP / f"json/{donor}.geojson"
        collection = json.loads(matches[0].read_text())
        # Correction E's outline is ours, not QuPath's, so it is merged in here rather than
        # expected in the export. Absent (no endocrine union built for this donor) the file is
        # the QuPath one unchanged and the viewer simply has one fewer class to draw.
        extra = S.OUT / f"geometry/nonendocrine_{donor}.geojson"
        if extra.exists():
            added = json.loads(extra.read_text())["features"]
            collection["features"] = list(collection["features"]) + added
        target.write_text(json.dumps(collection))
        with open(target, "rb") as source, gzip.open(APP / f"gson/{donor}.geojson.gz", "wb") as out:
            shutil.copyfileobj(source, out)
        classes = {}
        for feature in json.loads(target.read_text())["features"]:
            name = (feature.get("properties", {}).get("classification") or {}).get("name")
            classes[name] = classes.get(name, 0) + 1
        print(f"  outline {donor}: {classes}", flush=True)


def write_cell_distance_parquet(donors: list[str]) -> None:
    """parquet/cell_distance/case_id=<donor>/ -- the DuckDB `cell_distance` view.

    Card C's right plot bins one signed distance per cell, so the view is that column and the
    three it groups by.  ``signed_dist_islet_um`` already carries the app's sign convention --
    negative inside the islet, zero at the border -- and is the same column correction F derives
    ``cell_region_20um`` from, so the region and the distance cannot disagree.

    The legacy view also carried ten ``<marker>_pos`` booleans from RESTORE.  This cohort has no
    per-marker positivity call: the subpanel classification decides a cell's identity from the
    markers that can answer its context's question, and never asks "is this cell CD31-positive"
    as a standalone question.  Writing the columns would mean inventing that call, so they are
    absent and the app's marker mode says so.
    """
    root = APP / "parquet/cell_distance"
    status = pd.read_csv(APP / "donor_metadata.csv", dtype={"case_id": str})
    status = status.set_index("case_id")["donor_status"]
    for donor in donors:
        frame = _cells(donor, ["phenotype", "cell_region_20um", "cell_class",
                               "signed_dist_islet_um"])
        frame = frame[frame["cell_class"] != "Excluded"].copy()
        frame["donor_status"] = status[donor]
        frame["cell_region"] = frame["cell_region_20um"].map(REGION_TO_APP)
        frame = frame.rename(columns={"signed_dist_islet_um": "dist_islet"})
        frame = frame[["donor_status", "cell_region", "phenotype", "dist_islet"]]
        frame = frame[frame["dist_islet"].notna()]
        target = root / f"case_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        pq.write_table(pa.Table.from_pandas(frame, preserve_index=False),
                       target / "part-0.parquet", compression="zstd")
        inside = int((frame["dist_islet"] <= 0).sum())
        print(f"  distance {donor}: {len(frame):>9,} cells, {inside:>8,} inside an islet",
              flush=True)


WRITERS = {"tissue": write_tissue_parquet,
           "cells": write_drilldown_cells, "outlines": write_outlines,
           "distance": write_cell_distance_parquet}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--only", nargs="*", choices=sorted(WRITERS), default=sorted(WRITERS))
    args = parser.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    APP.mkdir(parents=True, exist_ok=True)
    for name in args.only:
        print(f"\n== {name}")
        WRITERS[name](donors)
    print(f"\n-> {APP}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
