#!/usr/bin/env python3
"""Which islet a cell is in, decided by geometry rather than by QuPath's `Parent` column.

WHY THIS EXISTS.

`context.py` decided islet membership from QuPath's `Parent`, on the stated reasoning that the
Parent "is what drew the islet in the first place". Measured 2026-09-15, the Parent disagrees with
the geometry: on 6539's `Islet_387` and four siblings, 1,073 of 1,678 cells fall inside an islet
CORE outline and the signed distance agrees at up to -74.6 um inside, while `Parent` labels all
1,678 `peri`.

The consequence was not a mislabelled region. `in_islet` false meant the cell was never asked the
islet question at all, so it fell through the chain to the complement and was called **Exocrine**:
11,556 islet-core cells cohort-wide (2.7 %; 35.3 % on 6539), of which 5,250 clear the hormone
reference level. That is how a 335-cell islet ends up 90 % Exocrine while 92.5 % of its cells are
hormone-positive, and it is why an islet could appear to contain no endocrine cell at all -- a
state the maintainer rules impossible, and it is: it was this bug.

WHAT IS AUTHORITATIVE HERE, AND WHY.

The islet OUTLINES are, because they answer both questions at once -- whether a cell is in an
islet, and WHICH islet -- from one test. The signed distance answers only the first. So
containment decides, and the signed distance is carried as an independent cross-check that is
reported every run rather than assumed: two sources derived from the same outlines by different
code, disagreeing, means one of them is broken and the run says so.

Note what is NOT claimed. This does not say the outlines are correct -- a classifier that draws a
sub-cell fragment still draws one, and 146 islets in this cohort contain no cell at all. It says
only that membership follows the outline the islet was drawn as, consistently, everywhere.

    python core/phenocycler/scripts/typing/islet_membership.py [--donors 6539] [--report]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from shapely import STRtree, points
from shapely.geometry import shape

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import build_islet_geometry as GEOM
import cohort as COHORT
import subpanel_io as S

DISTANCE = S.REPO / "data/development/kintsugi/distance"
OUT = S.OUT / "islet_membership"

def _manifest() -> tuple[dict, dict]:
    raw = json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    return ({r["donor_id"]: r["pixel_size_um_x"] for r in raw["images"]},
            {r["donor_id"]: r["image_id"].split(" - ")[0].replace(".qptiff", "")
             for r in raw["images"]})


def outlines(stem: str) -> dict:
    """`islet number -> actual core geometry`, in pixels.

    The file is chosen by ``build_islet_geometry.islet_export``, which searches the export
    directory RECURSIVELY, takes the newest matching annotation export and refuses a detection
    export with its reason printed. A flat, non-recursive glob was used here until 2026-09-17,
    and it could not see the re-exports that arrive in dated subdirectories -- so this lane and
    ``islet_geometry`` could read two different sets of outlines for the same donor and nothing
    would say so. One resolver, one answer.
    """
    core = {}
    for feature in json.loads(GEOM.islet_export(stem).read_text())["features"]:
        properties = feature.get("properties") or {}
        name = str(properties.get("name") or "")
        if not name.startswith("Islet_") or not feature.get("geometry"):
            continue
        geometry = shape(feature["geometry"])
        if not geometry.is_valid:
            geometry = geometry.buffer(0)
        if name.endswith("_exp20um"):
            raise ValueError(f"{stem}: stale expansion outline {name}; regenerate the no-expansion export")
        number = name[len("Islet_"):]
        core[number] = geometry
    if not core:
        raise ValueError(f"{stem}: no `Islet_N` core outlines in "
                         f"{GEOM.islet_export(stem).name}")
    return core


def _contain(store: dict, pts) -> tuple[np.ndarray, np.ndarray, int]:
    """(hit mask, islet number per point, how many points landed in more than one outline).

    Where outlines overlap, the SMALLEST containing one wins. A cell inside both a large islet and
    a small one drawn within it belongs to the small one; taking the first match instead would
    depend on the order features happen to sit in the file.
    """
    keys = list(store)
    geoms = [store[k] for k in keys]
    area = np.array([g.area for g in geoms])
    number = np.full(len(pts), None, dtype=object)
    best = np.full(len(pts), np.inf)
    hit = np.zeros(len(pts), dtype=bool)
    if not geoms:
        return hit, number, 0
    point_idx, geom_idx = STRtree(geoms).query(pts, predicate="within")
    multi = 0
    for p, g in zip(point_idx, geom_idx):
        if hit[p]:
            multi += 1
        if area[g] < best[p]:
            best[p] = area[g]
            number[p] = keys[g]
            hit[p] = True
    return hit, number, multi


def membership(donor: str, x_um: np.ndarray, y_um: np.ndarray,
               pixel_size_um: float, stem: str) -> pd.DataFrame:
    """Per cell: whether its centroid is inside an endocrine-object core and which object.

    Only true core membership is categorical. Signed distance remains continuous in the cell
    table; it is not converted into a peri-object band here.
    """
    core = outlines(stem)
    pts = points(np.column_stack([np.asarray(x_um) / pixel_size_um,
                                  np.asarray(y_um) / pixel_size_um]))
    in_core, core_num, multi_core = _contain(core, pts)
    number = np.where(in_core, core_num, None)
    metadata = {}
    for feature in json.loads(GEOM.islet_export(stem).read_text())["features"]:
        props = feature.get("properties") or {}
        key = str(props.get("name", "")).removeprefix("Islet_")
        metadata[key] = {"endocrine_object_class": props.get("size_class", "Islet"),
                         "endocrine_object_id": feature.get("id"),
                         "qupath_object_name": props.get("qupath_name", props.get("name"))}
    frame = pd.DataFrame({
        "in_islet": in_core,
        "islet_num": number,
    })
    for column in ("endocrine_object_class", "endocrine_object_id", "qupath_object_name"):
        frame[column] = [metadata.get(str(k), {}).get(column) for k in number]
    frame.attrs["n_core_outlines"] = len(core)
    frame.attrs["overlapping_core"] = multi_core
    return frame


#: Columns every cell source carries, whatever it is.
_CELL_COLUMNS = ["object_id", "X_centroid", "Y_centroid", "islet_num", "cell_region",
                 "qc_analysis_eligible"]


def cells_with_distance(donor: str) -> tuple[pd.DataFrame, str, str]:
    """(cells, which cell source, which distance source) for one donor.

    WHY THE SOURCE IS NOT FIXED. This lane read the historical run tree's expression partition
    (``data/runs/<run_id>/expression``, where the 20-donor cohort output still lives) and the
    2026-08-27 ``kintsugi/distance`` layer unconditionally, and on the seven re-exported
    donors both predate the re-export -- the annotations they were computed against are not the
    annotations `outlines()` now reads. Measured on 6450 (2026-09-17): against the run plus
    kintsugi the containment and the signed distance disagree on 12,232 cells at a median
    |d| of 265.79 um; against the revised export and its own re-run signed distance they
    disagree on 1,131 at 0.42 um, which is the regime every batch-1 donor sits in (0.24-0.46 um).
    That disagreement is the reason those donors were never built here, and it is a stale input
    rather than a missing one.

    So the cell geometry comes from ``S.cell_source`` -- the revised rebuild where one exists,
    the run otherwise -- and the signed distance comes from whichever program computed it
    against the SAME annotations: the revised export's own column, or kintsugi for the run.
    The cross-check stays independent either way: containment is decided by this file's own
    code from the outlines, the distance by another program entirely.
    """
    path, which = S.cell_source(donor)
    have = set(pq.read_schema(path).names)
    if "signed_dist_islet_um" in have:
        cells = pq.read_table(
            path, columns=_CELL_COLUMNS + ["signed_dist_islet_um"]).to_pandas()
        distance_source = f"{which} export"
    else:
        cells = pq.read_table(path, columns=_CELL_COLUMNS).to_pandas()
        distance = pq.read_table(DISTANCE / f"donor_id={donor}/data_0.parquet",
                                 columns=["object_id", "signed_dist_islet_um"]).to_pandas()
        cells = cells.merge(distance, on="object_id", how="left", validate="one_to_one")
        distance_source = "kintsugi"
    if cells["signed_dist_islet_um"].isna().any():
        raise AssertionError(
            f"{donor}: {int(cells['signed_dist_islet_um'].isna().sum())} cells have no signed "
            f"distance in {distance_source} -- the cross-check cannot run")
    return cells, which, distance_source


def build_donor(donor: str, pixel: dict, stem: dict) -> tuple[pd.DataFrame, dict]:
    run, cell_source, distance_source = cells_with_distance(donor)

    # QuPath's own answer is renamed, never dropped: the diff below is the whole point, and a
    # column called `islet_num` that silently means one thing here and another downstream is how
    # the original defect survived this long.
    run = run.rename(columns={"islet_num": "islet_num_parent"})
    geo = membership(donor, run["X_centroid"].to_numpy(), run["Y_centroid"].to_numpy(),
                     pixel[donor], stem[donor])
    frame = pd.concat([run.reset_index(drop=True), geo], axis=1)

    # --- the cross-check, reported every run ------------------------------------------------
    # Two routes from the same outlines: containment here, and the signed distance computed by
    # a different program. They must agree, and where they do not the run says by how much
    # rather than silently preferring one.
    by_distance = frame["signed_dist_islet_um"].to_numpy() <= 0
    agree = int((by_distance == frame["in_islet"]).sum())
    disagree = frame.loc[by_distance != frame["in_islet"]]

    # --- what changes against QuPath's Parent -----------------------------------------------
    parent_islet = run["islet_num_parent"].fillna("").astype(str).str.strip()
    was_in_islet = (parent_islet != "") & run["cell_region"].eq("core")
    eligible = frame["qc_analysis_eligible"]
    gained = int((frame["in_islet"] & ~was_in_islet & eligible).sum())
    lost = int((~frame["in_islet"] & was_in_islet & eligible).sum())
    renumbered = int((frame["in_islet"] & was_in_islet & eligible
                      & (frame["islet_num"].astype(str) != parent_islet)).sum())

    report = {
        "donor_id": donor,
        "cell_source": cell_source,
        "distance_source": distance_source,
        "cells": len(frame),
        "eligible": int(eligible.sum()),
        "core_outlines": geo.attrs["n_core_outlines"],
        # BOTH counts carry the same eligibility filter. They did not in the first version, and
        # the resulting table showed geometry finding FEWER islet cells than Parent on 12 of 13
        # donors -- an artefact of comparing an eligible count against a raw one.
        "in_islet_geometry": int((frame["in_islet"] & eligible).sum()),
        # Islets that actually hold an eligible cell -- not the outline count. 146 outlines in
        # this cohort contain no cell at all, so the two are different numbers and the app's
        # denominator is this one.
        "islets_occupied": int(frame.loc[frame["in_islet"] & eligible, "islet_num"].nunique()),
        "in_islet_parent": int((was_in_islet & eligible).sum()),
        "gained": gained,
        "lost": lost,
        "renumbered": renumbered,
        "distance_agrees": agree,
        "distance_disagrees": len(disagree),
        "disagree_median_um": float(disagree["signed_dist_islet_um"].abs().median())
        if len(disagree) else 0.0,
        "overlapping_core": geo.attrs["overlapping_core"],
    }
    keep = ["object_id", "in_islet", "endocrine_object_class", "endocrine_object_id", "qupath_object_name", "islet_num",
            "islet_num_parent", "signed_dist_islet_um"]
    return frame[keep], report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--report", action="store_true",
                        help="report the diff without writing the membership tables")
    args = parser.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    pixel, stem = _manifest()
    if not args.report:
        OUT.mkdir(parents=True, exist_ok=True)

    reports, refused = [], []
    for donor in donors:
        # A donor that cannot be built is recorded with the reason it could not, and the run goes
        # on to the rest. Coverage is the question this lane is asked; stopping at the first
        # failure answers it for one donor and leaves the others unmeasured. The exit status
        # still carries the failure, so nothing passes quietly.
        try:
            frame, report = build_donor(donor, pixel, stem)
        except Exception as error:                                   # noqa: BLE001
            refused.append({"donor_id": donor,
                            "reason": f"{type(error).__name__}: {error}"})
            print(f"{donor}: NOT BUILT -- {type(error).__name__}: {error}", flush=True)
            continue
        reports.append(report)
        print(f"{donor}: {report['in_islet_geometry']:>7,} islet cells by geometry vs "
              f"{report['in_islet_parent']:>7,} by Parent  "
              f"(+{report['gained']:,} / -{report['lost']:,} / renumbered {report['renumbered']:,})"
              f"   signed distance disagrees on {report['distance_disagrees']:,}"
              f" (median |d| {report['disagree_median_um']:.2f} um)", flush=True)
        if not args.report:
            target = OUT / f"donor_id={donor}"
            target.mkdir(parents=True, exist_ok=True)
            frame.to_parquet(target / "data.parquet", index=False)

    table = pd.DataFrame(reports)
    pd.set_option("display.width", 260)
    print("\n" + table.to_string(index=False))
    print(f"\ncovered {len(reports)} of {len(donors)} donors asked for")
    for row in refused:
        print(f"  not built: {row['donor_id']}  {row['reason']}")
    if not args.report:
        table.to_parquet(OUT / "membership_diff.parquet", index=False)
        if refused:
            pd.DataFrame(refused).to_json(OUT / "not_built.json",
                                          orient="records", indent=2)
        print(f"\n-> {OUT}")
    return 1 if refused else 0


if __name__ == "__main__":
    raise SystemExit(main())
