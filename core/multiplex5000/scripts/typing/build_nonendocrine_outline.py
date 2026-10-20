#!/usr/bin/env python3
"""The non-endocrine region's outline -- correction E.

Within an islet and its 20 um expansion, the part that is **not** endocrine cells: an outline the
app draws beside the islet and the expansion, and a region whose structure areas are measurable by
intersection like every other region's.  Until now it had neither.  ``build_islet_geometry.py`` set
``geometry["non_endocrine"] = None``, so all 23,346 ``_nonendo`` rows of ``Islet_Targets`` carried a
NaN structure area -- correctly Unavailable rather than falsely zero, but there was no shape to
intersect because none had been built.

    non_endocrine = (islet core U expansion) - union(the endocrine cells' own segmented outlines)

**The only per-cell polygons are the QuPath cell-outline exports the manifest declares**
(``cell_geojson.path``; read through ``subpanel_io.cell_outlines``).  For the v2 cohort those are
the 61 GB under ``data/redsea_scratch/geojson/``; for the v3 exports they are
``results/v3/export/cell_geojson/``.  No PhenoCycler stage output carries any -- neither the
current flat roots nor the historical ``data/runs/`` tree, whose ``cells``, ``expression`` and
``geometry_qc`` stages have 129, 69 and 15 columns and not one of them is geometry.  Those files
are keyed by the QuPath UUID, so they join to ``object_id`` directly.

The files are single-line FeatureCollections of several GB, so they are streamed rather than read.

    python core/multiplex5000/scripts/typing/build_nonendocrine_outline.py [--donors 6539]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import ijson
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from shapely.geometry import Polygon, shape
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cohort as COHORT
import subpanel_io as S

OUT = S.OUT / "geometry"


def endocrine_ids(donor: str) -> tuple[set[str], pd.DataFrame]:
    """The islet-associated cells of one donor, and which of them are endocrine."""
    cells = pq.read_table(
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        columns=["object_id", "cell_class", "cell_region_20um", "islet_key", "is_endocrine"],
    ).to_pandas()
    cells = cells[(cells["cell_class"] != "Excluded")
                  & cells["cell_region_20um"].isin(["core", "expansion"])
                  & cells["islet_key"].notna()]
    cells = cells.merge(S.cell_areas(donor), on="object_id", how="left", validate="one_to_one")
    return set(cells.loc[cells["is_endocrine"], "object_id"]), cells


def read_polygons(path, wanted: set[str]) -> dict[str, Polygon]:
    """Stream one multi-GB FeatureCollection, keeping only the wanted cells' outlines."""
    found: dict[str, Polygon] = {}
    started = time.time()
    with open(path, "rb") as handle:
        for feature in ijson.items(handle, "features.item", use_float=True):
            key = feature.get("id")
            if key in wanted:
                found[key] = shape(feature["geometry"])
                if len(found) == len(wanted):
                    break
    print(f"    {len(found):,} / {len(wanted):,} outlines in {time.time() - started:.0f} s",
          flush=True)
    return found


def is_fresh(donor: str, cell_geojson: Path) -> bool:
    """True when this donor's outline already postdates every input it is built from.

    The output is a function of exactly one thing that can move: WHICH cells are endocrine.
    The islet outlines are QuPath's, the cell polygons are the segmentation's, and neither is
    rewritten by this lane. So the inputs that matter are that donor's typing partition, its
    per-cell outline file, the QuPath manifest that supplies the pixel size, and this script.

    Why this is worth a guard rather than a comment: `read_polygons` streams a single-line
    FeatureCollection that reaches 5.5 GB, one feature at a time through `ijson`, and the loop
    did it unconditionally -- about a minute per donor, sixteen for the cohort, to re-derive a
    file that was already correct. The cost is JSON parsing, not geometry.

    mtimes, not hashes, to match `scripts/senior/build_layers.py --check`, which is the gate
    this lane is scheduled by. A no-op rewrite of an input therefore still forces a rebuild --
    that is the known limitation of the whole staleness scheme here, not something introduced
    at this call site.
    """
    out = OUT / f"endocrine_union_{donor}.parquet"
    if not out.exists():
        return False
    inputs = [
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        S.cell_area_source(donor),
        cell_geojson,
        S.REPO / "core/phenocycler/data/qupath_manifest.json",
        Path(__file__).resolve(),
    ]
    present = [i for i in inputs if i.exists()]
    if len(present) != len(inputs):
        # An input that is not there cannot be compared against, so nothing is claimed about
        # freshness: rebuild and let the build fail on the real missing file.
        return False
    return out.stat().st_mtime >= max(i.stat().st_mtime for i in present)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="rebuild a donor even when its outline is already newer than "
                             "every input it reads")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    # The scratch outlines are in PIXELS, like the islet outlines -- convert with each image's own
    # pixel size, exactly as build_islet_geometry.py does. An area left in px^2 looks plausible and
    # is wrong by a factor of ~4 here.
    pixel = {r["donor_id"]: r["pixel_size_um_x"]
             for r in json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json")
                                 .read_text())["images"]}
    # No islet geometry is read here. The flow is one-way: this script writes
    # endocrine_union_<donor>.parquet and build_islet_geometry.py reads it. A dead
    # `pd.read_parquet(S.OUT / "app/islet_geometry.parquet")` stood on this line until
    # 2026-09-18 -- assigned, never used -- making the script fail on a downstream output it
    # does not need and reading as though the outlines were joined to islet geometry.
    OUT.mkdir(parents=True, exist_ok=True)

    skipped = []
    for donor in donors:
        # the outline file the manifest declares for this donor, not a folder globbed for a name
        matches = [S.cell_outlines(donor)]
        if not matches[0].is_file():
            raise FileNotFoundError(f"{donor}: cell-outline file {matches[0]} is not on disk")

        if not args.force and is_fresh(donor, matches[0]):
            print(f"{donor}: fresh -- skipped (--force to rebuild)", flush=True)
            skipped.append(donor)
            continue

        print(f"{donor}:", flush=True)
        wanted, cells = endocrine_ids(donor)
        outlines = read_polygons(matches[0], wanted)
        missing = len(wanted) - len(outlines)
        if missing:
            print(f"    {missing:,} endocrine cells have no outline in the scratch export",
                  flush=True)

        area_um2 = pixel[donor] ** 2
        by_islet = cells[cells["is_endocrine"]].groupby("islet_key")["object_id"]
        records = []
        for islet, ids in by_islet:
            shapes = [outlines[i] for i in ids if i in outlines]
            if not shapes:
                continue
            endocrine = unary_union(shapes)
            records.append({"donor_id": donor, "islet_key": islet,
                            "endocrine_um2_from_outlines": float(endocrine.area) * area_um2,
                            "n_endocrine_outlines": len(shapes),
                            "wkb": endocrine.wkb})
        frame = pd.DataFrame(records)
        frame.to_parquet(OUT / f"endocrine_union_{donor}.parquet", index=False)

        # Cross-check against the SAME cells' QuPath-reported areas -- the endocrine ones only,
        # not the region's whole cell area, or the ratio absorbs the non-endocrine share and stops
        # being a check. The union is smaller wherever neighbouring cells touch, so it should sit
        # at or just below 1; far from it means the wrong polygons or the wrong units.
        reported = (cells[cells["is_endocrine"]].groupby("islet_key")["cell_area"].sum()
                    .rename("endocrine_reported_um2").reset_index())
        check = frame.merge(reported, on="islet_key", how="inner")
        ratio = (check["endocrine_um2_from_outlines"] / check["endocrine_reported_um2"]).replace(
            [np.inf, -np.inf], np.nan).dropna()
        print(f"    {len(frame):,} islets -> endocrine_union_{donor}.parquet  "
              f"(union/reported-area median {ratio.median():.3f})", flush=True)

    if skipped:
        print(f"skipped as fresh ({len(skipped)}): {', '.join(skipped)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
