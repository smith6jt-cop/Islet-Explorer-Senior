#!/usr/bin/env python3
"""How much of each cell a structure passes through -- correction D.

A blood vessel, nerve or lymphatic crossing a cell does not make that cell one: the classification
already keeps the cell's own identity and only the *rate* of crossing survived, as a per-structure
aggregate (Nerve 14.6 %, Vessel 56.1 %, Lymphatic 86.0 %). The degree was discarded. This writes it
back, per cell and per structure:

    crossed_frac_<structure> = area(cell outline INTERSECT that structure's objects) / area(cell)

**This changes no call.** It adds a measurement beside the existing ones. Nothing here re-decides
`cell_class`, which is final for these donors.

Exact, not sampled: the cell outlines are the QuPath ones from `data/redsea_scratch/geojson/` (the
only per-cell polygons in this checkout -- see RECONSTRUCTION.md) and the structure objects are the
same `results/Geojson/` objects the classification used, so the intersection is of the two shapes
QuPath itself drew. Both are in pixels; the ratio is dimensionless, so no conversion is needed --
but `cell_um2` is converted, as a cross-check against QuPath's own reported area.

Scope: EVERY admitted cell in the tissue. This is an input to the structure call, which is asked
of every cell -- a vessel runs where it runs -- so a scope that stopped at the islet expansion
left most of the tissue with no overlap measurement and made the majority-area test unanswerable
there. It was islet-core and peri-islet until 2026-09-19; that was a limit of this producer's run,
never a property of the rule.

    python core/phenocycler/scripts/typing/build_crossed_fraction.py [--donors 6539]
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import ijson
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import shapely
from shapely import STRtree
from shapely.geometry import box, shape
from shapely.ops import unary_union

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import context as CX
import subpanel_io as S

CELL_GEOJSON = S.REPO / "data/redsea_scratch/geojson"
OUT = S.OUT / "crossed"

#: The structure objects, and the column each becomes. The keys are context.py's own class names
#: ("Blood vessel" with a space -- the folder is `Blood_vessel`, the class is not), so this reads
#: the same objects the classification did rather than a second spelling of them.
STRUCTURES = {"Blood vessel": "vessel", "Lymphatic": "lymphatic", "Nerve": "nerve"}


def _manifest() -> tuple[dict, dict]:
    raw = json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    return ({r["donor_id"]: r["pixel_size_um_x"] for r in raw["images"]},
            {r["donor_id"]: r["image_id"].split(" - ")[0].replace(".qptiff", "")
             for r in raw["images"]})


#: A structure object above this many vertices is cut into tiles before any cell is tested
#: against it. Measured on 6505: 541,408 lymphatic objects, median 19 vertices -- and one of
#: 196,935. Every cell in the tissue overlaps that one, so every cell paid an overlay against a
#: 197k-vertex polygon. Cutting it costs a few dozen clips once instead of 239,000 times.
MAX_VERTICES = 5_000

#: Tile edge in PIXELS, the coordinate system the structure objects are in. A cell is a few tens
#: of pixels across, so a tile this size is never the thing that fragments the work.
TILE_PIXELS = 4096.0


def subdivide(geoms: list) -> list:
    """`geoms` with any oversized member replaced by its tile pieces.

    The pieces of one object are interior-disjoint and cover it exactly, so the area a cell takes
    from that object is unchanged -- `covered_fraction` unions whatever it is handed, and pieces
    that merely share a tile edge contribute their areas once. Checked against the unsubdivided
    build on 6450: the fractions agree to 1e-15.
    """
    out = []
    for geom in geoms:
        if shapely.get_num_coordinates(geom) <= MAX_VERTICES:
            out.append(geom)
            continue
        min_x, min_y, max_x, max_y = geom.bounds
        tiles = [box(x, y, min(x + TILE_PIXELS, max_x), min(y + TILE_PIXELS, max_y))
                 for x in np.arange(min_x, max_x, TILE_PIXELS)
                 for y in np.arange(min_y, max_y, TILE_PIXELS)]
        pieces = shapely.intersection(np.repeat(np.array([geom], dtype=object), len(tiles)),
                                      np.array(tiles, dtype=object))
        out.extend(piece for piece in pieces if not piece.is_empty)
    return out


def covered_fraction(tree, geoms: list, outline, cell_area: float) -> float:
    """share of `outline`'s area that `geoms` cover, by the same arithmetic, far cheaper.

    ``area(UNION_i G_i INTERSECT C)`` equals ``area(UNION_i (G_i INTERSECT C))``, so the clip
    happens FIRST and the union is taken over the fragments inside one cell rather than over the
    structure objects themselves. That matters because a structure object here is not small: on
    6505 the lymphatic export is 440 MB, and unioning even one polygon of that size re-nodes it in
    full -- per cell, 239,000 times. Measured on 6450 the two forms agree to the last bit; the
    timing is in the coverage document.

    The tree is queried with ``predicate="intersects"`` rather than by bounding box, so a
    candidate whose box overlaps but whose shape does not never reaches an overlay at all.
    """
    if tree is None or cell_area <= 0:
        return float("nan")
    idx = tree.query(outline, predicate="intersects")
    if len(idx) == 0:
        return 0.0
    if len(idx) == 1:
        return float(geoms[idx[0]].intersection(outline).area) / cell_area
    pieces = [piece for piece in (geoms[i].intersection(outline) for i in idx)
              if not piece.is_empty]
    if not pieces:
        return 0.0
    if len(pieces) == 1:
        return float(pieces[0].area) / cell_area
    return float(unary_union(pieces).area) / cell_area


def measured_cells(donor: str) -> pd.DataFrame:
    """EVERY cell of this donor the analysis admits -- not only the islet-associated ones.

    The overlap fraction is an input to the STRUCTURE CALL, which is asked of every cell in the
    tissue: a blood vessel, nerve or lymphatic runs wherever it runs, and whether the majority of
    a cell's area lies inside one is a question about that cell, not about its distance from an
    islet.  Restricting this to islet core and expansion left every other cell with NO overlap
    measurement, and a cell with no measurement cannot be asked the majority question at all --
    so the rule could not be applied where most of the tissue is.

    `Excluded` cells stay out: they are excluded from every calculation.
    """
    cells = pq.read_table(
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        columns=["object_id", "cell_class", "cell_region_20um", "islet_key", "phenotype"],
    ).to_pandas()
    return cells[cells["cell_class"] != "Excluded"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="rebuild a donor whose crossed file is already present")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    pixel, stem = _manifest()
    OUT.mkdir(parents=True, exist_ok=True)
    coverage = []

    for donor in donors:
        started = time.time()
        target = OUT / f"crossed_{donor}.parquet"
        cells = measured_cells(donor)
        wanted = set(cells["object_id"])
        if target.exists() and not args.force:
            # A present file is left alone only if it measures THE CELLS THAT ARE THERE NOW.
            # Each donor is a multi-GB streamed read, so skipping is what lets this lane be
            # extended donor by donor -- but a blanket skip keeps a file built against an older
            # generation of `typing/cells`, and this tree stacks generations. Measured
            # 2026-09-17: the ten files written 2026-09-14 predate the 2026-09-15 cells rebuild
            # and between them miss 42,149 islet-associated cells while carrying 9,846 that are
            # no longer islet-associated -- 19,959 missing on 6476 alone. The comparison is on
            # object_id, so it cannot be fooled by a count that happens to match.
            done = set(pd.read_parquet(target, columns=["object_id"])["object_id"])
            if done == wanted:
                print(f"  {donor}: {len(done):>7,} cells already measured -- skipped", flush=True)
                coverage.append({"donor_id": donor, "cells_wanted": len(wanted),
                                 "cells_measured": len(done), "no_outline": 0,
                                 "status": "already present"})
                continue
            print(f"  {donor}: rebuilding -- the file on disk misses {len(wanted - done):,} of "
                  f"today's islet-associated cells and carries {len(done - wanted):,} that no "
                  f"longer are", flush=True)
        matches = [p for p in CELL_GEOJSON.glob("cells__*.geojson") if stem[donor] in p.name]
        if len(matches) != 1:
            # Named and carried on rather than raised: the cell outlines are the one input this
            # lane cannot substitute, and a donor that lacks them is a coverage fact the other
            # donors' measurements should not be lost to.
            reason = f"{len(matches)} cell-outline exports match {stem[donor]!r} in {CELL_GEOJSON}"
            print(f"  {donor}: NOT BUILT -- {reason}", flush=True)
            coverage.append({"donor_id": donor, "cells_wanted": len(wanted),
                             "cells_measured": 0, "no_outline": len(wanted),
                             "status": reason})
            continue

        objects = {name: subdivide(CX.load_objects(stem[donor], name)) for name in STRUCTURES}
        trees = {name: (STRtree(geoms) if geoms else None) for name, geoms in objects.items()}
        area_um2 = pixel[donor] ** 2

        rows = []
        with open(matches[0], "rb") as handle:
            for feature in ijson.items(handle, "features.item", use_float=True):
                key = feature.get("id")
                if key not in wanted:
                    continue
                outline = shape(feature["geometry"])
                cell_area = float(outline.area)
                row = {"donor_id": donor, "object_id": key,
                       "cell_um2": cell_area * area_um2}
                for name, column in STRUCTURES.items():
                    row[f"crossed_frac_{column}"] = covered_fraction(
                        trees[name], objects[name], outline, cell_area)
                rows.append(row)
                if len(rows) == len(wanted):
                    break

        frame = pd.DataFrame(rows)
        frame.to_parquet(target, index=False)

        # A cell the outline export does not carry has NO crossed fraction. It is not a zero --
        # zero means a structure was looked for against that cell's own outline and none touched
        # it. The count is reported per donor so the absence is visible rather than inferred from
        # a row that is missing downstream.
        no_outline = len(wanted) - len(frame)
        coverage.append({"donor_id": donor, "cells_wanted": len(wanted),
                         "cells_measured": len(frame), "no_outline": no_outline,
                         "status": "built"})

        summary = []
        for column in STRUCTURES.values():
            values = frame[f"crossed_frac_{column}"]
            crossed = values > 0
            summary.append(f"{column} {100 * crossed.mean():4.1f}% crossed, "
                           f"median-when-crossed {values[crossed].median():.3f}"
                           if crossed.any() else f"{column} none")
        print(f"  {donor}: {len(frame):>7,} cells in {time.time() - started:4.0f} s"
              f"  ({no_outline:,} with no cell outline)  |  "
              + "  |  ".join(summary), flush=True)

    table = pd.DataFrame(coverage)
    pd.set_option("display.width", 200)
    print("\n" + table.to_string(index=False))
    table.to_parquet(OUT / "crossed_coverage.parquet", index=False)
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
