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

Scope: islet-core and peri-islet cells, the population the study question is about. Extending to
the whole tissue is the same code over ~13.4 M cells instead of 800 k.

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
from shapely import STRtree
from shapely.geometry import shape
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


def islet_cells(donor: str) -> pd.DataFrame:
    cells = pq.read_table(
        S.OUT / f"cells/donor_id={donor}/data.parquet",
        columns=["object_id", "cell_class", "cell_region_20um", "islet_key", "phenotype"],
    ).to_pandas()
    return cells[(cells["cell_class"] != "Excluded")
                 & cells["cell_region_20um"].isin(["core", "expansion"])
                 & cells["islet_key"].notna()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()
    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    pixel, stem = _manifest()
    OUT.mkdir(parents=True, exist_ok=True)

    for donor in donors:
        started = time.time()
        cells = islet_cells(donor)
        wanted = set(cells["object_id"])
        matches = [p for p in CELL_GEOJSON.glob("cells__*.geojson") if stem[donor] in p.name]
        if len(matches) != 1:
            raise FileNotFoundError(f"{donor}: {len(matches)} cell-outline files")

        objects = {name: CX.load_objects(stem[donor], name) for name in STRUCTURES}
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
                    tree = trees[name]
                    if tree is None or cell_area <= 0:
                        row[f"crossed_frac_{column}"] = np.nan
                        continue
                    hit = tree.query(outline)
                    if len(hit) == 0:
                        row[f"crossed_frac_{column}"] = 0.0
                        continue
                    overlap = unary_union([objects[name][i] for i in hit]).intersection(outline)
                    row[f"crossed_frac_{column}"] = float(overlap.area) / cell_area
                rows.append(row)
                if len(rows) == len(wanted):
                    break

        frame = pd.DataFrame(rows)
        frame.to_parquet(OUT / f"crossed_{donor}.parquet", index=False)

        summary = []
        for column in STRUCTURES.values():
            values = frame[f"crossed_frac_{column}"]
            crossed = values > 0
            summary.append(f"{column} {100 * crossed.mean():4.1f}% crossed, "
                           f"median-when-crossed {values[crossed].median():.3f}"
                           if crossed.any() else f"{column} none")
        print(f"  {donor}: {len(frame):>7,} cells in {time.time() - started:4.0f} s  |  "
              + "  |  ".join(summary), flush=True)

    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
