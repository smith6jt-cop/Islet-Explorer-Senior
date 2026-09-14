#!/usr/bin/env python3
"""Per-islet regions and their areas, recomputed so the denominator matches the cells.

Four regions per islet, not two:

  core           inside the islet outline
  expansion      the 20 um band outside it -- **recomputed** as a 20 um outward buffer of the
                 islet's own outline, not read from the annotation.  One donor's annotation is
                 drawn at 23.1 um while every other sits at 18.8-19.3, and all of them are named
                 ``_exp20um``, so reading the annotation would silently give that donor ~21% more
                 area than the cells assigned to it by the signed distance.
  non_endocrine  the part of core + expansion that is not endocrine cells.  Its area is the summed
                 segmented area of the non-endocrine cells in it -- a measurement, not the region
                 area minus an estimate.
  union          core + expansion

Structure areas come from intersecting the classified objects with each region, which is the same
evidence the classification used to decide identity -- measured as area here, as identity there.

Outlines are in pixels and centroids in micrometres, so every join uses that image's pixel size.

    python core/phenocycler/scripts/typing/build_islet_geometry.py [--donors ...]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from shapely.geometry import mapping, shape
from shapely import wkb
from shapely.ops import unary_union
from shapely.strtree import STRtree

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import cohort as COHORT
import context as CX
import subpanel_io as S

EXPANSION_UM = 20.0
STRUCTURES = ("Blood vessel", "Lymphatic", "Nerve")
#: The app's Microenvironment view names the vessel class "Capillary"; keep its vocabulary.
STRUCTURE_DISPLAY = {"Blood vessel": "Capillary", "Lymphatic": "Lymphatic", "Nerve": "Nerve"}


def _manifest() -> tuple[dict, dict]:
    raw = json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
    return ({r["donor_id"]: r["pixel_size_um_x"] for r in raw["images"]},
            {r["donor_id"]: r["image_id"].split(" - ")[0].replace(".qptiff", "")
             for r in raw["images"]})


def islet_outlines(stem: str) -> dict[str, object]:
    """`Islet_N` -> its outline, in pixels, for one image."""
    matches = [p for p in S.ISLETS.glob("*.geojson") if stem in p.name]
    if len(matches) != 1:
        raise FileNotFoundError(f"islet outlines for {stem!r}: found {len(matches)} files")
    out = {}
    for feature in json.loads(matches[0].read_text())["features"]:
        properties = feature.get("properties", {})
        label = (properties.get("classification") or {}).get("name")
        if label != "Islet" or not feature.get("geometry"):
            continue
        geometry = shape(feature["geometry"])
        out[properties["name"]] = geometry if geometry.is_valid else geometry.buffer(0)
    return out


def endocrine_unions(donor: str) -> dict[str, object]:
    """islet_key -> the union of that islet's endocrine cell outlines, in pixels.

    Written by ``build_nonendocrine_outline.py`` from the only per-cell polygons in this checkout.
    Absent, the non-endocrine region has no outline and its structure areas stay Unavailable --
    which is what they were before correction E had a source.
    """
    path = S.OUT / f"geometry/endocrine_union_{donor}.parquet"
    if not path.exists():
        return {}
    frame = pd.read_parquet(path)
    return {row.islet_key: wkb.loads(row.wkb) for row in frame.itertuples()}


def build_donor(donor: str, pixel: dict, stem: dict) -> pd.DataFrame:
    scale = pixel[donor]
    area_um2 = scale * scale                      # one pixel of outline area, in um^2
    outlines = islet_outlines(stem[donor])

    cells = pq.read_table(S.OUT / f"cells/donor_id={donor}/data.parquet",
                          columns=["islet_key", "cell_region_20um", "is_endocrine",
                                   "phenotype", "object_id"]).to_pandas()
    areas = pq.read_table(S.RUN / f"cells/donor_id={donor}/data_0.parquet",
                          columns=["object_id", "cell_area"]).to_pandas()
    cells = cells.merge(areas, on="object_id", how="left", validate="one_to_one")

    endocrine = endocrine_unions(donor)
    objects = {name: CX.load_objects(stem[donor], name) for name in STRUCTURES}
    trees = {name: (STRtree(geoms) if geoms else None) for name, geoms in objects.items()}

    rows, outlines_out = [], []
    for islet, core in outlines.items():
        expansion = core.buffer(EXPANSION_UM / scale).difference(core)
        union = unary_union([core, expansion])
        here = cells[cells["islet_key"] == islet]
        region_cells = {
            "core": here[here["cell_region_20um"] == "core"],
            "expansion": here[here["cell_region_20um"] == "expansion"],
        }
        region_cells["union"] = pd.concat(region_cells.values())
        region_cells["non_endocrine"] = region_cells["union"][~region_cells["union"]["is_endocrine"]]

        # Correction E: the non-endocrine region is core + expansion minus what the endocrine
        # cells actually occupy. With an outline it has a measurable area and its structure areas
        # come from intersection like every other region's; without one (no endocrine union on
        # disk for this donor) it stays None and those areas stay Unavailable.
        non_endocrine = None
        if islet in endocrine:
            non_endocrine = union.difference(endocrine[islet])
            if not non_endocrine.is_valid:
                non_endocrine = non_endocrine.buffer(0)
        geometry = {"core": core, "expansion": expansion, "union": union,
                    "non_endocrine": non_endocrine}
        for region, subset in region_cells.items():
            shape_here = geometry[region]
            # Every region's area is now its outline's; the non-endocrine one falls back to what
            # its cells occupy only where no endocrine union was built for that islet.
            region_area = (float(shape_here.area) * area_um2 if shape_here is not None
                           else float(subset["cell_area"].sum()))
            row = {"donor_id": donor, "islet_key": islet,
                   "combined_islet_id": f"{donor}_{islet}", "region": region,
                   "region_um2": region_area, "n_cells": int(len(subset)),
                   "cells_um2": float(subset["cell_area"].sum())}
            if region == "core":
                row["islet_diam_um"] = float(2.0 * np.sqrt(region_area / np.pi))
                row["centroid_x_um"] = float(core.centroid.x) * scale
                row["centroid_y_um"] = float(core.centroid.y) * scale
            for name in STRUCTURES:
                key = f"area_um2_{STRUCTURE_DISPLAY[name]}"
                if shape_here is None or trees[name] is None:
                    row[key] = np.nan            # no outline to intersect: unavailable, not zero
                    continue
                hit = trees[name].query(shape_here)
                if len(hit) == 0:
                    row[key] = 0.0
                    continue
                overlap = unary_union([objects[name][i] for i in hit]).intersection(shape_here)
                row[key] = float(overlap.area) * area_um2
            rows.append(row)

        # The outline itself, for the app to draw beside the islet and its expansion. Named and
        # classified the way the viewer matches: `<islet_key>_nonendo`, class IsletNonEndocrine.
        if non_endocrine is not None and not non_endocrine.is_empty:
            outlines_out.append({
                "type": "Feature",
                "geometry": mapping(non_endocrine),
                "properties": {"name": f"{islet}_nonendo", "objectType": "annotation",
                               "classification": {"name": "IsletNonEndocrine",
                                                  "color": [140, 107, 177]}},
            })
    return pd.DataFrame(rows), outlines_out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    args = parser.parse_args()

    donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)
    pixel, stem = _manifest()
    frames = []
    outline_dir = S.OUT / "geometry"
    outline_dir.mkdir(parents=True, exist_ok=True)
    for donor in donors:
        frame, nonendo = build_donor(donor, pixel, stem)
        frames.append(frame)
        (outline_dir / f"nonendocrine_{donor}.geojson").write_text(json.dumps(
            {"type": "FeatureCollection", "features": nonendo}))
        core = frame[frame.region == "core"]
        share = frame[frame.region == "non_endocrine"]
        print(f"  {donor}: {len(core):>5,} islets  median core {core.region_um2.median():>9,.0f} um^2"
              f"  median diam {core.islet_diam_um.median():>6.1f} um"
              f"  non-endocrine outlines {len(nonendo):>5,}"
              f"  median non-endo/union {(share.region_um2.median() / frame[frame.region == 'union'].region_um2.median()):.3f}",
              flush=True)

    out = pd.concat(frames, ignore_index=True)
    target = S.OUT / "app/islet_geometry.parquet"
    target.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(target, index=False)
    print(f"\n{out.combined_islet_id.nunique():,} islets x {out.region.nunique()} regions "
          f"-> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
