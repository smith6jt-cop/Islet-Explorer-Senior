#!/usr/bin/env python3
"""The viewer's region layer, built from the sources the numbers are measured from.

``data/development/typing/app/json/<donor>.geojson`` is what the Shiny drill-down reads through
``ieCore::domain_data_root()`` and what ``build_tissuumaps_project.py --regions`` copies into
``/data/tissuumaps/regions/``.  Three classes are drawn from it -- ``Islet``, ``IsletExpanded``
and ``IsletNonEndocrine`` -- and until now two of the three were not the thing they are named
after.  Measured 2026-09-17, on the shipped files:

  IsletNonEndocrine  was ``(core U 20 um expansion) - endocrine union`` -- the PRE-correction-E
                     region.  Against the corrected definition its area ran 15.67x too large on
                     6539 and 8.00x on 6436; it was larger than the whole drawn IsletExpanded
                     annotation on 67.9% / 99.7% of islets and larger than the islet core itself
                     on 97.7% / 92.8%.  ``build_islet_geometry.py`` fixed the definition on
                     2026-09-14 (``core.difference(endocrine)``) and rebuilt on 2026-09-17; the
                     shipped json was written 2026-09-14 11:11, between the two.

  IsletExpanded      was the drawn QuPath ``_exp20um`` annotation, which is not one shape.  On
                     6539 it is a filled dilated islet (contains the core on 57 of 60 sampled
                     islets); on 6436 and 6482 it is an annulus with the core punched out (58/60
                     and 60/60 carry an interior ring).  Its outward offset also varies: the 95th
                     percentile distance from the drawn boundary to the core is 19.9-20.5 um on
                     19 donors and 27.7 um on 6482.  The rule is that every expansion is 20 um and
                     comes from the SIGNED DISTANCE, so it is recomputed here as
                     ``core.buffer(20 um).difference(core)`` -- the same geometry
                     ``build_islet_geometry.py`` measures as the ``expansion`` region, so the
                     drawn region and its area are the same object.

  Islet              was already right: the QuPath outline, area identical to the ``core`` row of
                     ``islet_geometry.parquet`` to 1.0000 on every islet of both donors checked.
                     It is carried through unchanged, but resolved with
                     ``build_islet_geometry.islet_export`` -- recursive, newest-export-wins, and
                     refusing a detection export -- rather than a flat glob, so the outlines are
                     the ones the areas were measured from.

Coordinates stay in FULL-RESOLUTION PIXELS, the same frame the QuPath export uses and the frame
the app converts with ``PIXEL_SIZE_UM``.  Nothing here rescales: measured on the shipped files,
json Islet area x (0.5077663810243286)^2 equals ``islet_geometry`` ``core.region_um2`` at ratio
1.0000 for 100% of islets on 6539 and 6436, and the region bbox covers 94%/86% and 94%/93% of the
image at x1 -- x2 would put it at 189%/171% and 189%/186%, outside the image.

Ring orientation is set exterior-CCW / interior-CW, because TissUUmaps fills regions with the
nonzero winding rule (``glUtils._regionFillRule``, default ``"never"`` = outline only) and a hole
wound the same way as its exterior would fill in.

    python scripts/senior/build_tissuumaps_regions.py                 # every eligible donor
    python scripts/senior/build_tissuumaps_regions.py --donors 6539 6436
    python scripts/senior/build_tissuumaps_regions.py --check-only    # measure, write nothing
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from shapely.geometry import mapping, shape
from shapely.geometry.polygon import orient

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))

import build_islet_geometry as GEO  # noqa: E402
import cohort as COHORT  # noqa: E402
import subpanel_io as S  # noqa: E402

#: The expansion width. One number, the same one ``build_islet_geometry.py`` buffers by.
EXPANSION_UM = GEO.EXPANSION_UM

#: Classes this builder derives. Anything else in the QuPath export is carried through untouched
#: and counted, so a class added upstream shows up in the log instead of vanishing.
DERIVED = ("Islet", "IsletExpanded", "IsletNonEndocrine")

#: Same colours the app's drill-down uses (``segmentation_helpers.R``), so a region opened in
#: TissUUmaps and the same region in the app are not two different colours.
CLASS_COLOR = {
    "Islet": [0, 102, 204],
    "IsletExpanded": [0, 204, 204],
    "IsletNonEndocrine": [140, 107, 177],
}

#: Fraction by which a derived region's area may differ from the area
#: ``islet_geometry.parquet`` records for the same region before the build refuses to write.
#: They are the same computation on the same inputs, so this is a floating-point tolerance and
#: not a fit; anything larger means the two have diverged and the layer would misdescribe itself.
AREA_TOLERANCE = 1e-6


def feature(geometry, name: str, class_name: str, source_id: str | None = None) -> dict:
    """One QuPath-shaped annotation feature, oriented for the nonzero fill rule.

    ``source_id`` is the exported annotation's own QuPath UUID where there is one. It is carried
    through rather than regenerated, so a region in the viewer can still be traced back to the
    object in the project it came from.
    """
    if geometry.geom_type == "Polygon":
        geometry = orient(geometry, sign=1.0)
    elif geometry.geom_type == "MultiPolygon":
        geometry = type(geometry)([orient(p, sign=1.0) for p in geometry.geoms])
    feat = {
        "type": "Feature",
        "geometry": mapping(geometry),
        "properties": {
            "name": name,
            "objectType": "annotation",
            "classification": {"name": class_name, "color": CLASS_COLOR[class_name]},
        },
    }
    if source_id is not None:
        feat["id"] = source_id
    return feat


def nonendocrine_features(donor: str) -> tuple[list, dict]:
    """Correction E's outlines for one donor, and islet_key -> geometry beside them.

    Read from ``data/development/typing/geometry/nonendocrine_<donor>.geojson``, which
    ``build_islet_geometry.py`` writes as ``core.difference(endocrine union)``. Absent, the
    donor gets no non-endocrine class at all -- that region is then Unavailable, which is what it
    was before correction E had a source, and is not the same as drawing something else.
    """
    path = S.OUT / f"geometry/nonendocrine_{donor}.geojson"
    if not path.exists():
        return [], {}
    feats = json.loads(path.read_text())["features"]
    geoms = {}
    for feat in feats:
        name = (feat.get("properties") or {}).get("name", "")
        if feat.get("geometry"):
            geoms[name.removesuffix("_nonendo")] = shape(feat["geometry"])
    return feats, geoms


def measured_areas(donor: str) -> dict[str, pd.Series]:
    """islet_geometry's own area for each region, in um^2, keyed by islet."""
    frame = pd.read_parquet(S.OUT / "app/islet_geometry.parquet")
    frame = frame[frame["donor_id"] == donor]
    return {region: sub.set_index("islet_key")["region_um2"]
            for region, sub in frame.groupby("region")}


def _agreement(derived_px2: dict, measured_um2: pd.Series, px2: float) -> tuple[int, float]:
    """How many islets agree with the measured area, and the worst relative difference."""
    n, worst = 0, 0.0
    for key, area in derived_px2.items():
        want = measured_um2.get(key)
        if want is None or not np.isfinite(want) or want <= 0:
            continue
        rel = abs(area * px2 - want) / want
        worst = max(worst, rel)
        n += rel <= AREA_TOLERANCE
    return n, worst


def build_donor(donor: str, pixel: dict, stem: dict) -> tuple[list, dict]:
    scale = pixel[donor]
    px2 = scale * scale
    source = GEO.islet_export(stem[donor])
    export = json.loads(source.read_text())["features"]

    cores: dict[str, object] = {}
    core_ids: dict[str, str] = {}
    passthrough: list = []
    other: dict[str, int] = {}
    drawn_expansions = 0
    for feat in export:
        props = feat.get("properties") or {}
        label = (props.get("classification") or {}).get("name")
        if label == "Islet" and feat.get("geometry"):
            geometry = shape(feat["geometry"])
            cores[props["name"]] = geometry if geometry.is_valid else geometry.buffer(0)
            if feat.get("id") is not None:
                core_ids[props["name"]] = feat["id"]
        elif label == "IsletExpanded":
            drawn_expansions += 1          # replaced by the signed-distance expansion below
        elif label in DERIVED:
            continue                       # never carried from the export; derived here
        else:
            passthrough.append(feat)
            other[label] = other.get(label, 0) + 1

    ne_feats, ne_geoms = nonendocrine_features(donor)

    features: list = []
    core_px2, exp_px2 = {}, {}
    empty_expansions = 0
    for islet, core in cores.items():
        expansion = core.buffer(EXPANSION_UM / scale).difference(core)
        if not expansion.is_valid:
            expansion = expansion.buffer(0)
        features.append(feature(core, islet, "Islet", core_ids.get(islet)))
        core_px2[islet] = core.area
        if expansion.is_empty:
            empty_expansions += 1
            continue
        features.append(feature(expansion, f"{islet}_exp20um", "IsletExpanded"))
        exp_px2[islet] = expansion.area
    features.extend(ne_feats)
    features.extend(passthrough)

    measured = measured_areas(donor)
    ne_px2 = {k: g.area for k, g in ne_geoms.items()}
    checks = {
        "Islet": _agreement(core_px2, measured.get("core", pd.Series(dtype=float)), px2),
        "IsletExpanded": _agreement(exp_px2, measured.get("expansion", pd.Series(dtype=float)), px2),
        "IsletNonEndocrine": _agreement(
            ne_px2, measured.get("non_endocrine", pd.Series(dtype=float)), px2),
    }
    report = {
        "donor": donor, "export": source.name,
        "Islet": len(core_px2), "IsletExpanded": len(exp_px2),
        "IsletNonEndocrine": len(ne_feats),
        "drawn_expansions_replaced": drawn_expansions,
        "empty_expansions": empty_expansions,
        "other_classes": other,
        "checks": checks,
    }
    return features, report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--donors", nargs="*", default=None)
    # The default is CLEAN + MARGINAL, which is the 20 donors islet_geometry was rebuilt over.
    # `--clean-only` drops the three marginal donors; there is no flag that adds them, because a
    # `--include-marginal` that defaulted to true would be a switch that does nothing.
    parser.add_argument("--clean-only", action="store_true",
                        help="build only cohort.CLEAN, leaving out the marginal donors")
    parser.add_argument("--check-only", action="store_true",
                        help="measure and report, write nothing")
    parser.add_argument("--out-dir", type=pathlib.Path, default=S.OUT / "app/json")
    args = parser.parse_args()

    pixel, stem = GEO._manifest()
    eligible = set(COHORT.resolve(None, include_marginal=not args.clean_only))
    have = sorted(p.name[len("endocrine_union_"):-len(".parquet")]
                  for p in (S.OUT / "geometry").glob("endocrine_union_*.parquet"))
    donors = args.donors or [d for d in have if d in eligible]
    refused = sorted(set(have) - eligible)
    if refused:
        print(f"not built (not in the eligible cohort): {refused}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows, failed = [], []
    for donor in donors:
        features, report = build_donor(donor, pixel, stem)
        disagree = []
        for cls, (agree, worst) in report["checks"].items():
            expected = report[cls]
            if expected and agree < expected and worst > AREA_TOLERANCE:
                disagree.append(f"{donor}/{cls}: {expected - agree} of {expected} islets differ "
                                f"from islet_geometry.parquet, worst {worst:.3e}")
        failed.extend(disagree)
        target = args.out_dir / f"{donor}.geojson"
        # A layer that disagrees with the areas it is drawn beside is not shipped: the previous
        # generation of this file was wrong in exactly that way and nothing refused it.
        if not args.check_only and not disagree:
            tmp = target.with_suffix(".geojson.part")
            tmp.write_text(json.dumps({"type": "FeatureCollection", "features": features}))
            tmp.replace(target)
        size = target.stat().st_size / 1024 ** 2 if target.exists() else float("nan")
        print(f"  {donor:>5}  Islet {report['Islet']:>5,}  IsletExpanded "
              f"{report['IsletExpanded']:>5,}  IsletNonEndocrine "
              f"{report['IsletNonEndocrine']:>5,}   {size:7.1f} MB   "
              f"agree core/exp/nonendo "
              f"{report['checks']['Islet'][0]:,}/{report['checks']['IsletExpanded'][0]:,}/"
              f"{report['checks']['IsletNonEndocrine'][0]:,}"
              + (f"   other {report['other_classes']}" if report["other_classes"] else ""),
              flush=True)
        rows.append(report)

    print(f"\n{len(rows)} donors, "
          f"{sum(r['Islet'] for r in rows):,} Islet, "
          f"{sum(r['IsletExpanded'] for r in rows):,} IsletExpanded, "
          f"{sum(r['IsletNonEndocrine'] for r in rows):,} IsletNonEndocrine")
    if failed:
        print("\nREFUSED -- a derived region disagrees with the area islet_geometry records "
              "for it:")
        for line in failed:
            print(f"  {line}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
