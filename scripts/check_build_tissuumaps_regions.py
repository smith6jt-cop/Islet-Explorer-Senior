"""Unit tests for the viewer's region layer: how it is found, and how it is built.

Two modules. ``build_tissuumaps_project.py`` finds a donor's GeoJSON in a directory of exports,
reads it whether or not it is gzipped, and settles the coordinate scale into full-resolution
pixels. ``build_tissuumaps_regions.py`` decides WHAT each region is -- the join that names a
non-endocrine outline, the signed-distance expansion that replaces the drawn annotation, the ring
orientation the nonzero fill rule needs, and the area agreement that refuses a layer disagreeing
with the numbers it is drawn beside. No slide, no server, no pyvips -- both builders keep their
heavy imports function-local, so this runs on a bare interpreter in well under a second.

    python scripts/check_build_tissuumaps_regions.py
"""

from __future__ import annotations

import gzip
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "senior"))
import build_tissuumaps_project as b  # noqa: E402

_tmp = tempfile.TemporaryDirectory()
TMP = Path(_tmp.name)

fails = []


def check(label, got, want):
    ok = got == want
    print(f"{'ok  ' if ok else 'FAIL'}  {label}: {got!r}" + ("" if ok else f" != {want!r}"))
    if not ok:
        fails.append(label)


def fc(*polys):
    return {"type": "FeatureCollection", "features": [
        {"type": "Feature", "id": f"Islet_{i}",
         "properties": {"classification": {"name": "Islet"}},
         "geometry": {"type": "Polygon", "coordinates": [p]}}
        for i, p in enumerate(polys, 1)]}


def box(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]


IMG_W, IMG_H = 40000, 30000

# --- resolve_regions ------------------------------------------------------
(TMP / "d").mkdir(exist_ok=True)
plain = TMP / "d" / "results__6539_IO60panc2nd_resolution__1.geojson"
plain.write_text(json.dumps(fc(box(0, 0, 10, 10))))
gz = TMP / "d" / "results__115_IO60panc2nd_resolution__1.geojson.gz"
with gzip.open(gz, "wt") as fh:
    json.dump(fc(box(0, 0, 10, 10)), fh)

check("directory + case id", b.resolve_regions(TMP / "d", "6539"), plain)
check("directory + gz", b.resolve_regions(TMP / "d", "115"), gz)
check("explicit file", b.resolve_regions(plain, "6539"), plain)
check("no match", b.resolve_regions(TMP / "d", "9999"), None)
check("missing path", b.resolve_regions(TMP / "nope", "6539"), None)

# A case id that is a prefix of another donor's must not match it: a plain
# '*115*' substring glob picks up 1150 and silently attaches the wrong islets.
(TMP / "prefix").mkdir(exist_ok=True)
p1150 = TMP / "prefix" / "results__1150_resolution__1.geojson"
p1150.write_text(json.dumps(fc(box(0, 0, 10, 10))))
check("prefix id does not match longer id", b.resolve_regions(TMP / "prefix", "115"), None)
check("the longer id still matches itself",
      b.resolve_regions(TMP / "prefix", "1150"), p1150)
p115 = TMP / "prefix" / "results__115_resolution__1.geojson"
p115.write_text(json.dumps(fc(box(0, 0, 10, 10))))
check("both present, each resolves", b.resolve_regions(TMP / "prefix", "115"), p115)

# Genuine ambiguity fails closed rather than picking the first sorted hit:
# overlaying the wrong donor's outlines is worse than showing none.
(TMP / "ambig").mkdir(exist_ok=True)
for nm in ("a__6539_v1.geojson", "b__6539_v2.geojson"):
    (TMP / "ambig" / nm).write_text(json.dumps(fc(box(0, 0, 10, 10))))
check("ambiguous match attaches nothing", b.resolve_regions(TMP / "ambig", "6539"), None)

# .json / .json.gz are found too, and .geojson is not double-counted as .json.
(TMP / "ext").mkdir(exist_ok=True)
pj = TMP / "ext" / "islets_6539.json"
pj.write_text(json.dumps(fc(box(0, 0, 10, 10))))
check("plain .json found", b.resolve_regions(TMP / "ext", "6539"), pj)

# --- _load_geojson --------------------------------------------------------
check("gzip load", len(b.region_features(b._load_geojson(gz))), 1)
check("plain load", len(b.region_features(b._load_geojson(plain))), 1)

# --- infer_region_scale ---------------------------------------------------
NOMINAL = b.DEFAULT_REGION_UM_PER_PX / b.DEFAULT_UM_PER_PX   # 2.0
check("nominal is 2x", round(NOMINAL, 6), 2.0)

# Exported at resolution #1: tissue covers 70% of the half-res image.
half = fc(box(1000, 1000, int(0.70 * IMG_W / 2), int(0.70 * IMG_H / 2)))
check("res#1 export -> x2",
      b.infer_region_scale(b.region_features(half), IMG_W, IMG_H, NOMINAL), 2.0)

# Already full-res, tissue covers 70%: x2 would overflow -> caught.
full = fc(box(1000, 1000, int(0.70 * IMG_W), int(0.70 * IMG_H)))
check("full-res export -> x1 (falsified)",
      b.infer_region_scale(b.region_features(full), IMG_W, IMG_H, NOMINAL), 1.0)

# Full-res but tiny section: undetectable, documented limitation.
tiny = fc(box(0, 0, int(0.20 * IMG_W), int(0.20 * IMG_H)))
check("tiny section -> nominal (undetectable)",
      b.infer_region_scale(b.region_features(tiny), IMG_W, IMG_H, NOMINAL), 2.0)

check("no geometry -> nominal", b.infer_region_scale([], IMG_W, IMG_H, NOMINAL), 2.0)
check("zero image -> nominal",
      b.infer_region_scale(b.region_features(half), 0, 0, NOMINAL), 2.0)

# 4x-downsampled export against a nominal that overflows twice.
q = fc(box(0, 0, int(0.70 * IMG_W / 4), int(0.70 * IMG_H / 4)))
check("x8 nominal on quarter-res -> x4",
      b.infer_region_scale(b.region_features(q), IMG_W, IMG_H, 8.0), 4.0)

# --- copy_regions ---------------------------------------------------------
gj = b._load_geojson(plain)
out = TMP / "out.geojson"
n = b.copy_regions(gj, out, 2.0, layer_index=0)
written = json.loads(out.read_text())
f0 = written["features"][0]
check("feature count", n, 1)
check("coords doubled", f0["geometry"]["coordinates"][0][2], [20, 20])
check("collectionIndex stamped", f0["properties"]["collectionIndex"], 0)
check("name from id", f0["properties"]["name"], "Islet_1")
check("classification preserved", f0["properties"]["classification"]["name"], "Islet")

# scale 1.0 must not touch coordinates
n = b.copy_regions(b._load_geojson(plain), out, 1.0)
check("scale 1 leaves coords",
      json.loads(out.read_text())["features"][0]["geometry"]["coordinates"][0][2], [10, 10])

# bare list of features (single-object QuPath export)
bare = TMP / "bare.geojson"
bare.write_text(json.dumps(fc(box(0, 0, 5, 5))["features"]))
check("bare list", b.copy_regions(b._load_geojson(bare), out, 2.0), 1)

# --- build_tissuumaps_regions: what each region IS -------------------------
import math  # noqa: E402

import build_tissuumaps_regions as r  # noqa: E402
from shapely.geometry import Polygon, shape  # noqa: E402

# The join that the previous comparison failed on. The non-endocrine outline is named
# `<islet_key>_nonendo` and the expansion `<islet_key>_exp20um`; the islet_key is what is left
# when that suffix comes off, and it is what endocrine_union.islet_key carries. Stripping a fixed
# character count instead -- `name[:-9]` for an 8-character suffix -- silently joins Islet_17 to
# Islet_1 and produces a comparison of the wrong islets that still looks populated.
check("nonendo suffix strips to the islet key", "Islet_17_nonendo".removesuffix("_nonendo"), "Islet_17")
check("exp20um suffix strips to the islet key", "Islet_17_exp20um".removesuffix("_exp20um"), "Islet_17")
check("a bare islet key is unchanged", "Islet_17".removesuffix("_nonendo"), "Islet_17")

# Ring orientation: TissUUmaps fills with the nonzero winding rule, so an interior ring wound the
# same way as its exterior fills in rather than punching a hole.
square = Polygon(box(0, 0, 100, 100), [box(60, 60, 40, 40)])
oriented = shape(r.feature(square, "Islet_1_exp20um", "IsletExpanded")["geometry"])
check("exterior is counter-clockwise", oriented.exterior.is_ccw, True)
check("interior is clockwise", oriented.interiors[0].is_ccw, False)
check("orientation leaves the area alone", round(oriented.area, 6), round(square.area, 6))

f = r.feature(square, "Islet_1_exp20um", "IsletExpanded")
check("feature carries the app's name", f["properties"]["name"], "Islet_1_exp20um")
check("feature carries the class", f["properties"]["classification"]["name"], "IsletExpanded")
check("feature is an annotation", f["properties"]["objectType"], "annotation")

# The expansion is derived, not read: a 20 um band outside the islet outline, computed with that
# image's pixel size, and a RING -- the core is not part of it, because the core is its own region
# and islet_geometry measures the two separately.
PX = 0.5077663810243286
core = Polygon([(p[0] * 1.0, p[1] * 1.0) for p in box(0, 0, 200, 200)])
band = core.buffer(r.EXPANSION_UM / PX).difference(core)
check("the expansion excludes the core", band.intersection(core).area < 1e-9, True)
check("the expansion is 20 um wide",
      round(band.bounds[2] - core.bounds[2], 3), round(r.EXPANSION_UM / PX, 3))
check("20 um is the one number", r.EXPANSION_UM, 20.0)

# Area agreement is the refusal: a derived region is compared against the um^2 islet_geometry
# records for the same region, and a layer that disagrees is not written.
import pandas as _pd  # noqa: E402

derived = {"Islet_1": 1000.0, "Islet_2": 2000.0}
measured = _pd.Series({"Islet_1": 1000.0 * PX * PX, "Islet_2": 2000.0 * PX * PX})
n, worst = r._agreement(derived, measured, PX * PX)
check("matching areas agree", n, 2)
check("matching areas have no drift", worst < r.AREA_TOLERANCE, True)

skewed = _pd.Series({"Islet_1": 1000.0 * PX * PX, "Islet_2": 2400.0 * PX * PX})
n, worst = r._agreement(derived, skewed, PX * PX)
check("a 20% disagreement is caught", n, 1)
check("and is reported at its size", round(worst, 4), 0.1667)

# A region with no measured area is UNAVAILABLE, not a pass and not a failure: it is skipped.
absent = _pd.Series({"Islet_1": 1000.0 * PX * PX, "Islet_2": float("nan")})
n, worst = r._agreement(derived, absent, PX * PX)
check("an unavailable area is skipped, not counted", (n, worst < r.AREA_TOLERANCE), (1, True))
check("nothing measured means nothing to disagree with",
      r._agreement(derived, _pd.Series(dtype=float), PX * PX), (0, 0.0))

check("the classes this builder derives", r.DERIVED,
      ("Islet", "IsletExpanded", "IsletNonEndocrine"))
check("every derived class has a colour",
      sorted(r.CLASS_COLOR), sorted(r.DERIVED))

print()
print(f"{len(fails)} failed" if fails else "all checks passed")
sys.exit(1 if fails else 0)
