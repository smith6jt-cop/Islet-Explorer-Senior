"""Unit tests for build_tissuumaps_project.py's region handling.

Covers the parts the end-to-end smoke test only exercises down one path: finding
a donor's GeoJSON in a directory of exports, reading it whether or not it is
gzipped, and settling the coordinate scale into full-resolution pixels. No slide,
no server, no pyvips -- the builder keeps its heavy imports function-local, so
this runs on a bare interpreter in well under a second.

    python scripts/test_build_tissuumaps_regions.py
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

print()
print(f"{len(fails)} failed" if fails else "all checks passed")
sys.exit(1 if fails else 0)
