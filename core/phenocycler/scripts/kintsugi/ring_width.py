"""How wide each donor's peri-islet rings actually are, measured from the annotation geometry.

``cell_region == 'peri'`` is read off the QuPath ``Parent`` string, whose annotations are
named ``Islet_N_exp20um``.  The name asserts 20 um.  The geometry does not agree on every
donor, and nothing in ``expression`` can show it, because the categorical label survives any
width.

Measured over ``data/Islets/*.geojson``, which carries the ``Islet`` and ``IsletExpanded``
polygon for each islet: **19 of 20 run donors draw their rings at a median 20.1 um, and donor 6482 draws
its rings at a median 25.1 um.**  So a peri-islet comparison across the cohort puts 6482's ~25 um annulus beside
everyone else's 20 um one -- more area, more cells, and a density that is diluted by a
geometric fact rather than a biological one.

Two consequences, and only the second is a repair:

* any statistic taken over ``cell_region == 'peri'`` inherits this, including the existing
  cosine lane's islet refinement.  Recorded here so the number is not read as biology;
* the continuous signed distance does **not** inherit it.  A zone defined as ``0 < d <= 50``
  um is the same zone in every donor regardless of how the annotation was drawn, which is
  what KINTSUGI's peri-islet analysis uses.  This is the strongest argument for Step 0: it
  does not merely add resolution, it removes a cohort defect that the categorical label
  cannot express.

There is one ring per islet -- 67 of them in donor 6380, 1,304 in donor 115 -- so the
reported width is a MEDIAN OVER RINGS, not a single number the donor possesses.  The spread
matters and is returned beside it: rings whose expansions merge with a neighbour's read a
smaller width, which is why ``p5`` runs to 0 for several donors while ``p95`` sits at the
nominal.  Per ring, the width is estimated as the median distance from the expanded polygon's boundary vertices
to the islet polygon's boundary, over the first ``MAX_ISLETS`` matched pairs.  Median rather
than min or mean because expansions of neighbouring islets merge where islets are close,
which drives some vertices to zero distance -- a real property of the annotation, not an
error, and one that would wreck a mean.
"""
from __future__ import annotations

import glob
import json
import pathlib
import re
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "celltype"))
import panel                                        # noqa: E402

ISLETS = panel.REPO / "data" / "Islets"

#: Rings sampled per donor -- enough for a stable median without paying for all 1,304 of
#: donor 115's.  The sample is the first N by name, not a random draw, so it is reproducible.
MAX_ISLETS = 25
#: Vertices sampled per expanded polygon.
MAX_VERTICES = 150

#: What the annotation name asserts.  A donor whose measured width differs by more than
#: ``WIDTH_TOLERANCE_UM`` is reported, never silently rescaled.
NOMINAL_RING_UM = 20.0
WIDTH_TOLERANCE_UM = 1.0


def _geojson_by_donor() -> dict[str, pathlib.Path]:
    out = {}
    for path in glob.glob(str(ISLETS / "*.geojson")):
        match = re.search(r"(\d{3,4})", pathlib.Path(path).name)
        if match:
            out[match.group(1)] = pathlib.Path(path)
    return out


def measure(donor_id: str) -> dict:
    """Median expansion width in um for one donor, plus the spread it was taken over."""
    from shapely.geometry import Point, shape

    path = _geojson_by_donor().get(str(donor_id))
    if path is None:
        raise FileNotFoundError(f"no islet GeoJSON for donor {donor_id} under {ISLETS}")

    features = json.loads(path.read_text())
    features = features["features"] if isinstance(features, dict) else features
    pixel_um = _pixel_size_um(donor_id)

    core, expanded = {}, {}
    for feature in features:
        name = (feature["properties"].get("name") or "").replace("_exp20um", "")
        klass = (feature["properties"].get("classification") or {}).get("name")
        (expanded if klass == "IsletExpanded" else core)[name] = feature["geometry"]

    widths = []
    for name in sorted(core.keys() & expanded.keys())[:MAX_ISLETS]:
        inner, outer = shape(core[name]), shape(expanded[name])
        if not inner.is_valid:
            inner = inner.buffer(0)
        boundary = inner.boundary
        if outer.geom_type == "Polygon":
            vertices = np.asarray(outer.exterior.coords)
        else:
            vertices = np.vstack([np.asarray(p.exterior.coords) for p in outer.geoms])
        step = max(1, len(vertices) // MAX_VERTICES)
        distances = [boundary.distance(Point(v)) for v in vertices[::step]]
        if distances:
            widths.append(float(np.median(distances)) * pixel_um)

    n_islets_total, n_rings_total = len(core), len(expanded)
    if not widths:
        raise RuntimeError(f"donor {donor_id}: no matched Islet/IsletExpanded pair")

    widths = np.asarray(widths)
    median = float(np.median(widths))
    return {
        "donor_id": str(donor_id),
        "ring_width_um_median": median,
        "p5": float(np.percentile(widths, 5)),
        "p95": float(np.percentile(widths, 95)),
        "n_rings_measured": int(len(widths)),
        "n_rings_in_geojson": int(n_rings_total),
        "n_islets_in_geojson": int(n_islets_total),
        "matches_nominal": abs(median - NOMINAL_RING_UM) <= WIDTH_TOLERANCE_UM,
        "nominal_um": NOMINAL_RING_UM,
    }


def _pixel_size_um(donor_id: str) -> float:
    """From the QuPath manifest, per image -- never a constant."""
    manifest = json.loads((panel.PKG / "data" / "qupath_manifest.json").read_text())
    for image in manifest["images"]:
        if str(image["donor_id"]) == str(donor_id):
            x, y = image["pixel_size_um_x"], image["pixel_size_um_y"]
            if abs(x - y) > 1e-9:
                raise RuntimeError(f"donor {donor_id}: non-square pixel {x} x {y}")
            return float(x)
    raise KeyError(f"donor {donor_id} is not in the QuPath manifest")


def measure_cohort(donors: list[str] | None = None) -> dict[str, dict]:
    donors = donors or [str(d) for d in panel.cohort()["donor_id"]]
    return {d: measure(d) for d in donors}
