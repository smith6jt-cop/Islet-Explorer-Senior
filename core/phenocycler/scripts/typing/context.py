"""Where a cell is: which classified object its centroid falls in.

Region membership says *where* a cell is, not *what* it is -- the question a context asks is
answered by that context's own markers, not by the membership.  What membership does decide
is *which* question gets asked.

The objects are QuPath pixel-classifier output (`B3Tub_Ero0.5_Thr600`,
`PDPNmorphClos0-5_Thr300`, `CD34dil1430` / `Cavdil1200`) in PIXEL coordinates, while cell
centroids are in micrometres, so every join needs that image's own pixel size.  Assignment
validates against the markers the classifiers were built from with nothing tuned: on 6442,
B3TUBB is 13.2x enriched inside Nerve and Ker8-18 depleted to 0.13x, Caveolin 13.6x and
CD34 9.8x inside Blood vessel, Podoplanin 3.1x inside Lymphatic.
"""
from __future__ import annotations

import glob
import json

import numpy as np
import pandas as pd
from shapely import STRtree, points
from shapely.geometry import shape

from subpanel_io import GEOJSON

#: Folder under `results/Geojson/` for each context name.
OBJECT_FOLDER = {"Nerve": "Nerve", "Lymphatic": "Lymphatic", "Blood vessel": "Blood_vessel"}

#: Precedence when a centroid falls inside more than one object, strongest first.
#: Recovered from the classification on disk rather than assumed -- see `verify_against_calls`.
PRECEDENCE = ("Blood vessel", "Nerve", "Lymphatic")


def load_objects(image_stem: str, name: str) -> list:
    """Every valid outline of one classified class for one image, in pixel coordinates."""
    matches = [p for p in glob.glob(f"{GEOJSON}/{OBJECT_FOLDER[name]}/*.geojson")
               if image_stem in p]
    if len(matches) != 1:
        raise FileNotFoundError(
            f"{name} outlines for {image_stem!r}: expected 1 file, found {len(matches)}")
    with open(matches[0]) as handle:
        features = json.load(handle)["features"]
    out = []
    for feature in features:
        if not feature.get("geometry"):
            continue
        geom = shape(feature["geometry"])
        out.append(geom if geom.is_valid else geom.buffer(0))
    return out


def assign(x_um: np.ndarray, y_um: np.ndarray, in_islet: np.ndarray,
           image_stem: str, pixel_size_um: float) -> np.ndarray:
    """Context per cell: islet, one of the object classes, or tissue.

    Islet membership comes from QuPath's own Parent column rather than from a centroid test,
    because that is what drew the islet in the first place.
    """
    context = np.where(np.asarray(in_islet), "islet", "tissue").astype(object)
    candidate = ~np.asarray(in_islet)
    if not candidate.any():
        return context
    idx = np.flatnonzero(candidate)
    pts = points(np.column_stack([x_um[idx] / pixel_size_um, y_um[idx] / pixel_size_um]))

    # Weakest first, so a stronger class overwrites it -- one pass, no tie-breaking logic.
    for name in reversed(PRECEDENCE):
        outlines = load_objects(image_stem, name)
        if not outlines:
            continue
        hit = STRtree(outlines).query(pts, predicate="within")[0]
        context[idx[np.unique(hit)]] = name
    return context


def verify_against_calls(predicted: np.ndarray, recorded: pd.Series) -> pd.DataFrame:
    """Confusion of predicted context against the classification already on disk.

    The reconstruction has to reproduce the existing table before any correction is applied;
    otherwise the rebuild and the correction are confounded and neither is interpretable.
    """
    return pd.crosstab(pd.Series(predicted, name="predicted"),
                       pd.Series(np.asarray(recorded), name="recorded"))
