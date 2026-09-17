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

    `in_islet` is supplied by the caller and MUST come from `islet_membership.membership()`,
    which decides it by testing the centroid against the islet outlines.

    It used to come from QuPath's `Parent` column, on the reasoning that the Parent "is what drew
    the islet in the first place". That was measured wrong on 2026-09-15: the Parent labels cells
    `peri` that sit up to 74.6 um inside a core outline, and because a false `in_islet` means the
    cell is never asked the islet question at all, 11,556 islet-core cells cohort-wide fell
    through to Exocrine -- 5,250 of them hormone-positive. Do not reinstate it.
    """
    inside = memberships(x_um, y_um, image_stem, pixel_size_um)
    context = np.where(np.asarray(in_islet), "islet", "tissue").astype(object)
    candidate = ~np.asarray(in_islet)
    # Weakest first, so a stronger class overwrites it -- one pass, no tie-breaking logic.
    for name in reversed(PRECEDENCE):
        context[candidate & inside[name]] = name
    return context


def memberships(x_um: np.ndarray, y_um: np.ndarray, image_stem: str,
                pixel_size_um: float) -> dict:
    """class -> boolean mask, is each centroid inside an object of that class.

    Kept SEPARATE from :func:`assign`, which collapses these by precedence into one label, because
    the classification needs both and they are not interchangeable. A cell can sit inside a blood
    vessel object AND a lymphatic one; `assign` reports "Blood vessel" and the membership in the
    lymphatic is gone. But the decision rule asks each containing object's OWN markers whether
    they clear that object's own level, and precedence breaks ties only among the ones that do --
    so a cell in both, whose Caveolin/CD34 are at background and whose Podoplanin is not, is
    Lymphatic. Measured on 6442: 1,364 of 6,672 Lymphatic calls and 152 of 1,215 Nerve calls are
    exactly that case, and reading them off the collapsed label instead loses every one.
    """
    pts = points(np.column_stack([np.asarray(x_um) / pixel_size_um,
                                  np.asarray(y_um) / pixel_size_um]))
    out = {}
    for name in PRECEDENCE:
        mask = np.zeros(len(pts), dtype=bool)
        outlines = load_objects(image_stem, name)
        if outlines:
            mask[np.unique(STRtree(outlines).query(pts, predicate="within")[0])] = True
        out[name] = mask
    return out


def verify_against_calls(predicted: np.ndarray, recorded: pd.Series) -> pd.DataFrame:
    """Confusion of predicted context against the classification already on disk.

    The reconstruction has to reproduce the existing table before any correction is applied;
    otherwise the rebuild and the correction are confounded and neither is interpretable.
    """
    return pd.crosstab(pd.Series(predicted, name="predicted"),
                       pd.Series(np.asarray(recorded), name="recorded"))
