"""QuPath annotations exported from the post-Xenium IF, in microns.

The post-Xenium IF images the **same physical section** Xenium measured, so the correspondence
is exact by construction rather than inferred -- and the Xenium segmentation is already the
segmentation, which is why nothing here needs the re-stains segmented in QuPath.

The maintainer's older registration exports contain a ``Tissue`` outline plus thousands of
``Islet`` classifier polygons; newer downstream exports may be islet-only.  Both arrive in
**full-resolution qptiff pixels**, and this historical loader uses the cohort's recorded
``PixelSizeMicrons`` 0.25 (Vectra / PhenoImager Fusion 2.2.0, 7 pyramid levels).  The current
Explorer sidecar path in ``islet_annotations.py`` instead receives each inventory record's exact
qptiff calibration, so a future acquisition cannot silently inherit this historical constant.

**The ``Islet`` class needs a size filter and it is not optional.**  The classifier emits
thousands of sub-cell fragments alongside the real islets.  Measured on the two exported
donors:

    donor 6436: 4,671 objects, median area    27 um^2 (5.9 um across)
    donor 6442: 3,366 objects, median area    23 um^2 (5.4 um across)

At >= 50 um equivalent diameter that becomes **125** and **95** objects carrying **85%** and
**81%** of the total islet area -- which is an islet count in the range the PhenoCycler cohort
shows (67-1,304 per donor), and far more landmarks than a point-set registration needs.  The
median object is one cell wide, so an unfiltered centroid set is mostly noise.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

#: Recorded in the qptiff as ``PixelSizeMicrons``; asserted rather than assumed by
#: :func:`assert_scale_is_plausible`, because getting it wrong does not error -- it puts one
#: section's annotations a whole tissue away from its own cells.
PIXEL_UM = 0.25

#: Below this equivalent diameter an ``Islet`` object is a classifier fragment, not an islet.
#: Declared in microns rather than pixels so it does not silently change meaning if a future
#: export comes from a different pyramid level.
MIN_ISLET_DIAMETER_UM = 50.0

TISSUE, ISLET = "Tissue", "Islet"


class RestainGeoError(RuntimeError):
    """An export cannot be read as the geometry it claims to be -- never silently skipped."""


def load_annotations(path: Path) -> dict[str, list]:
    """Class name -> list of shapely geometries, in MICRONS.

    Reads ``.geojson`` or ``.geojson.gz``.  Classes other than Tissue/Islet are returned too
    rather than dropped, so a future export gaining a class is visible instead of ignored.
    """
    from shapely.affinity import scale as shp_scale
    from shapely.geometry import shape

    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        doc = json.load(fh)
    feats = doc["features"] if isinstance(doc, dict) and "features" in doc else doc
    if not feats:
        raise RestainGeoError(f"{path} carries no features")

    out: dict[str, list] = {}
    for ft in feats:
        cls = (ft.get("properties", {}).get("classification") or {})
        name = cls.get("name") if isinstance(cls, dict) else str(cls)
        geom = shp_scale(shape(ft["geometry"]), xfact=PIXEL_UM, yfact=PIXEL_UM, origin=(0, 0))
        out.setdefault(name or "", []).append(geom)
    return out


def has_tissue(ann: dict[str, list]) -> bool:
    """Whether annotations contain at least one ``Tissue`` geometry.

    A GeoJSON file being present is not equivalent to a tissue annotation being present: an
    islet-only export is valid input for islet windows but cannot define a tissue mask.
    """
    return bool(ann.get(TISSUE))


def tissue_polygon(ann: dict[str, list]):
    """The one tissue footprint, as a single geometry."""
    from shapely.ops import unary_union
    if not has_tissue(ann):
        raise RestainGeoError(f"no {TISSUE!r} annotation; classes present: {sorted(ann)}")
    return unary_union(ann[TISSUE])


def islets(ann: dict[str, list], *, min_diameter_um: float = MIN_ISLET_DIAMETER_UM) -> list:
    """The ``Islet`` objects big enough to be islets -- see the module docstring."""
    min_area = np.pi * (min_diameter_um / 2.0) ** 2
    return [g for g in ann.get(ISLET, []) if g.area >= min_area]


def centroids(geoms) -> np.ndarray:
    return np.array([[g.centroid.x, g.centroid.y] for g in geoms], dtype=float).reshape(-1, 2)


def assert_scale_is_plausible(tissue, xy_um: np.ndarray, *, tol: float = 0.25) -> dict:
    """Refuse a px->um conversion that puts the two modalities on different scales.

    The failure this exists for is silent: a wrong ``PIXEL_UM`` does not raise, it produces a
    registration that cannot converge and an overlay that looks merely misaligned.  Comparing
    the two footprints' extents catches an order-of-magnitude error immediately and a 2x error
    comfortably, while tolerating the few percent that a drawn outline genuinely differs from
    a cell footprint by.
    """
    bx0, by0, bx1, by1 = tissue.bounds
    ext_if = np.array([bx1 - bx0, by1 - by0])
    ext_xen = np.ptp(np.asarray(xy_um, dtype=float), axis=0)
    ratio = ext_if / np.where(ext_xen > 0, ext_xen, np.nan)
    if not np.all(np.abs(ratio - 1.0) <= tol):
        raise RestainGeoError(
            f"IF tissue extent {ext_if.round(0)} um vs Xenium {ext_xen.round(0)} um "
            f"(ratio {ratio.round(3)}) is outside +/-{tol:.0%}. Either PIXEL_UM={PIXEL_UM} is "
            "wrong for this export, or the geojson and the bundle are not the same section."
        )
    return {"extent_if_um": ext_if.tolist(), "extent_xenium_um": ext_xen.tolist(),
            "extent_ratio": ratio.tolist()}


def rasterize(geoms, grid) -> np.ndarray:
    """Binary mask of ``geoms`` (microns) on a :class:`GridSpec`.

    PIL rather than rasterio: the repo does not carry rasterio, a filled polygon is all this
    needs, and interior rings are honoured so a tissue outline with holes stays holed.
    """
    from PIL import Image, ImageDraw
    from shapely.affinity import affine_transform

    img = Image.new("1", (grid.width, grid.height), 0)
    draw = ImageDraw.Draw(img)
    s = 1.0 / grid.pixel_um
    for g in geoms:
        g = affine_transform(g, [s, 0, 0, s,
                                 -grid.origin_x_um / grid.pixel_um,
                                 -grid.origin_y_um / grid.pixel_um])
        for poly in ([g] if g.geom_type == "Polygon" else list(g.geoms)):
            draw.polygon([tuple(c) for c in poly.exterior.coords], fill=1)
            for ring in poly.interiors:
                draw.polygon([tuple(c) for c in ring.coords], fill=0)
    return np.array(img, dtype=bool)
