"""PhenoCycler QuPath annotations on the Xenium morphology grid, for Xenium Explorer.

This is the **serial-section** sibling of :mod:`islet_annotations`.  The post-Xenium IF images
the same physical section Xenium measured; the 55-plex PhenoCycler section does not.  So the
two lanes share the Xenium Explorer 4.1 output contract and nothing else about what the result
may claim.

Coordinate path, per annotation vertex::

    QuPath PhenoCycler pixels -> PhenoCycler microns -> Transform.inverse_apply
                              -> Xenium microns -> Xenium morphology level-0 pixels

The direction works out because ``[integration] fixed_modality = phenocycler``: S4 warps the
*moving* Xenium frame into the *fixed* PhenoCycler frame, so ``Transform.apply`` is
"Xenium microns -> PhenoCycler microns" -- the same stored direction the IF lane records, and
``inverse_apply`` is again the operation that carries an annotation onto the Xenium grid.

Three things this lane must handle that the same-slide lane does not.

**The Xenium section is a portion of the PhenoCycler section, and the part that matters is the
tissue, not the rectangle.**  Measured on the two donors registered so far, the Xenium extent is
39% (6436) and 54% (6442) of the PhenoCycler extent and lies entirely inside it, so annotations
outside the Xenium field have nowhere to land.  Until 2026-09-10 "the Xenium field" meant the
morphology image's bounding rectangle -- but the tissue is only **72.3%** of that rectangle on
*both* donors, so more than a quarter of the area an annotation could pass the test in was blank
background, and annotations were reaching Explorer sitting on it.  The region is now
:func:`xenium_tissue_outline`, and a polygon is kept only if it lies entirely inside it.
Straddlers are dropped whole and counted -- never clipped, which would invent a boundary that is
not an anatomical one and would misstate the annotation's area.

**The non-rigid field is not used here.**  ``reg_nonrigid`` is off for this lane, on measurement:
under the corrected transforms 15.0% (6436) and 6.1% (6442) of real annotation vertices fail the
round trip this module requires, costing 90 and 48 publishable islets, while affine-only
round-trips at 4e-12 um on the same vertices.  The earlier fields were worse than non-convergent
-- 6436's was *folded*, with a minimum Jacobian determinant of -1.657 over 0.85% of the field, so
the fixed-point inverse had no solution to find rather than converging slowly.
``use_displacement`` therefore defaults to False, and whichever ran is recorded.

**One QuPath class can exceed what Explorer will import.**  Explorer takes fewer than 10,000
polygons per file.  After dropping what falls outside the Xenium tissue, ``Blood_vessel`` still
has 76,159 (6436) and 88,578 (6442) parts, and no defensible diameter floor brings that under the
cap.  Rather than tune a filter to a display limit, a class over the cap is split into ordered
parts, largest object first, and every part is published.  ``min_diameter_um`` remains available,
defaults to 0, and is recorded when used.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import islet_annotations as IA
import numpy as np

#: Xenium Explorer imports fewer than 10,000 polygons from one third-party GeoJSON.
MAX_FEATURES_PER_FILE = IA.MAX_EXPLORER_POLYGONS

#: Affine-only inverses round-trip at ~5e-12 um, so the IF lane's strict bar is kept as the
#: default. It only needs relaxing if a caller asks for the displacement field back.
DEFAULT_ROUND_TRIP_LIMIT_UM = IA.MAX_TRANSFORM_ROUND_TRIP_UM


class PhenoAnnotationError(RuntimeError):
    """A PhenoCycler annotation cannot be placed on the Xenium grid and published."""


@dataclass(frozen=True)
class PhenoClass:
    """One publishable QuPath class, and the export directory it is drawn from."""

    name: str
    source_dir: str
    color: tuple[int, int, int]


#: The maintainer's QuPath exports, one directory per drawn class. ``Islets`` carries two
#: classes -- the islet core and its 20 um expansion -- matching ``Islet_N`` /
#: ``Islet_N_exp20um`` upstream. Colours are display metadata and carry no identity claim.
PHENO_CLASSES: tuple[PhenoClass, ...] = (
    PhenoClass("Islet", "Islets", (255, 0, 0)),
    PhenoClass("IsletExpanded", "Islets", (255, 145, 0)),
    PhenoClass("Tissue", "Tissue", (150, 150, 150)),
    PhenoClass("Blood_vessel", "Blood_vessel", (0, 90, 255)),
    PhenoClass("Lymphatic", "Lymphatic", (0, 200, 190)),
    PhenoClass("Nerve", "Nerve", (245, 215, 0)),
)

CLASSES_BY_NAME = {c.name: c for c in PHENO_CLASSES}


# --------------------------------------------------------------------------------------
# Source resolution
# --------------------------------------------------------------------------------------

def source_path(root: Path, cls: PhenoClass, image_stem: str) -> Path:
    """The QuPath export for one class of one section.

    ``image_stem`` is the qptiff-derived stem QuPath wrote into the filename, e.g.
    ``6436_Scan1.er.qptiff_-_resolution__1``. It is taken from the cohort manifest rather than
    rebuilt from the donor id, so a section whose filename does not encode its donor cannot be
    silently paired with another donor's annotations.
    """
    candidate = Path(root) / cls.source_dir / f"{cls.source_dir}__{image_stem}.geojson"
    if candidate.exists():
        return candidate
    gz = candidate.with_suffix(".geojson.gz")
    if gz.exists():
        return gz
    raise PhenoAnnotationError(
        f"no {cls.source_dir} export for image stem {image_stem!r} at {candidate}"
    )


def image_stem_for(manifest: dict[str, Any], donor: str) -> str:
    """QuPath's filename stem for a donor's PhenoCycler section, from the cohort manifest."""
    matches = [r for r in manifest["images"] if str(r.get("donor_id")) == str(donor)]
    if len(matches) != 1:
        raise PhenoAnnotationError(
            f"cohort manifest has {len(matches)} images for donor {donor}, expected exactly 1"
        )
    qptiff = matches[0].get("qptiff")
    path = qptiff.get("path") if isinstance(qptiff, dict) else qptiff
    if not path:
        raise PhenoAnnotationError(f"donor {donor} has no qptiff path in the cohort manifest")
    return f"{Path(str(path)).name}_-_resolution__1"


def pixel_size_for(manifest: dict[str, Any], donor: str) -> float:
    """Exact recorded qptiff calibration for a donor, refusing an anisotropic section."""
    matches = [r for r in manifest["images"] if str(r.get("donor_id")) == str(donor)]
    if len(matches) != 1:
        raise PhenoAnnotationError(
            f"cohort manifest has {len(matches)} images for donor {donor}, expected exactly 1"
        )
    sx = float(matches[0]["pixel_size_um_x"])
    sy = float(matches[0]["pixel_size_um_y"])
    if not np.isfinite(sx) or sx <= 0 or not np.isclose(sx, sy, rtol=0, atol=1e-9):
        raise PhenoAnnotationError(
            f"donor {donor} qptiff calibration is unusable for annotations: {sx} x {sy} um/px"
        )
    return sx


# --------------------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------------------

def _exterior_polygons(source: Path, class_name: str) -> tuple[list, dict[str, Any]]:
    """Every valid exterior polygon of one class, in source pixels, plus a full audit.

    Vectorised deliberately: ``Blood_vessel`` carries ~200k parts, and the per-polygon
    validity calls the same-slide writer makes would take minutes per class.
    """
    import shapely
    from shapely.geometry import shape

    doc = IA._read_geojson(Path(source), what="PhenoCycler annotation")
    features = doc["features"]
    if not features:
        raise PhenoAnnotationError(f"{source} carries no features")

    class_counts: Counter[str] = Counter()
    unsupported: Counter[str] = Counter()
    empty = 0
    matched_features = 0
    parts: list = []

    for index, feature in enumerate(features):
        if not isinstance(feature, dict) or feature.get("type") != "Feature":
            raise PhenoAnnotationError(f"{source}: feature {index + 1} is not a GeoJSON Feature")
        name = IA._class_name(feature)
        class_counts[name or "<unclassified>"] += 1
        if name != class_name:
            continue
        matched_features += 1
        geometry_doc = feature.get("geometry")
        if not isinstance(geometry_doc, dict):
            raise PhenoAnnotationError(
                f"{source}: {class_name} feature {index + 1} has no geometry object"
            )
        try:
            geometry = shape(geometry_doc)
        except Exception as exc:  # noqa: BLE001 - normalise malformed GeoJSON
            raise PhenoAnnotationError(
                f"{source}: cannot decode {class_name} feature {index + 1}: {exc}"
            ) from exc
        if geometry.is_empty:
            empty += 1
            continue
        if geometry.geom_type == "Polygon":
            parts.append(geometry)
        elif geometry.geom_type == "MultiPolygon":
            parts.extend(geometry.geoms)
        else:
            unsupported[geometry.geom_type] += 1

    if not parts:
        raise PhenoAnnotationError(f"{source}: no {class_name} polygon to export")

    array = np.array(parts, dtype=object)
    holes = shapely.get_num_interior_rings(array)
    # Finiteness per part, without a Python-level call per polygon: one bulk coordinate read,
    # then reduce the finite flag over each polygon's own slice.
    all_xy, owner = shapely.get_coordinates(array, return_index=True)
    finite_vertex = np.isfinite(all_xy).all(axis=1)
    finite = np.ones(len(parts), dtype=bool)
    np.logical_and.at(finite, owner, finite_vertex)
    ok = shapely.is_valid(array) & shapely.is_simple(array) & (shapely.area(array) > 0) & finite
    # A hole cannot be represented in Explorer's single-ring Polygon import, so the exterior is
    # published and the discarded interior is counted rather than silently absorbed.
    kept = array[ok]
    audit = {
        "total_source_features": int(len(features)),
        "source_classes": dict(sorted(class_counts.items())),
        "matched_source_features": int(matched_features),
        "polygon_parts": int(len(parts)),
        "excluded_empty_geometry": int(empty),
        "excluded_unsupported_geometry": int(sum(unsupported.values())),
        "unsupported_geometry_types": dict(sorted(unsupported.items())),
        "excluded_invalid_polygon": int((~ok).sum()),
        "parts_with_holes_exterior_published": int((holes[ok] > 0).sum()),
        "interior_rings_discarded": int(holes[ok].sum()),
    }
    return list(kept), audit


# --------------------------------------------------------------------------------------
# The Xenium tissue outline
# --------------------------------------------------------------------------------------

#: Defaults for the outline. ``close_um`` and ``min_area_um2`` are ``rasterize.tissue_mask``'s
#: own, kept identical so the region annotations are selected against is the same region the
#: registration QC scores Dice against. ``pixel_um`` is the raster this is built on: 10 um is
#: the grid ``register.evaluate`` already uses, and it is far finer than the ~500 um lobular
#: gaps the outline has to resolve while being coarse enough that one section is a few hundred
#: pixels on a side.
OUTLINE_PIXEL_UM = 10.0
OUTLINE_CLOSE_UM = 60.0
OUTLINE_MIN_AREA_UM2 = 50_000.0


def xenium_tissue_outline(
    xy_um: np.ndarray,
    *,
    pixel_um: float = OUTLINE_PIXEL_UM,
    close_um: float = OUTLINE_CLOSE_UM,
    min_area_um2: float = OUTLINE_MIN_AREA_UM2,
) -> tuple[Any, dict[str, Any]]:
    """Where the Xenium section actually has data, as geometry in Xenium microns.

    The exporter used to test annotations against the morphology image's bounding
    **rectangle**. That rectangle is the capture area, not the tissue: on donor 6436 it is
    7.9 x 14.5 mm of which well under half is section, so an annotation could pass the test
    while landing on blank background, which is exactly what the published figures show
    happening. This is the region that test should have been against.

    Built from Xenium **cell centroids** rather than the morphology image, for two reasons.
    Explorer draws cells, so "where Xenium looked" is honestly the cell field and not every
    pixel the camera covered. And it reuses ``rasterize.tissue_mask`` unchanged -- the same
    function ``register.evaluate`` scores tissue Dice with -- so the region this selects
    against and the region the registration is judged on cannot drift apart.

    Returns ``(geometry, provenance)``. The geometry is a shapely (Multi)Polygon in Xenium
    microns. Note that ``tissue_mask`` runs ``binary_fill_holes``, so on these sections the
    result is a filled silhouette and carries no interior rings -- an interlobular gap inside
    the section is a place Xenium imaged and found few cells, not a place it did not look, so
    filling it is the right call and an annotation there is not dropped for that reason. The
    hole handling below is kept anyway so a mask that does have them is not silently flattened.
    """
    import cv2
    import shapely

    from phenocycler.integration.rasterize import grid_for, tissue_mask

    xy = np.asarray(xy_um, dtype=np.float64)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) == 0:
        raise PhenoAnnotationError(f"need an (n, 2) array of cell centroids, got {xy.shape}")

    grid = grid_for(xy, float(pixel_um))
    mask = tissue_mask(xy, grid, close_um=float(close_um), min_area_um2=float(min_area_um2))
    if not mask.any():
        raise PhenoAnnotationError("Xenium cells produced an empty tissue mask")

    # RETR_CCOMP gives exteriors and their holes as a two-level hierarchy, which is what a
    # lobular gap is. CHAIN_APPROX_SIMPLE collapses the straight runs a rasterised boundary is
    # mostly made of; the remaining staircase is removed by the simplify below.
    contours, hierarchy = cv2.findContours(mask.astype(np.uint8), cv2.RETR_CCOMP,
                                           cv2.CHAIN_APPROX_SIMPLE)
    hierarchy = hierarchy[0] if hierarchy is not None else np.zeros((0, 4), int)

    def _ring(c: np.ndarray) -> np.ndarray:
        # findContours returns (x, y) = (col, row) pixel centres; GridSpec maps the same way.
        px = c.reshape(-1, 2).astype(np.float64)
        return grid.to_microns(px)

    shells, holes = [], []
    for i, contour in enumerate(contours):
        if len(contour) < 4:
            continue
        (shells if hierarchy[i][3] < 0 else holes).append(_ring(contour))

    if not shells:
        raise PhenoAnnotationError("tissue mask produced no exterior contour")

    outer = [shapely.Polygon(r) for r in shells]
    inner = [shapely.Polygon(r) for r in holes if len(r) >= 4]
    geom = shapely.union_all([g for g in outer if g.is_valid or g.buffer(0).is_valid])
    if inner:
        geom = geom.difference(shapely.union_all(
            [g if g.is_valid else g.buffer(0) for g in inner]))
    geom = geom.buffer(0).simplify(float(pixel_um))
    if geom.is_empty:
        raise PhenoAnnotationError("tissue outline collapsed to nothing after simplification")

    provenance = {
        "source": "xenium cell centroids",
        "n_cells": int(len(xy)),
        "pixel_um": float(pixel_um),
        "close_um": float(close_um),
        "min_area_um2": float(min_area_um2),
        "n_exterior_rings": int(len(shells)),
        "n_interior_rings": int(len(inner)),
        "area_um2": float(geom.area),
        "bounds_um": [float(v) for v in geom.bounds],
        "grid": {"origin_x_um": float(grid.origin_x_um), "origin_y_um": float(grid.origin_y_um),
                 "pixel_um": float(grid.pixel_um),
                 "width": int(grid.width), "height": int(grid.height)},
    }
    return geom, provenance


def outline_to_geojson(geom: Any, provenance: dict[str, Any]) -> dict[str, Any]:
    """The outline as a GeoJSON FeatureCollection in Xenium microns, plus its provenance."""
    import shapely

    parts = list(getattr(geom, "geoms", [geom]))
    features = []
    for i, part in enumerate(parts, start=1):
        features.append({
            "type": "Feature",
            "id": f"XeniumTissue-{i}",
            "geometry": shapely.geometry.mapping(part),
            "properties": {"objectType": "annotation", "name": f"Xenium tissue {i}",
                           "classification": {"name": "XeniumTissue", "color": [0, 200, 255]}},
        })
    return {"type": "FeatureCollection", "features": features,
            "properties": {"units": "xenium_microns", **provenance}}


def place_on_xenium_grid(
    polygons: list,
    transform,
    *,
    input_um_per_px: float,
    output_um_per_px: float,
    width: int,
    height: int,
    min_diameter_um: float = 0.0,
    round_trip_limit_um: float = DEFAULT_ROUND_TRIP_LIMIT_UM,
    tissue_outline: Any = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Carry source polygons onto the Xenium morphology grid, dropping what cannot land there.

    ``tissue_outline`` is the region from :func:`xenium_tissue_outline`, in Xenium microns. A
    polygon is kept only if it lies **entirely** inside it. Without it the only spatial test is
    the morphology image's bounding rectangle, and the tissue is 72% of that rectangle on both
    donors -- so nearly a third of the area an annotation could pass the old test in is blank
    background Xenium never imaged, which is where the stray outlines in the published figures
    came from.

    Entirely-inside, not clipped: cutting an annotation at the tissue edge would invent a
    boundary that is not an anatomical one and would misstate its area, which is the same rule
    this lane already applies at the grid edge. Straddlers are dropped whole and counted.

    Returns one record per surviving polygon (its ring in Xenium level-0 pixels and its area)
    plus the audit of everything excluded and why.
    """
    import shapely

    if not np.isfinite(input_um_per_px) or input_um_per_px <= 0:
        raise PhenoAnnotationError(f"source pixel size must be positive, got {input_um_per_px}")
    if not np.isfinite(output_um_per_px) or output_um_per_px <= 0:
        raise PhenoAnnotationError(f"Xenium pixel size must be positive, got {output_um_per_px}")
    if int(width) <= 0 or int(height) <= 0:
        raise PhenoAnnotationError(f"Xenium morphology grid is invalid: {width} x {height}")

    exteriors = np.array([p.exterior for p in polygons], dtype=object)
    coords, index = shapely.get_coordinates(exteriors, return_index=True)
    pheno_um = coords * float(input_um_per_px)

    xen_um = np.asarray(transform.inverse_apply(pheno_um, n_iter=10), dtype=np.float64)
    if xen_um.shape != pheno_um.shape or not np.isfinite(xen_um).all():
        raise PhenoAnnotationError("inverse transform produced malformed coordinates")
    round_trip = np.linalg.norm(
        np.asarray(transform.apply(xen_um), dtype=np.float64) - pheno_um, axis=1
    )
    if not np.isfinite(round_trip).all():
        raise PhenoAnnotationError("transform round trip is non-finite")

    xen_px = xen_um / float(output_um_per_px)

    n = len(polygons)
    start = np.searchsorted(index, np.arange(n))
    stop = np.searchsorted(index, np.arange(n), side="right")

    worst_round_trip = np.maximum.reduceat(round_trip, start)
    min_x = np.minimum.reduceat(xen_px[:, 0], start)
    max_x = np.maximum.reduceat(xen_px[:, 0], start)
    min_y = np.minimum.reduceat(xen_px[:, 1], start)
    max_y = np.maximum.reduceat(xen_px[:, 1], start)

    inside = (min_x >= 0.0) & (min_y >= 0.0) & (max_x <= float(width)) & (max_y <= float(height))
    converged = worst_round_trip <= float(round_trip_limit_um)

    # Containment against the tissue outline. Exact rather than a per-vertex approximation:
    # measured on donor 6436's 201,944 Blood_vessel parts the two agree on all but one, and the
    # exact test is the faster of the two, so there is nothing to buy by approximating.
    on_tissue = np.ones(n, dtype=bool)
    if tissue_outline is not None:
        outline_px = shapely.transform(tissue_outline, lambda c: c / float(output_um_per_px))
        shapely.prepare(outline_px)
        testable = (stop - start) >= 4
        if testable.any():
            order = np.full(n, -1)
            order[np.flatnonzero(testable)] = np.arange(int(testable.sum()))
            sel = testable[index]
            rings = shapely.linearrings(xen_px[sel], indices=order[index[sel]])
            on_tissue = np.zeros(n, dtype=bool)
            on_tissue[np.flatnonzero(testable)] = shapely.contains(
                outline_px, shapely.polygons(rings))
        else:
            on_tissue = np.zeros(n, dtype=bool)

    x, y = xen_px[:, 0], xen_px[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    contiguous = index == np.roll(index, -1)
    area_px2 = 0.5 * np.abs(np.add.reduceat(np.where(contiguous, cross, 0.0), start))
    area_um2 = area_px2 * float(output_um_per_px) ** 2
    diameter_um = 2.0 * np.sqrt(np.maximum(area_um2, 0.0) / np.pi)

    big_enough = diameter_um >= float(min_diameter_um) if min_diameter_um > 0 else np.ones(n, bool)
    keep = inside & on_tissue & converged & big_enough & (area_um2 > 0)

    records: list[dict[str, Any]] = []
    for i in np.flatnonzero(keep):
        ring = xen_px[start[i]:stop[i]]
        if len(ring) < 4:
            continue
        if not np.array_equal(ring[0], ring[-1]):
            ring = np.vstack([ring, ring[0]])
        records.append(
            {
                "ring_px": ring,
                "area_um2": float(area_um2[i]),
                "diameter_um": float(diameter_um[i]),
                "round_trip_um": float(worst_round_trip[i]),
                "source_part_index": int(i),
            }
        )

    audit = {
        "polygon_parts": int(n),
        "exported_polygons": int(len(records)),
        "excluded_outside_xenium_field": int((~inside).sum()),
        "excluded_outside_tissue_outline": int((inside & ~on_tissue).sum()),
        "excluded_inverse_did_not_converge": int((inside & on_tissue & ~converged).sum()),
        "excluded_below_minimum_size": int(
            (inside & on_tissue & converged & ~big_enough).sum()
        ),
        "excluded_degenerate_after_transform": int(
            (inside & on_tissue & converged & big_enough & (area_um2 <= 0)).sum()
        ),
        "minimum_equivalent_diameter_um": float(min_diameter_um),
        "round_trip_limit_um": float(round_trip_limit_um),
        "round_trip_max_um_exported": float(
            max((r["round_trip_um"] for r in records), default=0.0)
        ),
        "exported_area_um2": float(sum(r["area_um2"] for r in records)),
        "fraction_of_parts_inside_xenium_field": float(inside.mean()) if n else 0.0,
        "tissue_outline_applied": tissue_outline is not None,
        "fraction_of_parts_inside_tissue_outline": (
            float((inside & on_tissue).mean()) if n and tissue_outline is not None else None
        ),
    }
    return records, audit


def build_documents(
    records: list[dict[str, Any]],
    cls: PhenoClass,
    *,
    max_features: int = MAX_FEATURES_PER_FILE,
) -> list[dict[str, Any]]:
    """Explorer FeatureCollections for one class, largest object first, split under the cap.

    A class that fits is one document. A class that does not is split into ordered parts so
    every polygon is published; part 1 holds the largest objects, which is the layer worth
    loading first when a reviewer only wants the anatomy.
    """
    if max_features < 1:
        raise PhenoAnnotationError(f"max_features must be >= 1, got {max_features}")
    ordered = sorted(records, key=lambda r: r["area_um2"], reverse=True)
    chunks = [ordered[i:i + max_features] for i in range(0, len(ordered), max_features)] or [[]]
    documents: list[dict[str, Any]] = []
    ordinal = 0
    for part, chunk in enumerate(chunks, start=1):
        features = []
        for item in chunk:
            ordinal += 1
            features.append(
                {
                    "type": "Feature",
                    "id": f"{cls.name}-{ordinal}",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [[float(px), float(py)] for px, py in item["ring_px"]]
                        ],
                    },
                    "properties": {
                        "objectType": "annotation",
                        "name": f"{cls.name} {ordinal}",
                        "classification": {"name": cls.name, "color": list(cls.color)},
                    },
                }
            )
        documents.append(
            {
                "part": part,
                "n_parts": len(chunks),
                "document": {"type": "FeatureCollection", "features": features},
                "n_features": len(features),
                "area_um2": float(sum(i["area_um2"] for i in chunk)),
            }
        )
    return documents


def export_class(
    source: Path,
    destination_dir: Path,
    transform,
    cls: PhenoClass,
    *,
    serial: str,
    input_um_per_px: float,
    output_um_per_px: float,
    width: int,
    height: int,
    min_diameter_um: float = 0.0,
    round_trip_limit_um: float = DEFAULT_ROUND_TRIP_LIMIT_UM,
    max_features: int = MAX_FEATURES_PER_FILE,
    used_displacement: bool = False,
    tissue_outline: Any = None,
) -> dict[str, Any]:
    """Build, write, reopen and validate every sidecar for one class."""
    polygons, source_audit = _exterior_polygons(Path(source), cls.name)
    records, place_audit = place_on_xenium_grid(
        polygons,
        transform,
        input_um_per_px=input_um_per_px,
        output_um_per_px=output_um_per_px,
        width=width,
        height=height,
        min_diameter_um=min_diameter_um,
        round_trip_limit_um=round_trip_limit_um,
        tissue_outline=tissue_outline,
    )
    if not records:
        raise PhenoAnnotationError(
            f"{source}: no {cls.name} polygon survives placement on the Xenium grid "
            f"({place_audit['excluded_outside_xenium_field']} outside the field, "
            f"{place_audit['excluded_outside_tissue_outline']} off the tissue outline, "
            f"{place_audit['excluded_inverse_did_not_converge']} non-convergent)"
        )

    documents = build_documents(records, cls, max_features=max_features)
    destination_dir = Path(destination_dir)
    destination_dir.mkdir(parents=True, exist_ok=True)

    written: list[dict[str, Any]] = []
    for entry in documents:
        # "pheno" is load-bearing in the name: the same bundle already carries a same-slide
        # IF sidecar, and a reviewer must be able to tell a serial-section annotation from a
        # same-section one in Explorer's layer list without opening the QC file.
        if entry["n_parts"] == 1:
            filename = f"{serial}_pheno_{cls.name}_on_xenium_grid.geojson"
        else:
            filename = (
                f"{serial}_pheno_{cls.name}_on_xenium_grid_part{entry['part']:02d}"
                f"of{entry['n_parts']:02d}.geojson"
            )
        path = destination_dir / filename
        path.write_text(
            json.dumps(entry["document"], separators=(",", ":"), allow_nan=False) + "\n"
        )
        validation = IA.validate_islet_geojson(
            path,
            width=width,
            height=height,
            expected_features=entry["n_features"],
            class_name=cls.name,
            color=cls.color,
        )
        written.append(
            {
                "file": path.name,
                "part": entry["part"],
                "n_parts": entry["n_parts"],
                "n_features": entry["n_features"],
                "area_um2": entry["area_um2"],
                "sha256": validation["sha256"],
                "bytes": validation["bytes"],
                "validation": validation,
            }
        )

    return {
        "class": cls.name,
        "color_rgb": list(cls.color),
        "source_file": str(Path(source)),
        "source_sha256": IA._sha256(Path(source)),
        "source_coordinate_space": "PhenoCycler qptiff level-0 pixel edges",
        "source_um_per_px": float(input_um_per_px),
        "output_coordinate_space": "Xenium morphology level-0 pixel edges",
        "output_um_per_px": float(output_um_per_px),
        "morphology_width_px": int(width),
        "morphology_height_px": int(height),
        "transform_direction_stored": "Xenium microns -> PhenoCycler microns",
        "transform_operation_used": (
            "Transform.inverse_apply(PhenoCycler microns) -> Xenium microns"
        ),
        "includes_displacement": bool(used_displacement),
        "section_relationship": "serial",
        "source_audit": source_audit,
        "placement_audit": place_audit,
        "files": written,
        "n_files": len(written),
        "n_features_total": int(sum(f["n_features"] for f in written)),
        "import_route": (
            "Annotations -> Import GeoJSON; associate with the Xenium In Situ morphology image"
        ),
    }
