"""Transform post-Xenium IF islet annotations for Xenium Explorer.

QuPath writes annotation vertices in level-0 qptiff pixel coordinates.  Registration
``Transform`` objects in this lane map the other way -- Xenium microns to IF microns -- so each
annotation vertex follows this explicit path::

    QuPath IF pixels -> IF microns -> Transform.inverse_apply -> Xenium microns
                     -> Xenium morphology level-0 pixels

The last conversion matters.  Xenium Explorer 4.1 imports third-party GeoJSON relative to an
image, and these annotations are associated with the Xenium In Situ morphology image.  Emitting
microns instead would place them about 4.7-fold away from the registered OME-TIFF, whose pixels
already coincide with that morphology grid.

Only ``Islet`` polygons at least 50 microns in equivalent diameter are exported.  The QuPath
pixel classifier also emits sub-cell fragments (and occasionally line objects); their omission
is counted in the returned provenance rather than being silent.  The writer enforces Explorer's
documented Polygon-only limits before a file can be published.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import geo
import numpy as np

ISLET_CLASS = geo.ISLET
ISLET_COLOR = (255, 0, 0)
MIN_ISLET_DIAMETER_UM = geo.MIN_ISLET_DIAMETER_UM
MAX_EXPLORER_POLYGONS = 9_999
MAX_EXPLORER_VERTICES = 99_999
MAX_TRANSFORM_ROUND_TRIP_UM = 1e-4


class IsletAnnotationError(RuntimeError):
    """An islet sidecar cannot be imported faithfully into Xenium Explorer."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_geojson(path: Path, *, what: str) -> dict[str, Any]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt") as handle:
            doc = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise IsletAnnotationError(f"cannot read {what} {path}: {exc}") from exc
    if not isinstance(doc, dict) or doc.get("type") != "FeatureCollection":
        raise IsletAnnotationError(f"{what} {path} must be one GeoJSON FeatureCollection")
    if not isinstance(doc.get("features"), list):
        raise IsletAnnotationError(f"{what} {path} has no list-valued features member")
    return doc


def _class_name(feature: Mapping[str, Any]) -> str:
    properties = feature.get("properties")
    if not isinstance(properties, Mapping):
        return ""
    classification = properties.get("classification")
    if isinstance(classification, Mapping):
        name = classification.get("name")
    else:
        name = classification
    return "" if name is None else str(name)


def _valid_source_polygon(polygon, *, source: Path, feature_index: int, part_index: int) -> None:
    import shapely

    label = f"{source}: Islet feature {feature_index + 1}, part {part_index + 1}"
    coordinates = np.asarray(shapely.get_coordinates(polygon), dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2 or not np.isfinite(coordinates).all():
        raise IsletAnnotationError(f"{label} contains malformed or non-finite coordinates")
    if polygon.is_empty or polygon.area <= 0:
        raise IsletAnnotationError(f"{label} is empty or has zero area")
    if not polygon.is_valid or not polygon.is_simple:
        reason = shapely.is_valid_reason(polygon)
        raise IsletAnnotationError(f"{label} is not a simple valid polygon: {reason}")
    if len(polygon.interiors):
        raise IsletAnnotationError(
            f"{label} has {len(polygon.interiors)} hole(s), which Xenium Explorer cannot import"
        )


def _source_islet_polygons(
    source: Path,
    *,
    input_um_per_px: float,
    min_diameter_um: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Load, audit, scale, and size-filter source Islet polygons."""
    import shapely
    from shapely.geometry import shape

    source = Path(source)
    if not np.isfinite(input_um_per_px) or input_um_per_px <= 0:
        raise IsletAnnotationError(
            f"source qptiff pixel size must be positive and finite, got {input_um_per_px}"
        )
    if not np.isfinite(min_diameter_um) or min_diameter_um <= 0:
        raise IsletAnnotationError(
            f"minimum islet diameter must be positive and finite, got {min_diameter_um}"
        )

    doc = _read_geojson(source, what="source annotation")
    features = doc["features"]
    class_counts: Counter[str] = Counter()
    unsupported: Counter[str] = Counter()
    empty_geometries = 0
    islet_features = 0
    polygon_parts = 0
    below_size = 0
    kept: list[dict[str, Any]] = []
    minimum_area_um2 = float(np.pi * (min_diameter_um / 2.0) ** 2)

    for feature_index, feature in enumerate(features):
        if not isinstance(feature, Mapping) or feature.get("type") != "Feature":
            raise IsletAnnotationError(
                f"{source}: feature {feature_index + 1} is not a GeoJSON Feature"
            )
        class_name = _class_name(feature)
        class_counts[class_name or "<unclassified>"] += 1
        if class_name != ISLET_CLASS:
            continue
        islet_features += 1
        geometry_doc = feature.get("geometry")
        if not isinstance(geometry_doc, Mapping):
            raise IsletAnnotationError(
                f"{source}: Islet feature {feature_index + 1} has no geometry object"
            )
        try:
            geometry = shape(geometry_doc)
        except Exception as exc:  # noqa: BLE001 - normalize malformed GeoJSON geometry
            raise IsletAnnotationError(
                f"{source}: cannot decode Islet feature {feature_index + 1}: {exc}"
            ) from exc
        if geometry.is_empty:
            empty_geometries += 1
            continue
        if geometry.geom_type == "Polygon":
            parts = [geometry]
        elif geometry.geom_type == "MultiPolygon":
            parts = list(geometry.geoms)
        else:
            unsupported[geometry.geom_type] += 1
            continue

        for part_index, polygon_px in enumerate(parts):
            polygon_parts += 1
            _valid_source_polygon(
                polygon_px,
                source=source,
                feature_index=feature_index,
                part_index=part_index,
            )
            area_um2 = float(polygon_px.area * input_um_per_px**2)
            if area_um2 < minimum_area_um2:
                below_size += 1
                continue
            xy_if_um = (
                np.asarray(shapely.get_coordinates(polygon_px.exterior), dtype=np.float64)
                * input_um_per_px
            )
            kept.append(
                {
                    "source_feature_index": int(feature_index),
                    "source_feature_id": (
                        None if feature.get("id") is None else str(feature.get("id"))
                    ),
                    "source_part_index": int(part_index),
                    "source_area_um2": area_um2,
                    "xy_if_um": xy_if_um,
                }
            )

    if not features:
        raise IsletAnnotationError(f"source annotation {source} carries no features")
    if not kept:
        raise IsletAnnotationError(
            f"{source}: no Islet polygon meets the {min_diameter_um:g} um diameter filter"
        )
    if len(kept) > MAX_EXPLORER_POLYGONS:
        raise IsletAnnotationError(
            f"{source}: {len(kept)} retained polygons exceed Xenium Explorer's limit of "
            f"fewer than {MAX_EXPLORER_POLYGONS + 1}"
        )

    audit = {
        "total_source_features": int(len(features)),
        "source_classes": dict(sorted(class_counts.items())),
        "islet_source_features": int(islet_features),
        "islet_polygon_parts": int(polygon_parts),
        "exported_islet_polygons": int(len(kept)),
        "excluded_below_minimum_size": int(below_size),
        "excluded_unsupported_geometry": int(sum(unsupported.values())),
        "unsupported_geometry_types": dict(sorted(unsupported.items())),
        "excluded_empty_geometry": int(empty_geometries),
        "excluded_non_islet_class": int(len(features) - islet_features),
    }
    return kept, audit


def build_islet_geojson(
    source: Path,
    transform,
    *,
    input_um_per_px: float,
    output_um_per_px: float,
    width: int,
    height: int,
    min_diameter_um: float = MIN_ISLET_DIAMETER_UM,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build an Explorer GeoJSON document and complete transformation provenance."""
    from shapely.geometry import Polygon

    if not np.isfinite(output_um_per_px) or output_um_per_px <= 0:
        raise IsletAnnotationError(
            f"Xenium morphology pixel size must be positive and finite, got {output_um_per_px}"
        )
    if int(width) <= 0 or int(height) <= 0:
        raise IsletAnnotationError(f"Xenium morphology dimensions are invalid: {width} x {height}")

    retained, counts = _source_islet_polygons(
        source,
        input_um_per_px=float(input_um_per_px),
        min_diameter_um=float(min_diameter_um),
    )
    features: list[dict[str, Any]] = []
    max_round_trip_um = 0.0
    total_source_area_um2 = 0.0
    total_output_area_um2 = 0.0

    for ordinal, item in enumerate(retained, start=1):
        xy_if_um = item["xy_if_um"]
        try:
            xy_xenium_um = np.asarray(transform.inverse_apply(xy_if_um, n_iter=10), dtype=float)
            reconstructed_if_um = np.asarray(transform.apply(xy_xenium_um), dtype=float)
        except Exception as exc:  # noqa: BLE001 - normalize transform failures with source context
            raise IsletAnnotationError(
                f"cannot inverse-transform retained Islet {ordinal} from {source}: {exc}"
            ) from exc
        if xy_xenium_um.shape != xy_if_um.shape or not np.isfinite(xy_xenium_um).all():
            raise IsletAnnotationError(
                f"inverse transform produced malformed coordinates for retained Islet {ordinal}"
            )
        round_trip = np.linalg.norm(reconstructed_if_um - xy_if_um, axis=1)
        if not np.isfinite(round_trip).all():
            raise IsletAnnotationError(
                f"transform round trip is non-finite for retained Islet {ordinal}"
            )
        feature_round_trip_um = float(round_trip.max(initial=0.0))
        if feature_round_trip_um > MAX_TRANSFORM_ROUND_TRIP_UM:
            raise IsletAnnotationError(
                f"full inverse transform misses retained Islet {ordinal} by up to "
                f"{feature_round_trip_um:.6g} um after round trip; limit is "
                f"{MAX_TRANSFORM_ROUND_TRIP_UM:g} um"
            )
        max_round_trip_um = max(max_round_trip_um, feature_round_trip_um)

        xy_xenium_px = xy_xenium_um / float(output_um_per_px)
        polygon = Polygon(xy_xenium_px)
        _valid_source_polygon(
            polygon,
            source=Path(source),
            feature_index=int(item["source_feature_index"]),
            part_index=int(item["source_part_index"]),
        )
        ring = [[float(x), float(y)] for x, y in polygon.exterior.coords]
        if len(ring) > MAX_EXPLORER_VERTICES:
            raise IsletAnnotationError(
                f"retained Islet {ordinal} has {len(ring)} vertices; Xenium Explorer requires "
                f"fewer than {MAX_EXPLORER_VERTICES + 1}"
            )
        total_source_area_um2 += float(item["source_area_um2"])
        total_output_area_um2 += float(polygon.area * output_um_per_px**2)
        properties: dict[str, Any] = {
            "objectType": "annotation",
            "name": f"Islet {ordinal}",
            "classification": {"name": ISLET_CLASS, "color": list(ISLET_COLOR)},
        }
        features.append(
            {
                "type": "Feature",
                "id": f"Islet-{ordinal}",
                "geometry": {"type": "Polygon", "coordinates": [ring]},
                "properties": properties,
            }
        )

    document = {"type": "FeatureCollection", "features": features}
    provenance = {
        "status": "generated",
        "source_file": str(Path(source)),
        "source_sha256": _sha256(Path(source)),
        "source_coordinate_space": "post-Xenium IF qptiff level-0 pixel edges",
        "source_um_per_px": float(input_um_per_px),
        "class": ISLET_CLASS,
        "minimum_equivalent_diameter_um": float(min_diameter_um),
        "output_coordinate_space": "Xenium morphology level-0 pixel edges",
        "output_um_per_px": float(output_um_per_px),
        "morphology_width_px": int(width),
        "morphology_height_px": int(height),
        "transform_direction_stored": "Xenium microns -> IF microns",
        "transform_operation_used": "Transform.inverse_apply(IF microns) -> Xenium microns",
        "includes_displacement": bool(getattr(transform, "displacement", None) is not None),
        "transform_round_trip_max_um": float(max_round_trip_um),
        "transform_round_trip_limit_um": MAX_TRANSFORM_ROUND_TRIP_UM,
        "counts": counts,
        "retained_source_area_um2": float(total_source_area_um2),
        "transformed_area_um2": float(total_output_area_um2),
        "import_route": (
            "Annotations -> Import GeoJSON; associate with the Xenium In Situ morphology image"
        ),
        "alignment_file": None,
    }
    return document, provenance


def validate_islet_geojson(
    path: Path,
    *,
    width: int,
    height: int,
    expected_features: int | None = None,
) -> dict[str, Any]:
    """Reopen and enforce the Xenium Explorer 4.1 third-party annotation contract."""
    import shapely
    from shapely.geometry import Polygon

    path = Path(path)
    doc = _read_geojson(path, what="exported islet annotation")
    features = doc["features"]
    if not features:
        raise IsletAnnotationError(f"{path}: Xenium Explorer cannot import an empty annotation")
    if len(features) > MAX_EXPLORER_POLYGONS:
        raise IsletAnnotationError(
            f"{path}: {len(features)} polygons exceed Xenium Explorer's supported count"
        )
    if expected_features is not None and len(features) != int(expected_features):
        raise IsletAnnotationError(
            f"{path}: has {len(features)} features, expected {int(expected_features)}"
        )

    identifiers: set[str] = set()
    minimum = np.array([np.inf, np.inf], dtype=float)
    maximum = np.array([-np.inf, -np.inf], dtype=float)
    largest_vertex_count = 0
    total_vertices = 0
    for index, feature in enumerate(features):
        label = f"{path}: feature {index + 1}"
        if not isinstance(feature, Mapping) or feature.get("type") != "Feature":
            raise IsletAnnotationError(f"{label} is not a GeoJSON Feature")
        identifier = feature.get("id")
        if not isinstance(identifier, str) or not identifier or identifier in identifiers:
            raise IsletAnnotationError(f"{label} has an absent or duplicate string id")
        identifiers.add(identifier)
        geometry = feature.get("geometry")
        if not isinstance(geometry, Mapping) or geometry.get("type") != "Polygon":
            kind = geometry.get("type") if isinstance(geometry, Mapping) else None
            raise IsletAnnotationError(f"{label} geometry is {kind!r}, expected 'Polygon'")
        coordinates = geometry.get("coordinates")
        if not isinstance(coordinates, list) or len(coordinates) != 1:
            raise IsletAnnotationError(f"{label} must be one Polygon exterior with no holes")
        ring = coordinates[0]
        if not isinstance(ring, list) or not (4 <= len(ring) <= MAX_EXPLORER_VERTICES):
            raise IsletAnnotationError(
                f"{label} must have 4..{MAX_EXPLORER_VERTICES} exterior coordinates"
            )
        if ring[0] != ring[-1]:
            raise IsletAnnotationError(f"{label} is not explicitly closed")
        try:
            xy = np.asarray(ring, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise IsletAnnotationError(f"{label} has non-numeric coordinates") from exc
        if xy.shape != (len(ring), 2) or not np.isfinite(xy).all():
            raise IsletAnnotationError(f"{label} coordinates must be finite [x, y] pairs")
        if (
            np.any(xy[:, 0] < 0)
            or np.any(xy[:, 0] > int(width))
            or np.any(xy[:, 1] < 0)
            or np.any(xy[:, 1] > int(height))
        ):
            raise IsletAnnotationError(
                f"{label} lies outside the {int(width)} x {int(height)} morphology grid"
            )
        polygon = Polygon(xy)
        if polygon.is_empty or polygon.area <= 0 or not polygon.is_valid or not polygon.is_simple:
            raise IsletAnnotationError(
                f"{label} is not a simple valid polygon: {shapely.is_valid_reason(polygon)}"
            )
        properties = feature.get("properties")
        if not isinstance(properties, Mapping) or properties.get("objectType") != "annotation":
            raise IsletAnnotationError(f"{label} properties.objectType must be 'annotation'")
        classification = properties.get("classification")
        if not isinstance(classification, Mapping) or classification.get("name") != ISLET_CLASS:
            raise IsletAnnotationError(f"{label} classification.name must be {ISLET_CLASS!r}")
        color = classification.get("color")
        if (
            not isinstance(color, list)
            or len(color) != 3
            or any(type(value) is not int or not 0 <= value <= 255 for value in color)
            or tuple(color) != ISLET_COLOR
        ):
            raise IsletAnnotationError(
                f"{label} classification.color must be {list(ISLET_COLOR)}"
            )
        minimum = np.minimum(minimum, xy.min(axis=0))
        maximum = np.maximum(maximum, xy.max(axis=0))
        largest_vertex_count = max(largest_vertex_count, len(ring))
        total_vertices += len(ring)

    return {
        "status": "PASS",
        "checks": {
            "feature_collection": True,
            "polygon_only": True,
            "closed_rings": True,
            "no_holes": True,
            "simple_valid_polygons": True,
            "finite_xy_pairs": True,
            "inside_morphology_grid": True,
            "object_type_annotation": True,
            "classification": ISLET_CLASS,
            "classification_color_rgb": list(ISLET_COLOR),
            "feature_count_below_10000": True,
            "vertices_per_polygon_below_100000": True,
        },
        "n_features": int(len(features)),
        "n_vertices": int(total_vertices),
        "max_vertices_per_polygon": int(largest_vertex_count),
        "bounds_px": [float(minimum[0]), float(minimum[1]), float(maximum[0]), float(maximum[1])],
        "bytes": int(path.stat().st_size),
        "sha256": _sha256(path),
    }


def export_islet_geojson(
    source: Path,
    destination: Path,
    transform,
    *,
    input_um_per_px: float,
    output_um_per_px: float,
    width: int,
    height: int,
    min_diameter_um: float = MIN_ISLET_DIAMETER_UM,
) -> dict[str, Any]:
    """Build, write, reopen, and validate one Explorer annotation sidecar."""
    destination = Path(destination)
    document, provenance = build_islet_geojson(
        source,
        transform,
        input_um_per_px=input_um_per_px,
        output_um_per_px=output_um_per_px,
        width=width,
        height=height,
        min_diameter_um=min_diameter_um,
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(document, separators=(",", ":"), allow_nan=False) + "\n")
    validation = validate_islet_geojson(
        destination,
        width=width,
        height=height,
        expected_features=provenance["counts"]["exported_islet_polygons"],
    )
    return {
        **provenance,
        "file": str(destination),
        "output_sha256": validation["sha256"],
        "validation": validation,
    }
