"""Xenium Explorer annotation export contracts on planted geometry."""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import islet_annotations as IA  # noqa: E402

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402


def _feature(geometry, *, cls="Islet", identifier="source-1"):
    return {
        "type": "Feature",
        "id": identifier,
        "geometry": geometry,
        "properties": {
            "objectType": "annotation",
            "classification": {"name": cls, "color": [255, 0, 0]},
        },
    }


def _polygon(x0, y0, side):
    return {
        "type": "Polygon",
        "coordinates": [
            [
                [x0, y0],
                [x0 + side, y0],
                [x0 + side, y0 + side],
                [x0, y0 + side],
                [x0, y0],
            ]
        ],
    }


def _write_source(path, features):
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))
    return path


def test_qupath_pixels_are_inverse_transformed_to_xenium_grid_pixels(tmp_path):
    source = _write_source(tmp_path / "source.geojson", [_feature(_polygon(400, 800, 240))])
    # Stored direction is Xenium um -> IF um.  The source annotation therefore subtracts this
    # translation, and 0.25-um QuPath pixels become 0.2125-um Xenium morphology pixels.
    transform = Transform(affine=[[1, 0, 10], [0, 1, 20]])
    doc, rec = IA.build_islet_geojson(
        source,
        transform,
        input_um_per_px=0.25,
        output_um_per_px=0.2125,
        width=2_000,
        height=2_000,
    )
    ring = np.asarray(doc["features"][0]["geometry"]["coordinates"][0])
    assert ring[0] == pytest.approx([(400 * 0.25 - 10) / 0.2125, (800 * 0.25 - 20) / 0.2125])
    assert np.ptp(ring, axis=0) == pytest.approx([240 * 0.25 / 0.2125] * 2)
    assert rec["source_um_per_px"] == 0.25
    assert rec["output_um_per_px"] == 0.2125
    assert rec["transform_operation_used"].startswith("Transform.inverse_apply")
    assert rec["transform_round_trip_max_um"] < 1e-9


def test_full_inverse_retains_a_displacement_field(tmp_path):
    source = _write_source(tmp_path / "source.geojson", [_feature(_polygon(20, 20, 60))])
    displacement = np.zeros((2, 200, 200), dtype=np.float32)
    displacement[0] = 2.0
    displacement[1] = -1.0
    transform = Transform(
        displacement=displacement,
        grid=GridSpec(0.0, 0.0, 1.0, 200, 200),
    )
    doc, rec = IA.build_islet_geojson(
        source,
        transform,
        input_um_per_px=1.0,
        output_um_per_px=1.0,
        width=200,
        height=200,
    )
    ring = np.asarray(doc["features"][0]["geometry"]["coordinates"][0])
    assert ring[0] == pytest.approx([18.0, 21.0])
    assert rec["includes_displacement"] is True
    assert rec["transform_round_trip_max_um"] < 1e-9


def test_an_inaccurate_inverse_is_refused_by_round_trip(tmp_path):
    source = _write_source(tmp_path / "source.geojson", [_feature(_polygon(20, 20, 60))])

    class BadInverse:
        displacement = None

        @staticmethod
        def inverse_apply(xy, *, n_iter):
            return np.asarray(xy)

        @staticmethod
        def apply(xy):
            return np.asarray(xy) + [1.0, 0.0]

    with pytest.raises(IA.IsletAnnotationError, match="full inverse transform misses"):
        IA.build_islet_geojson(
            source,
            BadInverse(),
            input_um_per_px=1.0,
            output_um_per_px=1.0,
            width=200,
            height=200,
        )


def test_size_filter_and_unsupported_source_geometry_are_explicit(tmp_path):
    source = _write_source(
        tmp_path / "source.geojson",
        [
            _feature(_polygon(100, 100, 240), identifier="large"),
            _feature(_polygon(20, 20, 10), identifier="fragment"),
            _feature(
                {"type": "LineString", "coordinates": [[0, 0], [10, 10]]},
                identifier="line",
            ),
            _feature(_polygon(0, 0, 300), cls="Tissue", identifier="tissue"),
        ],
    )
    doc, rec = IA.build_islet_geojson(
        source,
        Transform(),
        input_um_per_px=0.25,
        output_um_per_px=0.25,
        width=1_000,
        height=1_000,
    )
    assert len(doc["features"]) == 1
    assert rec["counts"] == {
        "total_source_features": 4,
        "source_classes": {"Islet": 3, "Tissue": 1},
        "islet_source_features": 3,
        "islet_polygon_parts": 2,
        "exported_islet_polygons": 1,
        "excluded_below_minimum_size": 1,
        "excluded_unsupported_geometry": 1,
        "unsupported_geometry_types": {"LineString": 1},
        "excluded_empty_geometry": 0,
        "excluded_non_islet_class": 1,
    }


def test_writer_emits_and_reopens_the_explorer_4_1_polygon_contract(tmp_path):
    source = _write_source(tmp_path / "source.geojson", [_feature(_polygon(50, 60, 240))])
    destination = tmp_path / "islets_on_xenium_grid.geojson"
    rec = IA.export_islet_geojson(
        source,
        destination,
        Transform(),
        input_um_per_px=0.25,
        output_um_per_px=0.25,
        width=1_000,
        height=1_000,
    )
    doc = json.loads(destination.read_text())
    feature = doc["features"][0]
    assert doc["type"] == "FeatureCollection"
    assert feature["geometry"]["type"] == "Polygon"
    assert len(feature["geometry"]["coordinates"]) == 1
    assert feature["geometry"]["coordinates"][0][0] == feature["geometry"]["coordinates"][0][-1]
    assert feature["properties"]["objectType"] == "annotation"
    assert feature["properties"]["classification"] == {"name": "Islet", "color": [255, 0, 0]}
    assert "colorRGB" not in feature["properties"]["classification"]
    assert rec["output_sha256"] == IA._sha256(destination)
    assert rec["validation"]["status"] == "PASS"


def test_source_holes_are_refused_instead_of_publishing_an_unimportable_polygon(tmp_path):
    polygon_with_hole = {
        "type": "Polygon",
        "coordinates": [
            [[0, 0], [300, 0], [300, 300], [0, 300], [0, 0]],
            [[100, 100], [200, 100], [200, 200], [100, 200], [100, 100]],
        ],
    }
    source = _write_source(tmp_path / "source.geojson", [_feature(polygon_with_hole)])
    with pytest.raises(IA.IsletAnnotationError, match="hole.*cannot import"):
        IA.build_islet_geojson(
            source,
            Transform(),
            input_um_per_px=0.25,
            output_um_per_px=0.25,
            width=1_000,
            height=1_000,
        )


@pytest.mark.parametrize(
    "geometry,message",
    [
        ({"type": "LineString", "coordinates": [[0, 0], [1, 1]]}, "expected 'Polygon'"),
        (
            {
                "type": "Polygon",
                "coordinates": [[[0, 0], [20, 20], [0, 20], [20, 0], [0, 0]]],
            },
            "simple valid polygon",
        ),
        ({"type": "Polygon", "coordinates": [[[0, 0], [2, 0], [2, 2], [0, 0]]]}, None),
    ],
)
def test_reopen_validation_refuses_unsupported_or_out_of_bounds_output(
    tmp_path, geometry, message
):
    feature = _feature(geometry)
    if message is None:
        feature["geometry"]["coordinates"][0][1][0] = 2_000
        message = "outside"
    output = _write_source(tmp_path / "bad.geojson", [feature])
    with pytest.raises(IA.IsletAnnotationError, match=message):
        IA.validate_islet_geojson(output, width=1_000, height=1_000)
