"""Contract tests for the PhenoCycler -> Xenium-grid annotation sidecars.

Synthetic geometry only: the point is the coordinate path, the exclusions, and the Explorer
contract, none of which need a 33 GB qptiff to pin down.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2]))

import islet_annotations as IA  # noqa: E402
import shapely  # noqa: E402
import pheno_annotations as PA  # noqa: E402

from phenocycler.integration.rasterize import GridSpec  # noqa: E402
from phenocycler.integration.transform import Transform  # noqa: E402

PX_UM = 0.5  # source qptiff calibration
OUT_UM = 0.2125  # Xenium morphology calibration
WIDTH, HEIGHT = 12000, 15000  # 2550 x 3187 um at OUT_UM

ISLET = PA.CLASSES_BY_NAME["Islet"]


def _grid() -> GridSpec:
    return GridSpec(origin_x_um=0.0, origin_y_um=0.0, pixel_um=2.0, height=600, width=600)


def _identity() -> Transform:
    return Transform(affine=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]), grid=_grid())


def _shifted(dx_um: float, dy_um: float) -> Transform:
    """Xenium microns -> PhenoCycler microns, the direction S4 stores."""
    return Transform(affine=np.array([[1.0, 0.0, dx_um], [0.0, 1.0, dy_um]]), grid=_grid())


def _square(cx_px: float, cy_px: float, half_px: float, cls: str = "Islet") -> dict:
    ring = [
        [cx_px - half_px, cy_px - half_px],
        [cx_px + half_px, cy_px - half_px],
        [cx_px + half_px, cy_px + half_px],
        [cx_px - half_px, cy_px + half_px],
        [cx_px - half_px, cy_px - half_px],
    ]
    return {
        "type": "Feature",
        "geometry": {"type": "Polygon", "coordinates": [ring]},
        "properties": {"objectType": "annotation", "classification": {"name": cls}},
    }


def _write(tmp_path: Path, features: list[dict], name: str = "src.geojson") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))
    return path


def test_source_pixels_become_xenium_grid_pixels_through_the_inverse(tmp_path):
    """A known shift must appear as exactly that shift, in Xenium pixels, and no other."""
    source = _write(tmp_path, [_square(1000.0, 1200.0, 200.0)])
    polygons, _ = PA._exterior_polygons(source, "Islet")
    records, audit = PA.place_on_xenium_grid(
        polygons, _shifted(300.0, -100.0),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
    )
    assert audit["exported_polygons"] == 1
    ring = records[0]["ring_px"]
    # 800 px * 0.5 um/px = 400 um in PhenoCycler; minus the 300 um shift = 100 um in Xenium.
    expected_x = ((1000.0 - 200.0) * PX_UM - 300.0) / OUT_UM
    expected_y = ((1200.0 - 200.0) * PX_UM + 100.0) / OUT_UM
    assert ring[:, 0].min() == pytest.approx(expected_x)
    assert ring[:, 1].min() == pytest.approx(expected_y)


def test_annotations_outside_the_xenium_field_are_dropped_and_counted(tmp_path):
    """The Xenium section is a portion of the PhenoCycler one; the rest has nowhere to land."""
    inside = _square(1000.0, 1000.0, 100.0)
    outside = _square(50_000.0, 50_000.0, 100.0)  # far beyond the Xenium grid
    source = _write(tmp_path, [inside, outside])
    polygons, _ = PA._exterior_polygons(source, "Islet")
    records, audit = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
    )
    assert len(records) == 1
    assert audit["excluded_outside_xenium_field"] == 1
    assert audit["fraction_of_parts_inside_xenium_field"] == pytest.approx(0.5)


def test_a_partially_overlapping_annotation_is_dropped_not_clipped(tmp_path):
    """Clipping to the grid edge would invent an anatomical boundary that is not there."""
    straddling = _square(20.0, 1000.0, 200.0)  # crosses x = 0
    source = _write(tmp_path, [straddling])
    polygons, _ = PA._exterior_polygons(source, "Islet")
    records, audit = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
    )
    assert records == []
    assert audit["excluded_outside_xenium_field"] == 1


def test_a_non_convergent_inverse_excludes_only_the_affected_polygon(tmp_path):
    """A displacement field whose inverse misses must lose that object, not the whole class."""
    source = _write(tmp_path, [_square(600.0, 600.0, 100.0), _square(1000.0, 1000.0, 100.0)])
    polygons, _ = PA._exterior_polygons(source, "Islet")
    records, audit = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
        round_trip_limit_um=-1.0,  # nothing can satisfy this
    )
    assert records == []
    assert audit["excluded_inverse_did_not_converge"] == 2


def test_size_filter_is_off_by_default_and_recorded_when_used(tmp_path):
    big = _square(1000.0, 1000.0, 200.0)      # 200 um across at 0.5 um/px
    small = _square(2000.0, 2000.0, 5.0)      # 5 um across
    source = _write(tmp_path, [big, small])
    polygons, _ = PA._exterior_polygons(source, "Islet")

    kept_all, audit_all = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
    )
    assert len(kept_all) == 2
    assert audit_all["excluded_below_minimum_size"] == 0

    kept_filtered, audit_filtered = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
        min_diameter_um=50.0,
    )
    assert len(kept_filtered) == 1
    assert audit_filtered["excluded_below_minimum_size"] == 1
    assert audit_filtered["minimum_equivalent_diameter_um"] == 50.0


def test_a_class_over_the_cap_is_split_largest_first_and_nothing_is_lost(tmp_path):
    records = [
        {"ring_px": np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 0.0]]),
         "area_um2": float(i), "diameter_um": 1.0, "round_trip_um": 0.0, "source_part_index": i}
        for i in range(25)
    ]
    documents = PA.build_documents(records, ISLET, max_features=10)
    assert [d["n_features"] for d in documents] == [10, 10, 5]
    assert all(d["n_parts"] == 3 for d in documents)
    assert sum(d["n_features"] for d in documents) == 25
    # Largest object first, and ids continue across parts rather than restarting.
    assert documents[0]["area_um2"] > documents[-1]["area_um2"]
    first_ids = [f["id"] for f in documents[0]["document"]["features"]]
    last_ids = [f["id"] for f in documents[-1]["document"]["features"]]
    assert first_ids[0] == "Islet-1"
    assert last_ids[-1] == "Islet-25"
    assert not set(first_ids) & set(last_ids)


def test_written_sidecars_reopen_under_the_explorer_contract_for_their_own_class(tmp_path):
    features = [_square(500.0 + 400 * i, 700.0, 150.0, cls="Nerve") for i in range(4)]
    source = _write(tmp_path, features)
    nerve = PA.CLASSES_BY_NAME["Nerve"]
    result = PA.export_class(
        source, tmp_path / "out", _identity(), nerve,
        serial="0059947", input_um_per_px=PX_UM, output_um_per_px=OUT_UM,
        width=WIDTH, height=HEIGHT, max_features=3,
    )
    assert result["n_features_total"] == 4
    assert result["n_files"] == 2
    names = [f["file"] for f in result["files"]]
    assert names == [
        "0059947_pheno_Nerve_on_xenium_grid_part01of02.geojson",
        "0059947_pheno_Nerve_on_xenium_grid_part02of02.geojson",
    ]
    for entry in result["files"]:
        assert entry["validation"]["status"] == "PASS"
        assert entry["validation"]["checks"]["classification"] == "Nerve"
        assert entry["validation"]["checks"]["classification_color_rgb"] == list(nerve.color)
    assert result["transform_direction_stored"] == "Xenium microns -> PhenoCycler microns"
    assert result["section_relationship"] == "serial"
    assert result["includes_displacement"] is False


def test_the_islet_lane_contract_is_unchanged_by_the_shared_validator(tmp_path):
    """The same-slide sidecar still validates with its own class and colour by default."""
    document = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "id": "Islet-1",
                "geometry": {"type": "Polygon",
                             "coordinates": [[[1.0, 1.0], [9.0, 1.0], [9.0, 9.0], [1.0, 1.0]]]},
                "properties": {"objectType": "annotation", "name": "Islet 1",
                               "classification": {"name": "Islet", "color": [255, 0, 0]}},
            }
        ],
    }
    path = tmp_path / "islets.geojson"
    path.write_text(json.dumps(document))
    report = IA.validate_islet_geojson(path, width=100, height=100, expected_features=1)
    assert report["status"] == "PASS"
    assert report["checks"]["classification"] == "Islet"
    assert report["checks"]["classification_color_rgb"] == [255, 0, 0]

    # And a different class is refused against the islet default rather than waved through.
    with pytest.raises(IA.IsletAnnotationError, match="classification.name"):
        IA.validate_islet_geojson(path, width=100, height=100, class_name="Nerve")


def test_two_classes_in_one_file_are_separated(tmp_path):
    """`Islets` carries both the core and its 20 um expansion; they publish separately."""
    source = _write(
        tmp_path,
        [_square(1000.0, 1000.0, 200.0, cls="Islet"),
         _square(1000.0, 1000.0, 260.0, cls="IsletExpanded")],
    )
    core, _ = PA._exterior_polygons(source, "Islet")
    expanded, audit = PA._exterior_polygons(source, "IsletExpanded")
    assert len(core) == 1 and len(expanded) == 1
    assert audit["source_classes"] == {"Islet": 1, "IsletExpanded": 1}
    assert audit["matched_source_features"] == 1


def test_multipolygon_parts_are_exploded_rather_than_refused(tmp_path):
    multi = {
        "type": "Feature",
        "geometry": {
            "type": "MultiPolygon",
            "coordinates": [
                [[[100.0, 100.0], [300.0, 100.0], [300.0, 300.0], [100.0, 100.0]]],
                [[[900.0, 900.0], [1100.0, 900.0], [1100.0, 1100.0], [900.0, 900.0]]],
            ],
        },
        "properties": {"objectType": "annotation", "classification": {"name": "Tissue"}},
    }
    source = _write(tmp_path, [multi])
    polygons, audit = PA._exterior_polygons(source, "Tissue")
    assert len(polygons) == 2
    assert audit["polygon_parts"] == 2
    assert audit["matched_source_features"] == 1


# --------------------------------------------------------------------------------------
# The Xenium tissue outline
# --------------------------------------------------------------------------------------

def _outline_um(x0: float, y0: float, x1: float, y1: float):
    """A rectangular tissue outline in Xenium microns."""
    import shapely

    return shapely.Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])


def test_the_outline_is_the_region_annotations_are_selected_against(tmp_path):
    """Inside the morphology rectangle is not the same thing as on tissue.

    Both squares clear the grid bounds. Only one is on tissue, and the old rectangle test
    could not tell them apart -- which is why annotations reached Explorer sitting on blank
    background.
    """
    # Source pixels are 0.5 um, so these are 450-550 um and 1950-2050 um in Xenium microns.
    # The grid spans 0-2550 x 0-3187 um, so BOTH clear the rectangle; only one is on tissue.
    on_tissue = _square(1000.0, 1000.0, 100.0)
    off_tissue = _square(4000.0, 4000.0, 100.0)
    source = _write(tmp_path, [on_tissue, off_tissue])
    polygons, _ = PA._exterior_polygons(source, "Islet")

    records, audit = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
        tissue_outline=_outline_um(0.0, 0.0, 1000.0, 1000.0),
    )
    assert len(records) == 1
    assert audit["excluded_outside_xenium_field"] == 0   # both cleared the rectangle
    assert audit["excluded_outside_tissue_outline"] == 1
    assert audit["fraction_of_parts_inside_tissue_outline"] == pytest.approx(0.5)
    assert audit["tissue_outline_applied"] is True


def test_an_annotation_straddling_the_tissue_edge_is_dropped_not_clipped(tmp_path):
    """Cutting at the tissue edge would invent a boundary and misstate the annotation's area.

    The same rule the grid edge already follows, applied to the region that actually matters.
    """
    straddling = _square(2000.0, 2000.0, 500.0)  # 750-1250 um, crossing the edge at 1000 um
    source = _write(tmp_path, [straddling])
    polygons, _ = PA._exterior_polygons(source, "Islet")

    records, audit = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
        tissue_outline=_outline_um(0.0, 0.0, 1000.0, 1000.0),
    )
    assert audit["excluded_outside_xenium_field"] == 0   # it cleared the rectangle
    assert records == []
    assert audit["excluded_outside_tissue_outline"] == 1


def test_without_an_outline_the_behaviour_is_exactly_as_before(tmp_path):
    """The outline is opt-in; omitting it must not silently change what a caller publishes."""
    source = _write(tmp_path, [_square(1000.0, 1000.0, 100.0)])
    polygons, _ = PA._exterior_polygons(source, "Islet")
    records, audit = PA.place_on_xenium_grid(
        polygons, _identity(),
        input_um_per_px=PX_UM, output_um_per_px=OUT_UM, width=WIDTH, height=HEIGHT,
    )
    assert len(records) == 1
    assert audit["excluded_outside_tissue_outline"] == 0
    assert audit["tissue_outline_applied"] is False
    assert audit["fraction_of_parts_inside_tissue_outline"] is None


def test_the_outline_comes_from_the_cell_field_not_the_capture_rectangle():
    """A section occupying part of its capture area must yield an outline, not the rectangle."""
    rng = np.random.default_rng(0)
    # One 2 x 2 mm blob of cells inside a 6 x 6 mm capture area.
    xy = rng.uniform([1000.0, 1000.0], [3000.0, 3000.0], size=(40_000, 2))
    geom, prov = PA.xenium_tissue_outline(xy, pixel_um=20.0)

    assert prov["source"] == "xenium cell centroids"
    assert prov["n_cells"] == 40_000
    # ~4 mm^2 of tissue, and emphatically not the 36 mm^2 a capture rectangle would claim.
    assert 3.0e6 < geom.area < 5.0e6
    assert geom.contains(shapely.Point(2000.0, 2000.0))
    assert not geom.contains(shapely.Point(5000.0, 5000.0))
