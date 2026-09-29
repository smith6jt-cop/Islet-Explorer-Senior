"""Focused tests for content-derived input and stage/run contracts."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

import phenocycler.artifacts as artifacts
from phenocycler.artifacts import (
    ChannelMapping,
    ContractError,
    DuplicateObjectError,
    GeometryContractError,
    QuPathCohortManifest,
    QuPathImageManifest,
    hash_identifiers,
    image_id_to_donor_id,
)


def _geometry(x: float = 0.0):
    return {
        "type": "Polygon",
        "coordinates": [
            [
                [x, 0.0],
                [x + 1.0, 0.0],
                [x + 1.0, 1.0],
                [x, 0.0],
            ]
        ],
    }


def _write_geojson(
    path: Path,
    object_ids=("a", "b"),
    *,
    include_nucleus: bool = True,
    missing_cell_index: int | None = None,
    missing_nucleus_index: int | None = None,
):
    features = []
    for index, object_id in enumerate(object_ids):
        feature = {
            "type": "Feature",
            "id": object_id,
            "geometry": None if index == missing_cell_index else _geometry(float(index)),
            "properties": {"objectType": "cell"},
        }
        if include_nucleus and index != missing_nucleus_index:
            feature["nucleusGeometry"] = _geometry(float(index) + 0.2)
        features.append(feature)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def _qupath_inputs(tmp_path: Path):
    csv_path = tmp_path / "cells.csv"
    csv_path.write_text(
        "Image,Object ID,Cell: DAPI: Mean,Cell: E-cadherin: Mean\n"
        "image-1,a,10,20\nimage-1,b,11,21\n"
    )
    qptiff = tmp_path / "image.qptiff"
    qptiff.write_bytes(b"small synthetic qptiff")
    geojson = tmp_path / "cells.geojson"
    _write_geojson(geojson)
    return csv_path, qptiff, geojson


def _qupath_manifest(tmp_path: Path) -> QuPathImageManifest:
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    return QuPathImageManifest.create(
        donor_id="6374",
        image_id="6374_Scan1",
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=0.5,
        pixel_size_um_y=0.5,
        panel_id="panel-v1",
        channel_map=(
            ChannelMapping("DAPI", "DAPI", 0),
            ChannelMapping("E_cadherin", "E-cadherin", 1),
        ),
        segmentation_version="instanseg-0.1.1+script-v2",
    )


def _cohort_image(
    tmp_path: Path,
    *,
    donor_id: str,
    image_id: str,
    panel_id: str = "panel-v1",
    channel_map=(
        ChannelMapping("DAPI", "DAPI", 0),
        ChannelMapping("E_cadherin", "E-cadherin", 1),
    ),
) -> QuPathImageManifest:
    image_dir = tmp_path / f"{donor_id}_{image_id.replace('/', '_')}"
    image_dir.mkdir(parents=True, exist_ok=True)
    csv_path, qptiff, geojson = _qupath_inputs(image_dir)
    return QuPathImageManifest.create(
        donor_id=donor_id,
        image_id=image_id,
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=0.5,
        pixel_size_um_y=0.5,
        panel_id=panel_id,
        channel_map=channel_map,
        segmentation_version="instanseg-0.1.1+script-v2",
    )


def test_qupath_manifest_roundtrip_and_current_validation(tmp_path):
    manifest = _qupath_manifest(tmp_path)
    path = manifest.write_json(tmp_path / "qupath_manifest.json")
    loaded = QuPathImageManifest.read_json(path)

    assert loaded == manifest
    assert loaded.cell_geometry.feature_count == 2
    assert loaded.cell_geometry.nucleus_geometry_count == 2
    assert loaded.cell_geometry.object_id_sha256 == hash_identifiers(["a", "b"])


def test_qptiff_can_use_an_explicitly_labelled_sampled_fingerprint(tmp_path):
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    manifest = QuPathImageManifest.create(
        donor_id="1",
        image_id="image",
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=1.0,
        pixel_size_um_y=1.0,
        panel_id="panel",
        channel_map=(ChannelMapping("DAPI", "DAPI", 0),),
        segmentation_version="seg-v1",
        qptiff_hash_mode="sampled",
        qptiff_sample_size_bytes=4,
        qptiff_sample_count=3,
    )
    assert manifest.qptiff.hash_mode == "sampled"
    assert manifest.qptiff.sample_size_bytes == 4


def test_qupath_manifest_fails_on_missing_or_duplicate_geometry(tmp_path):
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    _write_geojson(geojson, include_nucleus=False)
    kwargs = dict(
        donor_id="1",
        image_id="image",
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=1.0,
        pixel_size_um_y=1.0,
        panel_id="panel",
        channel_map=(ChannelMapping("DAPI", "DAPI", 0),),
        segmentation_version="seg-v1",
    )
    with pytest.raises(GeometryContractError, match="nucleusGeometry"):
        QuPathImageManifest.create(**kwargs)

    _write_geojson(geojson, object_ids=("same", "same"))
    with pytest.raises(DuplicateObjectError, match="duplicate object_id"):
        QuPathImageManifest.create(**kwargs)

    _write_geojson(geojson, missing_cell_index=1)
    with pytest.raises(GeometryContractError, match="no usable cell geometry"):
        QuPathImageManifest.create(**kwargs)


def test_qupath_manifest_records_partially_missing_embedded_nuclei(tmp_path):
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    _write_geojson(geojson, missing_nucleus_index=0)
    manifest = QuPathImageManifest.create(
        donor_id="1",
        image_id="image",
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=1.0,
        pixel_size_um_y=1.0,
        panel_id="panel",
        channel_map=(ChannelMapping("DAPI", "DAPI", 0),),
        segmentation_version="seg-v1",
    )

    assert manifest.cell_geometry.feature_count == 2
    assert manifest.cell_geometry.nucleus_geometry_count == 1
    assert manifest.nucleus_geometry.primary_geometry_count == 1


def test_channel_map_and_measurement_contracts_are_explicit(tmp_path):
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    with pytest.raises(ContractError, match="duplicate markers"):
        QuPathImageManifest.create(
            donor_id="1",
            image_id="image",
            measurement_csv=csv_path,
            qptiff=qptiff,
            cell_geojson=geojson,
            nucleus_geojson=geojson,
            pixel_size_um_x=1.0,
            pixel_size_um_y=1.0,
            panel_id="panel",
            channel_map=(
                ChannelMapping("DAPI", "DAPI", 0),
                ChannelMapping("DAPI", "DAPI-2", 1),
            ),
            segmentation_version="seg-v1",
        )

    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("Image,Not Object ID\nimage,a\n")
    with pytest.raises(ContractError, match="Object ID"):
        QuPathImageManifest.create(
            donor_id="1",
            image_id="image",
            measurement_csv=bad_csv,
            qptiff=qptiff,
            cell_geojson=geojson,
            nucleus_geojson=geojson,
            pixel_size_um_x=1.0,
            pixel_size_um_y=1.0,
            panel_id="panel",
            channel_map=(ChannelMapping("DAPI", "DAPI", 0),),
            segmentation_version="seg-v1",
        )


def test_qupath_cohort_manifest_roundtrip_and_pure_ingest_mapping(tmp_path):
    images = (
        _cohort_image(tmp_path, donor_id="1", image_id="donor-1-image"),
        _cohort_image(tmp_path, donor_id="2", image_id="donor-2-image"),
    )
    cohort = QuPathCohortManifest.create(
        cohort_name="synthetic-cohort",
        expected_donors=("1", "2"),
        images=images,
    )
    path = cohort.write_json(tmp_path / "qupath_cohort.json")
    loaded = QuPathCohortManifest.read_json(path)

    assert loaded == cohort
    assert loaded.donor_ids == ("1", "2")
    assert loaded.image_ids == ("donor-1-image", "donor-2-image")
    assert len(loaded.panels) == 1
    assert loaded.panels[0].panel_id == "panel-v1"
    expected = {"donor-1-image": "1", "donor-2-image": "2"}
    assert image_id_to_donor_id(loaded) == expected
    assert image_id_to_donor_id(loaded.images) == expected
    assert loaded.image_to_donor() == expected


def test_qupath_cohort_requires_exactly_one_unique_image_per_expected_donor(tmp_path):
    first = _cohort_image(tmp_path, donor_id="1", image_id="image-1")
    second = _cohort_image(tmp_path, donor_id="2", image_id="image-2")
    with pytest.raises(ContractError, match="missing=.*3"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1", "2", "3"),
            images=(first, second),
        )
    with pytest.raises(ContractError, match="extra=.*2"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1",),
            images=(first, second),
        )

    duplicate_donor = _cohort_image(
        tmp_path, donor_id="1", image_id="another-image"
    )
    with pytest.raises(ContractError, match="duplicate donor_id"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1",),
            images=(first, duplicate_donor),
        )

    duplicate_image = _cohort_image(
        tmp_path, donor_id="2", image_id="image-1"
    )
    with pytest.raises(ContractError, match="duplicate image_id"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1", "2"),
            images=(first, duplicate_image),
        )


def test_qupath_cohort_enforces_one_exact_channel_map_per_panel(tmp_path):
    first = _cohort_image(tmp_path, donor_id="1", image_id="image-1")
    conflicting = _cohort_image(
        tmp_path,
        donor_id="2",
        image_id="image-2",
        panel_id="panel-v1",
        channel_map=(
            ChannelMapping("DAPI", "DAPI", 0),
            ChannelMapping("E_cadherin", "E-cadherin", 2),
        ),
    )
    with pytest.raises(ContractError, match="contradictory channel mappings"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1", "2"),
            images=(first, conflicting),
        )

    # A physical panel change is valid only when it has a distinct panel id.
    second_panel = _cohort_image(
        tmp_path,
        donor_id="2",
        image_id="image-2-panel2",
        panel_id="panel-v2",
        channel_map=(
            ChannelMapping("DAPI", "DAPI", 0),
            ChannelMapping("CD38", "CD38", 1),
        ),
    )
    cohort = QuPathCohortManifest.create(
        cohort_name="cohort",
        expected_donors=("1", "2"),
        images=(first, second_panel),
    )
    assert [panel.panel_id for panel in cohort.panels] == ["panel-v1", "panel-v2"]


