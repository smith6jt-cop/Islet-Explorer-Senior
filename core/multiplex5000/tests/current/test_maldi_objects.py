"""Object provenance, uncertainty, overlap and finite detection denominators."""
from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest
import shapely

from multiplex5000.integration.transform import Transform
from multiplex5000.stages import maldi_serial_integration as M
from multiplex5000.tasks import Input, Task, sha256_file


def _feature(identifier, geometry, *, kind="annotation"):
    return {"type": "Feature", "id": identifier, "geometry": shapely.geometry.mapping(geometry),
            "properties": {"objectType": kind, "classification": {"name": "Node"}}}


def _declared(tmp_path):
    aid = "phenocycler:exact-image.qptiff"
    export = tmp_path / "objects.geojson"
    export.write_text(json.dumps({"features": [_feature("polygon", shapely.box(0, 0, 200, 200)),
                                              _feature("line", shapely.LineString([(0, 0), (5, 5)]))]}))
    receipt = tmp_path / "source.json"
    receipt.write_text(json.dumps({"schema": "multiplex5000-object-source-observation/1",
        "status": "complete", "acquisition_id": aid, "frame_id": aid + ":um",
        "geometry_kind": "annotation", "classes": ["Node"], "frame_pixel_size_um": 0.5,
        "native_object_ids": ["polygon", "line"], "historical_producer": None,
        "outputs": [{"path": str(export), "sha256": sha256_file(export)}]}))
    entry = {"acquisition_id": aid, "frame_id": aid + ":um", "sources": [{
        "binding_kind": "observed_export", "role": "regions", "classes": ["Node"],
        "units": "pixel", "pixel_size_um": 0.5, "path": str(export), "sha256": sha256_file(export),
        "source_receipt": str(receipt), "source_receipt_sha256": sha256_file(receipt)}]}
    manifest = tmp_path / "sources.json"
    manifest.write_text(json.dumps({"schema": M.OBJECT_SCHEMA, "acquisitions": [entry]}))
    project = SimpleNamespace(source=lambda key: manifest if key == "maldi_object_sources" else None)
    row = {"tissue_id": "sample", "neighbour_acquisition": aid, "fixed_frame": aid + ":um"}
    declared = M._object_inputs(project, row)
    task = Task("08_integration", "candidate", tmp_path / "out", {i.name: i for i in declared}, project=project)
    transform = Transform(affine=np.array([[1.0, 0, 0], [0, 1.0, 0]]))
    transform.meta.update(fixed_acquisition_id=aid, fixed_frame=aid + ":um", fixed_geometry={"pixel_um": 0.5})
    return task, row, transform, entry


def test_exact_source_retains_nonarea_objects_as_unavailable(tmp_path):
    task, row, transform, entry = _declared(tmp_path)
    loaded = M._declared_objects(task, row, transform, entry)
    assert not loaded["islets"] and not loaded["structures"]
    polygon, line = loaded["regions"]
    assert polygon["geometry"].area == 10000
    assert line["geometry"] is None and "LineString" in line["unavailable_reason"]
    rows = M._object_table(loaded["regions"], np.array([[50.0, 50.0]]), np.array([[0.0]]),
                           ["feature"], 25.0, {})
    assert rows[0]["evaluable"] and rows[0]["detected_fraction__feature"] == 0
    assert not rows[1]["evaluable"] and np.isnan(rows[1]["area_um2"])
    assert np.isnan(rows[1]["detected_fraction__feature"])


@pytest.mark.parametrize("change", ["acquisition", "frame", "calibration", "source", "receipt"])
def test_wrong_frame_or_stale_source_rejected(tmp_path, change):
    task, row, transform, entry = _declared(tmp_path)
    if change == "acquisition":
        transform.meta["fixed_acquisition_id"] = "phenocycler:other-image.qptiff"
    elif change == "frame":
        transform.meta["fixed_frame"] = "whole-slide-pixels"
    elif change == "calibration":
        transform.meta["fixed_geometry"]["pixel_um"] = 0.51
    else:
        path = task.input("object_source:sample:0" if change == "source" else "object_receipt:sample:0")
        path.write_text(path.read_text() + " ")
    with pytest.raises(M.MaldiIntegrationError):
        M._declared_objects(task, row, transform, entry)


def test_cell_geometry_is_never_promoted_to_serial_object(tmp_path):
    path = tmp_path / "cells.json"
    path.write_text(json.dumps({"features": [_feature("cell", shapely.box(0, 0, 10, 10), kind="cell")]}))
    with pytest.raises(M.MaldiIntegrationError, match="cannot be cells"):
        M._objects(path, 1.0, strict=True)


def test_duplicate_native_object_ids_are_rejected(tmp_path):
    path = tmp_path / "duplicate.json"
    feature = _feature("same", shapely.box(0, 0, 10, 10))
    path.write_text(json.dumps({"features": [feature, feature]}))
    with pytest.raises(M.MaldiIntegrationError, match="duplicate native"):
        M._objects(path, 1.0, strict=True)


def test_overlap_is_explicit_and_independent_of_input_order():
    objects = [{"geometry": shapely.box(0, 0, 100, 100)}, {"geometry": shapely.box(50, 0, 150, 100)}]
    points = shapely.points([[25, 50], [75, 50], [125, 50], [200, 50]])
    np.testing.assert_array_equal(M._containing(points, objects), [0, -2, 1, -1])
    np.testing.assert_array_equal(M._containing(points, objects[::-1]), [1, -2, 0, -1])


def test_finite_denominators_and_core_summary_do_not_turn_missing_into_zero():
    objects = [{"object_id": "a", "object_class": "Node", "object_name": "", "size_class": "",
                "geometry": shapely.box(0, 0, 400, 400)}]
    xy = np.array([[10, 10], [100, 100], [200, 200], [300, 300]], float)
    values = np.array([[10.0], [0.0], [np.nan], [20.0]])
    result = M._object_table(objects, xy, values, ["feature"], 25.0, {})[0]
    assert result["finite_spots__feature"] == 3 and result["detected_spots__feature"] == 2
    assert result["detected_fraction__feature"] == 2 / 3
    assert result["core_finite_spots__feature"] == 2
    assert result["core_detected_fraction__feature"] == 0.5


def test_missing_previously_bound_displacement_is_not_ignored(tmp_path):
    root = tmp_path / "release"
    path = root / "07_registration/transforms/sample/transform_displacement.npz"
    task = Task("08_integration", "candidate", tmp_path / "out",
                {"displacement:sample": Input("displacement:sample", path, required=False)},
                project=SimpleNamespace(release_dir=root))
    receipt = {"outputs": [{"path": "transforms/sample/transform_displacement.npz", "sha256": "0" * 64}]}
    with pytest.raises(M.MaldiIntegrationError, match="missing"):
        M._bound_output(task, receipt, "07_registration", "displacement:sample")
