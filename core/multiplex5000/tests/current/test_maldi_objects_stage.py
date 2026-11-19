"""Stage08 serial-object contract exercised through the publishing task runner.

All acceptance records are made by the real gate API and bind synthetic input
receipts. These are test fixtures, not acceptance of any scientific acquisition.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import shapely

from multiplex5000 import gates, tasks
from multiplex5000.artifacts import sha256_file
from multiplex5000.integration.transform import Transform
from multiplex5000.stages import maldi_serial_integration as integration


def _csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _accept(folder: Path) -> None:
    """Real gate and receipt hashes, rebuilt only for intentional fixture edits."""
    gates.write(folder, gates.Gate(folder.name, "synthetic fixture", [
        gates.Check("synthetic fixture prepared", True, "true", True)
    ]), release="synthetic")
    outputs = [{"path": str(p.relative_to(folder)), "sha256": sha256_file(p)}
               for p in sorted(folder.rglob("*")) if p.is_file()
               and p.name not in {"receipt.json", "accepted.json"}]
    (folder / "receipt.json").write_text(json.dumps({
        "stage": folder.name, "release": "synthetic", "status": "complete",
        "outputs": outputs,
    }))
    gates.accept(folder, by="synthetic test fixture")


def _fixture(tmp_path: Path, *, platform: str = "phenocycler", with_line: bool = True):
    release = tmp_path / "release"
    tid = "synthetic-section"
    aid = ("phenocycler:exact-node.qptiff" if platform == "phenocycler"
           else "xenium:synthetic-slide:exact-region")
    fixed_frame, moving_frame = aid + ":um", f"maldi:{tid}:stage_um"
    rows = [{
        "tissue_id": tid, "donor_id": "synthetic-donor", "status": "registered",
        "neighbour_acquisition": aid, "neighbour_platform": platform,
        "neighbour_tissue": "pancreatic_lymph_node", "fixed_frame": fixed_frame,
        "moving_frame": moving_frame, "outline_rmse_um": 10.0, "reason": "",
    }]
    rows.append(dict(rows[0], tissue_id="unavailable-section", status="unavailable",
                     reason="synthetic missing outline"))
    registration = release / "07_registration"
    _csv(registration / "registration.csv", rows)
    transform_path = Transform(
        affine=np.array([[1.0, 0, 100.0], [0, 1.0, 200.0]]),
        meta={"fixed_acquisition_id": aid, "moving_acquisition_id": f"maldi:{tid}",
              "fixed_frame": fixed_frame, "moving_frame": moving_frame,
              "coordinate_units": "um", "fixed_geometry": {"pixel_um": 0.5}},
    ).save(registration / "transforms" / tid)
    _accept(registration)

    ingest = release / "04b_maldi_ingest"
    _csv(ingest / "features.csv", [{"feature": "glycan"}, {"feature": "glycogen"}])
    spots_path = ingest / f"data/spots/tissue_id={tid}/data.parquet"
    spots_path.parent.mkdir(parents=True)
    xy = np.array([[50, 50], [100, 100], [150, 150], [250, 150],
                   [350, 150], [500, 150], [100, 150], [150, 100]], float)
    pq.write_table(pa.table({
        "spot_uid": [f"maldi:{tid}:spot-{n}" for n in range(len(xy))],
        "x_um": xy[:, 0], "y_um": xy[:, 1],
        "on_tissue": [True] * 6 + [False, False],
        "glycan": [10.0, 0.0, np.nan, 20.0, 30.0, 40.0, 0.0, np.nan],
        "glycogen": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0, np.nan],
    }), spots_path)
    _accept(ingest)

    geometries = [("node-a", shapely.box(100, 200, 400, 500)),
                  ("node-b", shapely.box(300, 200, 500, 500)),
                  ("nonarea-line", shapely.LineString([(100, 200), (400, 500)]))][:3 if with_line else 2]
    export = tmp_path / "objects.geojson"
    export.write_text(json.dumps({"type": "FeatureCollection", "features": [
        {"type": "Feature", "id": identifier, "geometry": shapely.geometry.mapping(geometry),
         "properties": {"objectType": "annotation", "classification": {"name": "Node"}}}
        for identifier, geometry in geometries
    ]}))
    observation = tmp_path / "object-observation.json"
    observation.write_text(json.dumps({
        "schema": "multiplex5000-object-source-observation/1", "status": "complete",
        "acquisition_id": aid, "frame_id": fixed_frame, "geometry_kind": "annotation",
        "classes": ["Node"], "native_object_ids": [item[0] for item in geometries],
        "historical_producer": None,
        "outputs": [{"path": str(export), "sha256": sha256_file(export)}],
    }))
    manifest = tmp_path / "object-sources.json"
    manifest.write_text(json.dumps({"schema": integration.OBJECT_SCHEMA, "acquisitions": [{
        "acquisition_id": aid, "frame_id": fixed_frame, "sources": [{
            "binding_kind": "observed_export", "role": "regions", "classes": ["Node"],
            "units": "um", "path": str(export), "sha256": sha256_file(export),
            "source_receipt": str(observation), "source_receipt_sha256": sha256_file(observation),
        }],
    }]}))
    project = SimpleNamespace(
        release_dir=release, release="synthetic", path=tmp_path / "synthetic.yaml",
        settings={}, automatic_gates=(), previous=lambda stage: "07_registration",
        source=lambda key: manifest if key == "maldi_object_sources" else None,
    )
    return project, dict(tid=tid, acquisition=aid, xy=xy, transform=transform_path,
                         source=export, source_receipt=observation, manifest=manifest,
                         spots=spots_path)


@pytest.mark.parametrize("platform", ["phenocycler", "xenium"])
def test_stage08_publishes_serial_objects_with_finite_denominators(tmp_path, platform):
    project, fixture = _fixture(tmp_path, platform=platform)
    result = tasks.run_stage(integration, project)
    assert result.gate_passed and not result.accepted_automatically
    assert not (result.path / "gate/accepted.json").exists()
    assert (result.path / "gate/maldi_integration.png").is_file()

    relative = f"tissue_id={fixture['tid']}/data.parquet"
    spots = pq.read_table(result.path / "data/spots" / relative).to_pandas()
    assert len(spots) == 8 and spots.spot_uid.is_unique
    np.testing.assert_array_equal(spots[["x_neighbour_um", "y_neighbour_um"]],
                                  fixture["xy"] + [100, 200])
    assert spots.on_tissue.tolist() == [True] * 6 + [False, False]
    assert set(spots.relationship) == {"serial_section"}
    assert not any("cell" in name for name in spots.columns)
    assert spots.loc[3, "region_membership_state"] == "ambiguous_overlap"
    assert spots.loc[3, "region_object_id"] == ""
    assert spots.loc[0, "region_object_id"] == "node-a"
    assert spots.loc[4, "region_object_id"] == "node-b"
    assert spots.loc[5, "region_membership_state"] == "outside_or_unavailable"
    assert spots.signed_dist_islet_um.isna().all()
    assert set(spots.signed_dist_islet_status) == {"unavailable_no_islet_objects"}
    # The declared Node class holds a line that has no area: its cover is unavailable, never 0.
    assert spots["cover_fraction__regions__Node"].isna().all()
    assert set(spots.cover_status__regions) == {"unavailable_invalid_geometry"}
    assert set(spots.cover_status__islets) == {"unavailable_no_object_source"}
    assert "cover_fraction__islets__Islet" not in spots.columns
    assert set(spots.object_evidence_state) == {"partial_unavailable_objects"}

    objects = pq.read_table(result.path / "data/objects" / relative).to_pandas().set_index("object_id")
    assert set(objects.index) == {"node-a", "node-b", "nonarea-line"}
    a, b, line = (objects.loc[key] for key in ("node-a", "node-b", "nonarea-line"))
    # The overlapping spot belongs to both object summaries but has no unique spot label.
    assert a.n_spots == 4 and b.n_spots == 2
    assert a.evaluable and b.evaluable
    assert a.uncertainty_radius_um == 35  # recorded RMSE 10 + existing half-pitch 25
    assert a["finite_spots__glycan"] == 3 and a["detected_spots__glycan"] == 2
    assert a["detected_fraction__glycan"] == pytest.approx(2 / 3)
    assert a["median_detected__glycan"] == 15
    assert a["core_finite_spots__glycan"] == 3
    assert a["all_raster_n_spots"] == 6
    assert a["all_raster_finite_spots__glycan"] == 4
    assert a["all_raster_detected_fraction__glycan"] == 0.5
    assert a["all_raster_finite_spots__glycogen"] == 5
    assert a["all_raster_detected_fraction__glycogen"] == 0.8
    assert b["detected_fraction__glycan"] == 1
    assert not line.evaluable and "LineString" in line.reason
    assert np.isnan(line.area_um2) and np.isnan(line["detected_fraction__glycan"])

    receipt = json.loads((result.path / "receipt.json").read_text())
    summary = receipt["summary"]
    assert summary["sections"] == 1 and summary["declared_sections"] == 2
    assert summary["registration_unavailable"] == [
        {"tissue_id": "unavailable-section", "reason": "synthetic missing outline"}]
    assert summary["spots_on_tissue"] == 6
    assert summary["parameter_recording"]["status"] == "written"
    for record in receipt["outputs"]:
        assert sha256_file(result.path / record["path"]) == record["sha256"]
    parameters = json.loads((result.path / "resolved_parameters.json").read_text())
    assert parameters["context"]["relationship_type"] == "serial_section"
    assert parameters["automatic_application"] is False


def test_stage08_reports_each_spots_grid_square_cover(tmp_path):
    project, fixture = _fixture(tmp_path, with_line=False)
    result = tasks.run_stage(integration, project)
    spots = pq.read_table(result.path / f"data/spots/tissue_id={fixture['tid']}/data.parquet").to_pandas()
    # Squares carried by +(100, 200): spot 3's lies inside both nodes and is counted once; spot
    # 5's lies beyond them. Partial squares are tested on the helper itself.
    assert spots["cover_fraction__regions__Node"].tolist() == pytest.approx([1, 1, 1, 1, 1, 0, 1, 1])
    assert set(spots.cover_status__regions) == {"available"}
    sections = list(csv.DictReader((result.path / "sections.csv").open()))
    assert sections[0]["spots_touching_islets"] == ""


@pytest.mark.parametrize("changed", ["transform", "source", "source_receipt", "receipt", "gate", "spots"])
def test_stage08_refuses_stale_accepted_or_declared_content(tmp_path, changed):
    project, fixture = _fixture(tmp_path)
    path = (project.release_dir / "04b_maldi_ingest/receipt.json" if changed == "receipt" else
            project.release_dir / "04b_maldi_ingest/gate/gate.json" if changed == "gate" else
            fixture[changed])
    # Appended whitespace also changes a Parquet digest; no parser reads it before rejection.
    with path.open("ab") as handle:
        handle.write(b"\n")
    with pytest.raises(integration.MaldiIntegrationError, match="changed|acceptance does not bind"):
        tasks.run_stage(integration, project)
    assert not (project.release_dir / "08_integration/receipt.json").exists()


@pytest.mark.parametrize("field,value", [
    ("fixed_acquisition_id", "phenocycler:other-node.qptiff"),
    ("fixed_frame", "whole-slide-pixels"),
    ("moving_acquisition_id", "maldi:another-section"),
    ("coordinate_units", "pixel"),
])
def test_stage08_refuses_wrong_frame_even_when_receipt_hash_is_current(tmp_path, field, value):
    project, fixture = _fixture(tmp_path)
    record = json.loads(fixture["transform"].read_text())
    record["meta"][field] = value
    fixture["transform"].write_text(json.dumps(record))
    _accept(project.release_dir / "07_registration")
    with pytest.raises(integration.MaldiIntegrationError, match="exact acquisition/frame"):
        tasks.run_stage(integration, project)
    assert not (project.release_dir / "08_integration/receipt.json").exists()


def test_stage08_refuses_missing_frame_in_both_receipted_sources(tmp_path):
    project, fixture = _fixture(tmp_path)
    record = json.loads(fixture["transform"].read_text())
    del record["meta"]["moving_frame"]
    fixture["transform"].write_text(json.dumps(record))
    path = project.release_dir / "07_registration/registration.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        del row["moving_frame"]
    _csv(path, rows)
    _accept(project.release_dir / "07_registration")
    with pytest.raises(integration.MaldiIntegrationError, match="exact acquisition/frame"):
        tasks.run_stage(integration, project)


@pytest.mark.parametrize("document,field,value", [
    ("receipt", "stage", "wrong_stage"),
    ("receipt", "release", "wrong_release"),
    ("receipt", "status", "incomplete"),
    ("accepted", "stage", "wrong_stage"),
    ("accepted", "release", "wrong_release"),
    ("gate", "stage", "wrong_stage"),
    ("gate", "release", "wrong_release"),
    ("gate", "passed", False),
    ("gate", "checks", [{"name": "failed fixture", "blocking": True, "passed": False}]),
])
def test_stage08_refuses_incompatible_acceptance_metadata(tmp_path, document, field, value):
    project, _ = _fixture(tmp_path)
    folder = project.release_dir / "04b_maldi_ingest"
    paths = {"receipt": folder / "receipt.json", "gate": folder / "gate/gate.json",
             "accepted": folder / "gate/accepted.json"}
    records = {key: json.loads(path.read_text()) for key, path in paths.items()}
    records[document][field] = value
    for key in ("receipt", "gate"):
        paths[key].write_text(json.dumps(records[key]))
    # Hashes are current: this test rejects wrong semantics, independently of staleness.
    records["accepted"].update(receipt_sha256=sha256_file(paths["receipt"]),
                               gate_sha256=sha256_file(paths["gate"]))
    paths["accepted"].write_text(json.dumps(records["accepted"]))
    with pytest.raises(integration.MaldiIntegrationError, match="acceptance does not bind"):
        tasks.run_stage(integration, project)
    assert not (project.release_dir / "08_integration/receipt.json").exists()
