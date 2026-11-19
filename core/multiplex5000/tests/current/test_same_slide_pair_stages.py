"""A synthetic same-section pilot runs registration review before candidate links."""

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from multiplex5000 import gates, tasks
from multiplex5000.artifacts import sha256_file
from multiplex5000.integration.transform import Transform
from multiplex5000.same_slide_sources import SameSlidePair
from multiplex5000.stages import same_slide_integration as integration
from multiplex5000.stages import same_slide_registration as registration


def _table(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(pd.DataFrame(payload), preserve_index=False), path)


def _pilot(tmp_path, monkeypatch):
    release = tmp_path / "release"
    entry, acquisition, sample = "HDL090_PC67_0110618-S1", "xenium:0110618:090", "0110618__090"
    source = tmp_path / "output-XETG00298__0110618__090__20260826__224731"
    pair = SameSlidePair("spleen", "HDL090", entry, tmp_path / f"{entry}.csv.gz",
                         tmp_path / "0110618.qptiff", acquisition, source,
                         "same_section_after_phenocycler", "crop_local_um", "region_local_um")
    manifest = tmp_path / "same_slide.yaml"
    manifest.write_text("schema: test\n")
    monkeypatch.setattr(registration, "load_pairs", lambda path: (pair,))
    folder = release / "06_typing"
    folder.mkdir(parents=True)
    (folder / "receipt.json").write_text(json.dumps({"stage": "06_typing",
                                                     "status": "complete",
                                                     "inputs": [{
                                                         "name": "same_slide_source_manifest",
                                                         "sha256": sha256_file(manifest)}]}))
    gate = folder / "gate"
    gate.mkdir()
    (gate / "gate.json").write_text(json.dumps({"stage": "06_typing", "release": "pilot",
                                               "checks": [], "passed": True}))
    gates.accept(folder, by="test")
    regions = release / "04_rna_ingest/regions.csv"
    regions.parent.mkdir(parents=True)
    with regions.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["acquisition_id", "sample_key",
                                                   "donor_id", "tissue", "source_path"])
        writer.writeheader()
        writer.writerow({"acquisition_id": acquisition, "sample_key": sample,
                         "donor_id": "HDL090", "tissue": "spleen", "source_path": source})
    p_uid = [f"phenocycler:0110618:{entry}:digest:{i}" for i in (1, 2, 3)]
    x_uid = [f"{acquisition}:cells:{i}" for i in (1, 2, 3)]
    _table(release / f"03_protein_measurement/cells/{entry}.parquet", {
        "cell_uid": p_uid, "object_id": ["1", "2", "3"], "donor_id": ["HDL090"] * 3,
        "tissue": ["spleen"] * 3, "acquisition_id": ["phenocycler:0110618"] * 3,
        "x_um": [100.0, 110.0, 120.0], "y_um": [0.0] * 3,
        "protein_eligible": [True, True, False]})
    _table(release / f"04_rna_ingest/cells/{sample}.parquet", {
        "cell_uid": x_uid, "donor_id": ["HDL090"] * 3,
        "tissue": ["spleen"] * 3, "acquisition_id": [acquisition] * 3,
        "x_um": [0.5, 10.5, 20.5], "y_um": [0.0] * 3,
        "typing_eligible": [True, True, False]})
    _table(release / f"06_typing/protein/{entry}.parquet", {
        "cell_uid": p_uid, "object_id": ["1", "2", "3"],
        "candidate_phenotype": ["T cell", "", ""],
        "candidate_state": ["unique_candidate", "unavailable", "unavailable"]})
    _table(release / f"06_typing/rna/{sample}.parquet", {
        "cell_uid": x_uid, "typing_eligible": [True, True, False],
        "gene:CD3E": [1, -1, -1]})
    receipt = json.loads((folder / "receipt.json").read_text())
    receipt["inputs"] += [
        {"name": f"protein_cells:{entry}", "sha256": sha256_file(
            release / f"03_protein_measurement/cells/{entry}.parquet")},
        {"name": f"rna_cells:{sample}", "sha256": sha256_file(
            release / f"04_rna_ingest/cells/{sample}.parquet")}]
    receipt["outputs"] = [{"path": str(path.relative_to(folder)), "sha256": sha256_file(path)}
                          for path in sorted(folder.glob("*/*.parquet"))]
    (folder / "receipt.json").write_text(json.dumps(receipt))
    gates.accept(folder, by="test")
    transforms = tmp_path / "transforms"
    transform = Transform(
        affine=np.array([[1, 0, 100], [0, 1, 0]]),
        meta={"protein_entry": entry, "xenium_acquisition_id": acquisition,
              "moving_frame": "region_local_um", "fixed_frame": "crop_local_um",
              "reviewed": True})
    transform.save(transforms / entry)
    project = SimpleNamespace(
        path=tmp_path / "pilot.yaml", release="pilot", release_dir=release,
        settings={"same_slide_source_manifest": str(manifest),
                  "rna_acquisitions": [acquisition]}, automatic_gates=(),
        source=lambda name: {"same_slide_sources": manifest,
                             "same_slide_registrations": transforms}.get(name),
        previous=lambda stage: {"07_registration": "06_typing",
                                "08_integration": "07_registration"}[stage])
    return project, pair


def test_registration_then_candidate_link_stage(tmp_path, monkeypatch):
    project, pair = _pilot(tmp_path, monkeypatch)
    registered = tasks.run_stage(registration, project)
    assert registered.gate_passed
    assert (registered.path / "gate/registration_overlay.png").is_file()
    assert (registered.path / "gate/candidate_distance_map.png").is_file()
    assert not (registered.path / "gate/accepted.json").exists()
    with pytest.raises(gates.GateError, match="not been accepted"):
        tasks.run_stage(integration, project)
    gates.accept(registered.path, by="test")
    result = tasks.run_stage(integration, project)
    assert result.gate_passed
    links = pq.read_table(result.path / "candidate_links.parquet").to_pydict()
    assert links["physical_section"] == [pair.protein_entry] * 2
    assert links["gene:CD3E"] == [1, -1]
    assert links["protein_candidate_state"] == ["unique_candidate", "unavailable"]
    counts = json.loads((result.path / "pairing_counts.json").read_text())
    assert counts["candidate_links"] == 2
    assert counts["final_cell_classes"] == 0
    assert (result.path / "gate/candidate_links.png").is_file()
    assert len(pq.read_table(result.path / "spatial_candidates.parquet")) == 2
    assert len(pq.read_table(result.path / "cell_link_dispositions.parquet")) == 6
    metadata = json.loads((result.path / "linking_metadata.json").read_text())
    assert metadata["registration_uncertainty"]["status"] == "unavailable"
    parameters = json.loads((result.path / "resolved_parameters.json").read_text())
    assert parameters["automatic_application"] is False
    assert {row["name"] for row in parameters["input_bindings"]} == {
        "registration_receipt", "registration_acceptance", "typing_receipt", "typing_acceptance",
        "transform_json"}


def test_integration_refuses_transform_changed_after_registration_acceptance(tmp_path,
                                                                            monkeypatch):
    project, _ = _pilot(tmp_path, monkeypatch)
    registered = tasks.run_stage(registration, project)
    gates.accept(registered.path, by="test")
    transform = registered.path / "registration/transform.json"
    transform.write_text(transform.read_text() + "\n")
    with pytest.raises(integration.SameSlideIntegrationError, match="changed"):
        tasks.run_stage(integration, project)


def test_integration_refuses_source_manifest_changed_after_registration(tmp_path, monkeypatch):
    project, _ = _pilot(tmp_path, monkeypatch)
    registered = tasks.run_stage(registration, project)
    gates.accept(registered.path, by="test")
    manifest = project.source("same_slide_sources")
    manifest.write_text(manifest.read_text() + "# changed after review\n")
    with pytest.raises(integration.SameSlideIntegrationError,
                       match="source manifest changed"):
        tasks.run_stage(integration, project)


@pytest.mark.parametrize("name", ["protein_cells", "rna_cells", "protein_calls", "rna_calls",
                                  "rna_regions", "typing_receipt", "typing_acceptance"])
def test_current_cells_and_independent_calls_must_match_reviewed_receipts(tmp_path, monkeypatch, name):
    project, _ = _pilot(tmp_path, monkeypatch)
    registered = tasks.run_stage(registration, project)
    gates.accept(registered.path, by="test")
    item = next(item for item in integration.inputs(project) if item.name == name)
    if item.path.suffix == ".parquet":
        table = pq.read_table(item.path).to_pandas()
        table["post_review_change"] = True
        pq.write_table(pa.Table.from_pandas(table, preserve_index=False), item.path)
    else:
        item.path.write_text(item.path.read_text() + "\n")
    with pytest.raises(integration.SameSlideIntegrationError, match="changed|unbound"):
        tasks.run_stage(integration, project)


def test_registration_refuses_cells_changed_since_independent_typing(tmp_path, monkeypatch):
    project, pair = _pilot(tmp_path, monkeypatch)
    path = project.release_dir / f"03_protein_measurement/cells/{pair.protein_entry}.parquet"
    table = pq.read_table(path).to_pandas()
    table.loc[0, "x_um"] += 1
    pq.write_table(pa.Table.from_pandas(table, preserve_index=False), path)
    with pytest.raises(registration.SameSlideRegistrationError, match="protein_cells changed"):
        tasks.run_stage(registration, project)


def test_stage_preserves_and_rejects_mixed_source_frame_metadata(tmp_path, monkeypatch):
    project, pair = _pilot(tmp_path, monkeypatch)
    path = project.release_dir / f"03_protein_measurement/cells/{pair.protein_entry}.parquet"
    table = pq.read_table(path).to_pandas()
    table["coordinate_frame"] = ["crop_local_um", "full_slide_um", "crop_local_um"]
    pq.write_table(pa.Table.from_pandas(table, preserve_index=False), path)
    folder = project.release_dir / "06_typing"
    receipt = json.loads((folder / "receipt.json").read_text())
    next(row for row in receipt["inputs"] if row["name"].startswith("protein_cells:"))["sha256"] = sha256_file(path)
    (folder / "receipt.json").write_text(json.dumps(receipt))
    gates.accept(folder, by="test")
    with pytest.raises(ValueError, match="coordinate frame"):
        tasks.run_stage(registration, project)


def test_duplicate_receipt_bindings_do_not_choose_an_arbitrary_row(tmp_path, monkeypatch):
    project, _ = _pilot(tmp_path, monkeypatch)
    registered = tasks.run_stage(registration, project)
    receipt_path = registered.path / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["inputs"].append(receipt["inputs"][0])
    receipt_path.write_text(json.dumps(receipt))
    gates.accept(registered.path, by="test")
    with pytest.raises(integration.SameSlideIntegrationError, match="duplicate"):
        tasks.run_stage(integration, project)


def test_supplied_uncertainty_is_bound_by_reviewed_transform(tmp_path, monkeypatch):
    project, pair = _pilot(tmp_path, monkeypatch)
    folder = project.source("same_slide_registrations") / pair.protein_entry
    transform = Transform.load(folder)
    transform.meta["registration_uncertainty"] = {
        "bound_um": 0.2, "frame": "crop_local_um", "kind": "radial_displacement_bound",
        "evidence_ref": "synthetic-reviewed-bound"}
    transform.save(folder)
    registered = tasks.run_stage(registration, project)
    gates.accept(registered.path, by="test")
    result = tasks.run_stage(integration, project)
    links = pq.read_table(result.path / "candidate_links.parquet").to_pandas()
    assert links.stable_under_provided_bound.all()
    assert links.registration_uncertainty_um.tolist() == [0.2, 0.2]
    counts = json.loads((result.path / "pairing_counts.json").read_text())
    assert counts["registration_uncertainty"]["evidence_ref"] == "synthetic-reviewed-bound"


def test_stage_retains_actual_rna_calls_and_explicit_level_disagreement(tmp_path, monkeypatch):
    project, pair = _pilot(tmp_path, monkeypatch)
    folder = project.release_dir / "06_typing"
    path = folder / "rna/0110618__090.parquet"
    calls = pq.read_table(path).to_pandas().drop(columns="gene:CD3E")
    calls["cell_class"] = ["B cell", None, None]
    calls["call_state"] = ["resolved", "unavailable", "unavailable"]
    pq.write_table(pa.Table.from_pandas(calls, preserve_index=False), path)
    receipt = json.loads((folder / "receipt.json").read_text())
    next(row for row in receipt["outputs"] if row["path"] == "rna/0110618__090.parquet")["sha256"] = sha256_file(path)
    (folder / "receipt.json").write_text(json.dumps(receipt))
    gates.accept(folder, by="test")
    project.settings["same_slide_hierarchy_map"] = {
        "tissue": "spleen", "level": "lineage", "protein_column": "candidate_phenotype",
        "rna_column": "cell_class", "protein": {"T cell": "T"}, "rna": {"B cell": "B"}}
    registered = tasks.run_stage(registration, project)
    gates.accept(registered.path, by="test")
    result = tasks.run_stage(integration, project)
    links = pq.read_table(result.path / "candidate_links.parquet").to_pandas()
    assert links.rna_cell_class.iloc[0] == "B cell"
    assert links.protein_candidate_phenotype.iloc[0] == "T cell"
    assert links.identity_concordance.tolist() == ["disagree", "unavailable"]
    assert links.comparison_level.eq("lineage").all()
