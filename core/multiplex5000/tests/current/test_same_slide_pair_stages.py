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
