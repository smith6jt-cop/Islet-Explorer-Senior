"""The bounded IO60 Stage 06 route keeps unavailable channels and paths honest."""

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from multiplex5000 import project, tasks
from multiplex5000.stages import io60_mixed_protein_pilot as mixed
from multiplex5000.stages import typing
from multiplex5000.tasks import Input, Task


def test_parameter_producer_hashes_obey_actual_stage_read_audit(tmp_path, monkeypatch):
    """Producer hashing is a real file read and must be declared before execution."""
    monkeypatch.setattr(mixed, "_cohort_rows", lambda: [])
    monkeypatch.setattr(mixed, "_donors", lambda project: [])
    monkeypatch.setattr(mixed, "_source", lambda project, name: tmp_path / name)
    selected = SimpleNamespace(settings={})
    declared = mixed.inputs(selected)
    staging = tmp_path / "stage"
    staging.mkdir()
    task = Task(mixed.STAGE, "candidate", staging,
                {item.name: item for item in declared}, project=selected)
    opened = set()
    token = tasks._OPENED.set(opened)
    try:
        result = mixed.write_resolved_parameters(
            task, operation="protein_typing", context={}, parameters=[], acquisitions=[],
            producer_paths=(Path(mixed.__file__), mixed.RUNNER, mixed.RUNNER.parent / "refine.py"))
    finally:
        tasks._OPENED.reset(token)
    assert result["status"] == "written"
    assert tasks._undeclared(opened, declared, staging) == []


def test_tracked_mixed_project_selects_the_single_stage_module():
    selected = project.load("io60_mixed_protein_pilot")
    assert selected.stages == ("06_typing",)
    assert selected.automatic_gates == ()
    assert typing._route(selected) is mixed


def test_refused_channel_is_not_read_as_a_zero_scale(tmp_path, monkeypatch):
    path = tmp_path / "scales.csv"
    path.write_text("donor_id,marker,status\n6450,CD8,aligned\n"
                    "6450,CD57,refused\n6534,CD8,aligned\n6534,CD57,refused\n")
    monkeypatch.setattr(mixed, "_cohort_rows", lambda: [{"donor_id": "6450"},
                                                    {"donor_id": "6534"}])
    task = Task(mixed.STAGE, "pilot", tmp_path, {"normalization_scales": Input(
        "normalization_scales", path)})
    assert mixed._marker_roster(task) == ["CD8"]


def test_inner_log_names_published_outputs_but_preserves_execution_command(tmp_path):
    release = tmp_path / "release"
    task = Task(mixed.STAGE, "pilot", release / ".06_typing.incomplete-123", {},
                project=SimpleNamespace(release_dir=release))
    command = ["python", "run_calibration_arm.py", "--out", str(task.staging / "pancreas")]
    recorded = {"stages": [{"log": str(task.staging / "pancreas/subpanel.log"),
                            "command": command}]}
    result = mixed._published_inner_manifest(task, recorded)
    assert result["typing_root"] == str(release / "06_typing/pancreas")
    assert result["scale"] == str(release / "06_typing/pancreas_scale.parquet")
    assert result["stages"][0]["log"] == str(release / "06_typing/pancreas/subpanel.log")
    assert result["stages"][0]["command"] == command
    assert Path(result["stages"][0]["log"]).is_relative_to(release)


def test_refinement_provenance_paths_survive_stage_publication(tmp_path):
    release = tmp_path / "release"
    task = Task(mixed.STAGE, "candidate", release / ".06_typing.incomplete-123", {},
                project=SimpleNamespace(release_dir=release))
    folder = task.staging / "pancreas/refinement"
    folder.mkdir(parents=True)
    old = {"parameter_path": str(folder / "representation_params.parquet"),
           "scale_path": str(task.staging / "pancreas_scale.parquet"), "parameter_sha256": "retained"}
    (folder / "representation_provenance.json").write_text(json.dumps(old))
    (folder / "refine_log.json").write_text(json.dumps([{"measurement_policy": {"scale_path": old["scale_path"]}}]))
    mixed._publish_refinement_paths(task)
    result = json.loads((folder / "representation_provenance.json").read_text())
    assert result["parameter_sha256"] == "retained"
    assert result["execution_output_paths"]["parameter_path"] == old["parameter_path"]
    assert result["parameter_path"] == str(release / "06_typing/pancreas/refinement/representation_params.parquet")
    logs = json.loads((folder / "refine_log.json").read_text())
    assert logs[0]["measurement_policy"]["scale_path"] == str(release / "06_typing/pancreas_scale.parquet")


def test_marker_contract_is_bound_to_the_actual_producer_receipt(tmp_path):
    availability, scales = tmp_path / "availability.parquet", tmp_path / "scales.csv"
    availability.write_bytes(b"acquired marker contract")
    scales.write_bytes(b"alignment status contract")
    receipts = []
    for stage, path, rel in (("03_protein_measurement", availability, "data/availability/donor_id=0112/data.parquet"),
                             ("05_normalization", scales, "scales.csv")):
        receipt = tmp_path / f"{stage}.json"
        receipt.write_text(json.dumps({"stage": stage, "status": "complete", "outputs": [
            {"path": rel, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}]}))
        receipts.append(receipt)
    task = Task(mixed.STAGE, "candidate", tmp_path, {
        name: Input(name, path) for name, path in (("availability:0112", availability),
        ("normalization_scales", scales), ("measurement_receipt", receipts[0]),
        ("normalization_receipt", receipts[1]))}, settings={"donors": ["0112"]})
    mixed._verify_availability_inputs(task)
    availability.write_bytes(b"different marker contract")
    with pytest.raises(mixed.IO60MixedPilotError, match="producer receipt"):
        mixed._verify_availability_inputs(task)
