"""The joint Stage 06 publishes both assay routes with one gate and receipt."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from multiplex5000 import gates, tasks
from multiplex5000.stages import io60_joint_typing as joint
from multiplex5000.stages import typing


def _project(tmp_path: Path):
    release = tmp_path / "release"
    return SimpleNamespace(
        release_dir=release, release="joint-pilot", path=tmp_path / "joint.yaml",
        cohort="io60_v3", settings={"typing_route": "io60_joint_typing"},
        source=lambda name: None,
        previous=lambda stage: None, automatic_gates=())


def test_joint_stage_keeps_assays_separate_and_reports_both_gates(tmp_path, monkeypatch):
    project = _project(tmp_path)
    shared = tmp_path / "inventory.csv"
    shared.write_text("accepted fixture\n")
    monkeypatch.setattr(joint.protein, "settings", lambda project:
                        {"donors": ["6539"], "scope": "protein pilot"})
    monkeypatch.setattr(joint.rna, "settings", lambda project:
                        {"rna_acquisitions": ["xenium:sample:Panc"], "scope": "RNA pilot"})
    monkeypatch.setattr(joint.protein, "inputs", lambda project:
                        [tasks.Input("inventory", shared)])
    monkeypatch.setattr(joint.rna, "inputs", lambda project:
                        [tasks.Input("inventory", shared)])

    def protein_run(task):
        task.output("pancreas/calls.csv").write_text("protein cell\n")
        gates.figure_path(task, "protein.png").write_bytes(b"protein figure")
        return {"cells": 1}

    def rna_run(task):
        task.output("rna/calls.csv").write_text("RNA cell\n")
        gates.figure_path(task, "rna.png").write_bytes(b"RNA figure")
        return {"cells": 1}

    monkeypatch.setattr(joint.protein, "run", protein_run)
    monkeypatch.setattr(joint.rna, "run", rna_run)
    monkeypatch.setattr(joint.protein, "gate", lambda task, summary:
                        gates.Gate("06_typing", "protein", [gates.Check("cells", 1, "1", True)],
                                   ["protein.png"]))
    monkeypatch.setattr(joint.rna, "gate", lambda task, summary:
                        gates.Gate("06_typing", "RNA", [gates.Check("cells", 1, "1", True)],
                                   ["rna.png"]))

    result = tasks.run_stage(typing, project)
    assert result.gate_passed
    assert (result.path / "pancreas/calls.csv").read_text() == "protein cell\n"
    assert (result.path / "rna/calls.csv").read_text() == "RNA cell\n"
    receipt = json.loads((result.path / "receipt.json").read_text())
    assert [item["name"] for item in receipt["inputs"]] == ["inventory"]
    assert receipt["summary"]["paired_cell_links"] == 0
    gate = gates.load(result.path)
    assert [check["name"] for check in gate["checks"]] == [
        "protein: cells", "RNA: cells", "same-cell links"]
    assert gate["figures"] == ["protein.png", "rna.png"]


def test_joint_route_refuses_sources_from_another_release(tmp_path):
    project = _project(tmp_path)
    project.source = lambda name: tmp_path / "another-release"
    with pytest.raises(joint.IO60JointTypingError, match="this project's release"):
        joint.settings(project)


def test_joint_project_uses_one_release_for_both_assays():
    from multiplex5000 import project
    from multiplex5000.stages import io60_mixed_protein_pilot

    configured = project.load("io60_joint_typing_pilot")
    assert configured.settings["typing_donors"] == ["6539"]
    assert configured.settings["rna_acquisitions"] == [
        "xenium:0059865:Panc", "xenium:0059865:PanT"]
    assert io60_mixed_protein_pilot._upstream(configured) == configured.release_dir
