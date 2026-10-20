"""The task contract and the lifecycle repairs it carries (MultiPlex5000 D3).

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_tasks.py -q

Each repair the 2026-10-08 audit called for is held here on synthetic stages: a partial output is
never complete, --force really reruns (and the old output is moved aside, never deleted), every
input is declared and bound when the stage runs, and the next stage waits for the gate.
"""

from __future__ import annotations

import json
import os
import types
from pathlib import Path

import pytest

from multiplex5000 import gates, tasks
from multiplex5000 import project as P


def _project(tmp_path, stages=("s1", "s2"), automatic=()) -> P.Project:
    return P.Project(tmp_path / "p.yaml", "p", "r1", "io60_phenocycler", tuple(stages),
                     tuple(automatic), tmp_path / "releases")


def _stage(stage_id, *, reads=(), fail=False, passes=True, text="out", extra_read=None,
           during=None):
    module = types.SimpleNamespace(STAGE=stage_id)
    module.inputs = lambda project: [tasks.Input(name, path) for name, path in reads]

    def run(task):
        for name, _ in reads:
            task.input(name).read_text()
        if extra_read is not None:
            Path(extra_read).read_text()
        if during is not None:
            during()
        if fail:
            (task.staging / "half.txt").write_text("partial")
            raise RuntimeError("the stage crashed")
        task.output("out.txt").write_text(text)
        return {"rows": 1}

    module.run = run
    module.gate = lambda task, summary: gates.Gate(
        stage_id, "synthetic", [gates.Check("the check", 1, "1", passes)], [])
    return module


def _source(tmp_path, name="in.csv", text="a,b\n1,2\n") -> Path:
    path = tmp_path / name
    path.write_text(text)
    return path


def test_a_stage_publishes_its_outputs_gate_and_receipt(tmp_path):
    project = _project(tmp_path)
    source = _source(tmp_path)
    result = tasks.run_stage(_stage("s1", reads=[("table", source)]), project)
    assert result.status == "published" and result.gate_passed
    folder = project.release_dir / "s1"
    assert (folder / "out.txt").read_text() == "out"
    receipt = json.loads((folder / "receipt.json").read_text())
    assert receipt["status"] == "complete" and receipt["stage"] == "s1"
    assert receipt["inputs"][0]["sha256"] and receipt["inputs"][0]["name"] == "table"
    assert {o["path"] for o in receipt["outputs"]} == {"out.txt", "gate/gate.json",
                                                       "gate/gate.html"}
    assert not tasks.failed_attempts(project.release_dir)


def test_a_crashed_stage_leaves_no_output_and_is_never_complete(tmp_path):
    project = _project(tmp_path)
    with pytest.raises(RuntimeError, match="crashed"):
        tasks.run_stage(_stage("s1", fail=True), project)
    folder = project.release_dir / "s1"
    assert not folder.exists() and not tasks.is_complete(folder)
    [attempt] = tasks.failed_attempts(project.release_dir)
    failure = json.loads((attempt / "failure.json").read_text())
    assert failure["status"] == "failed" and "crashed" in failure["error"]
    assert (attempt / "half.txt").exists(), "the partial files stay where they are, unpublished"
    result = tasks.run_stage(_stage("s1"), project)
    assert result.status == "published" and tasks.is_complete(folder)


def test_without_force_a_complete_stage_is_left_alone(tmp_path):
    project = _project(tmp_path)
    tasks.run_stage(_stage("s1", text="first"), project)
    again = tasks.run_stage(_stage("s1", text="second"), project)
    assert again.status == "complete"
    assert (project.release_dir / "s1" / "out.txt").read_text() == "first"


def test_force_reruns_and_moves_the_old_output_and_its_acceptance_aside(tmp_path):
    project = _project(tmp_path)
    tasks.run_stage(_stage("s1", text="first"), project)
    gates.accept(project.release_dir / "s1", by="tester")
    result = tasks.run_stage(_stage("s1", text="second"), project, force=True)
    assert (project.release_dir / "s1" / "out.txt").read_text() == "second"
    assert result.superseded is not None and result.superseded.name.startswith(".s1.superseded-")
    assert (result.superseded / "out.txt").read_text() == "first", "moved aside, not deleted"
    assert (result.superseded / "gate" / "accepted.json").exists()
    assert gates.accepted(project.release_dir / "s1") is None, "a new result is accepted again"
    with pytest.raises(gates.GateError, match="not been accepted"):
        tasks.run_stage(_stage("s2"), project)


def test_an_empty_placeholder_folder_is_not_kept_as_an_old_output(tmp_path):
    """Snakemake makes the output's folder before a job runs; an empty one holds no output."""
    project = _project(tmp_path)
    (project.release_dir / "s1").mkdir(parents=True)
    result = tasks.run_stage(_stage("s1"), project)
    assert result.superseded is None
    assert not list(project.release_dir.glob(".s1.superseded-*"))


def test_reading_a_file_the_stage_did_not_declare_fails_it(tmp_path):
    project = _project(tmp_path)
    secret = _source(tmp_path, "undeclared.csv")
    with pytest.raises(tasks.TaskError, match="did not declare"):
        tasks.run_stage(_stage("s1", extra_read=secret), project)
    assert not (project.release_dir / "s1").exists()


def test_asking_for_an_undeclared_input_by_name_fails(tmp_path):
    task = tasks.Task("s1", "r1", tmp_path, {})
    with pytest.raises(tasks.TaskError, match="did not declare"):
        task.input("table")


def test_an_input_that_changes_while_the_stage_runs_fails_it(tmp_path):
    project = _project(tmp_path)
    source = _source(tmp_path)

    def touch():
        source.write_text("a,b\n1,2\n3,4\n")
        os.utime(source, ns=(1, 1))

    with pytest.raises(tasks.TaskError, match="changed while the stage ran"):
        tasks.run_stage(_stage("s1", reads=[("table", source)], during=touch), project)


def test_a_missing_input_is_refused_and_an_optional_one_is_recorded_absent(tmp_path):
    with pytest.raises(tasks.TaskError, match="missing"):
        tasks.bind(tasks.Input("table", tmp_path / "nope.csv"))
    record = tasks.bind(tasks.Input("copy", tmp_path / "nope.csv", required=False))
    assert record["exists"] is False
    raw = _source(tmp_path, "big.qptiff", "x" * 10)
    bound = tasks.bind(tasks.Input("raw", raw, hash=False, recorded_sha256="abc"))
    assert "sha256" not in bound and bound["bytes"] == 10 and bound["recorded_sha256"] == "abc"


def test_the_next_stage_waits_for_the_previous_gate(tmp_path):
    project = _project(tmp_path)
    with pytest.raises(gates.GateError, match="has not published"):
        tasks.run_stage(_stage("s2"), project)
    tasks.run_stage(_stage("s1"), project)
    with pytest.raises(gates.GateError, match="not been accepted"):
        tasks.run_stage(_stage("s2"), project)
    gates.accept(project.release_dir / "s1", by="tester")
    assert tasks.run_stage(_stage("s2"), project).status == "published"


def test_an_automatic_gate_accepts_itself_only_when_it_passes(tmp_path):
    project = _project(tmp_path, automatic=("s1",))
    result = tasks.run_stage(_stage("s1"), project)
    assert result.accepted_automatically
    record = gates.accepted(project.release_dir / "s1")
    assert record["accepted_by"] == "automatic (project file)"
    failing = _project(tmp_path / "other", automatic=("s1",))
    result = tasks.run_stage(_stage("s1", passes=False), failing)
    assert not result.accepted_automatically and not result.gate_passed
    assert gates.accepted(failing.release_dir / "s1") is None
