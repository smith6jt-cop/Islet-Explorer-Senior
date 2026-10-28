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

from multiplex5000 import gates, roots, tasks
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


def test_opted_in_stage_refuses_to_reuse_output_with_changed_settings(tmp_path):
    project = _project(tmp_path)
    stage = _stage("s1")
    stage.REQUIRE_SETTINGS_MATCH_ON_REUSE = True
    current = {"threads": 4}
    stage.settings = lambda _project: dict(current)
    tasks.run_stage(stage, project)
    assert tasks.run_stage(stage, project).status == "complete"
    current["threads"] = 8
    with pytest.raises(tasks.TaskError, match="use --force"):
        tasks.run_stage(stage, project)
    result = tasks.run_stage(stage, project, force=True)
    assert result.status == "published"
    receipt = json.loads((result.path / "receipt.json").read_text())
    assert receipt["settings"] == {"threads": 8}


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


def test_only_bytecode_of_an_exact_declared_python_source_is_audited_as_covered(tmp_path):
    package = tmp_path / "external" / "UniFORM"
    cache = package / "__pycache__"
    cache.mkdir(parents=True)
    declared = package / "registration.py"
    declared.write_text("value = 1\n")
    allowed = cache / "registration.cpython-312.pyc"
    undeclared = cache / "normalization.cpython-312.pyc"
    disguised = package / "registration.cpython-312.pyc"
    for path in (allowed, undeclared, disguised):
        path.write_bytes(b"bytecode fixture")
    opened = {("read", str(path)) for path in (allowed, undeclared, disguised)}
    stray = tasks._undeclared(opened, [tasks.Input("uniform_code", declared)],
                              tmp_path / "staging")
    assert str(allowed) not in stray
    assert set(stray) == {str(undeclared), str(disguised)}


def test_a_stage_s_temporary_files_go_under_the_site_s_tmp_dir_and_are_removed(tmp_path,
                                                                               monkeypatch):
    import tempfile

    root = tmp_path / "tmp"
    monkeypatch.setattr(roots, "tmp_dir", lambda: root)
    before = os.environ.get("TMPDIR")
    before_mpl = os.environ.get("MPLCONFIGDIR")
    seen = {}

    def during():
        seen["folder"] = Path(tempfile.gettempdir())
        seen["TMPDIR"] = os.environ.get("TMPDIR")
        seen["MPLCONFIGDIR"] = os.environ.get("MPLCONFIGDIR")
        with tempfile.NamedTemporaryFile(delete=False) as handle:
            seen["file"] = Path(handle.name)
        font_cache = Path(seen["MPLCONFIGDIR"]) / "fontlist-v3.11.0.json"
        font_cache.write_text("runtime cache")
        assert font_cache.read_text() == "runtime cache"
        seen["font_cache"] = font_cache

    result = tasks.run_stage(_stage("s1", during=during), _project(tmp_path))
    assert result.status == "published"
    # a folder of the run's own, under the site's tmp_dir, which pool workers inherit as TMPDIR
    assert seen["folder"].parent == root and seen["TMPDIR"] == str(seen["folder"])
    assert seen["file"].parent == seen["folder"]
    assert Path(seen["MPLCONFIGDIR"]).parent == seen["folder"]
    assert not seen["font_cache"].exists()
    assert not seen["folder"].exists(), "the run's temporary folder is removed when it ends"
    assert os.environ.get("TMPDIR") == before and tempfile.gettempdir() != str(seen["folder"])
    assert os.environ.get("MPLCONFIGDIR") == before_mpl
    receipt = json.loads((result.path / "receipt.json").read_text())
    assert receipt["environment"]["tmp_dir"] == str(root)


def test_a_file_read_in_a_pool_worker_is_checked_like_the_stage_s_own(tmp_path):
    """A worker is a fresh interpreter; ``tasks.audited`` hands what it opened back to the run."""
    from functools import partial

    from multiplex5000.parallel import map_donors

    project = _project(tmp_path)
    secrets = [_source(tmp_path, "undeclared.csv"), _source(tmp_path, "undeclared_too.csv")]

    def fan_out():
        for _, opened in map_donors(partial(tasks.audited, Path.read_text), secrets, n_jobs=2,
                                    ordered=True, on_error="raise"):
            tasks.record(opened)

    with pytest.raises(tasks.TaskError, match="did not declare"):
        tasks.run_stage(_stage("s1", during=fan_out), project)
    assert not (project.release_dir / "s1").exists()


def test_package_metadata_is_not_an_input_but_a_file_beside_it_is(tmp_path):
    """An editable install leaves ``<name>.egg-info`` beside its source, and scanning entry points
    reads it; that is the package describing itself, not data the stage depends on."""
    metadata = tmp_path / "multiplex5000.egg-info"
    metadata.mkdir()
    (metadata / "entry_points.txt").write_text("[console_scripts]\n")
    beside = _source(tmp_path, "beside.csv")
    opened = {("read", os.path.realpath(metadata / "entry_points.txt")),
              ("read", os.path.realpath(beside))}
    stray = tasks._undeclared(opened, [], tmp_path / "staging")
    assert stray == [os.path.realpath(beside)]


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
