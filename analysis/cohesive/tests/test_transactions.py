from pathlib import Path

import pytest

from analysis.cohesive.common import binding, read_json, run_task, verify_binding
from analysis.cohesive.driver import required_stages


def test_hash_is_required_even_when_size_and_mtime_are_unchanged(tmp_path):
    import os
    path = tmp_path / "source"
    path.write_text("first")
    record = binding(path)
    stamp = path.stat().st_mtime_ns
    path.write_text("other")
    os.utime(path, ns=(stamp, stamp))
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        verify_binding(record)


def test_complete_reuse_verifies_values_and_does_not_run_builder(tmp_path):
    source = tmp_path / "source"
    source.write_text("input")
    calls = []

    def build(directory):
        calls.append(directory)
        (directory / "values.csv").write_text("n\n3\n")
        return {"observations": 3}

    first, receipt = run_task(tmp_path / "run", "example", [binding(source)], {}, build)
    second, _ = run_task(tmp_path / "run", "example", [binding(source)], {}, build)
    assert first == second and len(calls) == 1
    assert receipt["execution_status"] == "complete"
    (first / "values.csv").write_text("n\n4\n")
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        run_task(tmp_path / "run", "example", [binding(source)], {}, build)


def test_failure_cannot_be_reused_and_successful_restart_preserves_evidence(tmp_path):
    source = tmp_path / "source"
    source.write_text("input")

    def fail(directory):
        (directory / "partial.csv").write_text("partial\n")
        raise RuntimeError("deliberate interruption")

    with pytest.raises(RuntimeError, match="deliberate interruption"):
        run_task(tmp_path / "run", "example", [binding(source)], {}, fail)
    failures = list((tmp_path / "run").rglob("failure.json"))
    assert len(failures) == 1
    assert read_json(failures[0])["partial_outputs_are_not_reusable"]

    def succeed(directory):
        (directory / "values.csv").write_text("n\n3\n")
        return {"observations": 3}

    directory, _ = run_task(tmp_path / "run", "example", [binding(source)], {}, succeed)
    assert (directory / "values.csv").is_file()
    assert failures[0].is_file()


def test_new_input_creates_a_distinct_task_without_overwriting(tmp_path):
    source = tmp_path / "source"
    source.write_text("old")

    def build(directory):
        (directory / "result").write_text(source.read_text())
        return {"value": source.read_text()}

    first, _ = run_task(tmp_path / "run", "example", [binding(source)], {}, build)
    source.write_text("new")
    second, _ = run_task(tmp_path / "run", "example", [binding(source)], {}, build)
    assert first != second
    assert (first / "result").read_text() == "old"
    assert (second / "result").read_text() == "new"


def test_missing_or_extra_outputs_refuse_reuse(tmp_path):
    def build(directory):
        (directory / "result").write_text("1")
        return {"value": 1}

    directory, _ = run_task(tmp_path, "example", [], {}, build)
    (directory / "undeclared").write_text("2")
    with pytest.raises(ValueError, match="Output set changed"):
        run_task(tmp_path, "example", [], {}, build)
    (directory / "undeclared").unlink()
    (directory / "result").unlink()
    with pytest.raises(ValueError, match="Output set changed"):
        run_task(tmp_path, "example", [], {}, build)


def test_h6_requires_h3_and_main_delivery_does_not_wait_for_external_owner():
    stages = required_stages("h6")
    assert stages.index("h3") < stages.index("h6")
    assert "h1" not in stages
    assert "same-slide" not in required_stages("verify")
    assert required_stages("atlas") == ["reconcile", "qualify", "atlas"]
