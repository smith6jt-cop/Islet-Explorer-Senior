"""Launch capture must bind the configured project without overwriting history."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SPEC = importlib.util.spec_from_file_location(
    "capture_protein_launch", Path(__file__).parents[1] / "capture_protein_launch.py")
CAPTURE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CAPTURE)


def fixture(tmp_path):
    for relative in CAPTURE.EXTRA_SOURCES:
        source = tmp_path / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text("frozen producer\n")
    base = tmp_path / "old.json"
    base.write_text('{"sources": []}\n')
    config = tmp_path / "project.yaml"
    config.write_text("id: fixture\n")
    measured = tmp_path / "cells.parquet"
    measured.write_bytes(b"measured values")
    def bind(item):
        return {"name": item.name, "path": str(item.path), "exists": True,
                "sha256": CAPTURE.sha(item.path), "bytes": item.path.stat().st_size}
    return {"repo": tmp_path, "project_path": config, "project_sha256": CAPTURE.sha(config),
            "declared_inputs": [SimpleNamespace(name="aligned:0112", path=measured)],
            "settings": {"donors": ["0112"]}, "bind_input": bind,
            "base_snapshot": base, "snapshot": tmp_path / "new.json", "archive": tmp_path / "archive",
            "bindings": tmp_path / "bindings.json", "attempt": 3,
            "code": {"branch": "fixture", "commit": "fixture", "working_changes": True}}


def test_capture_binds_exact_project_sources_and_inputs_without_starting_stage(tmp_path):
    args = fixture(tmp_path)
    result = CAPTURE.capture(**args)
    launch = json.loads(args["bindings"].read_text())
    assert result["stage_execution"] == "not_started_by_capture"
    assert launch["schema"] == "protein-prelaunch-bindings/1"
    assert launch["project"] == {"path": str(args["project_path"]), "sha256": args["project_sha256"]}
    assert launch["inputs"][0]["sha256"] == CAPTURE.sha(args["declared_inputs"][0].path)
    assert launch["producer_snapshot"]["sha256"] == CAPTURE.sha(args["snapshot"])
    for row in json.loads(args["snapshot"].read_text())["sources"]:
        assert CAPTURE.sha(args["archive"] / row["path"]) == row["sha256"]


def test_duplicate_input_names_fail_before_outputs_exist(tmp_path):
    args = fixture(tmp_path)
    args["declared_inputs"] *= 2
    with pytest.raises(ValueError, match="unique"):
        CAPTURE.capture(**args)
    assert not args["archive"].exists() and not args["snapshot"].exists()


@pytest.mark.parametrize("target", ["snapshot", "bindings", "archive"])
def test_existing_capture_outputs_are_never_overwritten(tmp_path, target):
    args = fixture(tmp_path)
    args[target].write_bytes(b"immutable prior evidence")
    with pytest.raises(FileExistsError, match="fresh"):
        CAPTURE.capture(**args)
    assert args[target].read_bytes() == b"immutable prior evidence"


def test_cli_target_resolution_does_not_hide_an_existing_dangling_symlink(tmp_path):
    target = tmp_path / "existing-link"
    target.symlink_to(tmp_path / "absent-target")
    with pytest.raises(FileExistsError, match="already exists"):
        CAPTURE.fresh_target(target)


def test_changed_project_identity_is_refused(tmp_path):
    args = fixture(tmp_path)
    args["project_path"].write_text("id: different\n")
    with pytest.raises(ValueError, match="project changed"):
        CAPTURE.capture(**args)


@pytest.mark.parametrize("mutation", ["config", "input_identity", "snapshot"])
def test_mutation_during_binding_does_not_publish_capture(tmp_path, mutation):
    args = fixture(tmp_path)
    original = args["bind_input"]
    def bind(item):
        row = original(item)
        if mutation == "config":
            args["project_path"].write_text("id: changed\n")
        elif mutation == "snapshot":
            args["base_snapshot"].write_text('{"sources": ["changed"]}\n')
        else:
            row["name"] = "another-input"
        return row
    args["bind_input"] = bind
    with pytest.raises(ValueError):
        CAPTURE.capture(**args)
    assert not args["bindings"].exists() and not args["snapshot"].exists()
