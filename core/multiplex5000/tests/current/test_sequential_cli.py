import json
import shlex
from pathlib import Path

import pytest

from multiplex5000.integration import sequential_cli as cli
from multiplex5000.integration import sequential_review as review


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


@pytest.fixture
def recipe(tmp_path):
    roots = {key: tmp_path / (key + " source with spaces") for key in ("base_bundle", "current_islets")}
    data = {"schema": "multiplex5000-sequential-recipe/1", "name": "explicit test project"}
    for key, root in roots.items():
        root.mkdir()
        artifact = root / "bound evidence.json"; write(artifact, {"generation": key})
        pin = {"path": artifact.name, "bytes": artifact.stat().st_size, "sha256": cli.digest(artifact)}
        write(root / "artifact_manifest.json", [pin])
        data[key] = {"path": str(root), "manifest_sha256": cli.digest(root / "artifact_manifest.json"), "pins": [pin]}
    native = tmp_path / "native construction"; native.mkdir()
    write(roots["base_bundle"] / "run_specs.json", [{"role": "RNA", "donor_id": "D1", "output": str(native)},
                                                   {"role": "protein", "donor_id": "D1", "output": str(native)}])
    path = tmp_path / "chosen recipe with spaces.json"; write(path, data)
    return path, data, native


def test_dry_run_preserves_explicit_recipe_and_output_from_neutral_directory(recipe, tmp_path, monkeypatch, capsys):
    path, data, native = recipe
    neutral = tmp_path / "neutral"; neutral.mkdir(); monkeypatch.chdir(neutral)
    output = tmp_path / "selected output location with spaces"
    monkeypatch.setattr(cli, "inspect_bundle", lambda p: {"missing_artifacts": [], "bound_bundle": str(p)})
    assert cli.main(["run", "--recipe", str(path), "--output-dir", str(output), "--dry-run"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["recipe"] == str(path)
    assert result["output_dir"] == str(output)
    assert result["corrected_multicell"]["bound_bundle"] == data["base_bundle"]["path"]
    assert result["side_effects"] is False and result["scientific_acceptance"] is False
    assert not output.exists()


def test_recipe_hash_drift_fails_before_output_creation(recipe, tmp_path, capsys):
    path, data, native = recipe
    (Path(data["current_islets"]["path"]) / "bound evidence.json").write_text("changed")
    output = tmp_path / "never created"
    assert cli.main(["run", "--recipe", str(path), "--output-dir", str(output), "--dry-run"]) == 2
    assert "Pinned source changed" in capsys.readouterr().err
    assert not output.exists()


def test_fresh_output_rejects_existing_and_native_source_directories(recipe, tmp_path):
    path, data, native = recipe
    existing = tmp_path / "existing review"; existing.mkdir(); (existing / "note").write_text("keep")
    with pytest.raises(FileExistsError):
        cli._fresh_output(existing, data)
    assert (existing / "note").read_text() == "keep"
    with pytest.raises(ValueError, match="separate"):
        cli._fresh_output(native / "new result", data)
    assert not (native / "new result").exists()


def test_followup_command_quotes_paths_and_preserves_interpreter_site(monkeypatch, tmp_path):
    site = tmp_path / "site chosen with spaces.json"
    recipe = tmp_path / "selected project's recipe.json"
    monkeypatch.setattr(cli, "runtime_context", lambda: {"site_file": str(site)})
    command = shlex.split(cli._command("verify", "--recipe", recipe))
    assert "MPLX5K_SITE_FILE=" + str(site) in command
    assert command[-2:] == ["--recipe", str(recipe)]
    assert cli.sys.executable in command
    assert "PYTHONPATH=" + str(Path(cli.__file__).resolve().parents[2]) in command
    assert command[command.index("-m") + 1:command.index("-m") + 3] == ["multiplex5000.cli", "sequential"]


def test_recipe_cannot_change_between_binding_and_verification(recipe, tmp_path, monkeypatch):
    path, data, native = recipe
    def changed_recipe(p):
        new = json.loads(path.read_text()); new["name"] = "changed while verification ran"; write(path, new)
        return {"status": "PASS", "recipe_sha256": cli.digest(path)}
    monkeypatch.setattr(cli, "verify_recipe", changed_recipe)
    monkeypatch.setattr(review, "build_workspace", lambda *a, **k: {"test": "must not build stale recipe"})
    output = tmp_path / "must not publish stale binding"
    with pytest.raises(ValueError, match="[Rr]ecipe changed"):
        cli.run_recipe(path, output)
    assert not output.exists()


@pytest.fixture
def progress(tmp_path):
    data = {"schema": review.SCHEMA, "dataset_fingerprint": "fixed-source-fingerprint",
            "cases": [{"id": "islet:A", "generation": "current_islets"},
                      {"id": "islet:B", "generation": "current_islets"},
                      {"id": "multicell:C", "generation": "corrected_multicell"}]}
    exported = {"schema": review.SCHEMA, "dataset_fingerprint": data["dataset_fingerprint"],
                "decisions": {"islet:A": {"status": "needs_reassessment", "reviewer": "JT",
                              "note": "Compare the competing neighboring structure", "updated_utc": "2026-11-07T12:00:00Z"}},
                "history": [], "import_history": [], "scientific_designations_changed": False, "production_activated": False}
    source = tmp_path / "source figure"; source.write_bytes(b"unchanged figure")
    workspace = tmp_path / "review with spaces"; workspace.mkdir()
    write(workspace / "review_data.json", data); (workspace / "REVIEW.html").write_text("review")
    write(workspace / "review_manifest.json", {"schema": review.SCHEMA, "dataset_fingerprint": data["dataset_fingerprint"],
          "files": {name: review.sha(workspace / name) for name in ("review_data.json", "REVIEW.html")},
          "source_artifacts": {str(source): {"sha256": review.sha(source), "bytes": source.stat().st_size}}})
    export = tmp_path / "saved progress.json"; write(export, exported)
    return data, exported, workspace, export, source


def test_progress_reports_notes_without_mutating_scientific_evidence(progress):
    data, exported, workspace, export, source = progress
    before = {p: review.sha(p) for p in (source, export, workspace / "review_data.json", workspace / "REVIEW.html")}
    result = review.review_status(workspace, export)
    assert result["counts"] == {"needs_reassessment": 1, "unreviewed": 2}
    assert result["by_generation"]["corrected_multicell"] == {"unreviewed": 1}
    assert result["scientific_designations_changed"] is False and result["production_activated"] is False
    assert result["linked_evidence_hashes_verified"] == 1
    assert before == {p: review.sha(p) for p in before}


@pytest.mark.parametrize("change", ["fingerprint", "unknown_case", "accepted_status"])
def test_foreign_or_acceptance_decisions_are_rejected(progress, change):
    data, exported, *_ = progress
    if change == "fingerprint": exported["dataset_fingerprint"] = "another result"
    elif change == "unknown_case": exported["decisions"]["islet:unknown"] = exported["decisions"].pop("islet:A")
    else: exported["decisions"]["islet:A"]["status"] = "accepted"
    with pytest.raises(ValueError):
        review.validate_progress(data, exported)


@pytest.mark.parametrize("malformed", [None, [], "not an object", 1])
def test_malformed_export_root_has_actionable_validation_error(progress, malformed):
    with pytest.raises(ValueError):
        review.validate_progress(progress[0], malformed)


@pytest.mark.parametrize("change", ["history", "import_history", "timestamp", "production", "designation"])
def test_malformed_history_time_and_approval_claims_are_rejected(progress, change):
    data, exported, *_ = progress
    if change in ("history", "import_history"): exported[change] = "not a history array"
    elif change == "timestamp": exported["decisions"]["islet:A"]["updated_utc"] = "not a timestamp"
    elif change == "production": exported["production_activated"] = True
    elif change == "designation": exported["scientific_designations_changed"] = True
    with pytest.raises(ValueError):
        review.validate_progress(data, exported)


def test_review_status_rejects_changed_source_evidence(progress):
    data, exported, workspace, export, source = progress
    source.write_bytes(b"changed figure")
    with pytest.raises(ValueError, match="Source evidence changed"):
        review.review_status(workspace, export)


def test_review_status_rejects_changed_workspace(progress):
    data, exported, workspace, export, source = progress
    (workspace / "REVIEW.html").write_text("changed review")
    with pytest.raises(ValueError, match="Changed review workspace"):
        review.review_status(workspace, export)
