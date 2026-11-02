"""Completed science and retained publication failures keep distinct bound evidence."""
import importlib.util
import json
import sys
from pathlib import Path

import pytest
from multiplex5000.parameter_observations import file_hash, validate


@pytest.fixture
def recorder(monkeypatch):
    directory = Path(__file__).parents[1]
    monkeypatch.syspath_prepend(str(directory))
    spec = importlib.util.spec_from_file_location("record_remediation_observations", directory / "record_remediation_observations.py")
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)
    return module


def fixture(tmp_path):
    source = tmp_path / "producer.py"
    source.write_text("# synthetic producer\n")
    path = tmp_path / "resolved_parameters.json"
    payload = {"schema": "multiplex5000-resolved-parameters/1", "stage": "06_typing",
               "release": "candidate", "stage_receipt": "receipt.json",
               "context": {"assay": "PhenoCycler", "tissue": "pancreatic_lymph_node"},
               "input_bindings": [], "producer_sources": [{"path": str(source), "sha256": file_hash(source)}],
               "parameters": [{"name": "protein_cutoff_fraction", "role": "fixed", "value": 0.1}],
               "acquisitions": [{"acquisition_id": "protein:S1", "markers": [{"marker": "CD38", "status": "unavailable", "fitted": None}]}]}
    path.write_text(json.dumps(payload))
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"status": "complete", "stage": "06_typing", "release": "candidate",
                                  "inputs": [], "outputs": [
        {"path": path.name, "sha256": file_hash(path)}]}))
    inventory = {"protein:S1": {"acquisition_id": "protein:S1", "panel": "test", "segmentation_id": "seg1"}}
    return path, source, receipt, inventory


def test_valid_observation_preserves_missing_evidence_and_cannot_approve(recorder, tmp_path):
    path, _, _, inventory = fixture(tmp_path)
    value, = recorder.typing_records(path, [], inventory)
    validate(value)
    assert value["status"] == "evaluated"
    assert value["settings"][-1]["value"]["markers"][0]["fitted"] is None
    assert value["context"]["segmentation_ids"] == ["seg1"]
    with pytest.raises(ValueError, match="cannot record"):
        validate({**value, "status": "approved"})


def test_changed_producer_is_not_relabelled_as_original(recorder, tmp_path):
    path, source, _, inventory = fixture(tmp_path)
    source.write_text("# changed\n")
    with pytest.raises(ValueError, match="source changed"):
        list(recorder.typing_records(path, [], inventory))


def test_revision_evidence_remains_available_after_source_changes(recorder, tmp_path):
    evidence = tmp_path / "diagnostics.json"
    evidence.write_text('{"result":1}\n')
    project = tmp_path / "project"
    reference = {"id": "diagnostics", "path": str(evidence), "sha256": file_hash(evidence)}
    value = {"evidence": [dict(reference)], "context": {"input_receipts": []}}
    saved = recorder.freeze_references(value, project)["evidence"][0]
    evidence.write_text('{"result":2}\n')
    assert Path(saved["path"]).read_text() == '{"result":1}\n'
    assert saved["id"] == reference["id"] and file_hash(Path(saved["path"])) == reference["sha256"]
    new = {"evidence": [{**reference, "sha256": file_hash(evidence)}], "context": {"input_receipts": []}}
    assert recorder.freeze_references(new, project)["evidence"][0]["path"] != saved["path"]


def test_mutated_evidence_is_rejected_before_recording(recorder, tmp_path):
    path = tmp_path / "diagnostics.json"
    path.write_text("original")
    value = {"evidence": [{"id": "diagnostics", "path": str(path), "sha256": file_hash(path)}],
             "context": {"input_receipts": []}}
    path.write_text("changed")
    with pytest.raises(ValueError, match="changed before recording"):
        recorder.freeze_references(value, tmp_path / "project")


@pytest.mark.parametrize("kind", ["incomplete", "duplicate", "changed_artifact", "wrong_stage", "wrong_release", "same_basename"])
def test_missing_or_ambiguous_receipt_binding_rejected(recorder, tmp_path, kind):
    path, _, receipt, inventory = fixture(tmp_path)
    data = json.loads(receipt.read_text())
    if kind == "incomplete":
        data["status"] = "failed"
    elif kind == "duplicate":
        data["outputs"] *= 2
    elif kind == "wrong_stage":
        data["stage"] = "03_protein_measurement"
    elif kind == "wrong_release":
        data["release"] = "another-candidate"
    elif kind == "same_basename":
        data["outputs"][0]["path"] = "unrelated/" + path.name
    else:
        path.write_text(path.read_text() + "\n")
    receipt.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="receipt"):
        list(recorder.typing_records(path, [], inventory))


def failed_publication_fixture(recorder, tmp_path, monkeypatch):
    from manifest_remediation_artifacts import FAILED_PROTEIN_ATTEMPT1_REL, RUN_REL

    monkeypatch.setattr(recorder, "BASE", tmp_path)
    monkeypatch.setattr(recorder, "RUN", tmp_path / RUN_REL)
    stage = tmp_path / FAILED_PROTEIN_ATTEMPT1_REL
    stage.mkdir(parents=True)
    donors = [f"donor-{i}" for i in range(25)]
    source = tmp_path / "core/multiplex5000/scripts/typing/refine.py"
    source.parent.mkdir(parents=True)
    source.write_text("current bytes deliberately differ from the old snapshot\n")
    old_source = {"path": str(source), "sha256": "1" * 64}
    paths = {"failure": stage / "failure.json", "resolved_parameters": stage / "resolved_parameters.json",
             "inner_run": stage / "pancreas/run.json", "failed_launch_snapshot": recorder.RUN / "protein/producer_snapshot.json",
             "project": recorder.RUN / "protein/protein_availability_remediation.yaml"}
    objects = {
        "failure": {"stage": "06_typing", "release": stage.parent.name, "status": "failed", "error_type": "TaskError",
                    "error": f"06_typing read file(s) it did not declare: {[str(source)]}"},
        "resolved_parameters": {"schema": "multiplex5000-resolved-parameters/1", "stage": "06_typing",
            "release": stage.parent.name, "producer_sources": [old_source],
            "context": {"modality": "protein", "assay": "PhenoCycler", "tissue": "pancreas_with_embedded_pLN", "units": "declared units"},
            "acquisitions": [{"donor_id": donor} for donor in donors],
            "parameters": [{"name": "fixed_setting", "role": "fixed", "value": 0.4},
                           {"name": "unaudited_fit", "role": "fitted", "value": 17}]},
        "inner_run": {"status": "complete", "donors": donors, "source_scripts": {str(source): old_source["sha256"]}},
        "failed_launch_snapshot": {"sources": [{"path": str(source.relative_to(tmp_path)), "sha256": old_source["sha256"]}]},
        "project": {"release": stage.parent.name, "settings": {"typing_donors": donors}},
    }
    for key, path in paths.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(objects[key]))
    assessment = {key: {"path": str(path), "sha256": file_hash(path)} for key, path in paths.items()}
    assessment.update({"decision": "no_provenance_safe_complete_recovery_import",
                       "framework_findings": {"persisted_outer_input_receipt": False},
                       "required_inputs_without_recovered_content_binding": 127,
                       "scientific_approval": False, "production_release_activated": False})
    assessment_path = recorder.RUN / "final_guard_review/protein_failed_attempt1_no_recovery_v1.json"
    assessment_path.parent.mkdir(parents=True)
    assessment_path.write_text(json.dumps(assessment))
    inventory = {f"acquisition-{i}": {"acquisition_id": f"acquisition-{i}", "donor_id": donor,
                 "cohort": "io60", "modality": "phenocycler", "tissue": "pancreas"} for i, donor in enumerate(donors)}
    return stage, paths, assessment_path, inventory


def test_actual_publication_failure_has_distinct_scope_and_archived_old_evidence(recorder, tmp_path, monkeypatch):
    _, paths, assessment, inventory = failed_publication_fixture(recorder, tmp_path, monkeypatch)
    value = recorder.failed_protein_publication(inventory)
    validate(value)
    assert value["observation_id"] == "remediation-20261102:protein-stage-publication-attempt1"
    assert value["status"] == "failed" and value["operation"] == "protein_stage_publication"
    assert value["context"]["input_receipts"] == []
    assert value["context"]["producer_version"] == "failed-launch-snapshot-sha256:" + file_hash(paths["failed_launch_snapshot"])
    assert len(value["context"]["acquisition_ids"]) == 25
    assert all(score is None for score in value["quality"].values())
    settings = {row["name"]: row["value"] for row in value["settings"]}
    assert "attempt1_unaudited_fit" not in settings
    assert settings["attempt1_publication_evidence_scope"]["complete_outer_receipt"] is None
    assert settings["attempt1_publication_evidence_scope"]["full_execution_input_closure"] is None
    assert settings["attempt1_publication_evidence_scope"]["scientific_parameter_failure"] == "not_established"
    assert settings["attempt1_publication_evidence_scope"]["retained_outputs"] == "unapproved_unpublished"
    frozen = recorder.freeze_references(value, tmp_path / "history")
    assert len(frozen["evidence"]) == 6
    assert {row["sha256"] for row in frozen["evidence"]} == {file_hash(path) for path in [*paths.values(), assessment]}
    assert all(Path(row["path"]).is_file() for row in frozen["evidence"])


@pytest.mark.parametrize("change", ["drift", "receipt", "error", "release", "missing_donor", "duplicate_donor",
                                    "project_scope", "missing_inventory", "ambiguous_inventory", "old_source",
                                    "assessment_approval", "closure", "missing_artifact"])
def test_publication_failure_rejects_changed_identity_scope_or_provenance(recorder, tmp_path, monkeypatch, change):
    stage, paths, assessment_path, inventory = failed_publication_fixture(recorder, tmp_path, monkeypatch)
    assessment = json.loads(assessment_path.read_text())
    if change == "drift":
        paths["failure"].write_text(paths["failure"].read_text() + "\n")
    elif change == "receipt":
        (stage / "receipt.json").write_text('{"status":"complete"}')
    elif change == "missing_artifact":
        paths["inner_run"].unlink()
    elif change == "missing_inventory":
        inventory.pop("acquisition-0")
    elif change == "ambiguous_inventory":
        inventory["other"] = dict(inventory["acquisition-0"])
    elif change == "assessment_approval":
        assessment["scientific_approval"] = True
    elif change == "closure":
        assessment["framework_findings"]["persisted_outer_input_receipt"] = True
    else:
        key = "failure" if change in {"error", "release"} else "project" if change == "project_scope" else "resolved_parameters"
        value = json.loads(paths[key].read_text())
        if change == "error":
            value["error"] = "different unrelated failure"
        elif change == "release":
            value["release"] = "another-attempt"
        elif change == "missing_donor":
            value["acquisitions"].pop()
        elif change == "duplicate_donor":
            value["acquisitions"][0] = value["acquisitions"][1]
        elif change == "project_scope":
            value["settings"]["typing_donors"][0] = "other-donor"
        else:
            value["producer_sources"][0]["sha256"] = "2" * 64
        paths[key].write_text(json.dumps(value))
        assessment[key]["sha256"] = file_hash(paths[key])
    assessment_path.write_text(json.dumps(assessment))
    with pytest.raises((ValueError, OSError)):
        recorder.failed_protein_publication(inventory)


def mixed_fixture(tmp_path):
    path, source, receipt, _ = fixture(tmp_path)
    payload = json.loads(path.read_text())
    payload["context"]["tissue"] = "pancreas_with_embedded_pLN"
    payload["acquisitions"] = [{"donor_id": "0112", "markers_used": ["INS"]}]
    path.write_text(json.dumps(payload))
    artifacts = {
        "pancreas/subpanel_audit.json": [{"donor_id": "0112", "n_cells": 10, "endocrine_quantile": 1.0,
            "islet_chain": True, "islet_levels": {"1": {"CD99": None}},
            "islet_levels_unguarded": {"1": {"CD99": 0.0}}, "per_stratum": {"1": {"cells": 10}}}],
        "pancreas/immune/coverage.json": {"donors_built": ["0112"], "donors_requested": ["0112"],
            "require_min": 0.1, "majority_term": 0.5,
            "not_comparable_criterion": {"high": 0.5, "low": 0.05, "permutations": 10000},
            "terms_not_cross_donor_comparable": {}, "supporting_only_notes": {"T-regulatory": "supporting marker"},
            "per_donor": {"0112": {"state": "built", "n_cells": 10, "n_immune": 2,
                "rule_markers_absent": ["CD38"], "panel_columns_entirely_null": [],
                "immune_class": {}, "immune_subtype": {}, "macrophage_profile": {},
                "term_positive_fraction": {"T-regulatory": {"FOXP3": 0.0}}, "unavailable_rules": {"B": ["CD38"]}}}},
        "pancreas/immune/unavailable_rules.json": {"0112": {"B": ["CD38"]}},
        "pancreas/immune/term_comparability.csv": "rule,marker,n_donors,cross_donor_comparable\nT,FOXP3,1,True\n",
        "pancreas/immune/per_donor_counts.csv": "donor_id,level,label,count,status,reason,count_without_majority_term\n0112,immune_class,B,,unavailable,missing CD38,\n",
    }
    for relative, value in artifacts.items():
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(value) if target.suffix == ".json" else value)
    value = json.loads(receipt.read_text())
    value["settings"] = {"donors": ["0112"]}
    value["outputs"] = [{"path": p, "sha256": file_hash(tmp_path / p)}
                        for p in [path.name, *artifacts]]
    receipt.write_text(json.dumps(value))
    inventory = {"protein:S1": {"acquisition_id": "protein:S1", "panel": "test", "segmentation_id": "seg1",
                               "modality": "phenocycler", "tissue": "pancreas", "donor_id": "0112"}}
    return path, source, receipt, inventory, artifacts


def test_mixed_record_captures_component_roles_missing_values_and_archived_outputs(recorder, tmp_path):
    path, _, _, inventory, artifacts = mixed_fixture(tmp_path)
    value, = recorder.typing_records(path, [], inventory, mixed_components=True)
    validate(value)
    settings = {row["name"]: row for row in value["settings"]}
    assert settings["primary_endocrine_quantile"]["kind"] == "fixed"
    assert settings["primary_islet_chain_selected"]["kind"] == "selected_method"
    assert settings["primary_islet_levels"]["kind"] == "fitted"
    assert settings["primary_islet_levels"]["value"]["1"]["CD99"] is None
    assert settings["primary_islet_levels_unguarded_diagnostics"]["value"]["1"]["CD99"] == 0.0
    assert settings["immune_require_min"]["value"] == 0.1
    assert settings["immune_majority_term"]["value"] == 0.5
    assert settings["immune_comparability_permutations"]["kind"] == "fixed"
    assert settings["immune_not_comparable_criterion"]["kind"] == "review_criterion"
    assert settings["immune_section_counts"]["value"][0]["count"] is None
    assert settings["immune_term_comparability_results"]["value"][0]["cross_donor_comparable"] is True
    assert value["quality"]["biological_acceptance"] is None
    assert {str(tmp_path / p) for p in artifacts} <= {r["path"] for r in value["evidence"]}
    frozen = recorder.freeze_references(value, tmp_path / "history")
    before = {r["id"]: Path(r["path"]).read_bytes() for r in frozen["evidence"]}
    (tmp_path / "pancreas/immune/coverage.json").write_text("changed after archive")
    assert all(Path(r["path"]).read_bytes() == before[r["id"]] for r in frozen["evidence"])
    with pytest.raises(ValueError, match="cannot record"):
        validate({**value, "status": "approved"})


@pytest.mark.parametrize("kind", ["missing_file", "missing_binding", "duplicate_binding", "stale", "wrong_relative",
                                  "wrong_stage", "primary_duplicate", "primary_missing", "coverage_missing",
                                  "coverage_duplicate_key", "missing_donor", "counts_duplicate", "counts_unknown",
                                  "unavailable_disagreement", "path_escape", "nonfinite_fixed", "wrong_roster",
                                  "primary_fixed", "immune_missing_fit", "invalid_comparability"])
def test_mixed_component_evidence_fails_closed(recorder, tmp_path, kind):
    path, _, receipt, inventory, _ = mixed_fixture(tmp_path)
    value = json.loads(receipt.read_text())
    relative = "pancreas/immune/coverage.json"
    coverage = tmp_path / relative
    primary = tmp_path / "pancreas/subpanel_audit.json"
    counts = tmp_path / "pancreas/immune/per_donor_counts.csv"
    if kind == "missing_file":
        coverage.unlink()
    elif kind == "missing_binding":
        value["outputs"] = [r for r in value["outputs"] if r["path"] != relative]
    elif kind == "duplicate_binding":
        value["outputs"].append(next(r for r in value["outputs"] if r["path"] == relative))
    elif kind == "stale":
        coverage.write_text(coverage.read_text() + "\n")
    elif kind == "wrong_relative":
        next(r for r in value["outputs"] if r["path"] == relative)["path"] = "other/coverage.json"
    elif kind == "wrong_stage":
        value["stage"] = "03_protein_measurement"
    elif kind in {"primary_duplicate", "primary_missing", "primary_fixed"}:
        data = json.loads(primary.read_text())
        if kind == "primary_duplicate":
            data *= 2
        elif kind == "primary_missing":
            del data[0]["per_stratum"]
        else:
            data[0]["endocrine_quantile"] = None
        primary.write_text(json.dumps(data))
    elif kind in {"coverage_missing", "missing_donor", "nonfinite_fixed", "immune_missing_fit"}:
        data = json.loads(coverage.read_text())
        if kind == "coverage_missing":
            del data["require_min"]
        elif kind == "missing_donor":
            data["per_donor"] = {}
        elif kind == "immune_missing_fit":
            del data["per_donor"]["0112"]["term_positive_fraction"]
        else:
            data["majority_term"] = float("nan")
        coverage.write_text(json.dumps(data))
    elif kind == "coverage_duplicate_key":
        coverage.write_text(coverage.read_text()[:-1] + ', "per_donor": {}}')
    elif kind == "counts_duplicate":
        counts.write_text(counts.read_text() + counts.read_text().splitlines()[1] + "\n")
    elif kind == "counts_unknown":
        counts.write_text(counts.read_text().replace("0112,", "unknown,"))
    elif kind == "invalid_comparability":
        (tmp_path / "pancreas/immune/term_comparability.csv").write_text("rule,marker,n_donors,cross_donor_comparable\nT,FOXP3,bad,True\n")
    elif kind == "unavailable_disagreement":
        (tmp_path / "pancreas/immune/unavailable_rules.json").write_text('{"0112": {}}')
    elif kind == "path_escape":
        elsewhere = tmp_path.parent / (tmp_path.name + "-external.json")
        elsewhere.write_text(coverage.read_text())
        coverage.unlink()
        coverage.symlink_to(elsewhere)
    else:
        value["settings"]["donors"] = ["9999"]
    # Rebind malformed payloads to distinguish schema rejection from stale-byte rejection.
    if kind not in {"missing_file", "missing_binding", "duplicate_binding", "stale", "wrong_relative"}:
        for row in value["outputs"]:
            row["sha256"] = file_hash(tmp_path / row["path"])
    receipt.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        list(recorder.typing_records(path, [], inventory, mixed_components=True))
