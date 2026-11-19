"""Delivery checks cannot hide pending cohorts, source drift or weak hash claims."""
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest


@pytest.fixture
def modules(monkeypatch):
    directory = Path(__file__).parents[1]
    monkeypatch.syspath_prepend(str(directory))
    found = []
    for name in ["manifest_remediation_artifacts", "verify_remediation_delivery"]:
        spec = importlib.util.spec_from_file_location(name, directory / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        monkeypatch.setitem(sys.modules, name, module)
        spec.loader.exec_module(module)
        found.append(module)
    return found


def write(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))


def ledger_fixture(tmp_path, verifier):
    original = [{"acquisition_id": f"acquisition-{i}", "disposition": "completed"} for i in range(201)]
    write(tmp_path / "docs/evaluations/dynamic_portability_20261102/acquisition_dispositions.json", original)
    folder = tmp_path / verifier.DOC_REL
    current = [{"acquisition_id": row["acquisition_id"], "baseline": dict(row),
                "current": {"disposition": "completed", "reason": "synthetic verified candidate",
                            "production_release_activated": False, "scientific_approval": "not_granted_by_evaluation"}}
               for row in original]
    write(folder / "acquisitions.json", {"acquisitions": current, "declared_acquisitions": 201,
          "disposition_counts": {"completed": 201}, "baseline_roster_unchanged": True})
    write(folder / "readiness.json", {"implementation_items": [{"id": "synthetic", "status": "completed_with_scope_limit"}],
          "routes": [{"axes": {key: "separate evidence" for key in verifier.AXES}}],
          "scientific_approvals_added": False, "production_release_activated": False})
    (folder / "README.md").write_text("201 completed candidates; no scientific approval.\n")
    return folder


def test_complete_ledger_counts_and_scope_are_checked(modules, tmp_path):
    _, verifier = modules
    ledger_fixture(tmp_path, verifier)
    assert verifier.check_accounting(tmp_path)["counts"] == {"completed": 201}


@pytest.mark.parametrize("change", ["pending", "missing", "duplicate", "baseline", "approval", "counts", "implementation"])
def test_incomplete_or_rewritten_ledger_cannot_pass(modules, tmp_path, change):
    _, verifier = modules
    folder = ledger_fixture(tmp_path, verifier)
    path = folder / "acquisitions.json"
    data = json.loads(path.read_text())
    if change == "pending":
        data["acquisitions"][0]["current"]["disposition"] = "pending"
    elif change == "missing":
        data["acquisitions"].pop()
    elif change == "duplicate":
        data["acquisitions"][0] = data["acquisitions"][1]
    elif change == "baseline":
        data["acquisitions"][0]["baseline"]["disposition"] = "blocked"
    elif change == "approval":
        data["acquisitions"][0]["current"]["scientific_approval"] = True
    elif change == "counts":
        data["disposition_counts"] = {"completed": 200}
    else:
        other = json.loads((folder / "readiness.json").read_text())
        other["implementation_items"][0]["status"] = "execution_pending"
        write(folder / "readiness.json", other)
    write(path, data)
    with pytest.raises(ValueError):
        verifier.check_accounting(tmp_path)


def test_current_source_drift_cannot_fall_back_to_git(modules, tmp_path):
    manifest, _ = modules
    source = tmp_path / "producer.py"
    source.write_text("# delivered current source\n")
    row = {"kind": "current_source", "verification": "fresh_sha256", **manifest.ref(source),
           "sourceGit": "a-historical-commit-is-not-current-code"}
    assert manifest.verify_artifact(row, tmp_path) == "fresh_sha256"
    source.write_text("# later different source\n")
    with pytest.raises(ValueError, match="current file changed"):
        manifest.verify_artifact(row, tmp_path)
    row["kind"] = "historical_git_fallback"
    with pytest.raises(ValueError, match="unknown artifact"):
        manifest.verify_artifact(row, tmp_path)


def audit_fixture(tmp_path, manifest):
    audit = tmp_path / "audit/receipt_checks.csv"
    receipt = audit.parent / "upstream/06_typing/receipt.json"
    release = tmp_path / "large-scientific-release"
    large = release / "06_typing/large.parquet"
    digest = "1" * 64
    write(receipt, {"stage": "06_typing", "status": "complete", "inputs": [],
                    "outputs": [{"path": "large.parquet", "sha256": digest, "bytes": 150_000_000_000}]})
    with audit.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["role", "path", "expected", "actual", "status"])
        writer.writeheader()
        writer.writerow({"role": "output", "path": str(large), "expected": digest, "actual": digest, "status": "pass"})
    return manifest.receipt_rows(audit, release)[0], large, receipt, audit


def test_large_data_is_explicitly_audit_bound_without_second_hash(modules, tmp_path, monkeypatch):
    manifest, _ = modules
    row, large, _, _ = audit_fixture(tmp_path, manifest)
    original = manifest.sha

    def guarded(path):
        assert path != large, "delivery verification must not reread150GB"
        return original(path)

    monkeypatch.setattr(manifest, "sha", guarded)
    assert manifest.verify_artifact(row, tmp_path) == "bound_to_prior_audit_not_rehashed"
    assert not large.exists()  # Remote/moved raw bytes are not represented as freshly checked.
    row["freshly_hashed_by_delivery_verifier"] = True
    with pytest.raises(ValueError, match="must not claim"):
        manifest.verify_artifact(row, tmp_path)


@pytest.mark.parametrize("change", ["receipt", "audit", "reported_digest", "missing_binding", "negative_index"])
def test_scientific_binding_needs_exact_successful_audit_and_receipt(modules, tmp_path, change):
    manifest, _ = modules
    row, _, receipt, audit = audit_fixture(tmp_path, manifest)
    if change == "receipt":
        receipt.write_text(receipt.read_text() + "\n")
    elif change == "audit":
        audit.write_text(audit.read_text().replace("pass", "fail"))
    elif change == "reported_digest":
        row["sha256"] = "2" * 64
    elif change == "negative_index":
        row["audit_row"] = -1
    else:
        data = json.loads(receipt.read_text())
        data["outputs"] = []
        write(receipt, data)
        row["receipt"] = manifest.ref(receipt)
    with pytest.raises(ValueError):
        manifest.verify_artifact(row, tmp_path)


def test_sampled_historical_hash_is_not_promoted_to_full_historical_bytes(modules, tmp_path):
    manifest, _ = modules
    audit = tmp_path / "verification.json"
    item = {"name": "qptiff:synthetic", "path": "/not-reread/large.qptiff", "sha256": "1" * 64,
            "bytes": 150_000_000_000, "recorded_hash_mode": "sampled", "status": "pass"}
    write(audit, {"complete": True, "status": "pass", "inputs": [item]})
    row = {"kind": "audited_recorded_input", "verification": "prior_full_current_bytes_with_declared_historical_algorithm",
           **{key: item[key] for key in ["path", "sha256", "bytes", "recorded_hash_mode"]},
           "input_name": item["name"], "audit": manifest.ref(audit), "freshly_hashed_by_delivery_verifier": False}
    assert manifest.verify_artifact(row, tmp_path) == "bound_to_prior_audit_not_rehashed"
    row["recorded_hash_mode"] = "full"
    with pytest.raises(ValueError, match="algorithm or binding"):
        manifest.verify_artifact(row, tmp_path)


def test_current_binding_checks_size_as_well_as_digest(modules, tmp_path):
    manifest, _ = modules
    path = tmp_path / "source.py"
    path.write_text("pass\n")
    row = {"kind": "current_source", "verification": "fresh_sha256", **manifest.ref(path), "bytes": 1}
    with pytest.raises(ValueError, match="file size differs"):
        manifest.verify_artifact(row, tmp_path)


def test_recovered_file_binding_retains_declared_origin_and_exact_hash(modules, tmp_path):
    manifest, _ = modules
    audit = tmp_path / "v20.json"
    digest = "1" * 64
    write(audit, {"source_release_modified": False, "checks": [{"path": "/recovered/table.parquet",
          "declared_path": "/original/table.parquet", "expected": digest, "actual": digest, "status": "pass"}]})
    row = {"kind": "audited_content_recovery", "verification": "prior_exact_release_content_recovery",
           "path": "/recovered/table.parquet", "declared_path": "/original/table.parquet", "sha256": digest,
           "audit": manifest.ref(audit), "freshly_hashed_by_delivery_verifier": False}
    assert manifest.verify_artifact(row, tmp_path) == "bound_to_prior_audit_not_rehashed"
    row["declared_path"] = "/other/table.parquet"
    with pytest.raises(ValueError, match="recovered content binding differs"):
        manifest.verify_artifact(row, tmp_path)


def test_empty_manifest_cannot_hide_missing_current_sources(modules, tmp_path, monkeypatch):
    manifest, verifier = modules
    path = tmp_path / "source.py"
    path.write_text("pass\n")
    monkeypatch.setattr(verifier, "catalog", lambda *_: ({path}, {path}))
    write(tmp_path / verifier.DOC_REL / "artifact_manifest.json", {"schema": manifest.SCHEMA, "artifacts": [], "counts": {}})
    with pytest.raises(ValueError, match="manifest omits"):
        verifier.check_manifest(tmp_path)


def test_manifest_includes_only_selected_small_failed_stage_evidence(modules, tmp_path, monkeypatch):
    manifest, _ = modules
    selected = {tmp_path / manifest.FAILED_PROTEIN_ATTEMPT1_REL / relative for relative in manifest.FAILED_PROTEIN_JSON}
    for path in selected:
        write(path, {"retained": "failed attempt evidence"})
    large = tmp_path / manifest.FAILED_PROTEIN_ATTEMPT1_REL / "pancreas/cells/data.parquet"
    large.parent.mkdir(parents=True)
    large.write_bytes(b"excluded scientific outputs are not rehashed by this manifest")
    write(tmp_path / manifest.RUN_REL / "codex/v20_final_verification/summary.json",
          {"source_release": "small-source.json", "recovery_manifest": "small-recovery.json"})
    monkeypatch.setattr(manifest.subprocess, "check_output", lambda *args, **kwargs: "")
    destination = tmp_path / manifest.DOC_REL / "artifact_manifest.json"
    files, sources = manifest.catalog(tmp_path, destination)
    assert selected <= files and not selected & sources and large not in files
    assert all(manifest._excluded(path.relative_to(tmp_path / manifest.RUN_REL)) for path in selected)
    next(iter(selected)).unlink()
    with pytest.raises(ValueError, match="failed-attempt evidence is missing"):
        manifest.catalog(tmp_path, destination)


def absence_fixture(tmp_path, manifest):
    from verify_protein_remediation import (
        verify_launch_bindings,
        verify_optional_declaration,
    )

    folder = tmp_path / "verification"
    audit = folder / "release_audit/receipt_checks.csv"
    receipt = audit.parent / "upstream/06_typing/receipt.json"
    release, upstream = tmp_path / "candidate", tmp_path / "upstream"
    absent = upstream / "05_normalization/alignment_block/donor_id=0112.parquet"
    item = {"name": "pln_blocks:0112", "path": str(absent), "exists": False, "recorded_sha256": None}
    project = tmp_path / "project.json"
    write(project, {"fixture": "synthetic provenance; no scientific run"})
    write(receipt, {"stage": "06_typing", "status": "complete", "inputs": [item], "outputs": [],
                    "project": str(project),
                    "settings": {"donors": ["0112"], "upstream_release": str(upstream)}})
    identity = audit.parent / "per_acquisition.json"
    write(identity, [{"donor_id": "0112", "identity_check": "pass", "embedded_pln": 0}])
    source = manifest.BASE / "core/multiplex5000/multiplex5000/stages/io60_mixed_protein_pilot.py"
    snapshot = tmp_path / "snapshot.json"
    write(snapshot, {"sources": [{"path": str(source.relative_to(manifest.BASE)), "sha256": manifest.sha(source)}]})
    launch = tmp_path / "attempt2_input_bindings.json"
    write(launch, {"schema": "protein-prelaunch-bindings/1", "attempt": 2,
                   "inputs": [item], "settings": manifest.read(receipt)["settings"],
                   "project": manifest.ref(project), "producer_snapshot": manifest.ref(snapshot)})
    report = folder / "remediation_checks.json"
    write(report, {"status": "pass", "receipt_required_unavailable_hashes": [], "receipt": manifest.ref(receipt),
          "verifier": manifest.ref(Path(__file__)), "generic_validator": manifest.ref(Path(__file__)),
          "optional_input_declaration": verify_optional_declaration(snapshot),
          "input_launch_bindings": verify_launch_bindings(launch, manifest.read(receipt), snapshot),
          "receipt_optional_absences": [{"input_name": "pln_blocks:0112", "path": str(absent), "exists": False,
              "donor_id": "0112", "embedded_pln_cells": 0, "status": "optional_absent_no_pln_cells",
              "receipt_input": item, "identity_audit": "release_audit/per_acquisition.json",
              "identity_audit_sha256": manifest.sha(identity), "evidence": "declared optional; audited zero pLN"}]})
    with audit.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["role", "path", "expected", "actual", "status"])
        writer.writeheader()
        writer.writerow({"role": "input", "path": str(absent), "expected": "", "actual": "", "status": "unavailable"})
    return audit, release, report, receipt, identity, absent


def test_declared_optional_absence_is_never_promoted_to_a_byte_hash(modules, tmp_path):
    manifest, _ = modules
    audit, release, report, _, _, _ = absence_fixture(tmp_path, manifest)
    row, = manifest.receipt_rows(audit, release, report)
    assert row["kind"] == "declared_optional_absence"
    assert row["sha256"] is row["bytes"] is None and row["exists"] is False
    assert manifest.verify_artifact(row, tmp_path) == "validated_declared_absence_not_byte_hash"
    row["kind"] = "audited_scientific_file"
    row["verification"] = "prior_full_byte_audit_and_receipt"
    with pytest.raises(ValueError, match="substituted"):
        manifest.verify_artifact(row, tmp_path)


@pytest.mark.parametrize("change", ["unclassified", "unknown", "required", "exists", "nonzero_pln", "identity_drift",
                                    "appeared", "dangling_link", "duplicate_name", "wrong_path", "forged_digest", "declaration"])
def test_unexplained_or_changed_absence_fails_closed(modules, tmp_path, change):
    manifest, _ = modules
    audit, release, report, receipt, identity, absent = absence_fixture(tmp_path, manifest)
    data, recorded = manifest.read(report), manifest.read(receipt)
    item = recorded["inputs"][0]
    entry = data["receipt_optional_absences"][0]
    if change == "unknown":
        item["name"] = entry["input_name"] = "required_measurement:0112"
    elif change == "required":
        item["required"] = True
    elif change == "exists":
        item["exists"] = True
    elif change == "nonzero_pln":
        write(identity, [{"donor_id": "0112", "identity_check": "pass", "embedded_pln": 1}])
        entry["identity_audit_sha256"] = manifest.sha(identity)
    elif change == "identity_drift":
        identity.write_text(identity.read_text() + "\n")
    elif change in {"appeared", "dangling_link"}:
        absent.parent.mkdir(parents=True)
        if change == "appeared":
            absent.write_bytes(b"new bytes after declared absence")
        else:
            absent.symlink_to(tmp_path / "missing-target")
    elif change == "duplicate_name":
        recorded["inputs"].append({**item, "path": str(tmp_path / "other.parquet")})
    elif change == "wrong_path":
        other = str(tmp_path / "wrong.parquet")
        audit.write_text(audit.read_text().replace(str(absent), other))
        item["path"] = entry["path"] = other
    elif change == "forged_digest":
        item["sha256"] = "1" * 64
    elif change == "declaration":
        data["optional_input_declaration"]["required"] = True
    entry["receipt_input"] = item
    write(receipt, recorded)
    data["receipt"] = manifest.ref(receipt)
    write(report, data)
    with pytest.raises(ValueError):
        manifest.receipt_rows(audit, release, None if change == "unclassified" else report)


def test_cached_audit_does_not_hide_new_optional_file(modules, tmp_path):
    manifest, _ = modules
    audit, release, report, _, _, absent = absence_fixture(tmp_path, manifest)
    row, = manifest.receipt_rows(audit, release, report)
    cache = {}
    manifest.verify_artifact(row, tmp_path, audit_cache=cache)
    absent.parent.mkdir(parents=True)
    absent.write_bytes(b"appeared since audit")
    with pytest.raises(ValueError, match="appeared or changed"):
        manifest.verify_artifact(row, tmp_path, audit_cache=cache)


def test_successful_protein_snapshot_comes_from_verified_report_not_original_default(modules, tmp_path):
    manifest, _ = modules
    _, _, report, _, _, _ = absence_fixture(tmp_path, manifest)
    original = tmp_path / "producer_snapshot.json"
    original.write_text("old failed attempt; not the successful run")
    chosen = manifest.verified_protein_snapshot(report)
    assert chosen == tmp_path / "snapshot.json" and chosen != original


@pytest.mark.parametrize("change", ["snapshot", "snapshot_digest", "receipt", "failed_report", "failed_receipt", "producer_binding"])
def test_protein_snapshot_requires_complete_report_and_exact_bindings(modules, tmp_path, change):
    manifest, _ = modules
    _, _, report, receipt, _, _ = absence_fixture(tmp_path, manifest)
    data = manifest.read(report)
    if change == "snapshot":
        path = Path(data["optional_input_declaration"]["snapshot_path"])
        path.write_text(path.read_text() + "\n")
    elif change == "snapshot_digest":
        data["optional_input_declaration"]["snapshot_sha256"] = "0" * 64
    elif change == "receipt":
        receipt.write_text(receipt.read_text() + "\n")
    elif change == "failed_report":
        data["status"] = "failed"
    elif change == "failed_receipt":
        value = manifest.read(receipt)
        value["status"] = "failed"
        write(receipt, value)
        data["receipt"] = manifest.ref(receipt)
    else:
        data["optional_input_declaration"]["sha256"] = "0" * 64
    write(report, data)
    with pytest.raises(ValueError):
        manifest.verified_protein_snapshot(report)


@pytest.mark.parametrize("change", ["missing", "artifact", "digest", "count", "boolean_count", "inputs",
                                    "duplicate_inputs", "settings", "project_bytes", "project_binding",
                                    "snapshot_binding", "snapshot_path", "receipt_input"])
def test_protein_run_requires_exact_prelaunch_input_closure(modules, tmp_path, change):
    manifest, _ = modules
    _, _, report, receipt, _, _ = absence_fixture(tmp_path, manifest)
    data = manifest.read(report)
    binding = data["input_launch_bindings"]
    path = Path(binding["path"])
    launch = manifest.read(path)
    if change == "missing":
        del data["input_launch_bindings"]
    elif change == "artifact":
        path.write_text(path.read_text() + "\n")
    elif change == "digest":
        binding["sha256"] = "0" * 64
    elif change == "count":
        binding["inputs"] += 1
    elif change == "boolean_count":
        binding["inputs"] = True
    elif change == "project_bytes":
        project = Path(launch["project"]["path"])
        project.write_text(project.read_text() + "\n")
    elif change == "receipt_input":
        consumed = manifest.read(receipt)
        consumed["inputs"][0]["exists"] = True
        write(receipt, consumed)
        data["receipt"] = manifest.ref(receipt)
    else:
        if change == "inputs":
            launch["inputs"][0]["path"] = str(tmp_path / "another.parquet")
        elif change == "duplicate_inputs":
            launch["inputs"].append(dict(launch["inputs"][0]))
        elif change == "settings":
            launch["settings"]["donors"].append("other")
        elif change == "project_binding":
            launch["project"]["sha256"] = "0" * 64
        elif change == "snapshot_binding":
            launch["producer_snapshot"]["sha256"] = "0" * 64
        else:
            launch["producer_snapshot"]["path"] = str(tmp_path / "earlier_snapshot.json")
        write(path, launch)
        binding["sha256"] = manifest.sha(path)
    write(report, data)
    with pytest.raises((ValueError, KeyError)):
        manifest.verified_protein_snapshot(report)


def source_stability_fixture(tmp_path, manifest):
    from verify_protein_remediation import verify_optional_declaration

    _, _, report_path, receipt_path, _, _ = absence_fixture(tmp_path, manifest)
    report = manifest.read(report_path)
    snapshot_path = Path(report["optional_input_declaration"]["snapshot_path"])
    snapshot = manifest.read(snapshot_path)
    stage_key = snapshot["sources"][0]["path"]
    stage_source = tmp_path / stage_key
    stage_source.parent.mkdir(parents=True)
    stage_source.write_text((manifest.BASE / stage_key).read_text())
    others = ["core/multiplex5000/multiplex5000/stages/typing.py",
              "core/multiplex5000/scripts/typing/refine.py",
              "core/multiplex5000/multiplex5000/tissues/helper.py"]
    for relative in others:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# Synthetic immutable producer fixture\n")
        snapshot["sources"].append({"path": relative, "sha256": manifest.sha(path)})
    write(snapshot_path, snapshot)
    report["optional_input_declaration"] = verify_optional_declaration(snapshot_path)
    launch_path = Path(report["input_launch_bindings"]["path"])
    launch = manifest.read(launch_path)
    launch["producer_snapshot"] = manifest.ref(snapshot_path)
    write(launch_path, launch)
    report["input_launch_bindings"]["sha256"] = manifest.sha(launch_path)
    producer = tmp_path / others[1]
    runner_path = receipt_path.parent / "pancreas/run.json"
    write(runner_path, {"status": "complete", "source_scripts": {str(producer): manifest.sha(producer)}})
    receipt = manifest.read(receipt_path)
    receipt["outputs"] = [{"path": "pancreas/run.json", "sha256": manifest.sha(runner_path)}]
    write(receipt_path, receipt)
    report["receipt"] = manifest.ref(receipt_path)
    write(report_path, report)
    checker = tmp_path / "scripts/evaluation/check_protein_producer_snapshot.py"
    checker.parent.mkdir(parents=True)
    checker.write_text("# Synthetic checker source identity; not a production execution claim\n")
    stability = {"status": "pass", "scope": "synthetic conservative static source closure", "limitations": "not a runtime trace",
                 "snapshot": {"path": str(snapshot_path.relative_to(tmp_path)), "sha256": manifest.sha(snapshot_path), "files": 4},
                 "runner_manifest": {"path": str(runner_path.relative_to(tmp_path)), "sha256": manifest.sha(runner_path)},
                 "closure": [{"path": item["path"], "launch_sha256": item["sha256"], "current_sha256": item["sha256"],
                              "unchanged": True, "status": "unchanged"} for item in snapshot["sources"]],
                 "runner_source_checks": [{"path": str(producer), "launch_sha256": manifest.sha(producer),
                                            "current_sha256": manifest.sha(producer), "unchanged": True}]}
    path = tmp_path / manifest.RUN_REL / "protein/producer_stability.json"
    write(path, stability)
    return report_path, path, receipt_path, runner_path, producer


def test_protein_stability_binds_successful_receipt_actual_snapshot_and_current_sources(modules, tmp_path):
    manifest, _ = modules
    report, stability, _, runner, _ = source_stability_fixture(tmp_path, manifest)
    result = manifest.verified_protein_source_stability(tmp_path, report)
    assert result["closure_sources"] == 4 and result["declared_runner_sources"] == 1
    assert result["audit"] == manifest.ref(stability) and result["runner_manifest"] == manifest.ref(runner)
    assert result["checker_current_source"] == manifest.ref(tmp_path / "scripts/evaluation/check_protein_producer_snapshot.py")
    assert "does not self-identify" in result["checker_binding_scope"]


@pytest.mark.parametrize("change", ["missing", "failed", "unavailable", "wrong_snapshot", "snapshot_hash", "snapshot_count",
                                    "running_manifest", "runner_hash", "runner_bytes", "duplicate_output", "source_drift",
                                    "unavailable_row", "duplicate_closure", "empty_closure", "missing_seed", "missing_runner",
                                    "duplicate_runner", "runner_digest", "outside_repo", "failed_receipt"])
def test_missing_stale_or_incomplete_protein_stability_fails_closed(modules, tmp_path, change):
    manifest, _ = modules
    report_path, path, receipt_path, runner_path, producer = source_stability_fixture(tmp_path, manifest)
    stability = manifest.read(path)
    if change == "missing":
        path.unlink()
    elif change == "failed":
        stability["status"] = "fail"
    elif change == "unavailable":
        stability["status"] = "pass_with_unavailable_static_candidates"
    elif change == "wrong_snapshot":
        stability["snapshot"]["path"] = "producer_snapshot_attempt1.json"
    elif change == "snapshot_hash":
        stability["snapshot"]["sha256"] = "0" * 64
    elif change == "snapshot_count":
        stability["snapshot"]["files"] += 1
    elif change == "running_manifest":
        stability["runner_manifest"]["path"] = str(tmp_path / ".06_typing.incomplete-other/pancreas/run.json")
    elif change == "runner_hash":
        stability["runner_manifest"]["sha256"] = "0" * 64
    elif change == "runner_bytes":
        runner_path.write_text(runner_path.read_text() + "\n")
    elif change in {"duplicate_output", "failed_receipt"}:
        receipt = manifest.read(receipt_path)
        if change == "duplicate_output":
            receipt["outputs"] *= 2
        else:
            receipt["status"] = "failed"
        write(receipt_path, receipt)
        report = manifest.read(report_path)
        report["receipt"] = manifest.ref(receipt_path)
        write(report_path, report)
    elif change == "source_drift":
        producer.write_text("# Changed after its source audit\n")
    elif change == "unavailable_row":
        stability["closure"][0].update(unchanged=None, status="not_snapshotted")
    elif change == "duplicate_closure":
        stability["closure"].append(stability["closure"][0])
    elif change == "empty_closure":
        stability["closure"] = []
    elif change == "missing_seed":
        stability["closure"] = [row for row in stability["closure"] if not row["path"].endswith("stages/typing.py")]
    elif change == "missing_runner":
        stability["runner_source_checks"] = []
    elif change == "duplicate_runner":
        stability["runner_source_checks"] *= 2
    elif change == "runner_digest":
        stability["runner_source_checks"][0]["launch_sha256"] = "0" * 64
    else:
        stability["closure"][0]["path"] = "../external.py"
    if change != "missing":
        write(path, stability)
    with pytest.raises((ValueError, OSError)):
        manifest.verified_protein_source_stability(tmp_path, report_path)
