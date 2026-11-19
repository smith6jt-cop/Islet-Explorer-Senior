"""A successful validator process must not conceal failed or incomplete evidence."""
import importlib.util
import json
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "verify_protein_remediation", Path(__file__).parents[1] / "verify_protein_remediation.py")
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


def audit(tmp_path, *, hashes=None, identities=None, rows=None, count=2):
    summary = {"hash_failures": hashes or [], "identity_failures": identities or [],
               "unavailable_hashes": [],
               "acquisitions": count}
    if rows is None:
        rows = [{"donor_id": donor, "identity_check": "pass"} for donor in ["0112", "6347"]]
    (tmp_path / "summary.json").write_text(json.dumps(summary))
    (tmp_path / "per_acquisition.json").write_text(json.dumps(rows))


@pytest.mark.parametrize("failure", ["hash", "identity"])
def test_zero_exit_with_recorded_failure_is_rejected(tmp_path, failure):
    # These are the artifacts produced even when the existing generic validator
    # returns 0. The actual recorded evidence must control acceptance.
    audit(tmp_path, hashes=[{"path": "changed.parquet"}] if failure == "hash" else [],
          identities=[{"donor_id": "0112"}] if failure == "identity" else [])
    with pytest.raises(ValueError, match="verification failed"):
        VERIFY.verify_audit_results(tmp_path, ["0112", "6347"], {"inputs": []})


@pytest.mark.parametrize("rows", [
    [{"donor_id": "0112", "identity_check": "pass"}],
    [{"donor_id": "0112", "identity_check": "pass"}] * 2,
    [{"donor_id": "0112", "identity_check": "pass"}, {"donor_id": "other", "identity_check": "pass"}],
])
def test_missing_duplicate_or_wrong_acquisition_cannot_pass(tmp_path, rows):
    audit(tmp_path, rows=rows)
    with pytest.raises(ValueError, match="identity roster"):
        VERIFY.verify_audit_results(tmp_path, ["0112", "6347"], {"inputs": []})


@pytest.mark.parametrize("state", ["blocked", "fail", None])
def test_every_declared_acquisition_requires_explicit_pass(tmp_path, state):
    audit(tmp_path, rows=[{"donor_id": "0112", "identity_check": "pass"},
                          {"donor_id": "6347", "identity_check": state}])
    with pytest.raises(ValueError, match="explicit passing identity"):
        VERIFY.verify_audit_results(tmp_path, ["0112", "6347"], {"inputs": []})


def optional_fixture(tmp_path, *, pLN_cells=0):
    audit(tmp_path, rows=[{"donor_id": "0112", "identity_check": "pass", "embedded_pln": pLN_cells},
                          {"donor_id": "6347", "identity_check": "pass", "embedded_pln": 0}])
    path = str(tmp_path / "05_normalization/alignment_block/donor_id=0112.parquet")
    Path(path).parent.mkdir(parents=True)
    summary = json.loads((tmp_path / "summary.json").read_text())
    summary["unavailable_hashes"] = [{"path": path, "role": "input", "expected": None, "actual": None}]
    (tmp_path / "summary.json").write_text(json.dumps(summary))
    return {"settings": {"upstream_release": str(tmp_path)}, "inputs": [{"name": "pln_blocks:0112", "path": path, "exists": False, "recorded_sha256": None}]}


def test_declared_absent_blocks_with_zero_pln_cells_are_distinct_from_hash_pass(tmp_path):
    receipt = optional_fixture(tmp_path)
    summary, rows, optional = VERIFY.verify_audit_results(tmp_path, ["0112", "6347"], receipt)
    assert len(rows) == 2 and len(summary["unavailable_hashes"]) == 1
    assert optional[0]["status"] == "optional_absent_no_pln_cells"
    assert optional[0]["exists"] is False
    assert optional[0]["receipt_input"] == receipt["inputs"][0]


@pytest.mark.parametrize("mutation", ["required", "unknown", "duplicate", "duplicate_name", "wrong_upstream", "exists", "sha", "file", "symlink", "output", "missing_count", "positive_count"])
def test_unexplained_or_scientifically_used_missing_binding_is_rejected(tmp_path, mutation):
    receipt = optional_fixture(tmp_path)
    item = receipt["inputs"][0]
    if mutation == "required":
        item["name"] = "aligned:0112"
    elif mutation == "unknown":
        receipt["inputs"] = []
    elif mutation == "duplicate":
        receipt["inputs"] = [item, item]
    elif mutation == "duplicate_name":
        receipt["inputs"].append({**item, "path": str(tmp_path / "other.parquet")})
    elif mutation == "wrong_upstream":
        receipt["settings"]["upstream_release"] = str(tmp_path / "wrong-release")
    elif mutation == "exists":
        item["exists"] = True
    elif mutation == "sha":
        item["sha256"] = "a" * 64
    elif mutation == "file":
        Path(item["path"]).write_text("new file")
    elif mutation == "symlink":
        Path(item["path"]).symlink_to(tmp_path / "absent-target")
    elif mutation == "output":
        path = tmp_path / "summary.json"
        summary = json.loads(path.read_text())
        summary["unavailable_hashes"][0]["role"] = "output"
        path.write_text(json.dumps(summary))
    else:
        path = tmp_path / "per_acquisition.json"
        rows = json.loads(path.read_text())
        if mutation == "missing_count":
            rows[0].pop("embedded_pln")
        else:
            rows[0]["embedded_pln"] = 1
        path.write_text(json.dumps(rows))
    with pytest.raises(ValueError, match="unclassified unavailable receipt binding"):
        VERIFY.verify_audit_results(tmp_path, ["0112", "6347"], receipt)


def test_complete_hashed_evidence_needs_no_optional_absence(tmp_path):
    audit(tmp_path)
    summary, rows, optional = VERIFY.verify_audit_results(tmp_path, ["0112", "6347"], {"inputs": []})
    assert len(rows) == 2 and not summary["unavailable_hashes"] and optional == []


@pytest.mark.parametrize("required", ["True", "None", "0"])
def test_optional_declaration_must_be_explicit_false(required):
    source = ("def inputs():\n    return Input(f'pln_blocks:{donor}', upstream / "
              "f'05_normalization/alignment_block/donor_id={donor}.parquet', required=" + required + ")")
    with pytest.raises(ValueError, match="optional pLN input contract"):
        VERIFY.check_optional_declaration(source)


def test_optional_declaration_requires_exact_path():
    source = "Input(f'pln_blocks:{donor}', upstream / f'other/{donor}.parquet', required=False)"
    with pytest.raises(ValueError, match="optional pLN input contract"):
        VERIFY.check_optional_declaration(source)


def test_exact_declared_optional_input_contract():
    VERIFY.check_optional_declaration("Input(f'pln_blocks:{donor}', upstream / f'05_normalization/alignment_block/donor_id={donor}.parquet', required=False)")


def launch_fixture(tmp_path):
    config, snapshot = tmp_path / "project.yaml", tmp_path / "snapshot.json"
    config.write_text("id: candidate\n")
    snapshot.write_text('{"sources": []}\n')
    receipt = {"project": str(config), "settings": {"donors": ["0112"]},
               "inputs": [{"name": "aligned:0112", "path": "/bound/input", "exists": True,
                           "sha256": "a" * 64, "bytes": 100, "mtime_ns": 123}]}
    launch = {"project": {"path": str(config), "sha256": VERIFY.sha(config)},
              "producer_snapshot": {"path": str(snapshot), "sha256": VERIFY.sha(snapshot)},
              "inputs": receipt["inputs"], "settings": receipt["settings"]}
    path = tmp_path / "launch.json"
    path.write_text(json.dumps(launch))
    return path, receipt, snapshot


@pytest.mark.parametrize("mutation", ["input_hash", "duplicate", "missing", "settings", "config", "snapshot"])
def test_execution_receipt_cannot_diverge_from_prelaunch_bindings(tmp_path, mutation):
    path, receipt, snapshot = launch_fixture(tmp_path)
    if mutation == "input_hash":
        receipt["inputs"][0]["sha256"] = "b" * 64
    elif mutation == "duplicate":
        receipt["inputs"] *= 2
    elif mutation == "missing":
        receipt["inputs"] = []
    elif mutation == "settings":
        receipt["settings"]["donors"] = ["other"]
    elif mutation == "config":
        Path(receipt["project"]).write_text("id: changed\n")
    else:
        snapshot.write_text('{"sources": ["changed"]}\n')
    with pytest.raises(ValueError):
        VERIFY.verify_launch_bindings(path, receipt, snapshot)


def test_exact_execution_and_prelaunch_bindings_pass(tmp_path):
    path, receipt, snapshot = launch_fixture(tmp_path)
    result = VERIFY.verify_launch_bindings(path, receipt, snapshot)
    assert result == {"path": str(path), "sha256": VERIFY.sha(path), "inputs": 1}
