#!/usr/bin/env python3
"""Typed bindings for current remediation sources and previously audited data.

Current source/evidence bytes are checked directly. Large scientific files retain
their prior full-byte audit and exact receipt bindings; this is explicitly not a
fresh hash of those files. Historical Git fallback belongs only to the separate
frozen-baseline verifier and is never used here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from collections import Counter
from importlib.metadata import version
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
RUN_REL = Path("outputs/blocker_remediation_20261102")
DOC_REL = Path("docs/evaluations/blocker_remediation_20261102")
FAILED_PROTEIN_ATTEMPT1_REL = RUN_REL / "protein/releases/protein-availability-candidate/.06_typing.incomplete-bac7a1c29fd6"
FAILED_PROTEIN_JSON = ("failure.json", "resolved_parameters.json", "pancreas/run.json")
SCHEMA = "portability-remediation-artifacts/1"
FRESH = {"current_source", "evidence_file"}
BOUND = {"audited_scientific_file", "audited_remote_image", "audited_recorded_input", "audited_content_recovery",
         "declared_optional_absence"}


def sha(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read(path: Path):
    return json.loads(path.read_text())


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def ref(path: Path) -> dict:
    return {"path": str(path.resolve()), "sha256": sha(path)}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def verified_protein_snapshot(report_path: Path) -> Path:
    """Bind the checked run's snapshot and full prelaunch inputs, never an earlier default."""
    from verify_protein_remediation import (
        verify_launch_bindings,
        verify_optional_declaration,
    )

    report = read(report_path)
    require(report["status"] == "pass", "protein report has not passed independent verification")
    receipt_ref = report["receipt"]
    receipt_path = Path(receipt_ref["path"])
    require(sha(receipt_path) == receipt_ref["sha256"], "protein report receipt changed")
    receipt = read(receipt_path)
    require(receipt.get("stage") == "06_typing" and receipt.get("status") == "complete",
            "protein snapshot requires an actual complete Stage06 receipt")
    declaration = report["optional_input_declaration"]
    path = Path(declaration["snapshot_path"])
    require(sha(path) == declaration["snapshot_sha256"], "protein producer snapshot changed")
    require(verify_optional_declaration(path) == declaration, "protein report/snapshot producer binding differs")
    launch = report["input_launch_bindings"]
    launch_path = Path(launch["path"])
    require(launch_path.is_absolute() and type(launch["inputs"]) is int and launch["inputs"] > 0,
            "protein prelaunch input verification needs an absolute artifact and exact input count")
    require(sha(launch_path) == launch["sha256"], "protein prelaunch input binding artifact changed")
    require(verify_launch_bindings(launch_path, receipt, path) == launch,
            "protein report/prelaunch input verification differs")
    return path


def verified_protein_source_stability(base: Path, report_path: Path) -> dict:
    """Require the successful attempt's producer audit, not a stale running report."""
    snapshot_path = verified_protein_snapshot(report_path)
    report = read(report_path)
    receipt_path = Path(report["receipt"]["path"])
    require(receipt_path.parent.name == "06_typing", "producer audit requires the published Stage06 directory")
    receipt = read(receipt_path)
    path = base / RUN_REL / "protein/producer_stability.json"
    stability = read(path)
    require(stability["status"] == "pass", "protein producer stability did not pass without unavailable sources")
    require(resolve(base, stability["snapshot"]["path"]).resolve() == snapshot_path.resolve()
            and stability["snapshot"]["sha256"] == sha(snapshot_path), "protein source audit uses a different launch snapshot")
    frozen_rows = read(snapshot_path)["sources"]
    frozen = {row["path"]: row["sha256"] for row in frozen_rows}
    require(len(frozen) == len(frozen_rows) == stability["snapshot"]["files"], "protein source snapshot count or identities differ")
    runner_path = receipt_path.parent / "pancreas/run.json"
    outputs = [row for row in receipt["outputs"] if row["path"] == "pancreas/run.json"]
    require(len(outputs) == 1 and outputs[0]["sha256"] == sha(runner_path)
            and runner_path.resolve().is_relative_to(receipt_path.parent.resolve()),
            "protein runner manifest is not uniquely bound by the successful receipt")
    require(resolve(base, stability["runner_manifest"]["path"]).resolve() == runner_path.resolve()
            and stability["runner_manifest"]["sha256"] == outputs[0]["sha256"],
            "protein source audit belongs to another or incomplete runner manifest")
    runner = read(runner_path)
    require(runner["status"] == "complete", "protein source audit runner is incomplete")
    closure = stability["closure"]
    paths = [row["path"] for row in closure]
    require(closure and len(set(paths)) == len(paths), "protein source audit has empty or duplicate closure identities")
    for row in closure:
        source = Path(row["path"])
        require(not source.is_absolute() and ".." not in source.parts
                and resolve(base, row["path"]).resolve().is_relative_to(base.resolve()),
                "protein closure source escapes the repository")
        require(row["status"] == "unchanged" and row["unchanged"] is True
                and row["launch_sha256"] == row["current_sha256"] == frozen.get(row["path"])
                and sha(base / source) == row["launch_sha256"],
                f"protein producer source is changed or unavailable: {row['path']}")
    declared = runner["source_scripts"]
    checks = stability["runner_source_checks"]
    names = [row["path"] for row in checks]
    require(declared and len(set(names)) == len(names) and set(names) == set(declared),
            "protein producer audit omits or duplicates declared runner sources")
    seeds = {"core/multiplex5000/multiplex5000/stages/io60_mixed_protein_pilot.py",
             "core/multiplex5000/multiplex5000/stages/typing.py"}
    for row in checks:
        source = Path(row["path"])
        require(source.is_absolute() and source.resolve().is_relative_to(base.resolve()),
                "protein runner source escapes the repository")
        relative = str(source.relative_to(base))
        seeds.add(relative)
        require(row["unchanged"] is True
                and row["launch_sha256"] == row["current_sha256"] == declared[row["path"]] == frozen.get(relative)
                and sha(source) == row["launch_sha256"], "protein declared producer differs from the launch snapshot")
    require(seeds <= set(paths), "protein conservative closure omits a required producer seed")
    checker = base / "scripts/evaluation/check_protein_producer_snapshot.py"
    return {"audit": ref(path), "snapshot": ref(snapshot_path), "runner_manifest": ref(runner_path),
            "closure_sources": len(closure), "declared_runner_sources": len(checks),
            "checker_current_source": ref(checker),
            "checker_binding_scope": "current delivery-time checker source; report schema does not self-identify its producer",
            "scope": stability["scope"], "limitations": stability["limitations"]}


def _classified_absences(classification: Path, audit: Path, receipt: Path) -> dict:
    """An audited zero-pLN section may lack only its declared optional block sidecar."""
    from verify_protein_remediation import verify_optional_declaration

    report = read(classification)
    require(report["status"] == "pass" and not report["receipt_required_unavailable_hashes"],
            "unexplained required/unclassified input absence")
    require(report["receipt"]["sha256"] == sha(receipt), "absence classification uses a different receipt")
    for key in ("verifier", "generic_validator"):
        value = report[key]
        require(sha(Path(value["path"])) == value["sha256"], "absence classifier source changed")
    declaration = report["optional_input_declaration"]
    require(verify_optional_declaration(Path(declaration["snapshot_path"])) == declaration,
            "optional input declaration or launch snapshot changed")
    identity_path = audit.parent / "per_acquisition.json"
    identities = read(identity_path)
    identity_sha = sha(identity_path)
    require(len({r["donor_id"] for r in identities}) == len(identities), "duplicate absence identity-audit donor")
    donors = {r["donor_id"]: r for r in identities}
    found = {}
    for row in report["receipt_optional_absences"]:
        donor, path = row["donor_id"], row["path"]
        require(path not in found, "duplicate classified optional absence")
        require(row["input_name"] == f"pln_blocks:{donor}" and row["exists"] is False
                and row["embedded_pln_cells"] == 0 and row["status"] == "optional_absent_no_pln_cells"
                and row.get("evidence"), "invalid classified optional absence")
        require(donor in donors and donors[donor]["identity_check"] == "pass"
                and donors[donor]["embedded_pln"] == 0, "missing block sidecar has no audited zero-pLN evidence")
        require(row["identity_audit"] == "release_audit/per_acquisition.json"
                and row["identity_audit_sha256"] == identity_sha, "unknown or changed absence identity audit")
        found[path] = row
    return found


def receipt_rows(audit: Path, release: Path, classification: Path | None = None) -> list[dict]:
    """Bind byte audits and explicitly classified optional absences without rereading data."""
    with audit.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    result = []
    absent = (_classified_absences(classification, audit, audit.parent / "upstream/06_typing/receipt.json")
              if classification else {})
    matched_absences = set()
    classification_ref = {"classification": ref(classification)} if classification else {}
    audit_ref = ref(audit)
    receipts = {}
    for index, row in enumerate(rows):
        require(row["role"] in {"input", "output"}, "scientific audit has an unknown role")
        stage = row.get("stage") or "06_typing"
        receipt_path = audit.parent / "upstream" / stage / "receipt.json"
        if receipt_path not in receipts:
            receipts[receipt_path] = read(receipt_path), ref(receipt_path)
        receipt, receipt_ref = receipts[receipt_path]
        require(receipt.get("status") == "complete" and receipt.get("stage") == stage,
                f"scientific receipt is incomplete: {receipt_path}")
        items = receipt["inputs" if row["role"] == "input" else "outputs"]
        matches = [item for item in items if str(Path(item["path"]) if row["role"] == "input"
                   else release / stage / item["path"]) == row["path"]]
        require(len(matches) == 1, f"audit is not uniquely bound by receipt: {row['path']}")
        item = matches[0]
        if row["status"] == "unavailable":
            known = absent.get(row["path"])
            require(known is not None and stage == "06_typing" and row["role"] == "input"
                    and not row["expected"] and not row["actual"], "unclassified scientific input absence")
            require(item == known["receipt_input"] and item.get("exists") is False
                    and item.get("sha256") is None and item.get("recorded_sha256") is None
                    and item.get("required", False) is False and item.get("bytes") is None,
                    "classified absence differs from exact declared optional receipt input")
            donor = known["donor_id"]
            require(item.get("name") == f"pln_blocks:{donor}" and donor in receipt["settings"]["donors"],
                    "optional absence is not a declared cohort block sidecar")
            require(sum(i.get("name") == item["name"] for i in receipt["inputs"]) == 1,
                    "optional absence has duplicate receipt names")
            upstream = Path(receipt["settings"]["upstream_release"])
            require(Path(row["path"]) == upstream / f"05_normalization/alignment_block/donor_id={donor}.parquet",
                    "optional absence has an incorrect declared path")
            path = Path(row["path"])
            require(not path.exists() and not path.is_symlink(), "declared absent optional input appeared or changed")
            matched_absences.add(row["path"])
            result.append({"kind": "declared_optional_absence", "verification": "receipt_declared_absence_and_zero_pln",
                           "path": row["path"], "sha256": None, "bytes": None, "exists": False,
                           "input_name": item["name"], "donor_id": donor, "embedded_pln_cells": 0,
                           "audit": audit_ref, "audit_row": index, "receipt": receipt_ref,
                           "release_root": str(release), "stage": stage, "role": row["role"],
                           **classification_ref, "freshly_hashed_by_delivery_verifier": False})
            continue
        require(row["status"] == "pass" and row["expected"] == row["actual"] == item.get("sha256") and row["expected"],
                f"scientific audit has an unsuccessful or unbound row: {audit}:{index}")
        result.append({"kind": "audited_scientific_file", "verification": "prior_full_byte_audit_and_receipt",
                       "path": row["path"], "sha256": row["expected"], "bytes": matches[0].get("bytes"),
                       "audit": audit_ref, "audit_row": index, "receipt": receipt_ref,
                       "release_root": str(release), "stage": stage, "role": row["role"],
                       **classification_ref, "freshly_hashed_by_delivery_verifier": False})
    require(matched_absences == set(absent), "classified absences do not match the complete original audit")
    return result


def verify_artifact(row: dict, base: Path = BASE, *, audit_cache: dict | None = None) -> str:
    """Verify the declared type strictly; no Git substitution for current bytes."""
    kind = row.get("kind")
    require(kind in FRESH | BOUND, f"unknown artifact binding type: {kind}")
    if kind in FRESH:
        require(row.get("verification") == "fresh_sha256", "fresh binding has incompatible verification mode")
        path = resolve(base, row["path"])
        require(path.is_file() and sha(path) == row["sha256"], f"current file changed or is missing: {path}")
        if "bytes" in row:
            require(path.stat().st_size == row["bytes"], f"current file size differs: {path}")
        return "fresh_sha256"
    require(row.get("freshly_hashed_by_delivery_verifier") is False,
            "prior scientific audit must not claim a fresh data hash")
    audit = row["audit"]
    audit_path = resolve(base, audit["path"])
    require(sha(audit_path) == audit["sha256"], f"scientific audit changed: {audit_path}")
    if kind in {"audited_scientific_file", "declared_optional_absence"}:
        mode = "prior_full_byte_audit_and_receipt" if kind == "audited_scientific_file" else "receipt_declared_absence_and_zero_pln"
        require(row["verification"] == mode, "wrong receipt binding mode")
        receipt = row["receipt"]
        require(sha(resolve(base, receipt["path"])) == receipt["sha256"], "scientific receipt changed")
        classification = row.get("classification")
        if classification:
            require(sha(resolve(base, classification["path"])) == classification["sha256"], "absence classification changed")
        key = (str(audit_path), audit["sha256"], row["release_root"], receipt["path"], receipt["sha256"],
               classification["sha256"] if classification else None)
        if audit_cache is None or key not in audit_cache:
            rows = receipt_rows(audit_path, Path(row["release_root"]),
                                resolve(base, classification["path"]) if classification else None)
            if audit_cache is not None:
                audit_cache[key] = rows
        else:
            rows = audit_cache[key]
        require(type(row["audit_row"]) is int and 0 <= row["audit_row"] < len(rows), "invalid scientific audit row index")
        expected = rows[row["audit_row"]]
        require(all(row.get(key) == expected.get(key) for key in
                    ["kind", "path", "sha256", "bytes", "stage", "role", "receipt", "classification",
                     "exists", "input_name", "donor_id", "embedded_pln_cells"]), "scientific audit row was substituted")
        if kind == "declared_optional_absence":
            path = Path(row["path"])
            require(not path.exists() and not path.is_symlink(), "declared absent optional input appeared or changed")
            return "validated_declared_absence_not_byte_hash"
    elif kind == "audited_remote_image":
        require(row["verification"] == "prior_full_image_audit", "wrong remote image binding mode")
        report = read(audit_path)
        matches = [item for item in report["acquisitions"] if item["acquisition_id"] == row["acquisition_id"]]
        require(len(matches) == 1 and report["complete"] and matches[0]["status"] == "validated",
                "remote image is not uniquely validated")
        require(matches[0]["image"] == row["path"] and matches[0]["image_sha256"] == row["sha256"]
                and matches[0]["image_bytes"] == row["bytes"], "remote image audit binding differs")
    elif kind == "audited_recorded_input":
        require(row["verification"] == "prior_full_current_bytes_with_declared_historical_algorithm",
                "wrong recorded-input binding mode")
        report = read(audit_path)
        require(report["complete"] and report["status"] == "pass", "recorded-input audit incomplete")
        matches = [item for item in report["inputs"] if item["name"] == row["input_name"]]
        require(len(matches) == 1 and matches[0]["status"] == "pass", "recorded input is not uniquely verified")
        require(all(matches[0][key] == row[key] for key in ["path", "sha256", "bytes", "recorded_hash_mode"]),
                "recorded-input algorithm or binding differs")
    else:
        require(row["verification"] == "prior_exact_release_content_recovery", "wrong recovery binding mode")
        report = read(audit_path)
        matches = [item for item in report["checks"] if item["path"] == row["path"]]
        require(len(matches) == 1 and not report["source_release_modified"], "recovered file is not uniquely bound")
        item = matches[0]
        require(item["status"] == "pass" and item["expected"] == item["actual"] == row["sha256"]
                and item["declared_path"] == row["declared_path"], "recovered content binding differs")
    return "bound_to_prior_audit_not_rehashed"


def _excluded(relative: Path) -> bool:
    parts = relative.parts
    return ("__pycache__" in parts or "releases" in parts or "v20_recovered" in parts
            or any(part.startswith(("pytest", "maldi_stage_tests")) for part in parts)
            or relative.suffix in {".pyc", ".lock", ".tmp"}
            or relative.name in {"observations.lock", "delivery_verification.json"})


def write_provenance(base: Path) -> dict:
    run, doc = base / RUN_REL, base / DOC_REL
    baseline = base / "docs/evaluations/dynamic_portability_20261102/pins.json"
    pins = read(baseline)
    checkout = base / "external/KINTSUGI"
    kintsugi = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
    require(kintsugi == pins["KINTSUGI"], "KINTSUGI pin changed from the original evaluation")
    pln = Path(read(run / "config/completion_summary.json")["scientific_output_root"])
    protein_report = run / "protein/verification/remediation_checks.json"
    protein_snapshot = verified_protein_snapshot(protein_report)
    protein_stability = verified_protein_source_stability(base, protein_report)
    runtime = {
        "CODEX_replay": run / "codex/exact_replay/replay.json",
        "MALDI_Xenium_sources": run / "hpg_maldi_xenium/results_relocated/source_validation.json",
        "protein": run / "protein/releases/protein-availability-candidate/06_typing/receipt.json",
        "protein_producer_snapshot": protein_snapshot,
        "protein_prelaunch_input_bindings": Path(read(protein_report)["input_launch_bindings"]["path"]),
        "protein_independent_verification": protein_report,
        "protein_source_stability": Path(protein_stability["audit"]["path"]),
        "dedicated_pLN": pln / "06_typing/receipt.json",
        "pLN_source_stability": run / "config/pln_source_stability.json",
    }
    value = {"schema": "portability-remediation-provenance/1",
             "launch_base": "610a6a1140427ebfd67b8bac86bbaa69b46ac739",
             "launch_branch": "fix/portability-blockers-20261102",
             "delivery_worktree_git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=base, text=True).strip(),
             "delivery_commit_is_execution_producer": False,
             "producer_policy": "Execution identities come from each actual receipt and immutable producer snapshot; later delivery commits are not execution producers.",
             "original_scientific_pins": {"binding": ref(baseline), "pins": pins},
             "KINTSUGI_current": kintsugi,
             "runtime_evidence": {name: ref(path) for name, path in runtime.items()},
             "protein_source_validation": protein_stability,
             "failed_attempt_evidence": {"protein_original_launch_snapshot": ref(run / "protein/producer_snapshot.json")},
             "delivery_verification_environment": {
                 "scope": "current delivery checks only; not substituted for producer environments",
                 "executable": sys.executable, "python": platform.python_version(), "platform": platform.platform(),
                 "packages": {name: version(name) for name in ["numpy", "pandas", "scipy", "pyarrow", "pytest"]}}}
    (doc / "provenance.json").write_text(json.dumps(value, indent=2) + "\n")
    return value


def catalog(base: Path, destination: Path) -> tuple[set[Path], set[Path]]:
    """Complete current source/evidence inventory, also required during verification."""
    run, doc = base / RUN_REL, base / DOC_REL
    files = {p for p in run.rglob("*") if p.is_file() and not _excluded(p.relative_to(run))}
    # Selected small failed-attempt records are evidence, never audited scientific outputs.
    for relative in FAILED_PROTEIN_JSON:
        path = base / FAILED_PROTEIN_ATTEMPT1_REL / relative
        require(path.is_file(), f"retained failed-attempt evidence is missing: {path}")
        files.add(path)
    files |= {p for p in doc.rglob("*") if p.is_file() and p != destination}
    files.add(base / "docs/PORTABILITY_REMEDIATION_USER_GUIDE.md")
    recovery_report = read(run / "codex/v20_final_verification/summary.json")
    files.add(resolve(base, recovery_report["source_release"]))
    files.add(resolve(base, recovery_report["recovery_manifest"]))
    names = subprocess.check_output(["git", "ls-files", "--cached", "--others", "--exclude-standard"],
                                    cwd=base, text=True).splitlines()
    sources = {base / name for name in names
               if name.startswith(("core/multiplex5000/", "core/panels/", "scripts/evaluation/"))
               and Path(name).suffix in {".py", ".md", ".yaml", ".yml", ".json", ".toml", ".csv"}}
    sources |= {base / path for path in [".gitignore", "README.md", "CLAUDE.md", "scripts/CLAUDE.md", "docs/PORTABILITY_EVALUATION_USER_GUIDE.md",
                                        "docs/PORTABILITY_REMEDIATION_USER_GUIDE.md"]}
    sources = {path for path in sources if path.is_file()}
    files |= sources
    return {p for p in files if p.is_file() and p != destination}, sources


def build(base: Path = BASE, destination: Path | None = None) -> dict:
    base = base.resolve()
    run, doc = base / RUN_REL, base / DOC_REL
    destination = destination or doc / "artifact_manifest.json"
    write_provenance(base)
    files, sources = catalog(base, destination)
    recovery_report = read(run / "codex/v20_final_verification/summary.json")
    artifacts = []
    for path in sorted(files):
        if not path.is_file() or path == destination:
            continue
        require(path.stat().st_size < 100_000_000,
                f"large file needs an explicit audited binding instead of a fresh delivery hash: {path}")
        artifacts.append({"kind": "current_source" if path in sources else "evidence_file",
                          "verification": "fresh_sha256", **ref(path), "bytes": path.stat().st_size})
    pln_root = Path(read(run / "config/completion_summary.json")["scientific_output_root"])
    artifacts += receipt_rows(run / "config/pln_full19_validation/receipt_checks.csv", pln_root)
    protein = run / "protein/verification/release_audit"
    protein_root = Path(read(protein / "summary.json")["root"])
    artifacts += receipt_rows(protein / "receipt_checks.csv", protein_root, protein.parent / "remediation_checks.json")
    audit = run / "hpg_maldi_xenium/results_relocated/source_validation.json"
    for row in read(audit)["acquisitions"]:
        if row["status"] == "validated":
            artifacts.append({"kind": "audited_remote_image", "verification": "prior_full_image_audit",
                              "path": row["image"], "sha256": row["image_sha256"], "bytes": row["image_bytes"],
                              "acquisition_id": row["acquisition_id"], "audit": ref(audit),
                              "freshly_hashed_by_delivery_verifier": False})
    audit = run / "hpg_maldi_xenium/spleen_algorithm_verification/verification.json"
    for row in read(audit)["inputs"]:
        artifacts.append({"kind": "audited_recorded_input",
                          "verification": "prior_full_current_bytes_with_declared_historical_algorithm",
                          "path": row["path"], "sha256": row["sha256"], "bytes": row["bytes"],
                          "recorded_hash_mode": row["recorded_hash_mode"], "input_name": row["name"],
                          "audit": ref(audit), "freshly_hashed_by_delivery_verifier": False})
    audit = run / "codex/v20_final_verification/summary.json"
    for row in recovery_report["checks"]:
        artifacts.append({"kind": "audited_content_recovery", "verification": "prior_exact_release_content_recovery",
                          "path": row["path"], "declared_path": row["declared_path"], "sha256": row["actual"],
                          "audit": ref(audit), "freshly_hashed_by_delivery_verifier": False})
    result = {"schema": SCHEMA, "repository": str(base),
              "current_source_git_base": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=base, text=True).strip(),
              "current_source_policy": "Exact current working-file SHA256 required; no historical Git fallback.",
              "baseline_policy": "Only verify_portability_delivery.py may use its pinned historical Git source fallback.",
              "artifacts": artifacts, "counts": dict(Counter(row["kind"] for row in artifacts)),
              "scope": "Current source and small evidence freshly hashed, including three selected failed-attempt metadata JSONs. Successful scientific files use completed byte audits and receipts; optional absences are separately classified. Large failed-attempt outputs are retained without a full audit or successful receipt coverage.",
              "exclusions": ["successful large scientific outputs are represented only by their explicit full-cohort audit/receipt bindings",
                             "large failed-attempt outputs remain in place but are not covered by successful cohort audits; only failure.json, resolved_parameters.json and pancreas/run.json are freshly hashed metadata evidence",
                             "unbound raw and scratch trees are not asserted to have full-cohort audit coverage",
                             "Python caches, locks and duplicate test-stage scratch trees",
                             "self and final verifier result avoid recursive hashes",
                             "v20 recovered scientific tables retain the independent 24-binding verification report"]}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE)
    args = parser.parse_args()
    result = build(args.base)
    print(json.dumps(result["counts"], indent=2))


if __name__ == "__main__":
    main()
