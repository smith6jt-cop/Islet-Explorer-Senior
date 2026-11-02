#!/usr/bin/env python3
"""Fail-closed final checks for the remediation delivery, separate from its baseline.

This does not recompute scientific typing, approve gates or reread large audited
images. It checks complete full-cohort reports, exact artifact/receipt bindings,
current sources, real-versus-synthetic scope and project-only observation history.
"""
from __future__ import annotations

import argparse
import json
import re
import sqlite3
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from manifest_remediation_artifacts import (
    BASE,
    DOC_REL,
    FAILED_PROTEIN_ATTEMPT1_REL,
    RUN_REL,
    build,
    catalog,
    read,
    receipt_rows,
    ref,
    require,
    resolve,
    sha,
    verified_protein_source_stability,
    verify_artifact,
)

AXES = {"code_availability", "valid_inputs", "scientific_specification", "executable_routing",
        "computational_validation", "biological_review"}
DISPOSITIONS = {"completed", "blocked", "excluded", "unavailable"}


def references(value):
    if isinstance(value, dict):
        if isinstance(value.get("path"), str) and isinstance(value.get("sha256"), str):
            yield value
        for child in value.values():
            yield from references(child)
    elif isinstance(value, list):
        for child in value:
            yield from references(child)


def check_refs(value, base: Path) -> int:
    found = {(item["path"], item["sha256"]) for item in references(value)}
    for path, expected in found:
        require(sha(resolve(base, path)) == expected, f"evidence reference changed: {path}")
    return len(found)


def indexed(rows: list[dict], key: str) -> dict:
    values = [row[key] for row in rows]
    require(len(values) == len(set(values)), f"duplicate {key} in cohort evidence")
    return {row[key]: row for row in rows}


def check_accounting(base: Path) -> dict:
    doc = base / DOC_REL
    baseline = read(base / "docs/evaluations/dynamic_portability_20261102/acquisition_dispositions.json")
    original = indexed(baseline, "acquisition_id")
    ledger, readiness = read(doc / "acquisitions.json"), read(doc / "readiness.json")
    current = indexed(ledger["acquisitions"], "acquisition_id")
    require(len(original) == len(current) == ledger["declared_acquisitions"] == 201 and set(current) == set(original),
            "final ledger does not account for the exact 201-acquisition roster")
    for aid, row in current.items():
        require(row["baseline"] == original[aid], f"frozen baseline row was rewritten: {aid}")
        value = row["current"]
        require(value["disposition"] in DISPOSITIONS and value.get("reason"), f"unfinished disposition: {aid}")
        require(not value.get("production_release_activated"), "delivery activated a production release")
        require(value.get("scientific_approval") in {None, False, "not_granted_by_evaluation"}, "scientific approval asserted")
    counts = dict(Counter(row["current"]["disposition"] for row in current.values()))
    require(counts == ledger["disposition_counts"], "current acquisition counts do not reconcile")
    require(ledger["baseline_roster_unchanged"] is True, "baseline preservation not declared")
    for row in readiness["implementation_items"]:
        require(row.get("status") and "pending" not in row["status"].lower(), f"implementation still pending: {row['id']}")
    require(all(set(row["axes"]) == AXES for row in readiness["routes"]), "readiness axes are incomplete")
    require(not readiness["scientific_approvals_added"] and not readiness["production_release_activated"],
            "readiness grants scientific or release approval")
    bindings = check_refs(ledger, base) + check_refs(readiness, base)
    text = (doc / "README.md").read_text()
    for status, count in counts.items():
        require(re.search(rf"\b{count}\s*(?:\*\*)?\s*{status}\b", text), f"README lacks reconciled {count} {status}")
    require(not re.search(r"\b25\s*(?:\*\*)?\s*pending\b", text, re.IGNORECASE), "README still reports pending full25 execution")
    return {"acquisitions": 201, "counts": counts, "evidence_references": bindings}


def check_pln(base: Path) -> dict:
    run = base / RUN_REL
    summary = read(run / "config/completion_summary.json")
    folder = run / "config/pln_full19_validation"
    sections = indexed(read(folder / "per_section.json"), "acquisition_id")
    baseline = read(base / "docs/evaluations/dynamic_portability_20261102/acquisition_dispositions.json")
    expected = {row["acquisition_id"] for row in baseline if row["cohort"] == "io60"
                and row["modality"] == "phenocycler" and row["tissue"] == "pancreatic_lymph_node"}
    require(summary["status"] == "complete" and summary["sections"] == len(sections) == 19
            and set(sections) == expected, "dedicated pLN full19 cohort incomplete")
    require(sum(row["cells"] for row in sections.values()) == summary["cells"] == 10_505_964, "pLN cell count mismatch")
    require(sum(row["eligible"] for row in sections.values()) == summary["eligible"] == 9_883_811, "pLN eligibility mismatch")
    require(all(row["identity_check"] == "pass" and row["nucleus_check"] == "exact raw admission; ineligible calls unavailable" for row in sections.values()),
            "pLN identity or nucleus verification failed")
    inspections = read(folder / "inspection_ledger.json")["sections"]
    require(set(indexed(inspections, "acquisition_id")) == expected, "pLN inspection roster incomplete")
    for row in inspections:
        require(row["visual_status"] == "inspected", "pLN section was not inspected")
        require(sha(resolve(base, row["plot"])) == row["plot_sha256"], "pLN plot changed")
        require(sha(resolve(base, row["viewed_in_contact_sheet"])) == row["contact_sheet_sha256"], "pLN viewed panel changed")
    samples = read(folder / "sampled_evidence_checks.json")
    require(samples["status"] == "pass" and samples["cells"] == 1900 and samples["marker_states_checked"] == 104500,
            "pLN independent sampled-state checks incomplete")
    require(len(samples["sections"]) == 19 and all(row["status"] == "pass" for row in samples["sections"]), "pLN sample cohort incomplete")
    release = Path(summary["scientific_output_root"])
    bindings = receipt_rows(folder / "receipt_checks.csv", release)
    require(len(bindings) == 204 and summary["receipt_hash_checks"] == {"fail": 0, "pass": 204, "unavailable": 0},
            "pLN complete input/output audit missing")
    require(not summary["typing_accepted"] and not (release / "06_typing/gate/accepted.json").exists(), "pLN scientific gate accepted")
    stability = read(run / "config/pln_source_stability.json")
    require(stability["status"] == "pass", "pLN execution sources changed")
    for item in stability["sources"]:
        require(item["unchanged"] and item["sha256"] == item["post_run_sha256"]
                and sha(Path(item["snapshot"])) == item["sha256"], "pLN runtime source snapshot mismatch")
    return {"sections": 19, "cells": summary["cells"], "audited_bindings_not_rehashed": len(bindings), "inspected": 19}


def check_protein(base: Path) -> dict:
    folder = base / RUN_REL / "protein/verification"
    report = read(folder / "remediation_checks.json")
    audit = read(folder / "release_audit/summary.json")
    require(report["status"] == "pass" and report["acquisitions"] == audit["acquisitions"] == 25,
            "full25 protein candidate verification incomplete")
    source_stability = verified_protein_source_stability(base, folder / "remediation_checks.json")
    rows = indexed(report["rows"], "donor_id")
    baseline = read(resolve(base, report["input_baseline"]["path"]))
    require(set(rows) == {row["donor_id"] for row in baseline["rows"]} and len(rows) == 25,
            "protein donor roster does not match exact baseline inputs")
    require(all(row["identity_check"] == "pass" and row["primary_calls_unchanged"] for row in rows.values()),
            "protein identities or primary calls changed")
    require(sum(row["primary_cells"] for row in rows.values()) == report["primary_cells"], "protein primary cell counts mismatch")
    require(not audit["hash_failures"] and not audit["identity_failures"]
            and not report["receipt_required_unavailable_hashes"], "protein full-release audit has unexplained failures")
    require(report["receipt_unavailable_hashes"] == audit["unavailable_hashes"], "raw unavailable audit rows were rewritten")
    check_refs({"receipt": report["receipt"], "baseline": report["input_baseline"]}, base)
    release = Path(audit["root"])
    require(not (release / "06_typing/gate/accepted.json").exists(), "protein Stage06 was accepted")
    bindings = receipt_rows(folder / "release_audit/receipt_checks.csv", release, folder / "remediation_checks.json")
    absent = [row for row in bindings if row["kind"] == "declared_optional_absence"]
    require(len(absent) == len(audit["unavailable_hashes"]) == len(report["receipt_optional_absences"]) == 15,
            "full25 optional-absence accounting differs from actual declared cohort")
    per = read(folder / "release_audit/per_acquisition.json")
    require(len(per) == 25 and all(row["identity_check"] == "pass" for row in per), "protein full per-acquisition audit missing")
    inspection = read(folder / "inspection_ledger.json")
    entries = inspection["sections"] if isinstance(inspection, dict) else inspection
    require(len(indexed(entries, "donor_id")) == 25 and {r["donor_id"] for r in entries} == set(rows),
            "protein spatial inspection does not cover all25 donors")
    for row in entries:
        require(row.get("visual_status", row.get("visual_review")) in {"inspected", "completed"}, "protein section uninspected")
        check_refs(row, base)
        if "plot" in row and "plot_sha256" in row:
            require(sha(resolve(base, row["plot"])) == row["plot_sha256"], "protein plot changed")
    return {"acquisitions": 25, "primary_cells": report["primary_cells"],
            "producer_source_stability": source_stability,
            "exact_prelaunch_input_bindings": report["input_launch_bindings"]["inputs"],
            "audited_byte_bindings_not_rehashed": len(bindings) - len(absent),
            "validated_declared_optional_absences_not_hashes": len(absent), "inspected": 25}


def check_maldi(base: Path) -> dict:
    run = base / RUN_REL
    source_path = run / "hpg_maldi_xenium/results_relocated/source_validation.json"
    source = read(source_path)
    rows = indexed(source["acquisitions"], "acquisition_id")
    baseline = read(base / "docs/evaluations/dynamic_portability_20261102/acquisition_dispositions.json")
    expected = {r["acquisition_id"] for r in baseline if r["cohort"] == "io60" and r["modality"] == "xenium"}
    require(source["complete"] and source["declared_acquisitions"] == len(rows) == 65 and set(rows) == expected,
            "Xenium morphology audit does not cover all65 exact acquisitions")
    require(Counter(r["status"] for r in rows.values()) == {"validated": 62, "unavailable": 3}, "unexpected morphology source dispositions")
    for row in rows.values():
        require(row["accepted_registration"] is False and row["neighbour_assignment"] is None, "source audit invented registration")
        if row["status"] == "unavailable":
            require(row["condition"] == "damaged" and row["reason"] == "recorded damaged Xenium tissue", "unexplained unavailable source")
        else:
            mask = source_path.parent / "masks" / Path(row["mask_artifact"]).name
            require(sha(mask) == row["mask_artifact_sha256"], "morphology mask artifact changed")
    for filename, key in [("producer_maldi_serial_registration.py", "producer_sha256"),
                          ("producer_validate_maldi_xenium_sources.py", "validation_script_sha256")]:
        require(sha(source_path.parent / filename) == source[key], "source audit producer snapshot changed")
    review = read(run / "maldi_xenium_review.json")
    require(review["source_report_sha256"] == sha(source_path) and review["visual_sections_inspected"] == 62,
            "morphology inspection is stale or incomplete")
    require(set(indexed(review["acquisitions"], "acquisition_id")) == expected, "morphology review roster differs")
    for row in review["acquisitions"]:
        require(not row["registration_accepted"] and not row["mask_accepted"] and row["neighbour_assignment"] is None,
                "morphology review granted scientific approval")
        if row["source_status"] == "validated":
            require(row["visual_review"] == "completed", "available morphology was not inspected")
            if "figure_sha256" in row:
                require(sha(resolve(base, row["figure"])) == row["figure_sha256"], "viewed morphology panel changed")
    guard = read(run / "final_guard_review/source_binding_review.json")
    require(guard["status"] == "pass" and guard["unique_experiment_bindings"] == guard["unique_source_index_bindings"] == 65
            and guard["prior_source_experiment_matches"] == guard["fresh_local_mask_hashes"] == 62,
            "final unique-binding guard not verified for full65")
    check_refs(guard, base)
    scaling_path = run / "maldi_scaling/manifest.json"
    scaling, inspected = read(scaling_path), read(run / "maldi_scaling/inspection.json")
    maldi_ids = {r["acquisition_id"] for r in baseline if r["modality"] == "maldi"}
    plotted, reviewed = indexed(scaling["sections"], "acquisition_id"), indexed(inspected["acquisitions"], "acquisition_id")
    require(set(plotted) == set(reviewed) == maldi_ids and len(plotted) == 32, "MALDI plots omit declared sections")
    require(scaling["spots"] == inspected["spots"] == 1_666_450 and inspected["manifest_sha256"] == sha(scaling_path), "scaling review is stale")
    for aid, row in plotted.items():
        require(sha(scaling_path.parent / row["plot"]) == row["sha256"] == reviewed[aid]["plot_sha256"], "MALDI scaling plot changed")
        view = reviewed[aid]
        require(view["visual_review"] == "completed" and view["registration_acceptance"] is None, "MALDI plot not reviewed or registration invented")
        require(sha(resolve(base, view["viewed_figure"])) == view["viewed_figure_sha256"], "viewed scaling figure changed")
    algorithm = read(run / "hpg_maldi_xenium/spleen_algorithm_verification/verification.json")
    require(algorithm["complete"] and algorithm["status"] == "pass" and algorithm["declared"] == algorithm["passed"] == 10,
            "full spleen current-byte/historical-algorithm audit incomplete")
    require(Counter(row["recorded_hash_mode"] for row in algorithm["inputs"]) == {"full": 8, "sampled": 2}, "historical hash algorithms were conflated")
    require(all(row["status"] == "pass" and row["hash_mode"] == "full" and row["stable_during_hash"]
                and row["size_matches_receipt"] and row["observed_recorded_mode_hash"] == row["recorded_sha256"]
                for row in algorithm["inputs"]), "recorded input verification failed")
    return {"declared_sources": 65, "validated_sources": 62, "damaged": 3, "scaling_plots_inspected": 32,
            "large_input_bytes_bound_not_rehashed": algorithm["bytes_hashed"]}


def check_geometry_codex(base: Path) -> dict:
    run = base / RUN_REL
    geom = read(run / "geom1/completion_summary.json")
    require(geom["tests"] == {"tests": 98, "failures": 0, "errors": 0, "skipped": 0}, "same-slide validation incomplete")
    real = read(run / "geom1/full_cohort_preflight/ledger.json")
    synthetic = read(run / "geom1/synthetic_batch/ledger/ledger.json")
    require(real["all_pairs_accounted"] and real["declared_pairs"] == real["accounted_pairs"] == 12
            and real["counts"] == {"blocked": 8, "unavailable": 4} and not real["full_cohort_candidate_completion"], "real pair blockers were hidden")
    require(synthetic["counts"] == {"completed": 2, "blocked": 1} and not synthetic["scientific_approval"], "synthetic batch result mismatched")
    require(all(not row.get("scientific_approval") for row in real["pairs"]), "real geometric approval invented")
    codex = read(run / "codex/completion_summary.json")
    replay = read(run / "codex/exact_replay/replay.json")
    require(codex["status"] == "complete" and codex["images"] == replay["expected_images"] == len(replay["images"]) == 20,
            "CODEX full20 replay incomplete")
    require(codex["exact_assertions_passed"] == replay["passed"] == 60 and replay["failed"] == codex["exact_assertions_failed"] == 0,
            "CODEX exact checks failed")
    require(all(len(row["checks"]) == 3 and all(c["status"] == "pass" for c in row["checks"]) for row in replay["images"]), "CODEX per-image assertions missing")
    require(not codex["tolerance_changed"] and not codex["original_fixture_or_historical_failures_modified"]
            and replay["unchanged_scientific_inputs_and_baselines"], "CODEX historical criteria or evidence changed")
    v20 = read(run / "codex/v20_final_verification/summary.json")
    require(len(v20["checks"]) == 24 and all(r["status"] == "pass" and r["actual"] == r["expected"] for r in v20["checks"]), "v20 content recovery incomplete")
    require(codex["v20_release_bindings_verified"] == 24 and not v20["source_release_modified"], "v20 source release changed")
    return {"real_pairs_accounted": 12, "real_candidate_pairs_completed": 0, "synthetic_completed_pairs": 2,
            "CODEX_exact_checks": 60, "v20_audited_bindings": 24}


def check_history(base: Path) -> dict:
    from multiplex5000.parameter_observations import digest, validate

    run, doc = base / RUN_REL, base / DOC_REL
    project = run / "parameter-history"
    exported = read(project / "parameter_observations.json")
    require(read(doc / "parameter_history.json") == exported, "delivered history differs from project export")
    require(exported["storage"] == "project_only", "history is not project-only")
    pins = read(base / "docs/evaluations/dynamic_portability_20261102/pins.json")
    require(exported["kintsugi_commit"] == pins["KINTSUGI"], "history producer differs from pinned KINTSUGI")
    require(not exported["automatic_application"] and not exported["recommendations_enabled"], "history enables tuning")
    revisions, aids = defaultdict(list), set()
    evidence = set()
    for row in exported["records"]:
        value = validate(row["observation"])
        require(row["sha256"] == digest(value), "observation digest mismatch")
        require(value["observation_id"].startswith("remediation-20261102:"), "new project history mixes unrelated records")
        revisions[value["observation_id"]].append(value["revision"])
        aids.update(value["context"]["acquisition_ids"])
        for item in value["evidence"] + value["context"]["input_receipts"]:
            path = Path(item["path"]).resolve()
            require(path.is_relative_to((project / "evidence").resolve()), "observation evidence is not preserved in its project archive")
            require(path.parent.name == item["sha256"], "observation archive does not use its content identity")
            evidence.add((path, item["sha256"]))
    require(revisions and all(values == list(range(1, len(values) + 1)) for values in revisions.values()), "history revisions are discontinuous")
    for path, expected in evidence:
        require(sha(path) == expected, f"archived observation evidence changed: {path}")
    baseline = read(base / "docs/evaluations/dynamic_portability_20261102/acquisition_dispositions.json")
    required = {row["acquisition_id"] for row in baseline if row["cohort"] == "io60" and row["modality"] in {"phenocycler", "xenium"}}
    required |= {row["acquisition_id"] for row in baseline if row["cohort"] == "codex_baseline"}
    require(required <= aids, "new full-cohort observation contexts omit remediated acquisitions")
    pairs = {r["observation"]["observation_id"] for r in exported["records"]
             if r["observation"]["operation"] == "same_slide_pair_preflight"}
    require(len(pairs) == 12, "new parameter history omits declared same-slide pair contexts")
    failed = [row["observation"] for row in exported["records"]
              if row["observation"]["observation_id"] == "remediation-20261102:protein-stage-publication-attempt1"]
    require(failed, "history omits the retained protein attempt1 publication failure")
    failure = read(base / FAILED_PROTEIN_ATTEMPT1_REL / "failure.json")
    for value in failed:
        settings = {item["name"]: item["value"] for item in value["settings"]}
        scope = settings["attempt1_publication_evidence_scope"]
        require(value["status"] == "failed" and value["operation"] == "protein_stage_publication"
                and not value["context"]["input_receipts"] and len(value["context"]["acquisition_ids"]) == 25
                and all(score is None for score in value["quality"].values())
                and settings["attempt1_publication_failure"] == failure
                and scope["complete_outer_receipt"] is None and scope["full_execution_input_closure"] is None
                and scope["retained_outputs"] == "unapproved_unpublished"
                and scope["scientific_parameter_failure"] == "not_established"
                and scope["large_scientific_outputs_audited"] is False and scope["recovered_receipt"] is False,
                "publication failure was omitted, promoted or conflated with a scientific parameter failure")
    recording, mcp = read(project / "recording_summary.json"), read(project / "scientific_mcp_check.json")
    require(recording["status"] == mcp["status"] == "pass" and recording["all_retries_duplicate"] and mcp["retry_duplicate"], "history/MCP idempotency failed")
    require(recording["records"] == mcp["records"] == len(exported["records"]) and recording["approved_records"] == 0, "record counts or approval state differ")
    require(mcp["export_sha256"] == sha(project / "parameter_observations.json"), "MCP export differs")
    require(recording["export"]["sha256"] == mcp["export_sha256"]
            and sha(project / "mcp_verified_export.json") == mcp["export_sha256"], "recorded and MCP export bytes differ")
    require(recording["global_before"] == recording["global_after"] == mcp["global_before"] == mcp["global_after"], "global history changed")
    global_db = Path.home() / ".kintsugi/parameter_learning.db"
    require((sha(global_db) if global_db.is_file() else None) == recording["global_after"], "global database differs from checked state")
    database = project / ".claude_workflow/parameter_learning.db"
    require(sha(database) == recording["database_sha256"], "project database changed after its recording check")
    with sqlite3.connect(f"file:{database}?mode=ro", uri=True) as connection:
        require(connection.execute("SELECT COUNT(*) FROM parameter_records").fetchone()[0] == len(exported["records"]), "SQLite/export counts differ")
        require(connection.execute("SELECT COUNT(*) FROM parameter_records WHERE user_approved != 0").fetchone()[0] == 0, "unapproved observations became approved")
    require(not list((project / "parameter_observation_outbox").glob("*.json")), "retryable recording failures remain")
    protocol = read(run / "mcp-protocol/protocol_check.json")
    require(protocol["status"] == "pass" and protocol["approval_rejected"] and protocol["idempotency"] == protocol["revision_history"] == "pass", "fresh MCP protocol checks failed")
    require(protocol["global_db_before"] == protocol["global_db_after"] and protocol["records"] == 2, "protocol check changed global history")
    require(protocol["server_sha256"] == sha(base / "core/multiplex5000/multiplex5000/parameter_observations.py"),
            "MCP protocol check used different server bytes")
    return {"observations": len(revisions), "records": len(exported["records"]), "archived_evidence_files": len(evidence), "approved": 0}


def check_manifest(base: Path) -> dict:
    path = base / DOC_REL / "artifact_manifest.json"
    manifest = read(path)
    require(manifest["schema"] == "portability-remediation-artifacts/1", "unknown final manifest schema")
    files, sources = catalog(base, path)
    fresh = [row for row in manifest["artifacts"] if row["kind"] in {"current_source", "evidence_file"}]
    require(len(fresh) == len(files) and {resolve(base, row["path"]) for row in fresh} == files,
            "manifest omits or duplicates current source/evidence files")
    require({resolve(base, row["path"]) for row in fresh if row["kind"] == "current_source"} == sources,
            "manifest current-source classifications differ from delivered source inventory")
    bound = [row for row in manifest["artifacts"] if row["kind"] not in {"current_source", "evidence_file"}]
    identities = [(row["kind"], row["audit"]["path"], row.get("audit_row", row.get("acquisition_id", row.get("input_name", row["path"]))))
                  for row in bound]
    require(len(set(identities)) == len(bound), "manifest duplicates prior audit bindings")
    run = base / RUN_REL
    protein_folder = run / "protein/verification/release_audit"
    protein_bindings = receipt_rows(protein_folder / "receipt_checks.csv", Path(read(protein_folder / "summary.json")["root"]),
                                   protein_folder.parent / "remediation_checks.json")
    require(dict(Counter(row["kind"] for row in bound)) == {"audited_scientific_file": 204 + len(protein_bindings) - 15,
            "declared_optional_absence": 15, "audited_remote_image": 62, "audited_recorded_input": 10, "audited_content_recovery": 24},
            "manifest scientific binding coverage is incomplete")
    audit_cache = {}
    checked = Counter(verify_artifact(row, base, audit_cache=audit_cache) for row in manifest["artifacts"])
    require(dict(Counter(row["kind"] for row in manifest["artifacts"])) == manifest["counts"], "manifest type counts differ")
    return dict(checked)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--build-manifest", action="store_true")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    base = args.base.resolve()
    checks = []
    for name, function in [("accounting", check_accounting), ("dedicated_pln", check_pln),
                           ("protein_full25", check_protein), ("maldi_sources_and_displays", check_maldi),
                           ("geometry_and_codex", check_geometry_codex), ("parameter_history", check_history)]:
        try:
            checks.append({"check": name, "status": "pass", "details": function(base)})
        except (OSError, ValueError, KeyError, TypeError, sqlite3.Error) as error:
            checks.append({"check": name, "status": "fail", "error": f"{type(error).__name__}: {error}"})
    baseline = subprocess.run([sys.executable, str(base / "scripts/evaluation/verify_portability_delivery.py")],
                              cwd=base, capture_output=True, text=True, check=False)
    try:
        baseline_result = json.loads(baseline.stdout) if baseline.returncode == 0 else None
        require(baseline_result is not None and baseline_result["status"] == "pass", "frozen baseline verification failed")
        checks.append({"check": "frozen_baseline", "status": "pass", "details": baseline_result})
    except (ValueError, KeyError) as error:
        checks.append({"check": "frozen_baseline", "status": "fail", "error": str(error), "stderr": baseline.stderr[-4000:]})
    if all(row["status"] == "pass" for row in checks):
        try:
            if args.build_manifest:
                build(base)
            checks.append({"check": "typed_artifact_manifest", "status": "pass", "details": check_manifest(base)})
        except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
            checks.append({"check": "typed_artifact_manifest", "status": "fail", "error": f"{type(error).__name__}: {error}"})
    result = {"schema": "portability-remediation-delivery-verification/1",
              "status": "pass" if all(row["status"] == "pass" for row in checks) else "fail",
              "checks": checks, "verifier_source": ref(Path(__file__)),
              "scientific_approval": False, "large_scientific_outputs_rehashed": False,
              "source_policy": "Current files require exact current bytes; historical Git fallback is confined to the frozen-baseline verifier."}
    out = args.out or base / RUN_REL / "delivery_verification.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
