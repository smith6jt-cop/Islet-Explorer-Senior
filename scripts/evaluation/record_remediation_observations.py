#!/usr/bin/env python3
"""Append remediation evidence to a separate project-only KINTSUGI history.

Reads completed receipt-bound artifacts and separately identified failed-attempt
evidence. Does not run, fit, approve, or apply scientific settings. Repeating the
command is idempotent; changed evidence adds a revision. The frozen original
evaluation history is never opened for writing.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sqlite3
import tempfile
from pathlib import Path

import yaml
from multiplex5000.parameter_observations import (
    ObservationHistory,
    canonical,
    file_hash,
    record_or_queue,
)
from record_portability_observations import (
    BASE,
    available,
    context,
    observation,
    read,
    ref,
    setting,
)

RUN = BASE / "outputs/blocker_remediation_20261102"
DOC = BASE / "docs/evaluations/blocker_remediation_20261102"
BASELINE = BASE / "docs/evaluations/dynamic_portability_20261102"


def record(name, status, operation, ctx, subject, settings, outcome, evidence):
    result = observation(name, status, operation, ctx, subject, settings, outcome, evidence)
    result["observation_id"] = "remediation-20261102:" + name
    return available(result)


def bind_sources(payload):
    """Do not record current files as the sources of earlier, different bytes."""
    for item in payload["input_bindings"] + payload["producer_sources"]:
        if file_hash(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"resolved parameter source changed: {item['path']}")


def freeze_references(value, project):
    """Keep referenced evidence bytes available after a later scientific revision."""
    for item in value["evidence"] + value["context"]["input_receipts"]:
        source = Path(item["path"])
        data = source.read_bytes()
        if hashlib.sha256(data).hexdigest() != item["sha256"]:
            raise ValueError(f"evidence changed before recording: {source}")
        target = (project / "evidence" / item["sha256"] / source.name).resolve()
        if not target.is_relative_to(project.resolve()):
            raise ValueError("evidence archive must remain inside project")
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if file_hash(target) != item["sha256"]:
                raise ValueError("immutable evidence archive changed")
        else:
            fd, temporary = tempfile.mkstemp(prefix=".evidence-", dir=target.parent)
            try:
                with os.fdopen(fd, "wb") as handle:
                    handle.write(data)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temporary, target)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
        # Original identity stays in id; path now locates immutable versioned bytes.
        item["path"] = str(target)
    return value


def bound_stage_output(stage, receipt, relative):
    """Read an exact Stage06 output identity; basename matching is insufficient."""
    if receipt.get("stage") != "06_typing" or receipt.get("status") != "complete":
        raise ValueError("component evidence requires a complete Stage 06 receipt")
    matches = [row for row in receipt["outputs"] if row.get("path") == relative]
    path = stage / relative
    if (len(matches) != 1 or not path.is_file() or not path.resolve().is_relative_to(stage.resolve())
            or matches[0].get("sha256") != file_hash(path)):
        raise ValueError(f"component output is missing, stale or not uniquely receipt-bound: {relative}")
    return path


def component_json(path):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate component JSON key: {key}")
            result[key] = value
        return result
    return json.loads(path.read_text(), object_pairs_hook=unique)


def component_csv(path, columns, keys):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if (reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames)
                or not columns <= set(reader.fieldnames)):
            raise ValueError(f"malformed component CSV header: {path}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"malformed component CSV row: {path}")
    identities = [tuple(row[key] for key in keys) for row in rows]
    if len(identities) != len(set(identities)):
        raise ValueError(f"duplicate component CSV identity: {path}")
    return rows


def mixed_component_settings(stage, receipt, donors):
    """Record emitted primary/immune fits and criteria, without fitting or approving."""
    roster = receipt.get("settings", {}).get("donors")
    if (not isinstance(roster, list) or not roster or len(set(roster)) != len(roster)
            or len(donors) != len(set(donors)) or set(roster) != set(donors)):
        raise ValueError("mixed component donor roster differs from the exact stage scope")
    relatives = ["pancreas/subpanel_audit.json", "pancreas/immune/coverage.json",
                 "pancreas/immune/term_comparability.csv", "pancreas/immune/unavailable_rules.json",
                 "pancreas/immune/per_donor_counts.csv"]
    paths = {name: bound_stage_output(stage, receipt, name) for name in relatives}
    primary = component_json(paths[relatives[0]])
    coverage = component_json(paths[relatives[1]])
    unavailable = component_json(paths[relatives[3]])
    if (not isinstance(primary, list) or any(not isinstance(row, dict) or not isinstance(row.get("donor_id"), str)
                                           for row in primary)
            or len(primary) != len(donors) or {row["donor_id"] for row in primary} != set(donors)):
        raise ValueError("missing or duplicate primary component donor audit")
    primary = {row["donor_id"]: row for row in primary}
    required = {"donors_built", "donors_requested", "per_donor", "require_min", "majority_term",
                "not_comparable_criterion", "terms_not_cross_donor_comparable", "supporting_only_notes"}
    if not isinstance(coverage, dict) or not required <= set(coverage):
        raise ValueError("malformed immune component coverage")
    for key in ("donors_built", "donors_requested"):
        if not isinstance(coverage[key], list) or len(coverage[key]) != len(donors) or set(coverage[key]) != set(donors):
            raise ValueError("immune component scope is incomplete or duplicated")
    if (not isinstance(coverage["per_donor"], dict) or set(coverage["per_donor"]) != set(donors)
            or not isinstance(unavailable, dict) or set(unavailable) != set(donors)):
        raise ValueError("immune component donor evidence is incomplete")
    criterion = coverage["not_comparable_criterion"]
    if not isinstance(criterion, dict) or set(criterion) != {"high", "low", "permutations"}:
        raise ValueError("malformed immune comparability criterion")
    for value in (coverage["require_min"], coverage["majority_term"], *criterion.values()):
        if not isinstance(value, (int, float)) or isinstance(value, bool) or not 0 <= value < float("inf"):
            raise ValueError("immune fixed/review setting is not a finite nonnegative number")
    comparison = component_csv(paths[relatives[2]], {"rule", "marker", "n_donors", "cross_donor_comparable"}, ("rule", "marker"))
    for row in comparison:
        if not row["rule"] or not row["marker"] or row["cross_donor_comparable"] not in {"True", "False"}:
            raise ValueError("malformed immune comparability result")
        row["cross_donor_comparable"] = row["cross_donor_comparable"] == "True"
        row["n_donors"] = int(row["n_donors"])
        if not 0 < row["n_donors"] <= len(donors):
            raise ValueError("invalid immune comparability donor count")
        for key in row:
            if key.startswith(("eta2_", "p_perm_")) or key in {"min_positive_fraction", "max_positive_fraction"}:
                row[key] = float(row[key]) if row[key] != "" else None
    counts = component_csv(paths[relatives[4]], {"donor_id", "level", "label", "count", "status", "reason"},
                           ("donor_id", "level", "label"))
    for row in counts:
        if row["donor_id"] not in donors or not row["level"] or not row["label"] or not row["status"]:
            raise ValueError("unknown or malformed immune count identity")
        for key in ("count", "count_without_majority_term"):
            if key in row:
                row[key] = int(row[key]) if row[key] != "" else None
                if row[key] is not None and row[key] < 0:
                    raise ValueError("negative immune count")
    result = {}
    for donor in donors:
        item = primary[donor]
        if not {"islet_levels", "islet_levels_unguarded", "per_stratum", "endocrine_quantile", "islet_chain"} <= set(item):
            raise ValueError("primary component lacks emitted settings or fits")
        if not all(isinstance(item[key], dict) for key in ("islet_levels", "islet_levels_unguarded", "per_stratum")):
            raise ValueError("malformed primary component fitted values")
        if (not isinstance(item["endocrine_quantile"], (int, float)) or isinstance(item["endocrine_quantile"], bool)
                or not 0 <= item["endocrine_quantile"] <= 100 or not isinstance(item["islet_chain"], bool)):
            raise ValueError("malformed primary fixed/selected settings")
        immune = coverage["per_donor"][donor]
        if not isinstance(immune, dict) or immune.get("state") not in {"built", "empty", "unavailable"}:
            raise ValueError("malformed immune per-donor component")
        if immune["state"] == "built" and not {"term_positive_fraction", "rule_markers_absent", "panel_columns_entirely_null",
                "immune_class", "immune_subtype", "macrophage_profile", "unavailable_rules"} <= set(immune):
            raise ValueError("built immune component lacks fitted/availability diagnostics")
        if "unavailable_rules" in immune and immune["unavailable_rules"] != unavailable[donor]:
            raise ValueError("immune unavailable-rule artifacts disagree")
        primary_source, immune_source = paths[relatives[0]], paths[relatives[1]]
        result[donor] = [
            setting("primary_endocrine_quantile", item["endocrine_quantile"], primary_source, "fixed", "percentile"),
            setting("primary_islet_chain_selected", item["islet_chain"], primary_source, "selected_method"),
            setting("primary_islet_levels", item["islet_levels"], primary_source, "fitted"),
            setting("primary_islet_levels_unguarded_diagnostics", item["islet_levels_unguarded"], primary_source, "fitted"),
            setting("primary_per_stratum_diagnostics", item["per_stratum"], primary_source, "fitted"),
            setting("primary_section_diagnostics", {key: value for key, value in item.items()
                    if key not in {"donor_id", "endocrine_quantile", "islet_chain", "islet_levels", "islet_levels_unguarded", "per_stratum"}},
                    primary_source, "fitted"),
            setting("immune_require_min", coverage["require_min"], immune_source, "fixed", "scaled intensity fraction"),
            setting("immune_majority_term", coverage["majority_term"], immune_source, "fixed", "positive fraction"),
            setting("immune_comparability_permutations", criterion["permutations"], immune_source, "fixed", "permutations"),
            setting("immune_not_comparable_criterion", criterion, immune_source, "review_criterion"),
            setting("immune_section_fits_and_availability", immune, immune_source, "fitted"),
            setting("immune_cohort_comparability_fits", coverage["terms_not_cross_donor_comparable"], immune_source, "fitted"),
            setting("immune_term_comparability_results", comparison, paths[relatives[2]], "fitted"),
            setting("immune_supporting_only_notes", coverage["supporting_only_notes"], immune_source, "review_criterion"),
            setting("immune_unavailable_rules", unavailable[donor], paths[relatives[3]], "fitted"),
            setting("immune_section_counts", [row for row in counts if row["donor_id"] == donor], paths[relatives[4]], "fitted", "cells"),
        ]
    return result, list(paths.values())


def typing_records(path, evidence, inventory, *, mixed_components=False):
    payload = read(path)
    bind_sources(payload)
    receipt = path.parent / "receipt.json"
    receipt_value = read(receipt)
    if (payload.get("schema") != "multiplex5000-resolved-parameters/1"
            or payload.get("stage") != "06_typing" or receipt_value.get("stage") != "06_typing"
            or not isinstance(payload.get("release"), str) or not payload["release"]
            or payload["release"] != receipt_value.get("release")
            or payload.get("stage_receipt") != "receipt.json"):
        raise ValueError("resolved parameters disagree with Stage 06 receipt identity")
    matches = [o for o in receipt_value["outputs"] if o["path"] == path.name]
    if len(matches) != 1 or matches[0]["sha256"] != file_hash(path):
        raise ValueError("resolved parameters are not uniquely bound by their stage receipt")
    if receipt_value.get("status") != "complete":
        raise ValueError("typing stage has no complete receipt")
    scope = payload["context"]
    producer = "receipt-sha256:" + file_hash(receipt)
    parameters = [setting(p["name"], p["value"], str(path) + "#" + p["name"], p["role"], p.get("units"))
                  for p in payload["parameters"]]
    components, component_evidence = ({}, [])
    if mixed_components:
        components, component_evidence = mixed_component_settings(
            path.parent, receipt_value, [row["donor_id"] for row in payload["acquisitions"]])
    for row in payload["acquisitions"]:
        if "acquisition_id" in row:
            inv = inventory[row["acquisition_id"]]
        else:
            # The primary producer is keyed by the exact inventory donor field;
            # this join identifies that input, never a cross-assay relationship.
            candidates = [v for v in inventory.values() if v["modality"] == "phenocycler"
                          and v["tissue"] == "pancreas" and v["donor_id"] == row["donor_id"]]
            if len(candidates) != 1:
                raise ValueError("primary producer donor key lacks one exact input acquisition")
            inv = candidates[0]
        aid = inv["acquisition_id"]
        units = scope.get("units") or "raw QuPath whole-cell means after cohort UniFORM shift; review candidate"
        ctx = context("protein", scope["assay"], scope["tissue"], inv["panel"],
                      "independent_acquisition", [aid], producer, units, receipt)
        ctx["segmentation_ids"] = [inv["segmentation_id"]] if inv.get("segmentation_id") else []
        ctx["registry_hashes"] = {p["name"]: p["sha256"] for p in payload["input_bindings"]
                                   if "rule" in p["name"] or "registry" in p["name"]}
        markers = [m["marker"] for m in row.get("markers", [])] or row.get("markers_used", [])
        yield record("typing:" + aid, "evaluated", "protein_typing", ctx,
                     {"kind": "markers", "identifiers": markers},
                     parameters + [setting("section_fits_and_availability", row, path, "fitted")]
                     + components.get(row.get("donor_id"), []),
                     "Completed computational candidate with explicit unavailable evidence. "
                     "No biological approval or release activation; original calls and review flags retained.",
                     [receipt, path, *evidence, *component_evidence])


def failed_protein_publication(inventory):
    """Record the real publication failure without promoting unreceipted scientific files."""
    from manifest_remediation_artifacts import FAILED_PROTEIN_ATTEMPT1_REL, require

    stage = BASE / FAILED_PROTEIN_ATTEMPT1_REL
    assessment_path = RUN / "final_guard_review/protein_failed_attempt1_no_recovery_v1.json"
    assessment = component_json(assessment_path)
    paths = {"failure": stage / "failure.json", "resolved_parameters": stage / "resolved_parameters.json",
             "inner_run": stage / "pancreas/run.json", "failed_launch_snapshot": RUN / "protein/producer_snapshot.json",
             "project": RUN / "protein/protein_availability_remediation.yaml"}
    for key, path in paths.items():
        require(assessment[key] == {"path": str(path), "sha256": file_hash(path)},
                f"failed attempt {key} differs from its retained assessment binding")
    require(assessment["decision"] == "no_provenance_safe_complete_recovery_import"
            and assessment["framework_findings"]["persisted_outer_input_receipt"] is False
            and assessment["required_inputs_without_recovered_content_binding"] > 0
            and assessment["scientific_approval"] is False and assessment["production_release_activated"] is False,
            "failed attempt assessment cannot establish a complete or approved candidate")
    receipt = stage / "receipt.json"
    require(not receipt.exists() and not receipt.is_symlink(), "failed attempt unexpectedly has an outer receipt")
    failure, payload, inner, snapshot = (component_json(paths[key]) for key in
                                       ("failure", "resolved_parameters", "inner_run", "failed_launch_snapshot"))
    project = yaml.safe_load(paths["project"].read_text())
    expected_error = f"06_typing read file(s) it did not declare: {[str(BASE / 'core/multiplex5000/scripts/typing/refine.py')]}"
    require(failure["status"] == "failed" and failure["stage"] == "06_typing"
            and failure["error_type"] == "TaskError" and failure["error"] == expected_error,
            "unexpected protein attempt1 publication failure identity")
    require(payload["schema"] == "multiplex5000-resolved-parameters/1" and payload["stage"] == failure["stage"]
            and payload["release"] == failure["release"] == project["release"] == stage.parent.name,
            "failed publication evidence has inconsistent stage/release identity")
    donors = inner["donors"]
    require(len(donors) == len(set(donors)) == 25, "failed attempt must preserve the declared full25 cohort")
    for other in (project["settings"]["typing_donors"], [row["donor_id"] for row in payload["acquisitions"]]):
        require(len(other) == len(donors) and set(other) == set(donors), "failed attempt donor scopes disagree")
    scope = payload["context"]
    require(scope["modality"] == "protein" and scope["assay"] == "PhenoCycler"
            and scope["tissue"] == "pancreas_with_embedded_pLN", "failed attempt assay/tissue differs")
    # Compare recorded old source identities only with the recorded old snapshot.
    # Current worktree bytes are deliberately not presented as this attempt's producer.
    recorded = {row["path"]: row["sha256"] for row in snapshot["sources"]}
    require(len(recorded) == len(snapshot["sources"]), "failed producer snapshot has duplicate source identities")
    old_sources = payload["producer_sources"] + [{"path": path, "sha256": digest} for path, digest in inner["source_scripts"].items()]
    for item in old_sources:
        require(recorded.get(str(Path(item["path"]).relative_to(BASE))) == item["sha256"],
                "failed-attempt recorded source differs from its original snapshot")
    acquisitions = []
    for donor in donors:
        candidates = [row for row in inventory.values() if row.get("cohort") == "io60"
                      and row.get("modality") == "phenocycler" and row.get("tissue") == "pancreas"
                      and row.get("donor_id") == donor]
        require(len(candidates) == 1, "failed attempt donor lacks one exact inventory acquisition")
        acquisitions.append(candidates[0]["acquisition_id"])
    require(len(set(acquisitions)) == 25, "failed attempt acquisition identities are not unique")
    ctx = context("protein", scope["assay"], scope["tissue"], None, "independent_acquisition", acquisitions,
                  "failed-launch-snapshot-sha256:" + file_hash(paths["failed_launch_snapshot"]), scope["units"])
    settings = [setting("attempt1_" + item["name"], item["value"], paths["resolved_parameters"], item["role"], item.get("units"))
                for item in payload["parameters"] if item["role"] in {"fixed", "selected_method", "review_criterion"}]
    settings += [setting("attempt1_declared_project_settings", project["settings"], paths["project"], "selected_method"),
                 setting("attempt1_declared_donors", donors, paths["inner_run"], "selected_method"),
                 setting("attempt1_publication_failure", failure, paths["failure"], "fitted"),
                 setting("attempt1_publication_evidence_scope", {
                     "attempt": 1, "failure_scope": "stage_publication_provenance_audit",
                     "scientific_parameter_failure": "not_established", "complete_outer_receipt": None,
                     "full_execution_input_closure": None, "retained_outputs": "unapproved_unpublished",
                     "large_scientific_outputs_audited": False, "recovered_receipt": False,
                     "required_inputs_without_recovered_content_binding": assessment["required_inputs_without_recovered_content_binding"],
                 }, assessment_path, "fitted")]
    return record("protein-stage-publication-attempt1", "failed", "protein_stage_publication", ctx,
                  {"kind": "acquisition", "identifiers": acquisitions}, settings,
                  "Attempt1 failed the Stage06 publication provenance audit: " + failure["error"]
                  + ". This is not evidence of failed scientific settings. No complete outer receipt or full input closure "
                    "is available; retained scientific files remain unpublished and unapproved, without a full audit or recovered receipt.",
                  [*paths.values(), assessment_path])


def build():
    from manifest_remediation_artifacts import verified_protein_snapshot
    from verify_remediation_delivery import check_protein

    inventory = {r["acquisition_id"]: r for r in read(BASELINE / "acquisition_dispositions.json")}
    records = [failed_protein_publication(inventory)]
    pln_summary = RUN / "config/completion_summary.json"
    pln = read(pln_summary)
    if pln["status"] != "complete":
        raise ValueError("full pLN cohort is incomplete")
    records.extend(typing_records(Path(pln["resolved_parameters"]), [pln_summary], inventory))
    protein = RUN / "protein/releases/protein-availability-candidate/06_typing/resolved_parameters.json"
    # This is intentionally required: a still-running cohort is not recorded as complete.
    check_protein(BASE)
    protein_report = RUN / "protein/verification/remediation_checks.json"
    producer_snapshot = verified_protein_snapshot(protein_report)
    records.extend(typing_records(protein, [
        producer_snapshot,
        RUN / "protein/producer_stability.json",
        Path(read(protein_report)["input_launch_bindings"]["path"]),
        protein_report,
        RUN / "protein/verification/release_audit/summary.json",
        RUN / "protein/verification/inspection_ledger.json",
    ], inventory, mixed_components=True))

    codex_path = RUN / "codex/completion_summary.json"
    codex = read(codex_path)
    ctx = context("protein", "CODEX", "pancreas", "codex-uniform32", "independent_acquisition",
                  [r["acquisition_id"] for r in inventory.values() if r["cohort"] == "codex_baseline"],
                  "contract-sha256:" + file_hash(BASE / "scripts/evaluation/codex_replay_contract.json"),
                  "whole-cell and membrane raw intensity")
    records.append(record("codex-exact-replay", "evaluated", "protein_parity", ctx,
                          {"kind": "feature_set", "identifiers": []},
                          [setting("replay_contract", read(BASE / "scripts/evaluation/codex_replay_contract.json"),
                                   "scripts/evaluation/codex_replay_contract.json", "selected_method"),
                           setting("exact_results", codex, codex_path, "fitted")],
                          "60 exact assertions pass across20 images;24 v20 bindings recovered exactly. "
                          "No numerical tolerance relaxation or v20 scientific approval.",
                          [codex_path, RUN / "codex/artifact_manifest.json"]))

    handoff = DOC / "codex_received_handoff.json"
    ctx = context("protein", "CODEX", "pancreas", None, "evaluation", [],
                  "package-sha256:" + read(handoff)["sha256sums_sha256"], None)
    records.append(record("received-codex-reference-context-fix", "evaluated", "source_provenance_review", ctx,
                          {"kind": "feature_set", "identifiers": []},
                          [setting("received_package_disposition", read(handoff), handoff, "review_criterion")],
                          "Package checksums and supplied evidence reviewed; recipient port not applied. "
                          "Frozen comparison and running cohort sources unchanged; no scientific acceptance.",
                          [handoff, DOC / "codex_received_handoff.md"]))
    order = DOC / "relationship_context.json"
    ctx = context("multimodal", "Xenium/PhenoCycler/MALDI", None, None, "serial_section", [],
                  "user-statement-sha256:" + file_hash(order), None)
    records.append(record("user-reported-typical-section-order", "proposed", "serial_section_relationship_context", ctx,
                          {"kind": "geometry", "identifiers": []},
                          [setting("usual_order_context", read(order), order, "review_criterion")],
                          "User-reported usual cutting order, with exceptions unenumerated. "
                          "No exact acquisition mapping, accepted adjacency or registration is added.", [order]))

    sources_path = RUN / "hpg_maldi_xenium/results_relocated/source_validation.json"
    sources = read(sources_path)
    if not sources["complete"] or sources["declared_acquisitions"] != 65:
        raise ValueError("full65 source audit is incomplete")
    for row in sources["acquisitions"]:
        inv = inventory[row["acquisition_id"]]
        ctx = context("rna", "Xenium", inv["tissue"], inv["panel"], "serial_section",
                      [row["acquisition_id"]], "source-sha256:" + sources["producer_sha256"], "um",
                      BASE / "outputs/dynamic_portability_20261102/remote/rna/upstream/00_inventory/receipt.json")
        parameters = row.get("mask_parameters") or {}
        source_settings = [setting("source_geometry_and_availability", row, sources_path, "fitted")]
        if parameters:
            source_settings.extend([
                setting("mask_method", parameters.get("method"), sources_path, "selected_method"),
                setting("mask_fixed_parameters", {k: v for k, v in parameters.items()
                                                  if k not in {"method", "threshold_log1p"}}, sources_path),
                setting("section_otsu_threshold", parameters.get("threshold_log1p"), sources_path, "fitted", "log1p DAPI intensity"),
                setting("grid_pixel_um", row["grid"]["pixel_um"], sources_path, units="um"),
            ])
        records.append(record("morphology-source:" + row["acquisition_id"],
                              "evaluated" if row["status"] == "validated" else "blocked",
                              "registration_source_validation", ctx,
                              {"kind": "geometry", "identifiers": []},
                              source_settings,
                              "Exact-region morphology source validation only; no MALDI neighbour selected, "
                              "transform estimated or registration accepted.",
                              [sources_path, RUN / "final_guard_review/source_binding_review.json"]))

    objects_path = RUN / "maldi_objects/validation.json"
    objects = read(objects_path)
    object_contracts = {r["acquisition_id"]: r for r in read(RUN / "maldi_objects/object_sources.json")["acquisitions"]}
    for row in objects["acquisitions"]:
        inv = inventory[row["acquisition_id"]]
        ctx = context("protein", "PhenoCycler", "pancreatic_lymph_node", inv["panel"], "serial_section",
                      [row["acquisition_id"]], "manifest-sha256:" + file_hash(RUN / "maldi_objects/object_sources.json"), "um")
        contract = object_contracts[row["acquisition_id"]]
        ctx["input_receipts"] = [ref(Path(item["source_receipt"])) for item in contract["sources"]]
        for item, bound in zip(contract["sources"], ctx["input_receipts"], strict=True):
            if bound["sha256"] != item["source_receipt_sha256"]:
                raise ValueError("object observation receipt changed")
        records.append(record("object-source:" + row["acquisition_id"], "evaluated", "maldi_object_source", ctx,
                              {"kind": "geometry", "identifiers": []},
                              [setting("source_geometry_availability", row, objects_path, "fitted"),
                               setting("object_frame_and_calibration_contract", contract, RUN / "maldi_objects/object_sources.json")],
                              "Observed export provenance retains native IDs and unavailable LineStrings. "
                              "No historical producer receipt, MALDI relationship or object association invented.",
                              [objects_path, RUN / "maldi_objects/object_sources.json"]))

    for name, path, operation in [
        ("same-slide-synthetic-checks", RUN / "geom2_revision1/completion_summary.json", "matching_validation"),
        ("maldi-object-synthetic-checks", RUN / "maldi_stage08_review.json", "maldi_object_validation"),
    ]:
        ctx = context("multimodal", None, None, None,
                      "same_slide" if name.startswith("same-slide") else "serial_section", [],
                      "evidence-sha256:" + file_hash(path), "um")
        records.append(record(name, "evaluated", operation, ctx, {"kind": "geometry", "identifiers": []},
                              [setting("synthetic_validation_results", read(path), path, "fitted")],
                              "Synthetic computational checks, not biological or cohort registration evidence.", [path]))
    for name, path, operation, outcome in [
        ("spleen-input-hashes", RUN / "hpg_maldi_xenium/spleen_algorithm_verification/verification.json",
         "input_provenance_validation", "Current full-byte hashes are separate from historical full/sampled manifest bindings; no new biological acceptance."),
        ("maldi-scaling", RUN / "maldi_scaling/manifest.json", "maldi_display_scaling",
         "All32 sections displayed in linear/log1p views; no scientific normalization or new mask."),
        ("morphology-inspection", RUN / "maldi_xenium_review.json", "registration_source_inspection",
         "All62 available exact-region masks inspected; concerns retained and no mask/registration accepted."),
        ("same-slide-cohort-preflight", RUN / "geom1/full_cohort_preflight/ledger.json", "same_slide_batch_preflight",
         "All declared pairs accounted independently. Missing inputs and scientific capabilities remain explicit; no gates accepted."),
    ]:
        ctx = context("multimodal", None, None, None, "evaluation", [],
                      "evidence-sha256:" + file_hash(path), None)
        records.append(record(name, "evaluated", operation, ctx, {"kind": "feature_set", "identifiers": []},
                              [setting("observed_results", read(path), path, "fitted")], outcome, [path]))
    batch_path = RUN / "geom1/full_cohort_preflight/ledger.json"
    batch = read(batch_path)
    for pair in batch["pairs"]:
        declared = pair["declared_pair"]
        ids = (declared.get("protein_acquisition_ids", []) + declared.get("rna_acquisition_ids", [])
               if "protein_acquisition_ids" in declared else
               [declared["protein_entry"], declared["xenium_acquisition_id"]])
        if any(aid not in inventory for aid in ids):
            raise ValueError("pair references an acquisition outside the frozen roster")
        ctx = context("protein/rna", "PhenoCycler/Xenium", declared["tissue"],
                      " / ".join(str(inventory[aid]["panel"]) for aid in ids), "same_slide", ids,
                      "source-sha256:" + batch["producer"]["sha256"], "um")
        ctx["segmentation_ids"] = sorted({inventory[aid]["segmentation_id"] for aid in ids
                                          if inventory[aid].get("segmentation_id")})
        status = "blocked" if pair["status"] in {"blocked", "unavailable"} else "failed" if pair["status"] == "failed" else "evaluated"
        records.append(record("pair-preflight:" + pair["pair_id"], status, "same_slide_pair_preflight", ctx,
                              {"kind": "acquisition_pair", "identifiers": ids},
                              [setting("required_evidence", pair["required_evidence"], batch_path, "review_criterion"),
                               setting("pair_disposition", pair, batch_path, "review_criterion")],
                              f"Exact declared pair preflight disposition: {pair['status']}. "
                              "This recording adds no cell links, uncertainty measurements or identity agreement.", [batch_path]))
    readiness = DOC / "readiness.json"
    for item in read(readiness)["implementation_items"]:
        if "pending" in item["status"]:
            raise ValueError("final recording requires dispositions for pending implementations")
        ctx = context("multimodal", None, None, None, "evaluation", [], "remediation-evaluation/1", None)
        status = "blocked" if ("required" in item["status"] or "blocked" in item["status"]) else "evaluated"
        records.append(record("dependency:" + item["id"], status, "portability_dependency", ctx,
                              {"kind": "feature_set", "identifiers": []},
                              [setting("implementation_disposition", item, readiness, "review_criterion")],
                              item["scope"], [readiness]))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=RUN / "parameter-history")
    args = parser.parse_args()
    project = args.project.resolve()
    if project == (BASE / "outputs/dynamic_portability_20261102/parameter-history").resolve():
        parser.error("the frozen baseline history cannot be modified")
    records = [freeze_references(value, project) for value in build()]
    global_db = Path.home() / ".kintsugi/parameter_learning.db"
    before = file_hash(global_db) if global_db.is_file() else None
    store = ObservationHistory(project, BASE / "external/KINTSUGI")
    for value in records:
        prior = store.history(value["observation_id"])
        if prior:
            previous = prior[-1]["observation"]
            value.update(revision=previous["revision"], supersedes=previous["supersedes"])
            if canonical(value) != canonical(previous):
                value.update(revision=previous["revision"] + 1, supersedes=previous["revision"])
        result = record_or_queue(project, BASE / "external/KINTSUGI", value)
        if result["status"] != "recorded":
            raise RuntimeError(f"retryable recording failure; scientific artifacts preserved: {result}")
    retries = [record_or_queue(project, BASE / "external/KINTSUGI", value) for value in records]
    assert all(result["duplicate"] for result in retries)
    exported = store.export()
    after = file_hash(global_db) if global_db.is_file() else None
    assert before == after
    with sqlite3.connect(store.db) as connection:
        approved = connection.execute("SELECT COUNT(*) FROM parameter_records WHERE user_approved = 1").fetchone()[0]
    assert approved == 0
    shutil.copyfile(exported["path"], DOC / "parameter_history.json")
    summary = {"status": "pass", "observations_submitted": len(records), "records": exported["records"],
               "approved_records": approved, "all_retries_duplicate": True,
               "global_before": before, "global_after": after, "export": exported,
               "automatic_application": False, "recommendations_enabled": False,
               "database_sha256": file_hash(store.db)}
    (project / "recording_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
