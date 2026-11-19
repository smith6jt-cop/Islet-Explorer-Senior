"""Coordinate exact same-slide projects without sharing their gates or outputs.

The batch ledger accounts for a hash-bound roster. Each supported pair still runs
through the ordinary Stage 07/08 contract in a distinct explicitly configured
project. This module never accepts a gate, invents a project or copies acceptance.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pyarrow.parquet as pq

from . import cohort, gates, tasks
from . import project as projects
from .artifacts import sha256_file
from .integration.paired_typing import (
    IDENTITY_COLUMNS,
    validate_pair_cells,
    validate_pair_transform,
)
from .integration.transform import Transform
from .stages import same_slide_integration as integration
from .stages import same_slide_registration as registration

SCHEMA = "multiplex5000-same-slide-batch/1"
SECTION_FIELDS = ("protein_entry", "xenium_acquisition_id", "donor_id", "tissue",
                  "protein_frame", "xenium_frame", "relationship", "xenium_output")


class BatchError(ValueError):
    """The batch cannot establish its exact roster or isolated execution contract."""


def _json(path: Path):
    return json.loads(path.read_text())


def _binding(path: Path) -> dict:
    return {"path": str(path.resolve()), "sha256": sha256_file(path)}


def _path(value: str, base: Path) -> Path:
    if not isinstance(value, str) or not value:
        raise BatchError("a configured path must be a nonempty string")
    path = Path(value)
    return (path if path.is_absolute() else base / path).resolve()


def load_config(path: Path) -> dict:
    path = path.resolve()
    spec = _json(path)
    if set(spec) != {"schema", "roster", "workspace_root", "pairs"} or spec["schema"] != SCHEMA:
        raise BatchError(f"expected {SCHEMA} with roster, workspace_root and pairs")
    if set(spec["roster"]) != {"path", "sha256"}:
        raise BatchError("roster needs an exact path and SHA256")
    roster_path = _path(spec["roster"]["path"], path.parent)
    if sha256_file(roster_path) != spec["roster"]["sha256"]:
        raise BatchError("roster changed from its configured SHA256")
    roster = _json(roster_path)
    rows = roster.get("explicit_pairs")
    if not isinstance(rows, list) or not rows:
        raise BatchError("roster needs nonempty explicit_pairs")
    identifiers = [row.get("pair_id") for row in rows]
    if (any(not isinstance(value, str) or not value for value in identifiers)
            or len(set(identifiers)) != len(identifiers)):
        raise BatchError("roster pair IDs must be present and unique")
    exact_sections = [(row["protein_entry"], row["xenium_acquisition_id"]) for row in rows
                      if row.get("protein_entry") and row.get("xenium_acquisition_id")]
    if len(set(exact_sections)) != len(exact_sections):
        raise BatchError("one physical pair cannot be duplicated under different pair IDs")
    configs = spec["pairs"]
    if not isinstance(configs, list):
        raise BatchError("pairs must be a list of explicit configurations")
    for row in configs:
        if set(row) != {"pair_id", "project", "required_evidence"}:
            raise BatchError("each configuration needs pair_id, project and required_evidence")
        if row["required_evidence"] not in {"independent_calls", "spatial_gene_evidence"}:
            raise BatchError("required_evidence must be independent_calls or spatial_gene_evidence")
        if row["project"] is not None and not isinstance(row["project"], str):
            raise BatchError("project must be an explicit path or null when unavailable")
    chosen = [row["pair_id"] for row in configs]
    if len(set(chosen)) != len(chosen) or set(chosen) != set(identifiers):
        raise BatchError("configurations must cover every roster pair exactly once")
    workspace = _path(spec["workspace_root"], path.parent)
    loaded, errors, roots = {}, {}, []
    for row in configs:
        if row["project"] is None:
            continue
        project_path = _path(row["project"], path.parent)
        try:
            project = projects.load(project_path)
            root = project.release_dir.resolve()
            # A shared or nested root can overwrite another pair's stages.
            for other in roots:
                if root == other or root.is_relative_to(other) or other.is_relative_to(root):
                    raise BatchError("pair projects must have distinct, nonoverlapping release roots")
            roots.append(root)
            loaded[row["pair_id"]] = project
        except (OSError, projects.ProjectError, ValueError) as error:
            if isinstance(error, BatchError):
                raise
            errors[row["pair_id"]] = str(error)
    return {"path": path, "binding": _binding(path), "roster_binding": _binding(roster_path),
            "roster": rows, "unresolved_relationships": roster.get("unresolved", []),
            "configs": {row["pair_id"]: row for row in configs}, "projects": loaded,
            "project_errors": errors, "workspace": workspace}


def _require(row: dict, code: str, detail: str, *, kind: str = "blocked", path: Path | None = None):
    item = {"code": code, "detail": detail, "disposition": kind}
    if path is not None:
        item["path"] = str(path)
    row["requirements"].append(item)


def _task(stage, project):
    declared = stage.inputs(project)
    if len({item.name for item in declared}) != len(declared):
        raise BatchError("stage declares duplicate input names")
    return tasks.Task(stage.STAGE, project.release, project.release_dir / ".batch-preflight",
                      {item.name: item for item in declared}, stage.settings(project), project)


def _acceptance(folder: Path, stage: str) -> dict:
    receipt_path = folder / "receipt.json"
    receipt = _json(receipt_path)
    accepted = gates.require_accepted(folder)
    if (receipt.get("status") != "complete" or receipt.get("stage") != stage
            or accepted.get("stage") != stage
            or accepted.get("receipt_sha256") != sha256_file(receipt_path)):
        raise BatchError(f"{stage} acceptance does not bind its complete current receipt")
    return receipt


def _current_stage(stage, project) -> dict | None:
    """Check inputs, settings and outputs before the task runner's reuse shortcut."""
    folder = project.release_dir / stage.STAGE
    if not folder.exists():
        return None
    receipt = _json(folder / "receipt.json")
    if (receipt.get("stage") != stage.STAGE or receipt.get("status") != "complete"
            or receipt.get("release") != project.release):
        raise BatchError(f"{stage.STAGE}: no complete receipt for this release")
    task = _task(stage, project)
    if receipt.get("settings") != task.settings:
        raise BatchError(f"{stage.STAGE}: current settings differ from the published receipt")
    bound = registration._receipt_hashes(receipt.get("inputs", []), "name")
    if set(bound) != set(task.inputs):
        raise BatchError(f"{stage.STAGE}: declared input roster changed")
    for name, item in task.inputs.items():
        path = Path(item.path)
        if path.is_file():
            if bound[name] != sha256_file(path):
                raise BatchError(f"{stage.STAGE}: input {name} changed")
        elif item.required or bound[name] is not None:
            raise BatchError(f"{stage.STAGE}: input {name} is missing")
    outputs = registration._receipt_hashes(receipt.get("outputs", []), "path")
    if not outputs:
        raise BatchError(f"{stage.STAGE}: no bound outputs")
    for relative, expected in outputs.items():
        path = folder / relative
        if not path.resolve().is_relative_to(folder.resolve()) or not path.is_file() or sha256_file(path) != expected:
            raise BatchError(f"{stage.STAGE}: output {relative} changed or is missing")
    return receipt


def _validate_project(project, declared: dict, workspace: Path):
    root = project.release_dir.resolve()
    if root == workspace or not root.is_relative_to(workspace):
        raise BatchError("project release must be inside the configured isolated workspace")
    for name in ["07_registration", "08_integration"]:
        if not (project.release_dir / name).resolve().is_relative_to(root):
            raise BatchError("stage output points outside its isolated pair release")
    if set(project.automatic_gates) & {"06_typing", "07_registration", "08_integration"}:
        raise BatchError("typing, registration and integration require explicit review; automatic gates refused")
    if project.previous("07_registration") != "06_typing" or project.previous("08_integration") != "07_registration":
        raise BatchError("pair project must retain the Stage 06 -> 07 -> 08 gate chain")
    if project.settings.get("registration_route") is not None or project.settings.get("integration_route") is not None:
        raise BatchError("project is not the supported same-slide section route")
    pair = registration._selected(project)
    for name in SECTION_FIELDS:
        actual = str(getattr(pair, name))
        if actual != str(declared[name]):
            raise BatchError(f"configured project differs from roster field {name}")
    for name in ["protein_export", "protein_image"]:
        if name in declared and str(getattr(pair, name)) != declared[name]:
            raise BatchError(f"configured project differs from roster field {name}")
    return pair


def _typing_evidence(project, pair, row: dict) -> bool:
    """Authenticate the section-keyed producer contract; never adapt donor tables."""
    folder = project.release_dir / "06_typing"
    receipt = _acceptance(folder, "06_typing")
    task = _task(registration, project)
    bound = registration._receipt_hashes(receipt.get("inputs", []), "name")
    if bound.get("same_slide_source_manifest") != sha256_file(task.input("same_slide_source_manifest")):
        raise BatchError("accepted typing source manifest changed or is unbound")
    registration._typed_cells_current(task, receipt)
    common = ["cell_uid", "donor_id", "tissue", "acquisition_id", "x_um", "y_um"]
    validate_pair_cells(pair, registration._cell_table(task, "protein_cells", common + ["protein_eligible"]),
                        registration._cell_table(task, "rna_cells", common + ["typing_eligible"]))
    sample = registration._sample(project, pair)
    outputs = registration._receipt_hashes(receipt.get("outputs", []), "path")
    present = {}
    for assay, name in [("protein", pair.protein_entry), ("rna", sample)]:
        relative = f"{assay}/{name}.parquet"
        path = folder / relative
        if not path.is_file() or outputs.get(relative) != sha256_file(path):
            raise BatchError(f"independent {assay} evidence is missing, changed or not receipt-bound")
        schema = set(pq.read_schema(path).names)
        columns = IDENTITY_COLUMNS & schema
        if (assay == "protein" and not columns) or (assay == "rna" and not columns and not any(col.startswith("gene:") for col in schema)):
            raise BatchError(f"independent {assay} table has no supported call/evidence columns")
        present[assay] = {"identity_columns": sorted(columns), "evidence": _binding(path),
                          "has_nonempty_identity_call": False}
        if columns:
            calls = pq.read_table(path, columns=sorted(columns)).to_pandas()
            present[assay]["has_nonempty_identity_call"] = bool(
                any((calls[col].notna() & calls[col].ne("")).any() for col in calls))
    row["independent_evidence"] = present
    row["evidence"].extend([_binding(folder / "receipt.json"), _binding(folder / "gate/accepted.json")])
    full = all(item["has_nonempty_identity_call"] for item in present.values())
    if row["required_evidence"] == "independent_calls" and not full:
        _require(row, "independent_identity_calls_unavailable",
                 "Gene evidence or entirely unavailable labels do not establish independent assay identity calls.",
                 kind="unavailable")
    if row["required_evidence"] == "independent_calls" and not project.settings.get("same_slide_hierarchy_map"):
        _require(row, "hierarchy_mapping_unavailable", "An explicit tissue/hierarchy comparison map is required.", kind="unavailable")
    return full


def _evaluate_pair(batch: dict, declared: dict, action: str) -> dict:
    identifier = declared["pair_id"]
    config = batch["configs"][identifier]
    row = {"pair_id": identifier, "declared_pair": declared, "required_evidence": config["required_evidence"],
           "status": "blocked", "requirements": [], "evidence": [],
           "registration": "not_evaluated", "integration": "not_evaluated",
           "scientific_approval": False, "identity_comparison": "unavailable"}
    donor = declared.get("donor_id")
    if donor and cohort.excluded_on(donor, "phenocycler"):
        row["status"] = "excluded"
        _require(row, "existing_protein_platform_exclusion", "The central cohort rule excludes this protein donor.", kind="excluded")
        row["evidence"].append(_binding(Path(cohort.__file__)))
        return row
    missing = [name for name in SECTION_FIELDS if not declared.get(name)]
    if missing:
        _require(row, "exact_section_adapter_unavailable",
                 f"The section-keyed route lacks {missing}; donor-level IO60 outputs cannot be substituted.", kind="unavailable")
        row["status"] = "unavailable"
        return row
    if identifier in batch["project_errors"]:
        _require(row, "invalid_project_configuration", batch["project_errors"][identifier])
        return row
    project = batch["projects"].get(identifier)
    if project is None:
        _require(row, "pair_project_unavailable", "Provide an exact isolated per-pair project with receipt-bound independent assay outputs.")
        _require(row, "accepted_registration_unavailable", "No configured per-pair project establishes a reviewed transform or accepted Stage 07 receipt.")
        if config["required_evidence"] == "independent_calls":
            _require(row, "independent_identity_calls_unavailable", "No configured project provides independent identity calls and an explicit hierarchy map.", kind="unavailable")
        return row
    row["project"] = str(project.path)
    row["release_dir"] = str(project.release_dir)
    project_binding = _binding(project.path)
    row["evidence"].append(project_binding)
    try:
        pair = _validate_project(project, declared, batch["workspace"])
    except (OSError, ValueError, KeyError) as error:
        _require(row, "invalid_pair_project", str(error))
        return row
    try:
        _typing_evidence(project, pair, row)
    except (OSError, ValueError, KeyError, gates.GateError, tasks.TaskError) as error:
        _require(row, "independent_typing_not_ready", str(error))
    try:
        transform_path = registration._transform_dir(project, pair) / "transform.json"
        transform = Transform.load(transform_path.parent)
        validate_pair_transform(pair, transform)
        row["evidence"].append(_binding(transform_path))
        row["registration"] = "source_geometry_available_not_accepted"
    except (OSError, ValueError, KeyError) as error:
        _require(row, "registration_source_not_ready", str(error))
    # Biology can remain unavailable while an independently supplied geometry is
    # prepared for review. Failed producer/geometry contracts cannot proceed.
    if any(item["disposition"] == "blocked" for item in row["requirements"]):
        return row
    try:
        registered = _current_stage(registration, project)
        if registered is None and action != "preflight":
            tasks.run_stage(registration, project)
            registered = _current_stage(registration, project)
        if registered is None:
            row["registration"] = "ready_to_run"
            _require(row, "registration_not_run", "Geometry and source contracts are ready for Stage 07 execution and review.")
            return row
        row["evidence"].append(_binding(project.release_dir / "07_registration/receipt.json"))
        row["registration"] = "published_awaiting_review"
        try:
            _acceptance(project.release_dir / "07_registration", "07_registration")
        except (OSError, ValueError, gates.GateError) as error:
            _require(row, "registration_review_required", str(error))
            return row
        row["registration"] = "accepted"
        row["evidence"].append(_binding(project.release_dir / "07_registration/gate/accepted.json"))
        integration._accepted_transform(_task(integration, project))
        if row["requirements"]:
            row["status"] = "unavailable"
            return row
        integrated = _current_stage(integration, project)
        if integrated is None and action == "integration":
            tasks.run_stage(integration, project)
            integrated = _current_stage(integration, project)
        if integrated is None:
            row["integration"] = "ready_to_run"
            _require(row, "integration_not_run", "Accepted geometry and required assay evidence are ready for Stage 08 execution.")
            return row
        row["evidence"].append(_binding(project.release_dir / "08_integration/receipt.json"))
        row["integration"] = "candidate_computation_complete"
        row["summary"] = integrated["summary"]
        if not integrated.get("gate", {}).get("passed"):
            _require(row, "integration_gate_failed", "Candidate output exists, but its computational gate failed.")
            return row
        row["status"] = "completed"
        row["identity_comparison"] = integrated["summary"].get("hierarchy_comparison", {}).get("status", "unavailable")
        # Completion is computational; final candidate review remains separate.
        try:
            _acceptance(project.release_dir / "08_integration", "08_integration")
            row["candidate_review"] = "accepted"
        except (OSError, ValueError, gates.GateError):
            row["candidate_review"] = "not_accepted"
    except (OSError, ValueError, KeyError, gates.GateError, tasks.TaskError) as error:
        _require(row, "stage_execution_or_reuse_failed", f"{type(error).__name__}: {error}")
    finally:
        if sha256_file(project.path) != project_binding["sha256"]:
            row["status"] = "blocked"
            _require(row, "project_changed_during_batch", "The configured project changed during evaluation.")
    return row


def coordinate(config: Path, out: Path, *, action: str = "preflight") -> dict:
    if action not in {"preflight", "registration", "integration"}:
        raise BatchError("action must be preflight, registration or integration")
    batch = load_config(config)
    out = out.resolve()
    if any(out.is_relative_to(project.release_dir.resolve()) for project in batch["projects"].values()):
        raise BatchError("ledger snapshots must be outside per-pair stage release roots")
    # Ledger snapshots are immutable. Stage candidates live in the configured
    # isolated roots and are reused only after validation, never force-replaced.
    out.mkdir(parents=True, exist_ok=False)
    ledger = {"schema": "multiplex5000-same-slide-batch-ledger/1", "action": action,
              "config": batch["binding"], "roster": batch["roster_binding"],
              "producer": _binding(Path(__file__)),
              "stage_producers": [_binding(Path(stage.__file__)) for stage in [registration, integration]],
              "workspace_root": str(batch["workspace"]),
              "declared_pairs": len(batch["roster"]), "pairs": [],
              "unresolved_relationships": batch["unresolved_relationships"],
              "scientific_approval": False, "automatic_application": False}
    for declared in batch["roster"]:
        try:
            row = _evaluate_pair(batch, declared, action)
        except Exception as error:
            row = {"pair_id": declared["pair_id"], "declared_pair": declared, "status": "blocked",
                   "requirements": [{"code": "pair_evaluation_failed", "disposition": "blocked",
                                     "detail": f"{type(error).__name__}: {error}"}]}
        ledger["pairs"].append(row)
        ledger["counts"] = dict(Counter(row["status"] for row in ledger["pairs"]))
        ledger["accounted_pairs"] = len(ledger["pairs"])
        ledger["all_pairs_accounted"] = ledger["accounted_pairs"] == ledger["declared_pairs"]
        ledger["full_cohort_candidate_completion"] = (
            ledger["all_pairs_accounted"] and not ledger["unresolved_relationships"]
            and all(row["status"] in {"completed", "excluded"} for row in ledger["pairs"]))
        (out / "ledger.json").write_text(json.dumps(ledger, indent=2, allow_nan=False) + "\n")
    (out / "artifact_manifest.json").write_text(json.dumps({
        "schema": "same-slide-batch-artifacts/1", "artifacts": [_binding(out / "ledger.json")]}, indent=2) + "\n")
    return ledger


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["preflight", "registration", "integration"])
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="new immutable ledger directory")
    args = parser.parse_args(argv)
    try:
        result = coordinate(args.config, args.out, action=args.action)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"same-slide batch: {error}\n")
    print(json.dumps({key: result[key] for key in ["declared_pairs", "accounted_pairs", "counts", "full_cohort_candidate_completion"]}))
    return 0 if result["full_cohort_candidate_completion"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
