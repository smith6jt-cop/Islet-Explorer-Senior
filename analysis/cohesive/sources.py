"""Select existing contracts and preserve the original, path-bound readers."""
from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
import importlib
from pathlib import Path
import shutil
import subprocess
import sys

import pyarrow.parquet as pq

from .common import Context, binding, read_json, require, unique_bindings, verify_binding, write_json


def plain(value):
    if isinstance(value, Mapping):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    return value


def recorded_bindings(value):
    """Extract declared file pins, without treating existence as verification."""
    if isinstance(value, dict):
        if isinstance(value.get("path"), str) and isinstance(value.get("sha256"), str):
            yield {key: value[key] for key in ("path", "sha256", "bytes") if key in value}
        for key, item in value.items():
            if isinstance(item, str) and isinstance(value.get(key + "_sha256"), str):
                yield {"path": item, "sha256": value[key + "_sha256"]}
            if key.endswith("_path") and isinstance(item, str):
                sha = value.get(key[:-5] + "_sha256")
                if isinstance(sha, str):
                    yield {"path": item, "sha256": sha}
            if isinstance(item, (dict, list)):
                yield from recorded_bindings(item)
    elif isinstance(value, list):
        for item in value:
            yield from recorded_bindings(item)


def audit_binding(context: Context, relative: str) -> dict:
    path = (context.inputs.parent / relative).resolve()
    candidates = [record for record in recorded_bindings(context.handoff["audit_reports"])
                  if Path(record["path"]).resolve() == path]
    require(bool(candidates), f"Artifact is not bound by the analytical handoff: {relative}")
    return unique_bindings(candidates)[0]


def original_reader(context: Context):
    api = Path(context.handoff["effective_annotation_reader"]["api"]).resolve()
    sys.path.insert(0, str(api.parent))
    module = importlib.import_module(api.stem)
    require(Path(module.__file__).resolve() == api, "Executing a foreign protein reader")
    return module.load_source_identity_amendment(
        context.handoff["effective_annotation_reader"]["contract"])


def input_bindings(context: Context, operation: str) -> list[dict]:
    require(operation == "reconcile", "Unknown source operation")
    records = list(context.handoff["contracts"].values())
    records.extend(audit_binding(context, relative) for relative in (
        "identity_review/selected_current_rna_inventory.json",
        "protein_composition/handoff.json",
        "protein_composition/qualified_broad_composition/completion.json",
    ))
    export = context.inputs.parent / "RNA_analysis_annotations/export_manifest.json"
    records.append(binding(export))
    export_manifest = read_json(export)
    records.extend(export_manifest[key] for key in ("reader", "input_inventory", "cohort_policy_binding"))
    records.extend(export_manifest["existing_matrix_binding_manifests"])
    amendment = read_json(context.handoff["effective_annotation_reader"]["contract"])
    records.extend(record for key, record in amendment.items()
                   if key.endswith("_api") or key == "builder")
    records.append(binding(context.inputs.parent / "selected_scope_inventory.json"))
    return unique_bindings(records)


def reconcile(context: Context, stage: Path) -> dict:
    handoff = context.handoff
    require(handoff["kind"] == "current_defensible_analytical_handoff", "Wrong input contract kind")
    require(handoff["production_activation"] is False, "This analysis cannot activate typing")
    view = original_reader(context)
    selected = read_json(context.inputs.parent / "selected_scope_inventory.json")
    require(len({row["scope_id"] for row in selected}) == len(selected), "Duplicate selected scopes")
    require(dict(Counter(row["lane"] for row in selected)) == handoff["scope_counts"],
            "Selected scope inventory disagrees with the input contract")

    from phenocycler import cohort
    cohort_pin = handoff["contracts"]["cohort_policy"]
    require(Path(cohort.__file__).resolve() == Path(cohort_pin["path"]).resolve(),
            "The interpreter resolves phenocycler.cohort to a different checkout")
    cohort.ensure_eligible_donors(row["donor_id"] for row in selected)

    scopes = []
    for row in selected:
        metadata = plain(view.scope_metadata(row["scope_id"]).fields)
        for key in ("lane", "donor_id", "path", "sha256", "cell_key_column", "claim_column",
                    "QC_column", "QC_kind", "expected_source_rows", "expected_QC_rows"):
            require(metadata[key] == row[key], f"Selected {key} differs for {row['scope_id']}")
        if "effective_relationship" in row:
            metadata["effective_relationship"] = row["effective_relationship"]
        scopes.append(metadata)
    rna_manifest_path = context.inputs.parent / "RNA_analysis_annotations/export_manifest.json"
    rna_manifest = read_json(rna_manifest_path)
    require(rna_manifest["biological_only"] is True and rna_manifest["named_lineage_only"] is False,
            "RNA export must preserve unresolved biological-QC observations")
    selected_rna = {row["scope_id"]: row for row in scopes if row["lane"].endswith("_RNA")}
    require(set(selected_rna) == {row["scope_id"] for row in rna_manifest["scopes"]},
            "RNA export has a different scope set")
    for row in rna_manifest["scopes"]:
        scope = selected_rna[row["scope_id"]]
        require(row["annotation_sidecar_binding"]["sha256"] == scope["sha256"],
                f"Export uses a different RNA generation: {row['scope_id']}")
        require(row["native_call_nucleus_binding"]["sha256"] == scope["native_source_sha256"],
                f"Export uses different nucleus evidence: {row['scope_id']}")
        require(row["generation"] == handoff["RNA_generations"][
            "pancreas" if scope["lane"] == "pancreas_RNA" else "pLN"],
            f"RNA generation mismatch: {row['scope_id']}")
    protein_handoff = read_json(context.inputs.parent / "protein_composition/handoff.json")
    protein_records = protein_handoff["per_donor_status"]
    require({record["donor_id"] for record in protein_records} ==
            {row["donor_id"] for row in scopes if row["lane"] == "pancreas_protein"},
            "Protein source mapping omits selected donors")
    compact = {
        "handoff": binding(context.inputs),
        "contracts": {key: verify_binding(value) for key, value in handoff["contracts"].items()},
        "scopes": scopes, "protein_records": protein_records,
        "rna_export_manifest": binding(rna_manifest_path), "rna_exports": rna_manifest["scopes"],
        "cohort_policy": verify_binding(cohort_pin),
        "reinstated_donors": sorted(cohort.REINSTATED),
        "production_activation": False, "across_scope_counts_are_additive": False,
        "source_reader_requires_original_location": True,
    }
    write_json(stage / "selected_inputs.json", compact)
    write_json(stage / "protein_handoff.json", protein_handoff)
    write_json(stage / "reader_contract_validation.json", plain(view.summary()))

    source_paths = [
        context.original_repo / "analysis/shared" / name for name in (
            "physical_source_identity_amendment_view.py", "typing_semantic_selection_view.py",
            "protein_semantic_view.py", "typing_release_view.py",
            "tests/test_physical_source_identity_amendment_view.py",
            "tests/test_protein_semantic_view.py",
            "tests/test_typing_semantic_selection_view.py",
            "tests/test_typing_release_view.py",
        )
    ]
    source_paths.extend(context.original_repo /
                        "results/v3/islet_immune_reassessment_20261006/scripts" / name for name in (
                            "evaluate_spatial_H1_H2.py", "evaluate_spatial_H3_H4.py",
                            "evaluate_spatial_RNA_H5_H6.py", "combine_spatial_hypotheses_v2.py",
                        ))
    source_paths.append(context.original_repo / "CLAUDE.md")
    snapshots = []
    for path in source_paths:
        record = binding(path)
        relative = path.relative_to(context.original_repo)
        target = stage / "source_snapshot" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        require(binding(target)["sha256"] == record["sha256"], f"Source snapshot differs: {path}")
        snapshots.append(dict(record, snapshot=str(target.relative_to(stage))))
    source_status = subprocess.check_output(
        ["git", "-C", str(context.original_repo), "status", "--short", "--",
         *[str(path.relative_to(context.original_repo)) for path in source_paths]], text=True)
    write_json(stage / "source_recovery.json", {
        "original_repo": str(context.original_repo),
        "original_HEAD": subprocess.check_output(
            ["git", "-C", str(context.original_repo), "rev-parse", "HEAD"], text=True).strip(),
        "controller_HEAD": subprocess.check_output(
            ["git", "-C", str(Path(__file__).resolve().parents[2]), "rev-parse", "HEAD"], text=True).strip(),
        "selected_source_status": source_status, "files": snapshots,
        "snapshots_are_provenance_not_alternate_executable_readers": True,
        "executing_interpreter": sys.executable, "cohort_module": str(Path(cohort.__file__).resolve()),
    })
    rows = []
    for scope in scopes:
        path = Path(scope["path"])
        pf = pq.ParquetFile(path)
        require(pf.metadata.num_rows == scope["expected_source_rows"],
                f"Scope source-row metadata mismatch: {scope['scope_id']}")
        rows.append({"scope_id": scope["scope_id"], "lane": scope["lane"],
                     "donor_id": scope["donor_id"], "source_rows": pf.metadata.num_rows,
                     "source_bytes": path.stat().st_size})
    import pandas as pd
    pd.DataFrame(rows).to_csv(stage / "scope_inventory.csv", index=False)
    return {"selected_scopes": len(scopes), "RNA_scopes": len(selected_rna),
            "pancreas_protein_scopes": len(protein_records),
            "source_workspace_preserved": True, "source_snapshots": len(snapshots),
            "qualification": "Compact contract and source metadata only; payload qualification is separate"}
