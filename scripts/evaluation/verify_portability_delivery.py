#!/usr/bin/env python3
"""Read-only consistency check of the completed dated review bundle."""
import gzip
import hashlib
import json
import sqlite3
import subprocess
from collections import defaultdict
from pathlib import Path

from multiplex5000.parameter_observations import digest, validate

BASE = Path(__file__).resolve().parents[2]
DOC = BASE / "docs/evaluations/dynamic_portability_20261102"
RUN = BASE / "outputs/dynamic_portability_20261102"
DELIVERY_REVISION = "610a6a1140427ebfd67b8bac86bbaa69b46ac739"


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def artifact_binding(artifact):
    """Historical code/docs can evolve; scientific evidence must retain its bytes.

    For an evolved tracked source, verify the exact delivered Git blob instead.
    This does not bless current code or permit replacement of history/data files.
    """
    path = Path(artifact["path"])
    if path.is_file() and sha(path) == artifact["sha256"]:
        return "working_file"
    relative = path.relative_to(BASE)
    if relative.parts[0] == "outputs" or path.is_relative_to(DOC):
        raise AssertionError(f"frozen scientific evidence changed or missing: {path}")
    if path.suffix not in {".py", ".md", ".toml"} and relative.as_posix() != ".gitignore":
        raise AssertionError(f"no source-version fallback for {path}")
    blob = subprocess.run(["git", "show", f"{DELIVERY_REVISION}:{relative.as_posix()}"],
                          cwd=BASE, check=True, capture_output=True).stdout
    assert hashlib.sha256(blob).hexdigest() == artifact["sha256"], path
    return "delivery_git_blob"


def main():
    acquisitions = read(DOC / "acquisition_dispositions.json")
    assert len(acquisitions) == 201
    assert len({r["acquisition_id"] for r in acquisitions}) == 201
    assert all(r["disposition"] in {"completed", "excluded", "unavailable", "blocked"} and r["reason"] and r["evidence"] for r in acquisitions)
    assert {s: sum(r["disposition"] == s for r in acquisitions) for s in ("completed", "excluded", "unavailable", "blocked")} == {"completed": 90, "excluded": 20, "unavailable": 5, "blocked": 86}
    assert len(read(DOC / "pair_dispositions.json")["explicit_pairs"]) == 12
    axes = {"code_availability", "valid_inputs", "scientific_specification", "executable_routing", "computational_validation", "biological_review"}
    assert all(set(r["axes"]) == axes for r in read(DOC / "readiness.json")["routes"])
    history = read(DOC / "parameter_history.json")
    assert not history["automatic_application"] and not history["recommendations_enabled"]
    recorded_acquisitions = {aid for r in history["records"] for aid in r["observation"]["context"]["acquisition_ids"]}
    assert {r["acquisition_id"] for r in acquisitions if r["disposition"] != "excluded"} <= recorded_acquisitions
    revisions, references = defaultdict(list), {}
    for record in history["records"]:
        observation = validate(record["observation"])
        assert record["sha256"] == digest(observation)
        revisions[observation["observation_id"]].append(observation["revision"])
        for ref in observation["evidence"] + observation["context"]["input_receipts"]:
            key = (ref["path"], ref["sha256"])
            references[key] = ref
    assert all(v == list(range(1, len(v) + 1)) for v in revisions.values())
    for path, expected in references:
        assert sha(Path(path)) == expected, f"observation evidence changed: {path}"
    project = RUN / "parameter-history"
    with sqlite3.connect(f"file:{project / '.claude_workflow/parameter_learning.db'}?mode=ro", uri=True) as con:
        assert con.execute("SELECT COUNT(*) FROM parameter_records").fetchone()[0] == len(history["records"])
        assert con.execute("SELECT COUNT(*) FROM parameter_records WHERE user_approved != 0").fetchone()[0] == 0
    assert not list((project / "parameter_observation_outbox").glob("*.json"))
    assert read(project / "scientific_mcp_check.json")["status"] == "pass"
    assert read(RUN / "remote/rna/state_reconciliation.json")["failures"] == 0
    for route, count in (("remote/rna", 65), ("remote/spleen", 12), ("pln_recheck", 2)):
        summary = read(RUN / route / "summary.json")
        assert summary["acquisitions"] == count and not summary["identity_failures"], route
    primary = read(RUN / "protein_primary/summary.json")
    assert primary["acquisitions"] == 25 and primary["completed_typing_acquisitions"] == 0
    assert primary["status"] == "blocked" and not primary["hash_failures"]
    assert not primary["identity_failures"]
    independent = read(RUN / "protein_primary_expected_universe.json")
    assert independent["acquisitions"] == 25 and not independent["identity_failures"] and not independent["hash_failures"]
    assert "0112" in primary["zero_measurable_acquisitions"]
    assert read(RUN / "protein_primary/failure.json")["status"] == "failed"
    nucleus = read(RUN / "protein_nucleus.json")
    assert nucleus["acquisitions"] == 25 and not nucleus["failures"] and not nucleus["unavailable"]
    coordinates = read(RUN / "protein_primary_coordinates.json")
    assert coordinates["acquisitions"] == 25 and coordinates["primary_cells"] == independent["primary_cells"]
    inspections = read(DOC / "inspection_ledger.json")
    assert len({(r["route"], r["section_key"]) for r in inspections}) == len(inspections)
    expected_panels = {(route, path.stem) for route in
                       ("remote/rna", "remote/spleen", "pln", "maldi", "remote/spleen_input_geometry", "protein_primary")
                       for path in (RUN / route / "overlays").glob("*.png")}
    assert len(expected_panels) == 144
    assert {(r["route"], r["section_key"]) for r in inspections} == expected_panels
    for entry in inspections:
        assert sha(BASE / entry["spatial_diagnostic"]) == entry["sha256"]
    source_review = read(DOC / "maldi_source_review.json")
    assert {r["acquisition_id"] for r in source_review["sections"]} == {r["acquisition_id"] for r in acquisitions if r["modality"] == "maldi"}
    assert all(r["visual_review"] == "completed" and r["accepted_transform"] is None and r["accepted_neighbour"] is None for r in source_review["sections"])
    for artifact in [source_review["source_figure"], *source_review["viewed_crops"]]:
        assert sha(Path(artifact["path"])) == artifact["sha256"]
    objects = read(DOC / "pln_object_availability.json")
    assert (objects["acquisitions"], objects["objects"], objects["area_objects"]) == (19, 47, 44)
    assert all(r["source_identity_check"] == "pass" and r["registration_acceptance"] is None for r in objects["rows"])
    for report in (read(DOC / "maldi_source_evidence.json"), objects):
        for source in report["sources"]:
            assert sha(Path(source["path"])) == source["sha256"]
    codex = read(DOC / "codex_completion_summary.json")
    assert codex["status"] == "complete" and codex["images"] == 20
    assert codex["unchanged_source_to_port_assertions"] == {"exact_pass": 60, "fail": 0}
    assert codex["unchanged_artifact_assertions"] == {"identity_exact_pass": 20, "levels_exact_fail": 20, "canary_exact_fail": 20}
    assert codex["controlled_table_precision_diagnostic"]["tables_with_any_difference_vs_roundtrip_csv"] == 0
    assert not codex["input_and_baseline_files_unchanged_during_comparisons"]["changed"]
    trace = json.loads(gzip.decompress((DOC / "parameter_traces.json.gz").read_bytes()))
    assert trace["summary"] == read(DOC / "parameter_trace_summary.json")["summary"]
    assert trace["summary"]["catalog_declarations"] == trace["summary"]["traced_declarations"] == 572
    assert trace["summary"]["registry_rows"] == 341
    assert trace["summary"]["remaining_literal_candidates"] == 11218
    requirements = read(DOC / "handoff_requirements.json")
    assert len(requirements["requirements"]) == 23 and not requirements["full_scientific_acceptance"]
    for row in requirements["requirements"]:
        assert all((DOC / p).exists() for p in row["evidence"]), row["id"]
    manifest = read(DOC / "artifact_manifest.json")
    source_replays = [artifact["id"] for artifact in manifest["artifacts"]
                      if artifact_binding(artifact) == "delivery_git_blob"]
    print(json.dumps({"status": "pass", "acquisitions": len(acquisitions),
                      "observations": len(revisions), "history_records": len(history["records"]),
                      "evidence_references_verified": len(references), "inspected_panels": len(inspections),
                      "maldi_source_comparisons": len(source_review["sections"]), "exact_codex_source_port_checks": 60,
                      "parameter_declarations_traced": 572, "handoff_requirements": 23,
                      "manifest_artifacts_verified": len(manifest["artifacts"]),
                      "historical_source_revision": DELIVERY_REVISION,
                      "source_artifacts_verified_from_git": source_replays,
                      "scope": "original delivered source bytes; current code is checked by the separate remediation suite"}, indent=2))


if __name__ == "__main__":
    main()
