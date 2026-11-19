#!/usr/bin/env python3
"""Hash the review bundle and retain scientific receipts without duplicating raw data."""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
RUN = BASE / "outputs/dynamic_portability_20261102"
DOC = BASE / "docs/evaluations/dynamic_portability_20261102"


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    destination = DOC / "artifact_manifest.json"
    files = set(DOC.rglob("*")) | set((BASE / "scripts/evaluation").glob("*.py"))
    files |= {BASE / "core/multiplex5000/multiplex5000/parameter_observations.py",
              BASE / "core/multiplex5000/tests/current/test_parameter_observations.py",
              BASE / "core/multiplex5000/pyproject.toml"}
    files |= {BASE / p for p in (".gitignore", "README.md", "CLAUDE.md",
                                "core/multiplex5000/README.md", "core/multiplex5000/docs/WORKFLOW.md",
                                "core/multiplex5000/docs/INTEGRATION.md", "docs/EXTERNAL_PACKAGES.md",
                                "docs/PLAN_2026-10-31_MALDI_SEQUENTIAL.md", "docs/PORTABILITY_EVALUATION_USER_GUIDE.md")}
    for root in (BASE / "outputs/handoff_completion_20261102", BASE / "outputs/maldi_completion_20261102"):
        for path in root.rglob("*"):
            if not any(part in {"__pycache__", "proposals"} for part in path.relative_to(root).parts):
                files.add(path)
    for path in RUN.rglob("*"):
        rel = path.relative_to(RUN)
        if (rel.parts[0] in {"typing_port", "releases"} or
                rel.parts[:2] == ("remote", "spleen-code") or "__pycache__" in rel.parts):
            continue
        if path.suffix in {".lock", ".tmp", ".pyc"} or path.name == "observations.lock":
            continue
        files.add(path)
    # Explicit failed-run evidence; the remaining incomplete tree is never a release.
    for stage in (RUN / "releases").glob("*/.06_typing.incomplete-*"):
        for name in ("failure.json", "pancreas/run.json", "pancreas/subpanel_audit.json",
                     "pancreas/build_cells_primary.log", "pancreas/subpanel.log", "pancreas/refine.log",
                     "pancreas_scale.parquet", "pancreas/refinement/representation_params.parquet",
                     "pancreas/refinement/representation_provenance.json"):
            files.add(stage / name)
    artifacts = []
    for path in sorted(files):
        if not path.is_file() or path == destination:
            continue
        artifacts.append({"id": str(path.relative_to(BASE)), "path": str(path),
                          "bytes": path.stat().st_size, "sha256": sha(path)})
    checks = []
    for path in sorted(RUN.rglob("receipt_checks.csv")):
        if "releases" in path.relative_to(RUN).parts:
            continue
        with path.open() as handle:
            for row in csv.DictReader(handle):
                checks.append({**row, "audit": str(path.relative_to(BASE))})
    manifest = {"schema": "dynamic-portability-artifacts/1",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "repository": str(BASE), "evaluation_root": str(RUN),
                "artifacts": artifacts, "scientific_receipt_checks": checks,
                "failed_stage_primary_bindings": json.loads((RUN / "protein_primary_expected_universe.json").read_text())["rows"],
                "coverage": "All review outputs including handoff follow-up and source inspections, workflow documentation, recorder, scripts and SQLite export/DB; selected failed-run diagnostics; scientific inputs/outputs retained by receipt checks",
                "exclusions": ["self hash (avoid recursion)", "Python caches and transient lock files",
                               "inactive, unvalidated scientific proposal patches are not deliverable implementation",
                               "raw scientific data and large release/scratch trees are not duplicated; follow receipt checks",
                               "source checkout trees are identified by Git pins and parameter-catalog hashes; later spleen archive is hashed"],
                "status_semantics": "unavailable expected hashes are not passes; original pLN mismatch has a separate exact historical-input recovery; no scientific approval is implied"}
    destination.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{len(artifacts)} review artifacts; {len(checks)} scientific receipt checks")


if __name__ == "__main__":
    main()
