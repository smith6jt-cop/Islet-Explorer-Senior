#!/usr/bin/env python3
"""Capture a fresh protein attempt's sources/config/input bindings; never run a stage.

The source roster starts from a named prior snapshot and includes the compatibility
modules and verifier sources added to the attempt-2 capture recipe. Use a new archive,
snapshot and bindings path for every attempt. This CLI did not create the historical
attempt-2 artifacts, which retain their original ad hoc capture provenance.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

EXTRA_SOURCES = (
    "core/multiplex5000/multiplex5000/_version.py",
    "core/multiplex5000/phenocycler/__init__.py",
    "scripts/evaluation/verify_protein_remediation.py",
    "scripts/evaluation/check_protein_producer_snapshot.py",
)


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fresh_target(path):
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"capture target already exists: {path}")
    return path.resolve()


def capture(*, repo, project_path, project_sha256, declared_inputs, settings,
            bind_input, base_snapshot, snapshot, archive, bindings, attempt, code):
    """Write only new artifacts after validating the declaration and source roster."""
    if attempt < 1:
        raise ValueError("attempt must be positive")
    targets = (snapshot, archive, bindings)
    if len(set(targets)) != len(targets) or any(path.exists() or path.is_symlink() for path in targets):
        raise FileExistsError("capture targets must be distinct, fresh paths; existing artifacts are immutable")
    if snapshot.is_relative_to(archive) or bindings.is_relative_to(archive):
        raise ValueError("snapshot and binding files must be outside the source archive")
    if sha(project_path) != project_sha256:
        raise ValueError("project changed while its configuration was loaded")
    names = [item.name for item in declared_inputs]
    if not names or len(set(names)) != len(names):
        raise ValueError("declared input names must be nonempty and unique")
    prior_digest = sha(base_snapshot)
    prior = json.loads(base_snapshot.read_text())
    paths = [row["path"] for row in prior["sources"]]
    if len(paths) != len(set(paths)):
        raise ValueError("base snapshot has duplicate source paths")
    sources = []
    for relative in sorted(set(paths) | set(EXTRA_SOURCES)):
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts or not (repo / path).resolve().is_relative_to(repo):
            raise ValueError(f"source path leaves the repository: {relative}")
        sources.append({"path": relative, "sha256": sha(repo / path)})
    bound = [bind_input(item) for item in declared_inputs]
    if [row.get("name") for row in bound] != names:
        raise ValueError("binder changed the declared input identity")
    if sha(project_path) != project_sha256:
        raise ValueError("project changed while input bindings were captured")
    if sha(base_snapshot) != prior_digest:
        raise ValueError("base source snapshot changed while input bindings were captured")
    if any(sha(repo / row["path"]) != row["sha256"] for row in sources):
        raise ValueError("producer source changed while input bindings were captured")
    creator = Path(__file__).resolve()
    producer = {**code, "attempt": attempt, "created_utc": datetime.now(timezone.utc).isoformat(),
                "prior_snapshot": {"path": str(base_snapshot), "sha256": prior_digest},
                "creator": {"path": str(creator), "sha256": sha(creator)}, "sources": sources}
    archive.mkdir(parents=True, exist_ok=False)
    for row in sources:
        target = archive / row["path"]
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / row["path"], target)
        if sha(target) != row["sha256"]:
            raise ValueError("producer source changed during archive copy; retain this partial capture and use new targets")
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    with snapshot.open("x") as handle:
        handle.write(json.dumps(producer, indent=2) + "\n")
    payload = {"schema": "protein-prelaunch-bindings/1", "attempt": attempt,
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "project": {"path": str(project_path), "sha256": project_sha256},
               "producer_snapshot": {"path": str(snapshot), "sha256": sha(snapshot)},
               "settings": settings, "inputs": bound}
    bindings.parent.mkdir(parents=True, exist_ok=True)
    with bindings.open("x") as handle:
        handle.write(json.dumps(payload, indent=2) + "\n")
    return {"snapshot": payload["producer_snapshot"],
            "input_bindings": {"path": str(bindings), "sha256": sha(bindings), "inputs": len(bound)},
            "archive": str(archive), "stage_execution": "not_started_by_capture"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--base-snapshot", required=True, type=Path)
    parser.add_argument("--snapshot", required=True, type=Path)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--bindings", required=True, type=Path)
    parser.add_argument("--attempt", required=True, type=int)
    args = parser.parse_args()
    from multiplex5000.project import load
    from multiplex5000.stages import io60_mixed_protein_pilot as stage
    from multiplex5000.tasks import bind

    repo = Path(__file__).resolve().parents[2]
    config = args.project.resolve()
    project_sha = sha(config)
    project = load(config)
    if project.settings.get("typing_route") != "io60_mixed_protein_pilot":
        raise ValueError("capture requires the explicit io60_mixed_protein_pilot evaluation route")
    def git(*command):
        return subprocess.check_output(["git", *command], cwd=repo, text=True).strip()
    result = capture(repo=repo, project_path=config, project_sha256=project_sha,
                     declared_inputs=stage.inputs(project), settings=stage.settings(project), bind_input=bind,
                     base_snapshot=args.base_snapshot.resolve(), snapshot=fresh_target(args.snapshot),
                     archive=fresh_target(args.archive), bindings=fresh_target(args.bindings), attempt=args.attempt,
                     code={"branch": git("branch", "--show-current"), "commit": git("rev-parse", "HEAD"),
                           "working_changes": bool(git("status", "--porcelain"))})
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
