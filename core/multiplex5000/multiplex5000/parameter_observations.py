"""Opt-in, project-only evaluation observations in KINTSUGI's SQLite history.

This adapter uses the pinned upstream schema, ParameterRecord and insert operation.
It deliberately does not expose the upstream success/recommendation MCP tools:
those default to global storage, inferred quality zeroes and user approval.
Scientific producers need not import this module. Recording is a separate step.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import sqlite3
import subprocess
import sys
import tempfile
import threading
from contextlib import closing
from pathlib import Path

from multiplex5000.artifacts import sha256_file as file_hash

KINTSUGI_COMMIT = "716d1ebd87db916ab04910c6374e04d2d9ea508e"
ENGINE_FILE = "src/kintsugi/claude/parameter_learning.py"
ENGINE_SHA256 = "9c29be884f5dd9d405033a9224b998778fe3905f99093036849a0b9436be6db5"
SCHEMA = "multiplex5000-parameter-observation/1"
CONTEXT_KEYS = {
    "modality", "assay", "tissue", "panel", "relationship_type", "acquisition_ids",
    "section_ids", "segmentation_ids", "measurement_units", "producer_version",
    "registry_hashes", "panel_hashes", "input_receipts",
}
SETTING_KINDS = {"fixed", "selected_method", "fitted", "review_criterion"}
_ENGINE_LOAD_LOCK = threading.Lock()


class ObservationConflict(ValueError):
    """An immutable history conflict needs a corrected revision, not a retry."""


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value) -> str:
    return hashlib.sha256(canonical(value).encode()).hexdigest()




def load_engine(checkout: Path):
    """Load only the pinned engine, without KINTSUGI image/GPU import side effects."""
    checkout = Path(checkout).resolve(strict=True)
    commit = subprocess.check_output(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
    source = checkout / ENGINE_FILE
    if commit != KINTSUGI_COMMIT or file_hash(source) != ENGINE_SHA256:
        raise ValueError("KINTSUGI engine differs from the recorded commit/content pin")
    name = "_multiplex5000_pinned_kintsugi_learning"
    # Dataclasses require publication before execution, but concurrent MCP calls
    # must not receive that module until its first import has finished.
    with _ENGINE_LOAD_LOCK:
        if name not in sys.modules:
            spec = importlib.util.spec_from_file_location(name, source)
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            try:
                spec.loader.exec_module(module)
            except BaseException:
                del sys.modules[name]
                raise
        return sys.modules[name]


def validate(observation: dict) -> dict:
    """Validate an observation; explicit null means unavailable, never measured zero."""
    value = json.loads(canonical(observation))
    required = {"schema", "observation_id", "revision", "supersedes", "status", "operation",
                "context", "subject", "settings", "quality", "outcome", "evidence"}
    if not isinstance(value, dict) or set(value) != required or value["schema"] != SCHEMA:
        raise ValueError("observation must contain exactly the versioned schema fields")
    for key in ("observation_id", "operation", "outcome"):
        if not isinstance(value[key], str) or not value[key].strip():
            raise ValueError(f"{key} must be nonempty text")
    revision = value["revision"]
    if type(revision) is not int or revision < 1:
        raise ValueError("revision must be a positive integer")
    if (value["supersedes"] is not None and
            (type(value["supersedes"]) is not int or value["supersedes"] != revision - 1)):
        raise ValueError("supersedes must name the preceding revision")
    if (revision == 1) != (value["supersedes"] is None):
        raise ValueError("only the first revision has no predecessor")
    if value["status"] not in {"proposed", "evaluated", "failed", "blocked"}:
        raise ValueError("evaluation cannot record successful/user-approved settings")
    context = value["context"]
    if not isinstance(context, dict) or set(context) != CONTEXT_KEYS:
        raise ValueError("context must explicitly include every compatibility/provenance field")
    for key in CONTEXT_KEYS - {"acquisition_ids", "section_ids", "segmentation_ids",
                               "registry_hashes", "panel_hashes", "input_receipts"}:
        if context[key] is not None and (not isinstance(context[key], str) or not context[key]):
            raise ValueError(f"context.{key} must be text or null (unavailable)")
    for key in ("acquisition_ids", "section_ids", "segmentation_ids"):
        items = context[key]
        if not isinstance(items, list) or any(not isinstance(s, str) or not s for s in items):
            raise ValueError(f"{key} must be a list of explicit identifiers")
        if len(items) != len(set(items)):
            raise ValueError(f"{key} contains repeated identifiers")
    for key in ("registry_hashes", "panel_hashes"):
        if not isinstance(context[key], dict):
            raise ValueError(f"{key} must be an explicit hash mapping (empty if unavailable)")
        for sha in context[key].values():
            _sha(sha)
    for refs in (value["evidence"], context["input_receipts"]):
        if not isinstance(refs, list):
            raise ValueError("evidence/receipts must be lists of artifact references")
        for ref in refs:
            if (not isinstance(ref, dict) or set(ref) != {"id", "path", "sha256"}
                    or not isinstance(ref["id"], str) or not ref["id"]
                    or not isinstance(ref["path"], str) or not ref["path"]):
                raise ValueError("artifact reference needs id, path and sha256")
            _sha(ref["sha256"])
    subject = value["subject"]
    if (not isinstance(subject, dict) or set(subject) != {"kind", "identifiers"}
            or subject["kind"] not in {"markers", "genes", "feature_set", "acquisition_pair",
                                       "geometry", "acquisition"}
            or not isinstance(subject["identifiers"], list)
            or any(not isinstance(s, str) or not s for s in subject["identifiers"])):
        raise ValueError("subject needs a typed identifier list; geometry needs no marker")
    if not isinstance(value["settings"], list):
        raise ValueError("settings must be a list")
    names = set()
    for setting in value["settings"]:
        if (not isinstance(setting, dict) or set(setting) != {"name", "kind", "value", "units", "source"}
                or setting["kind"] not in SETTING_KINDS
                or not isinstance(setting["name"], str) or not setting["name"]
                or not isinstance(setting["source"], str) or not setting["source"]
                or (setting["units"] is not None and not isinstance(setting["units"], str))):
            raise ValueError("setting needs name, kind, value, units and evidence source")
        if setting["name"] in names:
            raise ValueError("duplicate setting name")
        names.add(setting["name"])
    if not isinstance(value["quality"], dict):
        raise ValueError("quality must be a mapping")
    for score in value["quality"].values():
        if score is not None and (type(score) not in (int, float) or not math.isfinite(score)):
            raise ValueError("quality must be finite numbers or null")
    return value


def _sha(value):
    if (not isinstance(value, str) or len(value) != 64
            or any(c not in "0123456789abcdef" for c in value)):
        raise ValueError("expected lowercase SHA256")


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".observation-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class ObservationHistory:
    def __init__(self, project_path: Path, checkout: Path, *, storage: str = "project_only"):
        if storage != "project_only" or not project_path or not Path(project_path).is_absolute():
            raise ValueError("explicit absolute project_path and project_only storage required")
        self.project = Path(project_path).resolve()
        self.module = load_engine(checkout)
        self.project.mkdir(parents=True, exist_ok=True)
        # Resolve all storage paths before upstream initialization (including symlinks).
        self.db = self._inside(".claude_workflow/parameter_learning.db")
        self.engine = self.module.ParameterLearningEngine(str(self.project), global_db=False)
        assert self.engine.global_db_path is None
        with closing(sqlite3.connect(self.db)) as con, con:
            con.execute("""CREATE UNIQUE INDEX IF NOT EXISTS evaluation_observation_revision
                ON parameter_records(json_extract(parameters_json, '$.observation.observation_id'),
                                     json_extract(parameters_json, '$.observation.revision'))
                WHERE json_extract(parameters_json, '$.observation.schema') =
                      'multiplex5000-parameter-observation/1'""")
            # Even another upstream client cannot promote or rewrite an observation.
            for event in ("UPDATE", "DELETE"):
                con.execute(f"""CREATE TRIGGER IF NOT EXISTS evaluation_no_{event.lower()}
                    BEFORE {event} ON parameter_records
                    WHEN json_extract(OLD.parameters_json, '$.observation.schema') = '{SCHEMA}'
                    BEGIN SELECT RAISE(ABORT, 'evaluation observations are immutable'); END""")
            con.execute(f"""CREATE TRIGGER IF NOT EXISTS evaluation_no_approval
                BEFORE INSERT ON parameter_records
                WHEN json_extract(NEW.parameters_json, '$.observation.schema') = '{SCHEMA}'
                     AND NEW.user_approved != 0
                BEGIN SELECT RAISE(ABORT, 'evaluation observations cannot be approved'); END""")

    def _inside(self, relative: str) -> Path:
        path = (self.project / relative).resolve()
        if path == self.project or not path.is_relative_to(self.project):
            raise ValueError("observation storage/export must remain inside the project")
        return path

    def history(self, observation_id: str | None = None) -> list[dict]:
        self._inside(".claude_workflow/parameter_learning.db")
        with closing(sqlite3.connect(self.db)) as con:
            rows = con.execute("""SELECT id, timestamp, parameters_json FROM parameter_records
                WHERE json_extract(parameters_json, '$.observation.schema') = ? ORDER BY id""",
                               (SCHEMA,)).fetchall()
        return [{"record_id": rid, "recorded_at": when, **json.loads(payload)}
                for rid, when, payload in rows
                if observation_id is None or json.loads(payload)["observation"]["observation_id"]
                == observation_id]

    def record(self, observation: dict) -> dict:
        value = validate(observation)
        fingerprint = digest(value)
        lock = self._inside(".claude_workflow/observations.lock")
        with lock.open("a") as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            previous = self.history(value["observation_id"])
            for old in previous:
                if old["observation"]["revision"] == value["revision"]:
                    if old["sha256"] != fingerprint:
                        raise ObservationConflict("id/revision already records different content; append a revision")
                    return {"status": "recorded", "record_id": old["record_id"],
                            "sha256": fingerprint, "duplicate": True}
            if value["revision"] != len(previous) + 1:
                raise ObservationConflict("revision must append to the current history without gaps")
            # Record directly through upstream's insertion API: record_parameters subtracts
            # qualities and normalizes marker strings, so cannot represent missing evidence.
            record = self.module.ParameterRecord(
                tissue_type=value["context"]["tissue"],
                tissue_type_normalized=value["context"]["tissue"],
                marker_name=None, marker_name_normalized=None,
                operation=value["operation"],
                parameters={"observation": value, "sha256": fingerprint},
                quality_before=value["quality"].get("before"),
                quality_after=value["quality"].get("after"), quality_improvement=None,
                user_approved=False, project_path=str(self.project))
            # No image characteristics were measured by this recorder.
            record.characteristics = _UnavailableCharacteristics()
            record.algorithm_version = value["context"]["producer_version"]
            rid = self.engine._insert_record(self.db, record)
        return {"status": "recorded", "record_id": rid, "sha256": fingerprint, "duplicate": False}

    def export(self, relative_path: str = "parameter_observations.json") -> dict:
        path = self._inside(relative_path)
        # Never overwrite SQLite or its lock/journal files with an export.
        if ".claude_workflow" in path.relative_to(self.project).parts or path.suffix != ".json":
            raise ValueError("export needs a .json path outside the database directory")
        rows = self.history()
        atomic_json(path, {"schema": SCHEMA, "storage": "project_only",
                           "automatic_application": False, "recommendations_enabled": False,
                           "kintsugi_commit": KINTSUGI_COMMIT, "records": rows})
        return {"path": str(path), "sha256": file_hash(path), "records": len(rows)}


class _UnavailableCharacteristics:
    def to_dict(self):
        return None


def record_or_queue(project_path: Path, checkout: Path, observation: dict) -> dict:
    """Durable outbox before database access; failures never mutate scientific outputs.

    Invalid observations fail explicitly. Operational failures return a retryable status.
    If even the outbox cannot be written, the response says so; the caller retains its input.
    """
    if not project_path or not Path(project_path).is_absolute():
        raise ValueError("explicit absolute project_path required")
    value = validate(observation)
    project = Path(project_path).resolve()
    pending = project / "parameter_observation_outbox" / f"{digest(value)}.json"
    if not pending.resolve().is_relative_to(project):
        raise ValueError("outbox must remain inside project")
    try:
        atomic_json(pending, value)
        result = ObservationHistory(project, checkout).record(value)
        # Another identical retry may already have removed this shared outbox item.
        pending.unlink(missing_ok=True)
        return result
    except ObservationConflict:
        pending.unlink(missing_ok=True)
        raise
    except (OSError, sqlite3.Error, subprocess.SubprocessError, ValueError) as exc:
        return {"status": "retryable", "error": f"{type(exc).__name__}: {exc}",
                "pending_path": str(pending) if pending.is_file() else None,
                "sha256": digest(value)}


def make_server(project_path: Path, checkout: Path):
    """Project is fixed at launch; MCP clients cannot redirect storage or enable learning."""
    from mcp.server.fastmcp import FastMCP

    server = FastMCP("KINTSUGI evaluation observations")

    @server.tool()
    def record_parameter_observation(observation: dict) -> dict:
        """Record proposed/evaluated/failed/blocked parameters, never approve or apply them."""
        return record_or_queue(project_path, checkout, observation)

    @server.tool()
    def get_parameter_observation_history(observation_id: str | None = None) -> list[dict]:
        """Inspect the project's append-only evaluation history."""
        return ObservationHistory(project_path, checkout).history(observation_id)

    @server.tool()
    def export_parameter_observations(relative_path: str = "parameter_observations.json") -> dict:
        """Export evaluation history inside the project, including unavailable measurements."""
        return ObservationHistory(project_path, checkout).export(relative_path)

    return server


def retry_pending(project_path: Path, checkout: Path) -> list[dict]:
    """Replay predecessors before revisions, independent of hashed outbox filenames."""
    if not project_path or not Path(project_path).is_absolute():
        raise ValueError("explicit absolute project_path required")
    project = Path(project_path).resolve()
    folder = (project / "parameter_observation_outbox").resolve()
    if not folder.is_relative_to(project):
        raise ValueError("outbox must remain inside project")
    pending = []
    for path in folder.glob("*.json"):
        if not path.resolve().is_relative_to(project):
            raise ValueError("outbox entry must remain inside project")
        pending.append((path, validate(json.loads(path.read_text()))))
    pending.sort(key=lambda entry: (entry[1]["observation_id"], entry[1]["revision"]))
    results, blocked = [], set()
    for path, value in pending:
        if value["observation_id"] in blocked:
            result = {"status": "retryable", "error": "predecessor is still queued",
                      "pending_path": str(path), "sha256": digest(value)}
        else:
            result = record_or_queue(project, checkout, value)
        results.append(result)
        if result["status"] == "retryable":
            blocked.add(value["observation_id"])
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--kintsugi", required=True, type=Path)
    parser.add_argument("action", choices=("record", "retry", "history", "export", "serve"))
    parser.add_argument("--input", type=Path)
    args = parser.parse_args()
    if not args.project.is_absolute():
        parser.error("--project must be absolute")
    if args.action == "serve":
        make_server(args.project, args.kintsugi).run(transport="stdio")
        return
    if args.action == "record":
        if args.input is None:
            parser.error("record requires --input")
        result = record_or_queue(args.project, args.kintsugi, json.loads(args.input.read_text()))
    elif args.action == "retry":
        result = retry_pending(args.project, args.kintsugi)
    else:
        history = ObservationHistory(args.project, args.kintsugi)
        result = history.history() if args.action == "history" else history.export()
    print(json.dumps(result, indent=2, allow_nan=False))
    results = result if isinstance(result, list) else [result]
    if any(r.get("status") == "retryable" for r in results):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
