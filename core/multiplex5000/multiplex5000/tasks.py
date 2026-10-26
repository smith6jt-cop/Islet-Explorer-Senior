"""Tasks: a stage runs in a staging folder, is published in one rename, and leaves a receipt.

Three of the four lifecycle repairs of the MultiPlex5000 plan live here (the fourth, the
wrong-checkout import, was repaired in D1):

* **A partial output never counts as complete.** A stage writes into
  ``.<stage>.incomplete-<id>`` beside its destination and is renamed into ``<stage>/`` only once
  its outputs, its gate and its receipt are written. A stage is complete when
  ``<stage>/receipt.json`` says so -- never because a folder has something in it, which is what
  the 2026-10-08 audit found the old pipeline doing.
* **--force really reruns.** It recomputes in a new staging folder and replaces the published
  stage in one step; the previous output is moved aside to ``.<stage>.superseded-<time>``, never
  deleted. Without --force a complete stage is left as it is and reported.
* **Every input is declared.** A stage names its inputs; each is bound in the receipt at the
  moment the stage runs, a change while it runs fails it, and a file the stage opens that it did
  not declare fails it too (an audit hook sees every ``open``).

Receipts record; they never stop anyone reading or re-running. Only gates stop progression.
MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import contextvars
import fcntl
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import uuid
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import roots
from .artifacts import sha256_file

RECEIPT = "receipt.json"
FAILURE = "failure.json"
RECEIPT_SCHEMA = "multiplex5000-receipt/1"


class TaskError(RuntimeError):
    """A stage that did not complete, and why."""


@dataclass(frozen=True)
class Input:
    """One declared input of a stage.

    ``hash=False`` binds a large raw file by path, size, modification time and the digest a source
    already records for it, instead of reading every byte again; tables and manifests are hashed.
    """

    name: str
    path: Path
    hash: bool = True
    recorded_sha256: Optional[str] = None
    required: bool = True


def bind(item: Input) -> dict:
    """What a declared input was when the stage read it."""
    path = Path(item.path)
    record: dict[str, Any] = {"name": item.name, "path": str(path), "exists": path.exists(),
                              "recorded_sha256": item.recorded_sha256}
    if not record["exists"]:
        if item.required:
            raise TaskError(f"declared input {item.name} is missing: {path}")
        return record
    if path.is_dir():
        raise TaskError(f"declared input {item.name} is a directory; declare its files: {path}")
    stat = path.stat()
    record.update(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if item.hash:
        record["sha256"] = sha256_file(path)
    return record


def _unchanged(before: dict, item: Input) -> bool:
    path = Path(item.path)
    if not before["exists"]:
        return not path.exists()
    stat = path.stat()
    return (stat.st_size, stat.st_mtime_ns) == (before["bytes"], before["mtime_ns"])


# --------------------------------------------------------------------- what the stage reads ----
_OPENED: contextvars.ContextVar[Optional[set]] = contextvars.ContextVar("mplx5k_opened",
                                                                       default=None)


def _audit(event: str, args: tuple) -> None:
    if event != "open":
        return
    opened = _OPENED.get()
    if opened is None:
        return
    target, mode, flags = (tuple(args) + (None, None, None))[:3]
    if not isinstance(target, (str, bytes, os.PathLike)):
        return
    path = os.path.realpath(os.fsdecode(target))
    if isinstance(mode, str):
        creating = any(c in mode for c in "wax")
    else:
        creating = isinstance(flags, int) and bool(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT))
    # A file the stage itself creates or writes (its outputs, a library's scratch file) may be read
    # back; only a file it found already there is an input.
    opened.add(("written" if creating else "read", path))


sys.addaudithook(_audit)


def audited(func, *args):
    """``func(*args)`` and every file it opened -- for a stage's work done in a pool worker.

    A worker is a fresh interpreter with no record of the running stage, so what it opens would go
    unchecked. The stage hands each worker's record back to the run with :func:`record`."""
    opened: set = set()
    token = _OPENED.set(opened)
    try:
        return func(*args), opened
    finally:
        _OPENED.reset(token)


def record(opened: set) -> None:
    """Add what a worker opened to the running stage's record, so it is checked like its own."""
    current = _OPENED.get()
    if current is not None:
        current.update(opened)


def _always_readable() -> list[str]:
    """Where reading needs no declaration: the interpreter and its libraries, this package's own
    code and knowledge (versioned in the receipt as the commit), system files and caches."""
    home = Path.home()
    places = {sys.prefix, sys.base_prefix, sys.exec_prefix, str(Path(__file__).resolve().parent),
              str(roots.REPO / "core/multiplex5000/profiles"), "/usr", "/lib", "/lib64", "/etc",
              "/proc", "/sys", "/dev", "/var/cache/fontconfig",
              str(home / ".cache"), str(home / ".config/matplotlib"), str(home / ".matplotlib"),
              str(home / ".local"), str(home / ".config/multiplex5000")}
    places.update(p for p in sys.path if p and os.path.isdir(p) and "site-packages" in p)
    return sorted(os.path.realpath(p) for p in places)


def _undeclared(opened: set, inputs: list[Input], staging: Path) -> list[str]:
    declared = {os.path.realpath(i.path) for i in inputs}
    written = {path for kind, path in opened if kind == "written"}
    free = _always_readable() + [os.path.realpath(staging)]
    out = []
    for path in sorted({path for kind, path in opened if kind == "read"}):
        if path in declared or path in written or not os.path.isfile(path):
            continue
        if any(path == root or path.startswith(root.rstrip("/") + "/") for root in free):
            continue
        if _is_package_metadata(path):
            continue
        out.append(path)
    return out


def _is_package_metadata(path: str) -> bool:
    """A distribution's own metadata (``*.dist-info``, or the ``*.egg-info`` an editable install
    leaves beside its source), which importlib reads whenever entry points are scanned."""
    return any(part.endswith((".dist-info", ".egg-info")) for part in Path(path).parts[:-1])


# ----------------------------------------------------------------------------- the task ----
@dataclass
class Task:
    """What a stage's ``run`` and ``gate`` get: its declared inputs, its settings and the staging
    folder it writes into."""

    stage: str
    release: str
    staging: Path
    inputs: dict[str, Input]
    settings: dict = field(default_factory=dict)
    project: Any = None

    def input(self, name: str) -> Path:
        if name not in self.inputs:
            raise TaskError(f"{self.stage} reads {name!r}, which it did not declare")
        return Path(self.inputs[name].path)

    def output(self, relative: str) -> Path:
        path = self.staging / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        return path


@dataclass(frozen=True)
class Result:
    stage: str
    status: str            # "published" or "complete" (already, left as it was)
    path: Path
    superseded: Optional[Path] = None
    gate_passed: Optional[bool] = None
    accepted_automatically: bool = False


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_receipt(stage_dir: Path) -> Optional[dict]:
    path = Path(stage_dir) / RECEIPT
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except ValueError:
        return None


def is_complete(stage_dir: Path) -> bool:
    receipt = read_receipt(stage_dir)
    return bool(receipt) and receipt.get("status") == "complete"


def code_identity() -> dict:
    """The package version and the commit the code came from, and whether its tree was edited."""
    from . import __version__

    out = {"version": __version__, "commit": None, "dirty": None}
    bundled = roots.REPO / "_build_provenance.json"
    if bundled.is_file():
        record = json.loads(bundled.read_text())
        out["commit"], out["dirty"] = record["commit"], record["dirty"]
        return out
    try:
        out["commit"] = subprocess.run(["git", "-C", str(roots.REPO), "rev-parse", "HEAD"],
                                       capture_output=True, text=True, check=True,
                                       timeout=30).stdout.strip()
        status = subprocess.run(["git", "-C", str(roots.REPO), "status", "--porcelain", "--",
                                 "core/multiplex5000", "core/panels/data",
                                 "core/phenocycler/data/integration"],
                                capture_output=True, text=True, check=True, timeout=60).stdout
        out["dirty"] = bool(status.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return out


def _outputs(staging: Path) -> list[dict]:
    files = sorted(p for p in staging.rglob("*") if p.is_file() and p.name != RECEIPT)
    for path in files:
        if not path.resolve().is_relative_to(staging.resolve()):
            raise TaskError(f"an output points outside the stage: {path}")
    return [{"path": str(p.relative_to(staging)), "bytes": p.stat().st_size,
             "sha256": sha256_file(p)} for p in files]


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")


@contextmanager
def _temporary_files(stage_id: str, run_id: str):
    """While a stage runs, Python's temporary files go to a folder of their own under the site's
    ``tmp_dir`` -- and DuckDB's overflow with them (``cells_parquet``) -- so a machine whose home
    disk is small never holds them. Pool workers inherit ``TMPDIR``. The folder is removed when
    the run ends; with no ``tmp_dir`` they stay where Python puts them."""
    root = roots.tmp_dir()
    if root is None:
        yield None
        return
    folder = Path(root) / f"{stage_id}.{run_id}"
    folder.mkdir(parents=True)
    before = os.environ.get("TMPDIR")
    os.environ["TMPDIR"] = str(folder)
    tempfile.tempdir = None              # gettempdir() reads TMPDIR again
    try:
        yield folder
    finally:
        if before is None:
            os.environ.pop("TMPDIR", None)
        else:
            os.environ["TMPDIR"] = before
        tempfile.tempdir = None
        shutil.rmtree(folder, ignore_errors=True)


def _publish(staging: Path, destination: Path) -> Optional[Path]:
    """Rename the staging folder into place; a previous output is moved aside, never deleted."""
    release = destination.parent
    with (release / ".publish.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        aside = None
        if destination.is_dir() and not destination.is_symlink() and not any(destination.iterdir()):
            destination.rmdir()          # an empty placeholder (Snakemake makes one); no output
        if destination.exists() or destination.is_symlink():
            stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
            aside = release / f".{destination.name}.superseded-{stamp}"
            destination.rename(aside)
        staging.rename(destination)
    return aside


def run_stage(stage, project, *, force: bool = False) -> Result:
    """Run one stage of a project's release through the task contract.

    ``stage`` is the stage module (``inputs``, ``run``, ``gate``, and optionally ``settings``);
    ``project`` the loaded project file. The previous stage's gate must have been accepted.
    """
    from . import gates

    stage_id = stage.STAGE
    release_dir = Path(project.release_dir)
    destination = release_dir / stage_id
    if is_complete(destination) and not force:
        return Result(stage_id, "complete", destination)
    previous = project.previous(stage_id)
    if previous is not None:
        gates.require_accepted(release_dir / previous)
    inputs = list(stage.inputs(project))
    settings = dict(stage.settings(project)) if hasattr(stage, "settings") else {}
    release_dir.mkdir(parents=True, exist_ok=True)
    staging = release_dir / f".{stage_id}.incomplete-{uuid.uuid4().hex[:12]}"
    staging.mkdir()
    task = Task(stage_id, project.release, staging, {i.name: i for i in inputs}, settings, project)
    started = utc_now()
    try:
        bound = [bind(i) for i in inputs]
        opened: set = set()
        with _temporary_files(stage_id, staging.name.rsplit("-", 1)[-1]):
            token = _OPENED.set(opened)
            try:
                summary = stage.run(task)
                gate = stage.gate(task, summary)
            finally:
                _OPENED.reset(token)
        if not isinstance(summary, dict):
            raise TaskError(f"{stage_id}.run must return a summary dict")
        stray = _undeclared(opened, inputs, staging)
        if stray:
            raise TaskError(f"{stage_id} read file(s) it did not declare: {stray[:5]}"
                            f"{' ...' if len(stray) > 5 else ''}")
        changed = [i.name for i, b in zip(inputs, bound, strict=True) if not _unchanged(b, i)]
        if changed:
            raise TaskError(f"{stage_id}: input(s) changed while the stage ran: {changed}")
        gates.write(staging, gate, release=project.release)
        receipt = {
            "schema": RECEIPT_SCHEMA, "status": "complete", "stage": stage_id,
            "release": project.release, "project": str(project.path), "started": started,
            "completed": utc_now(), "inputs": bound, "settings": settings,
            "code": code_identity(),
            "environment": {"python": platform.python_version(), "platform": platform.platform(),
                            "tmp_dir": str(roots.tmp_dir()) if roots.tmp_dir() else None},
            "summary": summary, "gate": {"passed": gate.passed, "checks": len(gate.checks)},
            "outputs": _outputs(staging),
        }
        _write_json(staging / RECEIPT, receipt)
        superseded = _publish(staging, destination)
    except BaseException as error:
        if staging.exists():
            _write_json(staging / FAILURE, {
                "status": "failed", "stage": stage_id, "release": project.release,
                "started": started, "failed": utc_now(), "error_type": type(error).__name__,
                "error": str(error), "note": "a failed stage's files are never treated as output"})
        raise
    automatic = False
    if stage_id in project.automatic_gates and gate.passed:
        gates.accept(destination, by="automatic (project file)")
        automatic = True
    return Result(stage_id, "published", destination, superseded, gate.passed, automatic)


def failed_attempts(release_dir: Path) -> list[Path]:
    """Staging folders left by runs that did not publish, newest first."""
    found = [p for p in Path(release_dir).glob(".*.incomplete-*") if p.is_dir()]
    return sorted(found, key=lambda p: p.stat().st_mtime, reverse=True)

