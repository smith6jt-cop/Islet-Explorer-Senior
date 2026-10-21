"""Gates: every stage ends in checks and figures, and the next stage waits for its acceptance.

A stage's gate is ``<stage>/gate/``: ``gate.json`` (every check, its value, what it must be, and
whether it blocks), the stage's figures and ``gate.html`` to read them in. A person accepts it
with ``mplx5k gate accept <stage> --by <name>``, which writes ``gate/accepted.json`` with who,
when and the digest of the gate they accepted. Acceptance is refused while a blocking check
fails. The next stage, on either route, needs that file; publishing a stage again moves the old
gate and its acceptance aside with the old output, so a changed result is accepted again.

A stage the project file marks automatic accepts itself when every blocking check passes. None
is automatic unless the maintainer says so. Receipts record; gates are what stop progression.
MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import html
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from .artifacts import sha256_file

GATE_DIR = "gate"
GATE_FILE = "gate.json"
ACCEPTED_FILE = "accepted.json"
GATE_SCHEMA = "multiplex5000-gate/1"


class GateError(RuntimeError):
    """A gate that is missing, failing, unaccepted, or changed since it was accepted."""


@dataclass(frozen=True)
class Check:
    """One check: what was measured, what it must be, whether it passed, and whether it blocks.

    A reported (non-blocking) check is shown to the person accepting the gate but cannot stop it.
    """

    name: str
    value: Any
    expect: str
    passed: bool
    blocking: bool = True
    detail: str = ""


@dataclass
class Gate:
    stage: str
    summary: str
    checks: list[Check] = field(default_factory=list)
    figures: list[str] = field(default_factory=list)   # paths relative to the gate folder

    @property
    def passed(self) -> bool:
        return all(c.passed for c in self.checks if c.blocking)

    @property
    def failing(self) -> list[Check]:
        return [c for c in self.checks if c.blocking and not c.passed]


def write(stage_dir: Path, gate: Gate, *, release: str) -> Path:
    """Write ``gate.json`` and ``gate.html`` into ``<stage_dir>/gate``; figures are already there."""
    folder = Path(stage_dir) / GATE_DIR
    folder.mkdir(parents=True, exist_ok=True)
    for figure in gate.figures:
        if not (folder / figure).is_file():
            raise GateError(f"{gate.stage}: the gate names a figure that is not there: {figure}")
    record = {"schema": GATE_SCHEMA, "stage": gate.stage, "release": release,
              "summary": gate.summary, "passed": gate.passed,
              "checks": [asdict(c) for c in gate.checks], "figures": list(gate.figures)}
    path = folder / GATE_FILE
    path.write_text(json.dumps(record, indent=2, sort_keys=True, default=str) + "\n")
    (folder / "gate.html").write_text(_html(record))
    return path


def _html(record: dict) -> str:
    rows = []
    for c in record["checks"]:
        # a reported check is a note whatever its value; only a blocking one passes or fails
        state = ("pass" if c["passed"] else "FAIL") if c["blocking"] else "note"
        rows.append(f"<tr class='{state.lower()}'><td>{state}</td><td>{html.escape(c['name'])}</td>"
                    f"<td>{html.escape(json.dumps(c['value'], default=str))}</td>"
                    f"<td>{html.escape(c['expect'])}</td>"
                    f"<td>{'blocks' if c['blocking'] else 'reported'}</td>"
                    f"<td>{html.escape(c['detail'])}</td></tr>")
    figures = "".join(f"<figure><img src='{html.escape(f)}' alt='{html.escape(f)}'>"
                      f"<figcaption>{html.escape(f)}</figcaption></figure>"
                      for f in record["figures"])
    verdict = ("every blocking check passes" if record["passed"]
               else "a blocking check fails; it cannot be accepted")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Gate {html.escape(record['stage'])}</title>
<style>
:root {{ color-scheme: light dark; --fail: #b3261e; --pass: #2e7d32; --note: #6d6d6d; }}
body {{ font: 15px/1.5 system-ui, sans-serif; margin: 2rem auto; max-width: 70rem;
        padding: 0 1rem; background: Canvas; color: CanvasText; }}
table {{ border-collapse: collapse; width: 100%; }} td, th {{ padding: .35rem .5rem;
        border-bottom: 1px solid #8884; vertical-align: top; text-align: left; }}
tr.fail td:first-child {{ color: var(--fail); font-weight: 600; }}
tr.pass td:first-child {{ color: var(--pass); }} tr.note td:first-child {{ color: var(--note); }}
img {{ max-width: 100%; }} figure {{ margin: 1.5rem 0; }}
</style></head><body>
<h1>Gate: {html.escape(record['stage'])}</h1>
<p>Release <code>{html.escape(record['release'])}</code>. {html.escape(record['summary'])}</p>
<p><strong>{verdict}.</strong> Accept with
<code>mplx5k gate accept {html.escape(record['stage'])} --by &lt;name&gt;</code>.</p>
<table><tr><th></th><th>check</th><th>value</th><th>must be</th><th></th><th>detail</th></tr>
{''.join(rows)}</table>
{figures}
</body></html>
"""


def load(stage_dir: Path) -> dict:
    path = Path(stage_dir) / GATE_DIR / GATE_FILE
    if not path.is_file():
        raise GateError(f"{Path(stage_dir).name}: no gate at {path}; the stage has not published")
    return json.loads(path.read_text())


def digest(stage_dir: Path) -> str:
    return sha256_file(Path(stage_dir) / GATE_DIR / GATE_FILE)


def accept(stage_dir: Path, *, by: str, note: str = "") -> Path:
    """Record that ``by`` accepted this stage's gate. Refused while a blocking check fails."""
    if not by or not str(by).strip():
        raise GateError("an acceptance names who accepts it (--by)")
    record = load(stage_dir)
    failing = [c["name"] for c in record["checks"] if c["blocking"] and not c["passed"]]
    if failing:
        raise GateError(f"{record['stage']}: cannot accept while blocking check(s) fail: {failing}")
    receipt = Path(stage_dir) / "receipt.json"
    accepted = {"stage": record["stage"], "release": record["release"], "accepted_by": str(by),
                "accepted_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "gate_sha256": digest(stage_dir), "note": note,
                "receipt_sha256": sha256_file(receipt) if receipt.is_file() else None}
    path = Path(stage_dir) / GATE_DIR / ACCEPTED_FILE
    path.write_text(json.dumps(accepted, indent=2, sort_keys=True) + "\n")
    return path


def accepted(stage_dir: Path) -> Optional[dict]:
    """The acceptance, if there is one and it is of the gate that is there now."""
    path = Path(stage_dir) / GATE_DIR / ACCEPTED_FILE
    if not path.is_file():
        return None
    record = json.loads(path.read_text())
    if record.get("gate_sha256") != digest(stage_dir):
        raise GateError(f"{Path(stage_dir).name}: the gate changed after it was accepted; "
                        "accept the current one")
    return record


def require_accepted(stage_dir: Path) -> dict:
    """The previous stage's acceptance, or an error that says what is missing."""
    stage = Path(stage_dir).name
    record = accepted(stage_dir)
    if record is None:
        try:
            gate = load(stage_dir)
        except GateError:
            raise GateError(f"{stage} has not published; run it first") from None
        state = "passes" if gate["passed"] else "fails a blocking check"
        raise GateError(f"{stage}'s gate {state} and has not been accepted: read "
                        f"{Path(stage_dir) / GATE_DIR / 'gate.html'}, then "
                        f"`mplx5k gate accept {stage} --by <name>`")
    return record


def figure_path(task, name: str) -> Path:
    """Where a stage writes a gate figure: ``<staging>/gate/<name>``."""
    folder = Path(task.staging) / GATE_DIR
    folder.mkdir(parents=True, exist_ok=True)
    return folder / name
