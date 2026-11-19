"""Receipt-bound parameter artifacts without automatic tuning or global storage.

Stages provide actual consumer settings and emitted fits. This module adds
source/input identities and writes a JSON artifact for inspection or later
observation recording. A recording failure returns a retry payload and does not
discard the scientific outputs.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from .artifacts import sha256_file
from .tasks import Task, TaskError


def write_resolved_parameters(task: Task, *, operation: str, context: dict,
                              parameters: list[dict], acquisitions: list[dict],
                              input_names: Iterable[str] = (),
                              producer_paths: Iterable[Path] = (),
                              relative: str = "resolved_parameters.json") -> dict:
    """Write actual stage-boundary values, or return an explicit retryable status.

    Parameters use roles ``fixed``, ``selected_method``, ``fitted`` and
    ``review_criterion``. Unavailable values are null and retain an explanation
    in the supplied parameter/acquisition record. Source files and the selected
    declared inputs are hashed; large input tables can be bound by their upstream
    receipt instead of being hashed twice. The enclosing stage receipt binds
    this artifact, the complete input roster and the producer commit.
    """
    payload = {"schema": "multiplex5000-resolved-parameters/1", "stage": task.stage,
               "release": task.release, "operation": operation, "context": context,
               "parameters": parameters, "acquisitions": acquisitions,
               "stage_receipt": "receipt.json", "input_bindings": [], "producer_sources": [],
               "automatic_application": False, "scientific_approval": "not_granted_by_recording"}
    try:
        for item in parameters:
            if item.get("role") not in {"fixed", "selected_method", "fitted", "review_criterion"}:
                raise ValueError(f"unsupported parameter role for {item.get('name')!r}")
        for name in sorted(set(input_names)):
            path = task.input(name)
            payload["input_bindings"].append({"name": name, "path": str(path),
                                               "sha256": sha256_file(path)})
        for path in sorted({Path(p).resolve() for p in producer_paths} | {Path(__file__).resolve()}):
            payload["producer_sources"].append({"path": str(path), "sha256": sha256_file(path)})
        data = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
        target = task.output(relative)
        target.write_text(data)
        return {"status": "written", "path": relative, "sha256": sha256_file(target)}
    except (OSError, TypeError, ValueError, TaskError) as error:
        # The receipt retains enough context to retry recording after storage or
        # serialization is repaired; no scientific output is deleted or changed.
        return {"status": "retryable_failure", "path": relative,
                "error": f"{type(error).__name__}: {error}", "retry_payload": payload}
