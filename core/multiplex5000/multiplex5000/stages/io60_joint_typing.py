"""Run the packaged IO60 protein and RNA typing routes in one Stage 06.

The assays retain separate cell identities and outputs. A same-slide relationship
does not create a same-cell link here; registration and integration do that later.
"""

from __future__ import annotations

from pathlib import Path

from .. import gates
from ..tasks import Input, Task
from . import io60_mixed_protein_pilot as protein
from . import io60_rna_typing_pilot as rna

STAGE = "06_typing"


class IO60JointTypingError(ValueError):
    """The two assay routes cannot share one declared release."""


def _same_release(project) -> None:
    upstream = protein._upstream(project)
    if Path(upstream).resolve() != project.release_dir.resolve():
        raise IO60JointTypingError(
            "joint typing needs the protein route to read this project's release; "
            "Stage 00, 04 and 05 must describe the same selected inputs")


def settings(project) -> dict:
    _same_release(project)
    combined = {"scope": "io60_protein_and_rna_review_candidate",
                "paired_cell_links": 0}
    for route in (protein, rna):
        for name, value in route.settings(project).items():
            if name == "scope":
                continue
            if name in combined and combined[name] != value:
                raise IO60JointTypingError(f"assay routes disagree on setting {name}")
            combined[name] = value
    return combined


def inputs(project) -> list[Input]:
    _same_release(project)
    declared: dict[str, Input] = {}
    for route in (protein, rna):
        for item in route.inputs(project):
            prior = declared.get(item.name)
            if prior is not None and prior != item:
                raise IO60JointTypingError(
                    f"assay routes disagree on declared input {item.name}: "
                    f"{prior.path} versus {item.path}")
            declared[item.name] = item
    return list(declared.values())


def run(task: Task) -> dict:
    protein_summary = protein.run(task)
    rna_summary = rna.run(task)
    return {"protein": protein_summary, "rna": rna_summary,
            "paired_cell_links": 0,
            "scope": "protein and RNA cells typed separately in one release"}


def gate(task: Task, summary: dict) -> gates.Gate:
    protein_gate = protein.gate(task, summary["protein"])
    rna_gate = rna.gate(task, summary["rna"])
    figures = [*protein_gate.figures, *rna_gate.figures]
    if len(figures) != len(set(figures)):
        raise IO60JointTypingError("protein and RNA gate figures have colliding names")
    checks = [gates.Check(f"protein: {check.name}", check.value, check.expect,
                          check.passed, check.blocking, check.detail)
              for check in protein_gate.checks]
    checks += [gates.Check(f"RNA: {check.name}", check.value, check.expect,
                           check.passed, check.blocking, check.detail)
               for check in rna_gate.checks]
    checks.append(gates.Check("same-cell links", summary["paired_cell_links"],
                              "0 until registration and integration",
                              summary["paired_cell_links"] == 0))
    return gates.Gate(STAGE, "IO60 pancreas and pLN protein and RNA typing evidence "
                      "from one release. Assay cells remain separate; pLN protein labels "
                      "remain review candidates.", checks, figures)
