"""The notebook route: one notebook per stage, calling the same code the graph calls.

    from multiplex5000 import notebooks as nb
    ctx = nb.begin("00_inventory")      # refuses if the previous stage's gate is not accepted
    nb.run(ctx)                         # the same tasks.run_stage the graph runs
    nb.show(ctx)                        # the gate's checks and figures, inline
    nb.accept(ctx, by="<name>")         # or: mplx5k gate accept 00_inventory --by <name>

The project is ``io60`` unless ``$MPLX5K_PROJECT`` names another (an id or a path), and the
releases root is the site's unless ``$MPLX5K_RELEASES_ROOT`` says otherwise, so a notebook runs
unchanged on the workstation, on HiPerGator and in the tests.

MultiPlex5000 D3 (2026-10-20).
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from . import gates, stages, tasks
from . import project as projects


@dataclass(frozen=True)
class Context:
    stage: str
    project: projects.Project

    @property
    def folder(self) -> Path:
        return self.project.release_dir / self.stage


def begin(stage: str, project: str | None = None) -> Context:
    """Load the project, and refuse now if the previous stage's gate has not been accepted."""
    name = project or os.environ.get("MPLX5K_PROJECT") or "io60"
    root = os.environ.get("MPLX5K_RELEASES_ROOT")
    loaded = projects.load(name, releases_root=Path(root) if root else None)
    stages.get(stage)
    previous = loaded.previous(stage)
    if previous is not None:
        gates.require_accepted(loaded.release_dir / previous)
    return Context(stage, loaded)


def run(ctx: Context, *, force: bool = False) -> tasks.Result:
    return tasks.run_stage(stages.get(ctx.stage), ctx.project, force=force)


def show(ctx: Context) -> None:
    """The gate's checks as a table and its figures, inline."""
    from IPython.display import HTML, Image, display

    record = gates.load(ctx.folder)
    rows = "".join(
        f"<tr><td>{('ok' if c['passed'] else 'FAIL') if c['blocking'] else 'note'}</td>"
        f"<td>{c['name']}</td><td>{c['value']}</td><td>{c['detail']}</td></tr>"
        for c in record["checks"])
    display(HTML(f"<p>{record['summary']}</p><table>{rows}</table>"))
    for figure in record["figures"]:
        display(Image(filename=str(ctx.folder / gates.GATE_DIR / figure)))


def accept(ctx: Context, *, by: str, note: str = "") -> Path:
    return gates.accept(ctx.folder, by=by, note=note)
