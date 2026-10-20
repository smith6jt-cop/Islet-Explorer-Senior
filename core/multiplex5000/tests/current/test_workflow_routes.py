"""The two routes run the same code: the Snakemake graph and the notebook make the same outputs.

    make test-workflow

These need Snakemake in its own environment (``requirements-graph.txt``; named by
``$MPLX5K_SNAKEMAKE``, the site file's ``snakemake``, or PATH) and the ``workflow`` extra for the
notebook. They carry the ``workflow`` marker, so the default run leaves them out and the
Workflow CI job runs them -- where a missing Snakemake fails, it never skips.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from multiplex5000 import cli, gates

pytestmark = pytest.mark.workflow
NOTEBOOK = Path(__file__).resolve().parents[2] / "notebooks" / "00_inventory.ipynb"


def _snakemake() -> str:
    found = os.environ.get("MPLX5K_SNAKEMAKE") or cli._snakemake(None)
    if not found:
        pytest.fail("the workflow tests need Snakemake: set MPLX5K_SNAKEMAKE or the site's")
    return found


def _graph(project: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "multiplex5000.cli", "graph", "--project",
                           str(project), "--snakemake", _snakemake(), *extra],
                          capture_output=True, text=True, timeout=900)


def _outputs(folder: Path) -> dict:
    receipt = json.loads((folder / "receipt.json").read_text())
    return {o["path"]: o["sha256"] for o in receipt["outputs"]}


def test_the_graph_runs_stage_00_once(synthetic_project, tmp_path):
    project = synthetic_project()
    dry = _graph(project, "--dry-run")
    assert dry.returncode == 0, dry.stderr
    assert "stage_00_inventory" in dry.stdout + dry.stderr
    run = _graph(project)
    assert run.returncode == 0, run.stdout + run.stderr
    folder = tmp_path / "releases" / "r1" / "00_inventory"
    assert gates.load(folder)["passed"]
    again = _graph(project)
    assert "Nothing to be done" in again.stdout + again.stderr


def test_the_hipergator_profile_plans_the_same_jobs(synthetic_project):
    dry = _graph(synthetic_project(), "--profile", "hipergator", "--dry-run")
    assert dry.returncode == 0, dry.stderr
    assert "stage_00_inventory" in dry.stdout + dry.stderr


def test_the_notebook_makes_what_the_graph_makes(synthetic_project, tmp_path, monkeypatch):
    import nbformat
    from nbclient import NotebookClient

    graph_project = synthetic_project(release="by-graph")
    run = _graph(graph_project)
    assert run.returncode == 0, run.stdout + run.stderr
    notebook_project = synthetic_project(release="by-notebook")
    monkeypatch.setenv("MPLX5K_PROJECT", str(notebook_project))
    book = nbformat.read(NOTEBOOK, as_version=4)
    NotebookClient(book, timeout=600, kernel_name="python3").execute()
    errors = [o for c in book.cells if c.cell_type == "code" for o in c.get("outputs", [])
              if o.output_type == "error"]
    assert not errors, errors
    by_graph = _outputs(tmp_path / "releases" / "by-graph" / "00_inventory")
    by_notebook = _outputs(tmp_path / "releases" / "by-notebook" / "00_inventory")
    differ = sorted(p for p in by_graph if by_graph[p] != by_notebook.get(p)
                    and p != "gate/gate.json" and p != "gate/gate.html")
    assert set(by_graph) == set(by_notebook) and not differ, differ


def test_the_notebook_in_git_carries_no_outputs():
    import nbformat

    book = nbformat.read(NOTEBOOK, as_version=4)
    assert all(not c.get("outputs") and c.get("execution_count") is None
               for c in book.cells if c.cell_type == "code")
