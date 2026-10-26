"""Project files: the tracked io60 project loads, and anything the package cannot honour is refused.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_project.py -q
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from multiplex5000 import project as P
from multiplex5000 import roots, stages


def test_the_io60_project_loads_and_names_no_automatic_gate(tmp_path):
    project = P.load("io60", releases_root=tmp_path)
    assert project.release_dir == tmp_path / project.release
    assert project.stages[0] == "00_inventory" and project.previous("00_inventory") is None
    assert project.automatic_gates == (), "no gate is automatic until the maintainer says so"


def test_the_releases_root_is_the_argument_else_the_file_else_the_site(tmp_path, monkeypatch,
                                                                        synthetic_project):
    path = synthetic_project()
    assert P.load(path).releases_root == tmp_path / "releases"
    assert P.load(path, releases_root=tmp_path / "x").releases_root == tmp_path / "x"
    monkeypatch.setenv(roots.RELEASES_ROOT_ENV, str(tmp_path / "env"))
    assert P.load("io60").releases_root == tmp_path / "env"


def _write(tmp_path, **change) -> Path:
    body = {"schema": P.PROJECT_SCHEMA, "id": "t", "release": "r1", "cohort": "io60_phenocycler",
            "stages": ["00_inventory"], "automatic_gates": []}
    body.update(change)
    body = {k: v for k, v in body.items() if v is not None}
    path = tmp_path / "t.yaml"
    path.write_text(yaml.safe_dump(body))
    return path


@pytest.mark.parametrize("change, message", [
    ({"typo": 1}, "unknown key"),
    ({"cohort": None}, "missing key"),
    ({"release": "R 1"}, "release must be"),
    ({"stages": ["05_normalization"]}, "built in phase D5"),
    ({"stages": ["nonsense"]}, "no stage"),
    ({"stages": []}, "at least one"),
    ({"automatic_gates": ["03_protein_measurement"]}, "does not run"),
    ({"releases_root": "relative/path"}, "absolute"),
    ({"cohort": "nobody"}, "no cohort profile"),
    ({"id": "other"}, "file name"),
    ({"sources": {"donor_notes": "notes.txt"}}, "absolute"),
])
def test_a_project_the_package_cannot_honour_is_refused(tmp_path, change, message):
    with pytest.raises((P.ProjectError, stages.StageError), match=message):
        P.load(_write(tmp_path, **change))


def test_a_project_source_wins_over_the_sites(tmp_path, synthetic_project):
    project = P.load(synthetic_project())
    assert project.source("donor_notes") == tmp_path / "sources" / "notes.txt"
    assert project.source("not_a_source") == roots.source("not_a_source")


def test_only_built_stages_can_run_and_the_plan_keeps_its_order():
    assert stages.built() == ["00_inventory", "03_protein_measurement"]
    with pytest.raises(stages.StageError, match="phase D6"):
        stages.get("06_typing")
    order = [s for s, _ in stages.PLAN]
    assert order == sorted(order), "stage ids sort in the order they run"
