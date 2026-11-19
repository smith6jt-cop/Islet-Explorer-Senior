"""Project context must survive copying a suggested command to another directory."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pytest
import yaml

from multiplex5000 import cli, gates


def _invoke(capsys, *args):
    code = cli.main(list(args))
    output = capsys.readouterr()
    return code, output.out, output.err


def _hint(output, label):
    lines = [line.strip().split(label, 1)[1].strip() for line in output.splitlines()
             if line.strip().startswith(label)]
    assert len(lines) == 1, output
    return shlex.split(lines[0])


def _site(tmp_path, monkeypatch):
    # Pin the chosen site too: the copied command must survive a different shell
    # environment without reading that shell's invalid/default site file.
    site = tmp_path / "selected site 'quoted'.json"
    site.write_text(json.dumps({"schema": "multiplex5000-site/1", "site": "synthetic",
                                "data_root": str(tmp_path)}))
    monkeypatch.setenv("MPLX5K_SITE_FILE", str(site))
    monkeypatch.setattr(cli.roots, "SITE", cli.roots.read_site(site))
    monkeypatch.setattr(cli.roots, "_SITE_PATH", site)
    return site


@pytest.mark.parametrize("root_origin", ["argument", "project", "environment"])
def test_suggested_commands_keep_project_root_and_site_from_another_directory(
        capsys, monkeypatch, tmp_path, synthetic_project, root_origin):
    site = _site(tmp_path, monkeypatch)
    original = synthetic_project(release="selected-r1")
    body = yaml.safe_load(original.read_text())
    selected_root = tmp_path / "selected output 'quoted'"
    wrong_root = tmp_path / "unselected output"
    extra = []
    if root_origin == "argument":
        body["releases_root"] = str(wrong_root)
        extra = ["--releases-root", str(selected_root)]
    elif root_origin == "project":
        body["releases_root"] = str(selected_root)
    else:
        body.pop("releases_root")
        monkeypatch.setenv("MPLX5K_RELEASES_ROOT", str(selected_root))
    project = tmp_path / "named project 'quoted'" / "synthetic.yaml"
    project.parent.mkdir()
    project.write_text(yaml.safe_dump(body))
    argv = ["stage", "run", "00_inventory", "--project", str(project), *extra]
    code, output, error = _invoke(capsys, *argv)
    assert code == 0, error
    stage = selected_root / "selected-r1/00_inventory"
    assert gates.accepted(stage) is None
    gate_page = (stage / "gate/gate.html").read_text()
    assert "mplx5k gate accept" not in gate_page
    assert "project" in gate_page.lower()
    show = _hint(output, "show:")
    accept = _hint(output, "accept after review:")
    for command in (show, accept):
        assert command[command.index("--project") + 1] == str(project.resolve())
        assert command[command.index("--releases-root") + 1] == str(selected_root.resolve())
        assert sys.executable in command
        assert f"MPLX5K_SITE_FILE={site.resolve()}" in command
    # The reuse branch must expose the same review path, without rerunning it.
    code, repeated, error = _invoke(capsys, *argv)
    assert code == 0 and "already complete" in repeated, error
    assert _hint(repeated, "show:") == show
    assert _hint(repeated, "accept after review:") == accept
    neutral = tmp_path / "neutral directory"
    neutral.mkdir()
    other_code = tmp_path / "unselected code"
    (other_code / "multiplex5000").mkdir(parents=True)
    (other_code / "multiplex5000/__init__.py").write_text("raise RuntimeError('wrong producer')\n")
    env = dict(os.environ, MPLX5K_RELEASES_ROOT=str(wrong_root),
               MPLX5K_SITE_FILE=str(tmp_path / "missing-other-site.json"),
               PYTHONPATH=str(other_code),
               PYTHONDONTWRITEBYTECODE="1")
    result = subprocess.run(show, cwd=neutral, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "selected-r1" in result.stdout
    assert gates.accepted(stage) is None
    accept[accept.index("--by") + 1] = "synthetic reviewer"
    result = subprocess.run(accept, cwd=neutral, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert gates.accepted(stage)["accepted_by"] == "synthetic reviewer"
    assert not wrong_root.exists()


def test_context_is_read_only_and_reports_effective_named_project(
        capsys, monkeypatch, tmp_path, synthetic_project):
    site = _site(tmp_path, monkeypatch)
    project = synthetic_project(release="context-r1")
    output = tmp_path / "context output"
    neutral = tmp_path / "neutral"
    neutral.mkdir()
    monkeypatch.chdir(neutral)

    def forbidden(*args, **kwargs):
        pytest.fail("context must not run a stage or accept a gate")

    monkeypatch.setattr(cli.tasks, "run_stage", forbidden)
    monkeypatch.setattr(cli.gates, "accept", forbidden)
    code, text, error = _invoke(capsys, "context", "--project", str(project),
                               "--releases-root", str(output))
    assert code == 0, error
    context = json.loads(text)
    assert context["schema"] == "multiplex5000-context/1"
    assert context["side_effects"] is False
    assert context["receipt_freshness_checked"] is False
    assert context["project"]["id"] == "synthetic"
    assert context["project"]["path"] == str(project.resolve())
    assert context["project"]["release"] == "context-r1"
    assert context["project"]["releases_root"] == str(output.resolve())
    assert context["project"]["release_dir"] == str(output / "context-r1")
    assert context["runtime"]["executable"] == sys.executable
    assert context["runtime"]["site_file"] == str(site.resolve())
    assert context["runtime"]["package_file"]
    assert "revision" in context["runtime"]["git"]
    assert "dirty" in context["runtime"]["git"]
    assert not output.exists()
    assert not (tmp_path / "releases").exists()


def test_automatic_gate_does_not_offer_an_extra_acceptance(
        capsys, tmp_path, synthetic_project):
    project = synthetic_project(automatic=("00_inventory",))
    code, text, error = _invoke(capsys, "stage", "run", "00_inventory",
                               "--project", str(project))
    assert code == 0, error
    assert "accepted automatically" in text
    assert "accept after review:" not in text
    assert gates.accepted(tmp_path / "releases/r1/00_inventory") is not None


def test_runtime_context_reports_wheel_build_without_claiming_live_git(tmp_path, monkeypatch):
    from multiplex5000.runtime_context import runtime_context
    snapshot = tmp_path / "bundled source"
    snapshot.mkdir()
    provenance = {"commit": "fixture-commit", "dirty": True}
    (snapshot / "_build_provenance.json").write_text(json.dumps(provenance))
    monkeypatch.setattr(cli.roots, "REPO", snapshot)
    runtime = runtime_context()
    assert runtime["git"]["available"] is False
    assert runtime["git"]["revision"] is None
    assert runtime["build_provenance"]["record"] == provenance
    assert len(runtime["build_provenance"]["sha256"]) == 64
