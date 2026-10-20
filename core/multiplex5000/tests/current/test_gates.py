"""Gates: written with every check, accepted by a named person, refused while a blocking check
fails, and void once the gate they accepted is not the gate that is there (MultiPlex5000 D3).

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_gates.py -q
"""

from __future__ import annotations

import json

import pytest

from multiplex5000 import gates


def _write(tmp_path, *checks, figures=()):
    stage = tmp_path / "s1"
    (stage / gates.GATE_DIR).mkdir(parents=True)
    for figure in figures:
        (stage / gates.GATE_DIR / figure).write_bytes(b"png")
    gates.write(stage, gates.Gate("s1", "synthetic", list(checks), list(figures)), release="r1")
    return stage


PASS = gates.Check("ids are unique", [], "none", True)
FAIL = gates.Check("ids are unique", ["a"], "none", False)
NOTE = gates.Check("files not found", 3, "reported", False, blocking=False)


def test_a_gate_is_written_with_every_check_and_a_page(tmp_path):
    stage = _write(tmp_path, PASS, NOTE, figures=("coverage.png",))
    record = gates.load(stage)
    assert record["passed"] and record["release"] == "r1"
    assert [c["name"] for c in record["checks"]] == ["ids are unique", "files not found"]
    page = (stage / gates.GATE_DIR / "gate.html").read_text()
    assert "coverage.png" in page and "mplx5k gate accept s1" in page


def test_a_gate_that_names_a_missing_figure_is_refused(tmp_path):
    with pytest.raises(gates.GateError, match="figure"):
        gates.write(tmp_path / "s1", gates.Gate("s1", "x", [PASS], ["missing.png"]), release="r1")


def test_acceptance_is_refused_while_a_blocking_check_fails(tmp_path):
    stage = _write(tmp_path, FAIL, NOTE)
    with pytest.raises(gates.GateError, match="blocking"):
        gates.accept(stage, by="tester")
    assert gates.accepted(stage) is None


def test_a_reported_check_does_not_block(tmp_path):
    stage = _write(tmp_path, PASS, NOTE)
    gates.accept(stage, by="tester")
    assert gates.accepted(stage)["accepted_by"] == "tester"


def test_an_acceptance_names_who_and_binds_the_gate_it_accepted(tmp_path):
    stage = _write(tmp_path, PASS)
    with pytest.raises(gates.GateError, match="who"):
        gates.accept(stage, by=" ")
    gates.accept(stage, by="tester", note="looked at the figure")
    record = gates.accepted(stage)
    assert record["gate_sha256"] == gates.digest(stage) and record["note"]
    path = stage / gates.GATE_DIR / gates.GATE_FILE
    changed = json.loads(path.read_text())
    changed["summary"] = "edited after acceptance"
    path.write_text(json.dumps(changed))
    with pytest.raises(gates.GateError, match="changed after it was accepted"):
        gates.accepted(stage)


def test_require_accepted_says_what_is_missing(tmp_path):
    with pytest.raises(gates.GateError, match="has not published"):
        gates.require_accepted(tmp_path / "s0")
    stage = _write(tmp_path, PASS)
    with pytest.raises(gates.GateError, match="passes and has not been accepted"):
        gates.require_accepted(stage)
    gates.accept(stage, by="tester")
    assert gates.require_accepted(stage)["accepted_by"] == "tester"
