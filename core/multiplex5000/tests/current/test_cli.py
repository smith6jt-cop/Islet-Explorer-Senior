"""``mplx5k``: run a stage, see the release, read and accept a gate, describe the graph (D3).

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_cli.py -q
"""

from __future__ import annotations

import json

import pytest

from multiplex5000 import cli


def _mplx5k(capsys, *argv) -> tuple[int, str, str]:
    code = cli.main(list(argv))
    out = capsys.readouterr()
    return code, out.out, out.err


def test_version(capsys):
    with pytest.raises(SystemExit) as stop:
        cli.main(["--version"])
    assert stop.value.code == 0 and "multiplex5000" in capsys.readouterr().out


def test_run_status_show_and_accept(capsys, synthetic_project):
    path = str(synthetic_project())
    code, out, _ = _mplx5k(capsys, "stage", "run", "00_inventory", "--project", path)
    assert code == 0 and "its gate passes" in out and "gate accept 00_inventory" in out
    code, out, _ = _mplx5k(capsys, "stage", "run", "00_inventory", "--project", path)
    assert code == 0 and "already complete" in out
    code, out, _ = _mplx5k(capsys, "status", "--project", path)
    assert "00_inventory" in out and "complete; gate passes; not accepted" in out
    code, out, _ = _mplx5k(capsys, "gate", "show", "00_inventory", "--project", path)
    assert "acquisition ids are unique" in out and "gate.html" in out
    code, out, _ = _mplx5k(capsys, "gate", "accept", "00_inventory", "--by", "tester",
                           "--project", path)
    assert code == 0 and "accepted by tester" in out
    code, out, _ = _mplx5k(capsys, "status", "--project", path)
    assert "accepted by tester" in out


def test_status_keeps_research_comparisons_separate_from_release_gates(
        capsys, synthetic_project):
    path = str(synthetic_project(stages=("00_inventory", "04_rna_ingest", "06_typing")))
    code, out, _ = _mplx5k(capsys, "status", "--project", path)
    assert code == 0
    assert "release gates do not block separate analysis" in out
    assert "class-specific comparisons: need Stage 06 output" in out
    assert "available source measurements can be compared" in out
    assert "accepting Stage 06 does not run comparative analysis" in out


def test_accepting_a_failing_gate_is_refused_with_why(capsys, synthetic_project):
    path = str(synthetic_project(notes__txt="IO60panc2nd Notes\n6374: pancTail no pLN, sequential "
                                            "Xenium\n",
                                 notes_copy__md="IO60panc2nd Notes\n6374: pancTail no pLN, "
                                                "sequential Xenium\n"))
    _mplx5k(capsys, "stage", "run", "00_inventory", "--project", path)
    code, _, err = _mplx5k(capsys, "gate", "accept", "00_inventory", "--by", "tester",
                           "--project", path)
    assert code == 1 and "blocking" in err


def test_the_graph_spec_is_what_the_snakefile_reads(capsys, synthetic_project):
    path = str(synthetic_project())
    code, out, _ = _mplx5k(capsys, "graph-spec", "--project", path)
    spec = json.loads(out)
    [stage] = spec["stages"]
    assert stage["id"] == "00_inventory" and stage["receipt"].endswith("00_inventory/receipt.json")
    assert {"donor_notes", "pln_blocks", "cohort_spec"} <= set(stage["inputs"])
    assert "previous_gate_accepted" not in stage["inputs"], "stage 00 is the first"
    assert spec["target_receipt"] == stage["receipt"]


def test_the_graph_says_where_snakemake_must_come_from(capsys, monkeypatch, synthetic_project):
    monkeypatch.setattr(cli.roots, "SITE", None)
    monkeypatch.setattr(cli.shutil, "which", lambda name: None)
    code, _, err = _mplx5k(capsys, "graph", "--project", str(synthetic_project()))
    assert code == 2 and "requirements-graph.txt" in err
