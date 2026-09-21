"""`--run` must do less work when less has changed.

Until 2026-09-21 it did not. `run()` called `check()`, printed which nodes were stale, and
then rebuilt every scheduled node anyway -- including the two the same run had just reported
fresh. The staleness was computed and discarded, so any change to any input cost a full
26-node rebuild.

These tests drive `run()` against a synthetic graph whose "builders" are `touch` commands, so
they assert scheduling behaviour without running a builder.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from scripts.senior import build_layers as bl  # noqa: E402


def _node(name, command, inputs, outputs):
    # paths cross this API as strings, not Path
    return bl.Node(name=name, command=tuple(command), inputs=tuple(str(i) for i in inputs),
                   outputs=tuple(str(o) for o in outputs), why="", soft_inputs=())


@pytest.fixture
def graph(tmp_path):
    """upstream -> downstream, both already fresh, with touch as the builder."""
    src = tmp_path / "src.txt"
    src.write_text("x")
    up_out = tmp_path / "up.txt"
    down_out = tmp_path / "down.txt"
    up_out.write_text("u")
    down_out.write_text("d")
    # make both outputs newer than everything they read
    import os, time
    now = time.time()
    os.utime(src, (now - 100, now - 100))
    os.utime(up_out, (now - 50, now - 50))
    os.utime(down_out, (now - 10, now - 10))
    nodes = {
        "up": _node("up", [sys.executable, "-c", "open(r'%s','a').write('!')" % up_out],
                    [src], [up_out]),
        "down": _node("down", [sys.executable, "-c", "open(r'%s','a').write('!')" % down_out],
                      [up_out], [down_out]),
    }
    return nodes, tmp_path, up_out, down_out


def test_a_fresh_node_is_skipped_and_its_builder_never_runs(graph, capsys, tmp_path):
    nodes, repo, up_out, down_out = graph
    before = (up_out.read_text(), down_out.read_text())
    rc = bl.run(nodes, repo)
    out = capsys.readouterr().out
    assert rc == 0
    assert "SKIPPED -- fresh" in out
    assert (up_out.read_text(), down_out.read_text()) == before, (
        "a fresh node's builder ran: the outputs changed")


def test_only_forces_a_fresh_node(graph, capsys, tmp_path):
    """The escape hatch for a node whose inputs are not fully declared."""
    nodes, repo, up_out, _ = graph
    before = up_out.read_text()
    bl.run(nodes, repo, only=["up"])
    assert up_out.read_text() != before, "--only did not force the fresh node to rebuild"


def test_a_stale_node_still_runs(graph, capsys, tmp_path):
    nodes, repo, up_out, down_out = graph
    import os, time
    # make the source newer than up's output -> up is stale, down depends on it
    os.utime(Path(nodes["up"].inputs[0]), None)
    before = up_out.read_text()
    bl.run(nodes, repo)
    assert up_out.read_text() != before, "a stale node was skipped"


def test_downstream_of_a_skipped_fresh_node_is_not_refused(graph, capsys, tmp_path):
    """A skipped node is fresh, so anything reading it may proceed.

    The refusal exists for a node whose input is STALE and not being rebuilt. Treating a
    skip as 'not rebuilt' would refuse the whole graph below every fresh node.
    """
    nodes, repo, _, _ = graph
    bl.run(nodes, repo)
    assert "REFUSED" not in capsys.readouterr().out
