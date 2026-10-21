"""Every DuckDB read in this repository closes its result.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_duckdb_reads_close_their_results.py -q

``fetchone()`` and ``fetchmany()`` leave a result open until a later fetch returns no rows.
DuckDB 1.x let the next piece of work on the connection proceed anyway; DuckDB 2.0 refuses it --
measured on 2.0.0 (2026-10-21): an open result stops ``register()``, ``unregister()``,
``con.sql()`` and relation reads such as ``read_parquet()``, while a new ``execute()`` replaces
it. The distance-profile builder relied on the old behaviour by accident (a one-row count read
with ``fetchone()``, then ``register()`` of the next donor's cells), and its tests failed on 2.0.
A one-row read is written ``fetchall()[0]``, which closes the result on either version.
The stricter behaviour is the one this pipeline wants: an unfinished read is an error, not
something later work proceeds past.

The guard holds DuckDB code only: a tracked Python file that names DuckDB -- an import, a type,
or, for a helper that is handed a connection without importing the module
(``scripts/kintsugi/distance.py``), its text. Other code may use ``fetchone()`` as it likes. A
read that drains its result is allowed: the test of a ``while`` loop, which ends only when a
fetch returns no rows, or a call marked ``# drained`` where a later fetch reads to the end.
"""

from __future__ import annotations

import ast
import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
PARTIAL = {"fetchone", "fetchmany"}
DRAINED = "# drained"
_DUCKDB = re.compile(r"duckdb", re.IGNORECASE)


def is_duckdb_code(source: str) -> bool:
    return bool(_DUCKDB.search(source))


def partial_reads(source: str) -> list[int]:
    """The lines of every fetchone()/fetchmany() call that nothing drains."""
    tree = ast.parse(source)
    in_loop_test = {id(node) for loop in ast.walk(tree) if isinstance(loop, ast.While)
                    for node in ast.walk(loop.test)}
    lines = source.splitlines()
    found = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr in PARTIAL) or id(node) in in_loop_test:
            continue
        if any(DRAINED in line for line in lines[node.lineno - 1:node.end_lineno]):
            continue
        found.append(node.lineno)
    return sorted(found)


@pytest.mark.parametrize("snippet, flagged", [
    ("n = con.execute(q).fetchone()\n", [1]),
    ("rows = con.execute(q).fetchmany(5)\n", [1]),
    ("n, = con.execute(\n    q).fetchone()\n", [1]),
    ("n = con.execute(q).fetchall()[0]\n", []),
    ("while (row := con.fetchone()) is not None:\n    use(row)\n", []),
    ("while batch := con.fetchmany(1000):\n    use(batch)\n", []),
    ("for row in iter(con.fetchone, None):\n    use(row)\n", []),
    ("first = con.execute(q).fetchone()  # drained: the rest is read below\n"
     "rest = con.fetchall()\n", []),
])
def test_what_counts_as_a_read_left_open(snippet, flagged):
    assert partial_reads(snippet) == flagged


def test_only_duckdb_code_is_held_to_it():
    assert is_duckdb_code("import duckdb\n")
    assert is_duckdb_code('"""The DuckDB table function for one export."""\n')
    assert not is_duckdb_code("import sqlite3\nrow = cur.execute(q).fetchone()\n")


def test_no_duckdb_code_leaves_a_result_open():
    listed = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", "--", "*.py"],
                            capture_output=True, check=True).stdout.decode().split("\0")
    here = Path(__file__).resolve()
    held, found = 0, []
    for name in filter(None, listed):
        path = REPO / name
        if path.resolve() == here or not path.is_file():
            continue
        source = path.read_text(encoding="utf-8")
        if not is_duckdb_code(source):
            continue
        held += 1
        lines = source.splitlines()
        found += [f"{name}:{n}: {lines[n - 1].strip()}" for n in partial_reads(source)]
    assert held, "no tracked Python names DuckDB; the guard is looking in the wrong place"
    assert not found, ("a partial fetch leaves the DuckDB result open; read one row with "
                       "fetchall()[0]:\n" + "\n".join(found))
