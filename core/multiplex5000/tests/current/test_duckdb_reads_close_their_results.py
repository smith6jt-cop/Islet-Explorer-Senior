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
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
PARTIAL = re.compile(r"\.fetch(one|many)\(")


def test_no_tracked_python_leaves_a_duckdb_result_open():
    listed = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", "--", "*.py"],
                            capture_output=True, check=True).stdout.decode().split("\0")
    here = Path(__file__).resolve()
    found = []
    for name in filter(None, listed):
        path = REPO / name
        if path.resolve() == here or not path.is_file():
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if PARTIAL.search(line):
                found.append(f"{name}:{number}: {line.strip()}")
    assert not found, ("a partial fetch leaves the result open; read one row with "
                       "fetchall()[0]:\n" + "\n".join(found))
