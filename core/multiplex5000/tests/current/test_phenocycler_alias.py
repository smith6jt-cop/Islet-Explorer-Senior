"""The 2026-10-20 rename (MultiPlex5000 D1): ``phenocycler`` is an alias, and the data did not move.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_phenocycler_alias.py -q

``core/phenocycler`` became ``core/multiplex5000`` and the package ``phenocycler`` became
``multiplex5000``. Two things keep everything that still says the old name working, and both are
held here:

* ``phenocycler`` (``core/multiplex5000/phenocycler/__init__.py``) returns the ``multiplex5000``
  modules themselves. A copy would split module state, class identity and monkeypatches between
  two objects that look the same, so these tests ask for identity, not equality.
* The data stays at ``core/phenocycler/data`` until phase H; ``core/multiplex5000/data`` is a
  tracked link to it, so a path built from the package root reaches the same files as before.

The import checks run in a fresh interpreter with this checkout's ``core/multiplex5000`` on
``PYTHONPATH``, so they test this tree whatever the environment has installed, and nothing they
import or reload leaks into the rest of the suite.
"""

from __future__ import annotations

import ast
import os
import subprocess
import sys
import textwrap
from pathlib import Path

PKG_ROOT = Path(__file__).resolve().parents[2]
REPO = Path(__file__).resolve().parents[4]
PACKAGE = PKG_ROOT / "multiplex5000"


def _run(code: str, *args: str) -> subprocess.CompletedProcess:
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    env["PYTHONPATH"] = str(PKG_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    command = [sys.executable, *args] if args else [sys.executable, "-c", textwrap.dedent(code)]
    return subprocess.run(command, cwd="/tmp", env=env, capture_output=True, text=True,
                          timeout=300)


def test_the_old_name_returns_the_same_objects():
    proc = _run("""
        import importlib, sys
        import multiplex5000, phenocycler
        assert phenocycler is multiplex5000
        for name in ("artifacts", "cohort", "roots", "integration", "integration.transform"):
            old = importlib.import_module("phenocycler." + name)
            new = importlib.import_module("multiplex5000." + name)
            assert old is new, name
            assert new.__spec__.name == "multiplex5000." + name, new.__spec__.name
        from phenocycler.integration.transform import Transform as A
        from multiplex5000.integration.transform import Transform as B
        assert A is B
        from phenocycler.cohort import DONOR_EXCLUSIONS
        assert DONOR_EXCLUSIONS is multiplex5000.cohort.DONOR_EXCLUSIONS
        twice = [n for n, m in sys.modules.items()
                 if n.startswith("phenocycler.") and getattr(m, "__name__", "").startswith("phenocycler")]
        assert not twice, f"loaded a second time under the old name: {twice}"
        print("ok")
    """)
    assert proc.returncode == 0 and proc.stdout.strip() == "ok", proc.stderr


def test_reload_through_the_old_name_reruns_the_module():
    proc = _run("""
        import importlib
        import phenocycler.cohort as cohort
        kept = cohort.DONOR_EXCLUSIONS
        cohort.DONOR_EXCLUSIONS = None
        again = importlib.reload(cohort)
        assert again is cohort and cohort.DONOR_EXCLUSIONS == kept
        print("ok")
    """)
    assert proc.returncode == 0 and proc.stdout.strip() == "ok", proc.stderr


def test_a_missing_module_is_still_missing():
    proc = _run("""
        try:
            import phenocycler.no_such_module
        except ModuleNotFoundError as exc:
            print(exc.name)
    """)
    assert proc.stdout.strip() == "phenocycler.no_such_module", proc.stderr


def test_python_dash_m_with_the_old_name_runs_the_new_module():
    proc = _run("", "-m", "phenocycler.pipeline", "run", "--help")
    assert proc.returncode == 0, proc.stderr
    assert "--through" in proc.stdout


def test_nothing_inside_multiplex5000_imports_the_alias():
    offenders = []
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module]
            if any(n == "phenocycler" or n.startswith("phenocycler.") for n in names):
                offenders.append(f"{path.relative_to(PKG_ROOT)}:{node.lineno}")
    assert not offenders, f"the package imports its own alias: {offenders}"


def test_the_package_data_is_a_link_to_where_the_data_stays():
    link = PKG_ROOT / "data"
    assert link.is_symlink(), "core/multiplex5000/data must be the tracked link, not a directory"
    assert os.readlink(link) == "../phenocycler/data"
    assert link.resolve() == (REPO / "core" / "phenocycler" / "data").resolve()
    assert (link / "integration" / "section_modes.csv").is_file()


def test_the_old_directory_holds_only_the_data():
    tracked = subprocess.run(["git", "ls-files", "core/phenocycler"], cwd=REPO,
                             capture_output=True, text=True, check=True).stdout.split()
    assert tracked, "core/phenocycler/data holds the tracked integration inputs"
    stray = [p for p in tracked
             if p != "core/phenocycler/.gitignore" and not p.startswith("core/phenocycler/data/")]
    assert not stray, f"code left behind in core/phenocycler: {stray}"
