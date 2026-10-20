#!/usr/bin/env python
"""Does this environment belong to THIS checkout, and where does it read data from?  `make env-check`.

Two things go wrong silently when several checkouts share one virtual environment, and both
happened here on 2026-09-21:

* an editable install records an ABSOLUTE path, so whichever checkout last ran
  ``pip install -e`` owns ``import phenocycler`` for every checkout sharing the .venv -- it was a
  deleted worktree, 56 commits stale, and ``from phenocycler.cohort import DONOR_EXCLUSIONS``
  returned a donor list that had been superseded five days earlier (and on 2026-10-03 it happened
  again: the shared .venv was re-installed from the v3 worktree);
* console scripts record the interpreter that installed them in their shebang, so six of them
  (``pip``, ``phenocycler``, ``tissuumaps_server``, ...) still named that deleted worktree and
  the TissUUmaps unit had been running from a binary that no longer existed.

The environment checked is the one running this script, so ``make env-check`` checks the repo
``.venv`` and ``make env-check PY=/data/envs/multiplex5000/bin/python`` checks the MultiPlex5000
environment (phase D1).

The import is probed from ``/tmp`` with an empty environment, which is how scripts, notebooks and
cron jobs actually import -- pytest's own ``pythonpath = ["."]`` would hide the fault.  A shebang
counts as this environment's when the directory holding its interpreter resolves to this
environment's ``bin/`` (so ``core/phenocycler/../../.venv/bin/python`` counts); resolving the whole
interpreter path instead would collapse every venv on the same base Python into one.  Binaries
(``ruff``) and non-python scripts are skipped.

Finally it reports where ``phenocycler.roots`` takes its data roots from: a site file
(``core/phenocycler/sites/``), or the checkout fallback. The fallback is reported, not failed --
it is right for CI and a fresh clone.
"""
from __future__ import annotations

import glob
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = Path(sys.executable)
BIN = Path(sys.prefix) / "bin"


def _probe(code: str) -> subprocess.CompletedProcess:
    env = {"HOME": os.environ.get("HOME", "/"), "PATH": "/usr/bin:/bin"}
    if os.environ.get("MPLX5K_SITE_FILE"):
        env["MPLX5K_SITE_FILE"] = os.environ["MPLX5K_SITE_FILE"]
    return subprocess.run([str(PY), "-c", code], cwd="/tmp", env=env,
                          capture_output=True, text=True)


def import_resolves_here() -> bool:
    expect = str(REPO / "core" / "phenocycler" / "phenocycler")
    out = _probe("import os, phenocycler; print(os.path.dirname(phenocycler.__file__))")
    got = out.stdout.strip()
    print(f"environment {sys.prefix}")
    print(f"phenocycler -> {got or out.stderr.strip()}")
    if got != expect:
        print(f"WRONG CHECKOUT: expected {expect}")
        return False
    return True


def shebangs_are_this_interpreter() -> bool:
    here = os.path.realpath(BIN)
    bad = []
    for f in sorted(glob.glob(str(BIN / "*"))):
        try:
            with open(f, "rb") as handle:
                first = handle.readline(1024)
        except OSError:
            continue
        if not first.startswith(b"#!") or b"python" not in first:
            continue
        interp = first[2:].decode(errors="replace").split()[0]
        if os.path.realpath(os.path.dirname(interp)) != here:
            bad.append(f"{f} -> {interp}")
    if bad:
        print("console scripts on another interpreter:")
        for b in bad:
            print(f"  {b}")
        return False
    print(f"console scripts: all on {BIN}")
    return True


def data_roots_report() -> bool:
    out = _probe("import warnings; warnings.simplefilter('ignore'); "
                 "import phenocycler.roots as r; "
                 "print(getattr(r, 'DATA_ROOT_SOURCE', 'this checkout predates site files'))")
    source = out.stdout.strip() or out.stderr.strip().splitlines()[-1]
    print(f"data roots -> {source}")
    if "checkout fallback" in source:
        print("  (no site file: the data roots follow the installed checkout)")
    return out.returncode == 0


def main() -> int:
    ok = import_resolves_here()
    ok = shebangs_are_this_interpreter() and ok
    ok = data_roots_report() and ok
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
