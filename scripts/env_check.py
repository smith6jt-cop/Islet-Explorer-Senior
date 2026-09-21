#!/usr/bin/env python
"""Does this .venv belong to THIS checkout?  `make env-check`.

Two things go wrong silently when several checkouts share one virtual environment, and both
happened here on 2026-09-21:

* an editable install records an ABSOLUTE path, so whichever checkout last ran
  ``pip install -e`` owns ``import phenocycler`` for every checkout sharing the .venv -- it was a
  deleted worktree, 56 commits stale, and ``from phenocycler.cohort import DONOR_EXCLUSIONS``
  returned a donor list that had been superseded five days earlier;
* console scripts record the interpreter that installed them in their shebang, so six of them
  (``pip``, ``phenocycler``, ``tissuumaps_server``, ...) still named that deleted worktree and
  the TissUUmaps unit had been running from a binary that no longer existed.

The import is probed from ``/tmp`` with an empty environment, which is how scripts, notebooks and
cron jobs actually import -- pytest's own ``pythonpath = ["."]`` would hide the fault.  Shebangs
are compared by RESOLVED interpreter path, so ``core/phenocycler/../../.venv/bin/python`` counts as
this interpreter and binaries (``ruff``) and non-python scripts are skipped.
"""
from __future__ import annotations

import glob
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY = REPO / ".venv" / "bin" / "python"


def import_resolves_here() -> bool:
    expect = str(REPO / "core" / "phenocycler" / "phenocycler")
    probe = (
        "import os, phenocycler; print(os.path.dirname(phenocycler.__file__))"
    )
    env = {"HOME": os.environ.get("HOME", "/"), "PATH": "/usr/bin:/bin"}
    out = subprocess.run([str(PY), "-c", probe], cwd="/tmp", env=env,
                         capture_output=True, text=True)
    got = out.stdout.strip()
    print(f"phenocycler -> {got or out.stderr.strip()}")
    if got != expect:
        print(f"WRONG CHECKOUT: expected {expect}")
        return False
    return True


def shebangs_are_this_interpreter() -> bool:
    py = os.path.realpath(PY)
    bad = []
    for f in sorted(glob.glob(str(REPO / ".venv" / "bin" / "*"))):
        try:
            with open(f, "rb") as handle:
                first = handle.readline(1024)
        except OSError:
            continue
        if not first.startswith(b"#!") or b"python" not in first:
            continue
        interp = first[2:].decode(errors="replace").split()[0]
        if os.path.realpath(interp) != py:
            bad.append(f"{f} -> {interp}")
    if bad:
        print("console scripts on another interpreter:")
        for b in bad:
            print(f"  {b}")
        return False
    print(f"console scripts: all on {py}")
    return True


def main() -> int:
    ok = import_resolves_here()
    ok = shebangs_are_this_interpreter() and ok
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
