"""No large data file is tracked at two paths.

On 2026-09-17 three were: ``panel_audit.csv`` (220,755 B), ``panel_audit_thymus.csv`` (215,647 B)
and ``XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv`` (1,091,372 B), each committed both under
``core/panels/`` and under ``analysis/xenium/``. They arrived that way from two separate subtree
absorptions -- ``core/panels`` on 2026-08-14 and ``analysis/xenium`` on 2026-08-28 -- and each copy
had exactly one commit, so neither was ever the maintained one.

The cost is not the disk. It is that a correction has to be applied twice and silently is not: the
HiPerGator session's ``990c5b2`` hand-edited KCNQ1, JAG1 and DMD in BOTH copies to keep them in
step, and the next such edit is one someone forgets. A reader then has two files with the same name
and no way to tell which the pipeline used.

Scoped to files at or above ``MIN_BYTES`` deliberately: small identical files (empty ``__init__.py``,
a shared ``.gitignore``, a stub README) are duplicated for good reasons and are not this defect.
What this catches is a data artifact maintained in two places.
"""
from __future__ import annotations

import hashlib
import subprocess
from collections import defaultdict
from pathlib import Path

#: Asked of git rather than counted with ``parents[n]``. Counting is how this file was wrong the
#: first time -- ``parents[3]`` lands on ``core/``, not the repo root, and a test rooted there finds
#: no duplicates because it cannot see ``analysis/``. It passed, including against a duplicate
#: deliberately re-added to prove it fired. A test that silently inspects nothing is worse than none.
REPO = Path(subprocess.run(["git", "rev-parse", "--show-toplevel"],
                           cwd=Path(__file__).resolve().parent,
                           capture_output=True, text=True, check=True).stdout.strip())

#: Below this, identical content at two paths is ordinary and not a maintenance hazard.
MIN_BYTES = 64 * 1024

#: Paths whose duplication is deliberate and recorded. Empty on purpose -- add an entry only with
#: the reason the second copy must exist, never to silence this test.
ALLOWED: dict[str, str] = {}


def test_the_repo_root_is_the_repo_root():
    """Guards the guard: everything below is scoped by REPO, so a wrong REPO passes vacuously."""
    assert (REPO / ".git").exists(), f"{REPO} is not a git root"
    for marker in ("core/phenocycler", "analysis/xenium", "apps/islet", "CLAUDE.md"):
        assert (REPO / marker).exists(), (
            f"{REPO} does not look like the monorepo root -- {marker} is missing, so a "
            f"duplicate-file scan rooted here would not see the trees that matter")


def _tracked() -> list[Path]:
    out = subprocess.run(["git", "ls-files", "-z"], cwd=REPO, capture_output=True, check=True)
    return [REPO / p for p in out.stdout.decode().split("\0") if p]


def test_no_large_data_file_is_tracked_at_two_paths():
    by_digest: dict[tuple[str, int], list[str]] = defaultdict(list)
    for path in _tracked():
        try:
            size = path.stat().st_size
        except OSError:
            continue                      # a tracked path absent from the worktree is not our business
        if size < MIN_BYTES:
            continue
        rel = path.relative_to(REPO).as_posix()
        if rel in ALLOWED:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        by_digest[(digest, size)].append(rel)

    duplicated = {k: sorted(v) for k, v in by_digest.items() if len(v) > 1}
    assert not duplicated, (
        "these files are tracked at more than one path with identical content; keep one and make "
        "the readers resolve it, so a correction cannot be applied to one copy and not the other:\n"
        + "\n".join(f"  {size:>9,} B  {' == '.join(paths)}"
                    for (_, size), paths in sorted(duplicated.items(), key=lambda kv: -kv[0][1]))
    )
