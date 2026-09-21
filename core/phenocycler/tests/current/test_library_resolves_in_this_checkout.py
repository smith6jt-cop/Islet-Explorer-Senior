"""``import phenocycler`` must read THIS checkout, not another one on the machine.

This repository is worked in through git worktrees, and every worktree's ``.venv`` is a
symlink to one shared environment. That environment carries exactly one editable install of
``phenocycler``, and an editable install records an ABSOLUTE path. So whichever checkout ran
``pip install -e core/phenocycler`` last owns ``import phenocycler`` for every worktree at
once, and nothing announces it.

Measured 2026-09-21, before this test existed: the mapping pointed at
``/home/smith6jt/wt-Q-ci/core/phenocycler/phenocycler``, a worktree 56 commits behind. Its
``cohort.py`` still carried ``6579`` in ``DONOR_EXCLUSIONS`` and had no ``REINSTATED`` dict
at all, five days after the maintainer lifted that exclusion. A plain
``from phenocycler.cohort import filter_eligible_donors`` therefore returned a 20-donor
cohort where the tree on disk says 21 -- silently, with the right-looking answer.

WHY THIS TEST USES A SUBPROCESS, and why the obvious version of it is useless.

The first draft asserted on ``phenocycler.__file__`` in the test process. It passed with the
bug deliberately reintroduced. ``core/phenocycler/pyproject.toml`` sets
``pythonpath = ["."]``, so pytest puts this checkout's ``core/phenocycler`` on ``sys.path``
before the editable finder is ever consulted -- the in-process import is shadowed by the
test runner's own configuration and can never see the misconfiguration. That is also why
CI has never caught this: the CI job runs pytest exactly the same way.

So the condition has to be measured where it actually bites: a fresh interpreter, a working
directory outside the repository, and no inherited ``PYTHONPATH``. That is the environment
every script, notebook and cron job on this machine imports from.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]

PROBE = textwrap.dedent(
    """
    import json, sys
    out = {}
    try:
        import phenocycler
        out["phenocycler"] = phenocycler.__file__
        import phenocycler.integration as _i
        out["integration"] = _i.__file__
        from phenocycler import cohort
        out["excluded"] = sorted(cohort.DONOR_EXCLUSIONS)
        out["has_reinstated"] = hasattr(cohort, "REINSTATED")
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    print(json.dumps(out))
    """
)


@pytest.fixture(scope="module")
def probe(tmp_path_factory) -> dict:
    """What a plain `import phenocycler` resolves to outside this test runner."""
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    neutral = tmp_path_factory.mktemp("neutral")
    proc = subprocess.run(
        [sys.executable, "-c", PROBE],
        cwd=str(neutral),
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert proc.returncode == 0, f"probe failed:\n{proc.stdout}\n{proc.stderr}"
    return json.loads(proc.stdout.strip().splitlines()[-1])


def _under_repo(path: str) -> bool:
    return REPO.resolve() in Path(path).resolve().parents


def test_phenocycler_resolves_inside_the_checkout_this_test_lives_in(probe):
    assert "error" not in probe, probe["error"]
    assert _under_repo(probe["phenocycler"]), (
        f"`import phenocycler`, from a neutral directory with no PYTHONPATH, resolved to\n"
        f"    {probe['phenocycler']}\n"
        f"which is NOT under the checkout this test lives in:\n"
        f"    {REPO}\n\n"
        "Every script that imports phenocycler without a sys.path workaround is reading "
        "another checkout's code -- including its donor exclusions. Repoint it:\n"
        "    .venv/bin/python -m pip install -e core/phenocycler --no-deps\n"
        "or edit MAPPING in "
        ".venv/lib/python*/site-packages/__editable___phenocycler_analysis_*_finder.py"
    )


def test_the_integration_subpackage_resolves_to_the_same_checkout(probe):
    """A split mapping is worse than a wrong one: two halves of one package, two trees."""
    assert "error" not in probe, probe["error"]
    assert _under_repo(probe["integration"]), probe["integration"]
    assert (
        Path(probe["phenocycler"]).resolve().parent
        == Path(probe["integration"]).resolve().parent.parent
    ), (
        "phenocycler and phenocycler.integration resolved to DIFFERENT trees:\n"
        f"  phenocycler             {probe['phenocycler']}\n"
        f"  phenocycler.integration {probe['integration']}"
    )


def test_the_donor_exclusions_are_the_ones_in_this_working_tree(probe):
    """The consequence the resolution bug actually had, asserted directly.

    Read from the file rather than compared to a literal: hardcoding the donor list here
    would make a cohort decision look like a test regression, which is the mistake
    ``test_cohort.py`` already records having made once.
    """
    assert "error" not in probe, probe["error"]
    source = (REPO / "core/phenocycler/phenocycler/cohort.py").read_text()
    for donor in probe["excluded"]:
        assert f'"{donor}"' in source, (
            f"donor {donor!r} is excluded by the imported module but does not appear in "
            f"{REPO}/core/phenocycler/phenocycler/cohort.py -- the import is coming from "
            "another checkout"
        )
    assert probe["has_reinstated"], (
        "the imported cohort module has no REINSTATED mapping. That attribute was added "
        "when 6579's exclusion was lifted on 2026-09-17; a module without it predates the "
        "decision and is not this working tree's."
    )


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v", "-p", "no:cacheprovider"]))
