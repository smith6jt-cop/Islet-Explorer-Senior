"""
Shared pytest fixtures and collection helpers.

The full Xenium Analysis test suite mixes three tiers:

1. **Fast unit tests** (``tests/test_xenium_explorer_export.py``,
   ``tests/test_repo_health.py``) -- no scanpy required, always safe in CI.
2. **Structural notebook tests** (``tests/test_preprocessing_refactor.py``,
   structural parts of ``tests/test_deg_fixes.py``) -- parse notebook JSON
   only, but depend on specific notebooks being present.  They auto-skip
   when the notebook is missing from this checkout.
3. **Functional tests** (``@pytest.mark.slow``) -- load real data and run
   scanpy pipelines.  Skipped by default so the suite runs anywhere, and
   explicitly enabled on HiPerGator with ``pytest --run-slow``.

This conftest wires up all three.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

# Expose the repo root on sys.path so ``import utils.xenium_explorer_export``
# works from any test module.
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run slow/functional tests that need real Xenium data.",
    )
    parser.addoption(
        "--run-hipergator",
        action="store_true",
        default=False,
        help="Run HiPerGator-specific integration tests.",
    )


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    if not config.getoption("--run-slow"):
        skip_slow = pytest.mark.skip(reason="use --run-slow to execute")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

    if not config.getoption("--run-hipergator"):
        skip_hpg = pytest.mark.skip(reason="use --run-hipergator to execute")
        for item in items:
            if "hipergator" in item.keywords:
                item.add_marker(skip_hpg)


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO


@pytest.fixture(scope="session")
def notebooks_dir(repo_root: Path) -> Path:
    return repo_root / "notebooks"
