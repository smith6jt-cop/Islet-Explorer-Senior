"""The data roots, spelled once in phenocycler/roots.py.

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_roots.py -q
"""

from __future__ import annotations

from pathlib import Path

from phenocycler import roots

REPO = Path(__file__).resolve().parents[4]


def test_the_roots_are_this_checkouts():
    """The editable install is this checkout (scripts/env_check.py holds the repo to that), so
    the roots derived from the package are the ones beside this test."""
    assert roots.REPO == REPO
    assert roots.INTEGRATION == REPO / "core/phenocycler/data/integration"
    assert roots.DEVELOPMENT == REPO / "data/development"
    assert roots.TYPING == REPO / "data/development/typing"
    assert (roots.INTEGRATION / "vocab_crosswalk.csv").is_file(), "a tracked file is missing"


def test_the_layers_dir_is_read_when_asked_not_at_import(monkeypatch, tmp_path):
    monkeypatch.delenv("ISLET_LAYERS_DIR", raising=False)
    assert roots.layers_dir() == Path("/data/islet_layers")
    monkeypatch.setenv("ISLET_LAYERS_DIR", str(tmp_path))
    assert roots.layers_dir() == tmp_path
