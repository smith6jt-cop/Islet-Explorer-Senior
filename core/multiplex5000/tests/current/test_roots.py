"""The roots, spelled once in phenocycler/roots.py: tracked inventories follow the checkout, data
roots follow the site file.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_roots.py -q
"""

from __future__ import annotations

import importlib
import json
import warnings
from pathlib import Path

import pytest

from phenocycler import roots

REPO = Path(__file__).resolve().parents[4]


@pytest.fixture
def reload_roots(monkeypatch):
    """Reload ``roots`` under a test's environment, then once more under the original one."""
    def _reload():
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            module = importlib.reload(roots)
        return module, caught
    yield _reload
    monkeypatch.undo()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        importlib.reload(roots)


def _site(tmp_path, **fields) -> Path:
    path = tmp_path / "site.json"
    path.write_text(json.dumps({"schema": roots.SITE_SCHEMA, **fields}))
    return path


def test_tracked_inventories_follow_the_checkout():
    """INTEGRATION is versioned with the code, so it is the checkout's, whatever the site says."""
    assert roots.REPO == REPO
    assert roots.INTEGRATION == REPO / "core/phenocycler/data/integration"
    assert (roots.INTEGRATION / "vocab_crosswalk.csv").is_file(), "a tracked file is missing"


def test_data_roots_come_from_the_named_site_file(monkeypatch, tmp_path, reload_roots):
    data = tmp_path / "data"
    data.mkdir()
    monkeypatch.setenv(roots.SITE_FILE_ENV, str(_site(tmp_path, site="test", data_root=str(data))))
    module, caught = reload_roots()
    assert module.DATA_ROOT == data
    assert module.DEVELOPMENT == data / "development"
    assert module.TYPING == data / "development/typing"
    assert module.DATA_ROOT_SOURCE.startswith("site test")
    assert not [w for w in caught if issubclass(w.category, module.SiteFallbackWarning)]
    # The tracked inventories still follow the code, not the site.
    assert module.INTEGRATION == REPO / "core/phenocycler/data/integration"


def test_without_a_site_file_the_checkout_is_used_and_said_so(monkeypatch, tmp_path, reload_roots):
    monkeypatch.delenv(roots.SITE_FILE_ENV, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))  # no ~/.config/multiplex5000/site.json here
    module, caught = reload_roots()
    assert module.SITE is None
    assert module.DATA_ROOT == REPO / "data"
    assert module.DEVELOPMENT == REPO / "data/development"
    assert "checkout fallback" in module.DATA_ROOT_SOURCE
    assert [w for w in caught if issubclass(w.category, module.SiteFallbackWarning)], (
        "the fallback must announce itself")


def test_the_per_user_site_file_is_used_when_no_variable_names_one(monkeypatch, tmp_path,
                                                                   reload_roots):
    data = tmp_path / "data"
    data.mkdir()
    config = tmp_path / ".config" / "multiplex5000"
    config.mkdir(parents=True)
    (config / "site.json").write_text(json.dumps(
        {"schema": roots.SITE_SCHEMA, "site": "per-user", "data_root": str(data)}))
    monkeypatch.delenv(roots.SITE_FILE_ENV, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    module, _ = reload_roots()
    assert module.DATA_ROOT == data and module.SITE["site"] == "per-user"


@pytest.mark.parametrize("fields, message", [
    ({"site": "x", "data_root": "relative/data"}, "absolute"),
    ({"site": "x", "data_root": "/no/such/folder/for/this/test"}, "not a directory"),
    ({"site": "x"}, "data_root"),
    ({"data_root": "/tmp"}, "site"),
    ({"site": "x", "data_root": "/tmp", "tmp_dir": "relative/tmp"}, "tmp_dir"),
])
def test_a_site_file_that_cannot_be_honoured_is_refused(tmp_path, fields, message):
    with pytest.raises(roots.SiteError, match=message):
        roots.read_site(_site(tmp_path, **fields))


def test_a_wrong_schema_or_unreadable_file_is_refused(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"schema": "something-else/9", "site": "x", "data_root": "/tmp"}))
    with pytest.raises(roots.SiteError, match="schema"):
        roots.read_site(bad)
    bad.write_text("{not json")
    with pytest.raises(roots.SiteError, match="JSON"):
        roots.read_site(bad)


def test_a_named_site_file_that_does_not_exist_is_refused(monkeypatch, tmp_path):
    monkeypatch.setenv(roots.SITE_FILE_ENV, str(tmp_path / "missing.json"))
    with pytest.raises(roots.SiteError, match="not a file"):
        roots.site_file()


def test_the_layers_dir_is_read_when_asked_not_at_import(monkeypatch, tmp_path, reload_roots):
    data = tmp_path / "data"
    data.mkdir()
    layers = tmp_path / "layers"
    monkeypatch.setenv(roots.SITE_FILE_ENV, str(_site(tmp_path, site="t", data_root=str(data),
                                                      layers_dir=str(layers))))
    module, _ = reload_roots()
    monkeypatch.delenv("ISLET_LAYERS_DIR", raising=False)
    assert module.layers_dir() == layers  # the site's, without the variable
    monkeypatch.setenv("ISLET_LAYERS_DIR", str(tmp_path / "env"))
    assert module.layers_dir() == tmp_path / "env"  # the variable wins
    monkeypatch.delenv(roots.SITE_FILE_ENV)
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    module, _ = reload_roots()
    monkeypatch.delenv("ISLET_LAYERS_DIR", raising=False)
    assert module.layers_dir() == Path("/data/islet_layers")  # neither: the default


def test_a_site_names_the_recorded_sources_that_live_outside_the_repo(monkeypatch, tmp_path,
                                                                       reload_roots):
    data = tmp_path / "data"
    data.mkdir()
    notes = tmp_path / "notes.txt"
    monkeypatch.setenv(roots.SITE_FILE_ENV, str(_site(
        tmp_path, site="test", data_root=str(data), sources={"donor_notes": str(notes)})))
    module, _ = reload_roots()
    assert module.source("donor_notes") == notes
    assert module.source("pln_blocks") is None, "a source the site does not name is unknown"


def test_a_site_names_where_stage_runs_keep_temporary_files(monkeypatch, tmp_path):
    site = roots.read_site(_site(tmp_path, site="x", data_root="/tmp", tmp_dir="/data/tmp/x"))
    assert site["tmp_dir"] == Path("/data/tmp/x")
    monkeypatch.setattr(roots, "SITE", site)
    assert roots.tmp_dir() == Path("/data/tmp/x")
    monkeypatch.setattr(roots, "SITE", {**site, "tmp_dir": None})
    assert roots.tmp_dir() is None, "no tmp_dir: Python's own temporary folder"
    monkeypatch.setattr(roots, "SITE", None)
    assert roots.tmp_dir() is None


def test_a_source_given_as_a_relative_path_is_refused(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    path = _site(tmp_path, site="test", data_root=str(data), sources={"donor_notes": "notes.txt"})
    with pytest.raises(roots.SiteError, match="sources"):
        roots.read_site(path)
