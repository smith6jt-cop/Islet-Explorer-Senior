"""Evaluation infrastructure contracts, using the actual pinned KINTSUGI engine."""
import copy
import json
import sqlite3
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from multiplex5000 import parameter_observations as recording


@pytest.fixture
def observation():
    return {
        "schema": recording.SCHEMA, "observation_id": "registration:pair-1", "revision": 1,
        "supersedes": None, "status": "blocked", "operation": "registration",
        "context": {**dict.fromkeys(recording.CONTEXT_KEYS),
                    "modality": "protein/rna", "assay": "PhenoCycler/Xenium",
                    "tissue": "spleen", "panel": "PC67/5K",
                    "relationship_type": "same_section_after_phenocycler",
                    "acquisition_ids": ["protein:1", "xenium:1"], "section_ids": ["section:1"],
                    "segmentation_ids": [], "measurement_units": "um",
                    "producer_version": "test/1", "registry_hashes": {}, "panel_hashes": {},
                    "input_receipts": []},
        "subject": {"kind": "geometry", "identifiers": []},
        "settings": [{"name": "max_dist", "kind": "fixed", "value": 5,
                      "units": "um", "source": "same_slide_registration.MAX_DIST_UM"}],
        "quality": {"before": None, "after": None},
        "outcome": "No accepted registration receipt", "evidence": [],
    }


@pytest.fixture
def checkout():
    path = Path(__file__).resolve().parents[4] / "external/KINTSUGI"
    if not path.exists():
        pytest.skip("requires external/KINTSUGI at the pinned commit (see recording README)")
    return path


@pytest.mark.parametrize("status", ["proposed", "evaluated", "failed", "blocked"])
def test_observation_states_do_not_imply_approval(observation, status):
    observation["status"] = status
    assert recording.validate(observation)["status"] == status


@pytest.mark.parametrize("mutation", [
    {"status": "approved"}, {"user_approved": True}, {"revision": True},
    {"revision": 3, "supersedes": 1}, {"quality": {"after": float("nan")}},
    {"context": {"tissue": "pancreas"}}, {"subject": {"kind": "marker", "identifiers": []}},
])
def test_rejects_implicit_or_invalid_context(observation, mutation):
    observation.update(mutation)
    with pytest.raises(ValueError):
        recording.validate(observation)


@pytest.mark.data
def test_actual_engine_roundtrip_revisions_and_global_isolation(tmp_path, monkeypatch, checkout, observation):
    home = tmp_path / "home"
    global_db = home / ".kintsugi/parameter_learning.db"
    global_db.parent.mkdir(parents=True)
    global_db.write_bytes(b"global history sentinel")
    monkeypatch.setenv("HOME", str(home))
    project = tmp_path / "evaluation"
    store = recording.ObservationHistory(project, checkout)
    first = store.record(observation)
    assert store.record(observation)["duplicate"]
    changed = copy.deepcopy(observation)
    changed["status"] = "evaluated"
    with pytest.raises(ValueError, match="different content"):
        store.record(changed)
    changed.update(revision=2, supersedes=1)
    second = store.record(changed)
    assert second["record_id"] != first["record_id"]
    restored = recording.ObservationHistory(project, checkout)
    assert len(restored.history()) == 2
    exported = json.loads(Path(restored.export()["path"]).read_text())
    assert [r["observation"] for r in exported["records"]] == [observation, changed]
    with sqlite3.connect(store.db) as con:
        assert con.execute("SELECT marker_name, quality_before, quality_after, characteristics_json, "
                           "user_approved FROM parameter_records").fetchall() == [
                               (None, None, None, "null", 0), (None, None, None, "null", 0)]
        with pytest.raises(sqlite3.IntegrityError, match="immutable"):
            con.execute("UPDATE parameter_records SET user_approved=1")
        with pytest.raises(sqlite3.IntegrityError, match="immutable"):
            con.execute("DELETE FROM parameter_records")
    # Existing KINTSUGI inspection works; its approved-only export stays empty.
    assert store.engine.get_statistics()["project_records"] == 2
    approved = store.engine.export_learned_parameters(project / "upstream-approved.json")
    assert json.loads(approved.read_text())["total_records"] == 0
    assert store.engine._query_similar_records(store.db, "spleen", "CD3e", "registration", 0, True) == []
    assert global_db.read_bytes() == b"global history sentinel"


@pytest.mark.data
def test_concurrent_cold_start_loads_complete_engine(checkout):
    # A fresh interpreter matters: the other engine tests already warm the cache.
    result = subprocess.run([
        sys.executable, "-c", """
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
sys.path.insert(0, sys.argv[1])
from multiplex5000.parameter_observations import load_engine
barrier = Barrier(8)
def load(_):
    barrier.wait()
    engine = load_engine(Path(sys.argv[2]))
    return engine.ParameterLearningEngine.__name__
with ThreadPoolExecutor(max_workers=8) as pool:
    assert list(pool.map(load, range(8))) == ["ParameterLearningEngine"] * 8
""", str(Path(recording.__file__).resolve().parent.parent), str(checkout)],
        capture_output=True, text=True, timeout=30, check=False)
    assert result.returncode == 0, result.stderr


@pytest.mark.data
def test_concurrent_retries_insert_once(tmp_path, checkout, observation):
    store = recording.ObservationHistory(tmp_path, checkout)
    with ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(lambda _: store.record(observation), range(18)))
    assert len({r["record_id"] for r in results}) == 1
    assert len(store.history()) == 1
    with ThreadPoolExecutor(max_workers=6) as pool:
        queued = list(pool.map(
            lambda _: recording.record_or_queue(tmp_path, checkout, observation), range(18)))
    assert all(r["status"] == "recorded" and r["duplicate"] for r in queued)
    assert not list((tmp_path / "parameter_observation_outbox").glob("*.json"))


@pytest.mark.data
def test_recording_failure_is_durable_and_retryable(tmp_path, checkout, observation, monkeypatch):
    output = tmp_path / "scientific.parquet"
    output.write_bytes(b"scientific output remains intact")
    original = recording.ObservationHistory.record
    def broken(*args):
        raise sqlite3.OperationalError("database is locked")
    monkeypatch.setattr(recording.ObservationHistory, "record", broken)
    failed = recording.record_or_queue(tmp_path, checkout, observation)
    assert failed["status"] == "retryable"
    pending = Path(failed["pending_path"])
    assert json.loads(pending.read_text()) == observation
    monkeypatch.setattr(recording.ObservationHistory, "record", original)
    assert recording.record_or_queue(tmp_path, checkout, observation)["status"] == "recorded"
    assert not pending.exists()
    assert output.read_bytes() == b"scientific output remains intact"


@pytest.mark.data
def test_history_conflicts_are_not_queued_as_transient_failures(tmp_path, checkout, observation):
    assert recording.record_or_queue(tmp_path, checkout, observation)["status"] == "recorded"
    observation["outcome"] = "changed content requires a new revision"
    with pytest.raises(recording.ObservationConflict, match="append a revision"):
        recording.record_or_queue(tmp_path, checkout, observation)
    assert not list((tmp_path / "parameter_observation_outbox").glob("*.json"))


@pytest.mark.data
def test_outbox_replays_predecessors_before_later_revisions(tmp_path, checkout, observation, monkeypatch):
    original = recording.ObservationHistory.record
    def unavailable(*args):
        raise sqlite3.OperationalError("temporarily unavailable")
    monkeypatch.setattr(recording.ObservationHistory, "record", unavailable)
    first = copy.deepcopy(observation)
    second = {**observation, "revision": 2, "supersedes": 1, "outcome": "later review"}
    for value in (second, first):
        assert recording.record_or_queue(tmp_path, checkout, value)["status"] == "retryable"
    assert all(r["status"] == "retryable" for r in recording.retry_pending(tmp_path, checkout))
    assert len(list((tmp_path / "parameter_observation_outbox").glob("*.json"))) == 2
    monkeypatch.setattr(recording.ObservationHistory, "record", original)
    assert all(r["status"] == "recorded" for r in recording.retry_pending(tmp_path, checkout))
    history = recording.ObservationHistory(tmp_path, checkout).history()
    assert [r["observation"] for r in history] == [first, second]
    assert recording.retry_pending(tmp_path, checkout) == []


@pytest.mark.data
def test_storage_and_export_cannot_escape_project(tmp_path, checkout, observation):
    with pytest.raises(ValueError):
        recording.ObservationHistory(tmp_path, checkout, storage="global")
    project = tmp_path / "project"
    store = recording.ObservationHistory(project, checkout)
    for target in ("../outside.json", str(tmp_path / "outside.json"),
                   ".claude_workflow/parameter_learning.db", ".claude_workflow/lock.json"):
        with pytest.raises(ValueError):
            store.export(target)
    another = tmp_path / "another"
    another.mkdir()
    (another / ".claude_workflow").symlink_to(project / ".claude_workflow")
    with pytest.raises(ValueError):
        recording.ObservationHistory(another, checkout)


@pytest.mark.data
def test_pin_drift_and_missing_predecessor_refused(tmp_path, checkout, observation):
    store = recording.ObservationHistory(tmp_path, checkout)
    observation.update(revision=2, supersedes=1)
    with pytest.raises(ValueError, match="without gaps"):
        store.record(observation)
    with pytest.raises((ValueError, FileNotFoundError)):
        recording.load_engine(tmp_path / "not-a-checkout")
