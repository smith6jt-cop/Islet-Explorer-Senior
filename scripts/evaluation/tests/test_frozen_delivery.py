"""Evolved source may be replayed; scientific artifacts cannot be substituted."""
import hashlib
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

spec = importlib.util.spec_from_file_location("frozen_delivery", Path(__file__).parents[1] / "verify_portability_delivery.py")
delivery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(delivery)


def test_only_exact_delivered_source_blob_can_replace_changed_working_source(tmp_path, monkeypatch):
    monkeypatch.setattr(delivery, "BASE", tmp_path)
    monkeypatch.setattr(delivery, "DOC", tmp_path / "docs/evaluation")
    source = tmp_path / "runner.py"
    source.write_bytes(b"new code")
    artifact = {"path": str(source), "sha256": hashlib.sha256(b"original code").hexdigest()}
    calls = []

    def git(args, **kwargs):
        calls.append(args)
        return SimpleNamespace(stdout=b"original code")

    monkeypatch.setattr(delivery.subprocess, "run", git)
    assert delivery.artifact_binding(artifact) == "delivery_git_blob"
    assert calls == [["git", "show", f"{delivery.DELIVERY_REVISION}:runner.py"]]
    artifact["sha256"] = hashlib.sha256(b"different source").hexdigest()
    with pytest.raises(AssertionError):
        delivery.artifact_binding(artifact)


@pytest.mark.parametrize("relative", ["outputs/science.json", "docs/evaluation/parameter_history.json"])
def test_scientific_evidence_never_falls_back_to_git(tmp_path, monkeypatch, relative):
    monkeypatch.setattr(delivery, "BASE", tmp_path)
    monkeypatch.setattr(delivery, "DOC", tmp_path / "docs/evaluation")
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    path.write_bytes(b"changed")
    monkeypatch.setattr(delivery.subprocess, "run", lambda *a, **kw: pytest.fail("Git must not substitute scientific evidence"))
    with pytest.raises(AssertionError, match="frozen scientific evidence"):
        delivery.artifact_binding({"path": str(path), "sha256": hashlib.sha256(b"original").hexdigest()})
