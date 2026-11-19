"""Sampled historical fingerprints must never be compared as full-file hashes."""
import importlib.util
import json
import sys
from pathlib import Path

import pytest
from multiplex5000.artifacts import sha256_file, sha256_sampled_file

SPEC = importlib.util.spec_from_file_location("verify_recorded_inputs", Path(__file__).parents[1] / "verify_recorded_inputs.py")
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


def fixture(tmp_path, mode):
    source = tmp_path / "source.bin"
    source.write_bytes(bytes(range(256)) * 4)
    expected = sha256_file(source) if mode == "full" else sha256_sampled_file(source, sample_size_bytes=16, sample_count=3)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"images": [{"qptiff": {
        "path": str(source), "hash_mode": mode, "sha256": expected,
        "sample_size_bytes": 16, "sample_count": 3}}]}))
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"inputs": [
        {"name": "qupath_manifest", "path": str(manifest), "sha256": sha256_file(manifest)},
        {"name": "qptiff:S1", "path": str(source), "recorded_sha256": expected, "bytes": source.stat().st_size}]}))
    return source, manifest, receipt


@pytest.mark.parametrize("mode", ["full", "sampled"])
def test_declared_algorithm_passes_and_full_observation_is_separate(tmp_path, monkeypatch, mode):
    source, manifest, receipt = fixture(tmp_path, mode)
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["verify", "--receipt", str(receipt), "--manifest", str(manifest), "--out", str(out)])
    VERIFY.main()
    row, = json.loads((out / "verification.json").read_text())["inputs"]
    assert row["status"] == "pass" and row["hash_mode"] == "full"
    assert row["sha256"] == sha256_file(source)
    assert row["observed_recorded_mode_hash"] == row["recorded_sha256"]
    assert (row["sha256"] == row["recorded_sha256"]) == (mode == "full")


def test_manifest_cannot_redefine_hash_mode_after_receipt(tmp_path, monkeypatch):
    _, manifest, receipt = fixture(tmp_path, "sampled")
    manifest.write_text(manifest.read_text().replace('"sampled"', '"full"'))
    monkeypatch.setattr(sys, "argv", ["verify", "--receipt", str(receipt), "--manifest", str(manifest), "--out", str(tmp_path / "out")])
    with pytest.raises(ValueError, match="not uniquely bound"):
        VERIFY.main()


def test_changed_bytes_fail_even_with_same_size(tmp_path, monkeypatch):
    source, manifest, receipt = fixture(tmp_path, "full")
    source.write_bytes(b"x" * source.stat().st_size)
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["verify", "--receipt", str(receipt), "--manifest", str(manifest), "--out", str(out)])
    with pytest.raises(SystemExit):
        VERIFY.main()
    assert json.loads((out / "verification.json").read_text())["status"] == "fail"
