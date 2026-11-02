"""Parameter artifacts preserve provenance and leave scientific results intact on failure."""

import json
from pathlib import Path

from multiplex5000 import tasks
from multiplex5000.resolved_parameters import write_resolved_parameters


def test_write_records_hashes_null_unavailability_and_receipt_context(tmp_path):
    source = tmp_path / "input.json"
    source.write_text('{"source": "fixture"}')
    task = tasks.Task("06_typing", "test", tmp_path / "stage", {"upstream": tasks.Input("upstream", source)})
    result = write_resolved_parameters(
        task, operation="test", context={"modality": "RNA"},
        parameters=[{"name": "selected", "role": "selected_method", "value": "existing_method"}],
        acquisitions=[{"acquisition_id": "xenium:fixture:region", "fitted": None,
                       "unavailable_reason": "no measured signal"}],
        input_names=("upstream",), producer_paths=(Path(__file__),))
    assert result["status"] == "written"
    path = task.staging / result["path"]
    assert result["sha256"] == tasks.sha256_file(path)
    artifact = json.loads(path.read_text())
    assert artifact["input_bindings"][0]["sha256"] == tasks.sha256_file(source)
    assert artifact["acquisitions"][0]["fitted"] is None
    assert artifact["scientific_approval"] == "not_granted_by_recording"
    assert not artifact["automatic_application"]


def test_recording_failure_is_retryable_and_preserves_scientific_output(tmp_path, monkeypatch):
    task = tasks.Task("06_typing", "test", tmp_path, {})
    science = tmp_path / "calls.parquet"
    science.write_bytes(b"scientific-output-fixture")
    original = Path.write_text

    def failing_write(self, *args, **kwargs):
        if self.name == "resolved_parameters.json":
            raise OSError("storage fixture failure")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", failing_write)
    result = write_resolved_parameters(task, operation="test", context={}, parameters=[], acquisitions=[])
    assert result["status"] == "retryable_failure"
    assert "storage fixture failure" in result["error"]
    assert result["retry_payload"]["operation"] == "test"
    assert science.read_bytes() == b"scientific-output-fixture"
