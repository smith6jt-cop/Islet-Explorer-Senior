"""The bounded IO60 Stage 06 route keeps unavailable channels and paths honest."""

from pathlib import Path
from types import SimpleNamespace

from multiplex5000 import project
from multiplex5000.stages import io60_mixed_protein_pilot as mixed
from multiplex5000.stages import typing
from multiplex5000.tasks import Input, Task


def test_tracked_mixed_project_selects_the_single_stage_module():
    selected = project.load("io60_mixed_protein_pilot")
    assert selected.stages == ("06_typing",)
    assert selected.automatic_gates == ()
    assert typing._route(selected) is mixed


def test_refused_channel_is_not_read_as_a_zero_scale(tmp_path, monkeypatch):
    path = tmp_path / "scales.csv"
    path.write_text("donor_id,marker,status\n6450,CD8,aligned\n"
                    "6450,CD57,refused\n6534,CD8,aligned\n6534,CD57,refused\n")
    monkeypatch.setattr(mixed, "_cohort_rows", lambda: [{"donor_id": "6450"},
                                                    {"donor_id": "6534"}])
    task = Task(mixed.STAGE, "pilot", tmp_path, {"normalization_scales": Input(
        "normalization_scales", path)})
    assert mixed._marker_roster(task) == ["CD8"]


def test_inner_log_names_published_outputs_but_preserves_execution_command(tmp_path):
    release = tmp_path / "release"
    task = Task(mixed.STAGE, "pilot", release / ".06_typing.incomplete-123", {},
                project=SimpleNamespace(release_dir=release))
    command = ["python", "run_calibration_arm.py", "--out", str(task.staging / "pancreas")]
    recorded = {"stages": [{"log": str(task.staging / "pancreas/subpanel.log"),
                            "command": command}]}
    result = mixed._published_inner_manifest(task, recorded)
    assert result["typing_root"] == str(release / "06_typing/pancreas")
    assert result["scale"] == str(release / "06_typing/pancreas_scale.parquet")
    assert result["stages"][0]["log"] == str(release / "06_typing/pancreas/subpanel.log")
    assert result["stages"][0]["command"] == command
    assert Path(result["stages"][0]["log"]).is_relative_to(release)
