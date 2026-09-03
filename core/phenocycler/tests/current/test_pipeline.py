from __future__ import annotations

from types import SimpleNamespace

import pytest

from phenocycler.artifacts import ContractError, StaleArtifactError
from phenocycler.config import PipelineConfig
from phenocycler.pipeline import (
    DEFAULT_STAGE_ORDER,
    _assert_stages_available,
    _identity_config,
    _selected_stages,
    _validate_export,
)


def test_stage_selection_is_explicit_and_ordered():
    assert _selected_stages(None, None) == DEFAULT_STAGE_ORDER
    assert _selected_stages(None, "expression") == (
        "ingest",
        "geometry",
        "redsea",
        "expression",
    )
    assert _selected_stages(["type"], None) == ("type",)


def test_no_verified_typing_are_part_of_run_identity_config():
    # The upstream workflow must not require a placeholder RESTORE policy while
    # the registry declares that no pair roster has been verified.
    cfg = PipelineConfig()
    identity = _identity_config(cfg)
    assert identity["typing_model"] == {
        "mode": "unavailable",
        "reason": "NO_VERIFIED_TYPING",
    }
    assert identity["simultaneous_broad_typing"] == {
        "mode": "blocked",
        "reason": "NO_VERIFIED_TYPING",
    }
    assert identity["reference_controls"] == {
        "mode": "blocked",
        "reason": "NO_VERIFIED_TYPING",
    }
    assert identity["phenotyping"] == {
        "mode": "blocked",
        "reason": "NO_VERIFIED_TYPING",
    }


def test_production_typing_fails_closed_without_verified_pairs():
    import pandas as pd

    from phenocycler.marker_registry import load_registry
    from phenocycler.production_typing import type_production_cells

    with pytest.raises(ContractError, match="NO_VERIFIED_TYPING"):
        type_production_cells(pd.DataFrame(), marker_registry=load_registry())


def test_type_stage_has_only_the_simultaneous_interface(monkeypatch, tmp_path):
    from phenocycler import pipeline

    evidence = object()
    assignments = object()
    observed: dict[str, object] = {}
    monkeypatch.setattr(
        pipeline,
        "read_single_partition",
        lambda *_args, **_kwargs: evidence,
    )

    def fake_type(frame, **kwargs):
        observed.update(frame=frame, **kwargs)
        return assignments

    monkeypatch.setattr(pipeline, "type_calibrated_cells", fake_type)
    monkeypatch.setattr(
        pipeline,
        "write_donor_partition",
        lambda frame, *_args, **_kwargs: observed.update(written=frame),
    )
    context = SimpleNamespace(
        donors=("D1",),
        config=SimpleNamespace(
            marker_evidence_dir=tmp_path / "evidence",
            assignments_dir=tmp_path / "assignments",
        ),
        registry="markers",
        broad_gate_policy="gates",
    )

    pipeline._run_type(context)

    assert observed["frame"] is evidence
    assert observed["marker_registry"] == "markers"
    assert observed["broad_gate_policy"] == "gates"
    assert set(observed) == {
        "frame",
        "marker_registry",
        "broad_gate_policy",
        "written",
    }
    assert observed["written"] is assignments


def test_malformed_export_manifest_fails_with_a_contract_error(tmp_path):
    manifests = tmp_path / "manifests"
    manifests.mkdir()
    (manifests / "qupath_export.json").write_text("{}\n")
    context = SimpleNamespace(
        config=SimpleNamespace(manifests_dir=manifests),
        run_id="run",
    )

    with pytest.raises(ContractError, match="missing content_id"):
        _validate_export(context)


def test_stage_inputs_are_validated_before_the_runner(monkeypatch, tmp_path):
    from phenocycler import pipeline

    class StaleInput:
        def validate_current(self, *, mode):
            raise StaleArtifactError("upstream changed")

    context = SimpleNamespace(
        registry=SimpleNamespace(production_status="verified"),
        stage_manifest_path=lambda stage: tmp_path / f"{stage}.json",
        output_root=lambda stage: tmp_path / stage,
    )
    ran = False

    def runner(_context):
        nonlocal ran
        ran = True

    monkeypatch.setattr(pipeline, "_stage_inputs", lambda *_: (StaleInput(),))
    monkeypatch.setitem(pipeline.STAGE_RUNNERS, "ingest", runner)

    with pytest.raises(StaleArtifactError, match="upstream changed"):
        pipeline.run_stage(context, "ingest")
    assert not ran
    assert not (tmp_path / "ingest").exists()


def test_pipelined_cli_dispatches_dependency_prefix_and_worker_overrides(
    monkeypatch,
):
    from phenocycler import donor_pipeline, pipeline

    context = SimpleNamespace(
        registry=SimpleNamespace(production_status="verified")
    )
    loaded = {}
    dispatched = {}

    def fake_load(path, **overrides):
        loaded.update(overrides)
        return SimpleNamespace(marker_registry="markers")

    monkeypatch.setattr(pipeline, "load_config", fake_load)
    monkeypatch.setattr(pipeline, "load_registry", lambda _path: context.registry)
    monkeypatch.setattr(
        pipeline,
        "resolve_run_context",
        lambda _config: context,
    )
    monkeypatch.setattr(
        donor_pipeline,
        "run_pipelined_pipeline",
        lambda passed, *, stages, export: dispatched.update(
            context=passed,
            stages=stages,
            export=export,
        ),
    )

    assert (
        pipeline.main(
            [
                "run",
                "--pipelined",
                "--through",
                "geometry",
                "--jobs",
                "4",
                "--geometry-workers",
                "2",
                "--redsea-workers",
                "1",
                "--downstream-workers",
                "3",
                "--no-export",
            ]
        )
        == 0
    )
    assert loaded == {
        "n_jobs": 4,
        "pipeline_geometry_workers": 2,
        "pipeline_redsea_workers": 1,
        "pipeline_downstream_workers": 3,
    }
    assert dispatched == {
        "context": context,
        "stages": ("ingest", "geometry"),
        "export": False,
    }


def test_pair_dependent_stages_fail_before_execution():
    context = SimpleNamespace(
        registry=SimpleNamespace(
            production_status="development_no_verified_typing"
        )
    )
    with pytest.raises(ContractError, match="NO_VERIFIED_TYPING"):
        _assert_stages_available(context, ("controls",))
    with pytest.raises(ContractError, match="NO_VERIFIED_TYPING"):
        _assert_stages_available(context, ("type",), operation="QuPath export")
    _assert_stages_available(context, DEFAULT_STAGE_ORDER)


def test_no_pair_default_run_finalizes_at_expression(monkeypatch, tmp_path):
    from phenocycler import pipeline

    stage_root = tmp_path / "manifests"
    observed: dict[str, object] = {"ran": []}

    class FakeStageManifest:
        def __init__(self, stage):
            self.stage = stage

        def validate_current(self, **_kwargs):
            return None

    class FakeRunManifest:
        content_id = "upstream-run"

        def write_json(self, path):
            path.write_text("{}\n", encoding="utf-8")

    def fake_run_stage(_context, stage):
        observed["ran"].append(stage)
        path = stage_root / f"{stage}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}\n", encoding="utf-8")

    def fake_create(**kwargs):
        observed["completion_paths"] = tuple(kwargs["stage_manifest_paths"])
        return FakeRunManifest()

    monkeypatch.setattr(pipeline, "run_stage", fake_run_stage)
    monkeypatch.setattr(
        pipeline.StageManifest,
        "read_json",
        lambda path: FakeStageManifest(path.stem),
    )
    monkeypatch.setattr(pipeline, "_stage_config", lambda *_args: {})
    monkeypatch.setattr(pipeline.RunManifest, "create", fake_create)

    context = SimpleNamespace(
        run_id="upstream-only",
        run_root=tmp_path / "run",
        donors=("D1",),
        registry=SimpleNamespace(
            production_status="development_no_verified_typing",
            fingerprint="registry",
        ),
        broad_gate_policy=None,
        identity_config={},
        cohort=SimpleNamespace(content_id="cohort"),
        code="code",
        base_config=SimpleNamespace(runs_dir=tmp_path / "runs"),
        stage_manifest_path=lambda stage: stage_root / f"{stage}.json",
    )

    result = pipeline.run_pipeline(
        context,
        stages=DEFAULT_STAGE_ORDER,
        export=False,
    )

    assert result.content_id == "upstream-run"
    assert observed["ran"] == list(DEFAULT_STAGE_ORDER)
    assert [path.stem for path in observed["completion_paths"]] == list(
        DEFAULT_STAGE_ORDER
    )
    assert (context.run_root / "run.json").is_file()
    assert (context.base_config.runs_dir / "LATEST").read_text().strip() == (
        context.run_id
    )
    assert not any(
        context.stage_manifest_path(stage).exists()
        for stage in ("controls", "calibrate", "type", "states")
    )


@pytest.mark.parametrize(
    "argv",
    [
        ["run", "--only", "controls"],
        ["run", "--through", "type"],
        ["export"],
    ],
)
def test_cli_blocks_pair_work_before_resolving_cohort(monkeypatch, argv):
    from phenocycler import pipeline
    from phenocycler.marker_registry import load_registry

    config = SimpleNamespace(
        marker_registry=load_registry().fingerprint,
        qupath_manifest="must-not-be-read",
    )
    monkeypatch.setattr(pipeline, "load_config", lambda *_args, **_kwargs: config)
    monkeypatch.setattr(pipeline, "load_registry", lambda _path: load_registry())
    monkeypatch.setattr(
        pipeline,
        "resolve_run_context",
        lambda _config: pytest.fail("pair preflight must precede cohort resolution"),
    )

    with pytest.raises(ContractError, match="NO_VERIFIED_TYPING"):
        pipeline.main(argv)


def test_status_reports_donor_receipt_progress(capsys, tmp_path):
    from phenocycler.pipeline import status

    manifests = tmp_path / "manifests"
    receipt = (
        manifests
        / "donors"
        / "donor_id=1"
        / "geometry.json"
    )
    receipt.parent.mkdir(parents=True)
    receipt.write_text("{}\n")
    context = SimpleNamespace(
        donors=("1", "2"),
        run_id="run",
        registry=SimpleNamespace(
            production_status="development_no_verified_typing"
        ),
        config=SimpleNamespace(manifests_dir=manifests),
        stage_manifest_path=lambda stage: manifests / f"{stage}.json",
    )

    assert status(context) == 1
    output = capsys.readouterr().out
    assert "geometry    IN PROGRESS   1/2" in output
    assert "controls    BLOCKED NO_VERIFIED_TYPING" in output
    assert "calibrate   BLOCKED NO_VERIFIED_TYPING" in output
    assert "type        BLOCKED NO_VERIFIED_TYPING" in output
    assert "states      OPTIONAL (explicit-only)" in output
    assert "qupath      BLOCKED NO_VERIFIED_TYPING" in output
