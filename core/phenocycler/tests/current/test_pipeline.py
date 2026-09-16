from __future__ import annotations

from types import SimpleNamespace

import pytest

from phenocycler.artifacts import ContractError
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


def test_export_manifest_presence_is_not_an_identity_check(tmp_path):
    """A malformed export manifest reads as present, not as a contract violation.

    It used to recompute sha256_json over the payload, compare it to a recorded
    content_id, and additionally refuse a manifest whose run_id belonged to another
    run. Both checks are gone: the file is either there or it is not.
    """

    manifests = tmp_path / "manifests"
    manifests.mkdir()
    (manifests / "qupath_export.json").write_text("{}\n")
    context = SimpleNamespace(config=SimpleNamespace(manifests_dir=manifests))

    assert _validate_export(context) is True


def test_a_changed_input_does_not_block_the_runner(monkeypatch, tmp_path):
    """Inputs are no longer fingerprinted, so a changed one cannot refuse to run.

    This is the failure that motivated the removal: a QuPath export renamed on disk,
    or an edit to unrelated package source, raised StaleArtifactError and stopped a
    stage whose data was fine.
    """

    from phenocycler import pipeline

    context = SimpleNamespace(
        registry=SimpleNamespace(production_status="verified"),
        output_root=lambda stage: tmp_path / stage,
    )
    ran = False

    def runner(_context):
        nonlocal ran
        ran = True

    monkeypatch.setitem(pipeline.STAGE_RUNNERS, "ingest", runner)
    pipeline.run_stage(context, "ingest")

    assert ran


def test_present_output_is_skipped_unless_forced(monkeypatch, tmp_path):
    from phenocycler import pipeline

    context = SimpleNamespace(
        registry=SimpleNamespace(production_status="verified"),
        output_root=lambda stage: tmp_path / stage,
    )
    (tmp_path / "ingest").mkdir()
    (tmp_path / "ingest" / "data.parquet").write_bytes(b"x")
    runs = []
    monkeypatch.setitem(
        pipeline.STAGE_RUNNERS, "ingest", lambda _c: runs.append(1)
    )

    pipeline.run_stage(context, "ingest")
    assert runs == []

    pipeline.run_stage(context, "ingest", force=True)
    assert runs == [1]


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


def test_default_run_writes_no_run_manifest_and_no_latest_pointer(
    monkeypatch, tmp_path
):
    """The run completes at expression and leaves no run identity behind."""

    from phenocycler import pipeline

    ran: list[str] = []
    monkeypatch.setattr(
        pipeline,
        "run_stage",
        lambda _context, stage, force=False: ran.append(stage),
    )
    context = SimpleNamespace(
        run_root=tmp_path / "run",
        donors=("D1",),
        registry=SimpleNamespace(
            production_status="development_no_verified_typing",
            fingerprint="registry",
        ),
    )

    pipeline.run_pipeline(context, export=False)

    assert ran == list(pipeline.DEFAULT_STAGE_ORDER)
    assert ran[-1] == "expression"
    assert not (tmp_path / "run" / "run.json").exists()
    assert not (tmp_path / "runs" / "LATEST").exists()


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


def test_status_reports_partitions_on_disk_and_never_gates(capsys, tmp_path):
    """Status counts donor partitions and always returns 0.

    It used to return 1 for anything short of a validated manifest, which made a
    perfectly good half-finished cohort look like a failure to every caller.
    """

    from phenocycler.pipeline import status

    geometry = tmp_path / "geometry" / "donor_id=1"
    geometry.mkdir(parents=True)
    (geometry / "data.parquet").write_bytes(b"x")
    context = SimpleNamespace(
        donors=("1", "2"),
        registry=SimpleNamespace(
            production_status="development_no_verified_typing"
        ),
        config=SimpleNamespace(data_dir=tmp_path),
        output_root=lambda stage: tmp_path / stage,
    )

    assert status(context) == 0
    output = capsys.readouterr().out
    assert "geometry    PARTIAL" in output
    assert "missing=['2']" in output
    assert "ingest      ABSENT" in output
    assert "controls    BLOCKED NO_VERIFIED_TYPING" in output
    assert "states      OPTIONAL (explicit-only)" in output
    assert "qupath      BLOCKED NO_VERIFIED_TYPING" in output


