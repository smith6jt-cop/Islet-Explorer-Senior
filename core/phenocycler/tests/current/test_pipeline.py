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
        lambda passed, *, stages, export, force=False: dispatched.update(
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
            production_status="development_no_verified_typing",
            redsea_markers=("CD3e",),  # a registry that still routes a marker to REDSEA
        ),
        config=SimpleNamespace(data_dir=tmp_path),
        output_root=lambda stage: tmp_path / stage,
    )

    assert status(context) == 0
    output = capsys.readouterr().out
    assert "geometry    PARTIAL" in output
    assert "missing=['2']" in output
    assert "ingest      ABSENT" in output
    assert "redsea      ABSENT" in output
    assert "controls    BLOCKED NO_VERIFIED_TYPING" in output
    assert "states      OPTIONAL (explicit-only)" in output
    assert "qupath      BLOCKED NO_VERIFIED_TYPING" in output


def test_status_reports_redsea_as_a_noop_not_as_unfinished(capsys, tmp_path):
    """Donors a stopped REDSEA run reached are on disk, the rest are not; with no marker
    reading REDSEA that is not a stage waiting to be finished."""

    from phenocycler.marker_registry import load_registry
    from phenocycler.pipeline import status

    redsea = tmp_path / "redsea" / "donor_id=1"
    redsea.mkdir(parents=True)
    (redsea / "data_0.parquet").write_bytes(b"x")
    context = SimpleNamespace(
        donors=("1", "2"),
        registry=load_registry(),
        config=SimpleNamespace(data_dir=tmp_path),
        output_root=lambda stage: tmp_path / stage,
    )

    assert status(context) == 0
    output = capsys.readouterr().out
    assert "redsea      NO-OP   no registry marker reads REDSEA" in output
    assert "1 donor partition(s) on disk are not read" in output
    assert "redsea      PARTIAL" not in output


def _registry_with(**policies):
    import json

    from phenocycler.marker_registry import DEFAULT_REGISTRY_PATH, registry_from_dict

    payload = json.loads(DEFAULT_REGISTRY_PATH.read_text(encoding="utf-8"))
    for row in payload["markers"]:
        if row["name"] in policies:
            row["spillover_policy"] = policies[row["name"]]
    return registry_from_dict(payload)


def test_redsea_is_a_logged_noop_when_no_marker_reads_it(monkeypatch, capsys, tmp_path):
    """The stage keeps its place in DEFAULT_STAGE_ORDER and does nothing, even when forced."""

    from phenocycler import pipeline
    from phenocycler.marker_registry import load_registry

    registry = load_registry()
    assert pipeline.stage_noop_reason(registry, "redsea") is not None
    assert all(
        pipeline.stage_noop_reason(registry, stage) is None
        for stage in pipeline.STAGE_ORDER
        if stage != "redsea"
    )
    context = SimpleNamespace(
        registry=registry,
        output_root=lambda stage: tmp_path / stage,
    )
    monkeypatch.setitem(
        pipeline.STAGE_RUNNERS,
        "redsea",
        lambda _context: pytest.fail("REDSEA ran with no marker reading it"),
    )

    pipeline.run_stage(context, "redsea")
    pipeline.run_stage(context, "redsea", force=True)

    output = capsys.readouterr().out
    assert output.count("[no-op] redsea: no registry marker reads REDSEA") == 2
    assert "[run] redsea" not in output
    assert not (tmp_path / "redsea").exists()


def test_redsea_runs_again_once_a_marker_is_routed_to_it(monkeypatch, tmp_path):
    from phenocycler import pipeline

    registry = _registry_with(CD3e="redsea_corrected")
    assert pipeline.stage_noop_reason(registry, "redsea") is None
    context = SimpleNamespace(
        registry=registry,
        output_root=lambda stage: tmp_path / stage,
    )
    runs = []
    monkeypatch.setitem(pipeline.STAGE_RUNNERS, "redsea", lambda _c: runs.append(1))

    pipeline.run_stage(context, "redsea")
    assert runs == [1]


def _expression_inputs(tmp_path, donors=("1", "2")):
    """cells and geometry_qc partitions as ingest and geometry write them, and no redsea/."""

    import pandas as pd

    from phenocycler.expression import write_donor_partition
    from phenocycler.marker_registry import load_registry

    config = SimpleNamespace(
        data_dir=tmp_path,
        cells_dir=tmp_path / "cells",
        geometry_qc_dir=tmp_path / "geometry_qc",
        redsea_dir=tmp_path / "redsea",
        selected_expression_dir=tmp_path / "expression",
        audit_dir=tmp_path / "audit",
    )
    roots = {
        "ingest": config.cells_dir,
        "geometry": config.geometry_qc_dir,
        "redsea": config.redsea_dir,
        "expression": config.selected_expression_dir,
    }
    for offset, donor in enumerate(donors):
        write_donor_partition(
            pd.DataFrame(
                {
                    "object_id": ["b", "a"],
                    "image": [f"img{donor}"] * 2,
                    "X_centroid": [2.0, 1.0],
                    "Y_centroid": [4.0, 3.0],
                    "CD3e": [20.0 + offset, 10.0 + offset],
                    "Nucleus__CD3e": [-1.0, -1.0],
                    "Histone_H3_p_Ser28": [-1.0, -1.0],
                    "Nucleus__Histone_H3_p_Ser28": [7.0 + offset, 8.0 + offset],
                }
            ),
            config.cells_dir,
            donor_id=donor,
            filename="data_0.parquet",
        )
        write_donor_partition(
            pd.DataFrame(
                {
                    "object_id": ["a", "b"],
                    "analysis_eligible": [True, True],
                    "estimation_eligible": [True, False],
                }
            ),
            config.geometry_qc_dir,
            donor_id=donor,
        )
    channel_map = tuple(
        SimpleNamespace(marker=marker) for marker in ("CD3e", "Histone_H3_(p_Ser28)")
    )
    return SimpleNamespace(
        config=config,
        cohort=SimpleNamespace(
            images=tuple(
                SimpleNamespace(donor_id=donor, channel_map=channel_map)
                for donor in donors
            )
        ),
        donors=tuple(donors),
        registry=load_registry(),
        run_root=tmp_path,
        output_root=lambda stage: roots[stage],
    )


def test_expression_runs_through_the_stage_order_without_redsea_partitions(capsys, tmp_path):
    """`run --through expression` on donors REDSEA never reached: REDSEA is a no-op and
    expression reads the raw QuPath means from cells."""

    import pandas as pd

    from phenocycler import pipeline
    from phenocycler.expression import read_single_partition

    context = _expression_inputs(tmp_path)
    pipeline.run_pipeline(context, stages=("redsea", "expression"), export=False)

    assert not (tmp_path / "redsea").exists()
    output = capsys.readouterr().out
    assert "[no-op] redsea" in output and "[done] expression" in output
    for offset, donor in enumerate(context.donors):
        expression = read_single_partition(tmp_path / "expression", donor)
        assert expression["object_id"].tolist() == ["b", "a"]
        assert expression["CD3e"].tolist() == [20.0 + offset, 10.0 + offset]
        assert expression["Histone_H3_(p_Ser28)"].tolist() == [7.0 + offset, 8.0 + offset]
        assert expression["INS"].isna().all()  # not in this panel
    audit = pd.read_parquet(tmp_path / "audit" / "expression_availability.parquet")
    available = audit[audit["available"]]
    assert sorted(set(zip(available["marker"], available["source"], strict=True))) == [
        ("CD3e", "qupath:CD3e"),
        ("Histone_H3_(p_Ser28)", "qupath:Nucleus__Histone_H3_p_Ser28"),
    ]
    assert len(available) == 2 * len(context.donors)




def test_force_reaches_the_pipelined_path(monkeypatch, tmp_path):
    """--force was accepted and silently ignored by --pipelined until 2026-09-16.

    CLAUDE.md documents `run --config config.ini --pipelined --through expression` as THE
    invocation, and run_pipelined_pipeline took no force parameter, so the flag parsed, did
    nothing, and reported nothing. A flag that is accepted and ignored is worse than one that
    is rejected.
    """
    from phenocycler import donor_pipeline

    context = SimpleNamespace(
        donors=("1", "2"),
        registry=SimpleNamespace(production_status="verified"),
        output_root=lambda stage: tmp_path / stage,
    )
    for donor in context.donors:
        part = tmp_path / "geometry" / f"donor_id={donor}"
        part.mkdir(parents=True)
        (part / "data.parquet").write_bytes(b"x")

    # Without force, present partitions are complete and nothing is pending.
    assert donor_pipeline._initial_task_state(context, ("geometry",)) == {
        donor_pipeline.DonorTask("geometry", d) for d in context.donors
    }
    # With force, every donor is pending again.
    assert donor_pipeline._initial_task_state(context, ("geometry",), force=True) == set()
