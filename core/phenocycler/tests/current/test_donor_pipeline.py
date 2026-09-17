from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from phenocycler.artifacts import ContractError
from phenocycler.donor_pipeline import (
    DONOR_STAGES,
    DonorTask,
    ResourceLimits,
    ResourceState,
    _aggregate_donor_audit,
    _initial_task_state,
    _run_donor_stage_body,
    _terminate_executor_workers,
    execute_donor_stage,
    run_donor_queue,
    select_ready_task,
    task_ready,
    task_resources,
)
from phenocycler.expression import write_donor_partition

REPOSITORY = Path(__file__).resolve().parents[2]


def _frame(donor: str, objects=("a", "b")) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "donor_id": [donor] * len(objects),
            "object_id": list(objects),
            "value": range(len(objects)),
        }
    )


def test_task_readiness_follows_the_donor_local_dag():
    completed = {DonorTask("geometry", "1")}
    assert task_ready(DonorTask("redsea", "1"), completed)
    assert not task_ready(DonorTask("expression", "1"), completed)
    completed.add(DonorTask("redsea", "1"))
    assert task_ready(DonorTask("expression", "1"), completed)
    completed.add(DonorTask("expression", "1"))
    assert task_ready(DonorTask("controls", "1"), completed)
    assert task_ready(DonorTask("states", "1"), completed)


@pytest.mark.parametrize(
    "entrypoint",
    [execute_donor_stage],
)
def test_direct_donor_mutation_refuses_pair_stages_before_paths(entrypoint):
    context = SimpleNamespace(
        registry=SimpleNamespace(
            production_status="development_no_verified_typing"
        ),
        output_root=lambda _stage: pytest.fail("output path must not be resolved"),
    )

    with pytest.raises(ContractError, match="NO_VERIFIED_TYPING"):
        entrypoint(context, "controls", "D")


def test_resource_selection_keeps_one_gpu_redsea_job_and_feeds_geometry():
    context = SimpleNamespace(
        donors=("1", "2"),
        registry=SimpleNamespace(production_status="verified"),
        config=SimpleNamespace(
            n_jobs=4,
            use_gpu=True,
            pipeline_geometry_workers=2,
            pipeline_redsea_workers=4,
            pipeline_downstream_workers=4,
        ),
    )
    limits = ResourceLimits.from_context(context, task_count=10)
    assert limits.gpu_slots == 1
    assert limits.redsea_workers == 1
    resources = ResourceState()
    completed = {DonorTask("geometry", "1")}
    pending = {
        DonorTask("redsea", "1"),
        DonorTask("geometry", "2"),
    }

    first = select_ready_task(
        pending,
        completed,
        context=context,
        resources=resources,
        limits=limits,
    )
    assert first == DonorTask("redsea", "1")
    resources.reserve(task_resources(context, first))
    pending.remove(first)
    second = select_ready_task(
        pending,
        completed,
        context=context,
        resources=resources,
        limits=limits,
    )
    assert second == DonorTask("geometry", "2")
    assert not resources.can_start(
        task_resources(context, DonorTask("redsea", "2")),
        limits,
    )


def test_queue_abort_terminates_owned_workers_without_waiting():
    class Process:
        def __init__(self):
            self.alive = True
            self.terminated = False
            self.killed = False
            self.joins = []

        def is_alive(self):
            return self.alive

        def terminate(self):
            self.terminated = True

        def kill(self):
            self.killed = True
            self.alive = False

        def join(self, timeout):
            self.joins.append(timeout)

    class Executor:
        def __init__(self, process):
            self._processes = {1: process}
            self.shutdown_call = None

        def shutdown(self, *, wait, cancel_futures):
            self.shutdown_call = (wait, cancel_futures)

    process = Process()
    executor = Executor(process)
    _terminate_executor_workers({"geometry": executor})
    assert process.terminated
    assert process.killed
    assert process.joins == [2.0, 2.0]
    assert executor.shutdown_call == (False, True)


def test_serial_queue_obeys_dependencies_without_stage_barriers(
    monkeypatch,
):
    from phenocycler import donor_pipeline

    context = SimpleNamespace(
        donors=("1", "2"),
        registry=SimpleNamespace(production_status="verified"),
        config=SimpleNamespace(
            n_jobs=1,
            use_gpu=False,
            pipeline_geometry_workers=1,
            pipeline_redsea_workers=1,
            pipeline_downstream_workers=1,
        ),
    )
    events = []
    monkeypatch.setattr(
        donor_pipeline,
        "_initial_task_state",
        lambda *_, **__: set(),
    )
    monkeypatch.setattr(
        donor_pipeline,
        "execute_donor_stage",
        lambda _context, stage, donor: events.append((stage, donor)),
    )

    run_donor_queue(
        context,
        ("geometry", "redsea", "expression"),
    )
    positions = {
        event: index for index, event in enumerate(events)
    }
    for donor in context.donors:
        assert positions[("geometry", donor)] < positions[("redsea", donor)]
        assert positions[("redsea", donor)] < positions[("expression", donor)]
    assert positions[("redsea", "1")] < positions[("geometry", "2")]

    with pytest.raises(ContractError, match="dependency prefix"):
        run_donor_queue(context, ("geometry", "expression"))


def test_donor_audits_aggregate_in_cohort_order(tmp_path):
    context = SimpleNamespace(
        donors=("2", "1"),
        run_root=tmp_path,
        config=SimpleNamespace(audit_dir=tmp_path / "audit"),
    )
    for donor in context.donors:
        root = (
            context.config.audit_dir
            / "donors"
            / "expression_availability"
        )
        write_donor_partition(
            pd.DataFrame(
                {"donor_id": [donor], "marker": ["CD3e"]}
            ),
            root,
            donor_id=donor,
        )

    path = _aggregate_donor_audit(context, "expression")
    assert path == tmp_path / "audit" / "expression_availability.parquet"
    assert pd.read_parquet(path)["donor_id"].tolist() == ["2", "1"]
    assert _aggregate_donor_audit(context, "expression") == path


def test_all_declared_donor_stages_have_resource_requests():
    context = SimpleNamespace(
        config=SimpleNamespace(use_gpu=False)
    )
    assert {
        task_resources(context, DonorTask(stage, "1")).group
        for stage in DONOR_STAGES
    } == {"geometry", "redsea", "downstream"}


def test_simultaneous_donor_typing_needs_only_donor_calibration():
    context = SimpleNamespace(
        donors=("1",),
        config=SimpleNamespace(use_gpu=False),
    )
    task = DonorTask("type", "1")
    pending = {task}
    completed = {DonorTask("calibrate", "1")}
    resources = ResourceState()
    limits = ResourceLimits(
        total_cpu_slots=1,
        gpu_slots=0,
        geometry_workers=1,
        redsea_workers=1,
        downstream_workers=1,
    )

    assert select_ready_task(
        pending,
        completed,
        context=context,
        resources=resources,
        limits=limits,
    ) == task


def test_donor_type_body_uses_simultaneous_interface(monkeypatch, tmp_path):
    from phenocycler import donor_pipeline

    evidence = pd.DataFrame({"donor_id": ["1"], "object_id": ["a"]})
    assignments = evidence.assign(cell_type="Immune")
    observed = {}
    monkeypatch.setattr(
        donor_pipeline,
        "read_single_partition",
        lambda *_args, **_kwargs: evidence,
    )

    def fake_type(frame, **kwargs):
        observed.update(frame=frame, **kwargs)
        return assignments

    monkeypatch.setattr(donor_pipeline, "type_calibrated_cells", fake_type)
    monkeypatch.setattr(
        donor_pipeline,
        "write_donor_partition",
        lambda *_args, **_kwargs: tmp_path / "data.parquet",
    )
    context = SimpleNamespace(
        config=SimpleNamespace(
            marker_evidence_dir=tmp_path / "evidence",
            assignments_dir=tmp_path / "assignments",
        ),
        cohort=SimpleNamespace(images=()),
        registry="markers",
        broad_gate_policy="gates",
    )

    assert _run_donor_stage_body(context, "type", "1") == ()
    assert observed["frame"] is evidence
    assert observed["marker_registry"] == "markers"
    assert observed["broad_gate_policy"] == "gates"


def test_state_donor_body_writes_output_and_donor_audit(tmp_path):
    config = SimpleNamespace(
        selected_expression_dir=tmp_path / "expression",
        state_dir=tmp_path / "states",
        audit_dir=tmp_path / "audit",
    )
    context = SimpleNamespace(config=config, cohort=SimpleNamespace(images=()))
    expression = pd.DataFrame(
        {
            "donor_id": ["1", "1"],
            "object_id": ["a", "b"],
            "qc_analysis_eligible": [True, False],
            "qc_estimation_eligible": [True, False],
            "Ki67": [1.0, 2.0],
            "PCNA": [3.0, 4.0],
        }
    )
    write_donor_partition(
        expression,
        config.selected_expression_dir,
        donor_id="1",
    )

    sidecars = _run_donor_stage_body(context, "states", "1")
    assert len(sidecars) == 1
    assert sidecars[0].is_file()
    assert (
        config.state_dir / "donor_id=1" / "data.parquet"
    ).is_file()


def test_donor_output_without_a_receipt_counts_as_done(tmp_path):
    """Output present on disk IS the completion record.

    This used to raise PartialArtifactError: a donor whose stage had finished but
    whose receipt was missing -- a kill, a crash, or an edit to unrelated package
    source -- stopped the whole queue, and the only way back was a separate
    `recover` command. Now the partition speaks for itself and the queue skips it.
    """

    context = SimpleNamespace(
        donors=("1", "2"),
        config=SimpleNamespace(
            manifests_dir=tmp_path / "manifests",
            audit_dir=tmp_path / "audit",
        ),
        output_root=lambda stage: tmp_path / stage,
    )
    write_donor_partition(
        _frame("1"),
        context.output_root("geometry"),
        donor_id="1",
    )

    completed = _initial_task_state(context, ("geometry",))
    assert completed == {DonorTask("geometry", "1")}
