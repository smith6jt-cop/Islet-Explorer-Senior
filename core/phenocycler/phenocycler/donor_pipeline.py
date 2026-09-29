"""A resource-aware cross-stage donor task queue.

Completion is read from the output tree: a donor whose partition has content is
done. There are no receipts, no fingerprints and no content ids -- output that
exists is output that counts.
"""

from __future__ import annotations

import multiprocessing
import os
import tempfile
from collections.abc import Mapping, Sequence
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq

from .artifacts import (
    ContractError,
    PartialArtifactError,
)
from .expression import (
    build_selected_expression,
    read_redsea_for_expression,
    read_single_partition,
    write_donor_partition,
)
from .geometry_qc import GeometryQCConfig
from .geometry_stage import run_geometry_qc_for_image
from .parallel import resolve_jobs, terminate_process_pool_workers
from .production_typing import type_production_cells as type_calibrated_cells
from .redsea import RedseaParams, process_donor
from .reference_controls import (
    build_reference_controls,
    reference_masks_to_wide,
)
from .restore_stage import RestoreStageConfig, restore_donor_from_paths
from .state_markers import annotate_proliferation

DONOR_STAGES = (
    "geometry",
    "redsea",
    "expression",
    "controls",
    "calibrate",
    "type",
    "states",
)
DONOR_DEPENDENCIES: dict[str, tuple[str, ...]] = {
    "geometry": (),
    "redsea": ("geometry",),
    "expression": ("geometry", "redsea"),
    "controls": ("expression",),
    "calibrate": ("expression", "controls"),
    "type": ("calibrate",),
    "states": ("expression",),
}
DONOR_AUDITS = {
    "expression": "expression_availability",
    "controls": "reference_control_models",
    "calibrate": "restore_pair_models",
    "states": "state_models",
}


def _atomic_write_frame(path: Path, frame: pd.DataFrame) -> Path:
    """Write the frame, replacing whatever is there.

    This used to read the existing file back, compare it row for row and raise
    StaleArtifactError when it differed, which made a rerun with changed inputs an
    error rather than an update. The write is still atomic via os.replace.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        frame.to_parquet(temporary, index=False)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def donor_audit_path(context, stage: str, donor_id: str) -> Path | None:
    name = DONOR_AUDITS.get(stage)
    if name is None:
        return None
    return (
        context.config.audit_dir
        / "donors"
        / name
        / f"donor_id={donor_id}"
        / "data.parquet"
    )


def _geometry_config(context, image) -> GeometryQCConfig:
    return GeometryQCConfig(
        pixel_size_um_x=image.pixel_size_um_x,
        pixel_size_um_y=image.pixel_size_um_y,
        min_cell_area_um2=context.config.cell_qc_min_cell_area,
        max_cell_area_median_multiple=(
            context.config.cell_qc_max_area_median_multiple
        ),
        min_cell_solidity=context.config.cell_qc_min_cell_solidity,
        min_raster_area_ratio=context.config.cell_qc_min_raster_area_ratio,
        max_raster_area_ratio=context.config.cell_qc_max_raster_area_ratio,
        no_cytoplasm_nc_ratio=context.config.cell_qc_no_cytoplasm_nc_ratio,
        duplicate_radius_um=context.config.cell_qc_duplicate_radius_um,
        duplicate_area_tolerance=context.config.cell_qc_duplicate_area_tol,
        missing_geometry_policy=(
            context.config.cell_qc_missing_geometry_policy
        ),
        duplicate_policy=context.config.cell_qc_duplicate_policy,
    )


def _run_donor_stage_body(context, stage: str, donor_id: str) -> tuple[Path, ...]:
    donor = str(donor_id)
    image_by_donor = {
        image.donor_id: image for image in context.cohort.images
    }
    if stage == "geometry":
        image = image_by_donor[donor]
        result = run_geometry_qc_for_image(
            read_single_partition(context.config.cells_dir, donor),
            image,
            config=_geometry_config(context, image),
        )
        write_donor_partition(
            result.to_frame(),
            context.config.geometry_qc_dir,
            donor_id=donor,
        )
        return ()
    if stage == "redsea":
        process_donor(
            donor,
            context.config,
            RedseaParams.from_config(context.config),
        )
        return ()
    if stage == "expression":
        image = image_by_donor[donor]
        expression, availability = build_selected_expression(
            read_single_partition(context.config.cells_dir, donor),
            read_redsea_for_expression(
                context.registry, context.config.redsea_dir, donor
            ),
            read_single_partition(context.config.geometry_qc_dir, donor),
            donor_id=donor,
            acquired_markers={
                channel.marker for channel in image.channel_map
            },
            registry=context.registry,
        )
        write_donor_partition(
            expression,
            context.config.selected_expression_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            availability,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    if stage == "controls":
        expression = read_single_partition(
            context.config.selected_expression_dir,
            donor,
        )
        result = build_reference_controls(
            expression,
            donor_id=donor,
            registry=context.registry,
            value_columns={
                marker.name: marker.name
                for marker in context.registry.markers
            },
        )
        write_donor_partition(
            reference_masks_to_wide(
                result.masks_long,
                registry=context.registry,
            ),
            context.config.reference_controls_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            result.audit,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    if stage == "calibrate":
        # The stage key remains ``calibrate`` for immutable run compatibility;
        # RESTORE is the production marker-calling implementation.
        result = restore_donor_from_paths(
            donor_id=donor,
            cells_root=context.config.cells_dir,
            redsea_root=context.config.redsea_dir,
            expression_root=context.config.selected_expression_dir,
            geometry_root=context.config.geometry_qc_dir,
            reference_controls_root=context.config.reference_controls_dir,
            registry=context.registry,
            config=RestoreStageConfig.from_pipeline(context.config),
            expression_reviews_path=context.config.restore_expression_reviews,
        )
        write_donor_partition(
            result.evidence_wide,
            context.config.marker_evidence_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            result.model_table,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    if stage == "type":
        assignments = type_calibrated_cells(
            read_single_partition(
                context.config.marker_evidence_dir,
                donor,
            ),
            marker_registry=context.registry,
            broad_gate_policy=context.broad_gate_policy,
        )
        write_donor_partition(
            assignments,
            context.config.assignments_dir,
            donor_id=donor,
        )
        return ()
    if stage == "states":
        states, models = annotate_proliferation(
            read_single_partition(
                context.config.selected_expression_dir,
                donor,
            ),
            donor_id=donor,
        )
        for row in models.itertuples(index=False):
            states[f"{row.marker}__threshold"] = row.threshold
            states[f"{row.marker}__model_status"] = row.status
            states[f"{row.marker}__method"] = row.method
        write_donor_partition(
            states,
            context.config.state_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            models,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    raise ContractError(f"unknown donor stage {stage!r}")


def _partition_has_content(root: Path, donor_id: str) -> bool:
    partition = root / f"donor_id={donor_id}"
    return partition.exists() and any(partition.iterdir())


def execute_donor_stage(context, stage: str, donor_id: str, *, force: bool = False) -> str:
    """Run one ready donor task. Present output is done; nothing is fingerprinted.

    A stage ``pipeline.stage_noop_reason`` names (REDSEA while no marker reads it)
    is logged and returns as done without running or writing anything, so the
    tasks that depend on it still become ready.
    """

    from .pipeline import _assert_stages_available, stage_noop_reason

    _assert_stages_available(
        context,
        (stage,),
        operation="direct donor-stage execution",
    )

    donor = str(donor_id)
    reason = stage_noop_reason(context.registry, stage)
    if reason is not None:
        print(f"[donor-noop] {stage}/{donor}: {reason}", flush=True)
        return donor
    output_root = context.output_root(stage)
    if _partition_has_content(output_root, donor) and not force:
        print(f"[donor-skip] {stage}/{donor}: output present (--force to rerun)", flush=True)
        return donor
    print(f"[donor-run] {stage}/{donor}", flush=True)
    _run_donor_stage_body(context, stage, donor)
    print(f"[donor-done] {stage}/{donor}", flush=True)
    return donor


@dataclass(frozen=True, order=True)
class DonorTask:
    stage: str
    donor_id: str


@dataclass(frozen=True)
class ResourceRequest:
    group: str
    cpu_slots: int = 1
    gpu_slots: int = 0


@dataclass(frozen=True)
class ResourceLimits:
    total_cpu_slots: int
    gpu_slots: int
    geometry_workers: int
    redsea_workers: int
    downstream_workers: int

    @classmethod
    def from_context(
        cls,
        context,
        *,
        task_count: int,
    ) -> ResourceLimits:
        total = resolve_jobs(context.config.n_jobs, task_count)
        redsea = min(context.config.pipeline_redsea_workers, total)
        if context.config.use_gpu:
            redsea = 1
        return cls(
            total_cpu_slots=total,
            gpu_slots=1 if context.config.use_gpu else 0,
            geometry_workers=min(
                context.config.pipeline_geometry_workers,
                total,
            ),
            redsea_workers=redsea,
            downstream_workers=min(
                context.config.pipeline_downstream_workers,
                total,
            ),
        )

    def group_limit(self, group: str) -> int:
        if group == "geometry":
            return self.geometry_workers
        if group == "redsea":
            return self.redsea_workers
        if group == "downstream":
            return self.downstream_workers
        raise ContractError(f"unknown resource group {group!r}")


@dataclass
class ResourceState:
    cpu_slots: int = 0
    gpu_slots: int = 0
    geometry_workers: int = 0
    redsea_workers: int = 0
    downstream_workers: int = 0

    def _group_value(self, group: str) -> int:
        return int(getattr(self, f"{group}_workers"))

    def can_start(
        self,
        request: ResourceRequest,
        limits: ResourceLimits,
    ) -> bool:
        return (
            self.cpu_slots + request.cpu_slots <= limits.total_cpu_slots
            and self.gpu_slots + request.gpu_slots <= limits.gpu_slots
            and self._group_value(request.group)
            < limits.group_limit(request.group)
        )

    def reserve(self, request: ResourceRequest) -> None:
        self.cpu_slots += request.cpu_slots
        self.gpu_slots += request.gpu_slots
        attribute = f"{request.group}_workers"
        setattr(self, attribute, getattr(self, attribute) + 1)

    def release(self, request: ResourceRequest) -> None:
        self.cpu_slots -= request.cpu_slots
        self.gpu_slots -= request.gpu_slots
        attribute = f"{request.group}_workers"
        setattr(self, attribute, getattr(self, attribute) - 1)


def task_resources(context, task: DonorTask) -> ResourceRequest:
    if task.stage == "geometry":
        return ResourceRequest(group="geometry")
    if task.stage == "redsea":
        return ResourceRequest(
            group="redsea",
            gpu_slots=1 if context.config.use_gpu else 0,
        )
    return ResourceRequest(group="downstream")


def task_ready(task: DonorTask, completed: set[DonorTask]) -> bool:
    return all(
        DonorTask(dependency, task.donor_id) in completed
        for dependency in DONOR_DEPENDENCIES[task.stage]
    )


def _cohort_barrier_ready(context, task: DonorTask) -> bool:
    """Donor-local dependencies are sufficient for simultaneous typing."""

    return True


def _task_priority(
    task: DonorTask,
    *,
    donor_order: Mapping[str, int],
) -> tuple[int, int, int]:
    stage_index = DONOR_STAGES.index(task.stage)
    if task.stage == "redsea":
        group_priority = 0
    elif task.stage == "geometry":
        group_priority = 1
    else:
        group_priority = 2
    downstream_priority = {
        "expression": 0,
        "controls": 1,
        "calibrate": 2,
        "type": 3,
        "states": 4,
    }
    return (
        group_priority,
        donor_order[task.donor_id],
        (
            downstream_priority[task.stage]
            if group_priority == 2
            else stage_index
        ),
    )


def select_ready_task(
    pending: set[DonorTask],
    completed: set[DonorTask],
    *,
    context,
    resources: ResourceState,
    limits: ResourceLimits,
) -> DonorTask | None:
    donor_order = {
        donor: index for index, donor in enumerate(context.donors)
    }
    candidates = sorted(
        (
            task
            for task in pending
            if task_ready(task, completed)
            and _cohort_barrier_ready(context, task)
            and resources.can_start(task_resources(context, task), limits)
        ),
        key=lambda task: _task_priority(
            task,
            donor_order=donor_order,
        ),
    )
    return candidates[0] if candidates else None


def _aggregate_donor_audit(context, stage: str) -> Path | None:
    from .pipeline import STAGE_SIDECARS

    name = DONOR_AUDITS.get(stage)
    if name is None:
        return None
    frames = []
    for donor in context.donors:
        path = donor_audit_path(context, stage, donor)
        if path is None or not path.is_file():
            raise PartialArtifactError(
                f"{stage}: donor audit is missing for {donor}"
            )
        frames.append(pd.read_parquet(path))
    combined = pd.concat(frames, ignore_index=True)
    relative_paths = STAGE_SIDECARS.get(stage, ())
    if len(relative_paths) != 1:
        raise ContractError(
            f"{stage}: expected exactly one aggregate audit sidecar"
        )
    target = context.run_root / relative_paths[0]
    return _atomic_write_frame(target, combined)


def _stage_is_complete(context, stage: str) -> bool:
    """True when every expected donor has a non-empty partition for this stage."""

    root = context.output_root(stage)
    return all(_partition_has_content(root, donor) for donor in context.donors)


def _initial_task_state(
    context,
    stages: Sequence[str],
    *,
    force: bool = False,
) -> set[DonorTask]:
    """Tasks already done, judged by output on disk rather than by a receipt.

    Output present without a receipt used to raise PartialArtifactError and stop the
    run. There are no receipts now, so partial output is simply the donors that are
    not done yet -- the queue picks up the rest.
    """

    if force:
        return set()
    completed: set[DonorTask] = set()
    for stage in stages:
        root = context.output_root(stage)
        for donor in context.donors:
            if _partition_has_content(root, donor):
                completed.add(DonorTask(stage, donor))
    return completed


def _run_tasks_serial(
    context,
    pending: set[DonorTask],
    completed: set[DonorTask],
    stages: Sequence[str],
    limits: ResourceLimits,
) -> None:
    resources = ResourceState()
    while pending:
        task = select_ready_task(
            pending,
            completed,
            context=context,
            resources=resources,
            limits=limits,
        )
        if task is None:
            blocked = sorted(
                f"{item.stage}/{item.donor_id}" for item in pending
            )
            raise PartialArtifactError(
                f"donor task queue is deadlocked: {blocked[:10]}"
            )
        request = task_resources(context, task)
        resources.reserve(request)
        try:
            execute_donor_stage(context, task.stage, task.donor_id)
        finally:
            resources.release(request)
        pending.remove(task)
        completed.add(task)


def _run_tasks_parallel(
    context,
    pending: set[DonorTask],
    completed: set[DonorTask],
    stages: Sequence[str],
    limits: ResourceLimits,
) -> None:
    resources = ResourceState()
    running = {}
    spawn_context = multiprocessing.get_context("spawn")
    executors = {
        group: ProcessPoolExecutor(
            max_workers=limits.group_limit(group),
            mp_context=spawn_context,
        )
        for group in ("geometry", "redsea", "downstream")
    }
    aborting = False
    try:
        while pending or running:
            while len(running) < limits.total_cpu_slots:
                task = select_ready_task(
                    pending,
                    completed,
                    context=context,
                    resources=resources,
                    limits=limits,
                )
                if task is None:
                    break
                request = task_resources(context, task)
                resources.reserve(request)
                pending.remove(task)
                future = executors[request.group].submit(
                    execute_donor_stage,
                    context,
                    task.stage,
                    task.donor_id,
                )
                running[future] = (task, request)
            if not running:
                blocked = sorted(
                    f"{item.stage}/{item.donor_id}" for item in pending
                )
                raise PartialArtifactError(
                    f"donor task queue is deadlocked: {blocked[:10]}"
                )
            finished, _ = wait(
                tuple(running),
                return_when=FIRST_COMPLETED,
            )
            for future in finished:
                task, request = running.pop(future)
                resources.release(request)
                try:
                    future.result()
                except Exception as exc:
                    raise RuntimeError(
                        f"donor task failed: {task.stage}/{task.donor_id}"
                    ) from exc
                completed.add(task)
    except BaseException:
        aborting = True
        for future in running:
            future.cancel()
        _terminate_executor_workers(executors)
        raise
    finally:
        if not aborting:
            for executor in executors.values():
                executor.shutdown(wait=True, cancel_futures=True)


def _terminate_executor_workers(
    executors: Mapping[str, ProcessPoolExecutor],
) -> None:
    """Stop active pool children so an interrupted queue cannot become orphaned.

    Python 3.12 does not expose ``terminate_workers`` on
    :class:`ProcessPoolExecutor`. The executor's process mapping is therefore
    captured before non-waiting shutdown, with TERM followed by a bounded KILL
    fallback. Only children owned by these three queue executors are targeted.
    """

    process_count = sum(
        len(getattr(executor, "_processes", None) or ())
        for executor in executors.values()
    )
    print(
        f"[queue] stopping {process_count} active worker processes",
        flush=True,
    )
    terminate_process_pool_workers(executors)


def run_donor_queue(
    context,
    stages: Sequence[str],
    *,
    force: bool = False,
) -> None:
    """Execute a donor-stage DAG under explicit CPU/GPU/group limits."""

    selected = tuple(stages)
    from .pipeline import _assert_stages_available

    _assert_stages_available(context, selected)
    unknown = sorted(set(selected).difference(DONOR_STAGES))
    if unknown:
        raise ContractError(f"unknown donor stages: {unknown}")
    expected_prefix = DONOR_STAGES[: len(selected)]
    if selected != expected_prefix:
        raise ContractError(
            "pipelined donor stages must be a dependency prefix"
        )
    completed = _initial_task_state(context, selected, force=force)
    all_tasks = {
        DonorTask(stage, donor)
        for stage in selected
        for donor in context.donors
    }
    pending = all_tasks.difference(completed)
    limits = ResourceLimits.from_context(
        context,
        task_count=max(1, len(pending)),
    )
    print(
        "[queue] "
        f"cpu={limits.total_cpu_slots}, gpu={limits.gpu_slots}, "
        f"geometry={limits.geometry_workers}, "
        f"redsea={limits.redsea_workers}, "
        f"downstream={limits.downstream_workers}, "
        f"pending={len(pending)}",
        flush=True,
    )
    if not pending:
        print("[queue] nothing to do: every donor partition is present", flush=True)
        return
    if limits.total_cpu_slots == 1:
        _run_tasks_serial(
            context,
            pending,
            completed,
            selected,
            limits,
        )
    else:
        _run_tasks_parallel(
            context,
            pending,
            completed,
            selected,
            limits,
        )


def run_pipelined_pipeline(
    context,
    *,
    stages: Sequence[str],
    export: bool,
    force: bool = False,
):
    """Run ingest as a barrier, then overlap ready donor stages.

    ``force`` reruns donors whose partition already has content. It was missing here until
    2026-09-16, so ``--force --pipelined`` -- the invocation CLAUDE.md documents -- accepted
    the flag and silently ignored it.
    """

    from .pipeline import (
        STAGE_ORDER,
        _assert_stages_available,
        export_qupath,
        run_pipeline,
        run_stage,
    )

    selected = tuple(stages)
    _assert_stages_available(context, selected)
    expected_prefix = STAGE_ORDER[: len(selected)]
    if selected != expected_prefix:
        raise ContractError(
            "pipelined execution requires a dependency-prefix stage selection"
        )
    context.run_root.mkdir(parents=True, exist_ok=True)
    if not selected:
        return run_pipeline(context, stages=(), export=False)
    run_stage(context, "ingest", force=force)
    donor_stages = tuple(
        stage for stage in selected if stage in DONOR_STAGES
    )
    if donor_stages:
        run_donor_queue(context, donor_stages, force=force)
    if export and "type" in selected:
        export_qupath(context)
    return run_pipeline(context, stages=(), export=False)


def donor_progress(context, stage: str) -> tuple[int, int]:
    """How many donors have output for this stage, counted from the output tree."""

    if stage not in DONOR_STAGES:
        return 0, len(context.donors)
    root = context.output_root(stage)
    count = sum(
        _partition_has_content(root, donor) for donor in context.donors
    )
    return count, len(context.donors)


__all__ = [
    "DONOR_DEPENDENCIES",
    "DONOR_STAGES",
    "DonorTask",
    "ResourceLimits",
    "ResourceRequest",
    "ResourceState",
    "donor_progress",
    "execute_donor_stage",
    "run_donor_queue",
    "run_pipelined_pipeline",
    "select_ready_task",
    "task_ready",
    "task_resources",
]
