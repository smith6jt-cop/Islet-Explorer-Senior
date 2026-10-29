#!/usr/bin/env python3
"""QuPath → ingest → geometry QC → REDSEA → expression workflow.

NOT CONTENT-ADDRESSED ANY MORE.  The content-addressed run framework was REMOVED on 2026-09-16.
There is no ``run_id``, no stage manifest, no donor receipt and no fingerprint validation.  Each
stage writes to a stable root (``data/cells``, ``data/geometry_qc``, ``data/redsea``,
``data/expression``) and a donor is done when its partition has content; a stage whose output is
already present is skipped and ``--force`` reruns it (see :func:`run_stage`).  The QuPath cohort
manifest still records a ``content_id`` and per-file sha256, but nothing validates them -- it is a
description, so it can be hand-edited when a source file moves.

Ordinary execution stops after ``expression``.  ``controls``, ``calibrate``, ``type`` and the
QuPath assignment export are blocked by ``NO_VERIFIED_TYPING``, and that latch has not moved.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import pandas as pd

from .artifacts import (
    ContractError,
    QuPathCohortManifest,
    sha256_file,
)
from .cells_parquet import build_cells_parquet
from .cohort import ensure_eligible_donors
from .config import PipelineConfig, load_config
from .expression import (
    METHOD_VERSION as EXPRESSION_VERSION,
)
from .expression import (
    build_selected_expression,
    read_single_partition,
    write_donor_partition,
)
from .geometry_qc import GeometryQCConfig
from .geometry_stage import (
    METHOD_VERSION as GEOMETRY_VERSION,
)
from .geometry_stage import (
    run_geometry_qc_for_image,
)
from .marker_registry import (
    NO_VERIFIED_TYPING_STATUS,
    MarkerRegistry,
    load_registry,
)
from .production_typing import (
    METHOD_VERSION as PRODUCTION_TYPING_VERSION,
)
from .production_typing import (
    type_production_cells as type_calibrated_cells,
)
from .qupath_export import export_assignment_frame
from .redsea import RedseaParams, run_redsea
from .reference_controls import (
    METHOD_VERSION as REFERENCE_VERSION,
)
from .reference_controls import (
    ReferenceControlConfig,
    build_reference_controls,
    reference_masks_to_wide,
)
from .restore_stage import (
    METHOD_VERSION as RESTORE_STAGE_VERSION,
)
from .restore_stage import RestoreStageConfig, restore_donor_from_paths
from .simultaneous_typing import (
    BroadGatePolicy,
    load_broad_gate_policy,
)
from .state_markers import annotate_proliferation

INGEST_VERSION = "qupath-ingest-union-v3-qupath-norm-inputs"
REDSEA_VERSION = "redsea-compartment-only-v1"
STATE_VERSION = "orthogonal-state-v1"
RUN_SCHEMA_VERSION = 1

STAGE_ORDER = (
    "ingest",
    "geometry",
    "redsea",
    "expression",
    "controls",
    "calibrate",
    "type",
    "states",
)
DEFAULT_STAGE_ORDER = STAGE_ORDER[: STAGE_ORDER.index("expression") + 1]
TYPING_DEPENDENT_STAGES = frozenset({"controls", "calibrate", "type"})
NO_VERIFIED_TYPING_REASON = "NO_VERIFIED_TYPING"
STAGE_METHODS = {
    "ingest": INGEST_VERSION,
    "geometry": GEOMETRY_VERSION,
    "redsea": REDSEA_VERSION,
    "expression": EXPRESSION_VERSION,
    "controls": REFERENCE_VERSION,
    # The historical stage key is retained so completed upstream runs remain
    # adoptable. Its scientific producer is RESTORE, not donor calibration.
    "calibrate": RESTORE_STAGE_VERSION,
    "type": PRODUCTION_TYPING_VERSION,
    "states": STATE_VERSION,
}
STAGE_OUTPUT_ATTRIBUTE = {
    "ingest": "cells_dir",
    "geometry": "geometry_qc_dir",
    "redsea": "redsea_dir",
    "expression": "selected_expression_dir",
    "controls": "reference_controls_dir",
    "calibrate": "marker_evidence_dir",
    "type": "assignments_dir",
    "states": "state_dir",
}
STAGE_DEPENDENCIES = {
    "ingest": (),
    "geometry": ("ingest",),
    "redsea": ("ingest", "geometry"),
    "expression": ("ingest", "geometry", "redsea"),
    "controls": ("expression",),
    "calibrate": ("expression", "controls"),
    "type": ("calibrate",),
    "states": ("expression",),
}
STAGE_SIDECARS = {
    "ingest": ("ingest_rejects.parquet", "panel_availability.json"),
    "expression": ("audit/expression_availability.parquet",),
    "controls": ("audit/reference_control_models.parquet",),
    "calibrate": ("audit/restore_pair_models.parquet",),
    "states": ("audit/state_models.parquet",),
}
STAGE_SOURCE_PATHS = {
    "ingest": (
        "multiplex5000/cells_parquet.py",
        "multiplex5000/artifacts.py",
        "multiplex5000/cohort.py",
        "multiplex5000/marker_registry.py",
        "multiplex5000/marker_registry.json",
    ),
    "geometry": (
        "multiplex5000/donor_pipeline.py",
        "multiplex5000/geometry_stage.py",
        "multiplex5000/geometry_qc.py",
        "multiplex5000/redsea.py",
        "multiplex5000/artifacts.py",
    ),
    "redsea": (
        "multiplex5000/donor_pipeline.py",
        "multiplex5000/redsea.py",
        "multiplex5000/gpu.py",
        "multiplex5000/parallel.py",
        "multiplex5000/marker_taxonomy.py",
        "multiplex5000/marker_registry.py",
    ),
    "expression": (
        "multiplex5000/donor_pipeline.py",
        "multiplex5000/expression.py",
        "multiplex5000/marker_registry.py",
    ),
    "controls": (
        "multiplex5000/donor_pipeline.py",
        "multiplex5000/reference_controls.py",
        "multiplex5000/marker_registry.py",
    ),
    "calibrate": (
        "multiplex5000/donor_pipeline.py",
        "multiplex5000/reference_controls.py",
        "multiplex5000/restore_seeded.py",
        "multiplex5000/restore_stage.py",
        "multiplex5000/restore_validation.py",
        "multiplex5000/marker_registry.py",
    ),
    "type": (
        "multiplex5000/donor_pipeline.py",
        "multiplex5000/production_typing.py",
        "multiplex5000/simultaneous_typing.py",
        "multiplex5000/broad_gates.json",
        "multiplex5000/marker_registry.py",
    ),
    "states": (
        "multiplex5000/donor_pipeline.py",
        "multiplex5000/state_markers.py",
    ),
}


@dataclass(frozen=True)
class RunContext:
    base_config: PipelineConfig
    config: PipelineConfig
    cohort: QuPathCohortManifest
    registry: MarkerRegistry
    broad_gate_policy: BroadGatePolicy | None
    donors: tuple[str, ...]
    repository: Path
    identity_config: dict

    @property
    def run_root(self) -> Path:
        return self.config.data_dir

    def output_root(self, stage: str) -> Path:
        return Path(getattr(self.config, STAGE_OUTPUT_ATTRIBUTE[stage]))


def _identity_config(
    cfg: PipelineConfig,
    *,
    registry: MarkerRegistry | None = None,
    broad_gate_policy: BroadGatePolicy | None = None,
) -> dict:
    """Configuration that changes scientific content, excluding output location."""

    active_registry = (
        load_registry(cfg.marker_registry) if registry is None else registry
    )
    verified_pairs = (
        active_registry.production_status != NO_VERIFIED_TYPING_STATUS
    )
    gate_policy = broad_gate_policy
    if verified_pairs and gate_policy is None:
        gate_policy = load_broad_gate_policy(
            cfg.broad_gate_rules,
            registry=active_registry,
        )
    verified_pair_config = None
    expression_reviews_content_id = None
    if verified_pairs:
        verified_pair_config = RestoreStageConfig.from_pipeline(cfg)
        if cfg.restore_expression_reviews is not None:
            expression_reviews_content_id = sha256_file(
                cfg.restore_expression_reviews
            )
    return {
        "run_schema_version": RUN_SCHEMA_VERSION,
        "typing_model": (
            {"mode": "restore_simultaneous"}
            if verified_pairs
            else {
                "mode": "unavailable",
                "reason": NO_VERIFIED_TYPING_REASON,
            }
        ),
        "simultaneous_broad_typing": (
            {
                "method_version": PRODUCTION_TYPING_VERSION,
                "rules_version": gate_policy.rules_version,
                "rules_fingerprint": gate_policy.fingerprint,
                "gates": {
                    compartment: list(gate_policy.gates[compartment])
                    for compartment in gate_policy.compartments
                },
            }
            if gate_policy is not None
            else {
                "mode": "blocked",
                "reason": NO_VERIFIED_TYPING_REASON,
            }
        ),
        "cells_min_cell_area": cfg.cells_min_cell_area,
        "geometry_qc": {
            "min_cell_area_um2": cfg.cell_qc_min_cell_area,
            "max_cell_area_median_multiple": cfg.cell_qc_max_area_median_multiple,
            "min_cell_solidity": cfg.cell_qc_min_cell_solidity,
            "min_raster_area_ratio": cfg.cell_qc_min_raster_area_ratio,
            "max_raster_area_ratio": cfg.cell_qc_max_raster_area_ratio,
            "no_cytoplasm_nc_ratio": cfg.cell_qc_no_cytoplasm_nc_ratio,
            "duplicate_radius_um": cfg.cell_qc_duplicate_radius_um,
            "duplicate_area_tolerance": cfg.cell_qc_duplicate_area_tol,
            "missing_geometry_policy": cfg.cell_qc_missing_geometry_policy,
            "duplicate_policy": cfg.cell_qc_duplicate_policy,
        },
        "redsea": {
            "downsample": cfg.redsea_downsample,
            "edge_radius": cfg.redsea_edge_radius,
            "comp_mode": cfg.redsea_comp_mode,
            "alpha": cfg.redsea_alpha,
            "gap_bridge": cfg.redsea_gap_bridge,
            # only when on, so stages built before the despeckle keep their identity
            **({"despeckle": cfg.redsea_despeckle} if cfg.redsea_despeckle > 1 else {}),
            "norm_form": cfg.redsea_norm_form,
            "exclude_no_spillover": cfg.redsea_exclude_no_spillover,
            "output_schema": "Nucleus/Cytoplasm/Membrane/Cell",
        },
        "reference_controls": (
            asdict(ReferenceControlConfig())
            if verified_pairs
            else {"mode": "blocked", "reason": NO_VERIFIED_TYPING_REASON}
        ),
        "phenotyping": (
            {
                "method": "RESTORE",
                "method_version": RESTORE_STAGE_VERSION,
                **asdict(verified_pair_config),
                "expression_reviews_content_id": expression_reviews_content_id,
            }
            if verified_pair_config is not None
            else {"mode": "blocked", "reason": NO_VERIFIED_TYPING_REASON}
        ),
    }


def resolve_run_context(cfg: PipelineConfig) -> RunContext:
    """Resolve the cohort and registry a run needs.

    Output goes to the stable stage roots under ``data/`` -- ``data/cells``,
    ``data/redsea`` and so on. There is no content-addressed run directory and no
    run id: a rerun overwrites in place, and two invocations that differ only in
    code write to the same tree. The immutable-run framework this replaced keyed
    the output directory on a hash of the cohort, the config AND the package
    source, so an unrelated edit orphaned finished stages and the only way to
    reuse them was a 688-line resume module that re-validated them first.
    """

    cohort = QuPathCohortManifest.read_json(cfg.qupath_manifest)
    donors = tuple(
        ensure_eligible_donors(
            cohort.expected_donors,
            context="Phenocycler run",
        )
    )
    if donors != cohort.expected_donors:
        raise ContractError(
            "cohort donor order/set changed after exclusion validation"
        )
    registry = load_registry(cfg.marker_registry)
    broad_gate_policy = (
        None
        if registry.production_status == NO_VERIFIED_TYPING_STATUS
        else load_broad_gate_policy(
            cfg.broad_gate_rules,
            registry=registry,
        )
    )
    repository = Path(__file__).resolve().parents[1]
    identity = _identity_config(
        cfg,
        registry=registry,
        broad_gate_policy=broad_gate_policy,
    )
    return RunContext(
        base_config=cfg,
        config=cfg,
        cohort=cohort,
        registry=registry,
        broad_gate_policy=broad_gate_policy,
        donors=donors,
        repository=repository,
        identity_config=identity,
    )


def _write_audit(context: RunContext, name: str, frame: pd.DataFrame) -> Path:
    context.config.audit_dir.mkdir(parents=True, exist_ok=True)
    path = context.config.audit_dir / f"{name}.parquet"
    if path.exists():
        raise FileExistsError(f"immutable audit artifact exists: {path}")
    frame.to_parquet(path, index=False)
    return path


def _run_ingest(context: RunContext) -> None:
    build_cells_parquet(context.config)


def _run_geometry(context: RunContext) -> None:
    image_by_donor = {image.donor_id: image for image in context.cohort.images}
    for donor in context.donors:
        cells = read_single_partition(context.config.cells_dir, donor)
        image = image_by_donor[donor]
        config = GeometryQCConfig(
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
            missing_geometry_policy=context.config.cell_qc_missing_geometry_policy,
            duplicate_policy=context.config.cell_qc_duplicate_policy,
        )
        result = run_geometry_qc_for_image(cells, image, config=config)
        write_donor_partition(
            result.to_frame(),
            context.config.geometry_qc_dir,
            donor_id=donor,
        )


def _run_redsea(context: RunContext) -> None:
    params = RedseaParams.from_config(context.config)
    run_redsea(
        context.config,
        list(context.donors),
        params,
        n_jobs=context.config.n_jobs,
    )


def _run_expression(context: RunContext) -> None:
    image_by_donor = {image.donor_id: image for image in context.cohort.images}
    availability = []
    for donor in context.donors:
        expression, donor_availability = build_selected_expression(
            read_single_partition(context.config.cells_dir, donor),
            read_single_partition(context.config.redsea_dir, donor),
            read_single_partition(context.config.geometry_qc_dir, donor),
            donor_id=donor,
            acquired_markers={
                channel.marker for channel in image_by_donor[donor].channel_map
            },
            registry=context.registry,
        )
        write_donor_partition(
            expression,
            context.config.selected_expression_dir,
            donor_id=donor,
        )
        availability.append(donor_availability)
    _write_audit(
        context,
        "expression_availability",
        pd.concat(availability, ignore_index=True),
    )


def _run_controls(context: RunContext) -> None:
    audits = []
    for donor in context.donors:
        expression = read_single_partition(
            context.config.selected_expression_dir, donor
        )
        result = build_reference_controls(
            expression,
            donor_id=donor,
            registry=context.registry,
            value_columns={
                marker.name: marker.name for marker in context.registry.markers
            },
        )
        write_donor_partition(
            reference_masks_to_wide(
                result.masks_long, registry=context.registry
            ),
            context.config.reference_controls_dir,
            donor_id=donor,
        )
        audits.append(result.audit)
    _write_audit(
        context,
        "reference_control_models",
        pd.concat(audits, ignore_index=True),
    )


def _run_calibrate(context: RunContext) -> None:
    """Run RESTORE while retaining the historical immutable stage key."""

    models = []
    restore_config = RestoreStageConfig.from_pipeline(context.config)
    for donor in context.donors:
        result = restore_donor_from_paths(
            donor_id=donor,
            cells_root=context.config.cells_dir,
            redsea_root=context.config.redsea_dir,
            expression_root=context.config.selected_expression_dir,
            geometry_root=context.config.geometry_qc_dir,
            reference_controls_root=context.config.reference_controls_dir,
            registry=context.registry,
            config=restore_config,
            expression_reviews_path=context.config.restore_expression_reviews,
        )
        write_donor_partition(
            result.evidence_wide,
            context.config.marker_evidence_dir,
            donor_id=donor,
        )
        models.append(result.model_table)
    _write_audit(
        context,
        "restore_pair_models",
        pd.concat(models, ignore_index=True),
    )


def _run_type(context: RunContext) -> None:
    for donor in context.donors:
        assignments = type_calibrated_cells(
            read_single_partition(context.config.marker_evidence_dir, donor),
            marker_registry=context.registry,
            broad_gate_policy=context.broad_gate_policy,
        )
        write_donor_partition(
            assignments,
            context.config.assignments_dir,
            donor_id=donor,
        )


def _run_states(context: RunContext) -> None:
    audits = []
    for donor in context.donors:
        states, models = annotate_proliferation(
            read_single_partition(
                context.config.selected_expression_dir, donor
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
        audits.append(models)
    _write_audit(
        context,
        "state_models",
        pd.concat(audits, ignore_index=True),
    )


STAGE_RUNNERS: dict[str, Callable[[RunContext], None]] = {
    "ingest": _run_ingest,
    "geometry": _run_geometry,
    "redsea": _run_redsea,
    "expression": _run_expression,
    "controls": _run_controls,
    "calibrate": _run_calibrate,
    "type": _run_type,
    "states": _run_states,
}


def run_stage(context: RunContext, stage: str, *, force: bool = False) -> None:
    """Run one stage into its stable output root.

    Present, non-empty output is treated as done and skipped -- `--force` reruns it.
    Nothing here validates a fingerprint or a manifest, so a renamed input, a touched
    file or an edit to unrelated package source can no longer refuse to run.
    """

    _assert_stages_available(context, (stage,))
    output = context.output_root(stage)
    if output.exists() and any(output.iterdir()):
        if not force:
            print(f"[skip] {stage}: output present at {output} (--force to rerun)")
            return
        print(f"[force] {stage}: rerunning over existing output at {output}")
    print(f"[run] {stage}")
    STAGE_RUNNERS[stage](context)
    print(f"[done] {stage}: {output}")


def _export_manifest_path(context: RunContext) -> Path:
    return context.config.manifests_dir / "qupath_export.json"


def _validate_export(context: RunContext) -> bool:
    """True when the QuPath export exists. Presence only -- no identity check."""

    return _export_manifest_path(context).exists()


def export_qupath(context: RunContext) -> None:
    """Write and fingerprint one uncertainty-preserving CSV per donor."""

    _assert_stages_available(context, ("type",), operation="QuPath export")

    if _validate_export(context):
        print("[current] qupath export")
        return
    type_output = context.output_root("type")
    if not type_output.is_dir() or not any(type_output.iterdir()):
        raise ContractError("QuPath export requires the type stage")
    paths = []
    for donor in context.donors:
        paths.extend(
            export_assignment_frame(
                read_single_partition(context.config.assignments_dir, donor),
                read_single_partition(context.config.cells_dir, donor).loc[
                    :, ["object_id", "image"]
                ],
                context.config.qupath_class_dir,
            )
        )
    payload = {
        "expected_donors": list(context.donors),
        "files": [str(path) for path in paths],
    }
    manifest_path = _export_manifest_path(context)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"[done] qupath export: {len(paths)} files")


def run_pipeline(
    context: RunContext,
    *,
    stages: Sequence[str] = DEFAULT_STAGE_ORDER,
    export: bool = True,
    force: bool = False,
) -> None:
    """Run the selected stages. No run manifest, no completion gate."""

    _assert_stages_available(context, stages)
    context.run_root.mkdir(parents=True, exist_ok=True)
    for stage in stages:
        run_stage(context, stage, force=force)
    if export and "type" in stages:
        export_qupath(context)
    print("[complete] " + ", ".join(stages))


def _completion_stage_order(context: RunContext) -> tuple[str, ...]:
    """Stages that define a complete content-addressed run in the current workflow."""

    if context.registry.production_status == NO_VERIFIED_TYPING_STATUS:
        return DEFAULT_STAGE_ORDER
    return STAGE_ORDER


def status(context: RunContext) -> int:
    """Print what each stage has on disk. Always returns 0: status reports, it does not gate."""

    print(f"data root {context.config.data_dir}")
    for stage in STAGE_ORDER:
        if (
            stage in TYPING_DEPENDENT_STAGES
            and context.registry.production_status == NO_VERIFIED_TYPING_STATUS
        ):
            print(f"  {stage:11s} BLOCKED {NO_VERIFIED_TYPING_REASON}")
            continue
        output = context.output_root(stage)
        present = sorted(
            path.name.split("=", 1)[1]
            for path in output.glob("donor_id=*")
        ) if output.is_dir() else []
        if not output.is_dir() or not any(output.iterdir()):
            if stage == "states":
                print(f"  {stage:11s} OPTIONAL (explicit-only)")
            else:
                print(f"  {stage:11s} ABSENT")
            continue
        expected = len(context.donors)
        missing = [donor for donor in context.donors if donor not in present]
        if present and missing:
            print(
                f"  {stage:11s} PARTIAL {len(present):>3}/{expected:<3} donors"
                f"  missing={missing}"
            )
        elif present:
            print(f"  {stage:11s} PRESENT {len(present):>3}/{expected:<3} donors")
        else:
            print(f"  {stage:11s} PRESENT {output}")
    if context.registry.production_status == NO_VERIFIED_TYPING_STATUS:
        print(f"  {'qupath':11s} BLOCKED {NO_VERIFIED_TYPING_REASON}")
    else:
        print(
            f"  {'qupath':11s} "
            f"{'PRESENT' if _validate_export(context) else 'ABSENT'}"
        )
    return 0


def _selected_stages(only: Sequence[str] | None, through: str | None) -> tuple[str, ...]:
    if only and through:
        raise ValueError("--only and --through cannot be combined")
    if only:
        return tuple(only)
    if through:
        return STAGE_ORDER[: STAGE_ORDER.index(through) + 1]
    return DEFAULT_STAGE_ORDER


def _assert_stages_available(
    context: RunContext,
    stages: Sequence[str],
    *,
    operation: str = "pipeline stage selection",
) -> None:
    """Fail before output creation when verified typing is unavailable."""

    _assert_registry_stages_available(
        context.registry,
        stages,
        operation=operation,
    )


def _assert_registry_stages_available(
    registry: MarkerRegistry,
    stages: Sequence[str],
    *,
    operation: str = "pipeline stage selection",
) -> None:
    """Registry-only preflight for CLIs that have not resolved cohort data."""

    requested = sorted(TYPING_DEPENDENT_STAGES.intersection(stages))
    if requested and registry.production_status == NO_VERIFIED_TYPING_STATUS:
        raise ContractError(
            f"{operation} is blocked by {NO_VERIFIED_TYPING_REASON}: "
            "the marker registry declares no verified production pairs; "
            f"requested={requested}"
        )


def main(argv=None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    if raw and raw[0] == "manifest":
        from .qupath_manifest import main as manifest_main

        return manifest_main(raw[1:])
    if not raw or raw[0].startswith("-"):
        raw.insert(0, "run")

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--config", type=Path)
    run.add_argument("--jobs", type=int)
    run.add_argument("--only", nargs="+", choices=STAGE_ORDER)
    run.add_argument("--through", choices=STAGE_ORDER)
    run.add_argument("--no-export", action="store_true")
    run.add_argument(
        "--force",
        action="store_true",
        help="rerun stages whose output is already present",
    )
    run.add_argument(
        "--pipelined",
        "--donor-pipeline",
        dest="pipelined",
        action="store_true",
    )
    run.add_argument("--geometry-workers", type=int)
    run.add_argument("--redsea-workers", type=int)
    run.add_argument("--downstream-workers", type=int)
    show = subparsers.add_parser("status")
    show.add_argument("--config", type=Path)
    export_parser = subparsers.add_parser("export")
    export_parser.add_argument("--config", type=Path)
    args = parser.parse_args(raw)

    overrides = {}
    if getattr(args, "jobs", None) is not None:
        overrides["n_jobs"] = args.jobs
    worker_overrides = {
        "geometry_workers": "pipeline_geometry_workers",
        "redsea_workers": "pipeline_redsea_workers",
        "downstream_workers": "pipeline_downstream_workers",
    }
    for argument, config_name in worker_overrides.items():
        value = getattr(args, argument, None)
        if value is not None:
            overrides[config_name] = value
    cfg = load_config(args.config, **overrides)
    if args.command == "status" and not Path(cfg.qupath_manifest).exists():
        # `status` is the first thing anyone runs, including on a checkout with no data at
        # all, so answering "there is no cohort manifest yet, here is how to make one" beats
        # a FileNotFoundError traceback out of the artifact reader. Every other subcommand
        # still hard-fails, because they cannot do anything useful without the manifest.
        print(f"no cohort manifest at {cfg.qupath_manifest}")
        print("  phenocycler manifest template --out data/qupath_manifest_spec.json")
        print("  phenocycler manifest create --spec data/qupath_manifest_spec.json "
              "--out data/qupath_manifest.json")
        return 0
    stages: tuple[str, ...] | None = None
    registry = load_registry(cfg.marker_registry)
    if args.command == "export":
        _assert_registry_stages_available(
            registry,
            ("type",),
            operation="QuPath export",
        )
    elif args.command == "run":
        stages = _selected_stages(args.only, args.through)
        if getattr(args, "pipelined", False) and args.only:
            parser.error("--pipelined supports --through, not --only")
        _assert_registry_stages_available(registry, stages)
    context = resolve_run_context(cfg)
    if args.command == "status":
        return status(context)
    if args.command == "export":
        export_qupath(context)
        return 0
    assert stages is not None
    _assert_stages_available(context, stages)
    if getattr(args, "pipelined", False):
        from .donor_pipeline import run_pipelined_pipeline

        run_pipelined_pipeline(
            context,
            stages=stages,
            export=not args.no_export,
            force=bool(getattr(args, "force", False)),
        )
        return 0
    run_pipeline(
        context,
        stages=stages,
        export=not args.no_export,
        force=bool(getattr(args, "force", False)),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
