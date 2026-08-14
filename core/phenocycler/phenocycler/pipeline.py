#!/usr/bin/env python3
"""Content-addressed QuPath → REDSEA → expression → cell-type workflow."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import pandas as pd

from .artifacts import (
    CodeMetadata,
    ContractError,
    FileFingerprint,
    InputArtifact,
    PartialArtifactError,
    QuPathCohortManifest,
    RunManifest,
    StageManifest,
    StaleArtifactError,
    sha256_json,
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
    NO_VERIFIED_PAIRS_STATUS,
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
PAIR_DEPENDENT_STAGES = frozenset({"controls", "calibrate", "type"})
NO_VERIFIED_PAIRS_REASON = "NO_VERIFIED_PAIRS"
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
        "phenocycler/cells_parquet.py",
        "phenocycler/artifacts.py",
        "phenocycler/cohort.py",
        "phenocycler/marker_registry.py",
        "phenocycler/marker_registry.json",
    ),
    "geometry": (
        "phenocycler/donor_pipeline.py",
        "phenocycler/geometry_stage.py",
        "phenocycler/geometry_qc.py",
        "phenocycler/redsea.py",
        "phenocycler/artifacts.py",
    ),
    "redsea": (
        "phenocycler/donor_pipeline.py",
        "phenocycler/redsea.py",
        "phenocycler/gpu.py",
        "phenocycler/parallel.py",
        "phenocycler/marker_taxonomy.py",
        "phenocycler/marker_registry.py",
    ),
    "expression": (
        "phenocycler/donor_pipeline.py",
        "phenocycler/expression.py",
        "phenocycler/marker_registry.py",
    ),
    "controls": (
        "phenocycler/donor_pipeline.py",
        "phenocycler/reference_controls.py",
        "phenocycler/marker_registry.py",
    ),
    "calibrate": (
        "phenocycler/donor_pipeline.py",
        "phenocycler/reference_controls.py",
        "phenocycler/restore_seeded.py",
        "phenocycler/restore_stage.py",
        "phenocycler/restore_validation.py",
        "phenocycler/marker_registry.py",
    ),
    "type": (
        "phenocycler/donor_pipeline.py",
        "phenocycler/production_typing.py",
        "phenocycler/simultaneous_typing.py",
        "phenocycler/broad_gates.json",
        "phenocycler/marker_registry.py",
    ),
    "states": (
        "phenocycler/donor_pipeline.py",
        "phenocycler/state_markers.py",
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
    code: CodeMetadata
    identity_config: dict
    run_id: str

    @property
    def run_root(self) -> Path:
        return self.config.data_dir

    def output_root(self, stage: str) -> Path:
        return Path(getattr(self.config, STAGE_OUTPUT_ATTRIBUTE[stage]))

    def stage_manifest_path(self, stage: str) -> Path:
        return self.config.manifests_dir / f"{stage}.json"


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
        active_registry.production_status != NO_VERIFIED_PAIRS_STATUS
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
            expression_reviews_content_id = InputArtifact.from_path(
                "restore_expression_reviews",
                cfg.restore_expression_reviews,
            ).content_sha256
    return {
        "run_schema_version": RUN_SCHEMA_VERSION,
        "typing_model": (
            {"mode": "restore_simultaneous"}
            if verified_pairs
            else {
                "mode": "unavailable",
                "reason": NO_VERIFIED_PAIRS_REASON,
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
                "reason": NO_VERIFIED_PAIRS_REASON,
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
            "norm_form": cfg.redsea_norm_form,
            "exclude_no_spillover": cfg.redsea_exclude_no_spillover,
            "output_schema": "Nucleus/Cytoplasm/Membrane/Cell",
        },
        "reference_controls": (
            asdict(ReferenceControlConfig())
            if verified_pairs
            else {"mode": "blocked", "reason": NO_VERIFIED_PAIRS_REASON}
        ),
        "phenotyping": (
            {
                "method": "RESTORE",
                "method_version": RESTORE_STAGE_VERSION,
                **asdict(verified_pair_config),
                "expression_reviews_content_id": expression_reviews_content_id,
            }
            if verified_pair_config is not None
            else {"mode": "blocked", "reason": NO_VERIFIED_PAIRS_REASON}
        ),
    }


def resolve_run_context(cfg: PipelineConfig) -> RunContext:
    """Resolve and validate every source of run identity before writing output."""

    cohort = QuPathCohortManifest.read_json(cfg.qupath_manifest)
    cohort.validate_current(mode="fast")
    donors = tuple(
        ensure_eligible_donors(
            cohort.expected_donors,
            context="content-addressed Phenocycler run",
        )
    )
    if donors != cohort.expected_donors:
        raise ContractError(
            "cohort donor order/set changed after exclusion validation"
        )
    registry = load_registry(cfg.marker_registry)
    broad_gate_policy = (
        None
        if registry.production_status == NO_VERIFIED_PAIRS_STATUS
        else load_broad_gate_policy(
            cfg.broad_gate_rules,
            registry=registry,
        )
    )
    repository = Path(__file__).resolve().parents[1]
    code = CodeMetadata.capture(
        repository,
        source_paths=("phenocycler", "pyproject.toml"),
    )
    identity = _identity_config(
        cfg,
        registry=registry,
        broad_gate_policy=broad_gate_policy,
    )
    run_id = sha256_json(
        {
            "cohort": cohort.content_id,
            "registry": registry.fingerprint,
            "broad_gate_rules": (
                "" if broad_gate_policy is None else broad_gate_policy.fingerprint
            ),
            "configuration": identity,
            "code": code.source_sha256,
        }
    )[:20]
    run_root = cfg.runs_dir / run_id
    run_cfg = replace(cfg, data_dir=run_root)
    return RunContext(
        base_config=cfg,
        config=run_cfg,
        cohort=cohort,
        registry=registry,
        broad_gate_policy=broad_gate_policy,
        donors=donors,
        repository=repository,
        code=code,
        identity_config=identity,
        run_id=run_id,
    )


def _stage_config(context: RunContext, stage: str) -> dict:
    return {
        "run_id": context.run_id,
        "stage": stage,
        "method_version": STAGE_METHODS[stage],
        "configuration": context.identity_config,
        "registry_fingerprint": context.registry.fingerprint,
        "broad_gate_rules_fingerprint": (
            context.broad_gate_policy.fingerprint
            if context.broad_gate_policy is not None
            else ""
        ),
    }


def _stage_code(context: RunContext, stage: str) -> CodeMetadata:
    return CodeMetadata.capture(
        context.repository,
        source_paths=STAGE_SOURCE_PATHS[stage],
    )


def _stage_external_inputs(
    context: RunContext,
    stage: str,
) -> tuple[InputArtifact, ...]:
    """Capture non-stage inputs without requiring completed predecessors."""

    inputs: list[InputArtifact] = []
    if stage in {"ingest", "geometry", "redsea"}:
        inputs.append(
            InputArtifact.from_qupath_cohort_manifest(
                "qupath_cohort", context.base_config.qupath_manifest
            )
        )
    if stage in {"ingest", "redsea", "expression", "controls", "calibrate", "type"}:
        inputs.append(
            InputArtifact.from_path(
                "marker_registry", context.base_config.marker_registry
            )
        )
    if stage == "calibrate":
        if context.base_config.restore_expression_reviews is not None:
            inputs.append(
                InputArtifact.from_path(
                    "restore_expression_reviews",
                    context.base_config.restore_expression_reviews,
                )
            )
    if stage == "type":
        inputs.append(
            InputArtifact.from_path(
                "broad_gate_rules",
                context.base_config.broad_gate_rules,
            )
        )
    return tuple(inputs)


def _stage_inputs(context: RunContext, stage: str) -> tuple[InputArtifact, ...]:
    inputs: list[InputArtifact] = []
    for dependency in STAGE_DEPENDENCIES[stage]:
        path = context.stage_manifest_path(dependency)
        if not path.exists():
            raise PartialArtifactError(
                f"{stage} requires completed stage {dependency!r}"
            )
        inputs.append(InputArtifact.from_stage_manifest(dependency, path))
    inputs.extend(_stage_external_inputs(context, stage))
    return tuple(inputs)


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


def run_stage(context: RunContext, stage: str) -> StageManifest:
    """Run or validate one immutable stage."""

    _assert_stages_available(context, (stage,))

    manifest_path = context.stage_manifest_path(stage)
    if manifest_path.exists():
        manifest = StageManifest.read_json(manifest_path)
        manifest.validate_current(
            config=_stage_config(context, stage),
            validation_mode="fast",
        )
        print(f"[current] {stage}: {manifest.output.total_rows:,} rows")
        return manifest

    output = context.output_root(stage)
    if output.exists() and any(output.iterdir()):
        raise PartialArtifactError(
            f"{stage}: output exists without a valid stage manifest: {output}. "
            "This content-addressed run is incomplete; inspect it before "
            "explicitly removing only that stage directory."
        )
    inputs = _stage_inputs(context, stage)
    # Validate every prerequisite before the stage is allowed to create output.
    # StageManifest.create repeats this check after execution, closing the race
    # where an input changes while the runner is active.
    for input_artifact in inputs:
        input_artifact.validate_current(mode="fast")
    print(f"[run] {stage}")
    STAGE_RUNNERS[stage](context)
    manifest = StageManifest.create(
        stage=stage,
        method_version=STAGE_METHODS[stage],
        expected_donors=context.donors,
        inputs=inputs,
        output_root=output,
        sidecar_paths=tuple(
            context.run_root / relative
            for relative in STAGE_SIDECARS.get(stage, ())
        ),
        config=_stage_config(context, stage),
        code=_stage_code(context, stage),
    )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_json(manifest_path)
    print(
        f"[done] {stage}: {manifest.output.total_rows:,} rows, "
        f"manifest={manifest.content_id[:12]}"
    )
    return manifest


def _export_manifest_path(context: RunContext) -> Path:
    return context.config.manifests_dir / "qupath_export.json"


def _validate_export(context: RunContext) -> bool:
    path = _export_manifest_path(context)
    if not path.exists():
        return False
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ContractError("QuPath export manifest must be a JSON object")
    content_id = payload.pop("content_id", None)
    if not isinstance(content_id, str) or not content_id:
        raise ContractError("QuPath export manifest is missing content_id")
    if sha256_json(payload) != content_id:
        raise StaleArtifactError("QuPath export manifest identity is invalid")
    if payload.get("run_id") != context.run_id:
        raise StaleArtifactError("QuPath export belongs to another run")
    type_manifest_path = context.stage_manifest_path("type")
    if not type_manifest_path.exists():
        raise StaleArtifactError("QuPath export has no current type stage")
    type_manifest = StageManifest.read_json(type_manifest_path)
    type_manifest.validate_current(
        config=_stage_config(context, "type"),
        validation_mode="fast",
    )
    if payload.get("assignment_stage_content_id") != type_manifest.content_id:
        raise StaleArtifactError(
            "QuPath export does not match the current assignment stage"
        )
    if tuple(payload.get("expected_donors", ())) != context.donors:
        raise StaleArtifactError(
            "QuPath export donor contract differs from the run"
        )
    files = payload.get("files")
    if not isinstance(files, list):
        raise ContractError("QuPath export manifest files must be a list")
    if len(files) != len(context.donors):
        raise PartialArtifactError(
            "QuPath export does not contain exactly one file per donor"
        )
    for index, item in enumerate(files):
        try:
            fingerprint = FileFingerprint.from_dict(item)
        except (KeyError, TypeError, ValueError) as exc:
            raise ContractError(
                f"QuPath export file fingerprint {index} is malformed"
            ) from exc
        fingerprint.validate_current(mode="fast")
    return True


def export_qupath(context: RunContext) -> None:
    """Write and fingerprint one uncertainty-preserving CSV per donor."""

    _assert_stages_available(context, ("type",), operation="QuPath export")

    if _validate_export(context):
        print("[current] qupath export")
        return
    type_manifest = context.stage_manifest_path("type")
    if not type_manifest.exists():
        raise PartialArtifactError("QuPath export requires the type stage")
    current_type = StageManifest.read_json(type_manifest)
    current_type.validate_current(
        config=_stage_config(context, "type"),
        validation_mode="fast",
    )
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
        "run_id": context.run_id,
        "expected_donors": list(context.donors),
        "assignment_stage_content_id": current_type.content_id,
        "files": [asdict(FileFingerprint.capture(path)) for path in paths],
    }
    payload["content_id"] = sha256_json(payload)
    manifest_path = _export_manifest_path(context)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        raise FileExistsError(manifest_path)
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
) -> RunManifest | None:
    """Run selected stages and write a complete run manifest when possible."""

    _assert_stages_available(context, stages)
    context.run_root.mkdir(parents=True, exist_ok=True)
    for stage in stages:
        run_stage(context, stage)
    if export and "type" in stages:
        export_qupath(context)

    completion_stages = _completion_stage_order(context)
    all_paths = [context.stage_manifest_path(stage) for stage in completion_stages]
    if not all(path.exists() for path in all_paths):
        return None
    # Existence is not completion. Validate every required stage—including
    # stages not selected in this invocation—before writing run.json or
    # advancing LATEST. In development_no_verified_pairs, the content-addressed
    # run product intentionally ends at expression; states remains an optional
    # add-on and blocked pair stages cannot be completion requirements.
    for stage, path in zip(completion_stages, all_paths, strict=True):
        manifest = StageManifest.read_json(path)
        if manifest.stage != stage:
            raise ContractError(
                f"manifest at {path} declares stage {manifest.stage!r}, "
                f"expected {stage!r}"
            )
        manifest.validate_current(
            config=_stage_config(context, stage),
            validation_mode="fast",
        )
    run_config = {
        **context.identity_config,
        "cohort_content_id": context.cohort.content_id,
        "registry_fingerprint": context.registry.fingerprint,
        "broad_gate_rules_fingerprint": (
            ""
            if context.broad_gate_policy is None
            else context.broad_gate_policy.fingerprint
        ),
    }
    run_manifest = RunManifest.create(
        run_name=context.run_id,
        expected_donors=context.donors,
        stage_manifest_paths=all_paths,
        config=run_config,
        code=context.code,
    )
    run_manifest_path = context.run_root / "run.json"
    if run_manifest_path.exists():
        current = RunManifest.read_json(run_manifest_path)
        current.validate_current(
            config=run_config,
            validation_mode="fast",
        )
        if current.content_id != run_manifest.content_id:
            raise StaleArtifactError(
                "completed run manifest differs from the current stage set"
            )
        run_manifest = current
    else:
        run_manifest.write_json(run_manifest_path)
    context.base_config.runs_dir.mkdir(parents=True, exist_ok=True)
    (context.base_config.runs_dir / "LATEST").write_text(
        context.run_id + "\n", encoding="utf-8"
    )
    print(f"[complete] run {context.run_id}")
    return run_manifest


def _completion_stage_order(context: RunContext) -> tuple[str, ...]:
    """Stages that define a complete content-addressed run in the current workflow."""

    if context.registry.production_status == NO_VERIFIED_PAIRS_STATUS:
        return DEFAULT_STAGE_ORDER
    return STAGE_ORDER


def status(context: RunContext) -> int:
    """Print evidence-backed stage status; return nonzero for stale/partial state."""

    print(f"run_id {context.run_id}")
    healthy = True
    for stage in STAGE_ORDER:
        if (
            stage in PAIR_DEPENDENT_STAGES
            and context.registry.production_status == NO_VERIFIED_PAIRS_STATUS
        ):
            print(f"  {stage:11s} BLOCKED {NO_VERIFIED_PAIRS_REASON}")
            continue
        path = context.stage_manifest_path(stage)
        if not path.exists():
            if stage == "states":
                print(f"  {stage:11s} OPTIONAL (explicit-only)")
                continue
            if stage != "ingest":
                from .donor_pipeline import donor_progress

                completed, expected = donor_progress(context, stage)
            else:
                completed, expected = 0, len(context.donors)
            if completed:
                print(
                    f"  {stage:11s} IN PROGRESS "
                    f"{completed:>3}/{expected:<3} donor receipts"
                )
            else:
                print(f"  {stage:11s} MISSING")
            healthy = False
            continue
        try:
            manifest = StageManifest.read_json(path)
            manifest.validate_current(
                config=_stage_config(context, stage),
                validation_mode="fast",
            )
        except (ContractError, OSError, ValueError) as exc:
            print(f"  {stage:11s} STALE   {exc}")
            healthy = False
        else:
            print(
                f"  {stage:11s} CURRENT "
                f"{manifest.output.total_rows:>12,} rows "
                f"{manifest.content_id[:12]}"
            )
    if context.registry.production_status == NO_VERIFIED_PAIRS_STATUS:
        print(f"  {'qupath':11s} BLOCKED {NO_VERIFIED_PAIRS_REASON}")
    else:
        try:
            export_current = _validate_export(context)
        except (ContractError, OSError, ValueError) as exc:
            print(f"  {'qupath':11s} STALE   {exc}")
            healthy = False
        else:
            print(f"  {'qupath':11s} {'CURRENT' if export_current else 'MISSING'}")
            healthy &= export_current
    return 0 if healthy else 1


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
    """Fail before output creation when pair-derived production is unavailable."""

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

    requested = sorted(PAIR_DEPENDENT_STAGES.intersection(stages))
    if requested and registry.production_status == NO_VERIFIED_PAIRS_STATUS:
        raise ContractError(
            f"{operation} is blocked by {NO_VERIFIED_PAIRS_REASON}: "
            "the marker registry declares no verified production pairs; "
            f"requested={requested}"
        )


def main(argv=None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    if raw and raw[0] == "manifest":
        from .qupath_manifest import main as manifest_main

        return manifest_main(raw[1:])
    if raw and raw[0] == "resume":
        from .resume import main as resume_main

        return resume_main(raw[1:])
    if raw and raw[0] == "recover":
        from .recovery import main as recovery_main

        return recovery_main(raw[1:])
    if raw and raw[0] == "adopt":
        from .adoption import main as adoption_main

        return adoption_main(raw[1:])
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
        )
        return 0
    run_pipeline(context, stages=stages, export=not args.no_export)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
