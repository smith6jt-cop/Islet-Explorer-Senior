"""Pipeline-native RESTORE marker calling.

This module reconnects the validated RESTORE v10 pair fit to the current
content-addressed pipeline. RESTORE estimates one donor-local divisor for each
versioned target/reference pair, emits a tri-state marker call for every cell,
and presents those calls in the wide marker-evidence contract consumed by
ordered residual typing.

No reviewed cell labels or fitted object classifier enter this stage.  A
technical donor-marker failure becomes ``unavailable`` (or halts, when
explicitly configured); it is never converted to a biological negative.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from .artifacts import ContractError
from .marker_registry import MarkerRegistry, NoVerifiedTypingError
from .reference_controls import (
    reference_masks_as_arrays,
    reference_masks_from_wide,
    validate_reference_masks,
)
from .restore_seeded import (
    FEATURE_TRANSFORMS,
    evaluate_reference_seeded_pair,
)
from .restore_seeded import METHOD_VERSION as SEEDED_RESTORE_METHOD_VERSION
from .restore_validation import (
    ExpressionReview,
    PairReview,
    PairState,
    PairValidationConfig,
    load_expression_reviews,
    locked_pair_metrics,
)

# v4 consumes the registry-selected expression table for call intensities and
# routes NNMF compartments independently per marker. QuPath-normalized markers
# are left-censored at zero exactly once at this RESTORE boundary.
METHOD_VERSION = f"restore-marker-evidence-v4-selected-inputs+{SEEDED_RESTORE_METHOD_VERSION}"
CELL_POSITIVE = "positive"
CELL_NEGATIVE = "negative"
CELL_UNAVAILABLE = "unavailable"
CELL_STATES = (
    CELL_NEGATIVE,
    "indeterminate",
    CELL_POSITIVE,
    CELL_UNAVAILABLE,
)
CALLABLE_STATES = frozenset(
    {PairState.VALID, PairState.NO_TARGET_EXPRESSION_CONFIRMED}
)


@dataclass(frozen=True)
class RestoreStageConfig:
    """Scientific policy for the production RESTORE application."""

    # No defaults for the two scientific policy fields. Their historical defaults
    # were `p99` and ("Membrane", "Cell"), both of which docs/METHOD_HISTORY.md
    # records as measurably wrong -- p99 called 2.58x each pair's own anchor, and
    # a rank-2 fit on Membrane+Cell reproduces its own input to 5e-5 because
    # `Cell` is the sum of the others. A caller that forgets to pass them must
    # fail closed rather than silently revive the rejected settings, so these are
    # required. `from_pipeline` supplies them from config.ini.
    divisor_statistic: str
    feature_compartments: tuple[str, ...]
    control_definition: str = "reference_only"
    on_technical_failure: str = "unavailable"
    nnmf_fit_basis: str = "all"
    control_membership: str = "component"
    feature_transform: str = "linear_upper_scale"

    def __post_init__(self) -> None:
        # PairValidationConfig owns the detailed validation and the full frozen
        # specification. Constructing it here also catches drift immediately.
        self.pair_config()
        if self.on_technical_failure not in {"halt", "unavailable"}:
            raise ValueError(
                "on_technical_failure must be 'halt' or 'unavailable'"
            )
        if self.control_definition != "reference_only":
            raise ValueError(
                "production RESTORE requires control_definition='reference_only'"
            )
        if self.divisor_statistic != "empirical_alpha":
            raise ValueError(
                "production RESTORE requires divisor_statistic='empirical_alpha'; "
                "fixed quantiles and maxima are review-only"
            )
        if self.nnmf_fit_basis != "all":
            raise ValueError(
                "production RESTORE requires nnmf_fit_basis='all'; "
                "raw-floor/balanced-arm fitting is review-only"
            )
        if self.control_membership != "component":
            raise ValueError(
                "production RESTORE requires control_membership='component'"
            )
        if self.feature_transform not in FEATURE_TRANSFORMS:
            raise ValueError(
                "feature_transform must be one of "
                f"{sorted(FEATURE_TRANSFORMS)}"
            )
        if self.feature_transform != "linear_upper_scale":
            raise ValueError(
                "production RESTORE requires feature_transform="
                "'linear_upper_scale'; shifted and nonlinear transforms are "
                "review-only"
            )

    @classmethod
    def from_pipeline(cls, config: object) -> RestoreStageConfig:
        # These two are `None` on an unconfigured PipelineConfig, and `str(None)`
        # is the string "None" -- which would reach PairValidationConfig as a
        # feature compartment named "None" and fail there with an error that says
        # nothing about the missing configuration. Name the real problem here.
        missing = [
            key
            for key in ("restore_divisor_statistic", "restore_feature_compartments")
            if getattr(config, key, None) is None
        ]
        if missing:
            raise ValueError(
                "RESTORE scientific policy is unset: "
                + ", ".join(sorted(missing))
                + ". These have no defaults; supply them from config.ini's "
                "[restore] section."
            )
        compartments = tuple(
            item.strip()
            for item in str(
                config.restore_feature_compartments
            ).split(",")
            if item.strip()
        )
        return cls(
            control_definition=str(
                config.restore_control_definition
            ),
            divisor_statistic=str(
                config.restore_divisor_statistic
            ),
            feature_compartments=compartments,
            on_technical_failure=str(
                config.restore_on_technical_failure
            ),
            nnmf_fit_basis=str(config.restore_nnmf_fit_basis),
            control_membership=str(config.restore_control_membership),
            feature_transform=str(config.restore_feature_transform),
        )

    def pair_config(self) -> PairValidationConfig:
        return PairValidationConfig(
            control_definition=self.control_definition,
            divisor_statistic=self.divisor_statistic,
            feature_compartments=self.feature_compartments,
            nnmf_fit_basis=self.nnmf_fit_basis,
            control_membership=self.control_membership,
        )


@dataclass(frozen=True)
class RestoreStageResult:
    donor_id: str
    model_table: pd.DataFrame
    evidence_wide: pd.DataFrame


def production_pairs(
    registry: MarkerRegistry,
) -> tuple[tuple[str, str], ...]:
    """Refuse production marker calling until a pair roster is verified."""

    del registry
    raise NoVerifiedTypingError(
        "no verified production marker pairs exist; use the mutable fundamental "
        "QuPath input diagnostics while pair design is under development"
    )


def _keyed(
    frame: pd.DataFrame,
    *,
    name: str,
) -> pd.DataFrame:
    if "object_id" not in frame:
        raise ContractError(f"{name} is missing object_id")
    result = frame.copy()
    if result["object_id"].isna().any():
        raise ContractError(f"{name} contains null object_id")
    result["object_id"] = result["object_id"].astype(str)
    if result["object_id"].duplicated().any():
        examples = result.loc[
            result["object_id"].duplicated(keep=False), "object_id"
        ].head(5)
        raise ContractError(
            f"{name} contains duplicate object_id; examples={examples.tolist()}"
        )
    return result.set_index("object_id", drop=False)


def _align_exact(
    frame: pd.DataFrame,
    object_ids: pd.Index,
    *,
    name: str,
) -> pd.DataFrame:
    keyed = _keyed(frame, name=name)
    observed = keyed.index
    missing = object_ids.difference(observed)
    extra = observed.difference(object_ids)
    if len(missing) or len(extra):
        raise ContractError(
            f"{name} object universe differs from canonical cells "
            f"(missing={len(missing):,}, extra={len(extra):,}, "
            f"missing_examples={missing[:5].tolist()}, "
            f"extra_examples={extra[:5].tolist()})"
        )
    return keyed.reindex(object_ids)


def build_restore_input_frame(
    cells: pd.DataFrame,
    redsea: pd.DataFrame,
    selected_expression: pd.DataFrame,
    geometry: pd.DataFrame,
    *,
    donor_id: str,
    pairs: Sequence[tuple[str, str]],
    registry: MarkerRegistry,
    feature_compartments: Sequence[str] = ("Membrane", "Cell"),
) -> pd.DataFrame:
    """Build the aligned production RESTORE input contract.

    The selected-expression artifact is the sole threshold/call scale. NNMF
    compartments are routed independently by the registry: exact QuPath
    ``Mean (norm)`` measurements for ``qupath_norm`` markers, and REDSEA
    compartments for all other markers. Signed QuPath-normalized values are
    preserved upstream and finite negative values are left-censored to zero
    here because NNMF requires a nonnegative matrix. Missing values stay
    missing. No raw, REDSEA, or RESTORE-output fallback can satisfy a missing
    QuPath-normalized input.
    """

    donor = str(donor_id)
    marker_order = tuple(
        dict.fromkeys(marker for pair in pairs for marker in pair)
    )
    required_cells = {
        "object_id",
        "image",
        "X_centroid",
        "Y_centroid",
    }
    missing_cells = sorted(required_cells.difference(cells.columns))
    if missing_cells:
        raise ContractError(
            f"donor {donor}: canonical cells are missing RESTORE inputs "
            f"{missing_cells}"
        )

    canonical = _keyed(cells, name=f"donor {donor} cells")
    object_ids = canonical.index
    redsea_aligned = _align_exact(
        redsea, object_ids, name=f"donor {donor} REDSEA"
    )
    expression_aligned = _align_exact(
        selected_expression,
        object_ids,
        name=f"donor {donor} selected expression",
    )
    geometry_aligned = _align_exact(
        geometry, object_ids, name=f"donor {donor} geometry QC"
    )
    if "analysis_eligible" not in geometry_aligned:
        raise ContractError(
            f"donor {donor}: geometry QC is missing analysis_eligible"
        )
    if "estimation_eligible" not in geometry_aligned:
        raise ContractError(
            f"donor {donor}: geometry QC is missing estimation_eligible"
        )
    eligibility = geometry_aligned["analysis_eligible"]
    if eligibility.isna().any() or not eligibility.isin(
        [True, False, 0, 1]
    ).all():
        raise ContractError(
            f"donor {donor}: analysis_eligible must contain explicit booleans"
        )
    analysis_eligible = eligibility.astype(bool).to_numpy()
    estimation = geometry_aligned["estimation_eligible"]
    if estimation.isna().any() or not estimation.isin(
        [True, False, 0, 1]
    ).all():
        raise ContractError(
            f"donor {donor}: estimation_eligible must contain explicit booleans"
        )
    estimation_eligible = estimation.astype(bool).to_numpy()
    if np.any(estimation_eligible & ~analysis_eligible):
        raise ContractError(
            f"donor {donor}: estimation-eligible cells must also be analysis eligible"
        )

    columns: dict[str, object] = {
        "donor_id": np.full(len(canonical), donor, dtype=object),
        "object_id": object_ids.to_numpy(dtype=object),
        "image": canonical["image"].to_numpy(),
        "x_um": pd.to_numeric(
            canonical["X_centroid"], errors="coerce"
        ).to_numpy(),
        "y_um": pd.to_numeric(
            canonical["Y_centroid"], errors="coerce"
        ).to_numpy(),
        "qc_analysis_eligible": analysis_eligible,
        "qc_estimation_eligible": estimation_eligible,
    }
    measurement_sources: dict[str, str] = {}
    input_transforms: dict[str, str] = {}
    left_censored_counts: dict[str, int] = {}

    def restore_values(
        source: pd.Series,
        *,
        marker: str,
        measurement_source: str,
    ) -> np.ndarray:
        values = np.array(
            pd.to_numeric(source, errors="coerce").to_numpy(dtype=float),
            copy=True,
        )
        if measurement_source != "qupath_norm":
            return values
        finite_negative = np.isfinite(values) & (values < 0)
        if marker not in left_censored_counts:
            left_censored_counts[marker] = int(np.count_nonzero(finite_negative))
        values[finite_negative] = 0.0
        return values

    for marker in marker_order:
        try:
            spec = registry.marker(marker)
        except KeyError as exc:
            raise ContractError(
                f"donor {donor}: RESTORE marker {marker!r} is absent from the registry"
            ) from exc
        measurement_source = str(spec.measurement_source)
        measurement_sources[marker] = measurement_source
        input_transforms[marker] = (
            "qupath_norm_left_censor_zero"
            if measurement_source == "qupath_norm"
            else "none"
        )
        if marker not in expression_aligned:
            raise ContractError(
                f"donor {donor}: selected expression is missing marker {marker!r}"
            )
        if measurement_source == "qupath_norm":
            cell_source_name = f"QuPathNorm__Cell__{marker}"
            if cell_source_name not in canonical:
                raise ContractError(
                    f"donor {donor}: QuPath-normalized RESTORE marker {marker!r} "
                    f"requires exact cells column {cell_source_name!r}; raw, "
                    "REDSEA, and RESTORE-normalized fallbacks are forbidden"
                )
            selected_signed = pd.to_numeric(
                expression_aligned[marker], errors="coerce"
            ).to_numpy(dtype=float)
            canonical_signed = pd.to_numeric(
                canonical[cell_source_name], errors="coerce"
            ).to_numpy(dtype=float)
            if not np.array_equal(
                selected_signed, canonical_signed, equal_nan=True
            ):
                raise ContractError(
                    f"donor {donor}: selected expression for QuPath-normalized "
                    f"marker {marker!r} differs from exact cells column "
                    f"{cell_source_name!r}"
                )
        columns[f"selected__{marker}"] = restore_values(
            expression_aligned[marker],
            marker=marker,
            measurement_source=measurement_source,
        )

        # NNMF features are marker-source routed, not pair-source routed.
        #
        # Chang et al. 2020, Methods: "If we consider a 2-dimensional feature
        # set (i.e., mean intensity profile), there exists a trivial solution,
        # i.e., W = Y and H = I." Measured here: with `Cell` alone the rank-2
        # fit reproduces its input to a relative error of 5e-5 to 2e-3 -- it is
        # that exact trivial solution and carries no information the two input
        # columns did not already have. With Nucleus+Cytoplasm+Membrane the
        # error is 0.21-0.25, i.e. a genuine approximation.
        #
        # `redsea.compensate_compartments` emits Nucleus/Cytoplasm/Membrane/Cell
        # for EVERY channel, so every accepted pair can have a >2-dimensional
        # feature set. The QuPath compartment export covered 12 markers, which
        # is why 279 of 479 pairs were falling back to whole cell alone.
        for compartment in feature_compartments:
            name = f"{compartment}__{marker}"
            if measurement_source == "qupath_norm":
                source_name = f"QuPathNorm__{compartment}__{marker}"
                if source_name not in canonical:
                    raise ContractError(
                        f"donor {donor}: QuPath-normalized RESTORE marker "
                        f"{marker!r} requires exact cells column {source_name!r}; "
                        "raw, REDSEA, and RESTORE-normalized fallbacks are forbidden"
                    )
                source = canonical[source_name]
            else:
                source_name = name
                if source_name not in redsea_aligned:
                    raise ContractError(
                        f"donor {donor}: RESTORE marker {marker!r} requires "
                        f"REDSEA compartment column {source_name!r}"
                    )
                source = redsea_aligned[source_name]
            columns[name] = restore_values(
                source,
                marker=marker,
                measurement_source=measurement_source,
            )

    result = pd.DataFrame(columns)
    measurement_columns = [
        column
        for column in result.columns
        if column.startswith("selected__")
        or any(
            column.startswith(f"{compartment}__")
            for compartment in feature_compartments
        )
    ]
    result.loc[~analysis_eligible, measurement_columns] = np.nan
    result.attrs["analysis_eligible_n"] = int(analysis_eligible.sum())
    result.attrs["estimation_eligible_n"] = int(estimation_eligible.sum())
    result.attrs["measurement_source_by_marker"] = measurement_sources
    result.attrs["input_transform_by_marker"] = input_transforms
    result.attrs["left_censored_n_by_marker"] = left_censored_counts
    return result


def cell_calls_from_evaluation(
    evaluation: Mapping[str, object],
    *,
    analysis_eligible: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map one RESTORE evaluation to aligned positive/norm/state arrays."""

    n = int(evaluation["n_input"])
    pair_state = PairState(str(evaluation["state"]))
    positive = np.zeros(n, dtype=bool)
    normalized = np.full(n, np.nan, dtype=float)
    state = np.full(n, CELL_UNAVAILABLE, dtype="<U13")

    if pair_state in CALLABLE_STATES:
        # ``valid_idx`` is the set of cells that carried finite, nonnegative
        # measurements for both markers.  A cell outside it has no measurement,
        # and a reviewed absence of target expression is a statement about the
        # donor's biology, not evidence about a cell that was never measured.
        # Both callable states therefore call only within ``valid_idx``.
        # The CALL universe, not the fit universe: a cell needs a finite
        # corrected intensity to be called, not a complete NNMF feature row.
        # Falls back to valid_idx for an evaluation produced before call_idx
        # existed (and for the compact dicts used in tests).
        available = np.zeros(n, dtype=bool)
        index = evaluation.get("call_idx", evaluation["valid_idx"])
        available[np.asarray(index, dtype=int)] = True
        positive = np.asarray(
            evaluation["full_lineage_positive"], dtype=bool
        ).copy()
        normalized = np.asarray(
            evaluation["full_normalized"], dtype=float
        ).copy()
        state[available & positive] = CELL_POSITIVE
        state[available & ~positive] = CELL_NEGATIVE

    if analysis_eligible is not None:
        eligible = np.asarray(analysis_eligible, dtype=bool)
        if eligible.shape != (n,):
            raise ValueError("analysis_eligible has the wrong shape")
        positive[~eligible] = False
        normalized[~eligible] = np.nan
        state[~eligible] = CELL_UNAVAILABLE
    return positive, normalized, state


def _technical_row(
    *,
    donor: str,
    target: str,
    reference: str,
    n_input: int,
    reason: str,
    config: RestoreStageConfig,
) -> dict[str, object]:
    return {
        "method_version": SEEDED_RESTORE_METHOD_VERSION,
        "pipeline_stage_method_version": METHOD_VERSION,
        "donor": donor,
        "target": target,
        "reference": reference,
        "expression_review": ExpressionReview.UNREVIEWED.value,
        "pair_review": PairReview.ACCEPTED.value,
        "state": PairState.TECHNICAL_FAILURE.value,
        "state_reason": reason,
        "n_input": int(n_input),
        "n_assigned": 0,
        "n_unavailable": int(n_input),
        "canonical_divisor": np.nan,
        "reference_n": 0,
        "divisor_statistic": config.divisor_statistic,
        "control_definition": config.control_definition,
        "feature_compartments": "",
    }


def restore_marker_evidence(
    donor_frame: pd.DataFrame,
    *,
    donor_id: str,
    registry: MarkerRegistry,
    config: RestoreStageConfig | None = None,
    expression_reviews: Mapping[
        tuple[str, str], ExpressionReview
    ] | None = None,
    reference_controls: Mapping[str, np.ndarray] | None = None,
    evaluator: Callable[..., dict] | None = None,
) -> RestoreStageResult:
    """Fit every production RESTORE pair and emit wide tri-state evidence."""

    donor = str(donor_id)
    if config is None:
        raise ContractError(
            "restore_marker_evidence requires an explicit RestoreStageConfig; "
            "there is no safe default for divisor_statistic or "
            "feature_compartments (see the class docstring)"
        )
    pairs = production_pairs(registry)
    pair_config = config.pair_config()
    reviews = {} if expression_reviews is None else dict(expression_reviews)
    if "qc_analysis_eligible" not in donor_frame:
        raise ContractError("RESTORE donor frame is missing qc_analysis_eligible")
    analysis_eligible = donor_frame["qc_analysis_eligible"].to_numpy(
        dtype=bool
    )
    n = len(donor_frame)

    evidence: dict[str, object] = {
        "donor_id": np.full(n, donor, dtype=object),
        "object_id": donor_frame["object_id"].astype(str).to_numpy(),
        "qc_analysis_eligible": analysis_eligible,
    }
    audit_rows: list[dict[str, object]] = []
    for target, reference in pairs:
        evaluation: dict | None = None
        try:
            if evaluator is None:
                if reference_controls is None or reference not in reference_controls:
                    raise ContractError(
                        f"donor {donor}: no validated reference-control mask "
                        f"for {reference}"
                    )
                evaluation = evaluate_reference_seeded_pair(
                    donor_frame,
                    donor,
                    target,
                    reference,
                    config=pair_config,
                    marker_alpha=registry.marker_alpha(target),
                    reference_positive=reference_controls[reference],
                    feature_transform=config.feature_transform,
                    expression_review=reviews.get(
                        (donor, target), ExpressionReview.UNREVIEWED
                    ),
                    pair_review=PairReview.UNREVIEWED,
                )
            else:
                # Explicit injection remains available for focused tests and
                # historical review tooling. Production never supplies it.
                evaluation = evaluator(
                    donor_frame,
                    donor,
                    target,
                    reference,
                    config=pair_config,
                    expression_review=reviews.get(
                        (donor, target), ExpressionReview.UNREVIEWED
                    ),
                    pair_review=PairReview.UNREVIEWED,
                )
            pair_state = PairState(evaluation["state"])
            if pair_state not in CALLABLE_STATES:
                reason = (
                    f"{pair_state.value}: "
                    f"{evaluation.get('state_reason', '')}"
                )
                if config.on_technical_failure == "halt":
                    raise RuntimeError(
                        f"donor {donor} {target} <- {reference}: {reason}"
                    )
            row = dict(
                evaluation.get("model_row")
                or locked_pair_metrics(evaluation)
            )
            positive, normalized, states = cell_calls_from_evaluation(
                evaluation, analysis_eligible=analysis_eligible
            )
        except Exception as exc:
            if config.on_technical_failure == "halt":
                raise
            reason = f"{type(exc).__name__}: {exc}"
            row = _technical_row(
                donor=donor,
                target=target,
                reference=reference,
                n_input=n,
                reason=reason,
                config=config,
            )
            positive = np.zeros(n, dtype=bool)
            normalized = np.full(n, np.nan, dtype=float)
            states = np.full(n, CELL_UNAVAILABLE, dtype="<U13")

        finite_norm = np.isfinite(normalized) & (normalized >= 0)
        log2_ratio = np.full(n, np.nan, dtype=float)
        log2_ratio[finite_norm] = np.log2(
            np.maximum(normalized[finite_norm], np.finfo(float).tiny)
        )
        model_valid = (
            np.isin(states, [CELL_POSITIVE, CELL_NEGATIVE])
            & analysis_eligible
        )
        selected = pd.to_numeric(
            donor_frame[f"selected__{target}"], errors="coerce"
        ).to_numpy(dtype=float)
        prefix = f"{target}__"
        evidence[prefix + "selected_intensity"] = selected
        evidence[prefix + "state"] = pd.Categorical(
            states, categories=CELL_STATES
        )
        evidence[prefix + "log2_threshold_ratio"] = log2_ratio
        evidence[prefix + "model_valid"] = model_valid
        evidence[prefix + "restore_positive"] = positive
        evidence[prefix + "restore_normalized"] = normalized
        evidence[prefix + "empirical_tail_probability"] = np.asarray(
            evaluation.get(
                "full_empirical_tail_probability",
                np.full(n, np.nan, dtype=float),
            )
            if evaluation is not None
            else np.full(n, np.nan, dtype=float),
            dtype=float,
        )
        evidence[prefix + "censored_below"] = np.asarray(
            evaluation.get(
                "full_censored_below", np.zeros(n, dtype=bool)
            )
            if evaluation is not None
            else np.zeros(n, dtype=bool),
            dtype=bool,
        )
        evidence[prefix + "calibration_status"] = np.full(
            n, str(row["state"]), dtype=object
        )
        evidence[prefix + "reference"] = np.full(
            n, reference, dtype=object
        )
        evidence[prefix + "threshold"] = np.full(
            n,
            float(row.get("canonical_divisor", np.nan) or np.nan),
            dtype=float,
        )
        evidence[prefix + "n_controls"] = np.full(
            n, int(row.get("reference_n", 0) or 0), dtype=np.int32
        )
        evidence[prefix + "method_version"] = np.full(
            n, METHOD_VERSION, dtype=object
        )

        row["called_positive_n"] = int(
            np.count_nonzero(states == CELL_POSITIVE)
        )
        row.setdefault("donor", donor)
        row.setdefault("target", target)
        row.setdefault("reference", reference)
        row.setdefault("method_version", SEEDED_RESTORE_METHOD_VERSION)
        row["called_negative_n"] = int(
            np.count_nonzero(states == CELL_NEGATIVE)
        )
        row["called_unavailable_n"] = int(
            np.count_nonzero(states == CELL_UNAVAILABLE)
        )
        row["analysis_eligible_n"] = int(analysis_eligible.sum())
        sources = donor_frame.attrs.get("measurement_source_by_marker", {})
        transforms = donor_frame.attrs.get("input_transform_by_marker", {})
        censored_counts = donor_frame.attrs.get(
            "left_censored_n_by_marker", {}
        )
        row["target_measurement_source"] = str(sources.get(target, "unknown"))
        row["reference_measurement_source"] = str(
            sources.get(reference, "unknown")
        )
        row["target_input_transform"] = str(transforms.get(target, "unknown"))
        row["reference_input_transform"] = str(
            transforms.get(reference, "unknown")
        )
        row["target_left_censored_n"] = int(censored_counts.get(target, 0))
        row["reference_left_censored_n"] = int(
            censored_counts.get(reference, 0)
        )
        row["pipeline_stage_method_version"] = METHOD_VERSION
        audit_rows.append(row)

    evidence_wide = pd.DataFrame(evidence)
    model_table = pd.DataFrame.from_records(audit_rows).sort_values(
        ["donor", "target"], kind="stable", ignore_index=True
    )
    return RestoreStageResult(
        donor_id=donor,
        model_table=model_table,
        evidence_wide=evidence_wide,
    )


def _single_partition(root: Path, donor_id: str) -> Path:
    partition = Path(root) / f"donor_id={donor_id}"
    files = sorted(partition.glob("*.parquet"))
    if len(files) != 1:
        raise ContractError(
            f"expected exactly one parquet for donor {donor_id} in {root}; "
            f"found {len(files)}"
        )
    return files[0]


def _read_columns(
    root: Path,
    donor_id: str,
    columns: Sequence[str],
    *,
    allow_missing: bool = False,
) -> pd.DataFrame:
    path = _single_partition(root, donor_id)
    available = set(pq.ParquetFile(path).schema_arrow.names)
    missing = sorted(set(columns).difference(available))
    if missing and not allow_missing:
        raise ContractError(f"{path} is missing required columns {missing}")
    selected = [column for column in columns if column in available]
    return pd.read_parquet(path, columns=selected)


def restore_input_frame_from_paths(
    *,
    donor_id: str,
    cells_root: Path,
    redsea_root: Path,
    expression_root: Path,
    geometry_root: Path,
    registry: MarkerRegistry,
    config: RestoreStageConfig,
    pairs: Sequence[tuple[str, str]] | None = None,
) -> pd.DataFrame:
    """Read and assemble one donor's source-routed RESTORE input frame."""

    donor = str(donor_id)
    selected_pairs = production_pairs(registry) if pairs is None else tuple(pairs)
    markers = tuple(
        dict.fromkeys(marker for pair in selected_pairs for marker in pair)
    )
    qupath_norm_markers = tuple(
        marker
        for marker in markers
        if registry.marker(marker).measurement_source == "qupath_norm"
    )
    redsea_markers = tuple(
        marker for marker in markers if marker not in qupath_norm_markers
    )
    cells = _read_columns(
        cells_root,
        donor,
        [
            "object_id",
            "image",
            "X_centroid",
            "Y_centroid",
            *(
                f"QuPathNorm__{compartment}__{marker}"
                for compartment in dict.fromkeys(
                    ("Cell", *config.feature_compartments)
                )
                for marker in qupath_norm_markers
            ),
        ],
    )
    expression = _read_columns(
        expression_root,
        donor,
        ["object_id", *markers],
    )
    redsea = _read_columns(
        redsea_root,
        donor,
        [
            "object_id",
            *(
                f"{compartment}__{marker}"
                for compartment in config.feature_compartments
                for marker in redsea_markers
            ),
        ],
    )
    geometry = _read_columns(
        geometry_root,
        donor,
        ["object_id", "analysis_eligible", "estimation_eligible"],
    )
    return build_restore_input_frame(
        cells,
        redsea,
        expression,
        geometry,
        donor_id=donor,
        pairs=selected_pairs,
        registry=registry,
        feature_compartments=config.feature_compartments,
    )


def restore_donor_from_paths(
    *,
    donor_id: str,
    cells_root: Path,
    redsea_root: Path,
    expression_root: Path,
    geometry_root: Path,
    reference_controls_root: Path,
    registry: MarkerRegistry,
    config: RestoreStageConfig,
    expression_reviews_path: Path | None = None,
) -> RestoreStageResult:
    """Read one donor with column projection, then run RESTORE."""

    donor = str(donor_id)
    pairs = production_pairs(registry)
    review_map: Mapping[tuple[str, str], ExpressionReview] = {}
    if expression_reviews_path is not None:
        review_map, _ = load_expression_reviews(expression_reviews_path)
    donor_frame = restore_input_frame_from_paths(
        donor_id=donor,
        cells_root=cells_root,
        redsea_root=redsea_root,
        expression_root=expression_root,
        geometry_root=geometry_root,
        registry=registry,
        config=config,
        pairs=pairs,
    )
    controls_wide = pd.read_parquet(
        _single_partition(reference_controls_root, donor)
    )
    controls_long = reference_masks_from_wide(
        controls_wide, registry=registry
    )
    validated_controls = validate_reference_masks(
        controls_long,
        donor_id=donor,
        object_ids=donor_frame["object_id"].to_numpy(),
        registry=registry,
    )
    control_arrays = reference_masks_as_arrays(
        validated_controls.masks_long,
        object_ids=donor_frame["object_id"].to_numpy(),
        donor_id=donor,
    )
    return restore_marker_evidence(
        donor_frame,
        donor_id=donor,
        registry=registry,
        config=config,
        expression_reviews=review_map,
        reference_controls=control_arrays,
    )


__all__ = [
    "CALLABLE_STATES",
    "CELL_NEGATIVE",
    "CELL_POSITIVE",
    "CELL_UNAVAILABLE",
    "METHOD_VERSION",
    "RestoreStageConfig",
    "RestoreStageResult",
    "build_restore_input_frame",
    "cell_calls_from_evaluation",
    "production_pairs",
    "restore_donor_from_paths",
    "restore_input_frame_from_paths",
    "restore_marker_evidence",
]
