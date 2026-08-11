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

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from .artifacts import ContractError
from .marker_registry import MarkerRegistry
from .restore_validation import (
    ACCEPTED_PAIRS,
    CANDIDATE_PAIRS,
    METHOD_VERSION as RESTORE_METHOD_VERSION,
    ExpressionReview,
    PairReview,
    PairState,
    PairValidationConfig,
    evaluate_locked_pair,
    load_expression_reviews,
    locked_pair_metrics,
)


# v2: a reviewed absence of target expression calls only cells that were actually measured.  The
# NO_TARGET_EXPRESSION_CONFIRMED branch used to write `negative` over every row regardless of
# ``valid_idx``, so a cell with a missing measurement became a confident biological negative.
METHOD_VERSION = f"restore-marker-evidence-v2+{RESTORE_METHOD_VERSION}"
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

    def __post_init__(self) -> None:
        # PairValidationConfig owns the detailed validation and the full frozen
        # specification. Constructing it here also catches drift immediately.
        self.pair_config()
        if self.on_technical_failure not in {"halt", "unavailable"}:
            raise ValueError(
                "on_technical_failure must be 'halt' or 'unavailable'"
            )

    @classmethod
    def from_pipeline(cls, config: object) -> "RestoreStageConfig":
        compartments = tuple(
            item.strip()
            for item in str(
                getattr(config, "restore_feature_compartments")
            ).split(",")
            if item.strip()
        )
        return cls(
            control_definition=str(
                getattr(config, "restore_control_definition")
            ),
            divisor_statistic=str(
                getattr(config, "restore_divisor_statistic")
            ),
            feature_compartments=compartments,
            on_technical_failure=str(
                getattr(config, "restore_on_technical_failure")
            ),
        )

    def pair_config(self) -> PairValidationConfig:
        return PairValidationConfig(
            control_definition=self.control_definition,
            divisor_statistic=self.divisor_statistic,
            feature_compartments=self.feature_compartments,
        )


@dataclass(frozen=True)
class RestoreStageResult:
    donor_id: str
    model_table: pd.DataFrame
    evidence_wide: pd.DataFrame


def production_pairs(
    registry: MarkerRegistry,
) -> tuple[tuple[str, str], ...]:
    """Return the frozen, accepted RESTORE v10 pair web.

    The registry remains the machine-readable target/reference contract used
    by the rest of the pipeline, but it may not silently drift away from the
    pair web that passed RESTORE's validation gates.
    """

    pairs: list[tuple[str, str]] = []
    for marker in registry.active_markers:
        if marker.reference is None:
            raise ContractError(
                f"active RESTORE target {marker.name!r} has no reference"
            )
        pair = (str(marker.name), str(marker.reference))
        if pair not in CANDIDATE_PAIRS:
            raise ContractError(
                f"active registry pair {pair[0]} <- {pair[1]} was not "
                "validated in the RESTORE candidate specification"
            )
        pairs.append(pair)
    targets = [target for target, _ in pairs]
    if len(targets) != len(set(targets)):
        raise ContractError("RESTORE requires exactly one reference per target")
    observed = {target: reference for target, reference in pairs}
    accepted = {target: reference for target, reference in ACCEPTED_PAIRS}
    if observed != accepted:
        mismatches = {
            target: {
                "accepted": accepted.get(target),
                "registry": observed.get(target),
            }
            for target in sorted(set(observed) | set(accepted))
            if observed.get(target) != accepted.get(target)
        }
        raise ContractError(
            "active registry pairs differ from the accepted RESTORE v10 web: "
            f"{mismatches}"
        )
    # Preserve the reviewed order in model audits and deterministic outputs.
    return tuple(ACCEPTED_PAIRS)


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
    geometry: pd.DataFrame,
    *,
    donor_id: str,
    pairs: Sequence[tuple[str, str]],
    feature_compartments: Sequence[str] = ("Membrane", "Cell"),
    compartment_sample: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the exact column contract consumed by ``evaluate_locked_pair``.

    Canonical raw whole-cell intensities define RESTORE's input floors. REDSEA
    whole-cell intensities define divisors and calls. Optional raw QuPath
    compartment means enrich NNMF group inference. Geometry-ineligible cells
    remain in the output universe but all fit inputs are nulled, ensuring they
    cannot estimate a divisor or receive a marker call.
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
        *marker_order,
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
    geometry_aligned = _align_exact(
        geometry, object_ids, name=f"donor {donor} geometry QC"
    )
    if "analysis_eligible" not in geometry_aligned:
        raise ContractError(
            f"donor {donor}: geometry QC is missing analysis_eligible"
        )
    eligibility = geometry_aligned["analysis_eligible"]
    if eligibility.isna().any() or not eligibility.isin(
        [True, False, 0, 1]
    ).all():
        raise ContractError(
            f"donor {donor}: analysis_eligible must contain explicit booleans"
        )
    analysis_eligible = eligibility.astype(bool).to_numpy()

    sample_aligned: pd.DataFrame | None = None
    sample_coverage = 0
    if compartment_sample is not None:
        sample_keyed = _keyed(
            compartment_sample,
            name=f"donor {donor} RESTORE compartment sample",
        )
        extras = sample_keyed.index.difference(object_ids)
        if len(extras):
            raise ContractError(
                f"donor {donor}: RESTORE compartment sample has "
                f"{len(extras):,} noncanonical cells; examples={extras[:5].tolist()}"
            )
        sample_coverage = int(
            object_ids[analysis_eligible].isin(sample_keyed.index).sum()
        )
        sample_aligned = sample_keyed.reindex(object_ids)

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
    }
    for marker in marker_order:
        raw = pd.to_numeric(canonical[marker], errors="coerce").to_numpy(
            dtype=float
        )
        redsea_column = f"Cell__{marker}"
        if redsea_column not in redsea_aligned:
            raise ContractError(
                f"donor {donor}: REDSEA is missing {redsea_column}"
            )
        corrected = pd.to_numeric(
            redsea_aligned[redsea_column], errors="coerce"
        ).to_numpy(dtype=float)
        # Whole-cell raw means are already canonical in cells and were proven
        # equivalent to the historical sample's Cell__ columns up to float32
        # round-trip. They are therefore available even for markers added after
        # the compartment extraction.
        columns[f"raw__{marker}"] = raw.copy()
        columns[f"redsea__{marker}"] = corrected.copy()

        # NNMF FEATURES COME FROM REDSEA'S COMPARTMENT DECOMPOSITION.
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
            source: pd.Series | None = None
            if name in redsea_aligned:
                source = redsea_aligned[name]
            elif sample_aligned is not None and name in sample_aligned:
                source = sample_aligned[name]
            elif compartment == "Cell":
                columns[name] = raw.copy()
                continue
            elif name in canonical:
                source = canonical[name]
            if source is not None:
                columns[name] = pd.to_numeric(
                    source, errors="coerce"
                ).to_numpy(dtype=float)

    result = pd.DataFrame(columns)
    measurement_columns = [
        column
        for column in result.columns
        if column.startswith(("raw__", "redsea__"))
        or any(
            column.startswith(f"{compartment}__")
            for compartment in feature_compartments
        )
    ]
    result.loc[~analysis_eligible, measurement_columns] = np.nan
    result.attrs["analysis_eligible_n"] = int(analysis_eligible.sum())
    result.attrs["compartment_sample_covered_n"] = sample_coverage
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
    state = np.full(n, CELL_UNAVAILABLE, dtype="<U11")

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
        "method_version": RESTORE_METHOD_VERSION,
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
    evaluator: Callable[..., dict] = evaluate_locked_pair,
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
            evaluation = evaluator(
                donor_frame,
                donor,
                target,
                reference,
                config=pair_config,
                expression_review=reviews.get(
                    (donor, target), ExpressionReview.UNREVIEWED
                ),
                pair_review=PairReview.ACCEPTED,
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
            row = dict(locked_pair_metrics(evaluation))
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
            states = np.full(n, CELL_UNAVAILABLE, dtype="<U11")

        finite_norm = np.isfinite(normalized) & (normalized >= 0)
        log2_ratio = np.full(n, np.nan, dtype=float)
        log2_ratio[finite_norm] = np.log2(
            np.maximum(normalized[finite_norm], np.finfo(float).tiny)
        )
        model_valid = (
            np.isin(states, [CELL_POSITIVE, CELL_NEGATIVE])
            & analysis_eligible
        )
        corrected = pd.to_numeric(
            donor_frame[f"redsea__{target}"], errors="coerce"
        ).to_numpy(dtype=float)
        prefix = f"{target}__"
        evidence[prefix + "corrected_intensity"] = corrected
        evidence[prefix + "state"] = pd.Categorical(
            states, categories=CELL_STATES
        )
        evidence[prefix + "log2_threshold_ratio"] = log2_ratio
        evidence[prefix + "model_valid"] = model_valid
        evidence[prefix + "restore_positive"] = positive
        evidence[prefix + "restore_normalized"] = normalized
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
        row.setdefault("method_version", RESTORE_METHOD_VERSION)
        row["called_negative_n"] = int(
            np.count_nonzero(states == CELL_NEGATIVE)
        )
        row["called_unavailable_n"] = int(
            np.count_nonzero(states == CELL_UNAVAILABLE)
        )
        row["analysis_eligible_n"] = int(analysis_eligible.sum())
        row["compartment_sample_covered_n"] = int(
            donor_frame.attrs.get("compartment_sample_covered_n", 0)
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


def restore_donor_from_paths(
    *,
    donor_id: str,
    cells_root: Path,
    redsea_root: Path,
    geometry_root: Path,
    registry: MarkerRegistry,
    config: RestoreStageConfig,
    compartment_sample_root: Path | None = None,
    expression_reviews_path: Path | None = None,
) -> RestoreStageResult:
    """Read one donor with column projection, then run RESTORE."""

    donor = str(donor_id)
    pairs = production_pairs(registry)
    markers = tuple(dict.fromkeys(marker for pair in pairs for marker in pair))
    cells_columns = [
        "object_id",
        "image",
        "X_centroid",
        "Y_centroid",
        *markers,
        *(
            f"{compartment}__{marker}"
            for compartment in config.feature_compartments
            if compartment != "Cell"
            for marker in markers
        ),
    ]
    cells = _read_columns(
        cells_root, donor, cells_columns, allow_missing=True
    )
    # Bare marker means and the key/context columns are mandatory. Only optional
    # raw compartment columns may be missing from canonical cells.
    required_cells = {
        "object_id",
        "image",
        "X_centroid",
        "Y_centroid",
        *markers,
    }
    missing_cells = sorted(required_cells.difference(cells.columns))
    if missing_cells:
        raise ContractError(
            f"donor {donor}: cells are missing required columns {missing_cells}"
        )
    # Cell__ is mandatory (it is the corrected whole-cell intensity behind
    # `redsea__`); the other compartments are the NNMF feature basis and are
    # requested from the same REDSEA output.
    redsea = _read_columns(
        redsea_root,
        donor,
        [
            "object_id",
            *(f"Cell__{marker}" for marker in markers),
            *(
                f"{compartment}__{marker}"
                for compartment in config.feature_compartments
                if compartment != "Cell"
                for marker in markers
            ),
        ],
        allow_missing=True,
    )
    missing_redsea = sorted(
        {f"Cell__{marker}" for marker in markers}.difference(redsea.columns)
    )
    if missing_redsea:
        raise ContractError(
            f"donor {donor}: REDSEA is missing {missing_redsea}"
        )
    geometry = _read_columns(
        geometry_root, donor, ["object_id", "analysis_eligible"]
    )
    sample = None
    if compartment_sample_root is not None:
        sample = _read_columns(
            compartment_sample_root,
            donor,
            [
                "object_id",
                *(
                    f"{compartment}__{marker}"
                    for compartment in config.feature_compartments
                    if compartment != "Cell"
                    for marker in markers
                ),
            ],
            allow_missing=True,
        )
    review_map: Mapping[tuple[str, str], ExpressionReview] = {}
    if expression_reviews_path is not None:
        review_map, _ = load_expression_reviews(expression_reviews_path)
    donor_frame = build_restore_input_frame(
        cells,
        redsea,
        geometry,
        donor_id=donor,
        pairs=pairs,
        feature_compartments=config.feature_compartments,
        compartment_sample=sample,
    )
    return restore_marker_evidence(
        donor_frame,
        donor_id=donor,
        registry=registry,
        config=config,
        expression_reviews=review_map,
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
    "restore_marker_evidence",
]
