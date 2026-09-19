"""Production RESTORE evaluation from explicit reference-positive seeds.

RETIRED -- NOT A FALLBACK.  Recorded here 2026-09-18 because the docstring below reads as a live
specification and is not one.  RESTORE was retired on 2026-08-14 and is on the never-revive list
(root ``CLAUDE.md``): the RESTORE pair webs, the Gate 2-5 apparatus and the GMM/Otsu gate searches
do not come back as fallbacks.  Nothing in this module produced the typing on disk -- that is the
subpanel classification under ``core/phenocycler/scripts/typing/`` -- and the pipeline stages that
would call it (``controls``, ``calibrate``, ``type``) are blocked by ``NO_VERIFIED_TYPING``, a
latch that has not moved.  Kept as a record of what was implemented and measured.

The historical RESTORE evaluator in :mod:`phenocycler.restore_validation`
remains the reproducible review implementation.  Its raw QuPath Otsu floors
both choose the fitting arms and help choose the final control, however, while
the production pipeline already computes a target-independent, spatially
balanced ``reference_controls`` artifact.  This module is the production
bridge: NNMF still discovers the pair's directional components, but the
negative control is the intersection of the reference component with those
versioned reference-positive seeds.  The call boundary is then a finite-sample
empirical upper-tail order statistic on that control.

The evaluator consumes only generic ``selected__`` call intensities and
source-routed compartment columns. It therefore cannot reach behind the
RESTORE boundary to raw QuPath, REDSEA, or a previous RESTORE output.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .empirical_tail import empirical_tail_probability
from .restore_validation import (
    ExpressionReview,
    PairReview,
    PairState,
    PairValidationConfig,
    assess_pair_state,
    directional_groups,
    divisor_reproducibility_across_seeds,
    fit_nnmf_predict,
    negative_control_statistics,
    robust_minmax_scale,
)

METHOD_VERSION = "restore-reference-seeded-v2-selected-inputs"
FEATURE_TRANSFORMS = frozenset(
    {
        "linear_robust_minmax",
        "log1p_robust_minmax",
        "asinh_robust_minmax",
        "linear_upper_scale",
    }
)
FIT_SAMPLE_STRATEGIES = frozenset({"uniform", "reference_stratified"})


def _fit_indices(
    n: int,
    cap: int,
    seed: int,
    *,
    pool_cap: int | None = None,
) -> np.ndarray:
    """Draw a deterministic uniform fit sample.

    ``pool_cap`` is review-only support for nested sample-size experiments.  A
    fixed outer pool is drawn first, then its prefix is used for the requested
    cap.  Thus 50k, 100k, and 200k experiments with the same ``pool_cap=200k``
    compare sample size rather than three unrelated random draws.  Production
    leaves it unset and therefore preserves the existing sampling exactly.
    """

    if n < 1 or cap < 1:
        raise ValueError("fit population and cap must be positive")
    if pool_cap is not None and pool_cap < cap:
        raise ValueError("fit sample pool cap must be at least the fit cap")
    if n <= cap:
        return np.arange(n, dtype=int)
    rng = np.random.default_rng(seed)
    if pool_cap is None:
        return np.sort(rng.choice(n, size=cap, replace=False))
    pool_n = min(n, pool_cap)
    pool = rng.choice(n, size=pool_n, replace=False)
    return np.sort(pool[:cap])


def _reference_stratified_fit_indices(
    reference_seed: np.ndarray,
    cap: int,
    seed: int,
) -> np.ndarray:
    """Include the fixed reference seeds before drawing the remaining fit rows.

    This review-only design uses no target value, target floor, desired
    frequency, phenotype, or disease label.  Up to half of the fit is reserved
    for the versioned, target-independent reference-positive mask; if that mask
    has fewer rows (the current artifact has 20k), every seed is retained and
    the remaining capacity is filled uniformly from its complement.
    """

    reference_seed = np.asarray(reference_seed, dtype=bool)
    n = len(reference_seed)
    if n < 2 or cap < 2:
        raise ValueError("reference-stratified fitting requires at least two rows")
    if n <= cap:
        return np.arange(n, dtype=int)
    reference_idx = np.flatnonzero(reference_seed)
    complement_idx = np.flatnonzero(~reference_seed)
    if not len(reference_idx) or not len(complement_idx):
        raise ValueError(
            "reference-stratified fitting requires nonempty reference and complement pools"
        )
    reference_n = min(len(reference_idx), cap // 2)
    complement_n = min(len(complement_idx), cap - reference_n)
    remaining = cap - reference_n - complement_n
    if remaining:
        extra_reference = min(remaining, len(reference_idx) - reference_n)
        reference_n += extra_reference
        remaining -= extra_reference
        complement_n += min(remaining, len(complement_idx) - complement_n)
    if reference_n + complement_n != cap:
        raise ValueError("reference-stratified pools cannot fill the requested fit cap")
    rng = np.random.default_rng(seed)
    reference_fit = rng.choice(reference_idx, size=reference_n, replace=False)
    complement_fit = rng.choice(complement_idx, size=complement_n, replace=False)
    return np.sort(np.concatenate((reference_fit, complement_fit)))


def preprocess_nnmf_features(
    features: np.ndarray,
    *,
    transform: str,
    lower_quantile: float,
    upper_quantile: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply one explicit, auditable nonnegative NNMF preprocessing rule.

    Bounds are learned from *all* estimation-eligible rows passed by the caller,
    never from Otsu-selected or frequency-balanced arms.  ``linear_upper_scale``
    is the no-shift comparator; the other options make their nonlinear display
    transform and robust affine scaling explicit.
    """

    values = np.asarray(features, dtype=float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("NNMF features must be a finite two-dimensional matrix")
    if np.any(values < 0):
        raise ValueError("NNMF features must be nonnegative")
    if transform not in FEATURE_TRANSFORMS:
        raise ValueError(
            f"feature_transform must be one of {sorted(FEATURE_TRANSFORMS)}, "
            f"got {transform!r}"
        )

    if transform.startswith("log1p"):
        transformed = np.log1p(values)
    elif transform.startswith("asinh"):
        transformed = np.arcsinh(values)
    else:
        transformed = values

    all_rows = np.arange(len(transformed), dtype=int)
    if transform == "linear_upper_scale":
        upper = np.quantile(transformed, upper_quantile, axis=0)
        if np.any(~np.isfinite(upper)) or np.any(upper <= 0):
            raise ValueError("upper-scale bounds must be finite and positive")
        lower = np.zeros(transformed.shape[1], dtype=float)
        scaled = np.clip(transformed / upper, 0.0, 1.0)
        return scaled, lower, upper

    return robust_minmax_scale(
        transformed,
        all_rows,
        lower_quantile=lower_quantile,
        upper_quantile=upper_quantile,
    )


def _jaccard(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=bool)
    second = np.asarray(second, dtype=bool)
    union = np.count_nonzero(first | second)
    return float(np.count_nonzero(first & second) / union) if union else 1.0


def _finite_mean_ratio(numerator: np.ndarray, denominator: np.ndarray) -> float:
    if not len(numerator) or not len(denominator):
        return float("nan")
    return float(
        np.mean(numerator)
        / max(float(np.mean(denominator)), np.finfo(float).tiny)
    )


def evaluate_reference_seeded_pair(
    donor_frame: pd.DataFrame,
    donor: str,
    target: str,
    reference: str,
    *,
    config: PairValidationConfig,
    marker_alpha: float,
    reference_positive: Sequence[bool] | np.ndarray,
    feature_transform: str = "linear_upper_scale",
    nnmf_max_iter: int = 2_000,
    fit_sample_strategy: str = "uniform",
    fit_sample_seeds: Sequence[int] | None = None,
    fit_sample_pool_cap: int | None = None,
    expression_review: ExpressionReview = ExpressionReview.UNREVIEWED,
    pair_review: PairReview = PairReview.UNREVIEWED,
) -> dict[str, object]:
    """Fit one research pair without presuming scientific acceptance."""

    n = len(donor_frame)
    if n == 0:
        raise ValueError("RESTORE donor frame is empty")
    if not 0 < float(marker_alpha) < 1:
        raise ValueError("marker_alpha must be between 0 and 1")
    if int(nnmf_max_iter) != nnmf_max_iter or nnmf_max_iter < 1:
        raise ValueError("nnmf_max_iter must be a positive integer")
    if fit_sample_strategy not in FIT_SAMPLE_STRATEGIES:
        raise ValueError(
            f"fit_sample_strategy must be one of {sorted(FIT_SAMPLE_STRATEGIES)}"
        )
    if fit_sample_strategy != "uniform" and fit_sample_pool_cap is not None:
        raise ValueError(
            "fit_sample_pool_cap is supported only for uniform nested samples"
        )
    if fit_sample_pool_cap is not None:
        if (
            int(fit_sample_pool_cap) != fit_sample_pool_cap
            or fit_sample_pool_cap < config.fit_cap
        ):
            raise ValueError(
                "fit_sample_pool_cap must be an integer at least as large as fit_cap"
            )
        fit_sample_pool_cap = int(fit_sample_pool_cap)
    if fit_sample_seeds is None:
        sample_seeds = tuple(int(seed) for seed in config.seeds)
        sample_seed_policy = "coupled_to_model_seed"
    else:
        sample_seeds = tuple(int(seed) for seed in fit_sample_seeds)
        if len(sample_seeds) != len(config.seeds):
            raise ValueError("fit_sample_seeds must align one-for-one with config.seeds")
        sample_seed_policy = "explicit"
    reference_seed = np.asarray(reference_positive, dtype=bool)
    if reference_seed.shape != (n,):
        raise ValueError("reference_positive must align with the donor frame")
    for column in ("qc_analysis_eligible", "qc_estimation_eligible"):
        if column not in donor_frame:
            raise ValueError(f"RESTORE donor frame is missing {column}")

    effective_compartments = tuple(
        compartment
        for compartment in config.feature_compartments
        if all(
            f"{compartment}__{marker}" in donor_frame
            for marker in (target, reference)
        )
    )
    if not effective_compartments:
        raise ValueError(
            f"{target} <- {reference}: no declared feature compartment is "
            "present for both markers"
        )
    feature_columns = [
        f"{compartment}__{marker}"
        for compartment in effective_compartments
        for marker in (target, reference)
    ]
    required = [f"selected__{target}", f"selected__{reference}", *feature_columns]
    missing = sorted(set(required).difference(donor_frame.columns))
    if missing:
        raise ValueError(f"{target} <- {reference}: missing columns {missing}")

    analysis_eligible = donor_frame["qc_analysis_eligible"].to_numpy(dtype=bool)
    estimation_eligible = donor_frame["qc_estimation_eligible"].to_numpy(dtype=bool)
    features_all = donor_frame[feature_columns].to_numpy(dtype=float)
    target_all = pd.to_numeric(
        donor_frame[f"selected__{target}"], errors="coerce"
    ).to_numpy(dtype=float)
    reference_all = pd.to_numeric(
        donor_frame[f"selected__{reference}"], errors="coerce"
    ).to_numpy(dtype=float)
    if (
        np.any(np.isfinite(features_all) & (features_all < 0))
        or np.any(np.isfinite(target_all) & (target_all < 0))
        or np.any(np.isfinite(reference_all) & (reference_all < 0))
    ):
        raise ValueError(
            f"{target} <- {reference}: selected RESTORE inputs must be "
            "nonnegative; source-specific censoring belongs at frame construction"
        )

    fit_valid = (
        estimation_eligible
        & np.isfinite(features_all).all(axis=1)
        & (features_all >= 0).all(axis=1)
        & np.isfinite(target_all)
        & (target_all >= 0)
        & np.isfinite(reference_all)
        & (reference_all >= 0)
    )
    valid_idx = np.flatnonzero(fit_valid)
    if len(valid_idx) < max(config.n_components, 2 * config.min_arm_n):
        raise ValueError(
            f"{target} <- {reference}: fewer than "
            f"{max(config.n_components, 2 * config.min_arm_n)} complete "
            "estimation-eligible cells"
        )

    features = features_all[valid_idx]
    target_values = target_all[valid_idx]
    reference_values = reference_all[valid_idx]
    seed_values = reference_seed[valid_idx]
    scaled, scale_lower, scale_upper = preprocess_nnmf_features(
        features,
        transform=feature_transform,
        lower_quantile=config.scale_lower_quantile,
        upper_quantile=config.scale_upper_quantile,
    )

    fits: list[dict[str, object]] = []
    for seed, sample_seed in zip(config.seeds, sample_seeds, strict=True):
        if fit_sample_strategy == "reference_stratified":
            fit_indices = _reference_stratified_fit_indices(
                seed_values,
                config.fit_cap,
                sample_seed,
            )
        else:
            fit_indices = _fit_indices(
                len(scaled),
                config.fit_cap,
                sample_seed,
                pool_cap=fit_sample_pool_cap,
            )
        labels, components, metadata = fit_nnmf_predict(
            scaled,
            n_components=config.n_components,
            seed=int(seed),
            fit_cap=None,
            fit_indices=fit_indices,
            max_iter=int(nnmf_max_iter),
            beta_loss=config.nnmf_beta_loss,
            solver=config.nnmf_solver,
            init=config.nnmf_init,
        )
        direction = directional_groups(target_values, reference_values, labels)
        reference_component = labels == direction.reference_group
        target_component = labels == direction.target_group
        # Exact zero is a left-censoring code at the selected-input boundary,
        # not an observed background value. It may participate in the sparse
        # NNMF representation, but it cannot estimate the empirical tail.
        control = reference_component & seed_values & (target_values > 0)
        target_population = target_component & ~seed_values
        statistics: Mapping[str, object] | None = None
        error = ""
        try:
            statistics = negative_control_statistics(
                target_values[control],
                "empirical_alpha",
                marker_alpha=float(marker_alpha),
            )
        except (ValueError, AssertionError) as exc:
            error = str(exc)
        fits.append(
            {
                "seed": int(seed),
                "sample_seed": int(sample_seed),
                "labels": labels,
                "components": components,
                "metadata": metadata,
                "direction": direction,
                "reference_component": reference_component,
                "target_component": target_component,
                "reference_control": control,
                "target_population": target_population,
                "statistics": statistics,
                "error": error,
                "fit_n": int(len(fit_indices)),
            }
        )

    primary = fits[0]
    primary_control = np.asarray(primary["reference_control"], dtype=bool)
    primary_target = np.asarray(primary["target_population"], dtype=bool)
    primary_direction = primary["direction"]
    primary_stats = primary["statistics"]
    stability_rows = []
    divisors = []
    for fit in fits:
        stats = fit["statistics"]
        divisor = (
            float(stats["divisor"])
            if stats is not None and np.isfinite(float(stats["divisor"]))
            else np.nan
        )
        divisors.append(divisor)
        direction = fit["direction"]
        stability_rows.append(
            {
                "seed": fit["seed"],
                "sample_seed": fit["sample_seed"],
                "converged": bool(fit["metadata"]["converged"]),
                "n_iter": int(fit["metadata"]["n_iter"]),
                "fit_n": int(fit["fit_n"]),
                "reconstruction_error": float(
                    fit["metadata"]["reconstruction_error"]
                ),
                "directional_mean": bool(direction.directional_mean),
                "directional_median": bool(direction.directional_median),
                "group_centers": json.dumps(
                    direction.groups.to_dict(orient="records"),
                    separators=(",", ":"),
                ),
                "component_jaccard": _jaccard(
                    primary["reference_component"], fit["reference_component"]
                ),
                "control_jaccard": _jaccard(
                    primary_control, fit["reference_control"]
                ),
                "control_divisor": divisor,
                "control_error": fit["error"],
            }
        )
    seed_divisor_n, divisor_reproducibility = divisor_reproducibility_across_seeds(
        divisors
    )
    converged = all(bool(row["converged"]) for row in stability_rows)
    stable = bool(
        all(
            bool(row["directional_mean"]) and bool(row["directional_median"])
            for row in stability_rows
        )
        and np.isfinite(divisor_reproducibility)
        and divisor_reproducibility >= config.min_divisor_reproducibility
    )

    target_fold = _finite_mean_ratio(
        target_values[primary_target], target_values[primary_control]
    )
    reference_fold = _finite_mean_ratio(
        reference_values[primary_control], reference_values[primary_target]
    )
    technical_error = str(primary["error"]) or None
    assessment = assess_pair_state(
        reference_n=int(primary_control.sum()),
        target_n=int(primary_target.sum()),
        distinct=primary_direction.reference_group != primary_direction.target_group,
        directional_mean=bool(primary_direction.directional_mean),
        directional_median=bool(primary_direction.directional_median),
        converged=converged,
        stable=stable,
        target_arm_fold=target_fold,
        reference_arm_fold=reference_fold,
        expression_review=expression_review,
        pair_review=pair_review,
        config=config,
        technical_error=technical_error,
    )

    canonical_divisor = (
        float(primary_stats["divisor"])
        if assessment.state is PairState.VALID and primary_stats is not None
        else None
    )
    measurable = (
        analysis_eligible & np.isfinite(target_all) & (target_all > 0)
    )
    call_idx = np.flatnonzero(measurable)
    full_normalized = np.full(n, np.nan, dtype=float)
    full_positive = np.zeros(n, dtype=bool)
    full_tail_probability = np.full(n, np.nan, dtype=float)
    if assessment.state is PairState.VALID:
        if canonical_divisor is None or primary_stats is None:
            raise AssertionError("VALID reference-seeded pair lacks a divisor")
        full_normalized[call_idx] = target_all[call_idx] / canonical_divisor
        full_positive[call_idx] = full_normalized[call_idx] > 1.0
        controls = np.asarray(primary_stats["sorted_clean_log1p"], dtype=float)
        full_tail_probability[call_idx] = empirical_tail_probability(
            np.log1p(target_all[call_idx]), controls
        )
    elif assessment.state is PairState.NO_TARGET_EXPRESSION_CONFIRMED:
        # Measured cells are confident negatives; censored or missing cells are
        # still unavailable and therefore absent from call_idx.
        full_normalized[call_idx] = np.nan

    full_reference_control = np.zeros(n, dtype=bool)
    full_reference_control[valid_idx] = primary_control
    full_target_population = np.zeros(n, dtype=bool)
    full_target_population[valid_idx] = primary_target
    full_labels = np.full(n, -1, dtype=np.int16)
    full_labels[valid_idx] = np.asarray(primary["labels"], dtype=np.int16)
    censored_below = analysis_eligible & np.isfinite(target_all) & (target_all == 0)

    control_stats = {} if primary_stats is None else dict(primary_stats)
    sources = donor_frame.attrs.get("measurement_source_by_marker", {})
    transforms = donor_frame.attrs.get("input_transform_by_marker", {})
    target_source = str(sources.get(target, "selected_expression"))
    reference_source = str(sources.get(reference, "selected_expression"))
    target_transform = str(transforms.get(target, "none"))
    reference_transform = str(transforms.get(reference, "none"))
    model_row = {
        "method_version": METHOD_VERSION,
        "pipeline_stage_method_version": METHOD_VERSION,
        "donor": str(donor),
        "target": target,
        "reference": reference,
        "expression_review": ExpressionReview(expression_review).value,
        "pair_review": PairReview(pair_review).value,
        "state": assessment.state.value,
        "state_reason": assessment.reason,
        "n_input": int(n),
        "n_assigned": int(len(valid_idx)),
        "n_unavailable": int(n - len(call_idx)),
        "analysis_eligible_n": int(analysis_eligible.sum()),
        "estimation_eligible_n": int(estimation_eligible.sum()),
        "censored_below_n": int(censored_below.sum()),
        "target_n": int(primary_target.sum()),
        "target_anchor_n": int(primary_target.sum()),
        "reference_n": int(primary_control.sum()),
        "reference_seed_n": int(reference_seed.sum()),
        "reference_component_n": int(
            np.count_nonzero(primary["reference_component"])
        ),
        "target_component_n": int(np.count_nonzero(primary["target_component"])),
        "target_fold": target_fold,
        "reference_fold": reference_fold,
        "target_low_sbr": bool(
            np.isfinite(target_fold) and target_fold < config.min_target_arm_fold
        ),
        "directional_mean": bool(primary_direction.directional_mean),
        "directional_median": bool(primary_direction.directional_median),
        "converged": converged,
        "stable": stable,
        "seed_divisor_n": int(seed_divisor_n),
        "divisor_reproducibility": divisor_reproducibility,
        "min_control_jaccard": float(
            min(row["control_jaccard"] for row in stability_rows)
        ),
        "canonical_divisor": canonical_divisor,
        "manuscript_divisor": control_stats.get("maximum", np.nan),
        "candidate_full_maximum": control_stats.get("maximum", np.nan),
        "divisor_statistic": "empirical_alpha",
        "marker_alpha": float(marker_alpha),
        "control_fraction_le_1": control_stats.get(
            "control_fraction_le_1", np.nan
        ),
        "control_clean_n": control_stats.get("clean_n", 0),
        "control_contaminated_n": control_stats.get("contaminated_n", 0),
        "control_contamination_fraction": control_stats.get(
            "contamination_fraction", np.nan
        ),
        "control_source": "reference_controls_intersect_nnmf_component",
        "call_intensity_source": "selected_expression",
        "call_intensity_scale": "source-routed nonnegative selected intensity",
        "target_measurement_source": target_source,
        "reference_measurement_source": reference_source,
        "target_input_transform": target_transform,
        "reference_input_transform": reference_transform,
        "nnmf_input_source": "source-routed selected compartments",
        "nnmf_feature_transform": feature_transform,
        "nnmf_fit_basis": "all_estimation_eligible",
        "nnmf_n_components": int(config.n_components),
        "nnmf_beta_loss": config.nnmf_beta_loss,
        "nnmf_solver": config.nnmf_solver,
        "nnmf_init": config.nnmf_init,
        "nnmf_max_iter": int(nnmf_max_iter),
        "nnmf_fit_cap": int(config.fit_cap),
        "nnmf_fit_sample_strategy": fit_sample_strategy,
        "nnmf_fit_sample_seed_policy": sample_seed_policy,
        "nnmf_fit_sample_seeds": json.dumps(sample_seeds, separators=(",", ":")),
        "nnmf_fit_sample_pool_cap": (
            np.nan if fit_sample_pool_cap is None else int(fit_sample_pool_cap)
        ),
        "nnmf_control_membership": "reference_seed_and_component",
        "feature_compartments": "+".join(effective_compartments),
        "feature_compartments_declared": "+".join(config.feature_compartments),
        "feature_scale_lower": json.dumps(scale_lower.tolist(), separators=(",", ":")),
        "feature_scale_upper": json.dumps(scale_upper.tolist(), separators=(",", ":")),
        # Historical fields retained as explicit non-applicability for pilot
        # comparison; no raw floor participates in this evaluator.
        "target_input_floor": np.nan,
        "target_input_floor_method": "not_used",
        "reference_input_floor": np.nan,
        "reference_input_floor_method": "not_used",
    }

    return {
        "method_version": METHOD_VERSION,
        "n_input": n,
        "state": assessment.state.value,
        "state_reason": assessment.reason,
        "valid_idx": valid_idx,
        "call_idx": call_idx,
        "full_lineage_positive": full_positive,
        "full_normalized": full_normalized,
        "full_empirical_tail_probability": full_tail_probability,
        "full_reference_control": full_reference_control,
        "full_target_population": full_target_population,
        "full_labels": full_labels,
        "full_censored_below": censored_below,
        "canonical_divisor": canonical_divisor,
        "reference_n": int(primary_control.sum()),
        "model_row": model_row,
        "stability_rows": stability_rows,
    }


__all__ = [
    "FEATURE_TRANSFORMS",
    "FIT_SAMPLE_STRATEGIES",
    "METHOD_VERSION",
    "evaluate_reference_seeded_pair",
    "preprocess_nnmf_features",
]
