#!/usr/bin/env python3
"""Render review figures from current fundamental-pair development artifacts.

This mutable development tool reads the narrow upstream QuPath projection,
compact survey and rank-overlap outputs, and raw-channel crops written by the
boundary/contact audit.  It has no input-refresh path and never reads RESTORE
outputs, marker calls, assignments, phenotypes, or disease labels.

The roster intentionally mixes raw QuPath cell means with signed QuPath
``Iba1: Mean (norm)``.  Native measurements are consequently drawn on separate,
explicitly labelled axes.  Pair overlap is drawn from donor-local ranks, where
the two measurements need not share a numeric scale.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from decimal import ROUND_CEILING, Decimal
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = REPO.parents[1]
sys.path.insert(0, str(REPO))

from scripts.review.fundamental_pair_survey import (  # noqa: E402
    METHOD_VERSION as SURVEY_METHOD_VERSION,
)
from scripts.review.fundamental_pair_survey import (  # noqa: E402
    ROSTER_PATH,
    DiagnosticRoster,
    load_diagnostic_roster,
)
from scripts.review.fundamental_pair_survey import (  # noqa: E402
    WORKSPACE_MARKER as SURVEY_WORKSPACE_MARKER,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    DEFAULT_OUTPUT_ROOT as DEFAULT_PROJECTION_ROOT,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    METHOD_VERSION as PROJECTION_METHOD_VERSION,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    WORKSPACE_MARKER as PROJECTION_WORKSPACE_MARKER,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    load_qupath_projection_donor,
    validate_qupath_projection,
)
from scripts.review.reference_overlap_audit import (  # noqa: E402
    METHOD_VERSION as OVERLAP_METHOD_VERSION,
)
from scripts.review.reference_overlap_audit import (  # noqa: E402
    WORKSPACE_MARKER as OVERLAP_WORKSPACE_MARKER,
)
from scripts.review.reference_overlap_audit import stable_descending_ranks  # noqa: E402

METHOD_VERSION = "fundamental-pair-review-figures-v8"
CLUSTER_METHOD_VERSION = "fundamental-pair-q70-fit-all-finite-rank-kmeans-review-v4"
CLUSTER_SELECTION_QUANTILE = 0.70
CLUSTER_RANDOM_STATE = 0
CLUSTER_N_INIT = 20
MARGINAL_DEPTH_INCHES = 0.55
MARGINAL_PAD_INCHES = 0.04
MARGINAL_DENSITY_HEADROOM = 1.08
MARGINAL_LINE_WIDTH = 1.25
JOINT_FIGURE_WIDTH_INCHES = 14.2
JOINT_GRID_LEFT = 0.085
JOINT_GRID_RIGHT = 0.89
JOINT_GRID_BOTTOM = 0.06
JOINT_GRID_TOP = 0.84
JOINT_GRID_WSPACE = 0.24
JOINT_GRID_HSPACE = 0.27
DEFAULT_DEVELOPMENT_ROOT = WORKSPACE / "data" / "development"
DEFAULT_SURVEY_ROOT = DEFAULT_DEVELOPMENT_ROOT / "fundamental_pair_diagnostics"
DEFAULT_OVERLAP_ROOT = DEFAULT_DEVELOPMENT_ROOT / "fundamental_pair_overlap"
DEFAULT_CROP_ROOT = DEFAULT_DEVELOPMENT_ROOT / "boundary_contact_diagnostics"
DEFAULT_LABEL = "fundamental_pair_figures"
WORKSPACE_MARKER = ".fundamental_pair_figures_workspace"
BOUNDARY_METHOD_VERSION = "boundary-contact-audit-development-v4"
BOUNDARY_WORKSPACE_MARKER = ".boundary_contact_development"
DEFAULT_CROP_DONOR = "6450"

PAIR_TABLE_NAME = "fundamental_pair_diagnostics.parquet"
FEATURE_TABLE_NAME = "fundamental_input_summary.parquet"
OVERLAP_TABLE_NAME = "reference_overlap_pairs.parquet"
SURVEY_METADATA_NAME = "development_summary.json"
PROJECTION_METADATA_NAME = "projection_summary.json"
CROP_METADATA_NAME = "development_summary.json"
CROP_MANIFEST_TEMPLATE = "crop_manifest/donor_id={donor}/data.parquet"
CONTACT_METRICS_TEMPLATE = "contact_metrics/donor_id={donor}/data.parquet"
FIGURE_STEMS = (
    "fundamental_pair_joint_distributions",
    "rank_tail_overlap",
    "crop_montage_6450",
)
FIGURE_NAMES = tuple(f"{stem}.{suffix}" for stem in FIGURE_STEMS for suffix in ("png", "pdf"))

PAIR_COLUMNS = (
    "method_version",
    "donor",
    "target",
    "target_column",
    "target_source_measurement",
    "target_source_kind",
    "reference",
    "reference_column",
    "reference_source_measurement",
    "reference_source_kind",
    "input_n",
    "measured_n",
    "unmeasured_n",
    "pearson",
    "spearman",
    "target_median",
    "reference_median",
)
FEATURE_COLUMNS = (
    "donor",
    "column",
    "source_measurement",
    "source_kind",
    "input_n",
    "finite_n",
    "missing_n",
    "negative_n",
    "zero_n",
    "minimum",
    "median",
    "q900",
    "q990",
    "maximum",
)
OVERLAP_COLUMNS = (
    "method_version",
    "donor",
    "target",
    "target_column",
    "target_source_measurement",
    "target_source_kind",
    "reference",
    "reference_column",
    "reference_source_measurement",
    "reference_source_kind",
    "finite_pair_n",
    "excluded_nonfinite_n",
    "spearman_score",
    "quantile",
    "target_tail_boundary_score",
    "reference_tail_boundary_score",
    "target_tail_n",
    "reference_tail_n",
    "both_tail_n",
    "both_of_target_tail",
    "both_of_reference_tail",
    "tail_jaccard",
    "tail_independence_lift",
)
CROP_COLUMNS = (
    "donor_id",
    "pair_id",
    "target",
    "reference",
    "target_selection_source_measurement",
    "reference_selection_source_measurement",
    "target_selection_source_kind",
    "reference_selection_source_kind",
    "target_display_marker",
    "reference_display_marker",
    "stratum",
    "stratum_requested_n",
    "stratum_selected_n",
    "crop_rank",
    "object_id",
    "selection_target_cell_value",
    "selection_reference_cell_value",
    "target_percentile_rank",
    "reference_percentile_rank",
    "render_path",
)
CROP_STRATA = ("target_bright", "dual_bright", "reference_bright")


@dataclass(frozen=True)
class FigureInputs:
    """Validated inputs used by the three review figures."""

    pair_diagnostics: pd.DataFrame
    input_summary: pd.DataFrame
    rank_overlap: pd.DataFrame
    crop_manifest: pd.DataFrame
    donors: tuple[str, ...]
    crop_donor: str
    projection_root: Path
    crop_root: Path
    source_files: tuple[Path, ...]
    source_method_versions: Mapping[str, tuple[str, ...]]


@dataclass(frozen=True)
class ClusterMarginals:
    """Bounded summary of a review-only two-cluster population assignment."""

    selection_quantile: float
    finite_n: int
    target_selection_n: int
    reference_selection_n: int
    selection_intersection_n: int
    selected_union_n: int
    fit_n: int
    fit_red_n: int
    fit_blue_n: int
    prediction_n: int
    predicted_outside_fit_n: int
    outside_fit_red_n: int
    outside_fit_blue_n: int
    red_n: int
    blue_n: int
    cluster_overlap_n: int
    red_target_displayed_n: int
    blue_target_displayed_n: int
    red_reference_displayed_n: int
    blue_reference_displayed_n: int
    red_percentile_center: tuple[float, float]
    blue_percentile_center: tuple[float, float]
    opposing: bool
    reference_high_color: str
    red_target_counts: np.ndarray
    blue_target_counts: np.ndarray
    red_reference_counts: np.ndarray
    blue_reference_counts: np.ndarray


@dataclass(frozen=True)
class JointDistribution:
    """Bounded-memory native-scale histogram and conditioned marginals."""

    density: np.ndarray
    finite_n: int
    displayed_n: int
    target_limits: tuple[float, float]
    reference_limits: tuple[float, float]
    cluster_marginals: ClusterMarginals | None = None


def _read_json_object(path: Path) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ValueError(f"missing required survey provenance: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON survey provenance: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"survey provenance must be a JSON object: {path}")
    return payload


def _require_exact_marker(path: Path, expected: str, *, description: str) -> Path:
    try:
        value = path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        raise ValueError(f"missing required {description} workspace marker: {path}") from exc
    if value != expected:
        raise ValueError(
            f"{description} workspace is stale; expected marker {expected!r}, found {value!r}"
        )
    return path.resolve()


def _read_columns(path: Path, columns: Sequence[str], *, description: str) -> pd.DataFrame:
    if not path.is_file():
        raise ValueError(f"missing {description}: {path}")
    try:
        frame = pd.read_parquet(path, columns=list(columns))
    except (KeyError, OSError, ValueError) as exc:
        raise ValueError(
            f"could not read the required {description} columns from {path}: {exc}"
        ) from exc
    if frame.empty:
        raise ValueError(f"{description} is empty: {path}")
    frame = frame.copy()
    frame["donor"] = frame["donor"].astype(str)
    return frame


def _pair_key_columns() -> tuple[str, ...]:
    return (
        "target",
        "target_column",
        "target_source_measurement",
        "target_source_kind",
        "reference",
        "reference_column",
        "reference_source_measurement",
        "reference_source_kind",
    )


def _expected_pair_keys(roster: DiagnosticRoster) -> set[tuple[str, ...]]:
    return {
        (
            pair.target,
            pair.target_column,
            pair.target_source_measurement,
            pair.target_source_kind,
            pair.reference,
            pair.reference_column,
            pair.reference_source_measurement,
            pair.reference_source_kind,
        )
        for pair in roster.pairs
    }


def _validate_method_versions(
    frame: pd.DataFrame,
    *,
    expected: str,
    description: str,
) -> tuple[str, ...]:
    versions = tuple(sorted(set(frame["method_version"].astype(str))))
    if versions != (expected,):
        raise ValueError(
            f"{description} is stale or is not the current upstream QuPath review artifact; "
            f"expected method {expected!r}, found {list(versions)!r}"
        )
    return versions


def _validate_pair_provenance(
    frame: pd.DataFrame,
    *,
    roster: DiagnosticRoster,
    description: str,
) -> None:
    columns = _pair_key_columns()
    actual = set(frame.loc[:, columns].itertuples(index=False, name=None))
    expected = _expected_pair_keys(roster)
    if actual != expected:
        raise ValueError(
            f"{description} pair/source contract does not match the diagnostic roster; "
            "RESTORE aliases and substituted measurements are not accepted"
        )


def _validate_source_provenance(frame: pd.DataFrame, roster: DiagnosticRoster) -> None:
    actual = set(
        frame.loc[:, ["column", "source_measurement", "source_kind"]].itertuples(
            index=False, name=None
        )
    )
    expected = {
        (column, measurement, kind)
        for column, (measurement, kind) in roster.source_contract.items()
    }
    if actual != expected:
        raise ValueError(
            "input-summary source contract does not match the diagnostic roster; "
            "RESTORE aliases and substituted measurements are not accepted"
        )


def _require_unique(frame: pd.DataFrame, columns: Sequence[str], description: str) -> None:
    duplicated = frame.duplicated(list(columns), keep=False)
    if duplicated.any():
        raise ValueError(f"{description} has duplicate rows for {list(columns)}")


def _require_numeric(
    frame: pd.DataFrame,
    columns: Sequence[str],
    *,
    description: str,
    finite: bool = True,
) -> None:
    for column in columns:
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if finite and not np.isfinite(values).all():
            raise ValueError(f"{description} column {column!r} must contain finite numeric values")
        frame[column] = values


def _resolve_donors(
    pair: pd.DataFrame,
    features: pd.DataFrame,
    overlap: pd.DataFrame,
    requested: Sequence[str] | None,
) -> tuple[str, ...]:
    available = [set(table["donor"]) for table in (pair, features, overlap)]
    if requested:
        donors = tuple(map(str, requested))
        if len(set(donors)) != len(donors):
            raise ValueError("requested donors must be unique")
        missing = {
            name: sorted(set(donors).difference(values))
            for name, values in zip(
                ("survey", "input summary", "rank overlap"), available, strict=True
            )
            if set(donors).difference(values)
        }
        if missing:
            raise ValueError(f"requested donors are missing from review artifacts: {missing}")
        return donors
    if not (available[0] == available[1] == available[2]):
        raise ValueError(
            "survey, input-summary, and rank-overlap donor sets differ; "
            "select an explicit complete subset with --donors"
        )
    return tuple(sorted(available[0]))


def _validate_complete_grid(
    frame: pd.DataFrame,
    *,
    donors: Sequence[str],
    pair_columns: Sequence[str],
    expected_pairs: Sequence[str],
    description: str,
) -> None:
    actual = set(frame.loc[:, pair_columns].itertuples(index=False, name=None))
    expected = {(donor, pair) for donor in donors for pair in expected_pairs}
    if actual != expected:
        raise ValueError(f"{description} does not contain one complete donor/pair grid")


def _validate_tables(
    pair: pd.DataFrame,
    features: pd.DataFrame,
    overlap: pd.DataFrame,
    *,
    roster: DiagnosticRoster,
    requested_donors: Sequence[str] | None,
) -> tuple[tuple[str, ...], dict[str, tuple[str, ...]]]:
    versions = {
        "survey": _validate_method_versions(
            pair, expected=SURVEY_METHOD_VERSION, description="pair survey"
        ),
        "rank_overlap": _validate_method_versions(
            overlap, expected=OVERLAP_METHOD_VERSION, description="rank-overlap audit"
        ),
    }
    _validate_pair_provenance(pair, roster=roster, description="pair survey")
    _validate_pair_provenance(overlap, roster=roster, description="rank-overlap audit")
    _validate_source_provenance(features, roster)

    donors = _resolve_donors(pair, features, overlap, requested_donors)
    selected = set(donors)
    pair.drop(pair.index[~pair["donor"].isin(selected)], inplace=True)
    features.drop(features.index[~features["donor"].isin(selected)], inplace=True)
    overlap.drop(overlap.index[~overlap["donor"].isin(selected)], inplace=True)

    _require_unique(pair, ("donor", "target_column", "reference_column"), "pair survey")
    _require_unique(features, ("donor", "column"), "input summary")
    _require_unique(
        overlap,
        ("donor", "target_column", "reference_column", "quantile"),
        "rank-overlap audit",
    )
    targets = [item.target_column for item in roster.pairs]
    _validate_complete_grid(
        pair,
        donors=donors,
        pair_columns=("donor", "target_column"),
        expected_pairs=targets,
        description="pair survey",
    )
    _validate_complete_grid(
        features,
        donors=donors,
        pair_columns=("donor", "column"),
        expected_pairs=list(roster.source_contract),
        description="input summary",
    )

    expected_overlap = {
        (donor, target, float(quantile))
        for donor in donors
        for target in targets
        for quantile in roster.required_tail_quantiles
    }
    actual_overlap = set(
        overlap.loc[:, ["donor", "target_column", "quantile"]].itertuples(index=False, name=None)
    )
    if actual_overlap != expected_overlap:
        raise ValueError("rank-overlap audit does not contain the roster's complete tail grid")

    _require_numeric(
        pair,
        (
            "input_n",
            "measured_n",
            "unmeasured_n",
            "pearson",
            "spearman",
            "target_median",
            "reference_median",
        ),
        description="pair survey",
    )
    _require_numeric(
        features,
        (
            "input_n",
            "finite_n",
            "missing_n",
            "negative_n",
            "zero_n",
            "minimum",
            "median",
            "q900",
            "q990",
            "maximum",
        ),
        description="input summary",
    )
    _require_numeric(
        overlap,
        (
            "finite_pair_n",
            "excluded_nonfinite_n",
            "spearman_score",
            "quantile",
            "target_tail_boundary_score",
            "reference_tail_boundary_score",
            "target_tail_n",
            "reference_tail_n",
            "both_tail_n",
            "both_of_target_tail",
            "both_of_reference_tail",
            "tail_jaccard",
            "tail_independence_lift",
        ),
        description="rank-overlap audit",
    )

    if ((pair[["pearson", "spearman"]].abs() > 1 + 1e-10).any()).any():
        raise ValueError("pair correlations must lie between -1 and 1")
    quantiles = features[["minimum", "median", "q900", "q990", "maximum"]].to_numpy()
    if (np.diff(quantiles, axis=1) < -1e-10).any():
        raise ValueError("input-summary native quantiles must be nondecreasing")
    for count in ("input_n", "finite_n", "missing_n", "negative_n", "zero_n"):
        if (features[count] < 0).any():
            raise ValueError(f"input-summary count {count!r} cannot be negative")
    if not np.array_equal(
        overlap["target_tail_n"].to_numpy(), overlap["reference_tail_n"].to_numpy()
    ):
        raise ValueError("rank-overlap figures require equal-size exact target/reference tails")
    if (overlap["both_tail_n"] < 0).any() or (
        overlap["both_tail_n"] > overlap["target_tail_n"]
    ).any():
        raise ValueError("rank-overlap intersection counts are inconsistent with tail sizes")
    bounded = overlap[["both_of_target_tail", "both_of_reference_tail", "tail_jaccard"]]
    if ((bounded < -1e-10) | (bounded > 1 + 1e-10)).any().any():
        raise ValueError("rank-overlap proportions must lie between zero and one")
    if (overlap["tail_independence_lift"] < 0).any():
        raise ValueError("rank-overlap independence lift cannot be negative")
    if not np.allclose(
        overlap["both_of_target_tail"],
        overlap["both_of_reference_tail"],
        rtol=1e-7,
        atol=1e-10,
    ):
        raise ValueError("equal-size rank tails must report symmetric overlap fractions")
    return donors, versions


def _same_resolved_path(value: object, expected: Path) -> bool:
    if not isinstance(value, str) or not value:
        return False
    return Path(value).expanduser().resolve() == expected


def _validate_survey_metadata(
    metadata: Mapping[str, object],
    *,
    projection_root: Path,
    roster: DiagnosticRoster,
) -> None:
    expected = {
        "status": "complete",
        "method_version": SURVEY_METHOD_VERSION,
        "method_status": "development_mutable",
        "roster_version": roster.roster_version,
        "restore_outputs_read": False,
        "creates_cell_calls": False,
        "production_outputs_written": False,
    }
    drift = [key for key, value in expected.items() if metadata.get(key) != value]
    if drift:
        raise ValueError(f"pair-survey provenance is stale or unsafe: {sorted(drift)}")
    if not _same_resolved_path(metadata.get("input_root"), projection_root):
        raise ValueError("pair-survey provenance does not point to the selected QuPath projection")


def _crop_pair_key(pair) -> tuple[str, ...]:
    return (
        f"{pair.target_column}:{pair.reference_column}",
        pair.target_column,
        pair.reference_column,
        pair.target_source_measurement,
        pair.reference_source_measurement,
        pair.target_source_kind,
        pair.reference_source_kind,
    )


def _validate_crop_metadata(
    metadata: Mapping[str, object],
    *,
    crop_donor: str,
    projection_root: Path,
    roster: DiagnosticRoster,
) -> None:
    if metadata.get("method") != BOUNDARY_METHOD_VERSION:
        raise ValueError(
            "boundary-crop artifact is stale; "
            f"expected {BOUNDARY_METHOD_VERSION!r}, found {metadata.get('method')!r}"
        )
    if metadata.get("review_only") is not True:
        raise ValueError("boundary-crop provenance must explicitly attest review_only=true")
    selection = metadata.get("selection")
    if not isinstance(selection, Mapping) or selection.get("restore_outputs_read") is not False:
        raise ValueError(
            "boundary-crop provenance must explicitly attest restore_outputs_read=false"
        )
    display = metadata.get("display")
    if (
        not isinstance(display, Mapping)
        or display.get("source") != "raw qptiff channels"
        or display.get("selection_values_used_for_display") is not False
    ):
        raise ValueError("boundary-crop display provenance must declare independent raw channels")
    if not _same_resolved_path(metadata.get("qupath_input_projection_root"), projection_root):
        raise ValueError(
            "boundary-crop provenance does not point to the selected QuPath projection"
        )

    expected_sources = {
        column: {"measurement": measurement, "kind": kind}
        for column, (measurement, kind) in roster.source_contract.items()
    }
    measurements = metadata.get("selection_source_measurements")
    kinds = metadata.get("selection_source_kinds")
    if not isinstance(measurements, Mapping) or not isinstance(kinds, Mapping):
        raise ValueError("boundary-crop provenance lacks exact selection sources")
    for column, source in expected_sources.items():
        if measurements.get(column) != source["measurement"] or kinds.get(column) != source["kind"]:
            raise ValueError(
                "boundary-crop selection sources do not match the roster; RESTORE aliases are refused"
            )
    actual_pairs = metadata.get("pairs")
    expected_pairs = [
        {"target": pair.target_column, "reference": pair.reference_column} for pair in roster.pairs
    ]
    if actual_pairs != expected_pairs:
        raise ValueError("boundary-crop roster does not match the current diagnostic roster")
    donor_rows = metadata.get("donors")
    if not isinstance(donor_rows, list) or not any(
        isinstance(row, Mapping)
        and str(row.get("donor_id")) == crop_donor
        and int(row.get("crop_n", -1)) > 0
        for row in donor_rows
    ):
        raise ValueError(f"boundary-crop provenance has no complete donor {crop_donor}")


def _load_crop_manifest(
    crop_root: Path,
    *,
    crop_donor: str,
    roster: DiagnosticRoster,
) -> tuple[pd.DataFrame, Path, Path]:
    path = crop_root / CROP_MANIFEST_TEMPLATE.format(donor=crop_donor)
    if not path.is_file():
        raise ValueError(f"missing boundary-crop manifest: {path}")
    try:
        frame = pd.read_parquet(path, columns=list(CROP_COLUMNS))
    except (KeyError, OSError, ValueError) as exc:
        raise ValueError(f"could not read current boundary-crop manifest {path}: {exc}") from exc
    if frame.empty:
        raise ValueError(f"boundary-crop manifest is empty: {path}")
    frame = frame.copy()
    frame["donor_id"] = frame["donor_id"].astype(str)
    if set(frame["donor_id"]) != {crop_donor}:
        raise ValueError("boundary-crop manifest contains the wrong internal donor")

    provenance_columns = (
        "pair_id",
        "target",
        "reference",
        "target_selection_source_measurement",
        "reference_selection_source_measurement",
        "target_selection_source_kind",
        "reference_selection_source_kind",
    )
    actual = set(frame.loc[:, provenance_columns].itertuples(index=False, name=None))
    expected = {_crop_pair_key(pair) for pair in roster.pairs}
    if actual != expected:
        raise ValueError(
            "boundary-crop pair/source contract does not match the roster; RESTORE aliases are refused"
        )
    _require_unique(frame, ("pair_id", "stratum", "crop_rank"), "boundary-crop manifest")
    _require_numeric(
        frame,
        (
            "stratum_requested_n",
            "stratum_selected_n",
            "crop_rank",
            "selection_target_cell_value",
            "selection_reference_cell_value",
            "target_percentile_rank",
            "reference_percentile_rank",
        ),
        description="boundary-crop manifest",
    )
    if set(frame["stratum"]) != set(CROP_STRATA):
        raise ValueError("boundary-crop manifest has an unexpected selection stratum")
    requested = set(frame["stratum_requested_n"].astype(int))
    selected = set(frame["stratum_selected_n"].astype(int))
    if len(requested) != 1 or requested != selected or next(iter(requested)) <= 0:
        raise ValueError("boundary-crop manifest has incomplete stratum selection")
    crop_n = next(iter(requested))
    expected_grid = {
        (f"{pair.target_column}:{pair.reference_column}", stratum, rank)
        for pair in roster.pairs
        for stratum in CROP_STRATA
        for rank in range(1, crop_n + 1)
    }
    actual_grid = set(
        frame.loc[:, ["pair_id", "stratum", "crop_rank"]].itertuples(index=False, name=None)
    )
    if actual_grid != expected_grid:
        raise ValueError("boundary-crop manifest does not contain a complete pair/stratum grid")
    percentiles = frame[["target_percentile_rank", "reference_percentile_rank"]]
    if ((percentiles < 0) | (percentiles > 1)).any().any():
        raise ValueError("boundary-crop percentile ranks must lie between zero and one")

    render_root = (crop_root / "crops" / f"donor_id={crop_donor}").resolve()
    resolved_paths: list[str] = []
    for value in frame["render_path"]:
        render = Path(str(value)).expanduser().resolve()
        try:
            render.relative_to(render_root)
        except ValueError as exc:
            raise ValueError(
                f"boundary crop is outside its declared raw-crop workspace: {render}"
            ) from exc
        if render.suffix.lower() != ".png" or not render.is_file():
            raise ValueError(f"missing rendered boundary crop: {render}")
        resolved_paths.append(str(render))
    frame["render_path"] = resolved_paths
    contact_path = crop_root / CONTACT_METRICS_TEMPLATE.format(donor=crop_donor)
    if not contact_path.is_file():
        raise ValueError(f"missing boundary-contact metrics: {contact_path}")
    try:
        contacts = pd.read_parquet(
            contact_path,
            columns=["object_id", "contact_neighbor_n", "contact_boundary_px"],
        )
    except (KeyError, OSError, ValueError) as exc:
        raise ValueError(f"could not read boundary-contact metrics {contact_path}: {exc}") from exc
    contacts = contacts.copy()
    contacts["object_id"] = contacts["object_id"].astype(str)
    _require_unique(contacts, ("object_id",), "boundary-contact metrics")
    _require_numeric(
        contacts,
        ("contact_neighbor_n", "contact_boundary_px"),
        description="boundary-contact metrics",
    )
    frame = frame.merge(contacts, on="object_id", how="left", validate="many_to_one")
    if frame[["contact_neighbor_n", "contact_boundary_px"]].isna().any().any():
        raise ValueError("selected boundary crops are absent from contact metrics")
    return frame, path.resolve(), contact_path.resolve()


def load_figure_inputs(
    projection_root: str | Path,
    survey_root: str | Path,
    overlap_root: str | Path,
    crop_root: str | Path,
    *,
    roster: DiagnosticRoster,
    donors: Sequence[str] | None = None,
    crop_donor: str = DEFAULT_CROP_DONOR,
) -> FigureInputs:
    """Version-gate the current upstream summaries, projection, and raw crops."""

    projection = Path(projection_root).expanduser().resolve()
    survey = Path(survey_root).expanduser().resolve()
    rank = Path(overlap_root).expanduser().resolve()
    crops = Path(crop_root).expanduser().resolve()
    projection_marker_path = _require_exact_marker(
        projection / PROJECTION_WORKSPACE_MARKER,
        PROJECTION_METHOD_VERSION,
        description="QuPath projection",
    )
    survey_marker_path = _require_exact_marker(
        survey / SURVEY_WORKSPACE_MARKER,
        SURVEY_METHOD_VERSION,
        description="pair survey",
    )
    overlap_marker_path = _require_exact_marker(
        rank / OVERLAP_WORKSPACE_MARKER,
        OVERLAP_METHOD_VERSION,
        description="rank overlap",
    )
    metadata_path = survey / SURVEY_METADATA_NAME
    metadata = _read_json_object(metadata_path)
    _validate_survey_metadata(metadata, projection_root=projection, roster=roster)

    pair_path = survey / PAIR_TABLE_NAME
    feature_path = survey / FEATURE_TABLE_NAME
    overlap_path = rank / OVERLAP_TABLE_NAME
    pair = _read_columns(pair_path, PAIR_COLUMNS, description="pair-survey table")
    features = _read_columns(feature_path, FEATURE_COLUMNS, description="input-summary table")
    overlap = _read_columns(overlap_path, OVERLAP_COLUMNS, description="rank-overlap table")
    selected, versions = _validate_tables(
        pair,
        features,
        overlap,
        roster=roster,
        requested_donors=donors,
    )
    crop_donor_id = str(crop_donor)
    projection_summary = validate_qupath_projection(
        projection, tuple(dict.fromkeys((*selected, crop_donor_id)))
    )
    if projection_summary.get("method_version") != PROJECTION_METHOD_VERSION:
        raise ValueError("QuPath projection is not the exact current method version")
    boundary_metadata_path = crops / CROP_METADATA_NAME
    boundary_marker_path = crops / BOUNDARY_WORKSPACE_MARKER
    boundary_marker = _read_json_object(boundary_marker_path)
    if boundary_marker != {"tool": "boundary_contact_audit", "status": "complete"}:
        raise ValueError("boundary-crop workspace marker is absent, incomplete, or invalid")
    boundary_metadata = _read_json_object(boundary_metadata_path)
    _validate_crop_metadata(
        boundary_metadata,
        crop_donor=crop_donor_id,
        projection_root=projection,
        roster=roster,
    )
    crop_manifest, crop_manifest_path, contact_path = _load_crop_manifest(
        crops,
        crop_donor=crop_donor_id,
        roster=roster,
    )
    projection_metadata_path = (projection / PROJECTION_METADATA_NAME).resolve()
    return FigureInputs(
        pair_diagnostics=pair.reset_index(drop=True),
        input_summary=features.reset_index(drop=True),
        rank_overlap=overlap.reset_index(drop=True),
        crop_manifest=crop_manifest.reset_index(drop=True),
        donors=selected,
        crop_donor=crop_donor_id,
        projection_root=projection,
        crop_root=crops,
        source_files=(
            projection_marker_path,
            projection_metadata_path,
            survey_marker_path,
            pair_path.resolve(),
            feature_path.resolve(),
            overlap_marker_path,
            overlap_path.resolve(),
            boundary_marker_path.resolve(),
            boundary_metadata_path.resolve(),
            crop_manifest_path,
            contact_path,
        ),
        source_method_versions={
            "projection": (PROJECTION_METHOD_VERSION,),
            **versions,
            "boundary_crops": (BOUNDARY_METHOD_VERSION,),
        },
    )


def prepare_development_workspace(root: str | Path, label: str) -> Path:
    """Replace exactly one figure-tool-owned mutable workspace."""

    if not label or Path(label).name != label or label in {".", ".."}:
        raise ValueError("label must be one nonempty path component")
    development_root = Path(root).expanduser().resolve()
    dangerous = {Path("/"), Path.home().resolve(), WORKSPACE.resolve(), REPO.resolve()}
    if development_root in dangerous:
        raise ValueError("development-root must be a dedicated subdirectory")
    development_root.mkdir(parents=True, exist_ok=True)
    output = (development_root / label).resolve()
    if output.parent != development_root:
        raise ValueError("development output must be directly below development-root")
    if output.exists():
        if not (output / WORKSPACE_MARKER).is_file():
            raise ValueError(f"refusing to replace a workspace not owned by this tool: {output}")
        shutil.rmtree(output)
    output.mkdir()
    (output / WORKSPACE_MARKER).write_text(METHOD_VERSION + "\n", encoding="utf-8")
    return output


@lru_cache(maxsize=1)
def _plot_modules():
    os.environ.setdefault(
        "MPLCONFIGDIR",
        str(Path(tempfile.gettempdir()) / "phenocycler-matplotlib-cache"),
    )
    try:
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise RuntimeError(
            "figure rendering requires matplotlib; install the notebook extra with "
            "python -m pip install -e '.[notebook]'"
        ) from exc
    return plt, Line2D


def _pair_title(pair) -> str:
    target_scale = (
        "signed Mean (norm)"
        if pair.target_source_kind == "qupath_acquisition_normalized_cell_mean"
        else "raw cell mean"
    )
    return f"{pair.target} vs {pair.reference}\n{target_scale} vs raw cell mean"


def _donor_style(donors: Sequence[str], plt) -> tuple[dict[str, object], dict[str, str]]:
    colors = plt.get_cmap("tab10").colors
    markers = ("o", "s", "^", "D", "P", "X", "v", "<", ">")
    return (
        {donor: colors[index % len(colors)] for index, donor in enumerate(donors)},
        {donor: markers[index % len(markers)] for index, donor in enumerate(donors)},
    )


def native_joint_histogram(
    *,
    target_values: np.ndarray,
    reference_values: np.ndarray,
    target_limits: tuple[float, float],
    reference_limits: tuple[float, float],
    bins: int,
) -> JointDistribution:
    """Bin one native-scale pair after filtering to an explicit display window.

    The full finite-pair count remains the denominator, so the histogram sum is
    exactly the displayed fraction. Values beyond the display window are
    excluded, never clamped into an edge bin. Only the bounded-size histogram
    survives this function.
    """

    if bins < 20:
        raise ValueError("joint-distribution bins must be at least 20")
    target = np.asarray(target_values, dtype=float)
    reference = np.asarray(reference_values, dtype=float)
    if target.ndim != 1 or reference.shape != target.shape:
        raise ValueError("joint-distribution inputs must be aligned one-dimensional arrays")
    if not all(np.isfinite((*target_limits, *reference_limits))):
        raise ValueError("joint-distribution display limits must be finite")
    if target_limits[0] >= target_limits[1] or reference_limits[0] >= reference_limits[1]:
        raise ValueError("joint-distribution display limits must be increasing")

    finite = np.isfinite(target) & np.isfinite(reference)
    target = target[finite]
    reference = reference[finite]
    finite_n = len(target)
    if finite_n < 2:
        raise ValueError("joint distribution requires at least two finite pairs")
    in_view = (
        (target >= target_limits[0])
        & (target <= target_limits[1])
        & (reference >= reference_limits[0])
        & (reference <= reference_limits[1])
    )
    displayed_n = int(in_view.sum())
    if displayed_n < 1:
        raise ValueError("joint-distribution display window contains no finite pairs")
    target_edges = np.linspace(*target_limits, bins + 1)
    reference_edges = np.linspace(*reference_limits, bins + 1)
    histogram, _, _ = np.histogram2d(
        reference[in_view],
        target[in_view],
        bins=(reference_edges, target_edges),
    )
    density = histogram.astype(float, copy=False) / finite_n
    return JointDistribution(
        density=density,
        finite_n=finite_n,
        displayed_n=displayed_n,
        target_limits=target_limits,
        reference_limits=reference_limits,
    )


def _exact_upper_tail_count(n: int, quantile: float) -> int:
    """Return the exact rank-tail size used by the overlap audit."""

    if n < 1 or not 0 < quantile < 1:
        raise ValueError("rank tail requires rows and a quantile between zero and one")
    count = int(
        ((Decimal("1") - Decimal(str(quantile))) * Decimal(n)).to_integral_value(
            rounding=ROUND_CEILING
        )
    )
    return min(n, max(1, count))


def _bright_percentile_ranks(
    values: np.ndarray, object_ids: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return stable descending ranks and continuous 0–1 bright percentiles."""

    ranks = stable_descending_ranks(values, object_ids)
    if len(ranks) < 2:
        raise ValueError("bright percentile ranks require at least two rows")
    percentiles = (len(ranks) - ranks).astype(float) / (len(ranks) - 1)
    return ranks, percentiles


def review_kmeans_cluster_marginals(
    *,
    object_ids: np.ndarray,
    target_values: np.ndarray,
    reference_values: np.ndarray,
    target_key: str,
    reference_key: str,
    target_limits: tuple[float, float],
    reference_limits: tuple[float, float],
    bins: int,
    selection_quantile: float = CLUSTER_SELECTION_QUANTILE,
    random_state: int = CLUSTER_RANDOM_STATE,
    n_init: int = CLUSTER_N_INIT,
) -> ClusterMarginals:
    """Fit a review-only two-cluster assignment on an exact upper-rank union.

    Both marker values become continuous donor-local bright-percentile ranks
    using the same stable descending-score/object-ID contract as the overlap
    audit. Exact upper-rank membership at ``selection_quantile`` removes only
    the double-low population before fitting. KMeans fits every selected-union
    cell in canonical marker-column and object-ID order, then predicts every
    finite donor/pair cell into exactly one of the two fitted clusters.

    Red is the cluster with the higher target-percentile center. Blue is the
    higher-reference cluster only when it is distinct; otherwise blue is the
    other cluster and ``opposing`` is false. This is an unverified visualization
    assignment, never a biological call or inferred threshold.
    """

    if bins < 20:
        raise ValueError("cluster-marginal bins must be at least 20")
    if not 0 < selection_quantile < 1:
        raise ValueError("cluster selection quantile must be between zero and one")
    if n_init < 1:
        raise ValueError("KMeans n_init must be positive")
    identifiers = np.asarray(object_ids)
    target = np.asarray(target_values, dtype=float)
    reference = np.asarray(reference_values, dtype=float)
    if (
        identifiers.ndim != 1
        or target.shape != identifiers.shape
        or reference.shape != identifiers.shape
    ):
        raise ValueError("cluster-marginal inputs must be aligned one-dimensional arrays")
    if not all(np.isfinite((*target_limits, *reference_limits))):
        raise ValueError("cluster-marginal display limits must be finite")
    if target_limits[0] >= target_limits[1] or reference_limits[0] >= reference_limits[1]:
        raise ValueError("cluster-marginal display limits must be increasing")
    target_name = str(target_key)
    reference_name = str(reference_key)
    if not target_name or not reference_name or target_name == reference_name:
        raise ValueError("cluster marginals require two distinct nonempty marker keys")

    finite = np.isfinite(target) & np.isfinite(reference)
    identifiers = identifiers[finite]
    target = target[finite]
    reference = reference[finite]
    if len(target) < 2:
        raise ValueError("cluster marginals require at least two finite pairs")
    target_ranks, target_percentiles = _bright_percentile_ranks(target, identifiers)
    reference_ranks, reference_percentiles = _bright_percentile_ranks(reference, identifiers)
    tail_n = _exact_upper_tail_count(len(target), selection_quantile)
    target_selected = target_ranks <= tail_n
    reference_selected = reference_ranks <= tail_n
    selected = target_selected | reference_selected
    selected_n = int(selected.sum())
    if selected_n < 2:
        raise ValueError("selected upper-rank union contains fewer than two cells")
    selected_ids = identifiers[selected]
    directed_features = np.column_stack((target_percentiles, reference_percentiles))
    target_first = target_name < reference_name
    canonical_features = directed_features if target_first else directed_features[:, ::-1]
    selected_features = canonical_features[selected]
    fit_order = np.argsort(selected_ids, kind="stable")
    if not np.any(np.any(selected_features != selected_features[0], axis=1)):
        raise ValueError("selected upper-rank union lacks two distinct fit rows for KMeans")

    from sklearn.cluster import KMeans

    model = KMeans(
        n_clusters=2,
        init="k-means++",
        n_init=n_init,
        max_iter=300,
        algorithm="lloyd",
        random_state=random_state,
    ).fit(selected_features[fit_order])
    labels = model.predict(canonical_features)
    fit_labels = labels[selected]
    if not np.any(fit_labels == 0) or not np.any(fit_labels == 1):
        raise ValueError("KMeans did not assign two nonempty clusters")
    canonical_centers = np.asarray(model.cluster_centers_, dtype=float)
    centers = canonical_centers if target_first else canonical_centers[:, ::-1]
    target_winner = int(np.argmax(centers[:, 0]))
    reference_winner = int(np.argmax(centers[:, 1]))
    opposing = target_winner != reference_winner
    red_cluster = target_winner
    blue_cluster = reference_winner if opposing else 1 - red_cluster
    red = labels == red_cluster
    blue = labels == blue_cluster
    if np.any(red & blue) or not np.all(red | blue):  # pragma: no cover - KMeans invariant
        raise AssertionError("two-cluster assignment must be disjoint and exhaustive")

    fit_red_n = int(np.count_nonzero(selected & red))
    fit_blue_n = int(np.count_nonzero(selected & blue))
    outside_fit = ~selected
    outside_fit_red_n = int(np.count_nonzero(outside_fit & red))
    outside_fit_blue_n = int(np.count_nonzero(outside_fit & blue))
    target_edges = np.linspace(*target_limits, bins + 1)
    reference_edges = np.linspace(*reference_limits, bins + 1)

    def counts(values: np.ndarray, group: np.ndarray, edges: np.ndarray) -> np.ndarray:
        return np.histogram(values[group], bins=edges)[0].astype(np.int64, copy=False)

    red_target_counts = counts(target, red, target_edges)
    blue_target_counts = counts(target, blue, target_edges)
    red_reference_counts = counts(reference, red, reference_edges)
    blue_reference_counts = counts(reference, blue, reference_edges)
    red_n = int(red.sum())
    blue_n = int(blue.sum())
    if fit_red_n + fit_blue_n != selected_n:
        raise AssertionError("fit-label counts must conserve the selected upper-rank union")
    if red_n + blue_n != len(target):
        raise AssertionError("cluster prediction counts must conserve all finite pairs")
    if outside_fit_red_n + outside_fit_blue_n != len(target) - selected_n:
        raise AssertionError("outside-fit prediction counts must conserve excluded fit rows")

    return ClusterMarginals(
        selection_quantile=float(selection_quantile),
        finite_n=int(len(target)),
        target_selection_n=int(target_selected.sum()),
        reference_selection_n=int(reference_selected.sum()),
        selection_intersection_n=int((target_selected & reference_selected).sum()),
        selected_union_n=selected_n,
        fit_n=selected_n,
        fit_red_n=fit_red_n,
        fit_blue_n=fit_blue_n,
        prediction_n=int(len(target)),
        predicted_outside_fit_n=int(len(target) - selected_n),
        outside_fit_red_n=outside_fit_red_n,
        outside_fit_blue_n=outside_fit_blue_n,
        red_n=red_n,
        blue_n=blue_n,
        cluster_overlap_n=0,
        red_target_displayed_n=int(red_target_counts.sum()),
        blue_target_displayed_n=int(blue_target_counts.sum()),
        red_reference_displayed_n=int(red_reference_counts.sum()),
        blue_reference_displayed_n=int(blue_reference_counts.sum()),
        red_percentile_center=tuple(map(float, centers[red_cluster])),
        blue_percentile_center=tuple(map(float, centers[blue_cluster])),
        opposing=opposing,
        reference_high_color="blue" if opposing else "red",
        red_target_counts=red_target_counts,
        blue_target_counts=blue_target_counts,
        red_reference_counts=red_reference_counts,
        blue_reference_counts=blue_reference_counts,
    )


def _unique_tail_row(
    inputs: FigureInputs,
    donor: str,
    target_column: str,
    quantile: float,
) -> pd.Series:
    selected = inputs.rank_overlap.loc[
        inputs.rank_overlap["donor"].eq(donor)
        & inputs.rank_overlap["target_column"].eq(target_column)
        & np.isclose(inputs.rank_overlap["quantile"], quantile)
    ]
    if len(selected) != 1:
        raise ValueError(
            f"missing unique q{100 * quantile:g} rank summary for {donor} {target_column}"
        )
    return selected.iloc[0]


def _native_display_limits(
    inputs: FigureInputs,
    *,
    donor: str,
    pair,
    target: np.ndarray,
    reference: np.ndarray,
    finite: np.ndarray,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return donor-native q99.5 windows without transforming plotted values."""

    signed_kind = "qupath_acquisition_normalized_cell_mean"
    q995 = _unique_tail_row(inputs, donor, pair.target_column, 0.995)
    target_lower = (
        float(np.quantile(target[finite], 0.005))
        if pair.target_source_kind == signed_kind
        else 0.0
    )
    target_limits = (target_lower, float(q995["target_tail_boundary_score"]))
    reference_limits = (0.0, float(q995["reference_tail_boundary_score"]))
    if target_limits[0] >= target_limits[1] or reference_limits[0] >= reference_limits[1]:
        raise ValueError(f"{donor} {pair.target_column}: q99.5 display window is not increasing")
    return target_limits, reference_limits


def build_joint_histograms(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    *,
    bins: int,
    selection_quantile: float = CLUSTER_SELECTION_QUANTILE,
) -> dict[tuple[str, str], JointDistribution]:
    """Stream one donor/pair projection and retain only fixed-size histograms."""

    signed_kind = "qupath_acquisition_normalized_cell_mean"
    histograms: dict[tuple[str, str], JointDistribution] = {}
    for pair in roster.pairs:
        for donor in inputs.donors:
            frame, _paths, _summary = load_qupath_projection_donor(
                inputs.projection_root,
                donor,
                (pair.target_column, pair.reference_column),
            )
            object_ids = frame["object_id"].to_numpy()
            target = frame[pair.target_column].to_numpy(dtype=float)
            reference = frame[pair.reference_column].to_numpy(dtype=float)
            finite = np.isfinite(target) & np.isfinite(reference)
            if (
                pair.target_source_kind != signed_kind
                and (target[finite] < 0).any()
            ):
                raise ValueError(f"raw QuPath target {pair.target_column!r} contains negatives")
            if pair.reference_source_kind != "qupath_raw_cell_mean":
                raise ValueError("joint-distribution reference must be a raw QuPath cell mean")
            if (reference[finite] < 0).any():
                raise ValueError(
                    f"raw QuPath reference {pair.reference_column!r} contains negatives"
                )

            target_limits, reference_limits = _native_display_limits(
                inputs,
                donor=donor,
                pair=pair,
                target=target,
                reference=reference,
                finite=finite,
            )
            result = native_joint_histogram(
                target_values=target,
                reference_values=reference,
                target_limits=target_limits,
                reference_limits=reference_limits,
                bins=bins,
            )
            clusters = review_kmeans_cluster_marginals(
                object_ids=object_ids,
                target_values=target,
                reference_values=reference,
                target_key=pair.target_column,
                reference_key=pair.reference_column,
                target_limits=target_limits,
                reference_limits=reference_limits,
                bins=bins,
                selection_quantile=selection_quantile,
            )
            result = replace(result, cluster_marginals=clusters)
            audit = inputs.rank_overlap.loc[
                inputs.rank_overlap["donor"].eq(donor)
                & inputs.rank_overlap["target_column"].eq(pair.target_column)
            ]
            expected_n = set(audit["finite_pair_n"].astype(int))
            if expected_n != {result.finite_n}:
                raise ValueError(
                    f"{donor} {pair.target_column}:{pair.reference_column}: projection finite "
                    f"row count {result.finite_n} differs from rank-overlap summary "
                    f"{sorted(expected_n)}"
                )
            histograms[(donor, pair.target_column)] = result
            del frame, object_ids, target, reference, result, clusters
    return histograms


def _save_static_figure(fig, paths: Sequence[Path], *, dpi: int) -> None:
    for path in paths:
        save_dpi = dpi if path.suffix.lower() == ".png" else 300
        fig.savefig(path, dpi=save_dpi, bbox_inches="tight", facecolor="white")


def _shared_peak_marginal_lines(
    red_counts: np.ndarray,
    red_n: int,
    blue_counts: np.ndarray,
    blue_n: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize two cluster marginals to one shared within-margin peak."""

    red = np.asarray(red_counts, dtype=float) / red_n
    blue = np.asarray(blue_counts, dtype=float) / blue_n
    peak = float(max(red.max(), blue.max()))
    if not np.isfinite(peak) or peak <= 0:
        raise ValueError("cluster marginal has no positive displayed density")
    return red / peak, blue / peak


def _draw_rotated_marginal(
    ax,
    *,
    coordinates: np.ndarray,
    red: np.ndarray,
    blue: np.ndarray,
    orientation: str,
    coordinate_limits: tuple[float, float],
    red_color: str,
    blue_color: str,
) -> None:
    """Draw one marginal with formatting shared under a 90-degree rotation."""

    line_style = {
        "linewidth": MARGINAL_LINE_WIDTH,
        "solid_capstyle": "round",
        "clip_on": True,
    }
    if orientation == "vertical":
        ax.plot(red, coordinates, color=red_color, **line_style)
        ax.plot(blue, coordinates, color=blue_color, **line_style)
        ax.set_xlim(MARGINAL_DENSITY_HEADROOM, 0)
        ax.set_ylim(coordinate_limits)
        ax.set_xticks([])
        ax.tick_params(
            axis="y",
            left=True,
            right=False,
            labelleft=True,
            labelright=False,
            labelsize=7,
            length=2,
            pad=2,
        )
        hidden_spines = ("top", "bottom", "left")
        joint_spine = "right"
    elif orientation == "horizontal":
        ax.plot(coordinates, red, color=red_color, **line_style)
        ax.plot(coordinates, blue, color=blue_color, **line_style)
        ax.set_xlim(coordinate_limits)
        ax.set_ylim(MARGINAL_DENSITY_HEADROOM, 0)
        ax.set_yticks([])
        ax.tick_params(
            axis="x",
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
            labelsize=7,
            length=2,
            pad=2,
        )
        hidden_spines = ("left", "right", "bottom")
        joint_spine = "top"
    else:  # pragma: no cover - internal call contract
        raise ValueError(f"unknown marginal orientation: {orientation!r}")
    ax.spines[list(hidden_spines)].set_visible(False)
    ax.spines[joint_spine].set_visible(True)
    ax.spines[joint_spine].set_color("black")
    ax.spines[joint_spine].set_linewidth(0.8)


def _joint_figure_size(pair_n: int, donor_n: int) -> tuple[float, float]:
    """Size the figure so each outer grid cell, then each joint box, is square."""

    if pair_n < 1 or donor_n < 1:
        raise ValueError("joint figure requires at least one pair and donor")
    grid_width = JOINT_FIGURE_WIDTH_INCHES * (JOINT_GRID_RIGHT - JOINT_GRID_LEFT)
    cell_width = grid_width / (donor_n + JOINT_GRID_WSPACE * (donor_n - 1))
    grid_height = cell_width * (pair_n + JOINT_GRID_HSPACE * (pair_n - 1))
    figure_height = grid_height / (JOINT_GRID_TOP - JOINT_GRID_BOTTOM)
    return JOINT_FIGURE_WIDTH_INCHES, figure_height


def _joint_outer_grid(fig, *, pair_n: int, donor_n: int):
    return fig.add_gridspec(
        pair_n,
        donor_n,
        left=JOINT_GRID_LEFT,
        right=JOINT_GRID_RIGHT,
        bottom=JOINT_GRID_BOTTOM,
        top=JOINT_GRID_TOP,
        wspace=JOINT_GRID_WSPACE,
        hspace=JOINT_GRID_HSPACE,
    )


def _append_joint_axes(fig, slot):
    """Create one joint box with equal fixed left/bottom marginal geometry."""

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    joint_ax = fig.add_subplot(slot)
    divider = make_axes_locatable(joint_ax)
    target_ax = divider.append_axes(
        "left",
        size=MARGINAL_DEPTH_INCHES,
        pad=MARGINAL_PAD_INCHES,
        sharey=joint_ax,
    )
    reference_ax = divider.append_axes(
        "bottom",
        size=MARGINAL_DEPTH_INCHES,
        pad=MARGINAL_PAD_INCHES,
        sharex=joint_ax,
    )
    return target_ax, joint_ax, reference_ax


def render_joint_distributions(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    paths: Sequence[Path],
    *,
    dpi: int,
    bins: int = 160,
    selection_quantile: float = CLUSTER_SELECTION_QUANTILE,
) -> tuple[dict[str, object], ...]:
    """Render Figure 7-style native joint densities with two marginals."""

    plt, Line2D = _plot_modules()
    from matplotlib.colors import LogNorm

    histograms = build_joint_histograms(
        inputs,
        roster,
        bins=bins,
        selection_quantile=selection_quantile,
    )
    positive = np.concatenate([item.density[item.density > 0] for item in histograms.values()])
    if not len(positive):
        raise ValueError("joint-distribution histograms contain no occupied bins")
    normalization = LogNorm(vmin=float(positive.min()), vmax=float(positive.max()))
    cmap = plt.get_cmap("magma").with_extremes(bad="white")
    target_color = "#d73027"
    reference_color = "#2166ac"
    warning_color = "#9c2f00"

    pair_n = len(roster.pairs)
    donor_n = len(inputs.donors)
    fig = plt.figure(figsize=_joint_figure_size(pair_n, donor_n))
    outer = _joint_outer_grid(
        fig,
        pair_n=pair_n,
        donor_n=donor_n,
    )
    image = None
    panel_states: list[dict[str, object]] = []
    selection_label = f"q{100 * selection_quantile:g}"
    for row_index, pair in enumerate(roster.pairs):
        for column_index, donor in enumerate(inputs.donors):
            item = histograms[(donor, pair.target_column)]
            target_ax, joint_ax, reference_ax = _append_joint_axes(
                fig,
                outer[row_index, column_index],
            )

            target_edges = np.linspace(*item.target_limits, bins + 1)
            reference_edges = np.linspace(*item.reference_limits, bins + 1)
            target_centers = (target_edges[:-1] + target_edges[1:]) / 2
            reference_centers = (reference_edges[:-1] + reference_edges[1:]) / 2
            density = np.ma.masked_less_equal(item.density.T, 0)
            image = joint_ax.pcolormesh(
                reference_edges,
                target_edges,
                density,
                cmap=cmap,
                norm=normalization,
                shading="flat",
                rasterized=True,
            )
            joint_ax.set_xlim(item.reference_limits)
            joint_ax.set_ylim(item.target_limits)

            clusters = item.cluster_marginals
            if clusters is None:  # pragma: no cover - build contract
                raise ValueError("joint distribution lacks review-only cluster marginals")
            target_red, target_blue = _shared_peak_marginal_lines(
                clusters.red_target_counts,
                clusters.red_n,
                clusters.blue_target_counts,
                clusters.blue_n,
            )
            reference_red, reference_blue = _shared_peak_marginal_lines(
                clusters.red_reference_counts,
                clusters.red_n,
                clusters.blue_reference_counts,
                clusters.blue_n,
            )
            _draw_rotated_marginal(
                target_ax,
                coordinates=target_centers,
                red=target_red,
                blue=target_blue,
                orientation="vertical",
                coordinate_limits=item.target_limits,
                red_color=target_color,
                blue_color=reference_color,
            )
            _draw_rotated_marginal(
                reference_ax,
                coordinates=reference_centers,
                red=reference_red,
                blue=reference_blue,
                orientation="horizontal",
                coordinate_limits=item.reference_limits,
                red_color=target_color,
                blue_color=reference_color,
            )
            joint_ax.tick_params(labelleft=False, labelbottom=False, length=2)

            row = _unique_tail_row(inputs, donor, pair.target_column, 0.99)
            survey_row = inputs.pair_diagnostics.loc[
                inputs.pair_diagnostics["donor"].eq(donor)
                & inputs.pair_diagnostics["target_column"].eq(pair.target_column)
            ]
            if len(survey_row) != 1:  # pragma: no cover - validated complete grid
                raise ValueError(f"missing unique pair survey for {donor} {pair.target_column}")
            pearson = float(survey_row.iloc[0]["pearson"])
            displayed_percent = 100 * item.displayed_n / item.finite_n
            joint_ax.text(
                0.025,
                0.975,
                f"n={item.finite_n:,}; shown={displayed_percent:.2f}%\n"
                f"r={pearson:.3f}; ρ={row['spearman_score']:.3f}\n"
                f"{selection_label} fit R/B={clusters.fit_red_n:,}/{clusters.fit_blue_n:,}; "
                f"all R/B={clusters.red_n:,}/{clusters.blue_n:,}\n"
                f"q99 both={int(row['both_tail_n']):,} "
                f"({100 * row['both_of_target_tail']:.2f}%)",
                transform=joint_ax.transAxes,
                ha="left",
                va="top",
                fontsize=6.8,
                color="black",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 2},
            )
            if not clusters.opposing:
                joint_ax.text(
                    0.975,
                    0.035,
                    "NON-OPPOSING\nreference-high center is red",
                    transform=joint_ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=6.8,
                    fontweight="bold",
                    color=warning_color,
                    bbox={
                        "facecolor": "#fff3e6",
                        "edgecolor": warning_color,
                        "linewidth": 0.65,
                        "alpha": 0.94,
                        "pad": 2,
                    },
                )
            panel_states.append(
                {
                    "donor": donor,
                    "target_column": pair.target_column,
                    "reference_column": pair.reference_column,
                    "state": "OPPOSING" if clusters.opposing else "NON_OPPOSING",
                    "reference_high_color": clusters.reference_high_color,
                    "finite_n": clusters.finite_n,
                    "selection_quantile": clusters.selection_quantile,
                    "target_selection_n": clusters.target_selection_n,
                    "target_selection_fraction": (
                        clusters.target_selection_n / clusters.finite_n
                    ),
                    "reference_selection_n": clusters.reference_selection_n,
                    "reference_selection_fraction": (
                        clusters.reference_selection_n / clusters.finite_n
                    ),
                    "selection_intersection_n": clusters.selection_intersection_n,
                    "selection_intersection_fraction": (
                        clusters.selection_intersection_n / clusters.finite_n
                    ),
                    "selected_union_n": clusters.selected_union_n,
                    "selected_union_fraction": clusters.selected_union_n / clusters.finite_n,
                    "fit_n": clusters.fit_n,
                    "fit_red_n": clusters.fit_red_n,
                    "fit_blue_n": clusters.fit_blue_n,
                    "prediction_n": clusters.prediction_n,
                    "predicted_outside_fit_n": clusters.predicted_outside_fit_n,
                    "predicted_outside_fit_fraction": (
                        clusters.predicted_outside_fit_n / clusters.finite_n
                    ),
                    "outside_fit_red_n": clusters.outside_fit_red_n,
                    "outside_fit_blue_n": clusters.outside_fit_blue_n,
                    "red_n": clusters.red_n,
                    "blue_n": clusters.blue_n,
                    "cluster_overlap_n": clusters.cluster_overlap_n,
                    "target_display_limits": list(item.target_limits),
                    "reference_display_limits": list(item.reference_limits),
                    "displayed_n": item.displayed_n,
                    "red_percentile_center": {
                        "target": clusters.red_percentile_center[0],
                        "reference": clusters.red_percentile_center[1],
                    },
                    "blue_percentile_center": {
                        "target": clusters.blue_percentile_center[0],
                        "reference": clusters.blue_percentile_center[1],
                    },
                    "red_target_displayed_n": clusters.red_target_displayed_n,
                    "blue_target_displayed_n": clusters.blue_target_displayed_n,
                    "red_reference_displayed_n": clusters.red_reference_displayed_n,
                    "blue_reference_displayed_n": clusters.blue_reference_displayed_n,
                }
            )
            if row_index == 0:
                joint_ax.set_title(f"donor {donor}", fontsize=11, pad=5)
            if column_index == 0:
                target_ax.set_ylabel(
                    f"{pair.target_column}\nnative QuPath value",
                    fontsize=8.4,
                    labelpad=4,
                )
            if row_index == len(roster.pairs) - 1:
                reference_ax.set_xlabel(
                    f"{pair.reference_column} native QuPath value",
                    fontsize=8.4,
                    labelpad=3,
                )
    if image is None:  # pragma: no cover - roster contract is nonempty
        raise ValueError("diagnostic roster contains no pairs")
    colorbar_axis = fig.add_axes((0.915, 0.17, 0.018, 0.62))
    colorbar = fig.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("fraction of all finite pairs per bin (log scale)")
    handles = (
        Line2D(
            [0],
            [0],
            color=target_color,
            linewidth=2,
            label="red: y-axis marker",
        ),
        Line2D(
            [0],
            [0],
            color=reference_color,
            linewidth=2,
            label="blue: x-axis marker",
        ),
    )
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=2)
    fig.suptitle(
        f"Upstream QuPath joint distributions for the {pair_n} rostered comparisons",
        y=0.99,
        fontsize=15,
    )
    _save_static_figure(fig, paths, dpi=dpi)
    plt.close(fig)
    return tuple(panel_states)


def render_pair_correlations(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    path: Path,
    *,
    dpi: int,
) -> None:
    """Plot native-value and rank correlations, separated by comparison."""

    plt, Line2D = _plot_modules()
    pairs = roster.pairs
    fig, axes = plt.subplots(
        1,
        len(pairs),
        figsize=(4.6 * len(pairs), 4.8),
        sharey=True,
        squeeze=False,
    )
    x = np.arange(len(inputs.donors), dtype=float)
    for index, pair in enumerate(pairs):
        ax = axes[0, index]
        rows = (
            inputs.pair_diagnostics.loc[
                inputs.pair_diagnostics["target_column"].eq(pair.target_column)
            ]
            .set_index("donor")
            .loc[list(inputs.donors)]
        )
        ax.scatter(x - 0.09, rows["pearson"], color="#31688e", marker="o", s=54, zorder=3)
        ax.scatter(x + 0.09, rows["spearman"], color="#d1495b", marker="s", s=48, zorder=3)
        for position, (_, row) in enumerate(rows.iterrows()):
            ax.plot(
                [position - 0.09, position + 0.09],
                [row["pearson"], row["spearman"]],
                color="#b8b8b8",
                linewidth=0.9,
                zorder=1,
            )
        ax.axhline(0, color="#555555", linewidth=0.8)
        ax.set_xticks(x, inputs.donors)
        ax.set_ylim(-1.04, 1.04)
        ax.grid(axis="y", alpha=0.22)
        ax.set_title(_pair_title(pair), fontsize=11)
        ax.set_xlabel("donor")
    axes[0, 0].set_ylabel("within-donor correlation")
    handles = [
        Line2D(
            [], [], color="#31688e", marker="o", linestyle="none", label="Pearson (native values)"
        ),
        Line2D([], [], color="#d1495b", marker="s", linestyle="none", label="Spearman (ranks)"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.96), ncol=2)
    fig.suptitle("Fundamental pair association by donor", fontsize=14, y=1.02)
    fig.text(
        0.5,
        0.01,
        "Correlations are dimensionless and each panel is a separate comparison. "
        "Iba1_norm and raw means do not share a numeric scale; these are descriptive inputs, "
        "not thresholds or marker calls.",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.88))
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _feature_order(roster: DiagnosticRoster) -> tuple[str, ...]:
    targets = [pair.target_column for pair in roster.pairs]
    references = [pair.reference_column for pair in roster.pairs]
    return tuple(dict.fromkeys([*targets, *references]))


def render_native_scales(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    path: Path,
    *,
    dpi: int,
) -> None:
    """Plot each input's quantiles on its own explicitly labelled native axis."""

    plt, Line2D = _plot_modules()
    columns = _feature_order(roster)
    ncols = min(2, len(columns))
    nrows = int(np.ceil(len(columns) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.1 * ncols, 3.35 * nrows), squeeze=False)
    colors, markers = _donor_style(inputs.donors, plt)

    for ax, column in zip(axes.ravel(), columns, strict=False):
        rows = (
            inputs.input_summary.loc[inputs.input_summary["column"].eq(column)]
            .set_index("donor")
            .loc[list(inputs.donors)]
        )
        kind = str(rows["source_kind"].iloc[0])
        measurement = str(rows["source_measurement"].iloc[0])
        y = np.arange(len(inputs.donors), dtype=float)
        for position, (donor, row) in enumerate(rows.iterrows()):
            color = colors[donor]
            ax.plot(
                [row["minimum"], row["maximum"]],
                [position, position],
                color=color,
                alpha=0.25,
                linewidth=2.0,
            )
            ax.plot(
                [row["q900"], row["q990"]],
                [position, position],
                color=color,
                linewidth=6.0,
                solid_capstyle="round",
            )
            ax.scatter(
                row["median"],
                position,
                color=color,
                marker=markers[donor],
                s=45,
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )
            if kind == "qupath_acquisition_normalized_cell_mean":
                negative = 100 * row["negative_n"] / row["finite_n"] if row["finite_n"] else 0
                ax.annotate(
                    f"{negative:.1f}% < 0",
                    (row["q990"], position),
                    xytext=(7, 0),
                    textcoords="offset points",
                    va="center",
                    fontsize=8,
                    color=color,
                )
        ax.set_yticks(y, inputs.donors)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.2)
        ax.set_title(f"{column}\n{measurement}", fontsize=10.5)
        if kind == "qupath_acquisition_normalized_cell_mean":
            ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
            ax.axvline(1, color="#8c6d31", linewidth=1.0, linestyle=":")
            ax.set_xlabel("signed QuPath Mean (norm), linear axis")
        else:
            ax.set_xscale("symlog", linthresh=1.0, linscale=0.5, base=10)
            ax.set_xlabel("raw QuPath cell mean, symlog axis")
    for ax in axes.ravel()[len(columns) :]:
        ax.axis("off")

    handles = [
        Line2D([], [], color="#777777", linewidth=2, alpha=0.35, label="minimum–maximum"),
        Line2D([], [], color="#777777", linewidth=6, label="q90–q99"),
        Line2D([], [], color="#777777", marker="o", linestyle="none", label="median"),
        Line2D([], [], color="#8c6d31", linestyle=":", label="Iba1 normalization anchor (1)"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=4)
    fig.suptitle("Upstream input summaries on independent native scales", fontsize=14, y=1.01)
    fig.text(
        0.5,
        0.012,
        "Axes are deliberately not shared across measurements. For Iba1_norm, 1 is the QuPath "
        "tissue-p99 normalization anchor—not a biological positivity threshold.",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.065, 1, 0.89))
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_rank_overlap(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    paths: Sequence[Path],
    *,
    dpi: int,
) -> None:
    """Plot exact-rank upper-tail intersection and independence lift."""

    plt, Line2D = _plot_modules()
    pairs = roster.pairs
    fig, axes = plt.subplots(
        len(pairs),
        2,
        figsize=(12.8, 3.2 * len(pairs)),
        sharex=True,
        squeeze=False,
    )
    colors, markers = _donor_style(inputs.donors, plt)
    tail_percent = 100 * (1 - np.asarray(roster.required_tail_quantiles, dtype=float))

    for row_index, pair in enumerate(pairs):
        overlap_ax, lift_ax = axes[row_index]
        pair_rows = inputs.rank_overlap.loc[
            inputs.rank_overlap["target_column"].eq(pair.target_column)
        ]
        for donor_index, donor in enumerate(inputs.donors):
            rows = pair_rows.loc[pair_rows["donor"].eq(donor)].sort_values(
                "quantile", ascending=True
            )
            x = 100 * (1 - rows["quantile"].to_numpy(dtype=float))
            lift = rows["tail_independence_lift"].to_numpy(dtype=float)
            overlap_ax.plot(
                x,
                100 * rows["both_of_target_tail"],
                color=colors[donor],
                marker=markers[donor],
                linewidth=1.7,
                markersize=5,
            )
            lift_ax.plot(
                x,
                lift,
                color=colors[donor],
                marker=markers[donor],
                linewidth=1.7,
                markersize=5,
            )
            for retained, value, both_n in zip(
                x,
                lift,
                rows["both_tail_n"].to_numpy(dtype=int),
                strict=True,
            ):
                if both_n <= 10:
                    lift_ax.annotate(
                        f"n={both_n}",
                        (retained, value),
                        xytext=(0, 5 + 6 * donor_index),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=6.4,
                        color=colors[donor],
                    )
        overlap_ax.plot(
            tail_percent,
            tail_percent,
            color="#666666",
            linestyle="--",
            linewidth=1.1,
            label="independent expectation",
        )
        overlap_ax.set_ylabel(
            f"{pair.target_column} vs {pair.reference_column}\nintersection (% of either tail)"
        )
        overlap_ax.set_ylim(bottom=-0.5)
        lift_ax.axhline(1, color="#666666", linestyle="--", linewidth=1.1)
        lift_ax.set_yscale("symlog", linthresh=1.0, linscale=0.5, base=10)
        lift_ax.set_ylabel("observed / independent\nexpectation (symlog)")
        for ax in (overlap_ax, lift_ax):
            ax.set_xscale("log")
            ax.set_xticks(tail_percent, [f"{value:g}" for value in tail_percent])
            ax.set_xlim(max(tail_percent) * 1.25, min(tail_percent) / 1.25)
            ax.grid(alpha=0.2)
    axes[0, 0].set_title("Exact-rank upper-tail overlap", fontsize=12)
    axes[0, 1].set_title("Enrichment relative to independent tails", fontsize=12)
    for ax in axes[-1]:
        ax.set_xlabel("upper tail retained per measurement (% of cells; stricter →)")

    donor_handles = [
        Line2D(
            [],
            [],
            color=colors[donor],
            marker=markers[donor],
            label=f"donor {donor}",
        )
        for donor in inputs.donors
    ]
    donor_handles.append(
        Line2D([], [], color="#666666", linestyle="--", label="independent expectation")
    )
    fig.legend(
        handles=donor_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=4,
    )
    fig.suptitle(
        f"Rank-only overlap of the {len(pairs)} rostered comparisons", fontsize=14, y=1.005
    )
    fig.text(
        0.5,
        0.006,
        "Each measurement is ranked independently within donor; no shared numeric cutoff is "
        "applied. Mixed raw and signed-normalized scales are therefore explicit and intentional. "
        "No threshold or cell call is created.",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))
    _save_static_figure(fig, paths, dpi=dpi)
    plt.close(fig)


def _montage_rows(inputs: FigureInputs, pair) -> pd.DataFrame:
    """Select extreme controls plus all four dual-bright examples for one pair."""

    pair_id = f"{pair.target_column}:{pair.reference_column}"
    rows = inputs.crop_manifest.loc[inputs.crop_manifest["pair_id"].eq(pair_id)].copy()
    selected = pd.concat(
        [
            rows.loc[rows["stratum"].eq("target_bright") & rows["crop_rank"].eq(1)],
            rows.loc[rows["stratum"].eq("dual_bright")].sort_values("crop_rank"),
            rows.loc[rows["stratum"].eq("reference_bright") & rows["crop_rank"].eq(1)],
        ],
        ignore_index=True,
    )
    if len(selected) != 6:
        raise ValueError(
            f"{pair_id}: montage requires target/reference-bright extremes and four dual crops"
        )
    return selected


def _crop_short_label(row: pd.Series) -> str:
    names = {
        "target_bright": "Target-bright extreme",
        "dual_bright": f"Dual {int(row['crop_rank'])}",
        "reference_bright": "Reference-bright extreme",
    }
    return names[str(row["stratum"])]


def _crop_annotation(row: pd.Series) -> str:
    return (
        f"pT={row['target_percentile_rank']:.3f}, "
        f"pR={row['reference_percentile_rank']:.3f}\n"
        f"neighbors={int(row['contact_neighbor_n'])}, "
        f"boundary={row['contact_boundary_px']:.0f}px"
    )


def _read_crop_image(path: object) -> np.ndarray:
    plt, _Line2D = _plot_modules()
    image = np.asarray(plt.imread(str(path)))
    if image.ndim not in (2, 3) or image.shape[1] < 8:
        raise ValueError(f"rendered boundary crop has an invalid image shape: {path}")
    return image


def render_crop_montage_png(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    path: Path,
    *,
    dpi: int,
) -> None:
    """Render one roster-pair-by-six RGB overview while retaining all source composites."""

    plt, _Line2D = _plot_modules()
    fig, axes = plt.subplots(
        len(roster.pairs),
        6,
        figsize=(18, 2.4 * len(roster.pairs) + 1.2),
        squeeze=False,
    )
    for row_index, pair in enumerate(roster.pairs):
        rows = _montage_rows(inputs, pair)
        for column_index, (_, row) in enumerate(rows.iterrows()):
            ax = axes[row_index, column_index]
            source = _read_crop_image(row["render_path"])
            rgb_quarter = source[:, : source.shape[1] // 4]
            ax.imshow(rgb_quarter)
            ax.set_axis_off()
            ax.set_title(
                f"{_crop_short_label(row)}\n{_crop_annotation(row)}",
                fontsize=7.2,
                pad=3,
            )
            if column_index == 0:
                display_target = str(row["target_display_marker"])
                label = f"{pair.target_column} vs {pair.reference_column}"
                if pair.target_column == "Iba1_norm":
                    label += "\nnorm selection; raw Iba1 shown"
                elif display_target != pair.target_column:
                    label += f"\nraw {display_target} shown"
                ax.text(
                    -0.06,
                    0.5,
                    label,
                    transform=ax.transAxes,
                    rotation=90,
                    ha="right",
                    va="center",
                    fontsize=8.3,
                    fontweight="bold",
                )
    fig.suptitle(
        f"Donor {inputs.crop_donor}: raw-channel boundary/contact crop overview",
        fontsize=14,
        y=0.99,
    )
    fig.text(
        0.5,
        0.008,
        "RGB = raw target / raw Ker8_18 / raw DAPI; gold outlines the selected cell. "
        "These deterministic rank-extreme crops are deliberately nonrepresentative qualitative "
        "checks, not labels or training data. Iba1 objects were selected by exact signed QuPath "
        "Mean (norm), while display uses the acquired raw Iba1 channel.",
        ha="center",
        va="bottom",
        fontsize=8.5,
    )
    fig.subplots_adjust(left=0.045, right=0.995, top=0.89, bottom=0.09, wspace=0.045, hspace=0.22)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_crop_montage_pdf(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    path: Path,
) -> None:
    """Write one pair per page, preserving each full RGB/target/reference/DAPI composite."""

    plt, _Line2D = _plot_modules()
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        for pair in roster.pairs:
            rows = _montage_rows(inputs, pair)
            rows = rows.iloc[[0, 5, 1, 2, 3, 4]].reset_index(drop=True)
            fig, axes = plt.subplots(3, 2, figsize=(16, 10.2), squeeze=False)
            for ax, (_, row) in zip(axes.ravel(), rows.iterrows(), strict=True):
                source = _read_crop_image(row["render_path"])
                ax.imshow(source)
                ax.set_axis_off()
                ax.set_title(
                    f"{_crop_short_label(row)} · {_crop_annotation(row).replace(chr(10), ' · ')}",
                    fontsize=8,
                    pad=3,
                )
            display_note = ""
            if pair.target_column == "Iba1_norm":
                display_note = (
                    " Iba1 selection uses exact signed QuPath Mean (norm); "
                    "the panels show raw Iba1."
                )
            fig.suptitle(
                f"Donor {inputs.crop_donor} · {pair.target_column}:{pair.reference_column} · "
                f"raw RGB / target / reference / DAPI composites",
                fontsize=13,
                y=0.985,
            )
            fig.text(
                0.5,
                0.012,
                "Rank-extreme, deterministic crops are nonrepresentative qualitative checks; "
                "gold outlines the selected cell, cyan nuclei, white other cells." + display_note,
                ha="center",
                va="bottom",
                fontsize=8.2,
            )
            fig.subplots_adjust(
                left=0.02, right=0.98, top=0.94, bottom=0.06, hspace=0.18, wspace=0.04
            )
            pdf.savefig(fig, dpi=300, bbox_inches="tight", facecolor="white")
            plt.close(fig)


def render_figures(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    output: Path,
    *,
    dpi: int = 160,
    selection_quantile: float = CLUSTER_SELECTION_QUANTILE,
) -> tuple[tuple[Path, ...], tuple[dict[str, object], ...]]:
    """Render all review figures into an existing owned workspace."""

    if dpi < 72:
        raise ValueError("dpi must be at least 72")
    paths = {name: output / name for name in FIGURE_NAMES}
    cluster_panels = render_joint_distributions(
        inputs,
        roster,
        (
            paths["fundamental_pair_joint_distributions.png"],
            paths["fundamental_pair_joint_distributions.pdf"],
        ),
        dpi=dpi,
        selection_quantile=selection_quantile,
    )
    render_rank_overlap(
        inputs,
        roster,
        (paths["rank_tail_overlap.png"], paths["rank_tail_overlap.pdf"]),
        dpi=dpi,
    )
    render_crop_montage_png(
        inputs,
        roster,
        paths["crop_montage_6450.png"],
        dpi=dpi,
    )
    render_crop_montage_pdf(inputs, roster, paths["crop_montage_6450.pdf"])
    return tuple(paths[name] for name in FIGURE_NAMES), cluster_panels


def write_manifest(
    inputs: FigureInputs,
    roster: DiagnosticRoster,
    output: Path,
    figure_paths: Sequence[Path],
    cluster_panels: Sequence[Mapping[str, object]],
    *,
    selection_quantile: float = CLUSTER_SELECTION_QUANTILE,
) -> Path:
    uses_default_selection = bool(
        np.isclose(selection_quantile, CLUSTER_SELECTION_QUANTILE, rtol=0, atol=1e-12)
    )
    selection_label = f"q{100 * selection_quantile:g}"
    manifest = {
        "status": "complete",
        "method_version": METHOD_VERSION,
        "method_status": "development_mutable",
        "roster_version": roster.roster_version,
        "donors": list(inputs.donors),
        "crop_donor": inputs.crop_donor,
        "pairs": [
            {
                "target": pair.target,
                "target_column": pair.target_column,
                "reference": pair.reference,
                "reference_column": pair.reference_column,
                "target_source_measurement": pair.target_source_measurement,
                "target_source_kind": pair.target_source_kind,
                "reference_source_measurement": pair.reference_source_measurement,
                "reference_source_kind": pair.reference_source_kind,
            }
            for pair in roster.pairs
        ],
        "source_files": [str(path) for path in inputs.source_files],
        "source_method_versions": {
            key: list(values) for key, values in inputs.source_method_versions.items()
        },
        "figures": [path.name for path in figure_paths],
        "source_semantics": (
            "current narrow upstream QuPath projection, compact survey and donor-local "
            "exact-rank overlap summaries, plus boundary-audit raw-channel crops"
        ),
        "mixed_scale_policy": (
            "joint distributions use separate donor-specific native-linear axes, while pair "
            "overlap ranks each measurement independently within donor; no common native-value "
            "boundary is used"
        ),
        "fundamental_pair_joint_distributions": {
            "bins_per_axis": 160,
            "axes": "donor-specific native-linear QuPath values",
            "raw_display_window": "0 through donor-local exact-rank q99.5 boundary",
            "signed_iba1_display_window": (
                "finite-pair empirical q0.5 through donor-local exact-rank q99.5 boundary"
            ),
            "outside_window_handling": "filtered, never clipped into edge bins",
            "statistics_population": "all finite pairs, including cells outside the display window",
            "density": (
                "fraction of all finite donor/pair cells per native-value bin; one common "
                "logarithmic color normalization"
            ),
            "marginals": (
                "left target and bottom reference each overlay both mutually exclusive review-only "
                "KMeans clusters predicted over all finite cells as lines; no unconditional fill"
            ),
            "marginal_geometry": {
                "depth_inches": MARGINAL_DEPTH_INCHES,
                "pad_inches": MARGINAL_PAD_INCHES,
                "density_headroom": MARGINAL_DENSITY_HEADROOM,
                "line_width_points": MARGINAL_LINE_WIDTH,
                "rotation_contract": (
                    "bottom reference marginal is the exact 90-degree formatting counterpart "
                    "of the left target marginal"
                ),
                "joint_box": "physically square after equal fixed left/bottom margins",
            },
            "legend": ["red: y-axis marker", "blue: x-axis marker"],
            "bottom_explanatory_text": "none; scientific caveats remain in this manifest",
            "cluster_populations": {
                "method_version": CLUSTER_METHOD_VERSION,
                "method_status": "review_only_unverified",
                "selection": (
                    f"exact donor-local {selection_label} upper-tail union of either marker"
                ),
                "selection_quantile": selection_quantile,
                "default_selection_quantile": CLUSTER_SELECTION_QUANTILE,
                "selection_configuration": (
                    "manuscript_anchored_default"
                    if uses_default_selection
                    else "explicit_review_sensitivity_variant"
                ),
                "selection_rationale": (
                    "q70 was pinned as the closest reciprocal longitudinal analogue of the "
                    "RESTORE Figure 7 grouping; it was not optimized against these panels"
                ),
                "selection_dependence": (
                    "marginal-line separation and overlap depend on the clusters fitted from the "
                    "selected upper-rank union and can change with the selection quantile"
                ),
                "selection_tie_break": "object_id_ascending",
                "preprocessing": (
                    "continuous donor-local 0–1 bright-percentile ranks over all finite pairs; "
                    "descending native value with object_id ascending ties; no native-value "
                    "scale enters KMeans"
                ),
                "feature_range": [0.0, 1.0],
                "brightest_feature_value": 1.0,
                "algorithm": "sklearn KMeans",
                "n_clusters": 2,
                "init": "k-means++",
                "n_init": CLUSTER_N_INIT,
                "max_iter": 300,
                "algorithm_variant": "lloyd",
                "random_state": CLUSTER_RANDOM_STATE,
                "fit_population": f"every cell in the {selection_label} union",
                "fit_row_order": "object_id ascending",
                "feature_order": "canonical marker-column order during fit and prediction",
                "prediction_population": "every finite donor/pair cell",
                "outside_fit_policy": (
                    f"cells outside the {selection_label} fit union are predicted into one of the "
                    "two fitted clusters and included in both native-axis marginal curves"
                ),
                "third_population": (
                    "none; the manuscript's two-group construction assigns its double-low Y3 "
                    "population to one of the two groups rather than plotting a third curve"
                ),
                "overlap_policy": "impossible; each finite cell receives exactly one cluster",
                "red_mapping": "cluster with higher target bright-percentile center",
                "blue_mapping": (
                    "cluster with higher reference bright-percentile center when distinct; "
                    "otherwise other cluster and panel flagged NON-OPPOSING"
                ),
                "projection": (
                    "both all-finite predicted clusters projected without tail truncation onto "
                    "both native-value margins"
                ),
                "density": (
                    "displayed bin count divided by full cluster count; red and blue share one "
                    "peak normalization within each margin"
                ),
                "interpretation": (
                    "exploratory assignment only; not a verified population, inferred baseline, "
                    "threshold, or cell call"
                ),
                "panels": [dict(panel) for panel in cluster_panels],
            },
            "q99_guides": "none",
            "sampling": (
                "none; joint density uses all finite cells, KMeans fits the full selected union, "
                "and prediction covers all finite cells"
            ),
            "smoothing": "none",
            "retained_cell_level_arrays": False,
        },
        "crop_montage": {
            "selection": (
                "deterministic rank extremes; one target-bright and one reference-bright "
                "extreme plus four dual-bright examples per pair; deliberately nonrepresentative"
            ),
            "display": "raw target/raw Ker8_18/raw DAPI qptiff channels",
            "iba1_selection": "exact signed Cell: Iba1: Mean (norm)",
            "iba1_display": "raw acquired Iba1 channel",
            "pdf_pages": len(roster.pairs),
        },
        "restore_outputs_read": False,
        "creates_cell_calls": False,
        "creates_thresholds": False,
        "mutable_workspace": True,
        "production_outputs_written": False,
    }
    path = output / "figure_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projection-root", type=Path, default=DEFAULT_PROJECTION_ROOT)
    parser.add_argument("--survey-root", type=Path, default=DEFAULT_SURVEY_ROOT)
    parser.add_argument("--overlap-root", type=Path, default=DEFAULT_OVERLAP_ROOT)
    parser.add_argument("--crop-root", type=Path, default=DEFAULT_CROP_ROOT)
    parser.add_argument("--development-root", type=Path, default=DEFAULT_DEVELOPMENT_ROOT)
    parser.add_argument("--label", default=DEFAULT_LABEL)
    parser.add_argument("--roster", type=Path, default=ROSTER_PATH)
    parser.add_argument("--donors", nargs="+")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument(
        "--cluster-selection-quantile",
        type=float,
        default=CLUSTER_SELECTION_QUANTILE,
        help=(
            "review-only upper-rank union quantile (default: 0.70); non-default "
            "sensitivity variants require an explicit non-default --label"
        ),
    )
    args = parser.parse_args(argv)

    try:
        if not np.isfinite(args.cluster_selection_quantile) or not (
            0 < args.cluster_selection_quantile < 1
        ):
            raise ValueError("cluster selection quantile must be finite and between zero and one")
        uses_default_selection = bool(
            np.isclose(
                args.cluster_selection_quantile,
                CLUSTER_SELECTION_QUANTILE,
                rtol=0,
                atol=1e-12,
            )
        )
        if not uses_default_selection and args.label == DEFAULT_LABEL:
            raise ValueError(
                "non-default cluster-selection sensitivity variants require a non-default "
                "--label and cannot replace the pinned q70 workspace"
            )
        roster = load_diagnostic_roster(args.roster)
        inputs = load_figure_inputs(
            args.projection_root,
            args.survey_root,
            args.overlap_root,
            args.crop_root,
            roster=roster,
            donors=args.donors,
        )
        _plot_modules()  # Fail before replacing an existing owned workspace.
        candidate_output = (
            Path(args.development_root).expanduser().resolve() / args.label
        ).resolve()
        input_roots = {
            Path(args.projection_root).expanduser().resolve(),
            Path(args.survey_root).expanduser().resolve(),
            Path(args.overlap_root).expanduser().resolve(),
            Path(args.crop_root).expanduser().resolve(),
        }
        if candidate_output in input_roots:
            raise ValueError("figure output must not replace an input workspace")
        output = prepare_development_workspace(args.development_root, args.label)
        figures, cluster_panels = render_figures(
            inputs,
            roster,
            output,
            dpi=args.dpi,
            selection_quantile=args.cluster_selection_quantile,
        )
        write_manifest(
            inputs,
            roster,
            output,
            figures,
            cluster_panels,
            selection_quantile=args.cluster_selection_quantile,
        )
        shutil.copyfile(args.roster, output / "roster_used.json")
    except (RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    print(f"rendered {len(figures)} review figures in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
