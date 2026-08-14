#!/usr/bin/env python3
"""Summarize the three fundamental upstream QuPath input measurements.

The active development roster has two unique comparisons: raw QuPath CD3e and
QuPath acquisition-normalized Iba1, each against the raw QuPath Ker8_18 cell
mean. CD79a is not a pair candidate because its expression is too sparse in
normal tissue. These are inputs to RESTORE, not RESTORE outputs.

The three measurements do not share a call scale.  In particular, QuPath's
signed Iba1 normalization places the tissue median near zero and tissue p99
near one; one is not a biological positivity boundary.  This tool therefore
reports measurement coverage and association only.  It never assigns a cell
type or creates a positive call.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
WORKSPACE = REPO.parents[1]
sys.path.insert(0, str(REPO))

from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    DEFAULT_OUTPUT_ROOT as DEFAULT_INPUT_ROOT,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    DEFAULT_SPEC,
    DIAGNOSTIC_SOURCE_KINDS,
    SOURCE_MEASUREMENTS,
    load_qupath_projection_donor,
    refresh_qupath_projection,
)

METHOD_VERSION = "fundamental-qupath-input-survey-v3"
ROSTER_PATH = REPO / "phenocycler" / "fundamental_pair_diagnostics.json"
DEFAULT_DEVELOPMENT_ROOT = WORKSPACE / "data" / "development"
DEFAULT_LABEL = "fundamental_pair_diagnostics"
WORKSPACE_MARKER = ".fundamental_pair_survey_workspace"
ROSTER_SCHEMA_VERSION = 3


@dataclass(frozen=True)
class DiagnosticPair:
    target: str
    target_column: str
    target_source_measurement: str
    target_source_kind: str
    target_class: str
    reference: str
    reference_column: str
    reference_source_measurement: str
    reference_source_kind: str
    reference_class: str

    @property
    def specification(self) -> tuple[str, ...]:
        return (
            self.target,
            self.target_column,
            self.target_source_measurement,
            self.target_source_kind,
            self.target_class,
            self.reference,
            self.reference_column,
            self.reference_source_measurement,
            self.reference_source_kind,
            self.reference_class,
        )


@dataclass(frozen=True)
class DiagnosticRoster:
    schema_version: int
    roster_version: str
    method_status: str
    description: str
    required_tail_quantiles: tuple[float, ...]
    pairs: tuple[DiagnosticPair, ...]

    @property
    def source_columns(self) -> tuple[str, ...]:
        return tuple(
            dict.fromkeys(
                column
                for pair in self.pairs
                for column in (pair.target_column, pair.reference_column)
            )
        )

    @property
    def source_contract(self) -> dict[str, tuple[str, str]]:
        contract: dict[str, tuple[str, str]] = {}
        for pair in self.pairs:
            for column, measurement, kind in (
                (
                    pair.target_column,
                    pair.target_source_measurement,
                    pair.target_source_kind,
                ),
                (
                    pair.reference_column,
                    pair.reference_source_measurement,
                    pair.reference_source_kind,
                ),
            ):
                prior = contract.setdefault(column, (measurement, kind))
                if prior != (measurement, kind):
                    raise ValueError(f"source column {column} has inconsistent provenance")
        return contract


PAIR_FIELDS = {
    "target",
    "target_column",
    "target_source_measurement",
    "target_source_kind",
    "target_class",
    "reference",
    "reference_column",
    "reference_source_measurement",
    "reference_source_kind",
    "reference_class",
}


def roster_from_dict(payload: Mapping[str, object]) -> DiagnosticRoster:
    """Parse the editable roster and pin every alias to an exact QuPath header."""

    if int(payload.get("schema_version", -1)) != ROSTER_SCHEMA_VERSION:
        raise ValueError(f"fundamental roster schema_version must be {ROSTER_SCHEMA_VERSION}")
    if payload.get("method_status") != "development_mutable":
        raise ValueError("fundamental roster must be development_mutable")
    quantiles = tuple(float(value) for value in payload.get("required_tail_quantiles", ()))
    if (
        not quantiles
        or tuple(sorted(set(quantiles))) != quantiles
        or any(not 0 < value < 1 for value in quantiles)
    ):
        raise ValueError("required_tail_quantiles must be unique, increasing, and in (0,1)")
    raw_pairs = payload.get("pairs")
    if not isinstance(raw_pairs, list) or not raw_pairs:
        raise ValueError("fundamental roster pairs must be a nonempty list")

    pairs: list[DiagnosticPair] = []
    for index, row in enumerate(raw_pairs):
        if not isinstance(row, Mapping):
            raise ValueError(f"fundamental pair row {index} must be an object")
        missing = sorted(PAIR_FIELDS.difference(row))
        if missing:
            raise ValueError(f"fundamental pair row {index} is missing {missing}")
        pair = DiagnosticPair(**{field: str(row[field]) for field in PAIR_FIELDS})
        if pair.target == pair.reference or pair.target_column == pair.reference_column:
            raise ValueError(f"circular diagnostic pair {pair.target} <- itself")
        if pair.target_class == pair.reference_class:
            raise ValueError(f"{pair.target} versus {pair.reference}: classes must differ")
        for role, column, measurement, kind in (
            (
                "target",
                pair.target_column,
                pair.target_source_measurement,
                pair.target_source_kind,
            ),
            (
                "reference",
                pair.reference_column,
                pair.reference_source_measurement,
                pair.reference_source_kind,
            ),
        ):
            expected_measurement = SOURCE_MEASUREMENTS.get(column)
            expected_kind = DIAGNOSTIC_SOURCE_KINDS.get(column)
            if expected_measurement is None or expected_kind is None:
                raise ValueError(f"{role} column {column!r} is not in the QuPath input projection")
            if measurement != expected_measurement:
                raise ValueError(
                    f"{role} column {column!r} must use exact QuPath measurement "
                    f"{expected_measurement!r}"
                )
            if kind != expected_kind:
                raise ValueError(
                    f"{role} column {column!r} must declare source kind {expected_kind!r}"
                )
        pairs.append(pair)

    directed = [(pair.target_column, pair.reference_column) for pair in pairs]
    if len(directed) != len(set(directed)):
        raise ValueError("fundamental roster contains duplicate pairs")
    undirected = [frozenset(pair) for pair in directed]
    if len(undirected) != len(set(undirected)):
        raise ValueError("fundamental roster contains a redundant reciprocal pair")
    if len({pair.reference_column for pair in pairs}) != 1:
        raise ValueError("fundamental diagnostics require one shared reference column")
    roster = DiagnosticRoster(
        schema_version=ROSTER_SCHEMA_VERSION,
        roster_version=str(payload.get("roster_version", "")),
        method_status="development_mutable",
        description=str(payload.get("description", "")),
        required_tail_quantiles=quantiles,
        pairs=tuple(pairs),
    )
    _ = roster.source_contract
    return roster


def load_diagnostic_roster(path: str | Path = ROSTER_PATH) -> DiagnosticRoster:
    return roster_from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


ROSTER = load_diagnostic_roster()


def _correlation(left: np.ndarray, right: np.ndarray, *, rank: bool) -> float:
    if len(left) < 2:
        return float("nan")
    if rank:
        left = pd.Series(left).rank(method="average").to_numpy(dtype=float)
        right = pd.Series(right).rank(method="average").to_numpy(dtype=float)
    if np.ptp(left) == 0 or np.ptp(right) == 0:
        return float("nan")
    return float(np.corrcoef(left, right)[0, 1])


def qupath_pair_diagnostics(
    frame: pd.DataFrame,
    *,
    donor: str,
    pair: DiagnosticPair,
) -> dict[str, object]:
    """Summarize one mixed-scale QuPath pair without making a call."""

    missing = sorted({pair.target_column, pair.reference_column}.difference(frame.columns))
    if missing:
        raise ValueError(f"{donor}: QuPath pair is missing columns {missing}")
    target = pd.to_numeric(frame[pair.target_column], errors="coerce").to_numpy(dtype=float)
    reference = pd.to_numeric(frame[pair.reference_column], errors="coerce").to_numpy(dtype=float)
    measured = np.isfinite(target) & np.isfinite(reference)
    target_values = target[measured]
    reference_values = reference[measured]
    n = int(measured.sum())
    return {
        "method_version": METHOD_VERSION,
        "donor": str(donor),
        "target": pair.target,
        "target_column": pair.target_column,
        "target_source_measurement": pair.target_source_measurement,
        "target_source_kind": pair.target_source_kind,
        "target_class": pair.target_class,
        "reference": pair.reference,
        "reference_column": pair.reference_column,
        "reference_source_measurement": pair.reference_source_measurement,
        "reference_source_kind": pair.reference_source_kind,
        "reference_class": pair.reference_class,
        "input_n": int(len(frame)),
        "measured_n": n,
        "unmeasured_n": int(len(frame) - n),
        "target_negative_n": int(np.count_nonzero(target_values < 0)),
        "reference_negative_n": int(np.count_nonzero(reference_values < 0)),
        "target_zero_n": int(np.count_nonzero(target_values == 0)),
        "reference_zero_n": int(np.count_nonzero(reference_values == 0)),
        "pearson": _correlation(target_values, reference_values, rank=False),
        "spearman": _correlation(target_values, reference_values, rank=True),
        "target_median": float(np.median(target_values)) if n else float("nan"),
        "reference_median": float(np.median(reference_values)) if n else float("nan"),
    }


def summarize_input_features(
    frame: pd.DataFrame,
    *,
    donor: str,
    roster: DiagnosticRoster,
) -> pd.DataFrame:
    """Report signed-value and missingness behavior for each unique input."""

    rows: list[dict[str, object]] = []
    for column, (measurement, kind) in roster.source_contract.items():
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        row: dict[str, object] = {
            "donor": str(donor),
            "column": column,
            "source_measurement": measurement,
            "source_kind": kind,
            "input_n": int(len(values)),
            "finite_n": int(len(finite)),
            "missing_n": int(len(values) - len(finite)),
            "negative_n": int(np.count_nonzero(finite < 0)),
            "zero_n": int(np.count_nonzero(finite == 0)),
        }
        for name, quantile in (
            ("minimum", 0.0),
            ("median", 0.5),
            ("q900", 0.9),
            ("q990", 0.99),
            ("maximum", 1.0),
        ):
            row[name] = float(np.quantile(finite, quantile)) if len(finite) else float("nan")
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def load_qupath_donor(
    input_root: str | Path,
    donor: str,
    roster: DiagnosticRoster,
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    """Read only the projected upstream QuPath columns needed by the roster."""

    frame, paths, _summary = load_qupath_projection_donor(
        input_root,
        str(donor),
        roster.source_columns,
    )
    frame = frame[["object_id", *roster.source_columns]].copy()
    for column in roster.source_columns:
        values = frame[column]
        if not values.notna().any():
            raise ValueError(f"{donor}: {column} has no finite measurements")
    return frame, tuple(str(path.resolve()) for path in paths)


def prepare_development_workspace(root: str | Path, label: str) -> Path:
    """Replace exactly one survey-owned mutable workspace."""

    if not label or Path(label).name != label or label in {".", ".."}:
        raise ValueError("label must be one nonempty path component")
    development_root = Path(root).expanduser().resolve()
    dangerous_roots = {Path("/"), Path.home().resolve(), WORKSPACE.resolve(), REPO.resolve()}
    if development_root in dangerous_roots:
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


def _atomic_write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="+", required=True)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--refresh-inputs", action="store_true")
    parser.add_argument("--cohort-spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--development-root", type=Path, default=DEFAULT_DEVELOPMENT_ROOT)
    parser.add_argument("--label", default=DEFAULT_LABEL)
    parser.add_argument("--roster", type=Path, default=ROSTER_PATH)
    args = parser.parse_args(argv)
    donors = tuple(map(str, args.donors))
    if len(set(donors)) != len(donors):
        parser.error("--donors must be unique")
    roster = load_diagnostic_roster(args.roster)
    if args.refresh_inputs:
        try:
            refresh_qupath_projection(args.cohort_spec, args.input_root, donors)
        except ValueError as exc:
            parser.error(str(exc))
    try:
        output = prepare_development_workspace(args.development_root, args.label)
    except ValueError as exc:
        parser.error(str(exc))

    started = time.time()
    pair_rows: list[dict[str, object]] = []
    feature_tables: list[pd.DataFrame] = []
    sources: dict[str, tuple[str, ...]] = {}
    for donor in donors:
        frame, paths = load_qupath_donor(args.input_root, donor, roster)
        sources[donor] = paths
        pair_rows.extend(
            qupath_pair_diagnostics(frame, donor=donor, pair=pair) for pair in roster.pairs
        )
        feature_tables.append(summarize_input_features(frame, donor=donor, roster=roster))
        print(f"[{donor}] {len(frame):,} QuPath cells; {len(roster.pairs)} comparisons")

    pair_table = pd.DataFrame.from_records(pair_rows)
    feature_table = pd.concat(feature_tables, ignore_index=True)
    _atomic_write_parquet(pair_table, output / "fundamental_pair_diagnostics.parquet")
    _atomic_write_parquet(feature_table, output / "fundamental_input_summary.parquet")
    shutil.copyfile(args.roster, output / "roster_used.json")
    metadata = {
        "status": "complete",
        "method_version": METHOD_VERSION,
        "method_status": "development_mutable",
        "input_root": str(Path(args.input_root).expanduser().resolve()),
        "source_partitions": {donor: list(paths) for donor, paths in sources.items()},
        "source_semantics": (
            "upstream QuPath cell measurements; Iba1_norm is exact "
            "Cell: Iba1: Mean (norm); no RESTORE output is read"
        ),
        "roster_version": roster.roster_version,
        "donors": list(donors),
        "pair_n_per_donor": len(roster.pairs),
        "creates_cell_calls": False,
        "restore_outputs_read": False,
        "elapsed_seconds": time.time() - started,
        "output_root": str(output),
        "mutable_workspace": True,
        "production_outputs_written": False,
    }
    (output / "development_summary.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
