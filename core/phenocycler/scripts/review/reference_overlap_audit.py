#!/usr/bin/env python3
"""Refresh rank-only overlap diagnostics for upstream QuPath inputs.

The two comparisons intentionally mix raw QuPath cell means with QuPath's
signed acquisition-normalized Iba1 measurement.  Exact donor-local ranks are
comparable across those scales; a shared numeric positivity threshold is not.
No RESTORE table, divisor, call, classification, phenotype, or disease label is
read by this tool.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections.abc import Sequence
from decimal import ROUND_CEILING, Decimal
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from scripts.review.fundamental_pair_survey import (  # noqa: E402
    DEFAULT_INPUT_ROOT,
    DiagnosticPair,
    DiagnosticRoster,
    load_diagnostic_roster,
    load_qupath_donor,
)
from scripts.review.qupath_fundamental_inputs import (  # noqa: E402
    DEFAULT_SPEC,
    refresh_qupath_projection,
)

DEFAULT_ROSTER = REPO / "phenocycler" / "fundamental_pair_diagnostics.json"
DEFAULT_OUTPUT_ROOT = REPO.parents[1] / "data" / "development"
METHOD_VERSION = "qupath-input-rank-overlap-development-v4"
WORKSPACE_MARKER = ".reference_overlap_workspace"
TAIL_QUANTILES = (0.9, 0.95, 0.99, 0.995, 0.999)


def _validate_vectors(
    values: np.ndarray,
    object_ids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    numeric = np.asarray(values, dtype=float)
    identifiers = np.asarray(object_ids).astype(str)
    if numeric.ndim != 1 or identifiers.ndim != 1 or numeric.shape != identifiers.shape:
        raise ValueError("values and object_ids must be aligned one-dimensional arrays")
    if not len(numeric) or not np.isfinite(numeric).all():
        raise ValueError("rank values must be nonempty and finite")
    if len(np.unique(identifiers)) != len(identifiers):
        raise ValueError("object_ids must be unique")
    return numeric, identifiers


def stable_descending_ranks(values: np.ndarray, object_ids: np.ndarray) -> np.ndarray:
    """Rank descending scores, resolving equal scores by ascending object ID."""

    numeric, identifiers = _validate_vectors(values, object_ids)
    order = np.lexsort((identifiers, -numeric))
    ranks = np.empty(len(order), dtype=np.int64)
    ranks[order] = np.arange(1, len(order) + 1, dtype=np.int64)
    return ranks


def exact_upper_tail(
    values: np.ndarray,
    object_ids: np.ndarray,
    quantile: float,
) -> tuple[np.ndarray, dict[str, object]]:
    """Select exactly ``ceil((1-q)N)`` rows with deterministic score ties."""

    numeric, identifiers = _validate_vectors(values, object_ids)
    ranks = stable_descending_ranks(numeric, identifiers)
    return _upper_tail_from_ranks(numeric, ranks, quantile)


def _upper_tail_from_ranks(
    values: np.ndarray,
    ranks: np.ndarray,
    quantile: float,
) -> tuple[np.ndarray, dict[str, object]]:
    numeric = np.asarray(values, dtype=float)
    stable_ranks = np.asarray(ranks, dtype=np.int64)
    if numeric.ndim != 1 or numeric.shape != stable_ranks.shape:
        raise ValueError("values and ranks must be aligned one-dimensional arrays")
    q = float(quantile)
    if not len(numeric) or not 0 < q < 1:
        raise ValueError("rank tail requires rows and a quantile between zero and one")
    k = int(
        ((Decimal("1") - Decimal(str(q))) * Decimal(len(numeric))).to_integral_value(
            rounding=ROUND_CEILING
        )
    )
    k = min(len(numeric), max(1, k))
    mask = stable_ranks <= k
    boundary = float(numeric[stable_ranks == k][0])
    return mask, {
        "quantile": q,
        "tail_n": k,
        "boundary_score": boundary,
        "boundary_tie_n": int(np.count_nonzero(numeric == boundary)),
        "boundary_tie_selected_n": int(np.count_nonzero(mask & (numeric == boundary))),
        "tie_break": "object_id_ascending",
    }


def _safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator > 0 else float("nan")


def pairwise_overlap_metrics(
    target_mask: np.ndarray,
    reference_mask: np.ndarray,
    *,
    prefix: str,
) -> dict[str, object]:
    """Return symmetric overlap, Jaccard, and independence-lift statistics."""

    target = np.asarray(target_mask, dtype=bool)
    reference = np.asarray(reference_mask, dtype=bool)
    if target.ndim != 1 or target.shape != reference.shape:
        raise ValueError("overlap masks must be aligned one-dimensional arrays")
    n = len(target)
    target_n = int(target.sum())
    reference_n = int(reference.sum())
    both_n = int(np.count_nonzero(target & reference))
    union_n = int(np.count_nonzero(target | reference))
    expected_n = target_n * reference_n / n if n else float("nan")
    return {
        f"target_{prefix}_n": target_n,
        f"reference_{prefix}_n": reference_n,
        f"both_{prefix}_n": both_n,
        f"both_of_target_{prefix}": _safe_ratio(both_n, target_n),
        f"both_of_reference_{prefix}": _safe_ratio(both_n, reference_n),
        f"{prefix}_jaccard": _safe_ratio(both_n, union_n),
        f"{prefix}_independence_expected_n": float(expected_n),
        f"{prefix}_independence_lift": _safe_ratio(both_n, expected_n),
    }


def _rank_correlation(target: np.ndarray, reference: np.ndarray) -> float:
    if np.ptp(target) == 0 or np.ptp(reference) == 0:
        return float("nan")
    return float(pd.Series(target).corr(pd.Series(reference), method="spearman"))


def audit_pair_scores(
    *,
    donor: str,
    pair: DiagnosticPair,
    object_ids: np.ndarray,
    target_values: np.ndarray,
    reference_values: np.ndarray,
    quantiles: Sequence[float] = TAIL_QUANTILES,
) -> list[dict[str, object]]:
    """Audit one mixed-scale pair using only deterministic donor-local ranks."""

    identifiers = np.asarray(object_ids).astype(str)
    target = np.asarray(target_values, dtype=float)
    reference = np.asarray(reference_values, dtype=float)
    if not (identifiers.shape == target.shape == reference.shape) or target.ndim != 1:
        raise ValueError("pair scores and object_ids must be aligned one-dimensional arrays")
    finite = np.isfinite(target) & np.isfinite(reference)
    if int(finite.sum()) < 2:
        raise ValueError(f"{donor} {pair.target} vs {pair.reference}: fewer than two finite rows")
    identifiers = identifiers[finite]
    target = target[finite]
    reference = reference[finite]
    if len(np.unique(identifiers)) != len(identifiers):
        raise ValueError(f"donor {donor} contains duplicate object_id")

    spearman = _rank_correlation(target, reference)
    target_ranks = stable_descending_ranks(target, identifiers)
    reference_ranks = stable_descending_ranks(reference, identifiers)
    pair_fields: dict[str, object] = {
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
        "input_n": int(len(finite)),
        "finite_pair_n": int(finite.sum()),
        "excluded_nonfinite_n": int((~finite).sum()),
        "target_negative_n": int(np.count_nonzero(target < 0)),
        "reference_negative_n": int(np.count_nonzero(reference < 0)),
        "spearman_score": spearman,
    }
    rows: list[dict[str, object]] = []
    for quantile in quantiles:
        target_tail, target_audit = _upper_tail_from_ranks(target, target_ranks, quantile)
        reference_tail, reference_audit = _upper_tail_from_ranks(
            reference, reference_ranks, quantile
        )
        rows.append(
            {
                **pair_fields,
                "quantile": float(quantile),
                "target_tail_boundary_score": target_audit["boundary_score"],
                "target_boundary_tie_n": target_audit["boundary_tie_n"],
                "target_boundary_tie_selected_n": target_audit["boundary_tie_selected_n"],
                "reference_tail_boundary_score": reference_audit["boundary_score"],
                "reference_boundary_tie_n": reference_audit["boundary_tie_n"],
                "reference_boundary_tie_selected_n": reference_audit["boundary_tie_selected_n"],
                "tail_tie_break": "object_id_ascending",
                **pairwise_overlap_metrics(target_tail, reference_tail, prefix="tail"),
            }
        )
    return rows


def load_donor_scores(
    input_root: Path | str,
    donor: str,
    roster: DiagnosticRoster,
) -> pd.DataFrame:
    """Read the shared upstream QuPath projection for one donor."""

    frame, _paths = load_qupath_donor(input_root, donor, roster)
    return frame


def _atomic_write_parquet(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_parquet(temporary, index=False)
    temporary.replace(path)


def prepare_development_workspace(root: Path | str, label: str) -> Path:
    """Replace only a workspace previously created by this development tool."""

    if not label or Path(label).name != label or label in {".", ".."}:
        raise ValueError("label must be one nonempty path component")
    development_root = Path(root).expanduser().resolve()
    dangerous_roots = {Path("/"), Path.home().resolve(), REPO.resolve(), REPO.parent.resolve()}
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


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="+", required=True)
    parser.add_argument("--label", default="fundamental_pair_overlap")
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--refresh-inputs", action="store_true")
    parser.add_argument("--cohort-spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument(
        "--development-root",
        "--output-root",
        dest="development_root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
    )
    parser.add_argument("--roster", type=Path, default=DEFAULT_ROSTER)
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

    rows: list[dict[str, object]] = []
    for donor in donors:
        frame = load_donor_scores(args.input_root, donor, roster)
        for pair in roster.pairs:
            rows.extend(
                audit_pair_scores(
                    donor=donor,
                    pair=pair,
                    object_ids=frame["object_id"].to_numpy(str),
                    target_values=frame[pair.target_column].to_numpy(float),
                    reference_values=frame[pair.reference_column].to_numpy(float),
                    quantiles=roster.required_tail_quantiles,
                )
            )
        print(f"[{donor}] {len(frame):,} rows; {len(roster.pairs)} QuPath input pairs")

    try:
        output = prepare_development_workspace(args.development_root, args.label)
    except ValueError as exc:
        parser.error(str(exc))
    table = pd.DataFrame.from_records(rows)
    path = output / "reference_overlap_pairs.parquet"
    _atomic_write_parquet(table, path)
    print(f"refreshed {path} ({len(table)} pair-tail rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
