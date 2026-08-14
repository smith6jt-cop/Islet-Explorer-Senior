from __future__ import annotations

import json
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.review.qupath_fundamental_inputs import (
    METHOD_VERSION as PROJECTION_METHOD_VERSION,
)
from scripts.review.qupath_fundamental_inputs import (
    SOURCE_MEASUREMENTS,
)
from scripts.review.qupath_fundamental_inputs import (
    WORKSPACE_MARKER as PROJECTION_WORKSPACE_MARKER,
)
from scripts.review.reference_overlap_audit import (
    DEFAULT_ROSTER,
    TAIL_QUANTILES,
    DiagnosticPair,
    audit_pair_scores,
    exact_upper_tail,
    load_diagnostic_roster,
    load_donor_scores,
    main,
    pairwise_overlap_metrics,
    prepare_development_workspace,
)


def test_direct_script_help_loads_from_repository_root():
    script = Path(__file__).parents[2] / "scripts" / "review" / "reference_overlap_audit.py"
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=script.parents[2],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr


def test_roster_is_exactly_two_upstream_qupath_comparisons_without_sparse_cd79a():
    roster = load_diagnostic_roster(DEFAULT_ROSTER)

    assert roster.schema_version == 3
    assert roster.method_status == "development_mutable"
    assert roster.required_tail_quantiles == TAIL_QUANTILES
    assert [pair.target_column for pair in roster.pairs] == ["CD3e", "Iba1_norm"]
    assert {pair.reference_column for pair in roster.pairs} == {"Ker8_18"}
    assert roster.source_columns == ("CD3e", "Ker8_18", "Iba1_norm")
    assert roster.pairs[-1].target_source_measurement == "Cell: Iba1: Mean (norm)"
    assert not hasattr(roster, "normalization_threshold")


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda value: value.update(schema_version=2), "schema_version"),
        (lambda value: value.update(method_status="review_only"), "development_mutable"),
        (
            lambda value: value.update(required_tail_quantiles=[0.99, 0.9]),
            "unique, increasing",
        ),
        (
            lambda value: value["pairs"][1].update(target_source_measurement="Cell: Iba1: Mean"),
            r"Cell: Iba1: Mean \(norm\)",
        ),
    ],
)
def test_roster_rejects_contract_drift(tmp_path, mutation, message):
    payload = json.loads(DEFAULT_ROSTER.read_text(encoding="utf-8"))
    mutation(payload)
    path = tmp_path / "roster.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        load_diagnostic_roster(path)


def test_roster_rejects_reciprocal_and_accepts_reduction(tmp_path):
    payload = json.loads(DEFAULT_ROSTER.read_text(encoding="utf-8"))
    reciprocal = deepcopy(payload["pairs"][0])
    for stem in ("", "_column", "_source_measurement", "_source_kind", "_class"):
        target = f"target{stem}"
        reference = f"reference{stem}"
        reciprocal[target], reciprocal[reference] = reciprocal[reference], reciprocal[target]
    payload["pairs"].append(reciprocal)
    path = tmp_path / "reciprocal.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="redundant reciprocal"):
        load_diagnostic_roster(path)

    payload = json.loads(DEFAULT_ROSTER.read_text(encoding="utf-8"))
    payload["pairs"] = payload["pairs"][:2]
    payload["required_tail_quantiles"] = [0.95, 0.99]
    path = tmp_path / "edited.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    edited = load_diagnostic_roster(path)
    assert len(edited.pairs) == 2
    assert edited.required_tail_quantiles == (0.95, 0.99)


def test_exact_upper_tail_has_decimal_size_and_row_invariant_id_ties():
    identifiers = np.array([f"id-{index:03d}" for index in range(100)])
    tied = np.ones(100)
    mask, audit = exact_upper_tail(tied, identifiers, 0.95)

    assert mask.sum() == 5
    assert set(identifiers[mask]) == set(sorted(identifiers)[:5])
    assert audit == {
        "quantile": 0.95,
        "tail_n": 5,
        "boundary_score": 1.0,
        "boundary_tie_n": 100,
        "boundary_tie_selected_n": 5,
        "tie_break": "object_id_ascending",
    }
    permutation = np.random.default_rng(7).permutation(100)
    shuffled, _ = exact_upper_tail(tied[permutation], identifiers[permutation], 0.95)
    assert set(identifiers[permutation][shuffled]) == set(identifiers[mask])


def test_pairwise_overlap_reports_symmetric_independence_metrics():
    target = np.array([1, 1, 0, 0, 0, 0, 0, 0], dtype=bool)
    reference = np.array([1, 0, 1, 0, 0, 0, 0, 0], dtype=bool)
    metrics = pairwise_overlap_metrics(target, reference, prefix="tail")

    assert metrics["both_tail_n"] == 1
    assert metrics["both_of_target_tail"] == pytest.approx(0.5)
    assert metrics["both_of_reference_tail"] == pytest.approx(0.5)
    assert metrics["tail_jaccard"] == pytest.approx(1 / 3)
    assert metrics["tail_independence_lift"] == pytest.approx(2.0)


def _pair() -> DiagnosticPair:
    return load_diagnostic_roster(DEFAULT_ROSTER).pairs[-1]


def test_pair_audit_uses_signed_values_and_rank_grid_without_calls():
    n = 1_000
    identifiers = np.array([f"cell-{index:04d}" for index in range(n)])
    target = np.linspace(-0.5, 2.0, n)
    reference = np.roll(np.linspace(0, 5_000, n), 100)
    target[4] = np.nan

    rows = audit_pair_scores(
        donor="D",
        pair=_pair(),
        object_ids=identifiers,
        target_values=target,
        reference_values=reference,
    )

    assert [row["quantile"] for row in rows] == list(TAIL_QUANTILES)
    assert all(row["finite_pair_n"] == 999 for row in rows)
    assert all(row["excluded_nonfinite_n"] == 1 for row in rows)
    assert all(row["target_negative_n"] > 0 for row in rows)
    assert rows[-1]["target_tail_n"] == rows[-1]["reference_tail_n"] == 1
    assert all("tail_independence_lift" in row for row in rows)
    assert all(
        not any(fragment in key for fragment in ("positive", "threshold"))
        for row in rows
        for key in row
    )
    assert all("pearson_score" not in row for row in rows)


def test_rank_overlap_is_invariant_to_nonlinear_monotone_partner_transform():
    n = 500
    identifiers = np.array([f"cell-{index:04d}" for index in range(n)])
    target = np.sin(np.linspace(0, 8, n))
    reference = np.roll(np.arange(1, n + 1, dtype=float), 73)
    baseline = audit_pair_scores(
        donor="D",
        pair=_pair(),
        object_ids=identifiers,
        target_values=target,
        reference_values=reference,
    )
    transformed = audit_pair_scores(
        donor="D",
        pair=_pair(),
        object_ids=identifiers,
        target_values=target,
        reference_values=reference**2,
    )

    for original, changed in zip(baseline, transformed, strict=True):
        for field in (
            "reference_tail_n",
            "both_tail_n",
            "tail_jaccard",
            "tail_independence_lift",
            "spearman_score",
        ):
            assert changed[field] == pytest.approx(original[field], nan_ok=True)


@pytest.mark.parametrize("dtype", [np.uint16, np.int32, np.float32, np.float64])
@pytest.mark.parametrize(("scale", "offset"), [(1.0, 0.0), (1e-3, 7.0), (1e3, 11.0)])
def test_rank_overlap_is_invariant_to_partner_dtype_and_positive_affine_scale(dtype, scale, offset):
    n = 1_000
    identifiers = np.array([f"cell-{index:04d}" for index in range(n)])
    target = np.linspace(-0.5, 2.0, n)
    raw_partner = np.roll(np.arange(1, n + 1, dtype=np.uint16), 137)
    baseline = audit_pair_scores(
        donor="D",
        pair=_pair(),
        object_ids=identifiers,
        target_values=target,
        reference_values=raw_partner,
    )
    transformed = audit_pair_scores(
        donor="D",
        pair=_pair(),
        object_ids=identifiers,
        target_values=target,
        reference_values=raw_partner.astype(dtype) * scale + offset,
    )

    invariant_fields = {
        "target_tail_n",
        "reference_tail_n",
        "both_tail_n",
        "both_of_target_tail",
        "both_of_reference_tail",
        "tail_jaccard",
        "tail_independence_expected_n",
        "tail_independence_lift",
        "spearman_score",
    }
    for original, changed in zip(baseline, transformed, strict=True):
        for field in invariant_fields:
            assert changed[field] == pytest.approx(original[field], nan_ok=True)
        assert changed["reference_tail_boundary_score"] == pytest.approx(
            original["reference_tail_boundary_score"] * scale + offset
        )


def _write_parts(root, donor: str, frame: pd.DataFrame) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / PROJECTION_WORKSPACE_MARKER).write_text(
        PROJECTION_METHOD_VERSION + "\n", encoding="utf-8"
    )
    (root / "projection_summary.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "method_version": PROJECTION_METHOD_VERSION,
                "source_kind": "upstream QuPath full-cell measurement export",
                "source_measurements": SOURCE_MEASUREMENTS,
                "restore_outputs_read": False,
                "donors": {str(donor): len(frame)},
                "donor_images": {str(donor): "image"},
            }
        ),
        encoding="utf-8",
    )
    partition = root / f"donor_id={donor}"
    partition.mkdir(parents=True, exist_ok=True)
    midpoint = len(frame) // 2
    frame.iloc[:midpoint].to_parquet(partition / "part-0.parquet", index=False)
    frame.iloc[midpoint:].to_parquet(partition / "part-1.parquet", index=False)


def _input_frame(n: int = 100) -> pd.DataFrame:
    values = np.linspace(0.0, 2_000.0, n)
    return pd.DataFrame(
        {
            "donor_id": ["D"] * n,
            "image": ["image"] * n,
            "object_id": [f"cell-{index:04d}" for index in range(n)],
            "CD3e": values,
            "Iba1_norm": np.roll(np.linspace(-0.5, 1.5, n), 20),
            "Ker8_18": values[::-1],
            "Iba1": 9999.0,
            "Iba1__restore_normalized": 8888.0,
        }
    )


def test_donor_loader_accepts_parts_and_reads_only_qupath_contract(tmp_path):
    roster = load_diagnostic_roster(DEFAULT_ROSTER)
    frame = _input_frame()
    _write_parts(tmp_path, "D", frame)

    loaded = load_donor_scores(tmp_path, "D", roster)

    assert len(loaded) == len(frame)
    assert tuple(loaded.columns) == ("object_id", *roster.source_columns)
    assert "Iba1" not in loaded
    assert "Iba1__restore_normalized" not in loaded

    broken_root = tmp_path / "broken"
    _write_parts(broken_root, "D", frame.drop(columns="Iba1_norm"))
    with pytest.raises(ValueError, match="Iba1_norm"):
        load_donor_scores(broken_root, "D", roster)


def test_same_label_refresh_replaces_one_rank_only_table(tmp_path):
    input_root = tmp_path / "inputs"
    output_root = tmp_path / "development"
    frame = _input_frame()
    _write_parts(input_root, "D", frame)
    argv = [
        "--donors",
        "D",
        "--label",
        "scratch",
        "--input-root",
        str(input_root),
        "--development-root",
        str(output_root),
    ]

    assert main(argv) == 0
    path = output_root / "scratch" / "reference_overlap_pairs.parquet"
    first = pd.read_parquet(path)
    assert len(first) == len(load_diagnostic_roster(DEFAULT_ROSTER).pairs) * len(TAIL_QUANTILES)
    assert not any(
        fragment in column
        for column in first.columns
        for fragment in ("positive", "threshold", "restore")
    )

    changed = frame.copy()
    changed["CD3e"] = changed["Ker8_18"]
    _write_parts(input_root, "D", changed)
    assert main(argv) == 0
    refreshed = pd.read_parquet(path)
    before = first.loc[first["target"].eq("CD3e") & first["quantile"].eq(0.9), "both_tail_n"].iloc[
        0
    ]
    after = refreshed.loc[
        refreshed["target"].eq("CD3e") & refreshed["quantile"].eq(0.9), "both_tail_n"
    ].iloc[0]
    assert after > before


def test_workspace_refuses_to_replace_unowned_directory(tmp_path):
    root = tmp_path / "development"
    unowned = root / "diagnostic"
    unowned.mkdir(parents=True)
    (unowned / "user.txt").write_text("keep", encoding="utf-8")

    with pytest.raises(ValueError, match="not owned"):
        prepare_development_workspace(root, "diagnostic")

    assert (unowned / "user.txt").read_text(encoding="utf-8") == "keep"
