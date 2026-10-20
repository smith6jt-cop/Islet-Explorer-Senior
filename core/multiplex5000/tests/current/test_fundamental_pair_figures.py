from __future__ import annotations

import importlib.util
import inspect
import json
import struct
import subprocess
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scripts.review import fundamental_pair_figures as figures
from scripts.review.fundamental_pair_survey import METHOD_VERSION as SURVEY_METHOD_VERSION
from scripts.review.fundamental_pair_survey import ROSTER
from scripts.review.qupath_fundamental_inputs import (
    METHOD_VERSION as PROJECTION_METHOD_VERSION,
)
from scripts.review.qupath_fundamental_inputs import SOURCE_MEASUREMENTS
from scripts.review.reference_overlap_audit import METHOD_VERSION as OVERLAP_METHOD_VERSION
from scripts.review.reference_overlap_audit import exact_upper_tail


def _png_chunk(kind: bytes, data: bytes) -> bytes:
    return (
        struct.pack(">I", len(data))
        + kind
        + data
        + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)
    )


def _write_test_png(path: Path, *, width: int = 80, height: int = 20) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for y in range(height):
        pixels = b"".join(
            bytes(((x * 7) % 256, (y * 13) % 256, ((x + y) * 5) % 256)) for x in range(width)
        )
        rows.append(b"\x00" + pixels)
    payload = (
        b"\x89PNG\r\n\x1a\n"
        + _png_chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + _png_chunk(b"IDAT", zlib.compress(b"".join(rows)))
        + _png_chunk(b"IEND", b"")
    )
    path.write_bytes(payload)


def _write_review_artifacts(
    root: Path,
    donors: tuple[str, ...] = ("D1", "D2"),
    *,
    n: int = 120,
):
    projection_root = root / "projection"
    survey_root = root / "survey"
    overlap_root = root / "overlap"
    crop_root = root / "crops"
    for path in (projection_root, survey_root, overlap_root, crop_root):
        path.mkdir(parents=True)

    projection_donors = tuple(dict.fromkeys((*donors, figures.DEFAULT_CROP_DONOR)))
    (projection_root / ".qupath_fundamental_inputs_workspace").write_text(
        PROJECTION_METHOD_VERSION + "\n", encoding="utf-8"
    )
    (projection_root / figures.PROJECTION_METADATA_NAME).write_text(
        json.dumps(
            {
                "status": "complete",
                "method_version": PROJECTION_METHOD_VERSION,
                "source_kind": "upstream QuPath full-cell measurement export",
                "source_measurements": SOURCE_MEASUREMENTS,
                "restore_outputs_read": False,
                "donors": {donor: n for donor in projection_donors},
                "donor_images": {donor: f"image-{donor}" for donor in projection_donors},
            }
        ),
        encoding="utf-8",
    )
    base = np.linspace(0.0, 1.0, n)
    for donor_index, donor in enumerate(projection_donors):
        partition = projection_root / f"donor_id={donor}"
        partition.mkdir()
        pd.DataFrame(
            {
                "donor_id": [donor] * n,
                "image": [f"image-{donor}"] * n,
                "object_id": [f"{donor}-cell-{index:04d}" for index in range(n)],
                "CD3e": np.roll(base, donor_index + 3) * 1_000,
                "CD79a": np.roll(base, donor_index + 17) * 800,
                "Iba1_norm": np.roll(base, donor_index + 29) * 2 - 0.5,
                "Ker8_18": base[::-1] * 4_000,
            }
        ).to_parquet(partition / "data.parquet", index=False)

    (survey_root / figures.SURVEY_METADATA_NAME).write_text(
        json.dumps(
            {
                "status": "complete",
                "method_version": SURVEY_METHOD_VERSION,
                "method_status": "development_mutable",
                "roster_version": ROSTER.roster_version,
                "input_root": str(projection_root.resolve()),
                "restore_outputs_read": False,
                "creates_cell_calls": False,
                "production_outputs_written": False,
            }
        ),
        encoding="utf-8",
    )
    (survey_root / ".fundamental_pair_survey_workspace").write_text(
        SURVEY_METHOD_VERSION + "\n", encoding="utf-8"
    )
    (overlap_root / ".reference_overlap_workspace").write_text(
        OVERLAP_METHOD_VERSION + "\n", encoding="utf-8"
    )

    pair_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    for donor_index, donor in enumerate(donors):
        object_ids = np.array([f"{donor}-cell-{index:04d}" for index in range(n)])
        reference_values = base[::-1] * 4_000
        for pair_index, pair in enumerate(ROSTER.pairs):
            target_values = {
                "CD3e": np.roll(base, donor_index + 3) * 1_000,
                "CD79a": np.roll(base, donor_index + 17) * 800,
                "Iba1_norm": np.roll(base, donor_index + 29) * 2 - 0.5,
            }[pair.target_column]
            pair_rows.append(
                {
                    "method_version": SURVEY_METHOD_VERSION,
                    "donor": donor,
                    "target": pair.target,
                    "target_column": pair.target_column,
                    "target_source_measurement": pair.target_source_measurement,
                    "target_source_kind": pair.target_source_kind,
                    "reference": pair.reference,
                    "reference_column": pair.reference_column,
                    "reference_source_measurement": pair.reference_source_measurement,
                    "reference_source_kind": pair.reference_source_kind,
                    "input_n": n,
                    "measured_n": n,
                    "unmeasured_n": 0,
                    "pearson": -0.2 + 0.15 * pair_index + 0.02 * donor_index,
                    "spearman": -0.1 + 0.2 * pair_index + 0.02 * donor_index,
                    "target_median": 0.02 if pair.target_column == "Iba1_norm" else 100.0,
                    "reference_median": 900.0,
                }
            )
            for quantile in ROSTER.required_tail_quantiles:
                target_tail, target_audit = exact_upper_tail(target_values, object_ids, quantile)
                reference_tail, reference_audit = exact_upper_tail(
                    reference_values, object_ids, quantile
                )
                target_n = int(target_tail.sum())
                reference_n = int(reference_tail.sum())
                both_n = int(np.count_nonzero(target_tail & reference_tail))
                share = both_n / target_n
                expected_n = target_n * reference_n / n
                overlap_rows.append(
                    {
                        "method_version": OVERLAP_METHOD_VERSION,
                        "donor": donor,
                        "target": pair.target,
                        "target_column": pair.target_column,
                        "target_source_measurement": pair.target_source_measurement,
                        "target_source_kind": pair.target_source_kind,
                        "reference": pair.reference,
                        "reference_column": pair.reference_column,
                        "reference_source_measurement": pair.reference_source_measurement,
                        "reference_source_kind": pair.reference_source_kind,
                        "finite_pair_n": n,
                        "excluded_nonfinite_n": 0,
                        "spearman_score": -0.1 + 0.2 * pair_index,
                        "quantile": quantile,
                        "target_tail_boundary_score": target_audit["boundary_score"],
                        "reference_tail_boundary_score": reference_audit["boundary_score"],
                        "target_tail_n": target_n,
                        "reference_tail_n": reference_n,
                        "both_tail_n": both_n,
                        "both_of_target_tail": share,
                        "both_of_reference_tail": share,
                        "tail_jaccard": both_n / (target_n + reference_n - both_n),
                        "tail_independence_lift": both_n / expected_n,
                    }
                )
        for column_index, (column, (measurement, kind)) in enumerate(
            ROSTER.source_contract.items()
        ):
            normalized = kind == "qupath_acquisition_normalized_cell_mean"
            values = (
                (-0.5, 0.01, 0.8, 1.4, 1.5)
                if normalized
                else (
                    0.0,
                    100.0 * (column_index + 1),
                    300.0 * (column_index + 1),
                    900.0 * (column_index + 1),
                    4_000.0 * (column_index + 1),
                )
            )
            feature_rows.append(
                {
                    "donor": donor,
                    "column": column,
                    "source_measurement": measurement,
                    "source_kind": kind,
                    "input_n": n,
                    "finite_n": n,
                    "missing_n": 0,
                    "negative_n": 30 if normalized else 0,
                    "zero_n": 1,
                    "minimum": values[0],
                    "median": values[1],
                    "q900": values[2],
                    "q990": values[3],
                    "maximum": values[4],
                }
            )

    pd.DataFrame(pair_rows).to_parquet(survey_root / figures.PAIR_TABLE_NAME, index=False)
    pd.DataFrame(feature_rows).to_parquet(survey_root / figures.FEATURE_TABLE_NAME, index=False)
    pd.DataFrame(overlap_rows).to_parquet(overlap_root / figures.OVERLAP_TABLE_NAME, index=False)

    crop_rows: list[dict[str, object]] = []
    contact_rows: list[dict[str, object]] = []
    crop_directory = crop_root / "crops" / f"donor_id={figures.DEFAULT_CROP_DONOR}"
    for pair in ROSTER.pairs:
        for stratum in ("dual_bright", "target_bright", "reference_bright"):
            for rank in range(1, 5):
                object_id = f"{pair.target_column}-{stratum}-{rank}"
                render_path = crop_directory / f"{object_id}.png"
                _write_test_png(render_path)
                if stratum == "target_bright":
                    target_rank, reference_rank = 0.999, 0.005
                elif stratum == "reference_bright":
                    target_rank, reference_rank = 0.005, 0.999
                else:
                    target_rank, reference_rank = 0.99, 0.995
                crop_rows.append(
                    {
                        "donor_id": figures.DEFAULT_CROP_DONOR,
                        "pair_id": f"{pair.target_column}:{pair.reference_column}",
                        "target": pair.target_column,
                        "reference": pair.reference_column,
                        "target_selection_source_measurement": pair.target_source_measurement,
                        "reference_selection_source_measurement": pair.reference_source_measurement,
                        "target_selection_source_kind": pair.target_source_kind,
                        "reference_selection_source_kind": pair.reference_source_kind,
                        "target_display_marker": (
                            "Iba1" if pair.target_column == "Iba1_norm" else pair.target_column
                        ),
                        "reference_display_marker": pair.reference_column,
                        "stratum": stratum,
                        "stratum_requested_n": 4,
                        "stratum_selected_n": 4,
                        "crop_rank": rank,
                        "object_id": object_id,
                        "selection_target_cell_value": float(rank),
                        "selection_reference_cell_value": float(5 - rank),
                        "target_percentile_rank": target_rank,
                        "reference_percentile_rank": reference_rank,
                        "render_path": str(render_path.resolve()),
                    }
                )
                contact_rows.append(
                    {
                        "object_id": object_id,
                        "contact_neighbor_n": rank + 1,
                        "contact_boundary_px": 20.0 * rank,
                    }
                )
    crop_manifest_path = crop_root / figures.CROP_MANIFEST_TEMPLATE.format(
        donor=figures.DEFAULT_CROP_DONOR
    )
    crop_manifest_path.parent.mkdir(parents=True)
    pd.DataFrame(crop_rows).to_parquet(crop_manifest_path, index=False)
    contact_path = crop_root / figures.CONTACT_METRICS_TEMPLATE.format(
        donor=figures.DEFAULT_CROP_DONOR
    )
    contact_path.parent.mkdir(parents=True)
    pd.DataFrame(contact_rows).to_parquet(contact_path, index=False)
    (crop_root / figures.CROP_METADATA_NAME).write_text(
        json.dumps(
            {
                "method": figures.BOUNDARY_METHOD_VERSION,
                "review_only": True,
                "qupath_input_projection_root": str(projection_root.resolve()),
                "display": {
                    "source": "raw qptiff channels",
                    "selection_values_used_for_display": False,
                },
                "selection": {"restore_outputs_read": False},
                "selection_source_measurements": {
                    column: measurement
                    for column, (measurement, _kind) in ROSTER.source_contract.items()
                },
                "selection_source_kinds": {
                    column: kind for column, (_measurement, kind) in ROSTER.source_contract.items()
                },
                "pairs": [
                    {"target": pair.target_column, "reference": pair.reference_column}
                    for pair in ROSTER.pairs
                ],
                "donors": [
                    {
                        "donor_id": figures.DEFAULT_CROP_DONOR,
                        "crop_n": len(ROSTER.pairs) * len(figures.CROP_STRATA) * 4,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    (crop_root / figures.BOUNDARY_WORKSPACE_MARKER).write_text(
        json.dumps({"tool": "boundary_contact_audit", "status": "complete"}) + "\n",
        encoding="utf-8",
    )
    return projection_root, survey_root, overlap_root, crop_root


def _load(roots):
    return figures.load_figure_inputs(*roots, roster=ROSTER)


def test_direct_script_help_loads_from_repository_root():
    script = Path(__file__).parents[2] / "scripts" / "review" / "fundamental_pair_figures.py"
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=script.parents[2],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "--projection-root" in result.stdout
    assert "--crop-root" in result.stdout
    assert "--refresh-inputs" not in result.stdout


def test_loader_gates_current_projection_summaries_and_raw_crops(tmp_path):
    roots = _write_review_artifacts(tmp_path)

    loaded = _load(roots)

    assert loaded.donors == ("D1", "D2")
    assert loaded.crop_donor == "6450"
    assert loaded.pair_diagnostics["target_column"].tolist() == [
        pair.target_column for _donor in loaded.donors for pair in ROSTER.pairs
    ]
    assert "CD79a" not in set(loaded.pair_diagnostics["target_column"])
    assert len(loaded.crop_manifest) == len(ROSTER.pairs) * len(figures.CROP_STRATA) * 4
    assert (
        loaded.crop_manifest.loc[
            loaded.crop_manifest["target"].eq("Iba1_norm"), "target_display_marker"
        ]
        .eq("Iba1")
        .all()
    )
    assert set(loaded.source_method_versions) == {
        "projection",
        "survey",
        "rank_overlap",
        "boundary_crops",
    }
    assert not any("restore" in column.lower() for column in loaded.rank_overlap.columns)


@pytest.mark.parametrize("failure", ["metadata", "source", "stale", "stale_marker", "stale_crop"])
def test_loader_refuses_restore_or_stale_inputs(tmp_path, failure):
    roots = _write_review_artifacts(tmp_path)
    _projection_root, survey_root, overlap_root, crop_root = roots
    if failure == "metadata":
        metadata = json.loads(
            (survey_root / figures.SURVEY_METADATA_NAME).read_text(encoding="utf-8")
        )
        metadata["restore_outputs_read"] = True
        (survey_root / figures.SURVEY_METADATA_NAME).write_text(
            json.dumps(metadata), encoding="utf-8"
        )
        message = "stale or unsafe"
    elif failure == "source":
        table = pd.read_parquet(overlap_root / figures.OVERLAP_TABLE_NAME)
        table.loc[table["target_column"].eq("Iba1_norm"), "target_source_kind"] = (
            "restore_normalized"
        )
        table.to_parquet(overlap_root / figures.OVERLAP_TABLE_NAME, index=False)
        message = "RESTORE aliases"
    elif failure == "stale":
        table = pd.read_parquet(overlap_root / figures.OVERLAP_TABLE_NAME)
        table["method_version"] = "qupath-input-rank-overlap-development-v3"
        table.to_parquet(overlap_root / figures.OVERLAP_TABLE_NAME, index=False)
        message = "stale"
    elif failure == "stale_marker":
        (overlap_root / ".reference_overlap_workspace").write_text(
            "qupath-input-rank-overlap-development-v3\n", encoding="utf-8"
        )
        message = "workspace is stale"
    else:
        metadata = json.loads((crop_root / figures.CROP_METADATA_NAME).read_text(encoding="utf-8"))
        metadata["method"] = "boundary-contact-audit-development-v3"
        (crop_root / figures.CROP_METADATA_NAME).write_text(json.dumps(metadata), encoding="utf-8")
        message = "boundary-crop artifact is stale"

    with pytest.raises(ValueError, match=message):
        _load(roots)


def test_loader_refuses_incomplete_tail_grid(tmp_path):
    roots = _write_review_artifacts(tmp_path)
    overlap_root = roots[2]
    table = pd.read_parquet(overlap_root / figures.OVERLAP_TABLE_NAME)
    table = table.loc[
        ~(
            table["donor"].eq("D1")
            & table["target_column"].eq("Iba1_norm")
            & table["quantile"].eq(0.999)
        )
    ]
    table.to_parquet(overlap_root / figures.OVERLAP_TABLE_NAME, index=False)

    with pytest.raises(ValueError, match="complete tail grid"):
        _load(roots)


def test_native_joint_density_filters_without_clipping_and_is_bounded():
    target = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, np.nan])
    reference = np.arange(6, dtype=float)
    baseline = figures.native_joint_histogram(
        target_values=target,
        reference_values=reference,
        target_limits=(-1.0, 1.0),
        reference_limits=(0.0, 3.0),
        bins=20,
    )
    permutation = np.random.default_rng(12).permutation(len(target))
    shuffled = figures.native_joint_histogram(
        target_values=target[permutation],
        reference_values=reference[permutation],
        target_limits=(-1.0, 1.0),
        reference_limits=(0.0, 3.0),
        bins=20,
    )

    assert baseline.finite_n == shuffled.finite_n == 5
    assert baseline.displayed_n == shuffled.displayed_n == 3
    assert baseline.density.shape == (20, 20)
    assert baseline.density.sum() == pytest.approx(3 / 5)
    assert shuffled.density == pytest.approx(baseline.density)


def test_q995_native_display_windows_and_counts_use_declared_axes():
    rows = []
    for pair in ROSTER.pairs:
        for quantile, target_boundary, reference_boundary in (
            (0.995, 3.0 if pair.target_column == "Iba1_norm" else 4.0, 14.0),
            (0.999, 4.0 if pair.target_column == "Iba1_norm" else 5.0, 15.0),
        ):
            rows.append(
                {
                    "donor": "D",
                    "target_column": pair.target_column,
                    "quantile": quantile,
                    "target_tail_boundary_score": target_boundary,
                    "reference_tail_boundary_score": reference_boundary,
                }
            )
    inputs = figures.FigureInputs(
        pair_diagnostics=pd.DataFrame(),
        input_summary=pd.DataFrame(),
        rank_overlap=pd.DataFrame(rows),
        crop_manifest=pd.DataFrame(),
        donors=("D",),
        crop_donor="D",
        projection_root=Path("."),
        crop_root=Path("."),
        source_files=(),
        source_method_versions={},
    )
    reference = np.arange(10.0, 16.0)
    finite = np.ones(6, dtype=bool)
    for pair in ROSTER.pairs:
        target = np.arange(-1.0, 5.0) if pair.target_column == "Iba1_norm" else np.arange(6.0)
        target_limits, reference_limits = figures._native_display_limits(
            inputs,
            donor="D",
            pair=pair,
            target=target,
            reference=reference,
            finite=finite,
        )
        expected_lower = float(np.quantile(target, 0.005)) if target[0] < 0 else 0.0
        expected_upper = 3.0 if target[0] < 0 else 4.0
        assert target_limits == pytest.approx((expected_lower, expected_upper))
        assert reference_limits == pytest.approx((0.0, 14.0))
        result = figures.native_joint_histogram(
            target_values=target,
            reference_values=reference,
            target_limits=target_limits,
            reference_limits=reference_limits,
            bins=20,
        )
        direct = (
            (target >= target_limits[0])
            & (target <= target_limits[1])
            & (reference >= reference_limits[0])
            & (reference <= reference_limits[1])
        )
        assert result.displayed_n == int(direct.sum())
        assert result.density.sum() == pytest.approx(direct.mean())


def test_review_kmeans_clusters_are_reorder_stable_disjoint_and_conserve_counts():
    identifiers = np.array([f"cell-{index:03d}" for index in range(100)])
    target = np.arange(100, dtype=float)
    reference = target[::-1]
    kwargs = {
        "target_key": "target",
        "reference_key": "reference",
        "target_limits": (0.0, 99.0),
        "reference_limits": (0.0, 99.0),
        "bins": 20,
    }
    baseline = figures.review_kmeans_cluster_marginals(
        object_ids=identifiers,
        target_values=target,
        reference_values=reference,
        **kwargs,
    )
    permutation = np.random.default_rng(21).permutation(len(target))
    shuffled = figures.review_kmeans_cluster_marginals(
        object_ids=identifiers[permutation],
        target_values=target[permutation],
        reference_values=reference[permutation],
        **kwargs,
    )

    assert baseline.target_selection_n == baseline.reference_selection_n == 30
    assert baseline.selection_intersection_n == 0
    assert baseline.selected_union_n == 60
    assert baseline.fit_n == baseline.selected_union_n
    assert baseline.fit_red_n + baseline.fit_blue_n == baseline.fit_n
    assert baseline.prediction_n == baseline.finite_n == 100
    assert baseline.predicted_outside_fit_n == 40
    assert baseline.outside_fit_red_n + baseline.outside_fit_blue_n == 40
    assert baseline.red_n + baseline.blue_n == baseline.prediction_n
    assert baseline.cluster_overlap_n == 0
    assert baseline.opposing is True
    assert baseline.reference_high_color == "blue"
    assert baseline.red_percentile_center[0] > baseline.blue_percentile_center[0]
    assert baseline.blue_percentile_center[1] > baseline.red_percentile_center[1]
    assert baseline.red_target_counts.sum() == baseline.red_target_displayed_n
    assert baseline.blue_target_counts.sum() == baseline.blue_target_displayed_n
    assert baseline.red_reference_counts.sum() == baseline.red_reference_displayed_n
    assert baseline.blue_reference_counts.sum() == baseline.blue_reference_displayed_n
    edges = np.linspace(0.0, 99.0, 21)
    assert baseline.red_target_counts.tolist() == np.histogram(
        target[50:], bins=edges
    )[0].tolist()
    assert baseline.blue_target_counts.tolist() == np.histogram(
        target[:50], bins=edges
    )[0].tolist()
    assert (
        baseline.red_reference_counts.tolist()
        == np.histogram(reference[50:], bins=edges)[0].tolist()
    )
    assert (
        baseline.blue_reference_counts.tolist()
        == np.histogram(reference[:50], bins=edges)[0].tolist()
    )

    for field in (
        "red_target_counts",
        "blue_target_counts",
        "red_reference_counts",
        "blue_reference_counts",
    ):
        assert getattr(shuffled, field) == pytest.approx(getattr(baseline, field))
    assert (
        shuffled.red_n,
        shuffled.blue_n,
        shuffled.fit_red_n,
        shuffled.fit_blue_n,
        shuffled.selected_union_n,
        shuffled.opposing,
    ) == (
        baseline.red_n,
        baseline.blue_n,
        baseline.fit_red_n,
        baseline.fit_blue_n,
        baseline.selected_union_n,
        baseline.opposing,
    )


def test_all_finite_prediction_fills_support_absent_from_q70_fit_union():
    identifiers = np.array([f"cell-{index:03d}" for index in range(100)])
    target = np.concatenate((np.arange(200, 230), np.arange(50, 80), np.zeros(40)))
    reference = np.concatenate((np.arange(50, 80), np.arange(200, 230), np.arange(40)))
    target_ranks, _target_percentiles = figures._bright_percentile_ranks(target, identifiers)
    reference_ranks, _reference_percentiles = figures._bright_percentile_ranks(
        reference, identifiers
    )
    tail_n = figures._exact_upper_tail_count(len(target), 0.70)
    fit_union = (target_ranks <= tail_n) | (reference_ranks <= tail_n)
    result = figures.review_kmeans_cluster_marginals(
        object_ids=identifiers,
        target_values=target,
        reference_values=reference,
        target_key="target",
        reference_key="reference",
        target_limits=(0.0, 500.0),
        reference_limits=(0.0, 500.0),
        bins=20,
    )

    assert fit_union[:60].all() and not fit_union[60:].any()
    assert result.fit_red_n == result.fit_blue_n == 30
    assert not np.any(reference[fit_union] < 25.0)
    assert result.blue_reference_counts[0] > 0
    assert result.outside_fit_blue_n > 0
    assert result.fit_n == int(fit_union.sum())
    assert result.prediction_n == result.finite_n == len(target)
    assert result.red_n + result.blue_n == result.prediction_n
    assert result.predicted_outside_fit_n == int((~fit_union).sum())


def test_rank_cluster_histograms_are_exactly_transpose_invariant():
    identifiers = np.array([f"cell-{index:03d}" for index in range(100)])
    target = np.arange(100, dtype=float)
    reference = target[::-1]
    baseline = figures.review_kmeans_cluster_marginals(
        object_ids=identifiers,
        target_values=target,
        reference_values=reference,
        target_key="target",
        reference_key="reference",
        target_limits=(0.0, 99.0),
        reference_limits=(0.0, 99.0),
        bins=20,
    )
    transposed = figures.review_kmeans_cluster_marginals(
        object_ids=identifiers,
        target_values=reference,
        reference_values=target,
        target_key="reference",
        reference_key="target",
        target_limits=(0.0, 99.0),
        reference_limits=(0.0, 99.0),
        bins=20,
    )

    assert transposed.red_n == baseline.blue_n
    assert transposed.blue_n == baseline.red_n
    assert transposed.red_target_counts.tolist() == baseline.blue_reference_counts.tolist()
    assert transposed.blue_target_counts.tolist() == baseline.red_reference_counts.tolist()
    assert transposed.red_reference_counts.tolist() == baseline.blue_target_counts.tolist()
    assert transposed.blue_reference_counts.tolist() == baseline.red_target_counts.tolist()
    assert transposed.red_percentile_center == pytest.approx(baseline.blue_percentile_center[::-1])
    assert transposed.blue_percentile_center == pytest.approx(baseline.red_percentile_center[::-1])


def test_bright_percentile_ranks_use_stable_object_id_ties_and_exact_tail_count():
    ranks, percentiles = figures._bright_percentile_ranks(
        np.array([2.0, 2.0, 1.0]),
        np.array(["b", "a", "c"]),
    )

    assert ranks.tolist() == [2, 1, 3]
    assert percentiles.tolist() == [0.5, 1.0, 0.0]
    assert figures._exact_upper_tail_count(10, 0.7) == 3
    assert figures._exact_upper_tail_count(11, 0.7) == 4


def test_review_kmeans_flags_when_one_cluster_wins_both_markers():
    identifiers = np.array([f"cell-{index:03d}" for index in range(100)])
    target = np.arange(100, dtype=float)
    reference = target.copy()

    result = figures.review_kmeans_cluster_marginals(
        object_ids=identifiers,
        target_values=target,
        reference_values=reference,
        target_key="target",
        reference_key="reference",
        target_limits=(0.0, 99.0),
        reference_limits=(0.0, 99.0),
        bins=20,
    )

    assert result.opposing is False
    assert result.reference_high_color == "red"
    assert result.red_percentile_center[0] > result.blue_percentile_center[0]
    assert result.red_percentile_center[1] > result.blue_percentile_center[1]
    assert result.fit_red_n + result.fit_blue_n == result.selected_union_n
    assert result.red_n + result.blue_n == result.finite_n
    assert result.cluster_overlap_n == 0


def test_marginal_lines_use_one_shared_peak_normalization():
    red, blue = figures._shared_peak_marginal_lines(
        np.array([0, 4, 2]),
        4,
        np.array([0, 2, 4]),
        8,
    )

    assert red == pytest.approx([0, 1, 0.5])
    assert blue == pytest.approx([0, 0.25, 0.5])


def test_joint_box_and_marginal_bboxes_are_physically_rotation_invariant():
    if importlib.util.find_spec("matplotlib") is None:
        pytest.skip("matplotlib notebook extra is not installed")
    plt, _line2d = figures._plot_modules()
    dpi = 100
    fig = plt.figure(figsize=figures._joint_figure_size(2, 3), dpi=dpi)
    outer = figures._joint_outer_grid(fig, pair_n=2, donor_n=3)
    target_ax, joint_ax, reference_ax = figures._append_joint_axes(fig, outer[0, 0])
    fig.canvas.draw()

    target_box = target_ax.get_window_extent()
    joint_box = joint_ax.get_window_extent()
    reference_box = reference_ax.get_window_extent()
    assert joint_box.width == pytest.approx(joint_box.height, abs=1e-8)
    assert target_box.width == pytest.approx(figures.MARGINAL_DEPTH_INCHES * dpi)
    assert reference_box.height == pytest.approx(figures.MARGINAL_DEPTH_INCHES * dpi)
    assert joint_box.x0 - target_box.x1 == pytest.approx(figures.MARGINAL_PAD_INCHES * dpi)
    assert joint_box.y0 - reference_box.y1 == pytest.approx(figures.MARGINAL_PAD_INCHES * dpi)
    plt.close(fig)


def test_joint_renderer_has_requested_legend_and_no_bottom_explanatory_text():
    source = inspect.getsource(figures.render_joint_distributions)
    axes_source = inspect.getsource(figures._append_joint_axes)

    assert ".axvline(" not in source
    assert ".axhline(" not in source
    assert "descriptive q99" not in source
    assert 'label="red: y-axis marker"' in source
    assert 'label="blue: x-axis marker"' in source
    assert "fig.text(" not in source
    assert axes_source.count("size=MARGINAL_DEPTH_INCHES") == 2
    assert axes_source.count("pad=MARGINAL_PAD_INCHES") == 2
    assert source.count("_draw_rotated_marginal(") == 2
    assert "joint_ax.set_xscale(" not in source
    assert "joint_ax.set_yscale(" not in source


def test_nondefault_cluster_selection_cannot_replace_pinned_workspace(capsys):
    with pytest.raises(SystemExit) as error:
        figures.main(["--cluster-selection-quantile", "0.60"])

    assert error.value.code == 2
    assert "cannot replace the pinned q70 workspace" in capsys.readouterr().err


def test_cli_renders_owned_png_and_pdf_figures(tmp_path):
    if importlib.util.find_spec("matplotlib") is None:
        pytest.skip("matplotlib notebook extra is not installed")
    roots = _write_review_artifacts(tmp_path)
    development_root = tmp_path / "development"
    argv = [
        "--projection-root",
        str(roots[0]),
        "--survey-root",
        str(roots[1]),
        "--overlap-root",
        str(roots[2]),
        "--crop-root",
        str(roots[3]),
        "--development-root",
        str(development_root),
        "--label",
        "figures",
        "--dpi",
        "72",
    ]

    assert figures.main(argv) == 0
    output = development_root / "figures"
    stale = output / "stale.txt"
    stale.write_text("old", encoding="utf-8")
    sibling = development_root / "keep"
    sibling.mkdir()
    (sibling / "user.txt").write_text("keep", encoding="utf-8")
    assert figures.main(argv) == 0

    assert not stale.exists()
    assert not (output / "cell_rank_landscape.png").exists()
    assert (sibling / "user.txt").read_text(encoding="utf-8") == "keep"
    for name in figures.FIGURE_NAMES:
        path = output / name
        expected_magic = b"\x89PNG\r\n\x1a\n" if path.suffix == ".png" else b"%PDF-"
        assert path.read_bytes().startswith(expected_magic)
        assert path.stat().st_size > 5_000
    manifest = json.loads((output / "figure_manifest.json").read_text(encoding="utf-8"))
    assert manifest["figures"] == list(figures.FIGURE_NAMES)
    assert manifest["source_method_versions"]["projection"] == [PROJECTION_METHOD_VERSION]
    assert manifest["source_method_versions"]["boundary_crops"] == [figures.BOUNDARY_METHOD_VERSION]
    assert manifest["restore_outputs_read"] is False
    assert manifest["creates_cell_calls"] is False
    joint = manifest["fundamental_pair_joint_distributions"]
    assert joint["retained_cell_level_arrays"] is False
    assert joint["outside_window_handling"] == "filtered, never clipped into edge bins"
    assert joint["axes"] == "donor-specific native-linear QuPath values"
    assert joint["legend"] == ["red: y-axis marker", "blue: x-axis marker"]
    assert joint["bottom_explanatory_text"].startswith("none")
    assert joint["marginal_geometry"] == {
        "depth_inches": figures.MARGINAL_DEPTH_INCHES,
        "pad_inches": figures.MARGINAL_PAD_INCHES,
        "density_headroom": figures.MARGINAL_DENSITY_HEADROOM,
        "line_width_points": figures.MARGINAL_LINE_WIDTH,
        "rotation_contract": (
            "bottom reference marginal is the exact 90-degree formatting counterpart "
            "of the left target marginal"
        ),
        "joint_box": "physically square after equal fixed left/bottom margins",
    }
    clusters = joint["cluster_populations"]
    assert clusters["method_status"] == "review_only_unverified"
    assert clusters["selection_quantile"] == figures.CLUSTER_SELECTION_QUANTILE == 0.70
    assert clusters["default_selection_quantile"] == 0.70
    assert clusters["selection_configuration"] == "manuscript_anchored_default"
    assert "not optimized" in clusters["selection_rationale"]
    assert "depend" in clusters["selection_dependence"]
    assert clusters["selection_tie_break"] == "object_id_ascending"
    assert clusters["feature_range"] == [0.0, 1.0]
    assert clusters["brightest_feature_value"] == 1.0
    assert "bright-percentile ranks" in clusters["preprocessing"]
    assert clusters["fit_population"] == "every cell in the q70 union"
    assert clusters["fit_row_order"] == "object_id ascending"
    assert clusters["prediction_population"] == "every finite donor/pair cell"
    assert "included" in clusters["outside_fit_policy"]
    assert clusters["third_population"].startswith("none")
    assert clusters["overlap_policy"].startswith("impossible")
    assert "not a verified population" in clusters["interpretation"]
    panels = clusters["panels"]
    assert len(panels) == len(ROSTER.pairs) * 2
    assert {
        (panel["donor"], panel["target_column"], panel["reference_column"]) for panel in panels
    } == {
        (donor, pair.target_column, pair.reference_column)
        for donor in ("D1", "D2")
        for pair in ROSTER.pairs
    }
    for panel in panels:
        assert panel["state"] in {"OPPOSING", "NON_OPPOSING"}
        assert panel["reference_high_color"] in {"red", "blue"}
        assert panel["cluster_overlap_n"] == 0
        assert panel["selection_quantile"] == 0.70
        assert panel["fit_n"] == panel["selected_union_n"]
        assert panel["fit_red_n"] + panel["fit_blue_n"] == panel["fit_n"]
        assert panel["prediction_n"] == panel["finite_n"]
        assert panel["red_n"] + panel["blue_n"] == panel["prediction_n"]
        assert (
            panel["selected_union_n"] + panel["predicted_outside_fit_n"]
            == panel["finite_n"]
        )
        assert (
            panel["outside_fit_red_n"] + panel["outside_fit_blue_n"]
            == panel["predicted_outside_fit_n"]
        )
        assert panel["target_selection_fraction"] == pytest.approx(
            panel["target_selection_n"] / panel["finite_n"]
        )
        assert panel["reference_selection_fraction"] == pytest.approx(
            panel["reference_selection_n"] / panel["finite_n"]
        )
        assert panel["selection_intersection_fraction"] == pytest.approx(
            panel["selection_intersection_n"] / panel["finite_n"]
        )
        assert panel["selected_union_fraction"] == pytest.approx(
            panel["selected_union_n"] / panel["finite_n"]
        )
        assert panel["predicted_outside_fit_fraction"] == pytest.approx(
            panel["predicted_outside_fit_n"] / panel["finite_n"]
        )
        assert "excluded_double_low_n" not in panel
        assert panel["red_target_displayed_n"] <= panel["red_n"]
        assert panel["blue_target_displayed_n"] <= panel["blue_n"]
        assert panel["red_reference_displayed_n"] <= panel["red_n"]
        assert panel["blue_reference_displayed_n"] <= panel["blue_n"]
        assert "red_scaled_center" not in panel
        assert "target_q995_scale" not in panel
        frame = pd.read_parquet(roots[0] / f"donor_id={panel['donor']}" / "data.parquet")
        target = frame[panel["target_column"]].to_numpy(float)
        reference = frame[panel["reference_column"]].to_numpy(float)
        finite = np.isfinite(target) & np.isfinite(reference)
        target_limits = panel["target_display_limits"]
        reference_limits = panel["reference_display_limits"]
        direct = (
            finite
            & (target >= target_limits[0])
            & (target <= target_limits[1])
            & (reference >= reference_limits[0])
            & (reference <= reference_limits[1])
        )
        assert panel["displayed_n"] == int(direct.sum())
    assert joint["q99_guides"] == "none"
    assert "cell_rank_landscape" not in manifest
    assert "CD79a" not in {pair["target_column"] for pair in manifest["pairs"]}
    assert "CD79a" not in {panel["target_column"] for panel in panels}
    assert manifest["crop_montage"]["pdf_pages"] == len(ROSTER.pairs) == 2
