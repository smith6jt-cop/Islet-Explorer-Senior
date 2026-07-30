from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import build_bundle as build  # noqa: E402
import score  # noqa: E402


def _population() -> pd.DataFrame:
    rows = []
    for index in range(10):
        if index < 4:
            label, status, region = "Immune", "anchor", "tissue"
        elif index < 7:
            label, status, region = "Epithelial", "inferred", "core"
        else:
            label, status, region = "Other", "other", "tissue"
        rows.append(
            {
                "donor_id": "D1",
                "object_id": f"id-{index}",
                "image": "image",
                "X_centroid": float(index),
                "Y_centroid": float(index),
                "cell_region": region,
                "islet_num": "",
                "cell_area": 50.0,
                "cell_solidity": 0.9,
                "pipeline_broad_type": label,
                "pipeline_broad_assignment_status": status,
                "pipeline_best_broad_type": label,
                "pipeline_confidence": 0.9,
                "pipeline_ranked_broad_probabilities": "[]",
                "pipeline_broad_reason": "test",
                "pipeline_qc_analysis_eligible": True,
            }
        )
    return pd.DataFrame(rows)


def test_sampling_is_disjoint_deterministic_and_weighted():
    population = _population()
    first = build.sample_population(
        population, donor_id="D1", n=10, rng=np.random.default_rng(17)
    )
    second = build.sample_population(
        population, donor_id="D1", n=10, rng=np.random.default_rng(17)
    )
    pd.testing.assert_frame_equal(first, second)
    assert first["object_id"].is_unique
    assert first["stratum"].value_counts().to_dict() == {
        "immune_called": 4,
        "islet_proximal": 3,
        "uniform": 3,
    }
    reconstructed = first.groupby("stratum")["weight"].sum().sort_index()
    population_sizes = first.groupby("stratum")["stratum_N"].first().sort_index()
    np.testing.assert_allclose(reconstructed, population_sizes)

    blinded = first.copy()
    blinded.insert(0, "review_id", [f"GS{i:05d}" for i in range(len(first))])
    sheet = build.blinded_review_index(blinded)
    assert not any(
        "pipeline" in column or column in {"stratum", "weight"}
        for column in sheet.columns
    )
    assert sheet.columns.tolist() == [
        "review_id",
        "donor_id",
        "image",
        "X_centroid",
        "Y_centroid",
    ]
    with pytest.raises(ValueError, match="too small"):
        build.sample_population(
            population, donor_id="D1", n=2, rng=np.random.default_rng(17)
        )
    assert build._normalise_annotators(["Reviewer-A", "Reviewer_B"]) == (
        "Reviewer-A",
        "Reviewer_B",
    )
    with pytest.raises(ValueError, match="may contain only"):
        build._normalise_annotators(["Reviewer A", "Reviewer-B"])


def test_composite_key_mismatch_and_duplicate_fail_closed():
    population = _population()
    cells = population[
        ["donor_id", *build.CELL_COLUMNS]
    ].copy()
    assignments = population[
        ["donor_id", "object_id", *build.PREDICTION_RENAMES.values()]
    ].rename(columns={value: key for key, value in build.PREDICTION_RENAMES.items()})
    merged = build.merge_population(cells, assignments, donor_id="D1")
    assert len(merged) == len(population)

    duplicate = pd.concat([assignments, assignments.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        build.merge_population(cells, duplicate, donor_id="D1")

    with pytest.raises(ValueError, match="UUID universes differ"):
        build.merge_population(cells.iloc[:-1], assignments, donor_id="D1")


def test_frozen_abstentions_and_zero_f1_are_retained():
    key = pd.DataFrame(
        {
            "pipeline_broad_type": [
                "Immune",
                "Other",
                "Ambiguous",
                "Unavailable",
            ],
            "pipeline_broad_assignment_status": [
                "anchor",
                "other",
                "ambiguous",
                "unavailable",
            ],
        }
    )
    predicted = score.frozen_predictions(key, {"Immune", "Other"})
    assert predicted.tolist() == [
        "Immune",
        "Other",
        "Ambiguous",
        "Unavailable",
    ]

    metrics = score.rates(
        ["Immune", "Other"],
        ["Other", "Other"],
        [1.0, 1.0],
        ["Immune", "Other"],
    ).set_index("broad_type")
    assert metrics.loc["Immune", "recall"] == 0
    assert metrics.loc["Immune", "f1"] == 0

    bad = key.iloc[[0]].copy()
    bad.loc[:, "pipeline_broad_type"] = "Immune"
    bad.loc[:, "pipeline_broad_assignment_status"] = "other"
    with pytest.raises(ValueError, match="inconsistent"):
        score.frozen_predictions(bad, {"Immune", "Other"})


def test_consensus_reports_unsure_and_disagreement():
    annotations = pd.DataFrame(
        [
            ("A", "GS1", "id1", "Immune"),
            ("B", "GS1", "id1", "Immune"),
            ("A", "GS2", "id2", "Unsure"),
            ("B", "GS2", "id2", "Other"),
            ("A", "GS3", "id3", "Epithelial"),
            ("B", "GS3", "id3", "Other"),
        ],
        columns=["annotator", "review_id", "object_id", "compartment"],
    )
    agreed, diagnostics = score.consensus(annotations, ("A", "B"))
    assert agreed[["review_id", "truth"]].to_records(index=False).tolist() == [
        ("GS1", "Immune")
    ]
    assert diagnostics["n_unsure_either"] == 1
    assert diagnostics["n_disagree"] == 1


def test_score_reads_only_frozen_snapshot_and_writes_once(tmp_path):
    bundle = tmp_path / "bundle"
    completed = bundle / "completed"
    completed.mkdir(parents=True)
    key = pd.DataFrame(
        {
            "review_id": ["GS00000", "GS00001"],
            "donor_id": ["D1", "D1"],
            "object_id": ["id-0", "id-1"],
            "stratum": ["immune_called", "uniform"],
            "weight": [1.0, 1.0],
            "pipeline_broad_type": ["Immune", "Ambiguous"],
            "pipeline_broad_assignment_status": ["anchor", "ambiguous"],
        }
    )
    key.to_csv(bundle / "sample_key.csv", index=False)
    annotations = {
        "A": pd.DataFrame(
            {
                "annotator": ["A", "A"],
                "object_id": ["id-0", "id-1"],
                "gs_review_id": ["GS00000", "GS00001"],
                "compartment": ["Immune", "Epithelial"],
            }
        ),
        "B": pd.DataFrame(
            {
                "annotator": ["B", "B"],
                "object_id": ["id-0", "id-1"],
                "gs_review_id": ["GS00000", "GS00001"],
                "compartment": ["Immune", ""],
            }
        ),
    }
    for annotator, frame in annotations.items():
        frame.to_csv(completed / f"annotations_{annotator}_image.csv", index=False)

    manifest = {
        "schema_version": 2,
        "protocol_version": "gold-standard-broad-audit-v2",
        "annotators": ["A", "B"],
        "vocabulary": list(build.ANNOTATION_VOCABULARY),
        "files": [build._file_record(bundle / "sample_key.csv", root=bundle)],
    }
    manifest["content_id"] = score.sha256_json(manifest)
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match="must annotate every"):
        score.main(["--bundle", str(bundle)])
    assert not (bundle / "scores").exists()

    annotations["B"].loc[1, "compartment"] = "Epithelial"
    annotations["B"].to_csv(
        completed / "annotations_B_image.csv", index=False
    )
    assert score.main(["--bundle", str(bundle)]) == 0
    summary = json.loads((bundle / "scores" / "coverage_accuracy.json").read_text())
    assert summary["weighted_coverage"] == 0.5
    assert summary["n_ambiguous"] == 1
    score_manifest = json.loads(
        (bundle / "scores" / "score_manifest.json").read_text()
    )
    content_id = score_manifest.pop("content_id")
    assert score.sha256_json(score_manifest) == content_id
    assert len(score_manifest["annotation_files"]) == 2
    assert {item["path"] for item in score_manifest["output_files"]} == {
        "scores/broad_metrics_weighted.csv",
        "scores/confusion_raw.csv",
        "scores/confusion_weighted.csv",
        "scores/coverage_accuracy.json",
        "scores/inter_annotator.json",
        "scores/per_donor_coverage_accuracy.csv",
        "scores/scored_consensus.csv",
    }
    with pytest.raises(FileExistsError, match="will not be overwritten"):
        score.main(["--bundle", str(bundle)])
