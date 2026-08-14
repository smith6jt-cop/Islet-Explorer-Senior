#!/usr/bin/env python3
"""Score a frozen, blinded broad-cell-type audit bundle.

The scorer reads only the predictions captured in ``sample_key.csv`` when the
bundle was created. It never reruns typing and cannot change calibration,
typing rules, configuration, or release status.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
# `core/phenocycler/` since the subtree merge, so this runs without an editable install. Index 0
# pre-dates the move and is kept: the checkout wins over any installed copy.
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.artifacts import sha256_json  # noqa: E402


ABSTENTIONS = ("Ambiguous", "Unavailable")
ACCEPTED_STATUSES = frozenset({"anchor", "inferred", "other"})
STATUS_LABELS = {
    "ambiguous": "Ambiguous",
    "unavailable": "Unavailable",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path, *, root: Path) -> dict:
    resolved = path.resolve()
    root = root.resolve()
    name = (
        str(resolved.relative_to(root))
        if resolved.is_relative_to(root)
        else str(resolved)
    )
    return {
        "path": name,
        "size": resolved.stat().st_size,
        "sha256": _sha256(resolved),
    }


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def _write_json(path: Path, value) -> None:
    path.write_text(
        json.dumps(_json_safe(value), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def load_bundle(bundle: Path) -> tuple[dict, pd.DataFrame]:
    """Validate the immutable bundle contract and return its private key."""

    bundle = bundle.expanduser().resolve()
    manifest_path = bundle / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"audit manifest is missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    content_id = manifest.pop("content_id", None)
    if not isinstance(content_id, str) or sha256_json(manifest) != content_id:
        raise ValueError("audit manifest content_id is invalid")
    manifest["content_id"] = content_id
    if manifest.get("schema_version") != 2:
        raise ValueError(
            f"unsupported audit schema version {manifest.get('schema_version')!r}"
        )
    if manifest.get("protocol_version") != "gold-standard-broad-audit-v2":
        raise ValueError("bundle is not a broad-audit-v2 snapshot")

    for record in manifest.get("files", []):
        relative = Path(str(record["path"]))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"invalid bundle-relative path: {relative}")
        path = (bundle / relative).resolve()
        if not path.is_relative_to(bundle):
            raise ValueError(f"bundle file escapes its root: {relative}")
        if not path.is_file():
            raise FileNotFoundError(f"contracted bundle file is missing: {relative}")
        if (
            path.stat().st_size != int(record["size"])
            or _sha256(path) != record["sha256"]
        ):
            raise ValueError(f"contracted bundle file changed: {relative}")

    key_path = bundle / "sample_key.csv"
    key = pd.read_csv(
        key_path,
        dtype={"review_id": str, "donor_id": str, "object_id": str},
    )
    required = {
        "review_id",
        "donor_id",
        "object_id",
        "stratum",
        "weight",
        "pipeline_broad_type",
        "pipeline_broad_assignment_status",
    }
    missing = required.difference(key.columns)
    if missing:
        raise KeyError(f"sample key missing columns: {sorted(missing)}")
    if key.duplicated(["review_id"]).any():
        raise ValueError("sample key contains duplicate review IDs")
    if key.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("sample key contains duplicate donor/object keys")
    weights = pd.to_numeric(key["weight"], errors="coerce")
    if not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError("sample weights must be finite and positive")
    key["weight"] = weights.astype(float)
    return manifest, key


def frozen_predictions(key: pd.DataFrame, broad_labels: set[str]) -> pd.Series:
    """Translate the frozen label/status contract into scored categories."""

    predictions: list[str] = []
    for label, status in zip(
        key["pipeline_broad_type"].astype(str),
        key["pipeline_broad_assignment_status"].astype(str),
    ):
        if status in ACCEPTED_STATUSES:
            if label not in broad_labels:
                raise ValueError(
                    f"accepted prediction {label!r} is outside the frozen vocabulary"
                )
            if (status == "other") != (label == "Other"):
                raise ValueError(
                    f"inconsistent frozen prediction label/status: {label!r}/{status!r}"
                )
            predictions.append(label)
            continue
        expected = STATUS_LABELS.get(status)
        if expected is None or label != expected:
            raise ValueError(
                f"inconsistent frozen prediction label/status: {label!r}/{status!r}"
            )
        predictions.append(expected)
    return pd.Series(predictions, index=key.index, name="predicted")


def load_annotations(
    completed: Path,
    *,
    key: pd.DataFrame,
    annotators: tuple[str, str],
    vocabulary: set[str],
) -> pd.DataFrame:
    """Load exactly one annotation row per sampled cell and annotator."""

    files = sorted(completed.glob("annotations_*.csv"))
    if not files:
        raise FileNotFoundError(f"no annotations_*.csv under {completed}")
    frames = [
        pd.read_csv(
            path,
            dtype={
                "annotator": str,
                "object_id": str,
                "gs_review_id": str,
                "compartment": str,
            },
            keep_default_na=False,
        )
        for path in files
    ]
    annotations = pd.concat(frames, ignore_index=True)
    required = {"annotator", "object_id", "gs_review_id", "compartment"}
    missing = required.difference(annotations.columns)
    if missing:
        raise KeyError(f"annotation files missing columns: {sorted(missing)}")
    annotations = annotations.rename(columns={"gs_review_id": "review_id"})
    annotations["annotator"] = annotations["annotator"].astype(str)
    annotations["review_id"] = annotations["review_id"].astype(str)
    annotations["object_id"] = annotations["object_id"].astype(str)
    annotations["compartment"] = annotations["compartment"].astype(str).str.strip()

    unexpected_annotators = sorted(set(annotations["annotator"]) - set(annotators))
    if unexpected_annotators:
        raise ValueError(f"undeclared annotators: {unexpected_annotators}")
    unexpected_labels = sorted(set(annotations["compartment"]) - vocabulary - {""})
    if unexpected_labels:
        raise ValueError(f"labels outside the frozen vocabulary: {unexpected_labels}")
    if annotations.duplicated(["annotator", "review_id", "object_id"]).any():
        raise ValueError("annotation files contain duplicate annotator/review/object rows")

    expected_pairs = set(map(tuple, key[["review_id", "object_id"]].to_numpy()))
    for annotator in annotators:
        observed = annotations.loc[
            annotations["annotator"].eq(annotator), ["review_id", "object_id"]
        ]
        observed_pairs = set(map(tuple, observed.to_numpy()))
        missing_pairs = expected_pairs - observed_pairs
        extra_pairs = observed_pairs - expected_pairs
        if missing_pairs or extra_pairs or len(observed) != len(expected_pairs):
            raise ValueError(
                f"annotator {annotator!r} does not match the frozen sample: "
                f"missing={len(missing_pairs)}, extra={len(extra_pairs)}"
            )

    mapping = key[["review_id", "object_id"]]
    checked = annotations.merge(
        mapping,
        on=["review_id", "object_id"],
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    if not checked["_merge"].eq("both").all():
        raise ValueError("annotation review IDs do not map to their frozen UUIDs")
    return annotations


def _kappa(a: pd.Series, b: pd.Series) -> float:
    categories = sorted(set(a) | set(b))
    if not categories:
        return float("nan")
    matrix = pd.crosstab(
        pd.Categorical(a, categories=categories),
        pd.Categorical(b, categories=categories),
        dropna=False,
    ).to_numpy(float)
    n = float(matrix.sum())
    if n == 0:
        return float("nan")
    observed = float(np.trace(matrix) / n)
    expected = float(
        ((matrix.sum(axis=0) / n) * (matrix.sum(axis=1) / n)).sum()
    )
    return (observed - expected) / (1.0 - expected) if expected < 1 else float("nan")


def consensus(
    annotations: pd.DataFrame,
    annotators: tuple[str, str],
) -> tuple[pd.DataFrame, dict]:
    """Retain only agreed, non-unsure labels from the two blinded annotators."""

    wide = annotations.pivot(
        index=["review_id", "object_id"],
        columns="annotator",
        values="compartment",
    ).reset_index()
    first = wide[annotators[0]].fillna("")
    second = wide[annotators[1]].fillna("")
    both = first.ne("") & second.ne("")
    unsure = both & (first.eq("Unsure") | second.eq("Unsure"))
    scorable = both & ~unsure
    agree = scorable & first.eq(second)
    diagnostics = {
        "n_rows": int(len(wide)),
        "n_both_annotated": int(both.sum()),
        "n_unannotated_either": int((~both).sum()),
        "n_unsure_either": int(unsure.sum()),
        "n_scorable": int(scorable.sum()),
        "n_agree": int(agree.sum()),
        "n_disagree": int((scorable & ~agree).sum()),
        "raw_agreement": (
            float(agree.sum() / scorable.sum()) if scorable.any() else float("nan")
        ),
        "cohens_kappa": _kappa(first[scorable], second[scorable]),
        "annotators": list(annotators),
    }
    agreed = wide.loc[agree, ["review_id", "object_id"]].copy()
    agreed["truth"] = first[agree].to_numpy()
    return agreed, diagnostics


def rates(
    truth,
    predicted,
    weight,
    labels,
) -> pd.DataFrame:
    """Return weighted one-vs-rest metrics without hiding abstentions."""

    truth = np.asarray(truth)
    predicted = np.asarray(predicted)
    weight = np.asarray(weight, dtype=float)
    rows = []
    for label in labels:
        truth_mask = truth == label
        predicted_mask = predicted == label
        tp = float(weight[truth_mask & predicted_mask].sum())
        fp = float(weight[~truth_mask & predicted_mask].sum())
        fn = float(weight[truth_mask & ~predicted_mask].sum())
        precision = tp / (tp + fp) if tp + fp else np.nan
        recall = tp / (tp + fn) if tp + fn else np.nan
        if tp == 0 and (tp + fp + fn) > 0:
            f1 = 0.0
        elif np.isnan(precision) or np.isnan(recall):
            f1 = np.nan
        else:
            f1 = 2 * precision * recall / (precision + recall)
        rows.append(
            {
                "broad_type": label,
                "n_truth": int(truth_mask.sum()),
                "n_predicted": int(predicted_mask.sum()),
                "truth_weight": float(weight[truth_mask].sum()),
                "predicted_weight": float(weight[predicted_mask].sum()),
                "true_positive_weight": tp,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            }
        )
    return pd.DataFrame(rows)


def prediction_coverage(frame: pd.DataFrame) -> dict:
    """Weighted call coverage over the full frozen audit sample."""

    weight = frame["weight"].to_numpy(float)
    called = ~frame["predicted"].isin(ABSTENTIONS).to_numpy()
    total = float(weight.sum())
    called_weight = float(weight[called].sum())
    return {
        "n_frozen_sample": int(len(frame)),
        "weighted_sample_mass": total,
        "weighted_coverage": called_weight / total if total else np.nan,
        "n_ambiguous": int(frame["predicted"].eq("Ambiguous").sum()),
        "n_unavailable": int(frame["predicted"].eq("Unavailable").sum()),
    }


def accuracy_summary(scored: pd.DataFrame) -> dict:
    """Weighted accuracy over the agreed, non-unsure annotation subset."""

    weight = scored["weight"].to_numpy(float)
    called = ~scored["predicted"].isin(ABSTENTIONS).to_numpy()
    correct = scored["truth"].eq(scored["predicted"]).to_numpy()
    total = float(weight.sum())
    called_weight = float(weight[called].sum())
    return {
        "n_consensus": int(len(scored)),
        "weighted_consensus_mass": total,
        "weighted_consensus_accuracy_including_abstentions": (
            float(weight[correct].sum() / total) if total else np.nan
        ),
        "weighted_consensus_accuracy_when_called": (
            float(weight[correct & called].sum() / called_weight)
            if called_weight
            else np.nan
        ),
    }


def donor_summary(key: pd.DataFrame, scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor, full_frame in key.groupby("donor_id", sort=True):
        consensus_frame = scored.loc[scored["donor_id"].eq(donor)]
        rows.append(
            {
                "donor_id": donor,
                **prediction_coverage(full_frame),
                **accuracy_summary(consensus_frame),
            }
        )
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bundle",
        type=Path,
        required=True,
        help="exact immutable bundle directory containing manifest.json",
    )
    args = parser.parse_args(argv)
    bundle = args.bundle.expanduser().resolve()

    manifest, key = load_bundle(bundle)
    annotators_raw = tuple(map(str, manifest.get("annotators", ())))
    if len(annotators_raw) != 2 or len(set(annotators_raw)) != 2:
        raise ValueError("bundle must declare exactly two distinct annotators")
    annotators = (annotators_raw[0], annotators_raw[1])
    vocabulary = set(map(str, manifest.get("vocabulary", ())))
    if "Unsure" not in vocabulary:
        raise ValueError("bundle vocabulary must include Unsure")
    broad_labels = vocabulary - {"Unsure"}

    key = key.copy()
    key["predicted"] = frozen_predictions(key, broad_labels)
    annotations = load_annotations(
        bundle / "completed",
        key=key,
        annotators=annotators,
        vocabulary=vocabulary,
    )
    truth, agreement = consensus(annotations, annotators)
    if agreement["n_unannotated_either"]:
        raise ValueError(
            "both reviewers must annotate every frozen sample row before scoring; "
            f"{agreement['n_unannotated_either']} rows are incomplete"
        )
    if truth.empty:
        raise ValueError("no agreed, non-unsure annotations are available to score")
    scored = truth.merge(
        key[
            [
                "review_id",
                "object_id",
                "donor_id",
                "stratum",
                "weight",
                "predicted",
            ]
        ],
        on=["review_id", "object_id"],
        how="left",
        validate="one_to_one",
    )
    if scored["predicted"].isna().any():
        raise ValueError("a consensus cell has no frozen pipeline prediction")

    coverage = prediction_coverage(key)
    accuracy = accuracy_summary(scored)
    overall = {**coverage, **accuracy}
    per_donor = donor_summary(key, scored)
    weighted = rates(
        scored["truth"],
        scored["predicted"],
        scored["weight"],
        sorted(broad_labels),
    )
    confusion_raw = pd.crosstab(
        scored["truth"],
        scored["predicted"],
        rownames=["human_truth"],
        colnames=["pipeline_prediction"],
        dropna=False,
    )
    confusion_weighted = scored.pivot_table(
        index="truth",
        columns="predicted",
        values="weight",
        aggfunc="sum",
        fill_value=0.0,
        observed=False,
    )
    confusion_weighted.index.name = "human_truth"
    confusion_weighted.columns.name = "pipeline_prediction"

    print("=" * 72)
    print("1. BLINDED INTER-ANNOTATOR AGREEMENT")
    print("=" * 72)
    for name, value in agreement.items():
        print(f"  {name:36s} {value}")
    print("\n" + "=" * 72)
    print("2. FULL-SNAPSHOT COVERAGE AND CONSENSUS-CONDITIONAL ACCURACY")
    print("=" * 72)
    for name, value in overall.items():
        print(f"  {name:36s} {value}")
    print(
        "\nDesign-weighted estimates target only the selected donor population; "
        "reviewer exclusions may limit that interpretation."
    )
    print("\n" + weighted.round(4).to_string(index=False))
    print("\n" + "=" * 72)
    print("3. WEIGHTED CONFUSION (Ambiguous/Unavailable remain explicit)")
    print("=" * 72)
    print(confusion_weighted.round(3).to_string())

    out = bundle / "scores"
    if out.exists():
        raise FileExistsError(
            f"score output already exists and will not be overwritten: {out}"
        )
    out.mkdir()
    _write_json(out / "inter_annotator.json", agreement)
    _write_json(out / "coverage_accuracy.json", overall)
    weighted.to_csv(out / "broad_metrics_weighted.csv", index=False)
    per_donor.to_csv(out / "per_donor_coverage_accuracy.csv", index=False)
    confusion_raw.to_csv(out / "confusion_raw.csv")
    confusion_weighted.to_csv(out / "confusion_weighted.csv")
    scored.to_csv(out / "scored_consensus.csv", index=False)

    annotation_files = sorted((bundle / "completed").glob("annotations_*.csv"))
    output_files = sorted(
        path
        for path in out.iterdir()
        if path.is_file() and path.name != "score_manifest.json"
    )
    score_manifest = {
        "schema_version": 1,
        "protocol_version": "gold-standard-broad-score-v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "bundle_content_id": manifest["content_id"],
        "source_run_id": manifest.get("source", {}).get("run_id"),
        "annotation_files": [
            _file_record(path, root=bundle) for path in annotation_files
        ],
        "inter_annotator": _json_safe(agreement),
        "prediction_coverage": _json_safe(coverage),
        "consensus_accuracy": _json_safe(accuracy),
        "output_files": [_file_record(path, root=bundle) for path in output_files],
    }
    score_manifest["content_id"] = sha256_json(score_manifest)
    _write_json(out / "score_manifest.json", score_manifest)
    print(f"\n[done] immutable score outputs: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
