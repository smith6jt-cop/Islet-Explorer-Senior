#!/usr/bin/env python3
"""Build a blinded, audit-only broad-cell-type annotation bundle.

The bundle snapshots predictions from one validated content-addressed
Phenocycler run. Human annotations are never read here and never become an
input to calibration, typing, configuration, or release decisions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
# `core/phenocycler/` since the subtree merge; an installed copy takes precedence either way.
sys.path.insert(0, str(REPO / "core" / "phenocycler"))

from phenocycler.artifacts import RunManifest, StageManifest, sha256_json  # noqa: E402
from phenocycler.config import load_config  # noqa: E402
from phenocycler.expression import read_single_partition  # noqa: E402
from phenocycler.pipeline import resolve_run_context  # noqa: E402


DEFAULT_DONORS = ("6476", "6539", "6450", "6436", "6591")
N_PER_DONOR = 500
STRATA = {
    "immune_called": 0.40,
    "islet_proximal": 0.25,
    "uniform": 0.35,
}
BROAD_LABELS = (
    "Immune",
    "Endothelial",
    "Epithelial",
    "Neural",
    "Mesenchymal",
    "Other",
)
ANNOTATION_VOCABULARY = (*BROAD_LABELS, "Unsure")
ACCEPTED_IMMUNE_STATUSES = frozenset({"anchor", "inferred"})
SAFE_ANNOTATOR = re.compile(r"^[A-Za-z0-9._-]+$")

CELL_COLUMNS = (
    "object_id",
    "image",
    "X_centroid",
    "Y_centroid",
    "cell_region",
    "islet_num",
    "cell_area",
    "cell_solidity",
)
ASSIGNMENT_COLUMNS = (
    "donor_id",
    "object_id",
    "broad_type",
    "broad_assignment_status",
    "best_broad_type",
    "confidence",
    "ranked_broad_probabilities",
    "broad_reason",
    "qc_analysis_eligible",
)
PREDICTION_RENAMES = {
    "broad_type": "pipeline_broad_type",
    "broad_assignment_status": "pipeline_broad_assignment_status",
    "best_broad_type": "pipeline_best_broad_type",
    "confidence": "pipeline_confidence",
    "ranked_broad_probabilities": "pipeline_ranked_broad_probabilities",
    "broad_reason": "pipeline_broad_reason",
    "qc_analysis_eligible": "pipeline_qc_analysis_eligible",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file_record(path: Path, *, root: Path | None = None) -> dict:
    resolved = path.resolve()
    root_resolved = root.resolve() if root is not None else None
    name = (
        str(resolved.relative_to(root_resolved))
        if root_resolved is not None and resolved.is_relative_to(root_resolved)
        else str(resolved)
    )
    return {
        "path": name,
        "size": resolved.stat().st_size,
        "sha256": _sha256(resolved),
    }


def _allocate_counts(n: int, fractions: Mapping[str, float]) -> dict[str, int]:
    """Allocate exactly ``n`` rows using deterministic largest remainders."""

    if n <= 0:
        raise ValueError("n_per_donor must be positive")
    if not fractions or any(value < 0 for value in fractions.values()):
        raise ValueError("stratum fractions must be non-negative")
    total = float(sum(fractions.values()))
    if not np.isclose(total, 1.0):
        raise ValueError(f"stratum fractions must sum to one, got {total}")

    raw = {name: n * float(value) for name, value in fractions.items()}
    allocated = {name: int(np.floor(value)) for name, value in raw.items()}
    remainder = n - sum(allocated.values())
    order = sorted(
        fractions,
        key=lambda name: (-(raw[name] - allocated[name]), tuple(fractions).index(name)),
    )
    for name in order[:remainder]:
        allocated[name] += 1
    if sum(allocated.values()) != n:
        raise AssertionError("stratum allocation did not preserve sample size")
    return allocated


def merge_population(
    cells: pd.DataFrame,
    assignments: pd.DataFrame,
    *,
    donor_id: str,
) -> pd.DataFrame:
    """Join the canonical cell context to the exact broad-type snapshot."""

    donor = str(donor_id)
    cells = cells.copy()
    assignments = assignments.copy()
    if "donor_id" not in cells:
        cells.insert(0, "donor_id", donor)

    missing_cells = {"donor_id", *CELL_COLUMNS}.difference(cells.columns)
    missing_assignments = set(ASSIGNMENT_COLUMNS).difference(assignments.columns)
    if missing_cells:
        raise KeyError(f"cells missing required audit columns: {sorted(missing_cells)}")
    if missing_assignments:
        raise KeyError(
            "assignments missing required audit columns: "
            f"{sorted(missing_assignments)}"
        )

    cells = cells.loc[:, ["donor_id", *CELL_COLUMNS]]
    assignments = assignments.loc[:, ASSIGNMENT_COLUMNS]
    for frame, name in ((cells, "cells"), (assignments, "assignments")):
        frame["donor_id"] = frame["donor_id"].astype(str)
        frame["object_id"] = frame["object_id"].astype(str)
        if not frame["donor_id"].eq(donor).all():
            raise ValueError(f"{name} contains rows outside donor {donor}")
        if frame.duplicated(["donor_id", "object_id"]).any():
            raise ValueError(f"{name} contains duplicate donor/object keys")

    merged = cells.merge(
        assignments,
        on=["donor_id", "object_id"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    if not merged["_merge"].eq("both").all():
        counts = merged["_merge"].value_counts().to_dict()
        raise ValueError(
            f"cell and assignment UUID universes differ for donor {donor}: {counts}"
        )
    merged = merged.drop(columns="_merge").rename(columns=PREDICTION_RENAMES)
    if merged["pipeline_broad_type"].isna().any():
        raise ValueError("the frozen broad prediction contains null labels")
    return merged


def assign_strata(population: pd.DataFrame) -> pd.Series:
    """Return disjoint audit strata without changing any production decision."""

    required = {
        "pipeline_broad_type",
        "pipeline_broad_assignment_status",
        "cell_region",
    }
    missing = required.difference(population.columns)
    if missing:
        raise KeyError(f"population missing stratum columns: {sorted(missing)}")

    immune = population["pipeline_broad_type"].eq("Immune") & population[
        "pipeline_broad_assignment_status"
    ].isin(ACCEPTED_IMMUNE_STATUSES)
    islet = population["cell_region"].isin({"core", "peri"})
    values = np.full(len(population), "uniform", dtype=object)
    values[islet.to_numpy()] = "islet_proximal"
    values[immune.to_numpy()] = "immune_called"
    return pd.Series(values, index=population.index, name="stratum")


def sample_population(
    population: pd.DataFrame,
    *,
    donor_id: str,
    n: int,
    rng: np.random.Generator,
    fractions: Mapping[str, float] = STRATA,
) -> pd.DataFrame:
    """Draw one donor's deterministic stratified sample with inverse weights."""

    population = population.copy()
    population["stratum"] = assign_strata(population)
    allocation = _allocate_counts(n, fractions)
    empty_allocations = [
        name for name, fraction in fractions.items()
        if fraction > 0 and allocation[name] == 0
    ]
    if empty_allocations:
        raise ValueError(
            "n_per_donor is too small to sample every positive audit stratum: "
            f"{empty_allocations}"
        )
    picks: list[pd.DataFrame] = []
    for name in fractions:
        pool = population.index[population["stratum"].eq(name)].to_numpy()
        want = allocation[name]
        if len(pool) < want:
            raise ValueError(
                f"donor {donor_id}/{name}: requested {want} rows but only "
                f"{len(pool)} are available"
            )
        chosen = rng.choice(pool, size=want, replace=False)
        selected = population.loc[chosen].copy()
        selected["stratum_N"] = len(pool)
        selected["stratum_n"] = want
        selected["weight"] = len(pool) / want
        picks.append(selected)

    sample = pd.concat(picks, ignore_index=True)
    sample["donor_N"] = len(population)
    if len(sample) != n:
        raise AssertionError(f"sampled {len(sample)} rows instead of {n}")
    reconstructed = sample.groupby("stratum", observed=False)["weight"].sum()
    expected = sample.groupby("stratum", observed=False)["stratum_N"].first()
    if not np.allclose(
        reconstructed.reindex(expected.index), expected.astype(float), rtol=0, atol=1e-9
    ):
        raise AssertionError("inverse sampling weights do not reconstruct strata")
    return sample


def blinded_review_index(sample: pd.DataFrame) -> pd.DataFrame:
    """Create a read-only location index containing no prediction or stratum."""

    return sample[
        ["review_id", "donor_id", "image", "X_centroid", "Y_centroid"]
    ].copy()


def _validated_context(config_path: Path):
    """Resolve one current run and validate its complete and type manifests."""

    config = load_config(config_path)
    context = resolve_run_context(config)
    current_broad = {rule.name for rule in context.typing_registry.broad_rules}
    expected_broad = set(BROAD_LABELS) - {"Other"}
    if current_broad != expected_broad:
        raise ValueError(
            "the audit vocabulary does not match the current broad typing rules: "
            f"audit={sorted(expected_broad)}, current={sorted(current_broad)}"
        )
    run_path = context.run_root / "run.json"
    type_path = context.stage_manifest_path("type")
    if not run_path.is_file():
        raise FileNotFoundError(
            f"validated gold-standard sampling requires a complete run: {run_path}"
        )
    if not type_path.is_file():
        raise FileNotFoundError(f"validated type-stage manifest is missing: {type_path}")

    run_manifest = RunManifest.read_json(run_path)
    if run_manifest.run_name != context.run_id:
        raise ValueError(
            f"run manifest names {run_manifest.run_name}, expected {context.run_id}"
        )
    run_manifest.validate_current(validation_mode="fast")

    type_manifest = StageManifest.read_json(type_path)
    if type_manifest.stage != "type":
        raise ValueError(f"{type_path} declares stage {type_manifest.stage!r}")
    type_manifest.validate_current(validation_mode="fast")
    if type_manifest.expected_donors != context.donors:
        raise ValueError("type-stage donor contract differs from the resolved cohort")
    references = {reference.stage: reference for reference in run_manifest.stages}
    if (
        "type" not in references
        or references["type"].manifest_content_id != type_manifest.content_id
    ):
        raise ValueError("run manifest does not own the validated type stage")
    return config, context, run_manifest, type_manifest


def _normalise_donors(values: Iterable[str], expected: Iterable[str]) -> tuple[str, ...]:
    donors = tuple(map(str, values))
    if not donors or len(set(donors)) != len(donors):
        raise ValueError("audit donors must be a non-empty unique list")
    unknown = sorted(set(donors).difference(map(str, expected)))
    if unknown:
        raise ValueError(f"audit donors are absent from the current run: {unknown}")
    return donors


def _normalise_annotators(values: Iterable[str]) -> tuple[str, str]:
    annotators = tuple(map(str, values))
    if len(annotators) != 2 or len(set(annotators)) != 2:
        raise ValueError("exactly two distinct annotator identifiers are required")
    unsafe = [value for value in annotators if SAFE_ANNOTATOR.fullmatch(value) is None]
    if unsafe:
        raise ValueError(
            "annotator identifiers may contain only letters, numbers, '.', '_', "
            f"or '-': {unsafe}"
        )
    return annotators[0], annotators[1]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=REPO / "config.ini")
    parser.add_argument("--donor", action="append", default=None)
    parser.add_argument("--n-per-donor", type=int, default=N_PER_DONOR)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--annotators", nargs="+", default=["A", "B"])
    parser.add_argument("--seed", type=int, default=20260727)
    args = parser.parse_args(argv)

    config, context, run_manifest, type_manifest = _validated_context(args.config)
    donors = _normalise_donors(args.donor or DEFAULT_DONORS, context.donors)
    annotators = _normalise_annotators(args.annotators)

    out = (
        args.out
        if args.out is not None
        else config.data_dir / "gold_standard_v2" / context.run_id
    ).expanduser().resolve()
    if out.exists():
        raise FileExistsError(
            f"audit bundle path already exists and will not be overwritten: {out}"
        )

    rng = np.random.default_rng(args.seed)
    frames: list[pd.DataFrame] = []
    for donor in donors:
        cells = read_single_partition(context.config.cells_dir, donor)
        assignments = read_single_partition(context.config.assignments_dir, donor)
        population = merge_population(cells, assignments, donor_id=donor)
        frames.append(
            sample_population(
                population,
                donor_id=donor,
                n=args.n_per_donor,
                rng=rng,
            )
        )
    sample = pd.concat(frames, ignore_index=True)
    sample = sample.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
    sample.insert(0, "review_id", [f"GS{i:05d}" for i in range(len(sample))])

    out.mkdir(parents=True)
    private_columns = (
        "review_id",
        "donor_id",
        "object_id",
        "image",
        "X_centroid",
        "Y_centroid",
        "cell_region",
        "islet_num",
        "cell_area",
        "cell_solidity",
        "stratum",
        "stratum_N",
        "stratum_n",
        "weight",
        "donor_N",
        *PREDICTION_RENAMES.values(),
    )
    sample_key = out / "sample_key.csv"
    sample.loc[:, private_columns].to_csv(sample_key, index=False)

    blinded_review_index(sample).to_csv(
        out / "blinded_review_index.csv", index=False
    )

    qupath = out / "qupath"
    qupath.mkdir()
    for donor, donor_sample in sample.groupby("donor_id", sort=True):
        donor_sample[["object_id", "review_id", "image"]].to_csv(
            qupath / f"gold_standard_{donor}.csv", index=False
        )

    generated = sorted(
        path for path in out.rglob("*") if path.is_file() and path.name != "manifest.json"
    )
    source = {
        "config": _file_record(args.config.expanduser().resolve(), root=REPO),
        "run_manifest": _file_record(context.run_root / "run.json", root=REPO),
        "type_stage_manifest": _file_record(
            context.stage_manifest_path("type"), root=REPO
        ),
        "run_id": context.run_id,
        "run_manifest_content_id": run_manifest.content_id,
        "type_stage_content_id": type_manifest.content_id,
        "type_output_content_sha256": type_manifest.output.content_sha256,
        "cohort_content_id": context.cohort.content_id,
        "marker_registry_fingerprint": context.registry.fingerprint,
        "typing_rules_fingerprint": context.typing_registry.rules_fingerprint,
    }
    manifest = {
        "schema_version": 2,
        "protocol_version": "gold-standard-broad-audit-v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "audit_only": (
            "Human labels are blinded post hoc audit data. They are forbidden as "
            "calibration, typing, configuration, model-selection, or release-gate inputs."
        ),
        "source": source,
        "sampling": {
            "seed": args.seed,
            "donors": list(donors),
            "n_per_donor": args.n_per_donor,
            "n_total": int(len(sample)),
            "strata": STRATA,
            "immune_called": (
                "pipeline_broad_type == 'Immune' and broad status is anchor or inferred"
            ),
            "islet_proximal": (
                "non-immune-called cells with QuPath-derived cell_region core or peri"
            ),
            "weighting": (
                "weight=stratum_N/stratum_n; estimates describe the selected-donor "
                "canonical-cell population, not the full cohort"
            ),
        },
        "vocabulary": list(ANNOTATION_VOCABULARY),
        "prediction_abstentions": ["Ambiguous", "Unavailable"],
        "annotators": list(annotators),
        "blinding": (
            "Annotator-facing files contain location and review ID only; the private "
            "sample key contains frozen pipeline predictions and must not be shared."
        ),
        "files": [_file_record(path, root=out) for path in generated],
    }
    manifest["content_id"] = sha256_json(manifest)
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary = sample.groupby(["donor_id", "stratum"]).size().unstack(fill_value=0)
    print(summary.to_string())
    print(f"\nsampled {len(sample):,} cells from validated run {context.run_id}")
    print(f"[done] immutable audit bundle: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
