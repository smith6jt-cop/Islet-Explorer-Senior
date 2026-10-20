"""Strict, reviewable phenotype label spaces and identity releases.

The integration layer never imports a producer's Python vocabulary. A producing lane ships
JSON declarations, and this module validates them without mutable module-level state. An
undeclared label, mismatched expression content id, excluded donor, or retired input is a hard
error. Geometry remains usable when no release is configured.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional

import pandas as pd

from ..cohort import DONOR_EXCLUSIONS
from ..config import PipelineConfig


SCHEMA_VERSION = 1
COMMON_RESOLUTIONS = frozenset({"mapped", "ambiguous", "unavailable", "not_mappable"})
RETIRED_INPUT_TOKENS = (
    "restore_gated_redsea",
    "phenotype/broad",
    "production-default-restore-ordered-residual-typing-v1",
)


class PhenotypeSpaceError(ValueError):
    """A label-space or identity-release declaration is invalid."""


@dataclass(frozen=True)
class LabelRule:
    label: str
    parent: str = ""
    common_lineage: str = ""
    common_resolution: str = "mapped"
    inherit: bool = False
    reason: str = ""


@dataclass(frozen=True)
class PhenotypeSpace:
    schema_version: int
    space_id: str
    producer_lane: str
    rules: Mapping[str, LabelRule]
    source: Path

    def common_for(self, label: str) -> tuple[str, str]:
        """Resolve a producer label by walking declared parents, never by prefix."""

        current = str(label).strip()
        if not current:
            return "", "unavailable"
        seen: set[str] = set()
        while current:
            if current in seen:
                raise PhenotypeSpaceError(
                    f"{self.space_id}: parent cycle while resolving {label!r}: {sorted(seen)}"
                )
            seen.add(current)
            rule = self.rules.get(current)
            if rule is None:
                raise PhenotypeSpaceError(
                    f"{self.space_id}: undeclared producer label {current!r}"
                )
            if not rule.inherit:
                return rule.common_lineage, rule.common_resolution
            current = rule.parent
        raise PhenotypeSpaceError(
            f"{self.space_id}: label {label!r} inherits but has no declared parent mapping"
        )


@dataclass(frozen=True)
class PhenotypeRelease:
    schema_version: int
    release_id: str
    modality: str
    producer_lane: str
    producer_commit: str
    run_id: str
    expression_content_id: str
    cell_key: str
    label_column: str
    common_lineage_column: str
    evidence_column: str
    confidence_column: str
    fine_type_column: str
    calls_dir: Path
    donors: tuple[str, ...]
    space: PhenotypeSpace
    source: Path
    provenance: str = "unchecked"
    """Outcome of checking the declared expression input, never a content hash.

    ``verified:<path>``            the declared input exists and covers every donor claimed
    ``unverifiable:<reason>``      the declared input does not resolve, so nothing was checked

    This exists because the old run_id + content_id pin was checked against a stage manifest,
    and stage manifests were removed on 2026-09-16. A release naming a wrong content id would
    otherwise load with no error and no trace -- the failure mode being silent was the problem,
    not the failure mode being possible. A consumer that needs provenance reads this field;
    `integration --status` prints it.

    It deliberately cannot fire because package source changed: it reads only the filesystem
    the release itself names.
    """

    def calls_for(self, donor_id: str, roi: str = "") -> list[Path]:
        donor = str(donor_id)
        part = self.calls_dir / f"donor_id={donor}"
        if not part.is_dir():
            return []
        if roi:
            roi_part = part / f"roi={roi}"
            if roi_part.is_dir():
                return sorted(roi_part.rglob("*.parquet"))
            # Once a release partitions one donor by ROI, absence of the requested ROI is
            # absence of calls. Do not recurse through sibling ROIs and rely on a later cell
            # mismatch to catch the mistake.
            if any(part.glob("roi=*")):
                return []
        # A donor-level release follows the declared A.4 layout and stores shards directly
        # below donor_id=*. Nested ROI partitions were handled above.
        return sorted(part.glob("*.parquet"))


def _read_json(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise PhenotypeSpaceError(f"cannot read {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise PhenotypeSpaceError(f"{path} must contain one JSON object")
    return payload


def load_phenotype_space(path: str | Path) -> PhenotypeSpace:
    source = Path(path).resolve()
    payload = _read_json(source)
    version = int(payload.get("schema_version", -1))
    if version != SCHEMA_VERSION:
        raise PhenotypeSpaceError(
            f"{source}: unsupported label-space schema_version {version}; expected {SCHEMA_VERSION}"
        )
    space_id = str(payload.get("space_id", "")).strip()
    lane = str(payload.get("producer_lane", "")).strip()
    entries = payload.get("labels")
    if not space_id or not lane or not isinstance(entries, list):
        raise PhenotypeSpaceError(
            f"{source}: space_id, producer_lane and labels[] are required"
        )

    rules: dict[str, LabelRule] = {}
    for raw in entries:
        if not isinstance(raw, dict):
            raise PhenotypeSpaceError(f"{source}: every labels[] entry must be an object")
        label = str(raw.get("label", "")).strip()
        if not label or label in rules:
            raise PhenotypeSpaceError(f"{source}: blank or duplicate label {label!r}")
        resolution = str(raw.get("common_resolution", "mapped")).strip()
        if resolution not in COMMON_RESOLUTIONS:
            raise PhenotypeSpaceError(
                f"{source}: {label!r} has unknown common_resolution {resolution!r}"
            )
        rule = LabelRule(
            label=label,
            parent=str(raw.get("parent", "")).strip(),
            common_lineage=str(raw.get("common_lineage", "") or "").strip(),
            common_resolution=resolution,
            inherit=bool(raw.get("inherit", False)),
            reason=str(raw.get("reason", "")).strip(),
        )
        if rule.inherit and not rule.parent:
            raise PhenotypeSpaceError(f"{source}: {label!r} inherits but has no parent")
        if resolution == "mapped" and not rule.inherit and not rule.common_lineage:
            raise PhenotypeSpaceError(
                f"{source}: mapped label {label!r} has no common_lineage"
            )
        if resolution != "mapped" and rule.common_lineage:
            raise PhenotypeSpaceError(
                f"{source}: {resolution} label {label!r} cannot carry common_lineage"
            )
        rules[label] = rule

    for rule in rules.values():
        if rule.parent and rule.parent not in rules:
            raise PhenotypeSpaceError(
                f"{source}: {rule.label!r} names undeclared parent {rule.parent!r}"
            )
    space = PhenotypeSpace(version, space_id, lane, rules, source)
    for label in rules:
        space.common_for(label)
    return space


def _configured_path(cfg: PipelineConfig, modality: str) -> Optional[Path]:
    name = {
        "phenocycler": "phenocycler_identity_release",
        "xenium": "xenium_identity_release",
    }.get(modality)
    if name is None:
        raise PhenotypeSpaceError(f"unknown identity modality {modality!r}")
    value = getattr(cfg, name)
    return Path(value).expanduser().resolve() if value else None


def load_release(cfg: PipelineConfig, modality: str) -> PhenotypeRelease:
    configured = _configured_path(cfg, modality)
    if configured is None:
        raise PhenotypeSpaceError(f"no declared {modality} identity release is configured")
    if not configured.exists():
        # config.ini may name the release directory a producer lane WILL write (the
        # broadtype A.4a release, the Xenium typer release) before either exists. For the
        # strict loader that is still "no release": identity stages must fail closed on it.
        raise PhenotypeSpaceError(
            f"{modality} identity release is configured at {configured} but nothing exists "
            "there yet; the producer lane has not written its release"
        )
    source = configured / "release.json" if configured.is_dir() else configured
    payload = _read_json(source)
    version = int(payload.get("schema_version", -1))
    if version != SCHEMA_VERSION:
        raise PhenotypeSpaceError(
            f"{source}: unsupported release schema_version {version}; expected {SCHEMA_VERSION}"
        )

    release_id = str(payload.get("release_id", "")).strip()
    got_modality = str(payload.get("modality", "")).strip()
    producer = payload.get("producer", {})
    lane = str(producer.get("lane", "")).strip() if isinstance(producer, dict) else ""
    commit = str(producer.get("commit", "")).strip() if isinstance(producer, dict) else ""
    if not release_id or got_modality != modality or not lane or not commit:
        raise PhenotypeSpaceError(
            f"{source}: release_id, matching modality, producer.lane and producer.commit "
            "are required"
        )

    base = source.parent
    label_ref = Path(str(payload.get("label_space", "label_space.json")))
    calls_ref = Path(str(payload.get("calls_dir", "calls")))
    if label_ref.is_absolute() or calls_ref.is_absolute():
        raise PhenotypeSpaceError(
            f"{source}: label_space and calls_dir must stay inside the release directory"
        )
    space_path = (base / label_ref).resolve()
    calls_dir = (base / calls_ref).resolve()
    for declared_path, name in ((space_path, "label_space"), (calls_dir, "calls_dir")):
        try:
            declared_path.relative_to(base.resolve())
        except ValueError as exc:
            raise PhenotypeSpaceError(
                f"{source}: {name} escapes the release directory: {declared_path}"
            ) from exc
    space = load_phenotype_space(space_path)
    if space.producer_lane != lane:
        raise PhenotypeSpaceError(
            f"{source}: release lane {lane!r} != label-space lane {space.producer_lane!r}"
        )

    donors = tuple(str(d).strip() for d in payload.get("donors", []) if str(d).strip())
    if not donors or len(set(donors)) != len(donors):
        raise PhenotypeSpaceError(f"{source}: donors must be a non-empty, duplicate-free list")
    blocked = sorted(set(donors) & set(DONOR_EXCLUSIONS))
    if blocked:
        reasons = "; ".join(f"{d}: {DONOR_EXCLUSIONS[d]}" for d in blocked)
        raise PhenotypeSpaceError(f"{source}: release includes excluded donor(s): {reasons}")

    inputs = payload.get("inputs", [])
    if not isinstance(inputs, list):
        raise PhenotypeSpaceError(f"{source}: inputs must be a list")
    combined_provenance = " ".join(
        [lane, str(source), str(space_path), str(calls_dir)]
        + [str(item.get("path", "")) for item in inputs if isinstance(item, dict)]
    ).lower()
    if any(token.lower() in combined_provenance for token in RETIRED_INPUT_TOKENS):
        raise PhenotypeSpaceError(f"{source}: release names a retired identity input")
    if "/assignments/" in combined_provenance or "\\assignments\\" in combined_provenance:
        raise PhenotypeSpaceError(
            f"{source}: runs/<id>/assignments is retired and cannot be an identity release"
        )

    expression = next(
        (item for item in inputs if isinstance(item, dict) and item.get("kind") == "expression"),
        {},
    )
    run_id = str(expression.get("run_id", payload.get("run_id", ""))).strip()
    content_id = str(
        expression.get("content_id", payload.get("expression_content_id", ""))
    ).strip()
    # A PhenoCycler release used to pin an expression run_id + content_id, checked against
    # data/runs/<run_id>/manifests/expression.json. Stage manifests were removed on
    # 2026-09-16, so there is nothing left to compare a content id against.
    #
    # What replaces it is coverage, not identity: does the declared input exist, and does it
    # actually contain every donor this release claims to call? That catches the realistic
    # wrong-id release -- one pointing at a tree that is gone, or claiming donors the data
    # does not have -- without hashing anything, and it cannot fire because code changed.
    #
    # A release whose input does not resolve is NOT refused. It is recorded as unverifiable,
    # because refusing on a path that may simply be configured elsewhere is how the old check
    # ended up failing for reasons unrelated to provenance. Silence was the defect; a hard
    # gate is not the only cure for it.
    provenance = _check_declared_input(expression, donors, source)

    columns = payload.get("columns", {})
    if not isinstance(columns, dict):
        raise PhenotypeSpaceError(f"{source}: columns must be an object")
    cell_key = str(payload.get("cell_key", "object_id")).strip()
    expected_key = "object_id" if modality == "phenocycler" else "cell_id"
    if cell_key != expected_key:
        raise PhenotypeSpaceError(
            f"{source}: {modality} cell_key must be {expected_key!r}, got {cell_key!r}"
        )
    return PhenotypeRelease(
        schema_version=version,
        release_id=release_id,
        modality=got_modality,
        producer_lane=lane,
        producer_commit=commit,
        run_id=run_id,
        expression_content_id=content_id,
        cell_key=cell_key,
        label_column=str(columns.get("label", "lineage_native")).strip(),
        common_lineage_column=str(columns.get("common_lineage", "common_lineage")).strip(),
        evidence_column=str(columns.get("evidence_positive", "evidence_positive")).strip(),
        confidence_column=str(columns.get("confidence", "confidence")).strip(),
        fine_type_column=str(columns.get("fine_type", "celltype_fine")).strip(),
        calls_dir=calls_dir,
        donors=donors,
        space=space,
        source=source,
        provenance=provenance,
    )


# core/multiplex5000/multiplex5000/integration/this_file.py -> parents[4] is the repository root.
# parents[3] is core/, which resolves a repo-relative declared path only when the process
# happens to be running from the repo root -- so --status reported "unverifiable" for a
# release that is in fact verifiable.
_REPO_ROOT = Path(__file__).resolve().parents[4]


def _check_declared_input(
    expression: Mapping[str, object],
    donors: tuple[str, ...],
    source: Path,
) -> str:
    """Check the declared expression input by coverage. Never by content hash.

    Returns a ``verified:`` or ``unverifiable:`` string for
    :attr:`PhenotypeRelease.provenance`. Raises only when the input resolves AND
    contradicts the release -- a release that claims donors its own declared input does
    not contain is wrong in a way no configuration explains.
    """

    declared = str(expression.get("path", "")).strip()
    if not declared:
        return "unverifiable:release declares no expression input path"

    candidates = [Path(declared)]
    if not Path(declared).is_absolute():
        candidates.append(_REPO_ROOT / declared)
    root = next((c for c in candidates if c.is_dir()), None)
    if root is None:
        # Expected for a release whose input lives on another machine, and for the
        # synthetic releases in the test-suite. Recorded, not raised.
        return f"unverifiable:declared input {declared} does not resolve"

    present = {
        p.name.split("=", 1)[1]
        for p in root.glob("donor_id=*")
        if "=" in p.name
    }
    # donor_id is not normalised across trees: 0112 and 112 are one donor.
    normalised = {d.lstrip("0") or "0" for d in present}
    missing = [
        d for d in donors
        if d not in present and (d.lstrip("0") or "0") not in normalised
    ]
    if missing:
        raise PhenotypeSpaceError(
            f"{source}: release claims donors its declared expression input does not "
            f"contain: missing={sorted(missing)} under {root}. Either the release names "
            "the wrong input or it claims calls it cannot have made."
        )
    return f"verified:{root}"


def load_release_or_none(cfg: PipelineConfig, modality: str) -> Optional[PhenotypeRelease]:
    """Return None when no release is *available*; a configured invalid release still raises.

    Two absences are treated alike: no path configured, and a configured path that does not
    exist yet. The second is the normal state while a producer lane is still running --
    ``config.ini`` names the directory it will write so that geometry stages attach identity
    the moment it lands, without a config edit. A directory that exists but holds an invalid
    release is different: that is a broken declaration, and it raises.
    """

    configured = _configured_path(cfg, modality)
    if configured is None:
        return None
    if not configured.exists():
        print(f"[phenospace] {modality} identity release configured at {configured} does not "
              "exist yet; continuing geometry only", flush=True)
        return None
    return load_release(cfg, modality)


def identity_attached(cfg: PipelineConfig) -> bool:
    """True only when BOTH declared releases actually load.

    ``bool(cfg.phenocycler_identity_release and cfg.xenium_identity_release)`` was the
    previous test, and it became wrong the moment ``config.ini`` started naming release
    directories ahead of their producers: a configured path is a promise, not an attached
    identity. A configured release that fails validation is reported and counted as not
    attached, so a geometry-only report is never blocked by a broken declaration.
    """

    for modality in ("phenocycler", "xenium"):
        try:
            if load_release_or_none(cfg, modality) is None:
                return False
        except PhenotypeSpaceError as exc:
            print(f"[phenospace] {modality} identity release is configured but invalid; "
                  f"treating identity as not attached: {exc}", flush=True)
            return False
    return True


def require_declared_phenotype(cfg: PipelineConfig, modality: str) -> PhenotypeRelease:
    """Load the configured release or fail closed with one stable error type."""

    return load_release(cfg, modality)


def attach_release_identity(
    geometry: pd.DataFrame,
    release: PhenotypeRelease,
    donor_id: str,
) -> pd.DataFrame:
    """Attach a declared release by exact native cell id, preserving the geometry frame.

    This is modality-neutral.  A release may reorder calls, but it may not omit or add a
    cell.  Producer-supplied common labels are checked against the strict label-space
    declaration rather than trusted as a second vocabulary authority.
    """

    from .contract import (
        ContractError,
        coerce_cell_table,
        contract_profile,
        validate_cell_table,
    )

    validate_cell_table(geometry, profile="geometry")
    if contract_profile(geometry) != "geometry":
        raise PhenotypeSpaceError("identity attachment requires a geometry-profile table")
    modalities = set(geometry["modality"].astype(str))
    if modalities != {release.modality}:
        raise PhenotypeSpaceError(
            f"cannot attach {release.modality!r} release to table modality {sorted(modalities)}"
        )

    donor = str(donor_id)
    if release.donors and donor not in set(release.donors):
        raise PhenotypeSpaceError(
            f"identity release {release.release_id} does not declare donor {donor}"
        )
    roi_values = sorted(set(geometry["roi"].astype(str)))
    if len(roi_values) != 1:
        raise PhenotypeSpaceError(
            f"identity attachment requires one ROI, got {roi_values}"
        )
    roi = roi_values[0]
    files = release.calls_for(donor, roi)
    if not files:
        raise PhenotypeSpaceError(
            f"identity release {release.release_id} has no calls for {donor}/{roi} under "
            f"{release.calls_dir}"
        )
    calls = pd.concat([pd.read_parquet(path) for path in files], ignore_index=True)
    key = release.cell_key
    if key not in calls.columns:
        raise PhenotypeSpaceError(f"identity calls for donor {donor} have no cell_key {key!r}")
    call_ids = calls[key].astype(str)
    if call_ids.duplicated().any():
        raise PhenotypeSpaceError(
            f"identity calls for donor {donor} contain {int(call_ids.duplicated().sum())} "
            f"duplicate {key} values"
        )
    geometry_ids = geometry["cell_id"].astype(str)
    call_set = set(call_ids)
    geometry_set = set(geometry_ids)
    missing_calls = geometry_ids[~geometry_ids.isin(call_set)]
    extra_calls = call_ids[~call_ids.isin(geometry_set)]
    if len(missing_calls) or len(extra_calls):
        raise PhenotypeSpaceError(
            f"identity/geometry UUID mismatch (native cell-id mismatch) for donor "
            f"{donor}/{roi}: "
            f"{len(missing_calls)} geometry cells lack calls; {len(extra_calls)} calls lack "
            "geometry. Refusing to change the frame."
        )

    indexed = calls.assign(**{key: call_ids}).set_index(key).reindex(geometry_ids)
    label_col = release.label_column
    if label_col not in indexed.columns:
        raise PhenotypeSpaceError(
            f"identity calls have no declared label column {label_col!r}"
        )
    labels = indexed[label_col].fillna("").astype(str)
    mapped = [release.space.common_for(label) for label in labels]
    declared_common = [common for common, _resolution in mapped]
    resolutions = [resolution for _common, resolution in mapped]

    common_col = release.common_lineage_column
    if common_col and common_col in indexed.columns:
        supplied = indexed[common_col].fillna("").astype(str).tolist()
        mismatches = sum(a != b for a, b in zip(supplied, declared_common))
        if mismatches:
            raise PhenotypeSpaceError(
                f"identity release {release.release_id} has {mismatches} common-lineage "
                "values that disagree with its label_space.json"
            )

    out = geometry.copy()
    out["lineage_native"] = labels.to_numpy()
    out["lineage_common"] = declared_common
    out["common_resolution"] = resolutions
    for target, declared in (
        ("celltype_fine", release.fine_type_column),
        ("confidence", release.confidence_column),
        ("evidence_positive", release.evidence_column),
    ):
        if declared and declared in indexed.columns:
            out[target] = indexed[declared].to_numpy()
        elif target == "evidence_positive":
            out[target] = pd.Series(pd.NA, index=out.index, dtype="boolean")

    try:
        return coerce_cell_table(
            out,
            modality=release.modality,
            donor_id=donor,
            roi=roi,
            tissue=str(out["tissue"].iloc[0]) if len(out) else "",
            profile="identity",
            identity_release_id=release.release_id,
        )
    except ContractError as exc:
        raise PhenotypeSpaceError(
            f"identity release {release.release_id} violates the contract: {exc}"
        ) from exc


__all__ = [
    "LabelRule",
    "PhenotypeRelease",
    "PhenotypeSpace",
    "PhenotypeSpaceError",
    "attach_release_identity",
    "identity_attached",
    "load_phenotype_space",
    "load_release",
    "load_release_or_none",
    "require_declared_phenotype",
]
