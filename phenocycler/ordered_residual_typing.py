"""RESTORE broad phenotyping by the validated ordered residual hierarchy.

This is the pipeline-native form of the historical broad-lineage stage.  It
consumes RESTORE's categorical per-marker calls directly; it does not fit a
classifier, estimate probabilities, or require reviewed cell labels.

For each cell, gates are visited in their versioned priority order.  The first
gate containing a positive marker claims the cell.  If a gate has no positive
marker but at least one unavailable/indeterminate marker, that gate blocks all
lower-priority gates and the cell remains unavailable.  A cell with explicit
negative evidence at every gate is ``Other``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from .artifacts import ContractError, sha256_json
from .marker_registry import MarkerRegistry, load_registry


METHOD_VERSION = "restore-ordered-residual-typing-v1"
TYPING_MODE = "restore_ordered_residual"
SCORE_SEMANTICS = "categorical_restore_tristate"
DEFAULT_RULES_PATH = Path(__file__).with_name("ordered_residual_rules.json")

POSITIVE = "positive"
NEGATIVE = "negative"
UNAVAILABLE = "unavailable"
_UNAVAILABLE_INPUT_STATES = frozenset(
    {"", "unavailable", "indeterminate", "unknown", "missing"}
)


@dataclass(frozen=True)
class OrderedResidualPolicy:
    """Validated, content-fingerprinted gate order and marker membership."""

    schema_version: int
    rules_version: str
    description: str
    compartment_order: tuple[str, ...]
    gates: Mapping[str, tuple[str, ...]]
    other_label: str
    unavailable_label: str
    fingerprint: str
    source_path: Path

    @property
    def markers(self) -> tuple[str, ...]:
        return tuple(
            marker
            for compartment in self.compartment_order
            for marker in self.gates[compartment]
        )

    def validate_registry(self, registry: MarkerRegistry) -> None:
        """Require the gate web to be exactly the validated RESTORE targets."""

        active = {marker.name for marker in registry.active_markers}
        gated = set(self.markers)
        missing = sorted(active - gated)
        extra = sorted(gated - active)
        if missing or extra:
            raise ContractError(
                "ordered residual gates must equal the active RESTORE target "
                f"set (missing={missing}, extra={extra})"
            )

    @classmethod
    def read_json(cls, path: Path | str = DEFAULT_RULES_PATH) -> "OrderedResidualPolicy":
        source = Path(path).expanduser().resolve()
        try:
            raw = json.loads(source.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ContractError(
                f"cannot read ordered residual rules from {source}: {exc}"
            ) from exc
        if not isinstance(raw, dict):
            raise ContractError("ordered residual rules must be a JSON object")
        required = {
            "schema_version",
            "rules_version",
            "description",
            "compartment_order",
            "gates",
            "other_label",
            "unavailable_label",
        }
        missing_keys = sorted(required - set(raw))
        unknown_keys = sorted(set(raw) - required)
        if missing_keys or unknown_keys:
            raise ContractError(
                "ordered residual rule keys differ from the contract "
                f"(missing={missing_keys}, unknown={unknown_keys})"
            )
        if int(raw["schema_version"]) != 1:
            raise ContractError("ordered residual schema_version must be 1")
        order = tuple(str(value).strip() for value in raw["compartment_order"])
        if not order or any(not value for value in order):
            raise ContractError("compartment_order must contain nonempty labels")
        if len(order) != len(set(order)):
            raise ContractError("compartment_order contains duplicate labels")
        raw_gates = raw["gates"]
        if not isinstance(raw_gates, dict) or set(raw_gates) != set(order):
            raise ContractError(
                "ordered residual gate keys must exactly equal compartment_order"
            )
        gates: dict[str, tuple[str, ...]] = {}
        for compartment in order:
            values = tuple(str(value).strip() for value in raw_gates[compartment])
            if not values or any(not value for value in values):
                raise ContractError(
                    f"ordered residual gate {compartment!r} must be nonempty"
                )
            if len(values) != len(set(values)):
                raise ContractError(
                    f"ordered residual gate {compartment!r} has duplicate markers"
                )
            gates[compartment] = values
        markers = [marker for values in gates.values() for marker in values]
        if len(markers) != len(set(markers)):
            raise ContractError(
                "a RESTORE marker may occur in only one ordered residual gate"
            )
        other = str(raw["other_label"]).strip()
        unavailable = str(raw["unavailable_label"]).strip()
        if not other or not unavailable or other == unavailable:
            raise ContractError(
                "other_label and unavailable_label must be distinct and nonempty"
            )
        if other in order or unavailable in order:
            raise ContractError(
                "terminal ordered residual labels cannot also be compartments"
            )
        canonical = {
            "schema_version": 1,
            "rules_version": str(raw["rules_version"]).strip(),
            "description": str(raw["description"]).strip(),
            "compartment_order": order,
            "gates": gates,
            "other_label": other,
            "unavailable_label": unavailable,
        }
        if not canonical["rules_version"]:
            raise ContractError("ordered residual rules_version must be nonempty")
        return cls(
            **canonical,
            fingerprint=sha256_json(canonical),
            source_path=source,
        )


def load_ordered_residual_policy(
    path: Path | str = DEFAULT_RULES_PATH,
    *,
    registry: MarkerRegistry | None = None,
) -> OrderedResidualPolicy:
    policy = OrderedResidualPolicy.read_json(path)
    policy.validate_registry(load_registry() if registry is None else registry)
    return policy


def _canonical_state(series: pd.Series, *, marker: str) -> np.ndarray:
    raw = series.astype("string").fillna("").str.strip().str.lower()
    unknown = sorted(
        set(raw) - {POSITIVE, NEGATIVE, *_UNAVAILABLE_INPUT_STATES}
    )
    if unknown:
        raise ContractError(
            f"RESTORE marker {marker!r} has unknown cell states: {unknown}"
        )
    state = raw.to_numpy(dtype=object)
    unavailable = np.isin(state, tuple(_UNAVAILABLE_INPUT_STATES))
    state[unavailable] = UNAVAILABLE
    return state


def type_ordered_residual_cells(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None = None,
    policy: OrderedResidualPolicy | None = None,
) -> pd.DataFrame:
    """Assign the first positive RESTORE gate, failing closed on unavailable gates."""

    registry = load_registry() if marker_registry is None else marker_registry
    active_policy = (
        load_ordered_residual_policy(registry=registry)
        if policy is None
        else policy
    )
    active_policy.validate_registry(registry)
    required = {"donor_id", "object_id", "qc_analysis_eligible"}
    required.update(f"{marker}__state" for marker in active_policy.markers)
    missing = sorted(required - set(evidence_wide))
    if missing:
        raise ContractError(
            f"RESTORE ordered typing evidence is missing columns: {missing}"
        )
    keys = evidence_wide.loc[:, ["donor_id", "object_id"]].copy()
    if keys.isna().any().any():
        raise ContractError("RESTORE ordered typing keys cannot be null")
    keys["donor_id"] = keys["donor_id"].astype(str)
    keys["object_id"] = keys["object_id"].astype(str)
    if keys.duplicated(["donor_id", "object_id"]).any():
        raise ContractError("RESTORE ordered typing contains duplicate cells")
    eligibility = evidence_wide["qc_analysis_eligible"]
    if eligibility.isna().any() or not eligibility.isin([True, False, 0, 1]).all():
        raise ContractError("qc_analysis_eligible must contain explicit booleans")
    eligible = eligibility.astype(bool).to_numpy()
    n = len(evidence_wide)

    states = {
        marker: _canonical_state(
            evidence_wide[f"{marker}__state"], marker=marker
        )
        for marker in active_policy.markers
    }
    all_required_valid = eligible.copy()
    all_required_negative = eligible.copy()
    for state in states.values():
        all_required_valid &= state != UNAVAILABLE
        all_required_negative &= state == NEGATIVE

    label = np.full(n, active_policy.other_label, dtype=object)
    status = np.full(n, "other", dtype=object)
    reason = np.full(n, "all_negative", dtype=object)
    blocked_by = np.full(n, "", dtype=object)
    decisive_gate = np.full(n, "", dtype=object)
    authoritative_marker = np.full(n, "", dtype=object)
    best_candidate = np.full(n, "", dtype=object)
    decided = ~eligible

    for compartment in active_policy.compartment_order:
        gate_positive = np.zeros(n, dtype=bool)
        gate_unavailable = np.zeros(n, dtype=bool)
        first_positive = np.full(n, "", dtype=object)
        for marker in active_policy.gates[compartment]:
            state = states[marker]
            marker_positive = state == POSITIVE
            first_positive[(first_positive == "") & marker_positive] = marker
            gate_positive |= marker_positive
            gate_unavailable |= state == UNAVAILABLE
        # A positive sibling proves gate membership and dominates an unavailable
        # sibling exactly as in the validated historical implementation.
        gate_unavailable &= ~gate_positive
        take_positive = ~decided & gate_positive
        take_unavailable = ~decided & gate_unavailable
        label[take_positive] = compartment
        status[take_positive] = "anchor"
        reason[take_positive] = "positive"
        decisive_gate[take_positive] = compartment
        authoritative_marker[take_positive] = first_positive[take_positive]
        best_candidate[take_positive] = compartment
        label[take_unavailable] = active_policy.unavailable_label
        status[take_unavailable] = "unavailable"
        reason[take_unavailable] = "blocked_unavailable_gate"
        blocked_by[take_unavailable] = compartment
        decisive_gate[take_unavailable] = compartment
        best_candidate[take_unavailable] = compartment
        decided |= take_positive | take_unavailable

    label[~eligible] = active_policy.unavailable_label
    status[~eligible] = "unavailable"
    reason[~eligible] = "geometry_qc_ineligible"
    blocked_by[~eligible] = ""
    decisive_gate[~eligible] = ""
    authoritative_marker[~eligible] = ""
    best_candidate[~eligible] = ""

    deterministic = np.isin(status, ("anchor", "other"))
    confidence = np.full(n, np.nan, dtype=float)
    confidence[deterministic] = 1.0
    subtype_status = np.full(n, "not_applicable", dtype=object)
    subtype_status[~eligible] = "unavailable"
    subtype_reason = np.full(n, "broad_only_ordered_residual", dtype=object)
    subtype_reason[~eligible] = "geometry_qc_ineligible"
    empty = np.full(n, "", dtype=object)

    out = keys.copy()
    out["broad_label"] = label
    out["broad_assignment_status"] = status
    out["broad_best_candidate"] = best_candidate
    out["broad_best_probability"] = confidence
    out["broad_margin"] = confidence
    out["broad_stability"] = confidence
    out["broad_reason"] = reason
    out["broad_all_required_valid"] = all_required_valid
    out["broad_all_required_negative"] = all_required_negative
    out["broad_authoritative_markers"] = authoritative_marker
    out["broad_anchor_classes"] = np.where(status == "anchor", label, "")
    out["ranked_broad_probabilities"] = empty
    out["broad_type"] = label
    out["best_broad_type"] = best_candidate
    out["subtype_label"] = empty
    out["subtype_assignment_status"] = subtype_status
    out["subtype_best_candidate"] = empty
    out["subtype_best_probability"] = np.nan
    out["subtype_margin"] = np.nan
    out["subtype_stability"] = np.nan
    out["subtype_reason"] = subtype_reason
    out["ranked_subtype_probabilities"] = empty
    out["specific_type"] = label
    out["best_specific_type"] = np.where(best_candidate == "", label, best_candidate)
    out["cell_type"] = label
    out["assignment_status"] = status
    out["specific_assignment_status"] = status
    out["confidence"] = confidence
    out["ranked_specific_probabilities"] = empty
    out["ranked_probabilities"] = empty
    out["assignment_reason"] = reason
    out["blocked_by"] = blocked_by
    out["ordered_residual_gate"] = decisive_gate
    out["qc_analysis_eligible"] = eligible
    out["typing_rules_version"] = active_policy.rules_version
    out["typing_rules_fingerprint"] = active_policy.fingerprint
    out["marker_registry_fingerprint"] = registry.fingerprint
    out["typing_mode"] = TYPING_MODE
    out["typing_method_version"] = METHOD_VERSION
    out["typing_model_bundle_content_id"] = ""
    out["typing_model_bundle_version"] = ""
    out["typing_model_method_version"] = ""
    out["typing_feature_schema_version"] = "restore-tristate-v1"
    out["typing_score_semantics"] = SCORE_SEMANTICS
    out["typing_broad_score_semantics"] = SCORE_SEMANTICS
    out["typing_subtype_score_semantics"] = ""
    out["typing_subtype_model_source"] = "not_applicable"
    out["typing_broad_bootstrap_replicates"] = 0
    out["typing_broad_valid_bootstrap_replicates"] = 0
    out["typing_subtype_bootstrap_replicates"] = 0
    out["typing_subtype_valid_bootstrap_replicates"] = 0
    return out


__all__ = [
    "DEFAULT_RULES_PATH",
    "METHOD_VERSION",
    "OrderedResidualPolicy",
    "SCORE_SEMANTICS",
    "TYPING_MODE",
    "load_ordered_residual_policy",
    "type_ordered_residual_cells",
]
