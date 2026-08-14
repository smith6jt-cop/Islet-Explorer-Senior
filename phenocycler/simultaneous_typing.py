"""Simultaneous broad typing from RESTORE tri-state marker evidence.

Every broad gate is evaluated against the same evidence and no gate has
precedence over another.  This is the difference that matters: the ordered
residual hierarchy it replaces awarded a cell supported by several gates to
whichever gate was listed first, and on this cohort that single choice moved the
Immune fraction from 20.50% (Immune first) to 3.31% (Immune last) over the same
25,150,537 cells, because only 16.2% of Immune-gate-positive cells are positive
for no other gate.  A number that swings six-fold on a list ordering is set by
the ordering, not by the evidence.

Dispositions, per cell:

``anchor``
    Exactly one gate is supported.  That gate is the broad label.
``Ambiguous``
    More than one gate is supported.  The competing gates are recorded and no
    winner is chosen.  The margin between the two best-supported gates is
    reported so a later, separately versioned rule can decide these cells on
    evidence; this module does not invent that threshold.
``Other``
    No gate is supported and every gate marker is explicitly negative.
``Unavailable``
    The cell is geometry-QC ineligible, or no gate is supported and some gate
    marker has no measurement -- absence of evidence is not evidence of absence.

Support strength uses ``{marker}__restore_normalized``, which is each marker's
corrected intensity divided by its own donor-local divisor.  That makes markers
comparable across gates: 1.0 is the threshold for every marker by construction.
It is reported, never thresholded here.
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


METHOD_VERSION = "restore-simultaneous-broad-typing-v1"
TYPING_MODE = "restore_simultaneous"
SCORE_SEMANTICS = "restore_normalized_support"
DEFAULT_RULES_PATH = Path(__file__).with_name("broad_gates.json")

POSITIVE = "positive"
NEGATIVE = "negative"
UNAVAILABLE = "unavailable"
_UNAVAILABLE_INPUT_STATES = frozenset(
    {"", "unavailable", "indeterminate", "unknown", "missing"}
)


@dataclass(frozen=True)
class BroadGatePolicy:
    """Content-fingerprinted gate membership, with no precedence between gates."""

    schema_version: int
    rules_version: str
    description: str
    gates: Mapping[str, tuple[str, ...]]
    other_label: str
    ambiguous_label: str
    unavailable_label: str
    fingerprint: str
    source_path: Path

    @property
    def compartments(self) -> tuple[str, ...]:
        """Gate names in a fixed order, for stable column layout only.

        Sorted rather than file order, so that reordering the JSON cannot
        change any assignment even by accident.
        """

        return tuple(sorted(self.gates))

    @property
    def markers(self) -> tuple[str, ...]:
        return tuple(
            marker
            for compartment in self.compartments
            for marker in self.gates[compartment]
        )

    def validate_registry(self, registry: MarkerRegistry) -> None:
        """Require every gate marker to sit in its registered broad compartment.

        The ordered-residual policy this replaces checked only that the gate
        marker set equalled the active marker set.  That let SST -- declared
        ``uses: ["subtype"]`` with no ``broad_compartment`` -- sit in the
        Epithelial broad gate, contradicting the registry contract.
        """

        by_name = {marker.name: marker for marker in registry.markers}
        errors: list[str] = []
        for compartment in self.compartments:
            for marker in self.gates[compartment]:
                spec = by_name.get(marker)
                if spec is None:
                    errors.append(f"{marker}: not in the marker registry")
                    continue
                if spec.calibration_status != "active":
                    errors.append(
                        f"{marker}: gated as {compartment} but its calibration "
                        f"status is {spec.calibration_status!r}, not 'active'"
                    )
                    continue
                if "broad" not in tuple(spec.uses):
                    errors.append(
                        f"{marker}: gated as {compartment} but the registry "
                        f"declares uses={list(spec.uses)}, which excludes "
                        "'broad'"
                    )
                    continue
                declared = getattr(spec, "broad_compartment", None)
                if declared != compartment:
                    errors.append(
                        f"{marker}: gated as {compartment} but the registry "
                        f"declares broad_compartment={declared!r}"
                    )
        gated = set(self.markers)
        for spec in registry.active_markers:
            if "broad" in tuple(spec.uses) and spec.name not in gated:
                errors.append(
                    f"{spec.name}: registry declares it broad "
                    f"{spec.broad_compartment} evidence but no gate contains it"
                )
        if errors:
            raise ContractError(
                "broad gate membership contradicts the marker registry:\n- "
                + "\n- ".join(sorted(errors))
            )

    @classmethod
    def read_json(cls, path: Path | str = DEFAULT_RULES_PATH) -> "BroadGatePolicy":
        source = Path(path).expanduser().resolve()
        try:
            raw = json.loads(source.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ContractError(
                f"cannot read broad gate rules from {source}: {exc}"
            ) from exc
        if not isinstance(raw, dict):
            raise ContractError("broad gate rules must be a JSON object")
        required = {
            "schema_version",
            "rules_version",
            "description",
            "gates",
            "other_label",
            "ambiguous_label",
            "unavailable_label",
        }
        missing_keys = sorted(required - set(raw))
        unknown_keys = sorted(set(raw) - required)
        if missing_keys or unknown_keys:
            raise ContractError(
                "broad gate rule keys differ from the contract "
                f"(missing={missing_keys}, unknown={unknown_keys})"
            )
        if int(raw["schema_version"]) != 1:
            raise ContractError("broad gate schema_version must be 1")
        raw_gates = raw["gates"]
        if not isinstance(raw_gates, dict) or not raw_gates:
            raise ContractError("broad gate rules must declare a nonempty gates object")
        gates: dict[str, tuple[str, ...]] = {}
        for compartment, values in raw_gates.items():
            name = str(compartment).strip()
            if not name:
                raise ContractError("a broad gate name must be nonempty")
            members = tuple(str(value).strip() for value in values)
            if not members or any(not value for value in members):
                raise ContractError(f"broad gate {name!r} must be nonempty")
            if len(members) != len(set(members)):
                raise ContractError(f"broad gate {name!r} has duplicate markers")
            gates[name] = members
        markers = [marker for values in gates.values() for marker in values]
        if len(markers) != len(set(markers)):
            raise ContractError("a marker may occur in only one broad gate")
        labels = {
            key: str(raw[key]).strip()
            for key in ("other_label", "ambiguous_label", "unavailable_label")
        }
        if any(not value for value in labels.values()):
            raise ContractError("terminal broad labels must be nonempty")
        if len(set(labels.values())) != 3:
            raise ContractError("terminal broad labels must be distinct")
        if set(labels.values()) & set(gates):
            raise ContractError("terminal broad labels cannot also be gates")
        canonical = {
            "schema_version": 1,
            "rules_version": str(raw["rules_version"]).strip(),
            "description": str(raw["description"]).strip(),
            "gates": {key: gates[key] for key in sorted(gates)},
            **labels,
        }
        if not canonical["rules_version"]:
            raise ContractError("broad gate rules_version must be nonempty")
        return cls(
            **canonical,
            fingerprint=sha256_json(canonical),
            source_path=source,
        )


def load_broad_gate_policy(
    path: Path | str = DEFAULT_RULES_PATH,
    *,
    registry: MarkerRegistry | None = None,
) -> BroadGatePolicy:
    policy = BroadGatePolicy.read_json(path)
    policy.validate_registry(load_registry() if registry is None else registry)
    return policy


def _canonical_state(series: pd.Series, *, marker: str) -> np.ndarray:
    raw = series.astype("string").fillna("").str.strip().str.lower()
    unknown = sorted(set(raw) - {POSITIVE, NEGATIVE, *_UNAVAILABLE_INPUT_STATES})
    if unknown:
        raise ContractError(
            f"RESTORE marker {marker!r} has unknown cell states: {unknown}"
        )
    state = raw.to_numpy(dtype=object)
    state[np.isin(state, tuple(_UNAVAILABLE_INPUT_STATES))] = UNAVAILABLE
    return state


def type_simultaneous_cells(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None = None,
    policy: BroadGatePolicy | None = None,
) -> pd.DataFrame:
    """Assign a broad label from all gates at once, preserving genuine ambiguity."""

    registry = load_registry() if marker_registry is None else marker_registry
    active = (
        load_broad_gate_policy(registry=registry) if policy is None else policy
    )
    active.validate_registry(registry)

    required = {"donor_id", "object_id", "qc_analysis_eligible"}
    required.update(f"{marker}__state" for marker in active.markers)
    missing = sorted(required - set(evidence_wide))
    if missing:
        raise ContractError(
            f"simultaneous typing evidence is missing columns: {missing}"
        )
    keys = evidence_wide.loc[:, ["donor_id", "object_id"]].copy()
    if keys.isna().any().any():
        raise ContractError("simultaneous typing keys cannot be null")
    keys["donor_id"] = keys["donor_id"].astype(str)
    keys["object_id"] = keys["object_id"].astype(str)
    if keys.duplicated(["donor_id", "object_id"]).any():
        raise ContractError("simultaneous typing contains duplicate cells")
    eligibility = evidence_wide["qc_analysis_eligible"]
    if eligibility.isna().any() or not eligibility.isin([True, False, 0, 1]).all():
        raise ContractError("qc_analysis_eligible must contain explicit booleans")
    eligible = eligibility.astype(bool).to_numpy()
    n = len(evidence_wide)

    states = {
        marker: _canonical_state(evidence_wide[f"{marker}__state"], marker=marker)
        for marker in active.markers
    }
    compartments = active.compartments

    supported = np.zeros((n, len(compartments)), dtype=bool)
    incomplete = np.zeros((n, len(compartments)), dtype=bool)
    support = np.zeros((n, len(compartments)), dtype=float)
    first_marker = np.full((n, len(compartments)), "", dtype=object)

    for column, compartment in enumerate(compartments):
        for marker in active.gates[compartment]:
            state = states[marker]
            marker_positive = state == POSITIVE
            supported[:, column] |= marker_positive
            incomplete[:, column] |= state == UNAVAILABLE
            first_marker[(first_marker[:, column] == "") & marker_positive, column] = (
                marker
            )
            column_name = f"{marker}__restore_normalized"
            if column_name in evidence_wide:
                value = pd.to_numeric(
                    evidence_wide[column_name], errors="coerce"
                ).to_numpy(dtype=float)
                value = np.where(marker_positive & np.isfinite(value), value, 0.0)
                support[:, column] = np.maximum(support[:, column], value)
    # A positive sibling proves gate membership; an unmeasured sibling then
    # cannot make the gate uncertain.
    incomplete &= ~supported

    n_supported = supported.sum(axis=1)
    any_incomplete = incomplete.any(axis=1)

    ordered = np.argsort(-support, axis=1, kind="stable")
    best_column = ordered[:, 0]
    rows = np.arange(n)
    best_support = support[rows, best_column]
    second_support = (
        support[rows, ordered[:, 1]] if len(compartments) > 1 else np.zeros(n)
    )
    compartment_array = np.asarray(compartments, dtype=object)

    label = np.full(n, active.other_label, dtype=object)
    status = np.full(n, "other", dtype=object)
    reason = np.full(n, "all_negative", dtype=object)
    authoritative = np.full(n, "", dtype=object)
    candidates = np.full(n, "", dtype=object)
    margin = np.full(n, np.nan, dtype=float)

    single = eligible & (n_supported == 1)
    multiple = eligible & (n_supported > 1)
    blocked = eligible & (n_supported == 0) & any_incomplete

    single_column = supported[single].argmax(axis=1) if single.any() else np.array([], int)
    label[single] = compartment_array[single_column]
    status[single] = "anchor"
    reason[single] = "single_supported_gate"
    authoritative[single] = first_marker[single, single_column]
    candidates[single] = compartment_array[single_column]

    label[multiple] = active.ambiguous_label
    status[multiple] = "ambiguous"
    reason[multiple] = "multiple_supported_gates"
    if multiple.any():
        authoritative[multiple] = first_marker[multiple, best_column[multiple]]
        candidates[multiple] = [
            "|".join(sorted(compartment_array[row]))
            for row in supported[multiple]
        ]
        margin[multiple] = best_support[multiple] - second_support[multiple]

    label[blocked] = active.unavailable_label
    status[blocked] = "unavailable"
    reason[blocked] = "unmeasured_gate_marker"
    if blocked.any():
        candidates[blocked] = [
            "|".join(sorted(compartment_array[row])) for row in incomplete[blocked]
        ]

    label[~eligible] = active.unavailable_label
    status[~eligible] = "unavailable"
    reason[~eligible] = "geometry_qc_ineligible"
    authoritative[~eligible] = ""
    candidates[~eligible] = ""
    margin[~eligible] = np.nan

    all_required_valid = eligible.copy()
    all_required_negative = eligible.copy()
    for state in states.values():
        all_required_valid &= state != UNAVAILABLE
        all_required_negative &= state == NEGATIVE

    empty = np.full(n, "", dtype=object)
    nan = np.full(n, np.nan, dtype=float)

    out = keys.copy()
    out["broad_label"] = label
    out["broad_assignment_status"] = status
    out["broad_best_candidate"] = np.where(
        status == "anchor", label, np.where(status == "other", "", candidates)
    )
    # Categorical gate outcomes are not probabilities, margins or bootstrap
    # stabilities. The ordered lane stamped 1.0 into all three while also
    # recording zero bootstrap replicates; these stay NaN because they were
    # never estimated.
    out["broad_best_probability"] = nan
    out["broad_margin"] = nan
    out["broad_stability"] = nan
    out["broad_reason"] = reason
    out["broad_all_required_valid"] = all_required_valid
    out["broad_all_required_negative"] = all_required_negative
    out["broad_authoritative_markers"] = authoritative
    out["broad_anchor_classes"] = np.where(status == "anchor", label, "")
    out["broad_supported_gates"] = candidates
    out["broad_supported_gate_n"] = np.where(eligible, n_supported, 0)
    out["broad_support"] = np.where(eligible & (n_supported > 0), best_support, np.nan)
    out["broad_support_margin"] = margin
    out["ranked_broad_probabilities"] = np.full(n, "[]", dtype=object)
    out["broad_type"] = label
    out["best_broad_type"] = out["broad_best_candidate"]
    out["subtype_label"] = empty
    out["subtype_assignment_status"] = np.where(
        eligible, "not_applicable", "unavailable"
    )
    out["subtype_best_candidate"] = empty
    out["subtype_best_probability"] = nan
    out["subtype_margin"] = nan
    out["subtype_stability"] = nan
    out["subtype_reason"] = np.where(
        eligible, "broad_only_simultaneous", "geometry_qc_ineligible"
    )
    out["ranked_subtype_probabilities"] = np.full(n, "[]", dtype=object)
    out["specific_type"] = label
    out["best_specific_type"] = label
    out["cell_type"] = label
    out["assignment_status"] = status
    out["specific_assignment_status"] = status
    out["confidence"] = nan
    out["ranked_specific_probabilities"] = np.full(n, "[]", dtype=object)
    out["ranked_probabilities"] = np.full(n, "[]", dtype=object)
    out["assignment_reason"] = reason
    out["blocked_by"] = np.where(status == "unavailable", candidates, "")
    out["qc_analysis_eligible"] = eligible
    out["typing_rules_version"] = active.rules_version
    out["typing_rules_fingerprint"] = active.fingerprint
    out["marker_registry_fingerprint"] = registry.fingerprint
    out["typing_mode"] = TYPING_MODE
    out["typing_method_version"] = METHOD_VERSION
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
    "BroadGatePolicy",
    "SCORE_SEMANTICS",
    "TYPING_MODE",
    "load_broad_gate_policy",
    "type_simultaneous_cells",
]
