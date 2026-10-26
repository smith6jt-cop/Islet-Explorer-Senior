"""Evaluate a recorded protein hierarchy against explicit marker calls.

This is the shared Boolean step for the pLN and spleen rule sources. The caller
must establish per-cell marker positivity in the appropriate donor and assay
units first; this module neither chooses a cutoff nor assigns a final cell type.
Unavailable markers remain unavailable, including in inherited parent rules.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

from .lymph_node import RuleSource, RuleSourceError

RuleState = Literal["matched", "not_matched", "unavailable"]


@dataclass(frozen=True)
class RuleEvaluation:
    states: dict[str, RuleState]
    missing_markers: dict[str, tuple[str, ...]]


def _marker_calls(source: RuleSource, calls: Mapping[str, bool | None]) -> dict[str, bool]:
    """Resolve recorded aliases without treating a missing call as negative."""
    result: dict[str, bool] = {}
    for marker, value in calls.items():
        if value is None:
            continue
        if not isinstance(value, bool):
            raise RuleSourceError(f"{marker}: marker call must be true, false or unavailable")
        canonical = next((name for name, alias in source.panel_aliases.items()
                          if alias == marker), marker)
        if canonical in result and result[canonical] != value:
            raise RuleSourceError(f"{canonical}: alias and canonical marker calls disagree")
        result[canonical] = value
    return result


def evaluate_rules(source: RuleSource,
                   marker_calls: Mapping[str, bool | None]) -> RuleEvaluation:
    """Report each phenotype's rule result and exact missing inherited markers.

    ``pos`` and ``allpos`` require their markers to be positive; ``allneg``
    requires negatives; the ``anypos`` markers form one OR group. A child's
    rule can match only when its parent matches. Any absent marker in a rule's
    ancestry makes that phenotype unavailable, even if other clauses match.
    """
    calls = _marker_calls(source, marker_calls)
    missing = source.missing_by_phenotype(calls)
    rules = source.by_phenotype
    states: dict[str, RuleState] = {}

    def state(phenotype: str) -> RuleState:
        if phenotype in states:
            return states[phenotype]
        rule = rules[phenotype]
        if missing[phenotype]:
            result: RuleState = "unavailable"
        elif rule.parent in rules and state(rule.parent) != "matched":
            result = states[rule.parent]
        else:
            positives = [calls[marker] for marker, token in rule.constraints
                         if token in {"pos", "allpos"}]
            negatives = [calls[marker] for marker, token in rule.constraints
                         if token == "allneg"]
            any_positive = [calls[marker] for marker, token in rule.constraints
                            if token == "anypos"]
            ok = (all(positives) and not any(negatives)
                  and (not any_positive or any(any_positive)))
            result = "matched" if ok else "not_matched"
        states[phenotype] = result
        return result

    for phenotype in rules:
        state(phenotype)
    return RuleEvaluation(states, missing)
