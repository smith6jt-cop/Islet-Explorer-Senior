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

import numpy as np

from .lymph_node import RuleSource, RuleSourceError

RuleState = Literal["matched", "not_matched", "unavailable"]


@dataclass(frozen=True)
class RuleEvaluation:
    states: dict[str, RuleState]
    missing_markers: dict[str, tuple[str, ...]]


@dataclass(frozen=True)
class BatchRuleEvaluation:
    """One tri-state rule array per phenotype, in input cell order.

    Values are -1 (unavailable), 0 (not matched), and 1 (matched). These are
    rule matches, not a forced final cell class when leaves overlap.
    """

    states: dict[str, np.ndarray]
    unavailable_counts: dict[str, int]


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


def evaluate_rules_batch(source: RuleSource,
                         marker_calls: Mapping[str, np.ndarray]) -> BatchRuleEvaluation:
    """Run the same hierarchy on full cell columns without a Python loop per cell.

    An absent assay marker, or a -1 call in one cell, makes every rule that
    inherits that marker unavailable for that cell. Conflicting canonical and
    alias calls are refused. The caller retains cell identities separately.
    """
    if not marker_calls:
        raise RuleSourceError("batch rule evaluation needs at least one marker column")
    arrays = {marker: np.asarray(value) for marker, value in marker_calls.items()}
    if any(array.ndim != 1 for array in arrays.values()):
        raise RuleSourceError("marker columns must be one-dimensional")
    lengths = {len(array) for array in arrays.values()}
    if len(lengths) != 1:
        raise RuleSourceError("marker columns differ in cell count")
    n = lengths.pop()
    calls: dict[str, np.ndarray] = {}
    alias_to_marker = {alias: marker for marker, alias in source.panel_aliases.items()}
    for marker, array in arrays.items():
        if array.shape != (n,) or not np.issubdtype(array.dtype, np.integer):
            raise RuleSourceError(f"{marker}: batch calls must be one-dimensional integers")
        if not np.isin(array, (-1, 0, 1)).all():
            raise RuleSourceError(f"{marker}: batch call must be -1, 0 or 1")
        canonical = alias_to_marker.get(marker, marker)
        array = array.astype(np.int8, copy=False)
        if canonical in calls:
            old = calls[canonical]
            both = (old != -1) & (array != -1)
            if np.any(both & (old != array)):
                raise RuleSourceError(f"{canonical}: alias and canonical marker calls disagree")
            calls[canonical] = np.where(old == -1, array, old)
        else:
            calls[canonical] = array

    absent = np.full(n, -1, dtype=np.int8)
    rules = source.by_phenotype
    states: dict[str, np.ndarray] = {}

    def state(phenotype: str) -> np.ndarray:
        if phenotype in states:
            return states[phenotype]
        rule = rules[phenotype]
        inherited = source.inherited_markers(phenotype)
        unavailable = np.zeros(n, dtype=bool)
        for marker in inherited:
            unavailable |= calls.get(marker, absent) == -1
        parent = state(rule.parent) if rule.parent in rules else np.ones(n, dtype=np.int8)
        unavailable |= parent == -1
        matched = parent == 1
        any_positive = []
        for marker, token in rule.constraints:
            value = calls.get(marker, absent)
            if token in {"pos", "allpos"}:
                matched &= value == 1
            elif token == "allneg":
                matched &= value == 0
            elif token == "anypos":
                any_positive.append(value == 1)
        if any_positive:
            matched &= np.logical_or.reduce(any_positive)
        result = np.where(unavailable, -1, matched.astype(np.int8))
        states[phenotype] = result
        return result

    for phenotype in rules:
        state(phenotype)
    return BatchRuleEvaluation(states, {name: int((value == -1).sum())
                                        for name, value in states.items()})
