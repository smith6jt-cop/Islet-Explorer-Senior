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
CANDIDATE_STATES = ("unique_candidate", "ambiguous", "unavailable", "no_match")


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


@dataclass(frozen=True)
class BatchRuleCandidates:
    """Deepest uniquely matched lineage rule per cell, for review only.

    The selected 6551 manifest names lineage and state roots. A parent is superseded only by one
    of its matched descendants. Two matched lineage branches stay ambiguous;
    state annotations never compete with a lineage candidate. Missing lineage
    evidence is counted separately and no candidate is called verified.
    """

    phenotype: np.ndarray
    status: np.ndarray
    maximal_matches: np.ndarray
    unavailable_lineage_rules: np.ndarray


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


def resolve_unique_matched_rules(source: RuleSource,
                                 evaluation: BatchRuleEvaluation) -> BatchRuleCandidates:
    """Keep one deepest confirmed lineage, or report why no unique one exists.

    The sheet's state branches remain in the rule evidence, but are
    excluded from cell-identity competition. This does not choose among
    overlapping lineage leaves or turn a candidate into a served final class.
    When no lineage matches, availability is decided at the broad lineage
    roots; a missing subtype marker cannot erase a measured negative root.
    """
    rules = source.by_phenotype
    if set(evaluation.states) != set(rules):
        raise RuleSourceError("candidate resolution needs every recorded rule")
    arrays = {name: np.asarray(evaluation.states[name]) for name in rules}
    lengths = {len(value) for value in arrays.values() if value.ndim == 1}
    if (len(lengths) != 1 or any(value.ndim != 1 or not np.issubdtype(value.dtype, np.integer)
                                 or not np.isin(value, (-1, 0, 1)).all()
                                 for value in arrays.values())):
        raise RuleSourceError("candidate resolution needs equal-length tri-state rule columns")
    n = lengths.pop()
    children: dict[str, list[str]] = {name: [] for name in rules}
    for name, rule in rules.items():
        if rule.parent in rules:
            parent = arrays[rule.parent]
            if np.any((arrays[name] == 1) & (parent != 1)):
                raise RuleSourceError(f"{name}: matched child has no matched parent")
            children[rule.parent].append(name)

    declared_roots = set(source.candidate_lineage_roots)
    lineage_roots = {name for name, rule in rules.items()
                     if rule.parent == source.lineage_root_parent}
    if declared_roots and not declared_roots <= lineage_roots:
        raise RuleSourceError("candidate lineage root is absent from the hierarchy")

    def lineage_rule(name: str) -> bool:
        current, seen = name, set()
        while rules[current].parent in rules:
            if current in seen:
                raise RuleSourceError(f"cycle in rule parents at {current}")
            seen.add(current)
            current = rules[current].parent
        root_parent = rules[current].parent
        if root_parent not in {source.lineage_root_parent, source.state_root_parent}:
            raise RuleSourceError(f"{name}: unknown hierarchy root {root_parent}")
        return (root_parent == source.lineage_root_parent
                and (not source.candidate_lineage_roots
                     or current in source.candidate_lineage_roots))

    names = tuple(name for name in rules if lineage_rule(name))
    if not names:
        raise RuleSourceError("candidate resolution needs lineage rules")
    candidate_roots = (source.candidate_lineage_roots or
                       tuple(name for name in names
                             if rules[name].parent == source.lineage_root_parent))
    winner = np.full(n, -1, dtype=np.int16)
    maximal = np.zeros(n, dtype=np.uint16)
    unavailable = np.zeros(n, dtype=np.uint16)
    for index, name in enumerate(names):
        state = arrays[name]
        unavailable += (state == -1).astype(np.uint16)
        deepest = state == 1
        for child in children[name]:
            deepest &= arrays[child] != 1
        winner[deepest] = index
        maximal += deepest.astype(np.uint16)
    phenotype = np.asarray(names, dtype=object)[np.maximum(winner, 0)]
    phenotype[maximal != 1] = ""
    status = np.full(n, "no_match", dtype=object)
    root_unavailable = np.logical_or.reduce([arrays[name] == -1 for name in candidate_roots])
    status[(maximal == 0) & root_unavailable] = "unavailable"
    status[maximal > 1] = "ambiguous"
    status[maximal == 1] = "unique_candidate"
    return BatchRuleCandidates(phenotype, status, maximal, unavailable)
