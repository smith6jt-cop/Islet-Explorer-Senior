"""hierarchy.py -- name cells by a tissue's recorded protein hierarchy, from per-image positivity.

Maintainer, 2026-11-01 (spleen): cells are named by the tissue's selected hierarchy -- for spleen the 6551 sheet with
its spleen adaptation (``tissues.spleen.load_rule_source``), evaluated by ``tissues.evaluate`` -- and not by the CODEX
identity step (``identity.py`` and a labels file). The per-image method supplies the marker calls and the tiebreak:

- **Calls.** A cell's call for a marker is its value on the marker's positivity measurement against that section's
  level (``levels.choose``): positive above it, negative at or below it, unavailable where the measurement is missing or
  the marker has no level (a method that failed is not a negative).
- **Structure-only branches** (spleen: Endothelial) name a cell only inside a structure object of their class. Without
  an object layer a match there is unavailable, never negative (true AND unknown); with one, a match outside an object
  is not that class.
- **Ties.** A cell matching several lineage rules that are not one another's ancestors is settled by the CODEX canary
  rule (``identity.placed``): the groups are the branches diverging at the matches' lowest common ancestor, and each
  group's detecting markers are the positive-clause markers on its rules that have a typing role and a canary level.
  The canary basis comes first: the canary depth margin (CANARY_MARGIN, CANARY_MIN_DEPTH), then at equal depth the
  group that is both brighter (fold over level, on the differentiation reading) and covers at least
  CANARY_TIE_INTERIOR_LEAD more of the cell's interior. A tie only plain positivity can see is settled by the same
  rule on plain positivity. A tie that stays unresolved at either basis is not settled ("otherwise the cell stays
  ambiguous", maintainer 2026-11-01). CODEX itself places the top plain group below the canary level, so this route is
  stricter there. The winner's own subtype tie, if any, is settled the same way, one level down.

Classes (one class reached in two steps; only the final step is served):
- ``refined_class``: the deepest uniquely matched or tie-settled lineage rule. With an unresolved subtype tie it is
  the lineage the tied rules share (status subtype_ambiguous); with an unresolved tie between lineage roots it is empty
  (lineage_ambiguous, or doublet_suspected for three or more canary groups);
- ``cell_class``: the broad step. That is the refined class's root, or, where that root's children have subtypes of
  their own (spleen: Immune), the child on the way to the refined class (T_cell, B_cell, Myeloid ...).

The sheet's state roots (parent ``any``) never compete for identity; the ones a cell matches are listed in ``states``.
"""
from __future__ import annotations

import re
from dataclasses import replace
from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd

from ..tissues.evaluate import BatchRuleEvaluation, evaluate_rules_batch, resolve_unique_matched_rules
from ..tissues.lymph_node import RuleSource
from . import identity as ID
from . import registry as RG

STATUSES = ('unique', 'tie_settled', 'subtype_ambiguous', 'lineage_ambiguous', 'doublet_suspected', 'unavailable',
            'no_match', 'no_supported_subtype')


class HierarchyError(ValueError):
    """The calls, the levels or the rule source cannot be read together as declared."""


# -- names ------------------------------------------------------------------------------------------------------------
def _key(name: str) -> str:
    return re.sub(r'[^0-9a-z]', '', name.lower())


def sheet_names(source: RuleSource, markers: Iterable[str]) -> dict[str, str]:
    """Registry marker -> the sheet's marker, matched by name ignoring case, spaces, hyphens and underscores
    (Pan_Cytokeratin -> Pan-Cytokeratin, PD_L1 -> PD-L1, HLA_DR -> HLA-DR). Markers the sheet does not use are left
    out. Two sheet markers, or two registry markers, that reduce to one name are refused."""
    sheet: dict[str, str] = {}
    for m in sorted({m for rule in source.rules for m in rule.markers}):
        if _key(m) in sheet:
            raise HierarchyError(f'sheet markers {sheet[_key(m)]!r} and {m!r} cannot be told apart by name')
        sheet[_key(m)] = m
    out: dict[str, str] = {}
    for m in markers:
        if _key(m) in sheet:
            if sheet[_key(m)] in out.values():
                raise HierarchyError(f'two registry markers name the sheet marker {sheet[_key(m)]!r}')
            out[m] = sheet[_key(m)]
    return out


# -- calls ------------------------------------------------------------------------------------------------------------
def marker_calls(d: pd.DataFrame, reg: RG.Registry, levels: pd.DataFrame, markers: Iterable[str]) -> dict[str, np.ndarray]:
    """Tri-state calls (1 positive, 0 negative, -1 unavailable) of `markers` for the cells of `d`, from one section's
    levels (``levels`` indexed by marker: measurement, level_log). Positive is log1p(value) above the level, as
    ``identity`` and ``canary`` read it; a missing value or level is unavailable, never negative."""
    out = {}
    for m in markers:
        if m not in levels.index:
            raise HierarchyError(f'{m}: no level row for this section')
        v = d[reg.column(m, levels.loc[m, 'measurement'])].to_numpy(float)
        level = float(levels.loc[m, 'level_log'])
        with np.errstate(invalid='ignore'):
            call = (np.log1p(np.clip(v, 0, None)) > level).astype(np.int8)
        call[~np.isfinite(v)] = -1
        if not np.isfinite(level):
            call[:] = -1
        out[m] = call
    return out


# -- the hierarchy ------------------------------------------------------------------------------------------------------
def _chain(rules: Mapping, name: str) -> list[str]:
    """`name` and its ancestors, deepest first."""
    out = [name]
    while rules[out[-1]].parent in rules:
        out.append(rules[out[-1]].parent)
    return out


def evaluate(source: RuleSource, calls: Mapping[str, np.ndarray], structure_only: Iterable[str] = (),
             in_structure: Mapping[str, np.ndarray] | None = None) -> BatchRuleEvaluation:
    """The rule states of every cell (``evaluate_rules_batch``) with the structure-only branches applied: a match in
    such a branch is unavailable where no object layer of its class is given, and not matched outside its objects
    where one is (``in_structure[root]``: the cell lies in an object of that class)."""
    rules = source.by_phenotype
    branches = set(structure_only)
    if not branches <= set(rules):
        raise HierarchyError(f'structure-only branches absent from the rules: {sorted(branches - set(rules))}')
    ev = evaluate_rules_batch(source, calls)
    states = dict(ev.states)
    for name in rules:
        root = _chain(rules, name)[-1]
        if root not in branches:
            continue
        s = states[name]
        inside = (in_structure or {}).get(root)
        if inside is None:
            states[name] = np.where(s == 1, -1, s).astype(np.int8)
        else:
            states[name] = np.where((s == 1) & ~np.asarray(inside, bool), 0, s).astype(np.int8)
    return BatchRuleEvaluation(states, {k: int((v == -1).sum()) for k, v in states.items()})



def evaluate_profiles(source: RuleSource, calls: Mapping[str, np.ndarray],
                      structure_only: Iterable[str] = (),
                      in_structure: Mapping[str, np.ndarray] | None = None, *,
                      defer_exclusions: Mapping[str, Iterable[str]] | None = None) -> BatchRuleEvaluation:
    """Compare affirmative profiles before resolving recorded lineage exclusions.

    The tissue explicitly names exclusion markers whose positive signal is a
    competing profile, rather than a veto before that comparison. Original
    calls and rules are retained. Missing exclusion measurements still make
    their descendants unavailable; structure-only admission is unchanged.
    With no tissue policy this is the strict recorded Boolean evaluation.
    """
    if not defer_exclusions:
        return evaluate(source, calls, structure_only, in_structure)
    by_name = source.by_phenotype
    scopes = {name: set(markers) for name, markers in defer_exclusions.items()}
    if not set(scopes) <= set(by_name):
        raise HierarchyError('profile exclusions name an unknown hierarchy branch')
    known = {marker for rule in source.rules for marker in rule.markers}
    if any(not markers <= known for markers in scopes.values()):
        raise HierarchyError('profile exclusions name an unknown marker')
    profiled = []
    for rule in source.rules:
        ancestors = set(_chain(by_name, rule.phenotype))
        deferred = set().union(*(markers for branch, markers in scopes.items() if branch in ancestors))
        profiled.append(replace(rule, constraints=tuple(
            (marker, token) for marker, token in rule.constraints
            if not (token == 'allneg' and marker in deferred))))
    strict = evaluate(source, calls, structure_only, in_structure)
    candidates = evaluate(replace(source, rules=tuple(profiled)), calls, structure_only, in_structure)
    states = {name: np.where(strict.states[name] == -1, -1, value).astype(np.int8)
              for name, value in candidates.states.items()}
    return BatchRuleEvaluation(states, {name: int((value == -1).sum()) for name, value in states.items()})


def _lineage_names(source: RuleSource) -> list[str]:
    rules = source.by_phenotype
    roots = set(source.candidate_lineage_roots) or {n for n, r in rules.items() if r.parent == source.lineage_root_parent}
    return [n for n in rules if _chain(rules, n)[-1] in roots and rules[_chain(rules, n)[-1]].parent ==
            source.lineage_root_parent]


def broad_class(source: RuleSource, name: str) -> str:
    """The broad step of a refined class: its root, or the root's child on the way to it where the root's children
    have subtypes of their own (a compartment root such as Immune)."""
    if not name:
        return ''
    rules = source.by_phenotype
    chain = _chain(rules, name)
    root = chain[-1]
    children = [n for n, r in rules.items() if r.parent == root]
    compartment = any(r.parent in children for r in rules.values())
    return chain[-2] if compartment and len(chain) >= 2 else root


def name_cells(source: RuleSource, evaluation: BatchRuleEvaluation, evidence: Mapping, margin: float,
               min_depth: float, lead: float, *, require_specific: Iterable[str] = (),
               inherit_depth: bool = False) -> pd.DataFrame:
    """One row per cell: cell_class, refined_class, status, tie_basis, n_matched, states.

    ``evidence``: ``columns`` (sheet marker names, in matrix order), ``P`` and ``C`` (plain and canary positive, cells x
    markers), ``F`` (fold over level on the differentiation reading), ``COV`` (interior coverage), ``interior_valid``
    (cells with interior pixels), and ``typing`` (the sheet markers with a typing role and a canary level; only these
    detect a group)."""
    rules = source.by_phenotype
    required = set(require_specific)
    if not required <= set(rules):
        raise HierarchyError('specific-profile policy names an unknown rule')
    base = resolve_unique_matched_rules(source, evaluation)
    n = len(base.status)
    names = _lineage_names(source)
    states = {k: np.asarray(v) for k, v in evaluation.states.items()}
    children: dict[str, list[str]] = {k: [] for k in rules}
    for k, r in rules.items():
        if r.parent in rules:
            children[r.parent].append(k)
    deepest = np.column_stack([(states[k] == 1) & ~np.logical_or.reduce(
        [states[c] == 1 for c in children[k]] or [np.zeros(n, bool)]) for k in names])
    refined = np.where(base.status == 'unique_candidate', base.phenotype, '').astype(object)
    status = np.select([base.status == 'unique_candidate', base.status == 'unavailable', base.status == 'no_match'],
                       ['unique', 'unavailable', 'no_match'], '').astype(object)
    basis = np.full(n, '', dtype=object)

    col = {m: i for i, m in enumerate(evidence['columns'])}
    typing = set(evidence['typing'])
    P, C, F, COV = (np.asarray(evidence[k]) for k in ('P', 'C', 'F', 'COV'))
    valid = np.asarray(evidence['interior_valid'], bool)

    def detecting(group: str, members: list[str]) -> list[int]:
        """Positive-clause markers on the rules from the group's root down to its matched members."""
        out: dict[str, None] = {}
        for m in members:
            chain = _chain(rules, m)
            for r in chain[:chain.index(group) + 1]:
                for marker, token in rules[r].constraints:
                    if token in ('pos', 'allpos', 'anypos') and marker in typing and marker in col:
                        out[marker] = None
        return [col[x] for x in out]

    queue = []
    amb = np.flatnonzero(base.status == 'ambiguous')
    if len(amb):
        sig = pd.Series([tuple(np.asarray(names, dtype=object)[deepest[i]]) for i in amb], index=amb)
        queue = [(np.asarray(ix, int), list(key)) for key, ix in sig.groupby(sig).groups.items()]
    while queue:
        rows, current = queue.pop()
        chains = [_chain(rules, m) for m in current]
        common = [a for a in chains[0] if all(a in c for c in chains[1:])]
        lca = common[0] if common else ''
        groups: dict[str, list[str]] = {}
        for m, c in zip(current, chains):
            groups.setdefault(c[c.index(lca) - 1] if lca else c[-1], []).append(m)
        gnames = list(groups)
        det = [detecting(g, groups[g]) for g in gnames]
        # A subtype retains affirmative parent evidence for confidence.
        # OR-group alternatives do not independently support every branch.
        # Brightness and coverage still use only each branch's detecting markers.
        shared = set()
        if inherit_depth and lca:
            for ancestor in _chain(rules, lca):
                for marker, token in rules[ancestor].constraints:
                    if token in ('pos', 'allpos') and marker in typing and marker in col:
                        shared.add(col[marker])
        count = [sorted(set(markers) | shared) for markers in det]
        res = ID.placed(P[rows], C[rows], F[rows], COV[rows], det, count, margin, min_depth, lead, valid[rows])
        st = res['status'].astype(str)
        ok = np.where(res['basis'] == 'canary', np.isin(st, ID.PLACED),
                      (res['n_plain'] == 1) | (st == 'below_canary_resolved')) & (res['group'] >= 0)
        for i in np.unique(res['group'][ok]):
            sub = rows[ok & (res['group'] == i)]
            members = groups[gnames[i]]
            if len(members) == 1:
                refined[sub], status[sub] = members[0], 'tie_settled'
                basis[sub] = res['basis'][ok & (res['group'] == i)]
            else:
                queue.append((sub, members))
        left = rows[~ok]
        if len(left):
            if lca:
                refined[left], status[left] = lca, 'subtype_ambiguous'
            else:
                refined[left] = ''
                status[left] = np.where(st[~ok] == 'doublet_suspected', 'doublet_suspected', 'lineage_ambiguous')
            basis[left] = res['basis'][~ok]

    if (status == '').any():
        raise HierarchyError('a cell was left without a status')
    broad = [broad_class(source, name) for name in refined]
    for name in required:
        unresolved = refined == name
        unique = unresolved & np.isin(status, ('unique', 'tie_settled'))
        missing = np.logical_or.reduce([states[c] == -1 for c in children[name]] or [np.zeros(n, bool)])
        status[unique & missing] = 'unavailable'
        status[unique & ~missing] = 'no_supported_subtype'
        if rules[name].parent == source.lineage_root_parent:
            status[unresolved & (status == 'subtype_ambiguous')] = 'lineage_ambiguous'
        refined[unresolved] = ''
    state_roots = [k for k, r in rules.items() if _chain(rules, k)[-1] in rules
                   and rules[_chain(rules, k)[-1]].parent == source.state_root_parent]
    matched_states = np.full(n, '', dtype=object)
    for k in state_roots:
        hit = states[k] == 1
        matched_states[hit] = np.where(matched_states[hit] == '', k, matched_states[hit] + ';' + k)
    return pd.DataFrame({'cell_class': broad, 'refined_class': refined,
                         'status': status, 'tie_basis': basis, 'n_matched': base.maximal_matches.astype(int),
                         'states': matched_states})
