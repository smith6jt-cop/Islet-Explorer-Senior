"""identity.py -- canary evidence and identity for every cell (final typing plan of 2026-10-28: method Steps 3-4, order of
work Step 6).

Ported from Islet-Explorer scripts/uniform/typing/identity.py at 578767e (``positivity.SOURCE_COMMIT``). The computation
is unchanged. Changed in the move: cells, coverage and an established typing come from a source; the levels, canary and
selection tables are passed in (``tables``); the start class comes from the tissue hook (``start_class``; pancreas: the
endocrine object class, else Tissue, as here); ``classify_image`` returns the table and audit the script wrote.

Inputs, per image: the registry and its labels file (`registry/<panel>__<tissue>.csv`, `..__labels.csv`: marker roles,
compartments, lineages, label rules), every typing marker's positivity level (`typing/levels.csv`, levels.py) and canary
level (`typing/canary_levels.csv`, canary.py), the REDSEA value of each marker's one measurement, and the cell's context.
This module names compartments and lineages, never markers.

1 Start class: the parent object's class, Islet, EndocrineCluster, EndocrineSingle or Tissue (owner 2026-10-28).
2 Canary (Islet-Explorer-Senior Xenium notebook 02 cell 22, adapted as planned): a marker is canary-positive above its
  canary level. Per compartment, depth = the canary-positive markers among its defining markers and the shared markers it
  carries, counted only when a defining marker is canary-positive; groups detected = compartments with a canary-positive
  defining marker; status none / single / resolved (2 groups, top depth >= CANARY_MIN_DEPTH and >= CANARY_MARGIN above
  the runner-up) / lineage_ambiguous (2 otherwise) / doublet_suspected (>= 3). Groups rank by depth x 1000 + the summed
  log fold over level of their canary-positive markers; the margin is taken between compartments (plan: fixes Xenium's
  margin between lineages of one group); CANARY_MARGIN is 1 (owner 2026-10-29: this panel has 1-3 markers per group).
  Two groups of equal depth: the one with higher median calibrated brightness and the required exclusive interior area lead
  is placed, else the conflict stays (owner 2026-10-29). The same rule runs between the lineages of the cell's compartment: lineage-
  specific markers detect a lineage; markers on several lineages, the compartment's pan marker and its shared markers add
  depth to every detected lineage. A cell with no canary-positive defining marker is placed the same way by its positive
  markers (plain level) with status none: agent determination, recorded in the ledger (the plan does not cover a
  defining marker positive below its canary level; the owner's rule admits names from positive markers; Xenium types
  cells before its canary step and uses the canary only for conflicts).
3 Identity by the labels file. Endocrine (owner 2026-10-29): inside the endocrine objects the unique highest scaled hormone
  above its level names the cell, as Senior (subpanel step 5), no margin; an unresolved spatial tie keeps the start class
  (tied_highest_hormones). The established path below is kept for a labels row named_by 'established' (none now): owner
  2026-10-28: "The endocrine typing didn't need to change", "If the singles were already found and verified, don't undo
  that"): its decisions inside the islets, clusters and singles (Beta / Alpha / Delta, CHGA+ hormone-negative, unresolved
  subtype, and its exclusions hormone_unresolved_margin, exocrine_inside_islet_border, endocrine_single_*) are kept as
  they were, and no other cell is named endocrine: a cell positive for an endocrine marker is typed by its other evidence
  and tagged endocrine_marker_positive (endocrine_outside_object outside the objects). Exocrine: an epithelial cell
  without another lineage; inside an endocrine object Excluded (exocrine_inside_islet_border). A lineage named by its broad marker (T, B, Myeloid)
  needs that marker positive; T subtypes by positivity, both positive -> the higher fold over level (owner), neither ->
  the lineage label with t_subtype_status T_subtype_not_shown (interim, owner decision pending); B needs the registry's
  `requires` (CD45 positive, one hematopoietic lineage); Myeloid takes the profile label. A leukocyte without a named
  lineage is the pan label (Hematopoietic CD45+) when the pan marker is positive, else the pan-negative label
  (Hematopoietic CD45-; owner 2026-10-29). No cell is named by a marker (owner 2026-10-28: "Get
  rid of naming by top marker logic"): a cell its compartment and lineage do not name keeps its start class, with the
  status saying why: no_lineage_named (a compartment without a named lineage or pan marker), <lineage>_marker_only (a
  lineage barred where the cell is, e.g. an endocrine marker outside the endocrine objects), shared_marker_only (only
  shared markers positive), no_marker_above_background (no typing marker above its level). No marker below its level
  names a cell; markers without a typing role are never read.
4 Confirmed: compartment and lineage single or resolved, and a marker that confirms the label's lineage (registry
  `confirms`) and was not used to name it is canary-positive; NaN where the lineage has no such marker.

Readings (owner 2026-10-30: "Use membrane for differentiating, but then the whole cell still needs to be used for
positivity"): positivity and canary positivity are read on each marker's positivity measurement (whole cell; Ki67 the
nucleus); the fold over level that ranks and compares groups competing for a cell (groups' scores, the
conflict order, the tie rule's brightness, the structure step's top marker; fold__<marker> in the output) is read on the
marker's differentiation measurement where the registry gives one (typing/levels_differentiation.csv). The tie rule's
Native conflict resolution measures both groups on the same cell-interior pixels.
The median calibrated intensity must favor the group with the qualifying exclusive area lead.
Structure-defined compartments (owner 2026-10-30): a compartment whose labels are all structure rows (Neural; Endothelial,
whose blood-vessel and lymphatic cells come from their objects) is reached only through the structure step. Outside a
structure call its markers place no cell (they still add depth where the registry has another compartment carry them);
the overlap stays a feature (crossed_frac_*). A cell they would have placed is typed by its other evidence, else keeps
its start class with status structure_marker_only; tag <compartment>_marker_outside_structure.

``classify_image`` returns one row per analysis-eligible cell and an audit of counts.
"""
from __future__ import annotations

from collections.abc import Callable
from contextlib import ExitStack

import numpy as np
import pandas as pd

from . import coverage as AC
from .interior_signal import InteriorSignal
from . import registry as RG

STATUSES = ('none', 'single', 'resolved', 'lineage_ambiguous', 'doublet_suspected')
# the established endocrine typing: the run's typing state before the final run (v17's calls, never assembled; its
# endocrine chain unchanged since v15, decisions 98-99.5 % identical with v16's), and the exclusions it made
ESTABLISHED = 'superseded/typing-v17-20261028'
ESTABLISHED_REASONS = ('hormone_unresolved_margin', 'exocrine_inside_islet_border')
ESTABLISHED_REASON_PREFIX = 'endocrine_single_'
MARGIN_REASONS = ('hormone_unresolved_margin', 'endocrine_single_hormone_margin_below_2x')
PLACED = ('single', 'resolved')
COV_PREFIX = 'cov__'        # area_coverage.py
COV_INTERIOR_PREFIX = 'covint__'    # area_coverage.py, the cell's interior
# the cell context's columns for each structure class (cell_context.py): centroid in an object, area fraction inside
STRUCTURE_COLUMNS = {'Blood vessel': ('in_Blood vessel', 'crossed_frac_vessel'), 'Nerve': ('in_Nerve', 'crossed_frac_nerve'),
                     'Lymphatic': ('in_Lymphatic', 'crossed_frac_lymphatic')}


# -- the plan of one run: compartments, lineages and label rules from the registry ------------------------------------
def labels_table(reg: RG.Registry) -> pd.DataFrame:
    return pd.read_csv(reg.directory / f'{reg.name}__labels.csv', keep_default_na=False, dtype=str)


def used_markers(selection: pd.DataFrame, reg: RG.Registry) -> list[str]:
    """Typing markers the level selection kept (decision 'use'), in registry order."""
    sel = selection.set_index('marker')
    return [m for m in reg.typing_markers() if m in sel.index and sel.loc[m, 'decision'] == 'use']


def build_plan(reg: RG.Registry, markers: list[str]) -> dict:
    lab = labels_table(reg)
    ix = {m: i for i, m in enumerate(markers)}
    comps = [c for c in reg.compartments() if any(reg.table.loc[m, 'compartment'] == c for m in markers
                                                  if reg.role(m) == 'defining')]
    plan = {'markers': markers, 'ix': ix, 'compartments': comps, 'labels': lab, 'comp': {}}
    for c in comps:
        defining = [m for m in markers if reg.role(m) == 'defining' and reg.table.loc[m, 'compartment'] == c]
        carried = [m for m in markers if c in RG.split(reg.table.loc[m, 'shared_by']) and m not in defining]
        rows = lab[lab['compartment'] == c]
        pan_lin = set(rows.loc[rows['kind'] == 'pan', 'lineage'])
        pan = [m for m in defining if {lin for lin, _ in reg.lineage_entries(m)} <= pan_lin and reg.lineage_entries(m)]
        lineages = {}
        for m in defining:
            if m in pan:
                continue
            for lin, sub in reg.lineage_entries(m):
                L = lineages.setdefault(lin, {'specific': [], 'members': [], 'broad': [], 'subtypes': {}})
                L['members'].append(m)
                if {x for x, _ in reg.lineage_entries(m)} == {lin}:
                    L['specific'].append(m)
                if reg.lineage_entries(m) == [(lin, '')]:
                    L['broad'].append(m)
                if sub:
                    L['subtypes'][sub] = m
        lineages = {k: v for k, v in lineages.items() if v['specific']}
        residual = next(iter(rows.loc[rows['kind'] == 'residual', 'lineage']), '')
        kinds = set(rows['kind'])
        plan['comp'][c] = {'defining': defining, 'carried': carried, 'pan': pan, 'lineages': lineages,
                           'structure_only': bool(kinds) and kinds <= {'structure'},
                           'pan_label': next(iter(rows.loc[rows['kind'] == 'pan', 'label']), ''),
                           'pan_negative_label': next(iter(rows.loc[rows['kind'] == 'pan_negative', 'label']), ''),
                           'residual': residual,
                           'residual_members': [m for m in defining
                                                if residual in {lin for lin, _ in reg.lineage_entries(m)}]}
    return plan


def label_row(plan: dict, compartment: str, lineage: str) -> dict:
    lab = plan['labels']
    r = lab[(lab['compartment'] == compartment) & (lab['lineage'] == lineage) & ~lab['kind'].isin(['structure', 'pan_negative'])]
    return r.iloc[0].to_dict() if len(r) else {}


# -- canary status (Xenium cell 22) --------------------------------------------------------------------------------------
def group_status(detected: np.ndarray, depth: np.ndarray, score: np.ndarray, margin: float, min_depth: float,
                 bright: np.ndarray, cover: np.ndarray, lead: float = 0.0,
                 interior_valid: np.ndarray | None = None):
    """Per row: status, top group, runner-up group, their depths, the number of groups detected, and whether a tie in
    depth was settled. Two groups of equal depth (at least min_depth): the group whose detecting markers are brighter
    (median normalized interior signal) and win at least `lead` more of the cell's interior pixels is placed, resolved; otherwise the cell stays lineage_ambiguous (owner 2026-10-29: "both are
    brighter and cover more cell area"; later the same day "Interior coverage, clear lead", CANARY_TIE_INTERIOR_LEAD)."""
    n = len(detected)
    rows = np.arange(n)
    ng = detected.sum(axis=1)
    comp = np.where(detected, depth * 1000.0 + score, -np.inf)
    top = comp.argmax(axis=1)
    comp2 = comp.copy()
    comp2[rows, top] = -np.inf
    runner = comp2.argmax(axis=1)
    td = np.where(ng > 0, depth[rows, top], 0)
    rd = np.where(ng > 1, depth[rows, runner], 0)
    tie = (ng == 2) & (td == rd) & (td >= min_depth)
    if interior_valid is not None:
        tie &= interior_valid
    b_top, b_run, c_top, c_run = bright[rows, top], bright[rows, runner], cover[rows, top], cover[rows, runner]
    # A subtraction such as 0.35 - 0.20 can round just below an exact 0.15 lead.
    keep = tie & (b_top > b_run) & (c_top - c_run + 1e-12 >= lead) & (c_top > c_run)
    swap = tie & (b_run > b_top) & (c_run - c_top + 1e-12 >= lead) & (c_run > c_top)
    top, runner = np.where(swap, runner, top), np.where(swap, top, runner)
    settled = keep | swap
    status = np.select([ng == 0, ng == 1, ng >= 3, ((td - rd >= margin) & (td >= min_depth)) | settled],
                       ['none', 'single', 'doublet_suspected', 'resolved'], 'lineage_ambiguous').astype(object)
    return status, top, np.where(ng > 1, runner, -1), td, rd, ng, settled


def placed(P, C, F, detect, count, margin, min_depth, lead, pixels, source_rows, stage):
    """Use canary depth first, plain detection as fallback, and spatial evidence for eligible equal depths."""
    def assess(active, allowed, basis):
        detected = np.column_stack([active[:, d].any(axis=1) for d in detect])
        depth = np.column_stack([active[:, c].sum(axis=1) for c in count]) * detected
        score = np.column_stack([np.where(active[:, c], F[:, c], 0.).sum(axis=1) for c in count]) * detected
        eligible = (allowed & (detected.sum(axis=1) == 2) & (depth.max(axis=1) >= min_depth)
                    & ((depth == depth.max(axis=1)[:, None]).sum(axis=1) == 2))
        bright = np.full(score.shape, np.nan)
        cover = np.full(score.shape, np.nan)
        requests = []
        pair = {}
        for i in np.flatnonzero(eligible):
            a, b = np.flatnonzero(detected[i])
            row = source_rows[i]
            pair[row] = (i, a, b)
            requests.append((row, [m for m in detect[a] if active[i, m]],
                                  [m for m in detect[b] if active[i, m]]))
        measured = pixels.measure(requests, stage + '_' + basis) if pixels is not None else {}
        for row, (n, na, nb, ma, mb) in measured.items():
            i, a, b = pair[row]
            bright[i, [a, b]] = ma, mb
            cover[i, [a, b]] = na / n, nb / n
        return group_status(detected, depth, score, margin, min_depth, bright, cover, lead)

    status, top, runner, td, rd, ng, settled = assess(C, np.ones(len(C), bool), 'canary')
    statusp, topp, _, _, _, ngp, settledp = assess(P, ng == 0, 'plain')
    group = np.where(ng > 0, top, np.where(ngp > 0, topp, -1))
    basis = np.where(ng > 0, 'canary', np.where(ngp > 0, 'plain', '')).astype(object)
    # a cell placed below the canary level keeps its plain-level conflict visible (review 2026-10-28)
    status = np.where(ng > 0, status, np.where(ngp > 1, 'below_canary_' + statusp.astype(str), 'none')).astype(object)
    return {'group': group, 'basis': basis, 'status': status, 'n': ng, 'n_plain': ngp, 'depth': td,
            'runner': runner, 'runner_depth': rd, 'n_detected_at_basis': np.where(ng > 0, ng, ngp),
            'settled': np.where(ng > 0, settled, settledp)}


# -- one image -------------------------------------------------------------------------------------------------------------
def load_values(source, tag: str, reg: RG.Registry, markers: list[str], lv: pd.DataFrame,
                ld: pd.DataFrame | None = None) -> pd.DataFrame:
    cols = {m: reg.column(m, lv.loc[m, 'measurement']) for m in markers}
    diff_cols = {m: reg.column(m, ld.loc[m, 'measurement']) for m in markers if ld is not None and m in ld.index}
    d = source.identity_frame(tag, sorted(set(cols.values()) | set(diff_cols.values())),
                              [c for cols in STRUCTURE_COLUMNS.values() for c in cols])
    d.attrs['columns'] = cols
    d.attrs['diff_columns'] = diff_cols
    return d.reset_index(drop=True)


def established(calls: pd.DataFrame, tag: str, object_id: pd.Series, reg: RG.Registry, plan: dict, row) -> tuple:
    """The established typing's decisions in a lineage (row of the labels file), mapped to this run's labels: a subtype
    class (e.g. 'Beta') to the subtype's label, the lineage label as it is, its unresolved or ambiguous subtype classes
    ('<Lineage> (...)') to the lineage name, and its exclusions (reasons ESTABLISHED_REASONS) kept as Excluded.
    `calls`: the established typing's object_id, cell_class, excluded_reason (source.established_calls)."""
    c = pd.DataFrame({'object_id': object_id}).merge(calls, on='object_id', how='left', validate='one_to_one')
    cls, why = c['cell_class'].fillna('').to_numpy(dtype=object), c['excluded_reason'].fillna('').to_numpy(dtype=object)
    lin = row.lineage
    out = np.full(len(c), '', dtype=object)
    for sub in set(reg.subtypes(lin).values()):
        out[cls == sub] = label_row(plan, row.compartment, f'{lin}:{sub}')['label']
    out[cls == row.label] = row.label
    out[np.char.startswith(cls.astype(str), f'{lin} (')] = lin
    excl = (cls == 'Excluded') & np.isin(why, ESTABLISHED_REASONS) | \
        (cls == 'Excluded') & np.char.startswith(why.astype(str), ESTABLISHED_REASON_PREFIX)
    out[excl] = 'Excluded'
    keep = out != ''
    return keep, out, np.where(excl, why, '')


def default_start_class(d: pd.DataFrame) -> np.ndarray:
    """The CODEX start class: the parent endocrine object's class, else Tissue (the pancreas hook)."""
    return d['endocrine_object_class'].fillna('Tissue').to_numpy(dtype=object)


class IdentityError(ValueError):
    """A section cannot be typed until its levels or inputs are fixed."""


def classify_image(source, tag: str, reg_name: tuple, tables: dict, endocrine_from: str = ESTABLISHED,
                   start_class: Callable[[pd.DataFrame], np.ndarray] = default_start_class) -> tuple[pd.DataFrame, dict]:
    """`tables`: levels, levels_differentiation, selection (levels.choose) and canary_levels (canary.levels)."""
    with ExitStack() as stack:
        return _classify_image(source, tag, reg_name, tables, endocrine_from, start_class, stack)


def _classify_image(source, tag, reg_name, tables, endocrine_from, start_class, stack):
    reg = RG.load(*reg_name)
    markers = used_markers(tables['selection'], reg)
    state = reg.markers('state')
    lv_all = tables['levels']
    lv = lv_all[lv_all['image'] == tag].set_index('marker')
    cl = tables['canary_levels']
    cl = cl[cl['image'] == tag].set_index('marker')
    missing = [m for m in markers + state if not np.isfinite(lv.loc[m, 'level_log'])]
    missing += [m for m in markers if m not in cl.index or not np.isfinite(cl.loc[m, 'canary_log'])]
    # the differentiation readings (levels.py): the level of each marker compared on another measurement
    diff = [m for m in markers if reg.differentiation(m)]
    ld = None
    if diff:
        ld_all = tables['levels_differentiation']
        ld = ld_all[ld_all['image'] == tag].set_index('marker')
        missing += [m for m in diff if m not in ld.index or ld.loc[m, 'measurement'] != reg.differentiation(m)
                    or not np.isfinite(ld.loc[m, 'level_log'])]
    if missing:
        raise IdentityError(f'{tag}: no level or canary level for {sorted(set(missing))}: the method must be fixed first')
    plan = build_plan(reg, markers)
    ix = plan['ix']
    d = load_values(source, tag, reg, markers + state, lv, ld)
    n = len(d)
    V = np.column_stack([np.log1p(np.clip(np.nan_to_num(d[d.attrs['columns'][m]].to_numpy(float), nan=0.0), 0, None))
                         for m in markers])
    L = lv.loc[markers, 'level_log'].to_numpy(float)
    Lc = np.maximum(cl.loc[markers, 'canary_log'].to_numpy(float), L)
    P, C, F = V > L, V > Lc, V - L
    for m in diff:          # the fold competing groups are compared on: the differentiation reading
        vd = np.log1p(np.clip(np.nan_to_num(d[d.attrs['diff_columns'][m]].to_numpy(float), nan=0.0), 0, None))
        F[:, ix[m]] = vd - float(ld.loc[m, 'level_log'])
    # area coverage (area_coverage.py): proportion of the cell's pixels at or above the level, per marker
    cov_frame = source.coverage(tag)
    cov = AC.read_validated(cov_frame, source.redsea_ids(tag), tag, d['object_id'], markers)
    interior_count = cov[AC.INTERIOR_COUNT].to_numpy(np.int64)
    interior_valid = interior_count > 0
    lead = reg.const('CANARY_TIE_INTERIOR_LEAD')
    # the typing scale the hormone margin is read on (Senior islet_identity: top >= REFERENCE_MARGIN x runner-up on
    # scaled hormones), per image: (x - level) / (p99.9 - level), floored at 0 and not clipped above. Owner 2026-10-29:
    # ND_6232 ddb348e9 (GCG 2,054 / SST 118) was excluded when the folds over the levels were compared, because SST's
    # level is a ninth of GCG's; on the scale GCG is 3.5x SST, and ND_6232 af2ee599 (SST over INS) stays resolved.
    raw, Lraw = np.expm1(V), np.expm1(L)
    top = np.percentile(raw, reg.const('HORMONE_SCALE_TOP_QUANTILE'), axis=0)
    SC = np.clip((raw - Lraw) / np.where(top > Lraw, top - Lraw, np.nan), 0, None)

    window = source.native_window(tag) if hasattr(source, "native_window") else None
    pixels = (InteriorSignal(lambda: stack.enter_context(window), d["object_id"], interior_count, reg,
                             markers, Lraw, top) if window is not None else None)
    hormone_settled = np.zeros(n, dtype=bool)

    def highest_hormone(rows: np.ndarray, cols: list) -> tuple:
        """Senior subpanel step 5: the unique highest scaled hormone among those above their levels names the cell
        (owner 2026-10-29, no margin). Two positive hormones tied at the top can be resolved from native pixels."""
        Fs = np.where(P[np.ix_(rows, cols)], np.nan_to_num(SC[np.ix_(rows, cols)], nan=0.0), -np.inf)
        order = np.sort(Fs, axis=1)
        any_pos = np.isfinite(order[:, -1])
        unique = any_pos & (order[:, -1] > order[:, -2]) if len(cols) > 1 else any_pos
        winner = Fs.argmax(axis=1)
        tied = (Fs == order[:, -1, None]) & np.isfinite(Fs)
        requests = []
        positions = {}
        for i in np.flatnonzero(any_pos & (tied.sum(axis=1) == 2)):
            a, b = np.flatnonzero(tied[i])
            positions[rows[i]] = (i, a, b)
            requests.append((rows[i], [cols[a]], [cols[b]]))
        measured = pixels.measure(requests, 'hormone') if pixels is not None else {}
        for row, (size, na, nb, ma, mb) in measured.items():
            i, a, b = positions[row]
            if abs(na - nb) + 1e-12 < lead * size or (na - nb) * (ma - mb) <= 0:
                continue
            unique[i], winner[i], hormone_settled[row] = True, a if na > nb else b, True
        return any_pos, unique, winner
    start = np.asarray(start_class(d), dtype=object)
    margin, min_depth = reg.const('CANARY_MARGIN'), reg.const('CANARY_MIN_DEPTH')

    # compartments
    comps = plan['compartments']
    # a structure-defined compartment is reached only through the structure step: its markers detect nothing here
    detect = [[] if plan['comp'][c]['structure_only'] else [ix[m] for m in plan['comp'][c]['defining']] for c in comps]
    count = [[ix[m] for m in plan['comp'][c]['defining'] + plan['comp'][c]['carried']] for c in comps]
    cp = placed(P, C, F, detect, count, margin, min_depth, lead, pixels, np.arange(n), 'compartment')
    compartment = np.where(cp['group'] >= 0, np.asarray(comps, dtype=object)[np.maximum(cp['group'], 0)], '')
    comp_status = cp['status'].astype(object)

    # lineages inside the cell's compartment
    lineage = np.full(n, '', dtype=object)
    lin_basis = np.full(n, '', dtype=object)
    lin_status = np.full(n, '', dtype=object)
    lin_runner = np.full(n, '', dtype=object)
    lin_detected = np.zeros(n, int)
    lin_settled = np.zeros(n, bool)
    for c in comps:
        info = plan['comp'][c]
        lins = list(info['lineages'])
        rows = np.flatnonzero(compartment == c)
        if not lins or not len(rows):
            continue
        extra = info['pan'] + info['carried']
        det_cols = [[ix[m] for m in info['lineages'][lin]['specific']] for lin in lins]
        cnt_cols = [[ix[m] for m in dict.fromkeys(info['lineages'][lin]['members'] + extra)] for lin in lins]
        lp = placed(P[rows], C[rows], F[rows], det_cols, cnt_cols, margin, min_depth, lead,
                    pixels, rows, 'lineage_' + c)
        lin_settled[rows] = lp['settled']
        has = lp['group'] >= 0
        lineage[rows[has]] = np.asarray(lins, dtype=object)[lp['group'][has]]
        lin_basis[rows] = lp['basis']
        lin_status[rows] = np.where(has, lp['status'], '')
        lin_runner[rows] = np.where(lp['runner'] >= 0, np.asarray(lins, dtype=object)[np.maximum(lp['runner'], 0)], '')
        lin_detected[rows] = lp['n_detected_at_basis']

    # identity
    label = start.copy()
    reason = np.full(n, '', dtype=object)
    named_by = np.full(n, '', dtype=object)
    tags = [[] for _ in range(n)]
    for i in np.flatnonzero(cp['settled']):
        tags[i].append('compartment_tie_brightness_coverage')
    for i in np.flatnonzero(lin_settled):
        tags[i].append('lineage_tie_brightness_coverage')
    t_status = np.full(n, '', dtype=object)
    status = np.full(n, '', dtype=object)
    done = np.zeros(n, bool)
    names = np.asarray(markers, dtype=object)
    comp_basis = cp['basis'].copy()
    established_label = np.full(n, '', dtype=object)

    # the established typing of a lineage (labels file named_by 'established'; owner 2026-10-28: "The endocrine typing
    # didn't need to change", "If the singles were already found and verified, don't undo that"): its decisions are
    # kept as they were, and this hierarchy names no other cell in that lineage
    lab = plan['labels']
    # objects that are their endocrine cells (labels file defined_by_object; owner 2026-10-28: "the islets, clusters,
    # and single endocrine cells were already accurately identified", "If the singles were already found and verified,
    # don't undo that"): the object places the cell in its lineage; its markers name it, else it keeps the object class
    by_object = np.zeros(n, bool)
    for r in lab[lab['defined_by_object'] != ''].itertuples():
        rows = np.flatnonzero(np.isin(start, RG.split(r.defined_by_object)))
        for i in rows[(compartment[rows] != r.compartment) & (compartment[rows] != '')]:
            tags[i].append(f'canary_top_{compartment[i]}')
        compartment[rows], lineage[rows] = r.compartment, r.lineage
        comp_basis[rows], lin_basis[rows], lin_status[rows] = 'object', 'object', ''
        by_object[rows] = True
    # Structure objects (owner 2026-10-29, Senior's subpanel structure step): a cell whose centroid is in a vessel, nerve or
    # lymphatic object, with more than STRUCTURE_MAJORITY_AREA of its area inside that class's objects, the object's
    # marker positive and the highest fold over level among its positive typing markers, is that structure's cell type.
    # Asked weakest class first so a later class overrides (Blood vessel > Nerve > Lymphatic); before the canary
    # conflicts, which a structure call resolves. Object-defined endocrine cells (clusters, singles) are left as placed.
    # A cell in an object that is not called carries passed_through_<class>.
    Fpos = np.where(P, F, -np.inf)
    top_any = np.where(P.any(axis=1), Fpos.argmax(axis=1), -1)
    majority = reg.const('STRUCTURE_MAJORITY_AREA')
    structure_hit = np.zeros(n, bool)
    for r in list(lab[lab['kind'] == 'structure'].itertuples())[::-1]:
        inside_col, frac_col = STRUCTURE_COLUMNS[r.structure]
        inside = d[inside_col].fillna(False).to_numpy(bool)
        frac = d[frac_col].to_numpy(float)
        j = ix.get(r.named_by)
        if j is None:
            raise IdentityError(f'{tag}: structure marker {r.named_by} has no level')
        hit = inside & (np.nan_to_num(frac) > majority) & P[:, j] & (top_any == j) & ~by_object
        label[hit], named_by[hit] = r.label, r.named_by
        compartment[hit], lineage[hit] = r.compartment, r.lineage
        comp_basis[hit] = lin_basis[hit] = 'structure'
        lin_status[hit] = ''
        structure_hit |= hit
        for i in np.flatnonzero(inside & ~hit):
            tags[i].append(f"passed_through_{r.structure.replace(' ', '_')}")
    for i in np.flatnonzero(structure_hit):
        tags[i] = [t for t in tags[i] if not t.startswith('passed_through_')]
    done |= structure_hit
    # a structure-defined compartment's marker positive without that compartment's structure call names nothing
    # (owner 2026-10-30): the overlap stays a feature; tagged
    struct_cols = []
    for c in comps:
        if plan['comp'][c]['structure_only']:
            cols_c = [ix[m] for m in plan['comp'][c]['defining']]
            struct_cols += cols_c
            for i in np.flatnonzero(P[:, cols_c].any(axis=1) & ~(structure_hit & (compartment == c))):
                tags[i].append(f'{c.lower()}_marker_outside_structure')
    # A canary-level conflict between mutually exclusive compartments, or lineages of one compartment, that the margin
    # does not resolve names no cell type (owner 2026-10-28: "there is no such thing as a CD8+CD163+ cell"): the cell
    # is kept with the conflict as its class, never named as either side. Below the canary level a second group is not
    # a credible conflict (status below_canary_*). Cells their object places (clusters, singles) are left as placed.
    comp_runner = np.where(cp['runner'] >= 0, np.asarray(comps, dtype=object)[np.maximum(cp['runner'], 0)], '')
    unresolved = ('lineage_ambiguous', 'doublet_suspected')
    clash = ~by_object & ~done & np.isin(comp_status, unresolved)
    clash_lin = ~by_object & ~done & ~clash & np.isin(lin_status, unresolved)
    for mask, top_, run_, st in ((clash, compartment, comp_runner, comp_status), (clash_lin, lineage, lin_runner, lin_status)):
        rows = np.flatnonzero(mask)
        label[rows] = np.where(st[rows] == 'doublet_suspected', 'Doublet suspected',
                               'Conflict ' + top_[rows].astype(str) + '/' + run_[rows].astype(str))
        named_by[rows] = 'canary'
        done[rows] = True
    for r in lab[lab['named_by'] == 'established'].itertuples():
        keep, est_label, est_reason = established(source.established_calls(tag, endocrine_from), tag, d['object_id'],
                                                  reg, plan, r)
        for i in np.flatnonzero(keep & (compartment != r.compartment) & (compartment != '')):
            tags[i].append(f'canary_top_{compartment[i]}')
        label[keep], reason[keep] = est_label[keep], est_reason[keep]
        established_label[keep] = est_label[keep]
        named_by[keep] = 'established'
        compartment[keep], lineage[keep] = r.compartment, r.lineage
        comp_basis[keep], lin_basis[keep], lin_status[keep] = 'established', 'established', ''
        done[keep] = True
        # Which cells are endocrine stays established; which hormone names them is read again with the one measurement
        # per marker and each hormone's per-image level (owner 2026-10-28: an SST-high cell excluded "under 2x the
        # runner-up"; the established comparison used cohort-scaled cytoplasm values, the "two ways of typing" the
        # 10-28 ruling ended). The unique highest scaled hormone names the cell (highest_hormone; owner 2026-10-29), an
        # unresolved spatial tie keeps the start class (tied_highest_hormones); no hormone above its level -> the lineage label when
        # the broad marker is positive (CHGA+ hormone-negative), else the start class.
        info = plan['comp'][r.compartment]
        L_ = info['lineages'][r.lineage]
        subs = {sname: ix[m] for sname, m in L_['subtypes'].items()}
        reread = keep & (np.isin(est_label, [label_row(plan, r.compartment, f'{r.lineage}:{x}')['label'] for x in subs]
                                 + [r.label, r.lineage]) | np.isin(est_reason, MARGIN_REASONS))
        rows = np.flatnonzero(reread)
        if len(rows) and subs:
            cols = list(subs.values())
            anypos, clear, w = highest_hormone(rows, cols)
            win = np.asarray(list(subs), dtype=object)[w]
            new = np.full(len(rows), '', dtype=object)
            why = np.full(len(rows), '', dtype=object)
            for x in subs:
                new[anypos & clear & (win == x)] = label_row(plan, r.compartment, f'{r.lineage}:{x}')['label']
            new[anypos & ~clear] = start[rows[anypos & ~clear]]
            status[rows[anypos & ~clear]] = 'tied_highest_hormones'
            broad = P[np.ix_(rows, [ix[m] for m in L_['broad']])].any(axis=1) if L_['broad'] else np.zeros(len(rows), bool)
            new[~anypos & broad] = r.label
            nameless = ~anypos & ~broad
            new[nameless] = start[rows[nameless]]
            status[rows[nameless]] = 'no_hormone_above_level'
            label[rows], reason[rows] = new, why
            nb = np.asarray([markers[subs[x]] for x in win], dtype=object)
            named_by[rows] = np.where(anypos & clear, nb, np.where(~anypos & broad, ';'.join(L_['broad']), 'established'))
            for i, j in zip(rows, range(len(rows))):
                if new[j] != established_label[i]:
                    tags[i].append('hormone_reread')

    def name_lineage(rows: np.ndarray, c: str, lin: str, L_: dict, rule: dict) -> None:
        """Name `rows` by lineage `lin` of compartment `c` (its broad marker and the registry's requires, then its
        subtype rule); a cell that fails the requirements is left unnamed (or keeps its object's class)."""
        broad = [ix[m] for m in L_['broad']]
        subs = {s: ix[m] for s, m in L_['subtypes'].items()}
        if rule.get('named_by') == 'broad':
            ok = P[np.ix_(rows, broad)].any(axis=1) if broad else np.zeros(len(rows), bool)
        else:                                   # broad_or_subtype
            ok = P[np.ix_(rows, broad + list(subs.values()))].any(axis=1)
        for m in L_['broad']:
            for req in reg.requires(m):
                if req == 'single_lineage':
                    ok &= lin_detected[rows] == 1
                elif req in ix:
                    ok &= P[rows, ix[req]]
        kept = rows[~ok & by_object[rows]]          # their object names them: they keep its class
        done[kept] = True
        status[kept] = f'no_{lin.lower()}_marker_above_level'
        for i in rows[~ok & ~by_object[rows]]:
            tags[i].append(f'{lin}_not_named')
        rows = rows[ok]
        srule = rule.get('subtype_rule', '')
        if srule == 'unique_highest' and subs:
            any_pos, named, w = highest_hormone(rows, list(subs.values()))
            win = np.asarray(list(subs), dtype=object)[w]
            for s in subs:
                sel = rows[named & (win == s)]
                label[sel] = label_row(plan, c, f'{lin}:{s}')['label']
                named_by[sel] = markers[subs[s]]
            tied = rows[any_pos & ~named]
            label[tied] = start[tied]
            status[tied] = 'tied_highest_hormones'
            neg = rows[~any_pos]
            label[neg] = rule['label']
            named_by[neg] = ';'.join(L_['broad'])
        elif srule == 'higher' and subs:
            cols = list(subs.values())
            Fs = np.where(P[np.ix_(rows, cols)], F[np.ix_(rows, cols)], -np.inf)
            any_pos = np.isfinite(Fs.max(axis=1))
            win = np.asarray(list(subs), dtype=object)[Fs.argmax(axis=1)]
            for s in subs:
                sel = rows[any_pos & (win == s)]
                label[sel] = label_row(plan, c, f'{lin}:{s}')['label']
                named_by[sel] = markers[subs[s]]
            none_ = rows[~any_pos]
            label[none_] = rule['label']
            named_by[none_] = ';'.join(L_['broad'])
            t_status[none_] = f'{lin}_subtype_not_shown'
            both = rows[P[np.ix_(rows, cols)].sum(axis=1) > 1]
            t_status[both] = 'both_subtype_markers_higher_named'
        elif srule == 'profile':
            head, tail = (RG.split(x.replace(',', ';')) for x in rule['profile'].split('|'))
            sign = lambda m: np.where(P[rows, ix[m]], f'{m}+', f'{m}-') if m in ix else np.full(len(rows), f'{m}?')
            parts = [sign(m) for m in head] + [np.full(len(rows), f"{rule['label']}_")] + [sign(m) for m in tail]
            label[rows] = np.asarray([''.join(p) for p in zip(*parts)], dtype=object)
            named_by[rows] = ';'.join(head)
        else:
            label[rows] = rule['label']
            named_by[rows] = ';'.join(L_['broad'])
        done[rows] = True


    for c in comps:
        info = plan['comp'][c]
        in_c = (compartment == c) & ~done
        for lin, L_ in info['lineages'].items():
            rule = label_row(plan, c, lin)
            if not rule or rule.get('kind') == 'residual':       # the residual lineage is labeled below
                continue
            rows = np.flatnonzero(in_c & (lineage == lin))
            if not len(rows):
                continue
            only_in = set(RG.split(rule.get('only_in', '')))
            # a lineage named only by its established typing, or outside its objects, names no cell here: the cell is
            # typed by its other evidence of the compartment (the residual lineage's markers), else it keeps its start
            # class; tagged either way
            barred = np.ones(len(rows), bool) if rule.get('named_by') == 'established' else \
                (~np.isin(start[rows], list(only_in)) if only_in else np.zeros(len(rows), bool))
            out_rows = rows[barred]
            for i in out_rows:
                tags[i].append(f'{lin.lower()}_marker_positive' if start[i] in only_in else f'{lin.lower()}_outside_object')
            members = [ix[m] for m in info['residual_members']]
            other = P[np.ix_(out_rows, members)].any(axis=1) if members else np.zeros(len(out_rows), bool)
            lineage[out_rows[other]] = info['residual']
            # the lineage's own markers name nothing here (tiles 2026-10-28: 2-6 of 40 INS+/GCG+/SST+ cells showed the
            # hormone on the cell), and no marker names a cell (owner 2026-10-28): it keeps its start class
            kept = out_rows[~other]
            done[kept] = True
            status[kept] = f'{lin.lower()}_marker_only'
            rows = rows[~barred]
            if not len(rows):
                continue
            name_lineage(rows, c, lin, L_, rule)
        # the compartment's cells without a named lineage
        rows = np.flatnonzero((compartment == c) & ~done)
        if not len(rows):
            continue
        if info['residual']:
            rule = label_row(plan, c, info['residual'])
            excl = np.isin(start[rows], RG.split(rule.get('excluded_in', '')))
            label[rows[~excl]] = rule['label']
            named_by[rows[~excl]] = 'compartment'
            label[rows[excl]] = 'Excluded'
            reason[rows[excl]] = f"{info['residual'].lower()}_inside_islet_border"
            done[rows] = True
            continue
        if info['pan']:
            pan_pos = P[np.ix_(rows, [ix[m] for m in info['pan']])].any(axis=1)
            label[rows[pan_pos]] = info['pan_label']
            named_by[rows[pan_pos]] = ';'.join(info['pan'])
            done[rows[pan_pos]] = True
            rows = rows[~pan_pos]
        done[rows] = True                          # no marker names a cell (owner 2026-10-28): the start class stays,
        status[rows] = 'no_lineage_named'          # or the compartment's pan-negative label (owner 2026-10-29)
        if info['pan_negative_label']:
            named = rows[~by_object[rows]]
            label[named] = info['pan_negative_label']
            named_by[named] = 'compartment'
    # no compartment: the start class stays, with the status saying why
    left = ~done
    if (left & (compartment != '')).any():
        raise ValueError(f'{tag}: {int((left & (compartment != "")).sum())} cells of a compartment were not placed')
    structure_marker = P[:, struct_cols].any(axis=1) if struct_cols else np.zeros(n, bool)
    status[left & structure_marker] = 'structure_marker_only'
    status[left & ~structure_marker & P.any(axis=1)] = 'shared_marker_only'
    status[left & ~P.any(axis=1)] = 'no_marker_above_background'

    # confirmation
    confirms = {m: RG.split(reg.table.loc[m, 'confirms']) for m in markers} if 'confirms' in reg.table.columns else {}
    confirmed = np.full(n, np.nan)
    confirming = np.full(n, '', dtype=object)
    lab_lineage = lineage.copy()
    for c in comps:
        info = plan['comp'][c]
        if info['residual']:
            lab_lineage[(compartment == c) & (label == label_row(plan, c, info['residual'])['label'])] = info['residual']
    placed_ok = np.isin(comp_status, PLACED) & (np.isin(lin_status, PLACED) | (lin_status == ''))
    for lin in sorted(set(lab_lineage) - {''}):
        cand = [m for m, ls in confirms.items() if lin in ls]
        if not cand:
            continue
        rows = np.flatnonzero(lab_lineage == lin)
        for i in rows:
            used = set(str(named_by[i]).split(';'))
            use = [ix[m] for m in cand if m not in used]
            if not use:
                continue
            hit = C[i, use]
            confirmed[i] = float(placed_ok[i] and hit.any())
            if hit.any():
                confirming[i] = ';'.join(names[np.asarray(use)[hit]])

    for i in np.flatnonzero(hormone_settled):
        tags[i].append("hormone_tie_brightness_coverage")

    out = pd.DataFrame({
        'object_id': d['object_id'].to_numpy(), 'start_class': start, 'label': label, 'excluded_reason': reason,
        'named_by': named_by, 'compartment': compartment, 'compartment_basis': comp_basis,
        'compartment_status': comp_status, 'n_compartments_canary': cp['n'], 'n_compartments_plain': cp['n_plain'],
        'compartment_depth': cp['depth'], 'runner_up_compartment': np.where(cp['runner'] >= 0, np.asarray(
            comps, dtype=object)[np.maximum(cp['runner'], 0)], ''), 'runner_up_depth': cp['runner_depth'],
        'lineage': lineage, 'lineage_basis': lin_basis, 'lineage_status': lin_status,
        'subtype_status': t_status, 'status': status, 'tags': [';'.join(t) for t in tags],
        'confirmed': confirmed, 'confirming_markers': confirming, 'established_label': established_label,
        'n_positive': P.sum(axis=1).astype(np.int16),
        AC.INTERIOR_COUNT: interior_count,
        'n_canary': C.sum(axis=1).astype(np.int16)})
    out["spatial_conflict_evidence"] = pixels.audit(n) if pixels is not None else ""
    for m in state:
        x = np.log1p(np.clip(np.nan_to_num(d[d.attrs['columns'][m]].to_numpy(float), nan=0.0), 0, None))
        out[f'{m}_pos'] = (x > lv.loc[m, 'level_log']).astype(np.int8)
    for j, m in enumerate(markers):
        out[f'pos__{m}'] = P[:, j].astype(np.int8)
        out[f'canary__{m}'] = C[:, j].astype(np.int8)
        out[f'fold__{m}'] = F[:, j].astype(np.float32)
    audit = {'image': tag, 'cells': n, 'markers': markers,
             'readings': {m: {'positivity': str(lv.loc[m, 'measurement']),
                              **({'differentiation': str(ld.loc[m, 'measurement'])} if m in diff else {})}
                          for m in markers + state},
             'structure_only_compartments': [c for c in comps if plan['comp'][c]['structure_only']],
             'cells_without_interior': int((~interior_valid).sum()),
             'label': {str(k): int(v) for k, v in out['label'].value_counts().items()},
             'compartment_status': {str(k): int(v) for k, v in out['compartment_status'].value_counts().items()},
             'lineage_status': {str(k): int(v) for k, v in out['lineage_status'].value_counts().items()},
             'status': {str(k): int(v) for k, v in pd.Series(status[status != '']).value_counts().items()},
             'no_marker_above_background': int((status == 'no_marker_above_background').sum()),
             'excluded': {str(k): int(v) for k, v in out.loc[out['label'] == 'Excluded', 'excluded_reason']
                          .value_counts().items()}}
    return out, audit
