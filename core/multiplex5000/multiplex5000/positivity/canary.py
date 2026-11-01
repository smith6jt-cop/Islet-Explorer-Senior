"""canary.py -- each typing marker's canary detection level per image (final typing plan of 2026-10-28, Step 5).

Ported from Islet-Explorer scripts/uniform/typing/canary.py at 578767e (``positivity.SOURCE_COMMIT``). The computation is
unchanged; the levels tables are passed in (``tables``: the dict ``levels.choose`` returns) and cells come from a source.

A marker counts toward canary depth (identity.py) only above a stricter level than plain positivity, as the Xenium
canary step counts a gene only above a stricter level than its ordinary detection. The plan set that level at 0-3
background widths above the positivity level, chosen on the 7 ND images against a null (each marker's positive rate in
compartments that cannot carry it, in cells farther than 20 um from any structure object). Two parts of that failed on
these data (agent determination 2026-10-28, under the owner's rule "flag the methods you are using so they will work
with the images and markers"): a background width cannot be read for every marker and image (CD20 on 6 images,
E-cadherin on 2, CD68 on 2: the background is a plateau), and leukocytes farther than 20 um from any structure number
82-1,395 on the ND images, too few for a per-image percentile. The method used instead needs neither:

  canary level = max(the marker's positivity level, the CANARY_NULL_QUANTILE-th percentile of the marker among the
                 image's cells of every compartment that cannot carry it)

the known-negative rule the owner approved for SMA (Barlow et al. eLife 2024), applied to the canary step. Compartments
that can carry a marker: its registry compartment and shared_by. Reference cells (not isolated: a neighbor's
spillover is what the canary level must discount): Epithelial = endocrine-object hormone-positive cells and ordinary-
tissue cells positive for an Epithelial majority marker that are not leukocytes (levels.populations 'hormone' and
'exocrine'); Hematopoietic = leukocytes (the Leukocyte marker with a T or Myeloid marker). A marker carried by both has
its registry check population instead (LGALS3: endocrine-object hormone-positive cells). A reference needs KNOWN_NEGATIVE_MIN_ABOVE cells above the percentile, or the run stops.
Which cells are leukocytes or epithelial cells is read on each marker's differentiation reading (typing/
levels_differentiation.csv; owner 2026-10-30: "Use membrane for differentiating, but then the whole cell still needs
to be used for positivity"), as levels.py defines its reference cells; the percentile itself is read on the marker's
positivity measurement (typing/levels.csv).

``levels`` returns canary_levels: image, marker, level_log, canary_log, canary_fold_over_level, null_reference,
null_cells, null_rate_at_level (the plain level's rate in the null).
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

import numpy as np
import pandas as pd

from . import identity as ID
from . import levels as LV
from . import registry as RG

REFERENCE_COMPARTMENTS = ('Epithelial', 'Hematopoietic')


def image_levels(source, tag: str, reg_name: tuple, tables: dict) -> list[dict]:
    reg = RG.load(*reg_name)
    R = LV.roles(reg)
    markers = ID.used_markers(tables['selection'], reg)
    lv = tables['levels']
    lv = lv[lv['image'] == tag].set_index('marker')
    d = LV.load(source, tag, reg)
    xy = d[['x_um', 'y_um']].to_numpy(float)
    x = {m: np.clip(np.nan_to_num(d[reg.column(m, lv.loc[m, 'measurement'])].to_numpy(float), nan=0.0), 0, None)
         for m in markers}
    pos = {m: np.log1p(x[m]) > lv.loc[m, 'level_log'] for m in markers}
    # reference cells on the differentiation readings (levels.py)
    ld = tables['levels_differentiation']
    ld = ld[ld['image'] == tag].set_index('marker')
    ref_pos = dict(pos)
    for m in markers:
        if reg.differentiation(m):
            if m not in ld.index or ld.loc[m, 'measurement'] != reg.differentiation(m):
                raise LV.PositivityError(f'{tag} {m}: no level on its differentiation measurement')
            v = np.clip(np.nan_to_num(d[reg.column(m, ld.loc[m, 'measurement'])].to_numpy(float), nan=0.0), 0, None)
            ref_pos[m] = np.log1p(v) > ld.loc[m, 'level_log']
    q = reg.const('CANARY_NULL_QUANTILE')
    need = reg.const('KNOWN_NEGATIVE_MIN_ABOVE')
    cut = reg.const('ISOLATION_UM')
    rows = []
    hormone, _ = LV.cell_sets(d, ref_pos, R)
    epithelial_any = LV.any_of(ref_pos, [x for x in R['epithelial'] if x in ref_pos], len(d)) & d['ordinary'].to_numpy()
    hematopoietic = [x for x in markers if reg.role(x) == 'defining' and reg.table.loc[x, 'compartment'] == 'Hematopoietic']
    for m in markers:
        pops = LV.populations(d, reg, m, ref_pos, xy, cut, R)   # the same populations as the levels' references
        # Epithelial null: endocrine cells and ordinary-tissue epithelial cells, not isolated (the canary level must
        # discount a neighbor's spillover), leaving out cells positive for a hematopoietic marker other than m, so a
        # hematopoietic marker's null is never defined by its own negativity (review 2026-10-28)
        # leave out cells positive for a marker of a lineage mutually exclusive with m's (never the pan-leukocyte
        # marker every carrier also has, nor m's own lineage; for the pan marker itself, every lineage marker); a
        # second review (2026-10-28) found the first form still removed the carriers of every marker but CD45
        # a shared marker has no lineage of its own: its carriers (registry carrier_lineages) stand in, so cells of a
        # carrier lineage stay in its null (third review, 2026-10-28: HLADR, CD35)
        lin_m = {lin for lin, _ in reg.lineage_entries(m)} or set(RG.split(reg.table.loc[m, 'carrier_lineages']))
        pan = R['leukocyte']
        others = [x for x in hematopoietic if x not in pan and x != m and
                  (m in pan or not ({lin for lin, _ in reg.lineage_entries(x)} & lin_m))]
        refs = {'Epithelial': (hormone | epithelial_any) & ~LV.any_of(ref_pos, others, len(d)),
                'Hematopoietic': pops[reg.table.loc[m, 'check_reference']]
                if reg.table.loc[m, 'check_reference'] in ('t_cells', 'leukocyte') else pops['leukocyte']}
        carriers = {reg.table.loc[m, 'compartment']} | set(RG.split(reg.table.loc[m, 'shared_by']))
        nulls = [c for c in REFERENCE_COMPARTMENTS if c not in carriers]
        if nulls:
            ref = np.logical_or.reduce([refs[c] for c in nulls])
            name = ' + '.join(nulls)
        else:                                  # carried by both: the marker's own check population (registry)
            name = reg.table.loc[m, 'check_reference']
            if not name:
                raise LV.PositivityError(f'{m}: carried by every reference compartment and no check_reference in the '
                                         'registry')
            ref = pops[name]
        o = x[m][ref]
        if len(o) * (1 - q / 100) < need:
            raise LV.PositivityError(f'{tag} {m}: {len(o)} reference cells are too few for the {q:g}th percentile')
        level = float(lv.loc[m, 'level_log'])
        canary = max(level, float(np.log1p(np.percentile(o, q))))
        rows.append({'image': tag, 'marker': m, 'level_log': level, 'canary_log': canary,
                     'canary_fold_over_level': float(np.expm1(canary) / np.expm1(level)) if level > 0 else np.nan,
                     'null_reference': name, 'null_cells': int(ref.sum()),
                     'null_rate_at_level': float(pos[m][ref].mean())})
    return rows


def levels(source, reg_name: tuple, tags: list[str], tables: dict, workers: int = 1) -> pd.DataFrame:
    """canary_levels for every section (canary.py main(), returning the table it wrote)."""
    with ProcessPoolExecutor(max(1, min(workers, len(tags)))) as ex:
        rows = [r for rs in ex.map(image_levels, repeat(source), tags, repeat(reg_name), repeat(tables)) for r in rs]
    return pd.DataFrame(rows)
