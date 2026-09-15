# What the subpanel classification actually does

Recovered 2026-09-12 by reproducing `subpanel/subpanel_calls.parquet` from its inputs, because
nothing in the repository writes that file. Everything below was **checked against the calls on
disk**, not inferred from prose. Stopped short of bit-exactness deliberately: the QuPath export is
being regenerated, so the remaining residual would be chasing data that is about to change.

## Input chain — reproduced exactly

```
data/development/celltype/uniform/donor_id=<d>/data.parquet     UniFORM-aligned per-cell values
  -> apply_scale(floor, span)                                   broadtype/scratch/cosine_scale.parquet
  -> reference level = pooled p99 of max(subpanel) over cells with context == tissue
```

Pooled over all 20 donors this returns the four levels in `subpanel_log.json` to four decimals:

| subpanel | reproduced | recorded |
|---|---|---|
| hormone (INS/GCG/SST) | 0.1010 | 0.1010 |
| Nerve (B3TUBB/CD56) | 0.2272 | 0.2272 |
| Lymphatic (Podoplanin) | 0.1776 | 0.1776 |
| Blood vessel (Caveolin/CD34) | 0.3073 | 0.3073 |

The scale is **read, not recomputed** — `build_scale` draws a seeded 200,000-cell sample per donor,
so recomputing it would drift. One cohort floor (the aligned background mode) and one top (cohort
p99.9) per marker.

## Context — reproduced exactly, 100.0000% on 6442, 6436 and 115

Centroid-in-object against `results/Geojson/`, in pixel coordinates via each image's own pixel size.
Islet membership comes from the islet OUTLINES, by testing each centroid against them (`islet_membership.py`).
It came from QuPath's `Parent` column until 2026-09-15; that was wrong for 11,556 islet-core cells
(2.7 %; 35.3 % on 6539), which were never asked the islet question and fell through to Exocrine.

**Precedence, recovered from the calls rather than assumed:**

```
islet  >  Blood vessel  >  Nerve  >  Lymphatic  >  tissue
```

## The decision chain

**Inside an islet** — reproduced at 100.0000% on 6442 (19,675 cells), 1-2 cells off on 115 and 6436,
all of them `Excluded`:

1. centroid **inside** a structure object **and** that structure's markers clear its level -> that
   structure. Precedence as above.
2. else `max(INS, GCG, SST) >= 0.1010` -> argmax of the three names Beta / Alpha / Delta.
   **CD99 is not in the argmax** - three variants were tested and only this one fits.
3. else fall through to the tissue chain.

So the structure question **is** asked inside islets, with the object evidence. `context` simply
does not record it: the precedence overwrites the object with `islet`, which is why every one of
the 652,871 islet-core cells reads `context = islet`. The pass-through is live there:

```
donor  structure       in object   clears level   called   passed through
6442   Blood vessel        1,969          1,312    1,312        33.4 %
115    Blood vessel        6,108          1,659    1,659        72.8 %
6534   Lymphatic           6,206            108      100        98.3 %
```

**Outside an islet** — 99.52% from the upper tier alone, with **zero false positives**:

1. structure context -> its markers clear its level -> that structure; else fall through.
2. the immune definitions from `celltype/rules.requirements(..., with_negatives=True)`, evaluated at
   `require_min = 0.10` on the scaled values. Note the rule table, not the prose: `Imm-macrophage`
   is CD68 **or** Iba1 (an `any_of`), not both; `Imm-neutrophil` is CD66 **and** (CD11b or MPO).
3. a supporting tier over `bottom_tier_immune` at a level of **0.8715**, cohort-wide - the same
   value recovers independently on 6442 (0.8715) and 6436 (0.8716).
4. else Exocrine.

**Recovered 2026-09-15**: the supporting tier is eligible for a cell whose immune evidence
**clears** 0.8715 and whose exocrine evidence **does not** -- the same level on both sides, so the
tier never overrules a cell the exocrine markers speak for as loudly as the immune ones do.

All 35,332 cells the original table claims by this tier on 6436 carry exocrine evidence below
0.8715; none is above it. The two rejected alternatives are both measured: applying the tier to
every unclaimed cell wrongly claims 5,671 strongly exocrine cells on 6436 (exocrine max p50 0.98),
and restricting to cells with *no* exocrine evidence keeps only 19,142 of the 35,332.

**Still open, and small**: 1,223 cells on 6436 (0.08 %), 98.4 % of them driven by PD_1 on cells
whose exocrine evidence sits just under the level -- p50 0.843, against 0.044 for the cells the
tier should claim. Closing it means choosing a threshold because it fits the answer, which is
fitting rather than recovery, so it is left open.

## The gate, 2026-09-15

`subpanel.py --gate` rebuilds with `--membership parent` and compares row-for-row against the
table on disk. Every donor: identical row count, no cell only on one side.

| donor | cell_class | context | cells wrong |
|---|---|---|---|
| 6442 | 99.9995 % | 100.0000 % | 7 of 1,419,669 |
| 115 | 99.9990 % | 100.0000 % | 9 of 860,863 |
| 6436 | 99.9224 % | 100.0000 % | 1,226 of 1,579,440 |

All four reference levels reproduce to four decimals. The whole residual is the PD_1 set above.

## Why this matters after the re-export

None of the above is a fact about the current run. It is the method, and the method is what gets
re-run over the corrected export. The two modules beside this file - `subpanel_io.py` and
`context.py` - carry the parts that were verified exactly.

## Two marker facts about the clean batch-1 cohort

Measured while checking what the panel supports there; recorded because they bear on how any
immune number from these donors should be read, and neither is derivable from the code.

**CD38 is present on all ten.** Every clean-cohort donor is a `panel-cd38-59ch` donor, so the
`Imm-plasma-cell` rule runs as the table writes it. The CD79a+/CD20- stand-in that batch 2 needs is
not in force here.

**CD57 is effectively absent.** Its scaled p99 runs 0.0035-0.0266 across the ten against a 0.10
positivity level -- it never fires. So on this cohort `Imm-NK-cell` reduces to CD56 with Granzyme B.
That is an absence, not the b-Catenin1 reagent-change confound (which is a batch-2 property): within
batch 1 CD57 is *consistent*, just consistently at background.

## Not revisited

The donors here are already cell typed and that typing stands. Nothing in the app build re-derives
`cell_class`: `build_cells.py` reads it from `subpanel_calls.parquet` and maps it to a display name,
and the only cells it removes are the `exocrine_inside_islet_border` set the maintainer ruled out of
analysis. A script that re-derived immune subtypes from the rule table was written and removed
rather than left runnable, since the established workflow is not to be altered.
