# What the cell typing actually does

Two generations of calls are described here, and every statement says which one it is about,
because a reader in a month has to be able to tell them apart:

* **[ON DISK]** — about calls that already exist: `typing/subpanel/subpanel_calls.parquet`
  (25,150,537 rows, **superseded**) and `typing/cells/` (25,230,109 rows, **current**). Everything
  so marked was **recovered 2026-09-12 by reproducing those calls from their inputs** and checked
  against them, not inferred from prose.
* **[CHAIN NOW]** — the method as `subpanel.py` runs it after the maintainer's islet ruling of
  2026-09-18. Proven first on three donors (115, 6442, 6534) with `--membership geometry --levels
  per-donor`, then **run on all 20 on 2026-09-19**: `typing/subpanel_rebuilt/` (written
  09:41–09:46) and `typing/cells/` (09:55–09:59) carry it. Measured 2026-09-24: 25,230,109 rows
  over 20 donors in each, every partition with `islet_route` (cd99 408,950 · elongation 29,060 ·
  hormone 190,571) and `islet_stratum`, and no `Unresolved` anywhere.

  **Read the date on an `[ON DISK]` statement.** One written before 2026-09-19 describes the
  pre-chain build of 2026-09-15 that this rewrite replaced — the same 25,230,109 rows, but no
  `islet_cells`, `islet_stratum` or `islet_route` column, `resolved_by` carrying only
  `structure` / `definition` / `supporting` / `fallthrough` / `hormone`, and no CD99 question. A
  `[CHAIN NOW]` number quoted from the three-donor proving run is the proving run's, not the
  stage root's.

The refinement pass has its own file, [`REFINEMENT.md`](REFINEMENT.md), and a short section here.
It is summarised here on purpose: that pass ran on 2026-09-11, no producer for it existed anywhere
in the repository, this file did not mention it, and it survived only in a session transcript until
it was recovered on 2026-09-18. An omission from this file is how a method gets lost.

Nothing below is a verification target. The measured numbers are recorded because they are the
reason a rule exists or the size of what a rule moves — never as something a rebuild has to
reproduce.

## Who produces what

```
data/development/celltype/uniform/donor_id=<d>/data.parquet   UniFORM-aligned per-cell values
  -> subpanel.py     the classification                  -> data/development/typing/subpanel_rebuilt/
  -> build_cells.py  carries the producer's columns,      -> data/development/typing/cells/
                     adds regions, QC and display names
  -> refine.py       the refinement pass, a SIDECAR       -> <--out>/donor_id=<d>/data.parquet
```

`subpanel.py` is the producer that did not exist when this file was first written: the 25.15M-row
table had no script, so no correction to the typing could ever be applied. `refine.py` does not
rewrite `typing/cells/`; it writes a keyed sidecar beside it.

## Input chain — reproduced exactly  [ON DISK]

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

Pooled levels reproduce the old table and are **not the method**: `--levels per-donor` is, and it is
the default. The same reference population fitted per donor moves the hormone level 19.5x across the
cohort and separates almost perfectly by acquisition batch, which is aliased with disease status.

## Context — reproduced exactly, 100.0000% on 6442, 6436 and 115  [ON DISK]

Centroid-in-object against `results/Geojson/`, in pixel coordinates via each image's own pixel size.
Islet membership comes from the islet OUTLINES, by testing each centroid against them
(`islet_membership.py`). It came from QuPath's `Parent` column until 2026-09-15; that was wrong for
11,556 islet-core cells (2.7 %; 35.3 % on 6539), which were never asked the islet question and fell
through to Exocrine.

**Precedence, recovered from the calls rather than assumed:**

```
islet  >  Blood vessel  >  Nerve  >  Lymphatic  >  tissue
```

**The precedence is a loss, and it is now recorded rather than reconstructed later.** [CHAIN NOW]
Every islet-core cell reads `context = islet`, so the object answer disappears from the table —
6,116 of donor 115's islet cells sit inside a blood-vessel object, 55 in a lymphatic and 111 in a
nerve, and none of that was visible. `context.object_label()` returns the object answer with islets
not consulted, and `object_context` is written on **every** row, not only under `--emit-detail`. It
is not a second membership column; it already existed and was only being withheld by a flag.

## The islet chain  [CHAIN NOW]

The maintainer's ruling of 2026-09-18, in substance: CD99 is the key differentiator between
endocrine and non-endocrine inside an islet; the typical endocrine cell **in the larger islets** is
also elongated; the classified objects separate vessels and nerves; exocrine markers get into an
islet annotation only by mistake; and reference levels are **per donor AND per size stratum**.

> **Corrected 2026-09-18.** An earlier draft of this section read the ruling down to "per donor"
> for CD99, on the argument that a marker's reference population is ordinary tissue and ordinary
> tissue carries no islet size. That is not the ruling and it is not what `subpanel.py` does.
> The stratum applies to CD99 as well, because the population the endocrine level is fitted from
> is **not** ordinary tissue: it is that donor's hormone-confirmed cells *inside islets*, where
> size is a property of the thing being measured. Reading it down was a method decision, and
> method decisions are the maintainer's.

For a cell whose centroid is inside an islet outline:

1. **Structure.** Inside a classified object **and** clearing that object's level **and** the
   object's own marker is the cell's top-expressed identity marker. Three conditions, not one: an
   object is a region of tissue and most of what lies in it is whatever the region runs through.
2. **Exocrine inside an islet is an annotation error** — `Excluded`, with the reason on the row.
   Applied at the end of the chain, where the evidence is in.
3. **Endocrine or not, on CD99** — the key differentiator — at a level fitted **PER DONOR AND PER
   SIZE STRATUM**, from that donor-and-stratum's own hormone-confirmed cells, at
   `ENDOCRINE_QUANTILE` (the 1st percentile; see "The levels" below). A **FLOOR GUARD** refuses a
   level at or below the feature's own floor — `0.0` for a scaled marker, `1.0` for an elongation
   that is a ratio of diameters — and returns it UNAVAILABLE, because `apply_scale` clips at the
   measured background and a low quantile of a population sitting there returns `0.0000`, which
   admits every measured cell while looking fitted.

   The donor's ordinary-tissue CD99 reference level — the p99 of cells in no classified object
   and outside every islet, the same estimator the hormone, Blood vessel, Nerve and Lymphatic
   levels use — is **carried alongside** as `donor_reference_level` and **REPORTED**. It is never
   substituted for the per-stratum level: it is there so a reader can see both numbers where a
   stratum could not be fitted, and so any fallback to it is a maintainer's choice rather than a
   silent one. Measured 2026-09-18: 115 `0.0518`, 6442 `0.0496`, 6534 `0.0610`.

4. **Elongation** as a second ROUTE that can only ADD: a cell CD99 does not reach can still be
   endocrine if its nucleus is as elongated as that stratum's confirmed endocrine cells. Large
   stratum ONLY — islets of STRICTLY more than 10 admitted cells — and barred below it **twice**,
   masked by the stratum and with no level fitted there at all. It is `nucleus_elongation`, never
   a cell-boundary feature: cell circularity runs 0.514-0.523 across every class, a tessellation
   artefact. **Its level is not settled** — see "Elongation is not a level-able feature" below.
5. **Endocrine** -> the hormone argmax names Beta / Alpha / Delta. **CD99 is still not in the
   argmax** — it decides *whether*, never *which*; three variants of the naming step were tested
   against the calls on disk and only the three-hormone argmax fits. **At least one hormone must
   have been READ before anything may name a cell**: argmax over three NaNs filled to -1 returns
   index 0, "Beta", every time, and a shape is not a name.
6. **Not endocrine** -> the ordinary tissue chain, unchanged. That population — the non-endocrine
   cells inside the islet — is the analysis, not a remainder.

**Islet size is a property of a completed classification**, so `classify()` runs the chain twice:
pass A is the chain exactly as it stood before this change, which fixes the admitted set (islet
core, non-structural, not `Excluded`) and therefore the sizes; pass B is the islet chain. Pass A
doubles as the audit baseline, so "what moved" costs nothing extra, and `--no-islet-chain` returns
it. Strata are `1`, `2-3`, `4-10`, `>10`, with `>10` strictly greater — the same boundary
`islet_identity_mixture.py` already takes as `--min-islet-cells`.

### Why CD99 sits at a level and elongation at a shape  [measured, 20 donors, 2026-09-18]

`var_ratio` = spread(undecided) / spread(reference), where reference = unambiguous hormone and
undecided = max scaled hormone < 0.10, islet-core cells with structural classes removed. In
parentheses, how many of the 20 donors are above 1.

| feature | islet size 1 | 2-3 | 4-10 | >10 |
|---|---|---|---|---|
| nucleus_elongation | 1.519 (1/20) | 1.139 (5/20) | 1.657 (16/20) | **2.089 (20/20)** |
| nucleus_circularity | 0.988 (0/20) | 1.078 (3/20) | 1.329 (14/20) | 1.529 (20/20) |
| CD99 | 1.021 (1/20) | 0.483 (0/20) | 0.665 (4/20) | 0.460 (1/20) |
| CD99 reference/undecided MEDIAN ratio | 2.07x | 2.63x | 2.34x | 3.85x |

Read it as two different kinds of signal. Nucleus shape carries a **mixture** signature — the
undecided cells are WIDER than the reference, which is what two populations in one group look like —
and it is unanimous only above ten cells per islet and noise below four. CD99 carries **no** mixture
signature at any size (its undecided cells are narrower); what it carries is 2.1-3.9x **less of the
marker**, widest in the large islets. So CD99 is a **location discriminator** and belongs at a level;
elongation is a **shape discriminator** that is real only in the large stratum. The three-donor
proving run reproduces the direction: elongation var_ratio in the large stratum 1.832 / 2.058 /
1.758 against the cohort's 2.089.

Elongation is therefore barred below `>10` twice — masked by the stratum, and with no level fitted
there at all. `nucleus_elongation` is used rather than a cell-boundary feature because cell
circularity runs 0.514-0.523 across six donors and every class; real immune and real exocrine cells
are not the same shape, so that is a tessellation artefact, and nuclei are segmented from DAPI.

#### Elongation is not a level-able feature, and that is unresolved  [OPEN, 2026-09-18]

A level separates populations that differ in **location**. On CD99 the confirmed/undecided median
ratio is 2.07-3.85x and reaches 12.4x on 6534's large islets. On nucleus elongation the two medians
are **the same to within 5 %** — ratio 1.019 / 1.003 / 0.985 / 0.954 across the four strata. The
difference between the populations is in *spread*, not in location, so no threshold can separate
them, and the one currently fitted lands a few percent above `1.0`, the floor a max-over-min
diameter ratio cannot go below. On 6534 it calls 22,022 of 88,736 large-islet cells endocrine on
shape alone — more than CD99 itself moves, which inverts the ruling it implements.

`islet_identity_mixture.py` was written for exactly this and answers the right question: is the
undecided population one or two, tested against a null **sampled from the reference** rather than
simulated, with a disjoint control that must not reject.

**It has now been run** — 2026-09-18, the first time, 20 donors, output at
`data/development/typing/islet_identity/`. In the large stratum, of 20 donor rows: **10 displaced
subpopulation, 3 two populations, 1 one population, 5 UNFIT because the CONTROL rejected, 1 too
few cells** — 14 testable. Per feature, against a Bonferroni bar of 0.0083 over six features, the
large stratum puts CD99 at the strongest separation of anything tested and `nucleus_elongation`
well behind it, with `nucleus_solidity` firing nowhere.

**That does not close the question, and it must not be read as closing it.** The test asks one
threshold against a cohort-pooled null per stratum; the maintainer's standing point is that a
single stratification is not expected to hold across donors whose islets differ in shape by
disease. So the elongation route stays provisional, `n_by_elongation_only` is quoted separately
and never pooled into "endocrine", and **elongation's weight is open for the maintainer**.

### The levels, and where they cannot be fitted  [CHAIN NOW]

`ENDOCRINE_QUANTILE` is derived, not chosen again: a structure's level is the p99 of a population
that cannot contain it ("1 % of ordinary tissue reaches this"), and endocrine identity inside an
islet is that statement upside down, so the level is the **1st percentile of the confirmed endocrine
cells themselves**. Same 1 % convention, no second threshold to defend. The count guard
(`MIN_REFERENCE_CELLS = 500`) is `islet_identity`'s existing reference guard, not a new number.

Fitted CD99 levels on the three proving donors, with the confirmed-cell count behind each:

| donor | 1 | 2-3 | 4-10 | >10 |
|---|---|---|---|---|
| 115 | unavailable (n=1) | unavailable (40) | unavailable (384) | 0.0745 (8,927) |
| 6442 | unavailable (8) | unavailable (131) | 0.0324 (569) | 0.0508 (12,306) |
| 6534 | unavailable (43) | unavailable (335) | unavailable (1,134) | unavailable (12,174) |

**Where two strata both fit, they differ:** 6442's `4-10` level is 0.0324 against 0.0508 at `>10`,
1.57x. Stratifying is not cosmetic. But on these donors only one stratum per donor (two on 6442) can
supply 500 confirmed cells, so the ruling has practical effect in the large stratum and is
**unanswerable below it** — those cells fall back to the hormone question, and the fallback is
counted per stratum in the audit rather than hidden.

The per-donor spread is itself an argument for per-donor fitting: the CD99 reference/undecided
median ratio in the large stratum is 1.97x (115), 4.85x (6442) and 12.44x (6534) around a cohort
median of 3.85x.

### Three defects found while building this, and what each one does now

**A fitted level of 0.0000.** `apply_scale` clips at the measured background, so a marker at or
below background reads exactly 0.0. On 6534, 20.31 % of islet cells sit at that floor, which puts
the 1st percentile of its hormone-confirmed cells at 0.0000 in **all four strata**. `CD99 >= 0.0` is
true of every measured cell, so the first run of this chain called 91,831 of 6534's 100,612 islet
cells endocrine (up from 44,540) on a number that looked fitted and was not.
`subpanel_io.endocrine_level` now refuses a level at or below the feature's own floor — 0.0 for a
scaled marker, 1.0 for an elongation that is a ratio of diameters — and returns it as unavailable.
6534 now moves 0 cells and falls back to the hormone question, which is the honest answer.

**Naming a cell from three NaNs.** 25,472 of 6534's islet cells have no INS, GCG or SST value at
all. `np.argmax` over three NaNs filled to -1 returns index 0 — "Beta" — every time, and before the
guard elongation alone called 23,554 of them endocrine. `endocrine_decision` now requires that at
least one hormone was **read** before any route can fire: whatever calls a cell endocrine, the
argmax that names it has only those three channels, and a shape is not a name.

**An exclusion reason that is not true of the cells it is written on.** [ON DISK]
`exocrine_inside_islet_border` asserts that exocrine markers reached the cell. On 6534, 24,309 of
the 45,457 in-islet `Excluded` cells have no Ker8_18, E-cadherin or Pan-Cytokeratin value at all;
they reach Exocrine by fall-through with **no** evidence, not exocrine evidence. The count read as
"how often the islet annotation caught exocrine tissue" is inflated 2.15x on that donor (45,457
quoted against 21,148 real). Those rows now carry `no_marker_values`. **Only the reason moves, never
the class**, which is why the baseline pass still reproduces the table on disk at 100.0000 %.

**The root cause under all three is upstream of this lane and is not fixed here.**
`data/development/celltype/uniform/donor_id=6534/data.parquet` has 1,417,265 rows, of which 205,070
carry no value for **any** of the 57 lineage markers (155,666 carry no value in any of the 65
columns). Every islet number quoted for 6534 today rests on that.

### What the islet chain moves, on the three proving donors  [CHAIN NOW]

Against the identical baseline pass, so the difference is the chain and nothing else:

| donor | cells changed | immune -> endocrine | excluded -> endocrine | endocrine -> other |
|---|---|---|---|---|
| 115 | 3,290 | 468 | 2,794 | 28 |
| 6442 | 261 | 56 | 170 | 35 |
| 6534 | 0 | 0 | 0 | 0 |

`immune -> endocrine` is reported on its own line, never inside a confusion table: those are immune
cells inside islets being reclassified as Beta/Alpha/Delta, and the immune environment in and around
islets is the study question. The endocrine decision is asked before the immune definitions, so this
direction exists by construction and has to stay visible.

Endocrine routes (cd99 / elongation-only / hormone-fallback, and islet-core endocrine before -> after):

```
115    42,464 / 3,689 / 3,216     46,135 -> 49,369
6442   18,957 /   578 /   208     19,552 -> 19,743
6534        0 /     0 / 44,540    44,540 -> 44,540
```

**What elongation adds that CD99 would not.** In the large stratum, 115 has 49,259 cells, CD99 calls
42,464, and elongation rescues 3,689 of the 6,795 it does not — 54.3 %. On 6442: 19,507 cells, 18,148
called, 578 of 1,359 rescued — 42.5 %. The fitted elongation levels are 1.0890 / 1.0896 / 1.1182,
a few percent above the floor of 1.0 that a max/min diameter ratio cannot go below, so the criterion
is close to "the nucleus is not a perfect circle". That is a lot of weight on a near-degenerate
threshold, and it follows from the measurement rather than contradicting it: the undecided cells are
WIDER on elongation, so a low-tail cut catches only the round half of the mixture. If elongation is
to do real work, its level has to come from the mixture — which is what `islet_identity_mixture.py`
fits — and not from a tail quantile. **Open for the maintainer.**

**Elongation is implemented as a second route, not a requirement.** "A second criterion" reads two
ways. The additive reading is implemented: elongation can only ever ADD endocrine cells. The AND
reading is counted and not applied — cells CD99 calls that elongation would reject in the large
stratum: 447 (115), 172 (6442), 0 (6534), in the audit as `n_large_cd99_ok_but_round`. **Open for
the maintainer.**

**A stratum with no fittable CD99 level silences elongation too**, even where the shape level did
fit. The shape route is defined as "a cell the key marker does not reach can still be endocrine", so
with no key marker to ask there is no "did not reach". On 6534 this is the difference between 0 cells
moving and 84,079 of 88,669 large-islet cells being called endocrine on nucleus shape and nothing
else. Recorded here because it is a decision, not an implementation detail.

**CD99 is very nearly a superset of the old hormone gate.** Cells the hormone level kept and CD99
demotes: 28 on 115 (0 of them hormone-confirmed), 35 on 6442 (7 confirmed), 0 on 6534.

### The structure call, as the maintainer's rule states it  [CHAIN NOW, 2026-09-19]

Four steps, in order, none sufficient alone. The maintainer's words:

> FOR ANY CELL WITH ITS CENTROID EXACTLY WITHIN THE BOUNDS OF THE OBJECT, THE ENTIRE CELL GETS
> CHECKED AND IF THE MAJORITY OF THE CELL BOUNDARIES ARE WITHIN THE OBJECT IT IS THEN A CANDIDATE
> TO BE CLASSIFIED ACCORDING TO WHAT THE OBJECT IS. IF THE SAME MARKER THAT CREATED THE OBJECT IS
> SIMPLY HIGHER THAN ALL THE REST, IT GETS CLASSIFIED. IF NOT, IT GETS A "PASSTHROUGH" FEATURE
> INDICATING THE FRACTION OF THE OBJECT THAT OVERLAPS THE CELL AREA.

Two of the four were missing from the code until 2026-09-19:

- **(b) the majority of the cell's own area inside the object** did not exist at all. Centroid-inside
  was being treated as candidacy. Its input, `crossed_frac_*` (correction D), was also restricted to
  islet-core and peri-islet cells, so most of the tissue had no overlap measurement and the question
  could not be asked there; that producer now measures every admitted cell.
- **(d) the passthrough was a boolean.** It now carries the overlap FRACTION, NaN where unmeasured.

**A small number of structure calls is the expected result.** Most vessels, lymphatics and nerves
are thin processes running through and around cells rather than the cell bodies themselves, so most
of what an object contains is passthrough. Donor 6450, whole tissue: of cells a vessel touches at
all, 22 % are majority-inside; lymphatic 22 %; nerve 10 %. Do not loosen `MAJORITY_AREA` because a
class came out small.

**Step (c) reads `top_expressed_marker`** — the argmax of the aligned values BEFORE `apply_scale`.
It read `top_identity_marker` until 2026-09-19, which is the argmax of a per-marker-normalised value
and is close to the scaled head-to-head the approved plan records as already tested and rejected.

Measured on donor 6450, both steps in place, against the normalised reading:

| object | majority-inside | + clears level | called, normalised | called, the rule |
|---|---|---|---|---|
| Blood vessel | 42,423 | 23,478 | 21,595 | **23,181** |
| Nerve | 532 | 529 | 456 | **377** |
| Lymphatic | 35,313 | 12,330 | 1,616 | **39** |

Blood vessel is the rule working: marker withholding falls from 1,883 to 297. **Lymphatic does not
survive it** — 12,330 cell bodies clear the level and 39 are called. The expressed argmax on this
donor is led by Ker8_18 35.5 %, CD34 26.2 %, CD11b 18.9 %; raw intensity is dominated by a few bright
antibodies and Podoplanin is not one, so it loses on 12,291 of 12,330. That is a measurement about
the panel, reported and not worked around. OPEN, and the maintainer's.

### The top-marker requirement, isolated against `--no-top-marker`  [measured 2026-09-18, pre-(b)]

These numbers predate steps (b) and (d) and the move to the expressed argmax. They are what the
top-marker condition moved ON ITS OWN, kept because they are the only isolation of it.

From a `--no-top-marker` diff, so this is the tightening alone:

| donor | cells moved | of | inside an islet |
|---|---|---|---|
| 115 | 22,298 | 860,863 (2.59 %) | 1,246 |
| 6442 | 9,674 | 1,419,669 (0.68 %) | 669 |
| 6534 | 45,951 | 1,160,909 (3.96 %) | 2,028 |

Largest moves: Blood vessel -> Immune 7,675 and Blood vessel -> Exocrine 6,960 on 115; Blood vessel
-> Immune 4,117 on 6442; **Lymphatic -> Exocrine 37,102** on 6534. Inside islets: Blood vessel ->
Beta 615 / Alpha 430 / Delta 63 (115), Blood vessel -> Beta 600 (6442), Lymphatic -> Alpha 846 /
Excluded 533 / Delta 243 (6534).

Withheld per structure (inside the object and clearing its level, but the cell's top identity marker
is not one that drew the object), against cells still called:

```
115    Blood vessel 16,483 (1,142 in islet) vs 48,359 called
       Lymphatic     4,816 (   74)          vs  5,129
       Nerve         2,946 (   33)          vs    652
6442   Blood vessel  6,299 (  643)          vs 54,467
       Lymphatic     3,571 (   12)          vs  5,079
       Nerve           528 (    4)          vs  1,081
6534   Lymphatic    45,727 (1,789)          vs  4,923
       Nerve         1,993 (  157)          vs  8,300
       Blood vessel  1,597 (  155)          vs 59,570
```

On 6534, and on 115's Nerve, the tightening withholds **more** cells than it calls. That is the
tightening working on objects whose contents are mostly pass-through, but it is a large move and it
should be looked at before the cohort re-run rather than after.

## Outside an islet — the ordinary chain, unchanged  [ON DISK, still current]

99.52 % from the upper tier alone, with **zero false positives**:

1. structure context -> its markers clear its level -> that structure; else fall through.
2. the immune definitions from `celltype/rules.requirements(..., with_negatives=True)`, evaluated at
   `require_min = 0.10` on the scaled values. Note the rule table, not the prose: `Imm-macrophage`
   is CD68 **or** Iba1 (an `any_of`), not both; `Imm-neutrophil` is CD66 **and** (CD11b or MPO).
3. a supporting tier over `bottom_tier_immune` at a level of **0.8715**, cohort-wide — the same
   value recovers independently on 6442 (0.8715) and 6436 (0.8716).
4. else Exocrine.

**Recovered 2026-09-15**: the supporting tier is eligible for a cell whose immune evidence
**clears** 0.8715 and whose exocrine evidence **does not** — the same level on both sides, so the
tier never overrules a cell the exocrine markers speak for as loudly as the immune ones do.

All 35,332 cells the original table claims by this tier on 6436 carry exocrine evidence below
0.8715; none is above it. The two rejected alternatives are both measured: applying the tier to
every unclaimed cell wrongly claims 5,671 strongly exocrine cells on 6436 (exocrine max p50 0.98),
and restricting to cells with *no* exocrine evidence keeps only 19,142 of the 35,332.

**Still open, and small**: 1,223 cells on 6436 (0.08 %), 98.4 % of them driven by PD_1 on cells
whose exocrine evidence sits just under the level — p50 0.843, against 0.044 for the cells the
tier should claim. Closing it means choosing a threshold because it fits the answer, which is
fitting rather than recovery, so it is left open.

The structure question **is** asked inside islets, with the object evidence, and it always was.
The pass-through is live there  [ON DISK]:

```
donor  structure       in object   clears level   called   passed through
6442   Blood vessel        1,969          1,312    1,312        33.4 %
115    Blood vessel        6,108          1,659    1,659        72.8 %
6534   Lymphatic           6,206            108      100        98.3 %
```

## The refinement pass  [CHAIN NOW — RULED IN 2026-09-23, WIRED THROUGH `--refine`; ON DISK: NOT YET RUN; see REFINEMENT.md]

Recovered 2026-09-18 from the 2026-09-11 session and made maintainable as `refine.py`.

> **RULED IN 2026-09-23, WIRED THROUGH THE `--refine` JOIN.** The maintainer ruled that the pass
> "needs to be included in the cell typing for phenocycler data" and scoped it the same day to
> `refine.py --out data/development/typing/refine`, then `build_cells.py --include-marginal
> --refine data/development/typing/refine`: the pass's verdict is joined **beside**
> `cell_class`, and the islet chain above — `cell_class`, `is_endocrine`, `endocrine_route` — is
> not touched. **[ON DISK]** the join has not run until the DATA step: no cell carries
> `Unresolved (no lineage evidence)` and `typing/cells/`'s `refined_*` columns are empty, and
> every builder still reads `cell_class`.

It records a second verdict on a cell's class, beside the first rather than over it, and it is
**not** described in full here;
[`REFINEMENT.md`](REFINEMENT.md) carries the rule, the cohort-scale measurements and the open
decisions. What a reader of this file needs:

* **The rule.** `concordance(cell)` = share of its K=30 marker-space neighbours carrying the same
  label (`pynndescent`, `random_state=0`, self dropped, `LOW = 0.40`), measured on the matrix
  `to_Z` returns — standardised arcsinh, in SD of each marker, from a cohort-wide cofactor /
  centre / spread table. `to_Z` **refuses** a marker it has no parameters for and **refuses** a
  non-finite value rather than substituting one: a missing measurement is unavailable, and its
  cell is dropped and counted, never filled in. `bound_X` = 5th percentile of
  max(own markers of X) over the concordant members of X. A cell of X below `bound_X` is **demoted**:
  a demoted Exocrine cell becomes `Unresolved (no lineage evidence)`; anything else becomes `Immune`
  if its `immune_leaf` flag is true, else `Exocrine`.
* **Why.** Object membership and the islet prior are evidence of WHERE a cell is, not that it IS
  that thing when its own marker is absent.
* **The same reference-level ruling applies**: `bound_X` is fitted per donor **and** per size
  stratum, over `subpanel_io.SIZE_STRATA`, with no bound at all for a group of fewer than
  `MIN_CLEAR` = 50 concordant members — UNAVAILABLE, not borrowed from another stratum, not
  replaced by the donor's and not defaulted, because demotion needs positive evidence that a cell
  falls below its own class's clear members and an unfittable group supplies none. `MIN_CLEAR` is
  the recovered run's own guard and is deliberately **not** `subpanel_io`'s 500; the log reports
  `n_groups_under_500` as well, so a group that clears 50 but not 500 is visible rather than
  assumed adequate. Neighbours are drawn **within one donor**.
* **It is applied once**, and that is encoded rather than documented: there is no step or iteration
  argument, and re-refining already-refined labels raises.
* **It writes a sidecar keyed on `object_id`** (`donor_id=<d>/data.parquet` under `--out`, plus
  `refine_log.json`) and does not rewrite `typing/cells/`.
* **What it READS, so a reader can tell a live input from a refused one.** `typing/cells/` for the
  classification; the donor's aligned expression table; the frozen cohort scale; and
  `data/development/celltype/shifts.json`, from which `derive_markers()` **derives** the marker
  set — 54 markers, dropping CD38 and b_Catenin1 by name as panel variants and CD57 by its
  `refused_reagent_change` provenance, never by name. The one artifact it **REFUSES** is
  `data/development/typing/umap/refinement_log.json` — see the note below.
* **The 5th-percentile bound demotes about 5 % of every class by construction**, which is
  `REFINEMENT.md` section 6. The ruling of 2026-09-23 included the pass without revisiting that
  floor, so a demotion count is still not evidence about the cells it removed: read it against
  the floor every log reports.
* **`Unresolved (no lineage evidence)` is a new class value.** `build_cells.py`'s display map and
  the app's domain configuration were told on 2026-09-18. `trajectory_singlecell_features.py` was
  the last reader still taking its own-marker table from the 2026-09-11 `umap/refinement_log.json`;
  **corrected 2026-09-19** — it now loads `refine.OWN_MARKERS` by path, with no fallback, and
  importing `refine` installs the process-wide refusal of that orphan. Being told is still not
  being wired: no cell carries the value, because the pass has not been joined into a build.

### `immune_leaf` — the per-cell flag the routing turns on  [CHAIN NOW — in `refine.py` only]

The routing bullet above says "if its `immune_leaf` flag is true". That flag is a quantity of the
refinement and of nothing else: `subpanel.py` has no `immune_leaf`, no cell in `typing/cells/`
carries the column, and the classification on disk never asked the question.

**What it is.** A per-cell boolean, one per `object_id`: *does this cell satisfy AT LEAST ONE immune
leaf phenotype?* `refine.immune_leaf(donor, markers)` returns `object_id -> immune_leaf` and nothing
else; `refine_donor` joins it on `object_id`, never on position or row order, and a cell with no row
in the join is `False`.

**Which phenotypes.** `immune_leaves()` reads them off the rule table rather than listing them —
`rules.requirements(raw, with_negatives=True)`, then every row whose phenotype name starts `Imm-`,
`T-` or `Mac-` and survives into `requirements`. A table that gains a leaf gains it here. Derived
2026-09-19, that is **15 leaves**: `Imm-plasma-cell`, `Imm-neutrophil`, `Imm-NK-cell`,
`Imm-dendritic-cell`, `Imm-macrophage`, `Imm-B-cell`, `Imm-T-cell`, `T-regulatory`, `T-exhausted`,
`T-progenitor-memory`, `T-CD8-cytotoxic`, `T-CD8`, `T-CD4-helper`, `Mac-M1-like`, `Mac-M2-like`.

**What satisfying one means.** `score.meets_requirements(Xs, have, reqs, leaves, require_min)`, then
`.any(axis=1)`. A leaf is satisfied when every `all_of` marker sits at or above `require_min`, every
`any_of` group has at least one member there, **and every `none_of` marker sits strictly below it**.
The negatives are the reason `with_negatives=True` is passed: a cell bright in a forbidden marker
does not get the leaf. `require_min` is `run_step3.REQUIRE_MIN` = **0.10**, imported and not
restated, so it cannot drift from the value the classification's own immune tier uses — 0.10 of the
way above that marker's own measured background, which is a statement about background and not
about how common immune cells are.

**A marker that was not measured fails the leaf closed.** In `meets_requirements`, NaN is neither
hot nor cold: an unmeasured `all_of` marker makes the leaf unverifiable, and an unmeasured `none_of`
marker fails the cell rather than passing on a constraint nobody checked. A missing measurement is
unavailable, not negative.

**Scope.** Per donor, over `qc_analysis_eligible` cells only, against the **frozen** cohort scale
(`subpanel_io.load_scale`) — the scale these classes were called against, not a recomputed one.

**What it routes.** In `route`, only demoted **non-Exocrine** cells consult it: true → `Immune`,
false → `Exocrine`. A demoted Exocrine cell becomes `Unresolved (no lineage evidence)` whatever the
flag says, because the flag cannot supply a complement that is not there.

**The 5,131.**  [RECOVERED RUN, 2026-09-11 — 320,000-cell subsample, NOT the cohort]
`refine.py` and `REFINEMENT.md` both record that adding `immune_leaf` "moved 5,131 cells". It is
recomputable from the arithmetic the 2026-09-11 record carries, and was recomputed on 2026-09-19:

    step-1 input                        320,000 cells (the eight class counts sum to it exactly)
    demoted Exocrine                     14,191  -> Unresolved   (= the log's Unresolved n)
    demoted non-Exocrine                 15,060  -> routed by the flag
    Immune kept (26,433 - 3,193)         23,240
    Immune after (logged)                28,371
    => routed to Immune by immune_leaf    5,131
    Exocrine: 222,736 + (15,060 - 5,131) = 232,665  = the log's Exocrine n

So on that subsample the flag is what stood between 5,131 cells and Exocrine; without it every
demoted region cell went to Exocrine and none could fall back to Immune. It is a subsample figure
from the recovered run — it is **not** a cohort number, and no cohort run of the pass exists.

> **The record those figures come from is refused as an input.**
> `data/development/typing/umap/refinement_log.json` was written 2026-09-11 13:38, **has no
> producer** anywhere in this repository, carries a step 2 marked `"accepted": false`, and counts a
> 320,000-cell subsample against the cohort's 25,230,109. `refine.py` installs a process-wide guard
> that raises `SupersededRefinementLog` on any read of it. Nothing is deleted — it stays on disk as
> the only surviving record of that run and of the step-2 rejection. **The committed authority for
> the defining-marker table is `refine.OWN_MARKERS`**, which is identical to the table that file
> carries (checked 2026-09-19, both steps, all eight classes). Quote the record as history; never
> parameterise anything from it.

**The routing is a no-op for some of them.** `route` excludes nothing that was already `Immune`, so
an Immune cell demoted for carrying none of its own nine markers above its bound, but carrying a
leaf, is relabelled `Immune`. Those rows keep `demoted = True` with an unchanged `refined_class`, so
they are countable — 267 cells on 6539, 2,930 on 6450, 34,163 on 6505 (23 % of everything demoted
there). Whether such a cell should instead fall to `Unresolved` is an **OPEN** decision
(`REFINEMENT.md` section 7), not one this file settles.

**Prevalence varies 3x across donors** — 6.4 % of measured cells on 6450, 7.0 % on 6539, 20.9 % on
6505 — and so does the share of demoted region cells it rescues: 5.8 %, 37.9 %, 57.1 %. Read no
cohort-wide rescue rate off one donor.

### The defining-marker table, and how it is not the subpanels  [CHAIN NOW]

The refinement rests on a per-class table of **own defining markers** — `refine.OWN_MARKERS`. It is
verbatim the recovered run's: checked 2026-09-19, `refine.OWN_MARKERS` and the `own_markers` block
in `umap/refinement_log.json` are **equal**, key for key.

| class | own defining markers |
|---|---|
| Nerve | B3TUBB, CD56 |
| Lymphatic | Podoplanin |
| Blood vessel | **CD31**, Caveolin, CD34 |
| Beta / Alpha / Delta | INS / GCG / SST |
| Immune | CD3e, CD68, Iba1, CD20, CD79a, MPO, CD66, CD11c, CD209 |
| Exocrine | Ker8_18, E_cadherin, Pan_Cytokeratin |

`bound_X` is the 5th percentile of **max over that row** among the concordant members of X, and a
cell of X below it is demoted. It is a "does this cell carry any of its own evidence" test, not a
signature match: one marker reaching the bound keeps the cell.

**This is the third marker table in the lane, and the three answer three different questions.** Do
not substitute one for another:

| table | question |
|---|---|
| `subpanel.CLASSIFIER_MARKER` | what **drew the object** the cell sits in |
| `subpanel_io.SUBPANEL` | which markers **answer that context's question** once a cell is in it |
| `refine.OWN_MARKERS` | does the cell **carry its own evidence at all**, whatever put it in the class |

Two differences from `SUBPANEL` matter, and one of the two things commonly said about them is
wrong. Both were checked against the code on 2026-09-19:

* **Blood vessel gains CD31 — true.** `SUBPANEL["Blood vessel"]` is `("Caveolin", "CD34")` and
  deliberately **not** CD31: the objects were drawn by `Cavdil1200` / `CD34dil1430`, so CD31 never
  drew anything and cannot answer "is this object's cell what the object is made of".
  `OWN_MARKERS["Blood vessel"]` is `("CD31", "Caveolin", "CD34")`. The refinement is not asking what
  drew the object — it is asking whether the cell carries vessel evidence — and CD31 is vessel
  evidence whether or not a pixel classifier was built on it. The two tables disagree because the
  questions differ, not because one is stale.

* **Exocrine does NOT acquire defining markers it lacked — the markers are the same, the
  DIRECTION reverses.** `SUBPANEL["exocrine"]` is already
  `("Ker8_18", "E_cadherin", "Pan_Cytokeratin")` — set-identical to `OWN_MARKERS["Exocrine"]`,
  verified. What changes is how they are read:
  * **In the chain  [ON DISK]** those three are a **veto only**. At `subpanel.py:929-934` the
    supporting immune tier refuses a cell that `speaks_exocrine` (max over the three at or above
    `SUPPORTING_LEVEL = 0.8715`). The class itself is then assigned at `subpanel.py:939` as
    `cell_class[unclaimed] = "Exocrine"` with `resolved_by = "fallthrough"`. **No cell on disk is
    Exocrine because those markers are up**; it is Exocrine because nothing else claimed it.
  * **In the refinement  [CHAIN NOW]** the same three become a **positive keep-test**: an Exocrine
    cell whose max over them falls below `bound_Exocrine` is demoted to
    `Unresolved (no lineage evidence)`.

  That reversal is exactly why the Exocrine floor of `REFINEMENT.md` section 6 bites hardest here.
  Exocrine's members were never selected for showing these markers, so taking a 5th percentile of
  them removes about 5 % of a population that was never required to carry any — and Exocrine's
  measured excess over the floor is only 0.2-0.5 percentage points, meaning essentially the whole
  Exocrine demotion is the percentile relabelling a fixed share by construction. **OPEN, not
  decided.**

Two more properties of the table, so it is not misread:

* **CD99 is not a defining marker of anything.** The islet subpanel is `("INS","GCG","SST","CD99")`
  read together, but `OWN_MARKERS` splits the hormones one per class and CD99 appears in no row —
  consistent with the chain's own rule that **CD99 decides *whether*, never *which***. The
  refinement cannot demote a cell for want of CD99.
* **Immune has nine defining markers and no `SUBPANEL` entry at all.** There is no `"Immune"` key in
  `SUBPANEL`; the chain's immune tier comes from the rule table. The nine are the refinement's own,
  and they are not the 15 leaf phenotypes above — a cell can carry a leaf and still fall below
  `bound_Immune`, which is how the no-op demotion of section 7 arises.

### `scatter` — what step 1 was accepted on, and step 2 rejected on  [RECOVERED RUN, 2026-09-11]

A cell is **scattered** when its concordance is below `LOW = 0.40`: fewer than 40 % of its 30
marker-space neighbours carry its label. `scattered_pct_before` / `scattered_pct_after` in
`refine_log.json` are the share of measured cells below that line.

**Scatter is a diagnostic and was never the target.** From the recovered script's own docstring:

> Scatter is DETECTED in the embedding (that is where it is visible) but every cell is CHARACTERISED
> and acted on in marker space, so a rule is a statement about measurements and the tidier picture
> is a consequence, never the target.

`refine.py` builds **no embedding at all** — nothing downstream of it needs one. The only surviving
trace of the UMAP is the directory the recovered log sits in.

**Step 1 was accepted at 11.51 % → 10.44 %.** Sourcing, because the two numbers do not come from the
same place: the 2026-09-11 record carries step 1 as `total_scattered_pct: 10.44` with
`accepted: true`. It has **no step-0 entry**, so the 11.51 % pre-pass figure is not in it — that one
is the recovered run's own record, carried in `REFINEMENT.md` section 7. (Both are quoted here as
history; that file is refused as an input — see the note above.)

**Step 2 was run, rejected and rolled back at 10.32 %.** The record carries it with
`accepted: false` and the reason verbatim:

> Re-applying the step-1 rule to an already-cleaned population is a ratchet, not a new cause. The
> cells it removes carry 3+ SD of their own defining marker and lie above the Exocrine 95th
> percentile in 91-100% of cases (Nerve +3.00 SD 100%, Delta +3.52 100%, Alpha +3.31 100%, Beta
> +3.60 100%, Blood vessel +1.37 91%), for a 0.12 percentage-point gain in scatter (10.44% ->
> 10.32%). Only Lymphatic's removed cells were genuinely marginal (median Podoplanin +0.86 SD, 0%
> above Exocrine p95). Rolling back; the refinement converged at step 1.

Section 6's floor is the mechanism behind that record: a percentile bound cuts about 5 % of a class
every time it is applied, so iteration converges on nothing and keeps cutting. The module **encodes**
the single pass rather than documenting it — no `step`, `passes`, `iterations` or `prev` argument
exists on any public function, and `apply_once` / `refine_donor` raise `RefinementAlreadyApplied`
when handed labels already carrying `Unresolved (no lineage evidence)` or a frame already carrying
`refined_class` / `refine_pass`.

**Do not read 11.51 → 10.44 as a property the pass must reproduce.**  [measured 2026-09-18]
At cohort scale the total can rise: 4.90 → 4.87 % on 6539, but **4.05 → 4.68 %** on 6450 and
**5.05 → 5.69 %** on 6505. Per class on 6539 every class improves (Immune 29.9 → 21.3 %, Lymphatic
31.0 → 21.5 %, Blood vessel 10.6 → 6.1 %, Alpha 14.8 → 9.4 %) while Exocrine worsens (1.1 → 1.8 %)
because demoted cells from other classes land in it, and `Unresolved` is 26.6 % scattered by
construction — it is by definition the cells with the least evidence. A rising total is not a
failure of the rule, and a falling one is not a verification of it. **No result here is a
verification target.**

---

**Documenting the refinement is not applying it, and wiring it is not running it.** Everything in
these three subsections describes `refine.py` and the 2026-09-11 recovered run. Since 2026-09-23
the `build_cells.py --refine` join is the declared route into `typing/cells/`, but **[ON DISK]** it
has not run: no cell carries `Unresolved (no lineage evidence)`, the `refined_*` columns and
`immune_leaf` are empty on every cell, and no builder reads `refined_class`. What
is on disk is the 2026-09-19 islet-chain build: eight classes plus `Excluded`.

## The gate, 2026-09-15  [ON DISK]

`subpanel.py --gate --levels pooled` rebuilds with `--membership parent` and compares row-for-row
against the table on disk. Every donor: identical row count, no cell only on one side.

| donor | cell_class | context | cells wrong |
|---|---|---|---|
| 6442 | 99.9995 % | 100.0000 % | 7 of 1,419,669 |
| 115 | 99.9990 % | 100.0000 % | 9 of 860,863 |
| 6436 | 99.9224 % | 100.0000 % | 1,226 of 1,579,440 |

All four reference levels reproduce to four decimals. The whole residual is the PD_1 set above.

The gate is a **reconstruction** gate and only means that with pooled levels, parent membership and
no islet chain. `--gate` with per-donor levels or with the islet chain on is the correction's effect,
not a reproduction failure, and the tool says so before it prints.

## Two marker facts about the clean batch-1 cohort

Measured while checking what the panel supports there; recorded because they bear on how any
immune number from these donors should be read, and neither is derivable from the code.

**CD38 is present on all ten.** Every clean-cohort donor is a `panel-cd38-59ch` donor, so the
`Imm-plasma-cell` rule runs as the table writes it. The CD79a+/CD20- stand-in that batch 2 needs is
not in force here.

**CD57 is effectively absent.** Its scaled p99 runs 0.0035-0.0266 across the ten against a 0.10
positivity level — it never fires. So on this cohort `Imm-NK-cell` reduces to CD56 with Granzyme B.
That is an absence, not the b-Catenin1 reagent-change confound (which is a batch-2 property): within
batch 1 CD57 is *consistent*, just consistently at background.

## What is not claimed

**Never acinar, never ductal.** The 60-plex carries no amylase, trypsin, CPA1 or elastase, and no
CK19, CA19-9, SOX9, HNF1B, CFTR or MUC1. **Exocrine is one label**, and `Unresolved (no lineage
evidence)` is not a lineage — it is the absence of one.

A missing measurement is unavailable: not zero, not negative, and not a level of zero. A cell with
no nucleus cannot clear a shape level, and is out of every calculation that needs one. `object_id`
is carried verbatim throughout; no join here is by centroid, row order or filename.

## What this file used to say, and does not any more

Until 2026-09-18 this file ended with a section headed "Not revisited", which said that nothing
re-derives `cell_class`. Both halves of that have stopped being true, and the correction is the
point of this rewrite: `subpanel.py` re-derives the classification, and `refine.py` proposes a
refinement of it.

`build_cells.py` does not decide either, but it no longer "only maps a display name". Corrected
2026-09-18: it CARRIES — every column the producer writes reaches the output, asserted in
`tests/current/test_build_cells_carry.py` — and it reads `is_endocrine` off `resolved_by` rather
than re-deriving it from the class NAME, keeps `Excluded` as a class with its reason instead of
removing it, attaches the per-stratum level each cell was measured against from
`subpanel.py --audit-out`, and joins `refine.py`'s sidecar under `--refine`. It reads
`typing/subpanel_rebuilt/`, so what it carries is whatever the last run of the producer wrote —
today, the 2026-09-19 islet-chain build (written 09:41–09:46).
