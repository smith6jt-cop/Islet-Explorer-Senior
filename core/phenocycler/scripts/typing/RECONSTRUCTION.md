# What the cell typing actually does

Two generations of calls are described here, and every statement says which one it is about,
because a reader in a month has to be able to tell them apart:

* **[ON DISK]** — about calls that already exist: `typing/subpanel/subpanel_calls.parquet`
  (25,150,537 rows, **superseded**) and `typing/cells/` (25,230,109 rows, **current**). Everything
  so marked was **recovered 2026-09-12 by reproducing those calls from their inputs** and checked
  against them, not inferred from prose.
* **[CHAIN NOW]** — the method as `subpanel.py` runs it after the maintainer's islet ruling of
  2026-09-18. Proven on three donors (115, 6442, 6534) with `--membership geometry --levels
  per-donor`; **NO DONOR HAS BEEN RE-RUN**, so nothing marked this way is yet a fact about any
  table on disk.

  Measured 2026-09-18, because this is what stale text keeps getting wrong.
  `typing/subpanel_rebuilt/` (25,230,109 rows, 20 donors, written 2026-09-15 12:41) and
  `typing/cells/` (25,230,109 rows, 20 donors, written 2026-09-15 13:16) are **the 09-15 build**:
  nine `cell_class` values — the eight classes plus `Excluded` — no `islet_cells`,
  `islet_stratum` or `islet_route` column at all, `resolved_by` carrying only
  `structure` / `definition` / `supporting` / `fallthrough` / `hormone`, and no `Unresolved`
  anywhere. There is no CD99 question in either table. Every `[CHAIN NOW]` number below comes
  from a three-donor proving run, not from a stage root.

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

## The refinement pass  [CHAIN NOW — see REFINEMENT.md]

Recovered 2026-09-18 from the 2026-09-11 session and made maintainable as `refine.py`.

> **WRITTEN AND NOT WIRED, as of 2026-09-18.** `refined_class` has **zero consumers**: nothing
> reads it, no cell on disk carries `Unresolved (no lineage evidence)`, and `typing/cells/` has
> no `refined_*` column at all. `build_cells.py --refine <dir>` is the join that would put it in
> front of the app, and it has not been run. The pass is therefore a proposal, not a second
> decision that has been taken.

It would be the second thing that decides a cell's class in this lane, and it is **not** described
in full here;
[`REFINEMENT.md`](REFINEMENT.md) carries the rule, the cohort-scale measurements and the open
decisions. What a reader of this file needs:

* **The rule.** `concordance(cell)` = share of its K=30 marker-space neighbours carrying the same
  label (`pynndescent`, `random_state=0`, self dropped, `LOW = 0.40`). `bound_X` = 5th percentile of
  max(own markers of X) over the concordant members of X. A cell of X below `bound_X` is **demoted**:
  a demoted Exocrine cell becomes `Unresolved (no lineage evidence)`; anything else becomes `Immune`
  if it satisfies an immune leaf combination, else `Exocrine`.
* **Why.** Object membership and the islet prior are evidence of WHERE a cell is, not that it IS
  that thing when its own marker is absent.
* **The same reference-level ruling applies**: `bound_X` is fitted per donor **and** per size
  stratum, over `subpanel_io.SIZE_STRATA`, with no bound at all for a group of fewer than 50
  concordant members. Neighbours are drawn **within one donor**.
* **It is applied once**, and that is encoded rather than documented: there is no step or iteration
  argument, and re-refining already-refined labels raises.
* **It writes a sidecar keyed on `object_id`** and does not rewrite `typing/cells/`.
* **The 5th-percentile bound demotes about 5 % of every class by construction**, which is
  `REFINEMENT.md` section 6 and is an **OPEN question, not a decision**. Do not read a demotion
  count as evidence about the cells it removed until that is ruled on.
* **`Unresolved (no lineage evidence)` is a new class value.** `build_cells.py`'s display map and
  the app's domain configuration have been told (2026-09-18);
  `trajectory_singlecell_features.py` has **not** — it still reads the 2026-09-11
  `umap/refinement_log.json` for its own-marker table, not `refine.py`'s `refine_log.json`.

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
today, the 2026-09-15 build.
