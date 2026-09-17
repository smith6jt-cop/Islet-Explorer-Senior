# Islet membership and crossed fraction: what they now cover, and what stopped them

2026-09-17. Two layers of the chain that feeds the app were short of the cell-type call, which
covers every donor in `data/development/typing/cells/`. This records what the shortfall was — it
was measured, not guessed — what it is now, and the donor-by-donor numbers.

    layer                        before   after
    typing/islet_membership        13      20
    typing/crossed                 10      20

Both lanes resolve their donor list through `cohort.resolve(..., include_marginal=...)`. Nothing
here restates a donor list or a cohort size; every number in this document came out of a run.

## Why membership stopped at 13, measured

Not a missing input. Islet outlines exist for every donor, and `islet_membership.membership()` was
already being called per cell for all of them from inside `build_cells.py` — the *artifact* under
`typing/islet_membership/` was what stopped at 13, not the capability.

What stopped it was a **stale input on the seven re-exported donors**. `build_donor` read two
things unconditionally:

* the content-addressed run's `expression/donor_id=<d>/` partition, and
* the `data/development/kintsugi/distance/` layer, written 2026-08-27.

The islet annotations for those seven were re-exported on 2026-09-14 and 2026-09-15. Both of the
inputs above predate that, so the signed distance was computed against annotations that are not
the annotations `outlines()` reads. The cross-check this lane runs every time — containment from
the outlines against the signed distance from another program — then fails loudly:

| 6450, the same donor, two ways | cells where the two disagree | median abs distance |
|---|---|---|
| run expression + kintsugi distance (2026-08-27) | 12,232 | 265.79 um |
| revised export + its own re-run signed distance | 1,131 | 0.42 um |

0.42 um is the regime every batch-1 donor sits in (0.24–0.46 um): the disagreements are cells
sitting on the outline, where a centroid inside by a fraction of a micrometre and a distance of
zero are the same fact rounded differently. 265.79 um is not a boundary effect. It is two
different sets of islets.

### What changed in the script

1. The cell geometry now comes from `subpanel_io.cell_source()` — the revised rebuild where one
   exists, the run otherwise — and the signed distance from whichever program computed it against
   the *same* annotations: the revised export's own `signed_dist_islet_um` column, or kintsugi for
   a run source. This is the pattern `build_cells.py` already used; this lane was the one that had
   not adopted it. Both sources are recorded per donor in `membership_diff.parquet`.
2. The outline file is now resolved by `build_islet_geometry.islet_export()` — recursive, newest
   mtime wins, a detection export refused with its reason printed — instead of a flat, non-recursive
   glob. The flat glob could not see a re-export that lands in a dated subdirectory, so this lane
   and `islet_geometry` could read two different sets of outlines for the same donor and nothing
   would have said so. On 6436 the refusal fires and is printed: a 66.3 GiB per-cell detection
   export sits in the islets directory under a name that matches an islet export exactly.
3. A donor that cannot be built is now named with its reason and the run continues to the rest,
   writing `not_built.json`; the exit status still carries the failure.

The cross-check is unchanged in kind, and it stays independent: containment is decided here from
the outlines, the distance by a different program entirely. Neither is assumed correct — where
they disagree, the run says by how much.

## Membership, per donor

All 20 donors, from `data/development/typing/islet_membership/membership_diff.parquet`. `in_islet`
and `in_expansion` are geometric and analysis-eligible; `in_islet by Parent` is QuPath's own answer
under the same eligibility filter, kept for the audit and never used to decide anything.

| donor | cell source | cells | eligible | in_islet | in_expansion | islets w/ a cell | islet outlines | in_islet by Parent | gained | lost | distance disagrees | median abs d um |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 115 | run | 1,030,803 | 860,863 | 53,454 | 72,824 | 1,304 | 1,304 | 53,185 | 965 | 696 | 2,616 | 0.31 |
| 6374 | run | 2,529,819 | 2,132,181 | 53,004 | 65,076 | 1,002 | 1,003 | 52,985 | 47 | 28 | 1,327 | 0.29 |
| 6380 | run | 737,647 | 623,781 | 1,320 | 3,851 | 67 | 67 | 1,294 | 69 | 43 | 165 | 0.28 |
| 6436 | run | 1,932,797 | 1,579,440 | 46,188 | 42,066 | 599 | 600 | 45,576 | 1,040 | 428 | 2,056 | 0.35 |
| 6442 | run | 1,715,283 | 1,419,669 | 20,626 | 29,122 | 487 | 490 | 20,421 | 542 | 337 | 1,366 | 0.38 |
| 6450 | revised | 438,658 | 369,728 | 17,917 | 20,895 | 448 | 453 | 17,771 | 369 | 223 | 1,131 | 0.42 |
| 6461 | run | 1,877,403 | 1,445,027 | 35,004 | 50,015 | 928 | 933 | 34,559 | 1,003 | 558 | 2,543 | 0.37 |
| 6476 | run | 1,979,670 | 1,625,342 | 52,741 | 66,708 | 1,269 | 1,275 | 52,491 | 1,025 | 775 | 2,943 | 0.38 |
| 6482 | run | 1,366,260 | 1,135,669 | 29,331 | 39,928 | 488 | 489 | 29,036 | 314 | 19 | 766 | 0.24 |
| 6493 | run | 1,545,699 | 1,227,411 | 41,903 | 41,011 | 695 | 696 | 41,721 | 644 | 462 | 1,923 | 0.39 |
| 6505 | revised | 2,870,418 | 2,196,442 | 124,542 | 109,276 | 3,081 | 3,123 | 122,823 | 5,098 | 3,379 | 11,536 | 0.65 |
| 6521 | run | 1,399,888 | 1,128,971 | 26,047 | 22,498 | 395 | 403 | 26,048 | 318 | 319 | 1,149 | 0.46 |
| 6523 | revised | 1,017,490 | 814,739 | 18,800 | 20,747 | 363 | 364 | 18,620 | 604 | 424 | 1,428 | 0.50 |
| 6534 | revised | 1,384,743 | 1,160,909 | 100,612 | 101,631 | 2,021 | 2,030 | 99,727 | 2,999 | 2,114 | 10,338 | 0.48 |
| 6539 | run | 501,991 | 364,201 | 10,066 | 15,109 | 514 | 572 | 6,624 | 3,497 | 55 | 531 | 0.31 |
| 6551 | run | 1,845,970 | 1,379,108 | 26,351 | 21,945 | 425 | 427 | 26,318 | 394 | 361 | 1,220 | 0.41 |
| 6566 | revised | 2,198,514 | 1,790,922 | 28,340 | 44,678 | 881 | 881 | 27,784 | 1,783 | 1,227 | 3,862 | 0.62 |
| 6591 | revised | 1,349,175 | 1,129,638 | 74,084 | 83,852 | 1,171 | 1,172 | 73,358 | 1,879 | 1,153 | 5,852 | 0.41 |
| 6593 | run | 2,046,673 | 1,629,937 | 33,472 | 27,591 | 418 | 420 | 33,269 | 637 | 434 | 1,584 | 0.42 |
| 6623 | revised | 1,563,091 | 1,216,131 | 54,901 | 56,907 | 1,687 | 1,702 | 54,792 | 1,689 | 1,580 | 6,961 | 0.49 |

Cohort: 31,331,992 cells read, **25,230,109 analysis-eligible**, 848,703 in an islet core,
935,730 in the 20 um expansion, 18,243 islets holding at least one eligible cell out of
**18,404** core outlines. Against QuPath's Parent: 24,916 cells gained, 14,615 lost.

Two coherence checks, both passing and neither a target:

* The eligible column sums to 25,230,109 and the outline column to 18,404 — independently the row
  count of `typing/cells/` and the core-islet count of `typing/app/islet_geometry.parquet`. Three
  lanes counting the same cohort from different files agree.
* Per donor, `in_islet` and `in_expansion_outline` never both hold, `islet_num` is set exactly when
  one of them does, and `object_id` is unique. Checked on 115 and 6505.

The 13 donors that were already built come back **identical**, value for value, so the source
routing changed nothing where nothing should change.

### The renumbering on the seven, and what it means

On every re-exported donor the `renumbered` count is nearly the whole islet population — 119,444
of 124,542 on 6505. That is expected and it is why `islet_num_parent` is kept: the re-export
re-derived the islet annotations inside the new tissue boundary, so QuPath's numbers are not the
old numbers. **An islet id from a batch-2 donor is only meaningful against the export it came
from.** Anything that joins an islet id across the re-export boundary is joining two different
numberings.

### Stratify before pooling: the fragment pattern is not only 6539

Cells per islet core, from `typing/cells/`, core region only:

| donor | core cells | islets | median cells/islet | islets under 10 cells |
|---|---|---|---|---|
| 6539 | 10,133 | 512 | 3 | 73.6 % |
| 6505 | 122,819 | 3,073 | 5 | 66.8 % |
| 6623 | 54,799 | 1,686 | 8 | 57.8 % |
| 6450 | 17,772 | 447 | 8 | 55.2 % |
| 6566 | 27,827 | 881 | 10 | 49.7 % |
| 6380 | 1,287 | 67 | 13 | 46.3 % |
| 6523 | 18,644 | 363 | 12 | 45.7 % |
| 6534 | 99,743 | 2,021 | 11 | 44.1 % |
| 6442 | 20,654 | 488 | 13 | 43.6 % |
| 6461 | 34,936 | 930 | 12 | 42.9 % |
| 6551 | 26,338 | 424 | 14 | 40.3 % |
| 6436 | 45,490 | 595 | 18 | 39.3 % |
| 6521 | 26,153 | 395 | 18 | 38.7 % |
| 6476 | 52,869 | 1,271 | 16 | 38.6 % |
| 6493 | 41,966 | 694 | 17 | 38.2 % |
| 115 | 53,439 | 1,304 | 14 | 37.2 % |
| 6593 | 33,378 | 418 | 23.5 | 36.6 % |
| 6374 | 53,182 | 1,002 | 23 | 27.1 % |
| 6591 | 73,350 | 1,171 | 21 | 26.6 % |
| 6482 | 29,292 | 489 | 29 | 17.0 % |

6539's classifier fragments are known. The measurement says 6505 sits beside it — median 5 cells
per islet core, two thirds of its 3,073 islets under ten cells — and 6623 and 6450 are not far
behind. A per-islet quantity pooled over the cohort is therefore weighted toward a handful of
donors' fragments. Stratify by islet size, as the endocrine-shape work already does above about
ten cells per islet, before any cohort-pooled per-islet number is reported.

## Crossed fraction

`crossed_frac_{vessel,lymphatic,nerve}` is how much of a cell's own segmented outline a structure
covers. It is a measurement beside the call, not a competing label: a cell can be immune **and**
crossed, and nothing here re-decides `cell_class`.

**Every donor's inputs are present — checked, not assumed.** For all 20 donors, exactly one
per-cell outline export in `data/redsea_scratch/geojson/` and exactly one object export for each
of vessel, lymphatic and nerve in `results/Geojson/`. So nothing about the data holds this lane at
10; what held it there was the donor list it was last run over, on 2026-09-14, before the
re-exported donors were admitted. Object ids survive the re-export: on 6450, of the 33,204
islet-associated cells wanted, **0 had no cell outline**.

### The ten files on disk are stale, measured

`typing/cells` was rebuilt on 2026-09-15 with geometric membership. The crossed files were written
on 2026-09-14. They therefore measure a different population than the one that is there now:

| donor | islet-associated cells today | in the 2026-09-14 file | missing | no longer islet-associated |
|---|---|---|---|---|
| 115 | 122,125 | 116,991 | 6,597 | 1,463 |
| 6374 | 115,818 | 107,879 | 7,995 | 56 |
| 6436 | 85,116 | 86,824 | 958 | 2,666 |
| 6442 | 49,143 | 48,822 | 1,034 | 713 |
| 6461 | 83,769 | 83,021 | 1,860 | 1,112 |
| 6476 | 116,807 | 98,379 | 19,959 | 1,531 |
| 6482 | 62,344 | 61,669 | 731 | 56 |
| 6493 | 81,033 | 80,049 | 1,836 | 852 |
| 6539 | 23,971 | 24,059 | 411 | 499 |
| 6593 | 59,476 | 59,606 | 768 | 898 |
| 6450 (built today) | 33,204 | 33,204 | 0 | 0 |

42,149 islet-associated cells are missing across the ten and 9,846 carried that no longer belong —
17 % missing on 6476. A count alone would not have shown this: on 6436, 6539 and 6593 the file
holds *more* rows than there are cells to measure.

Changes to `build_crossed_fraction.py`:

* A donor whose file is already present is skipped **only when its `object_id` set equals the
  current islet-associated set**; otherwise it is rebuilt and the run says what it found. A
  blanket skip would have preserved the ten stale files above. `--force` redoes a donor outright.
* Cells the outline export does not carry are counted per donor and reported. A cell with no
  outline has **no** crossed fraction — not a zero. Zero means a structure was looked for against
  that cell's own outline and none touched it.
* A donor whose cell outlines cannot be resolved is named and skipped rather than raising, so one
  donor's absence does not cost the others their measurement. `crossed_coverage.parquet` records
  the per-donor outcome.

### Crossed coverage, complete

All 20 donors, from `data/development/typing/crossed/crossed_coverage.parquet` and the per-cell
files themselves. `cells` is the donor's islet-associated population as `typing/cells` defines it
today; `no outline` counts cells the per-cell outline export does not carry, and those would have
had **no** crossed fraction rather than a zero.

| donor | cells | no outline | vessel reach % | vessel depth | lymphatic reach % | lymphatic depth | nerve reach % | nerve depth |
|---|---|---|---|---|---|---|---|---|
| 115 | 122,125 | 0 | 53.2 | 0.103 | 12.5 | 0.063 | 2.7 | 0.084 |
| 6374 | 115,818 | 0 | 40.9 | 0.090 | 4.9 | 0.059 | 2.0 | 0.081 |
| 6380 | 5,119 | 0 | 37.9 | 0.110 | 11.7 | 0.120 | 4.9 | 0.117 |
| 6436 | 85,116 | 0 | 25.2 | 0.069 | 3.2 | 0.059 | 0.3 | 0.086 |
| 6442 | 49,143 | 0 | 35.2 | 0.075 | 3.4 | 0.067 | 0.3 | 0.060 |
| 6450 | 33,204 | 0 | 55.3 | 0.107 | 47.0 | 0.190 | 2.3 | 0.085 |
| 6461 | 83,769 | 0 | 30.1 | 0.068 | 0.8 | 0.045 | 0.8 | 0.071 |
| 6476 | 116,807 | 0 | 37.6 | 0.097 | 35.5 | 0.113 | 5.0 | 0.134 |
| 6482 | 62,344 | 0 | 34.4 | 0.063 | 6.0 | 0.051 | 0.7 | 0.085 |
| 6493 | 81,033 | 0 | 29.6 | 0.066 | 0.7 | 0.064 | 2.1 | 0.072 |
| 6505 | 180,664 | 0 | 20.8 | 0.081 | 60.3 | 0.193 | 1.3 | 0.068 |
| 6521 | 45,491 | 0 | 38.7 | 0.151 | 11.0 | 0.069 | 6.1 | 0.153 |
| 6523 | 33,783 | 0 | 41.8 | 0.080 | 3.0 | 0.054 | 1.2 | 0.085 |
| 6534 | 155,097 | 0 | 47.1 | 0.110 | 65.1 | 0.271 | 6.7 | 0.076 |
| 6539 | 23,971 | 0 | 43.4 | 0.079 | 3.5 | 0.045 | 1.2 | 0.060 |
| 6551 | 45,784 | 0 | 35.0 | 0.084 | 23.2 | 0.088 | 2.2 | 0.122 |
| 6566 | 60,482 | 0 | 25.8 | 0.087 | 1.9 | 0.072 | 1.6 | 0.089 |
| 6591 | 126,751 | 0 | 56.2 | 0.146 | 8.8 | 0.083 | 3.2 | 0.113 |
| 6593 | 59,476 | 0 | 39.3 | 0.074 | 12.2 | 0.069 | 0.9 | 0.054 |
| 6623 | 84,293 | 0 | 53.0 | 0.108 | 62.8 | 0.222 | 3.9 | 0.077 |

**1,570,270 islet-associated cells measured, 0 with no cell outline**, against 800,503 rows in the
files that were on disk this morning. Reach is the share of cells any object of that structure
touches; depth is the median covered fraction over the touched cells. The two are different
questions and a structure that grazes half the cells is not one that engulfs them: vessels reach
20.8–56.2 % of cells and cover a median 6.3–15.1 % of a touched cell.

The completeness was checked independently of the builder's own coverage table, by re-deriving
each donor's islet-associated set from `typing/cells` and comparing `object_id` sets against the
file: **0 missing and 0 extra on all 20**, and 0 NaN crossed fractions. A count alone would not
have shown this — on three donors the stale file held *more* rows than there were cells to measure.

### The islet-level rollup, rebuilt on 20

`summarise_crossed_fraction.py --include-marginal --no-merge` rebuilt
`crossed_islet_summary.parquet`, which had been written on 2026-09-14 from ten donors:

    core  18,063 islets,  5,834 below the 5-cell floor (Unavailable)
    peri  18,374 islets,     51 below the 5-cell floor
    18,384 islets total

    peri vessel      median reach 0.447   median depth 0.228
    peri lymphatic   median reach 0.292   median depth 0.206
    peri nerve       median reach 0.000   median depth 0.135

A peri nerve median reach of zero is a measurement, not an absence: most islets have no nerve
object touching any peri cell, and the depth column is conditional on the ones that do.
`--no-merge` is deliberate — without it the rollup writes into `typing/app/parquet/islets.parquet`,
which belongs to the app build and not to this lane.

### Acquisition batch: eta-squared with a permutation null

Crossed reach is a cross-donor quantity, so the batch question has to be asked of it before any
cohort-pooled number is read as biology. Per-islet peri reach, 18,323 islets over 20 donors, 1,000
label permutations per cell of the table:

| quantity | by donor | by stain_group (k=9) | by acq_run (k=13) | null median (donor) |
|---|---|---|---|---|
| peri vessel reach | 0.380 | 0.144 | 0.158 | 0.0010 |
| peri lymphatic reach | 0.819 | 0.588 | 0.811 | 0.0010 |
| peri nerve reach | 0.167 | 0.091 | 0.132 | 0.0010 |

Every observed eta-squared sits far above its permutation null (p = 0.001, the floor at 1,000
permutations; null 95th percentile 0.0009–0.0017 throughout), so none of these is noise. The number
that matters is the third column against the second: for **lymphatic reach, acq_run alone accounts
for 0.811 of the variance against 0.819 for donor identity**. Almost all of the between-donor
spread in lymphatic reach is between acquisition runs, which is what the 0.7–65.1 % range in the
per-donor table looks like from the other side. Vessel and nerve are much less exposed, though not
free of it.

Consequence for anyone reading these columns: a cohort-pooled lymphatic crossed fraction is largely
reporting which run a donor was acquired in. Correct on `stain_group` or `acq_run`, never on the
Fusion version, or stratify and say so. This is a property of the Podoplanin objects the
classification was given, not of the crossed-fraction arithmetic, which is exact.

### What it took, and why it is no longer days

Every cell was intersected against whole structure objects, and a structure object is not small:
6505 carries 541,408 lymphatic objects whose median is 19 vertices — and one of **196,935**. Every
cell in the tissue overlaps that one, and `unary_union` re-nodes its input in full, so each of
6505's 180,664 islet-associated cells paid a full overlay against a 197k-vertex polygon. Measured
throughput at that point: **0.4 MB/s** through the donor's 5.87 GB outline export.

Two changes, both arithmetically identical to what was there:

* **Clip first.** `area(UNION_i G_i INTERSECT C)` equals `area(UNION_i (G_i INTERSECT C))`, so the
  union is taken over the fragments inside one cell rather than over the structure objects.
* **Query by shape.** `STRtree.query(..., predicate="intersects")` instead of by bounding box, so a
  candidate whose box overlaps but whose shape does not never reaches an overlay.

Verified on 6450 by rebuilding it and comparing against the file built the old way, cell by cell:
the same 33,204 `object_id` values, `cell_um2` bit-identical, and the three crossed fractions
agreeing to **1e-15** — floating-point round-off on a quantity in [0, 1], not a difference. With
that, the ten remaining donors and the ten stale ones all rebuilt in one evening.

A third change, `subdivide()`, cuts an oversized object into 4096-pixel tiles before any cell is
tested against it. On 6505 it does nothing: all 103 objects above 5,000 vertices fit inside a
single tile — the 197k-vertex object is convoluted, not sprawling. It is kept because it is cheap
and because the next donor's export need not have that shape, but **it is not what made the
difference here, and this document should not be read as saying it did.**

Both lanes are idempotent and both check staleness by `object_id`, so re-running either picks up
exactly what still needs work:

```bash
.venv/bin/python core/phenocycler/scripts/typing/build_crossed_fraction.py --include-marginal
.venv/bin/python core/phenocycler/scripts/typing/summarise_crossed_fraction.py \
    --include-marginal --no-merge
```

## Donors this lane cannot reach at all, with the reason

* **6457** — centrally excluded in `core/phenocycler/phenocycler/cohort.py`; poor DAPI prevented
  accurate channel registration. Not a coverage question.
* **6579** — reinstated centrally on 2026-09-17, and it has islet outlines and a cell-outline
  export. It has **no expression partition** in the run this lane reads and no signed-distance
  partition, and no `typing/cells` partition, so there are no cells here to give a membership to.
  It becomes reachable when the typing lane covers it, not before.
* **6533** — the 21st donor, carried in `core/phenocycler/data/expression/` rather than the run,
  and outside this lane's cohort resolution. Its markers are not shared across the cohort, so
  nothing derived from it crosses donors anyway.

## What this does and does not unblock

Membership and crossed both now cover every donor `typing/cells` covers, so neither is a cap any
more. They do **not** by themselves make the app show 20 donors: the app's composition artifacts
are built from `typing/immune/`, and that layer is the binding constraint. It is not this lane's.

Nor does 20-donor coverage make a cohort-pooled number safe. Two measurements above say to
stratify first: islet size, where 6505 sits beside 6539 with a median of 5 cells per islet core,
and acquisition run, which accounts for 0.811 of the variance in peri lymphatic reach.
