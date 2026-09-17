# The app artifacts, rebuilt on the coverage that now exists — 2026-09-17

The app showed ten donors while the typing covered twenty. Nothing was broken: every path
resolved, every column existed, 772 R assertions passed. The six per-islet artifacts under
`data/development/typing/app/` had simply been built on 2026-09-14 against the donor list of that
day and never rebuilt, and no check compared one artifact's coverage against another's.

The binding constraint was `typing/immune` at ten donors. The membership, crossed and immune
tracks lifted their layers to twenty first; this rebuilds everything downstream of them.

## Verified on disk before rebuilding

The upstream tracks reported twenty. Checked rather than taken:

| lane | partitions | note |
|---|---|---|
| `typing/cells` | 20 | `islet_key`, `combined_islet_id`, `cell_region_20um` populated on all 20 |
| `typing/islet_membership` | 20 | |
| `typing/crossed` | 20 per-cell + `crossed_islet_summary` at 18,384 islets / 20 donors | |
| `typing/immune` | 20 | `coverage.json` `donors_built` = 20 |
| `celltype/uniform` | 20 | the aligned intensities the marker means average |
| `app/donor_metadata.csv` | 20 | |

## Before and after

Every file below is in `data/development/typing/app/`. "Donors" counts distinct `donor_id` /
`case_id` in the file itself.

| artifact | rows before | donors | rows after | donors |
|---|---|---|---|---|
| `islet_composition_long.parquet` | 111,823 | 10 | **255,753** | **20** |
| `islet_hormone_composition.parquet` | 7,634 | 10 | **18,063** | **20** |
| `islet_celltype_proportions.csv` | 7,780 | 10 | **18,384** | **20** |
| `islet_spatial_lookup.csv` | 7,782 | 10 | **18,404** | **20** |
| `islet_marker_means.parquet` | 7,781 | 10 | **18,384** | **20** |
| `neighborhood_metrics.csv` | 7,780 | 10 | **18,384** | **20** |
| `annotations.tsv` | 23,346 | 10 | **55,212** | **20** |
| `parquet/islets.parquet` | 7,780 × 104 | 10 | **18,384 × 104** | **20** |
| `parquet/uns_targets.parquet` | 93,384 | 10 | **220,848** | **20** |
| `parquet/uns_composition.parquet` | 7,634 | 10 | **18,063** | **20** |
| `master_results_run10.xlsx` (both sheets) | 7,634 / 93,384 | 10 | **18,063 / 220,848** | **20** |

`islet_geometry.parquet` was already 73,616 rows / 18,404 core islets / 20 donors and is
unchanged; it is the floor the rest now meet.

`islet_spatial_lookup.csv` carries 18,404 rows against 18,384 elsewhere on purpose: the lookup is
one row per core outline, and 20 outlines hold no analysis-eligible cell to compose. Similarly
18,063 of the 18,384 islets carry `prop_*`; the other 321 have an empty core after the
exocrine-inside-the-border set is dropped, and their proportions are **unavailable, not 0.0** — a
proportion over an empty region would state "0% of this islet's cells are Alpha" about an islet
with no cells to be 0% of, and drag every pooled median toward zero.

### What this changes about what the app can be asked

The ten donors it showed were all ND or T1D. The rebuilt table is the first that contains the
middle of the disease axis at all:

| | before | after |
|---|---|---|
| ND | 6,760 islets | 6,760 |
| Aab+ | **0 islets, 0 donors** | **3,966 islets, 3 donors** |
| T1D | 1,020 islets | 7,658 |

Confirmed through the app's own resolver: `domain_data_root()` under
`ISLET_DOMAIN_CONFIG=config/run20_default.yml` returns `../../data/development/typing/app`, and
`parquet/islets.parquet` read from there is 18,384 rows over 20 donors.

## The `--include-marginal` decision

**Used.** Every rebuild above ran `--include-marginal`, so the cohort is
`cohort.CLEAN` (17) **+** `cohort.MARGINAL` (3) = 20.

Not a default. A marginal donor is marginal for a measured reason, and the three it admits are,
per `core/phenocycler/scripts/typing/cohort.py`, each carrying its own measured defect:

| donor | measured defect | status |
|---|---|---|
| **6380** | 6 missed islets of 67 — **9.0%**, and **13% of its islet-core cells** affected | T1D |
| **6551** | 2–3 missed islets of 424 (0.5–0.7%) | T1D |
| **6521** | 1–2 missed islets of 395 (0.3–0.5%) | Aab+ |

6380 is the one to discount hardest: it is the smallest donor in the cohort (67 islets, against a
median of 648) and the defect touches an eighth of its islet-core cells, so any per-donor number
from 6380 is both the noisiest and the most biased in the set. 6551 and 6521 are sub-percent.

Taken because the layers underneath had already reached all twenty — `islet_geometry`,
`islet_membership`, `crossed` and `immune` all include these three — so building the app on 17
would have re-created the same mismatch one layer up, with the app outlining islets for donors it
could say nothing about. What excluding them would have cost: 897 core islets and one of the three
Aab+ donors. The other two Aab+ donors (6450, 6505) are in `CLEAN`, so the disease axis survives
either choice; the islets do not.

`6457` is centrally excluded in `phenocycler.cohort.DONOR_EXCLUSIONS` and appears in nothing.
`6579` (reinstated 2026-09-17) and `6533` have no partition in this lane, so they are absent
because there are no cells to build from — not because a list refused them.

## A wrong-islet join, found and fixed

`build_islet_marker_means.py` derived `islet_key` from the **aligned** per-cell table's own
`islet_num`. The batch-2 re-export renumbered the islets, and the aligned table was never rebuilt
against it. Measured 2026-09-17, on the cells the two tables share and both label:

| donor | islet numbering agreement, aligned vs rebuilt |
|---|---|
| 115, 6380, 6521 (batch 1, not re-exported) | 99.94%, 99.96%, 99.97% |
| **6505** | **0.15%** |
| **6450**, **6623** | **0.00%** |

So on the seven re-exported donors, every marker mean was filed under a **different islet** than
the one that key names in `islet_geometry`, `islet_celltype_proportions` and
`neighborhood_metrics` — a right-looking islet served the wrong cells' intensities, which is the
same failure mode the domain-root rule exists to prevent, one layer in. It also silently dropped
most of those donors' islets for want of a match: 6505 built 646 of its 3,117, 6450 built 70 of
448.

Fixed by taking islet membership from `typing/cells` joined on `object_id` — the geometric
membership every other app artifact already uses — and carrying `combined_islet_id` from that
table rather than re-spelling it. `object_id` is preserved verbatim; no centroid, row-order or
filename matching is involved.

Two consequences worth stating:

- **`Excluded` is now out of the marker means too.** `aggregate_islets.py` drops the
  exocrine-inside-an-islet-border set from the composition denominator; the marker means kept it.
  That is 2.7% of 115's islet cells but **21.9% of 6505's and 22.7% of 6534's**, so the two
  tables on one screen were answering the same question over different denominators, with those
  donors' means pulled toward exocrine. It also produced an islet the composition has no row for:
  `6539_Islet_536` holds exactly one eligible cell and that cell is `Excluded`. After the fix the
  two files agree on the islet set exactly — 18,384, with 0 on either side alone.
- **Batch-1 means moved as well**, and should have: 115 keeps 1,304 islets but only 532 have an
  unchanged `marker_mean_INS` (correlation 0.973), because the authoritative membership differs
  from the aligned table's by more than the numbering — it is the geometric assignment, which is
  also what corrects QuPath's Parent mislabel.
- `marker_mean_b_Catenin1` is new, and is not a regression: it exists only in donors that were
  outside the previous ten.

## Stratification: why no cohort-pooled per-islet number is quoted flat

Both stratifications are measured, from `islet_geometry`'s own `n_cells`:

**6539 is mostly one-cell classifier fragments — confirmed, and it is the extreme, not the rule.**

| donor | core outlines | median cells | single-cell | ≤3 cells | >10 cells |
|---|---|---|---|---|---|
| **6539** | 572 | **2** | **42.1%** | **62.4%** | 23.1% |
| 6505 | 3,123 | 5 | 8.1% | 35.8% | 31.6% |
| 6521 | 403 | 17 | 5.2% | 12.4% | 59.1% |
| 6623 | 1,702 | 7 | 4.1% | 20.6% | 38.1% |
| … | | | | | |
| 6482 | 489 | 29 | 0.2% | 1.4% | 80.8% |

Cohort: 18,404 core outlines, 4.1% single-cell, and only **50.9% above ten cells**.

**The endocrine-shape signal is real only above ~10 cells per islet**, so the >10-cell stratum is
the one a per-islet statistic belongs in. The difference is not cosmetic — on the per-islet core
immune proportion:

| stratum | n | mean | median |
|---|---|---|---|
| all core islets | 18,063 | 0.0474 | **0.0000** |
| >10 core cells | 9,368 | 0.0400 | 0.0049 |
| >10 cells, 6539 dropped | 9,236 | 0.0396 | 0.0049 |

The pooled median is **zero** — half the islets are too small to contain an immune cell at all, so
the unstratified figure measures islet size, not immune infiltration.

## Batch is aliased with disease status, and here it is larger

On the >10-cell stratum, per-islet core immune proportion, η² against a permutation null
(2,000 permutations of the labels, seed 20260917):

| factor | levels | η² | null mean | p |
|---|---|---|---|---|
| `donor_id` | 20 | 0.1881 | 0.0020 | < 0.0005 |
| `acq_run` | 13 | **0.1669** | 0.0013 | < 0.0005 |
| `stain_group` | 9 | **0.1112** | 0.0008 | < 0.0005 |
| `donor_status` | 3 | **0.0880** | 0.0002 | < 0.0005 |

The two finer acquisition variables each account for **more** variance than disease status does.
Any cross-donor claim off this table has to be corrected on `stain_group` or `acq_run`, or
stratified — never on the Fusion version, and never left uncorrected on the strength of the
disease term alone.

## Still ten donors, and named rather than filled

Three column families in `parquet/islets.parquet` come from producers outside this rebuild and
remain at the previous ten. They are `NA` — unavailable, not zero — for the other ten donors:

| family | cols | donors | islets with a value |
|---|---|---|---|
| `leiden_*` (islet-level Leiden + UMAP) | 6 | 10 | 4,371 / 18,384 |
| `sc_pseudotime*`, `sc_n_*` (single-cell ordering) | 7 | 10 | 1,269 / 18,384 |
| galectin `prop_Gal*` | 4 | 10 | 7,746 / 18,384 |

Unavailable on: 6380, 6450, 6505, 6521, 6523, 6534, 6551, 6566, 6591, 6623. Their sources
(`islet_leiden.parquet`, `typing/trajectory/singlecell_islet_summary.parquet`,
`islet_galectin.parquet`) are all still 10-donor files dated 2026-09-14. `crossed_*` (14 cols) and
`enrich_z_*` (8 cols) **are** at twenty.

`build_parquet_views.py` rewrites `islets.parquet` from scratch, which drops the 21 `crossed_*`
and `sc_*` columns that `summarise_crossed_fraction.py` and `summarise_singlecell_trajectory.py`
merge in afterwards. **Anyone re-running `build_parquet_views.py` must restore them, or the app
loses those columns silently** — the rebuild is not finished when that script exits.

Here they were restored by replaying only the merge block of each script — the same drop-stale,
left join on `(case_id, islet_key)`, `validate="one_to_one"` and row-count check those scripts
perform — rather than by running the scripts whole. `crossed_islet_summary.parquet` was already
correct at 18,384 islets over 20 donors, and re-deriving another track's artifact to obtain a
merge it had already earned would have risked leaving that file half-written for no gain. Only
`islets.parquet` was touched. Verified against the pre-rebuild schema: 104 columns before, 104
after, **0 lost and 0 gained**. Anyone reproducing from scratch should run the two scripts, as
below.

## The guard

`core/phenocycler/tests/current/test_app_artifacts_cover_the_cohort.py` (12 tests, marked `data`).
It compares each artifact's coverage to `islet_geometry`'s and to `cohort.resolve()`, **derives**
the expected donor set instead of naming it, refuses any centrally excluded donor, asserts the
four per-islet tables describe the same islet set, checks `combined_islet_id` agrees with
`case_id + islet_key` and resolves against the geometry, and asserts a proportion over an empty
region is `NA`.

Checked that it fails on the problem rather than merely passing now: run against the pre-rebuild
files, **7 of 7 artifacts fail**, each naming the ten missing donors. Against the rebuilt files,
12 of 12 pass.

## Test runs

- `pytest core/phenocycler/tests/current -q` — **586 passed**, 1 warning (a pre-existing
  `PerformanceWarning` in `test_grid_composition.py`).
- `testthat::test_dir("apps/islet/tests/testthat")` — **FAIL 0 | WARN 0 | SKIP 0 | PASS 772**.

### Tab smoke test

`NOT_CRAN=true ISLET_DOMAIN_CONFIG=config/run20_default.yml Rscript tests/tab_smoke.R /tmp/tabshots`

| tab | R errors | JS errors | drawn | verdict |
|---|---|---|---|---|
| Plot | 0 | 0 | 2 plotly drawn of 2 | ok |
| Trajectory | 0 | 0 | 3 plots, 2 plotly, all drawn | ok |
| Spatial | 0 | 0 | 2 plots, 2 plotly, all drawn | ok |
| Neighborhoods | 0 | 0 | 3 plotly drawn of 3 | ok |
| Statistics | 0 | 0 | 0 drawn; 2 tables present and **not drawn**, 10 empty html outputs | ok (see below) |
| CrossModal | 0 | 0 | 6 plots drawn of 6 | **FAIL** (see below) |
| LymphNode | 0 | 0 | 1 plot drawn of 1 | ok |
| Viewer | 0 | 0 | 0 | ok |

No R output errors and no JS console errors on any tab.

**Statistics** reports `tables 2 (drawn 0)` and 10 empty html outputs, and the harness still calls
it ok. That tab is entirely downstream of the Plot sidebar's current selection and shares its DOM,
so with no selection made its tables are empty by construction. Unchanged from before this
rebuild, and not evidence of a data problem — but it is the one tab this run does not positively
demonstrate drawing the new donors.

**CrossModal FAIL — not caused by this rebuild, and in a file outside this work's ownership.**
The harness expects 7 registry cards and saw 11: the 7 registered ones plus **4 elements carrying
the reserved `ie-card` class with no `data-card-id`**. Cause, located in the source:
`apps/islet/R/card_g1_grid.R:264` renders one `class = "ie-card"` div **per grid cell** — four
cells, pancreas/pLN × protein/RNA — inside the registered `g1_grid` card. The harness counts
those as unnamed cards and fails the tab. The counts line up exactly: of the 4 unnamed entries, 2
carry a drawn plot and 2 do not, matching the grid's measured and note-only cells. All 7 real
cards report `ok`.

This is the new-lanes track's file (untracked, along with the modification to
`mod_crossmodal_ui.R`), and none of the CrossModal cards read anything rebuilt here — g1_grid
reads `$ISLET_LAYERS_DIR/grid_composition`. **Needs the owner of `card_g1_grid.R`:** give the
inner grid cells their own class rather than reusing `ie-card`.

Two further pre-existing notes surfaced by the same run, both correctly named guards rather than
silent blanks: `x2_overlay` and `x3_identity` are gated off by `registration_gate_fail` (median
residual 7.20 µm ≥ 3.67 µm), and 5 screenshots failed with a Chromote
`Page.captureScreenshot` timeout — a capture failure, not a rendering one; those cards reported
their plots drawn.

## Reproducing

From the repo root, in dependency order. `build_workbook.py` reads what `aggregate_islets.py`
writes, and `build_parquet_views.py` reads the output of the two before it.

```bash
.venv/bin/python core/phenocycler/scripts/typing/aggregate_islets.py          --include-marginal
.venv/bin/python core/phenocycler/scripts/typing/build_workbook.py            --include-marginal
.venv/bin/python core/phenocycler/scripts/typing/neighborhood.py              --include-marginal
.venv/bin/python core/phenocycler/scripts/typing/build_islet_marker_means.py  --include-marginal
.venv/bin/python core/phenocycler/scripts/typing/build_parquet_views.py       --include-marginal
# then restore the two post-hoc merges build_parquet_views.py overwrites:
.venv/bin/python core/phenocycler/scripts/typing/summarise_crossed_fraction.py      --include-marginal
.venv/bin/python core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py
```

Nothing under `data/` was deleted; the pre-rebuild artifacts were copied to the session scratchpad
before the first builder ran.
