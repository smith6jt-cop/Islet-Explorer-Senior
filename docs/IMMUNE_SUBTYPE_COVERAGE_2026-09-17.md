# Immune subtyping: ten donors to twenty

Measured and written 2026-09-17. Producer: `core/phenocycler/scripts/typing/immune_subtype.py`.
Tests: `core/phenocycler/tests/current/test_immune_subtype.py` (24, all passing).
Artifacts: `data/development/typing/immune/` — one partition per donor plus `coverage.json`,
`per_donor_counts.csv`, `term_comparability.csv`, `unavailable_rules.json`.

```bash
.venv/bin/python core/phenocycler/scripts/typing/immune_subtype.py --include-marginal
```

## 1. Why it covered ten, per donor, measured

Not a missing input and not a marker the panels lack. **Every one of the twenty donors has both a
`typing/cells/` partition and a `celltype/uniform/` aligned partition. Thirteen carry all 37 markers
the fifteen immune rules name; the other seven carry 36, the one gap being CD38.** The ten were two
separate omissions:

| donors | why the run did not reach them | evidence |
|---|---|---|
| 6450, 6505, 6523, 6534, 6566, 6591, 6623 | the run predates their admission to the cohort | the immune partitions were written 2026-09-14 14:17, when `scripts/typing/cohort.CLEAN` was still the ten-donor tuple. Commit `4335647` (2026-09-15) added these seven after their rebuild from the revised QuPath export landed (`typing/revised_cells/`, written 2026-09-14 16:44 – 2026-09-15 12:12). `git log -p` on `cohort.py` shows the tuple before and after. |
| 6380, 6521, 6551 | `cohort.MARGINAL`, and the run was made without `--include-marginal` | `COHORT.resolve(None)` returns 17 donors today and `COHORT.resolve(None, include_marginal=True)` returns 20; the three are excluded from the default by their measured islet deficit, not by any immune input. |

Two further measurements make the point that nothing was missing:

* **The rule evaluation did not change.** Re-running the unchanged rules over donor 115 produced
  **0 disagreements on the 89,891 `object_id`s shared with the previous generation.** The row-count
  difference (92,270 → 93,357) is entirely the `typing/cells/` rebuild of 2026-09-15 13:12, which
  the 2026-09-14 immune run also predates. The old ten-donor output was stale against
  `typing/cells/` even for the ten it covered.
* **`unavailable_rules.json` was `{}` for all ten** — no rule had been refused for want of a
  marker, so nothing had been blocked by the panel.

## 2. Coverage now

**All 20 donors. 2,934,898 `Immune` cells over 25,230,109 cells — the full `typing/cells/` count,
to the row.**

| donor | acq_run | stain | cells | Immune | % | named | no rule fired | undecidable | markers | unavailable rule |
|---|---|---|---|---|---|---|---|---|---|---|
| 115 | 2025-09-16 | S2 | 860,863 | 93,357 | 10.8 | 90,088 | 3,269 | 0 | 57 | - |
| 6374 | 2026-04-01 | S3 | 2,132,181 | 159,383 | 7.5 | 154,435 | 4,948 | 0 | 57 | - |
| 6380 | 2025-09-13 | S1 | 623,781 | 102,384 | 16.4 | 101,393 | 991 | 0 | 57 | - |
| 6436 | 2026-04-07 | S5 | 1,579,440 | 224,244 | 14.2 | 200,174 | 24,070 | 0 | 57 | - |
| 6442 | 2026-04-05 | S4 | 1,419,669 | 225,150 | 15.9 | 218,800 | 6,350 | 0 | 57 | - |
| 6450 | 2026-06-18 | S7 | 369,728 | 23,850 | 6.5 | 21,597 | 1,982 | 271 | 57 | Imm-plasma-cell |
| 6461 | 2026-04-05 | S4 | 1,445,027 | 104,232 | 7.2 | 91,926 | 12,306 | 0 | 57 | - |
| 6476 | 2025-09-11 | S0 | 1,625,342 | 153,017 | 9.4 | 150,791 | 2,226 | 0 | 57 | - |
| 6482 | 2026-04-01 | S3 | 1,135,669 | 96,394 | 8.5 | 93,287 | 3,107 | 0 | 57 | - |
| 6493 | 2026-04-03 | S4 | 1,227,411 | 137,078 | 11.2 | 128,123 | 8,955 | 0 | 57 | - |
| 6505 | 2026-06-16 | S6 | 2,196,442 | 442,328 | 20.1 | 425,574 | 15,764 | 990 | 57 | Imm-plasma-cell |
| 6521 | 2025-09-10 | S0 | 1,128,971 | 80,859 | 7.2 | 77,834 | 3,025 | 0 | 57 | - |
| 6523 | 2026-06-21 | S8 | 814,739 | 52,823 | 6.5 | 50,508 | 2,276 | 39 | 56 | Imm-plasma-cell |
| 6534 | 2026-06-18 | S7 | 1,160,909 | 172,111 | 14.8 | 168,895 | 2,053 | 1,163 | 57 | Imm-plasma-cell |
| 6539 | 2026-04-03 | S4 | 364,201 | 37,497 | 10.3 | 27,426 | 10,071 | 0 | 57 | - |
| 6551 | 2025-09-13 | S1 | 1,379,108 | 103,639 | 7.5 | 98,112 | 5,527 | 0 | 57 | - |
| 6566 | 2026-06-19 | S7 | 1,790,922 | 307,715 | 17.2 | 301,043 | 6,591 | 81 | 56 | Imm-plasma-cell |
| 6591 | 2026-06-19 | S7 | 1,129,638 | 61,757 | 5.5 | 60,484 | 1,231 | 42 | 56 | Imm-plasma-cell |
| 6593 | 2026-04-04 | S4 | 1,629,937 | 290,894 | 17.8 | 257,634 | 33,260 | 0 | 57 | - |
| 6623 | 2026-06-16 | S6 | 1,216,131 | 66,186 | 5.4 | 61,250 | 4,822 | 114 | 57 | Imm-plasma-cell |
| **20 donors** | | | **25,230,109** | **2,934,898** | **11.6** | **2,779,374** | **152,824** | **2,700** | | |

### immune_class

| donor | Imm-plasma-cell | Imm-neutrophil | Imm-NK-cell | Imm-dendritic-cell | Imm-macrophage | Imm-B-cell | Imm-T-cell |
|---|---|---|---|---|---|---|---|
| 115 | 5743 | 706 | 864 | 1996 | 69484 | 9872 | 1423 |
| 6374 | 29202 | 15428 | 1580 | 3385 | 77560 | 22421 | 4859 |
| 6380 | 8150 | 3849 | 491 | 7122 | 64323 | 8774 | 8684 |
| 6436 | 17882 | 31108 | 766 | 10327 | 109666 | 24431 | 5994 |
| 6442 | 21587 | 6505 | 256 | 5813 | 151465 | 17419 | 15755 |
| 6450 | EMPTY | 116 | 815 | 894 | 14247 | 2503 | 3022 |
| 6461 | 27476 | 834 | 240 | 6763 | 37008 | 16936 | 2669 |
| 6476 | 14661 | 501 | 1985 | 3820 | 99459 | 25252 | 5113 |
| 6482 | 27244 | 3522 | 613 | 6623 | 40510 | 11432 | 3343 |
| 6493 | 23973 | 2714 | 558 | 4683 | 81555 | 9880 | 4760 |
| 6505 | EMPTY | 12883 | 19811 | 7831 | 328918 | 23632 | 32499 |
| 6521 | 18284 | 429 | 1364 | 7953 | 31514 | 13864 | 4426 |
| 6523 | EMPTY | 768 | 203 | 3455 | 33624 | 7922 | 4536 |
| 6534 | EMPTY | 2375 | 9387 | 4672 | 117197 | 17131 | 18133 |
| 6539 | 3927 | 80 | 234 | 2172 | 15042 | 4726 | 1245 |
| 6551 | 11715 | 131 | 545 | 7582 | 52372 | 13616 | 12151 |
| 6566 | EMPTY | 3052 | 596 | 32838 | 196154 | 37980 | 30423 |
| 6591 | EMPTY | 310 | 217 | 2440 | 39443 | 3319 | 14755 |
| 6593 | 16295 | 6553 | 643 | 13612 | 174823 | 20216 | 25492 |
| 6623 | EMPTY | 228 | 6298 | 2639 | 38853 | 9983 | 3249 |

### immune_subtype

| donor | T-regulatory | T-exhausted | T-progenitor-memory | T-CD8-cytotoxic | T-CD8 | T-CD4-helper | Mac-M1-like | Mac-M2-like |
|---|---|---|---|---|---|---|---|---|
| 115 | 36 | 288 | 256 | 209 | 104 | 280 | 188 | 39076 |
| 6374 | 20 | 2148 | 473 | 316 | 521 | 334 | 1155 | 11566 |
| 6380 | 365 | 2564 | 1339 | 504 | 847 | 1801 | 62 | 22867 |
| 6436 | 52 | 2460 | 881 | 596 | 1389 | 305 | 2007 | 19323 |
| 6442 | 254 | 3871 | 3068 | 1967 | 2216 | 2220 | 7028 | 5684 |
| 6450 | 9 | 1170 | 387 | 349 | 581 | 354 | 576 | 1834 |
| 6461 | 24 | 942 | 386 | 191 | 310 | 382 | 1421 | 9722 |
| 6476 | 107 | 1283 | 885 | 795 | 472 | 1003 | 470 | 50146 |
| 6482 | 39 | 1159 | 249 | 663 | 302 | 615 | 957 | 10947 |
| 6493 | 34 | 1641 | 503 | 537 | 755 | 824 | 2668 | 22232 |
| 6505 | 454 | 12164 | 5646 | 2194 | 3027 | 5120 | 6165 | 77766 |
| 6521 | 58 | 1440 | 493 | 406 | 712 | 470 | 601 | 14037 |
| 6523 | 25 | 707 | 725 | 643 | 849 | 473 | 204 | 11164 |
| 6534 | 64 | 6588 | 2976 | 1123 | 2599 | 2532 | 3080 | 12883 |
| 6539 | 5 | 268 | 201 | 324 | 235 | 127 | 890 | 3609 |
| 6551 | 41 | 4037 | 1647 | 1863 | 3034 | 700 | 830 | 32637 |
| 6566 | 189 | 6114 | 5144 | 2979 | 3344 | 6204 | 637 | 66426 |
| 6591 | 24 | 4008 | 2125 | 2169 | 3030 | 1387 | 486 | 15966 |
| 6593 | 373 | 8835 | 2622 | 3840 | 3662 | 2836 | 9328 | 39739 |
| 6623 | 5 | 1271 | 496 | 193 | 601 | 367 | 2123 | 10096 |

`EMPTY` is written, never `0`. `markers` is that donor's own acquired count, taken from that
donor's own aligned partition by two tests: the column is missing, or it is present and its
per-row-group null count equals the row count. A column exists for every marker in the cohort
union, so its presence is not evidence the slide carried it.

**On this cohort every gap is the first kind.** CD38 is column-absent on seven donors and
`b_Catenin1` on sixteen — the panel variant — and `panel_columns_entirely_null` is `[]` on all
twenty, so the all-null test fires on no donor here. It is carried and exercised in the suite
(`test_a_column_that_is_entirely_null_is_not_an_acquired_marker`) because it is the test that
catches the other shape of the same defect, not because it found anything in this build.
`b_Catenin1` is named by no immune rule, so only CD38 costs a rule.

Totals over the twenty: Imm-macrophage 1,773,217 · Imm-B-cell 301,309 · Imm-plasma-cell 226,139
(13 donors) · Imm-T-cell 202,531 · Imm-dendritic-cell 136,620 · Imm-neutrophil 92,092 ·
Imm-NK-cell 47,466 · uncalled-no-rule-fired 152,824 · uncalled-undecidable 2,700.
Subtypes: Mac-M2-like 477,720 · T-exhausted 62,958 · Mac-M1-like 40,876 · T-progenitor-memory
30,502 · T-CD8 28,590 · T-CD4-helper 28,334 · T-CD8-cytotoxic 21,861 · T-regulatory 2,178.

Internal coherence, asserted in the suite rather than claimed: named + undecidable +
no-rule-fired = the donor's Immune count on all twenty; no subtype row sits outside its own
parent (0 of 2,934,898); a subtype count never exceeds its parent.

## 3. What is EMPTY, and what "undecidable" means

**`Imm-plasma-cell` is EMPTY on 6450, 6505, 6523, 6534, 6566, 6591 and 6623.** The rule is
`CD38 + CD79a`, and those seven slides never carried CD38. An unasked question is not a negative
answer, so no number is written — not 0, and not a reduced rule that quietly answers a different
question. On the other thirteen the count is a measurement.

**2,700 cells are `undecidable`,** all on the seven donors rebuilt from the revised export
(6534 1,163 · 6505 990 · 6450 271 · 6623 114 · 6566 81 · 6591 42 · 6523 39). A rule naming a
marker that is NaN *on that cell* cannot be read on it. `value >= 0.10` is False for NaN, so this
population previously fell silently into "no rule fired" — a missing measurement reading as a
negative answer. It is now its own column (`immune_class_undecided`) and its own row in
`per_donor_counts.csv`, and it is never pooled with the 152,824 cells where every evaluable rule
*was* read and none fired.

## 4. The finding that changes how a count may be read: CD57

`Imm-NK-cell` is `CD56` **and** `(CD57 or Granzyme_B)` **and** one myeloid-group marker, with CD3e
negative. CD57 is one of the two corroborating terms, and it is not the same measurement across
the cohort:

| | CD57 above the level, among immune cells | raw CD57 median, immune cells | raw CD57 median, whole partition |
|---|---|---|---|
| 6450, 6505, 6534, 6623 | 86.03 %, 91.10 %, 95.95 %, 93.81 % | 435, 358, 430, 533 | 573, 496, 520, 625 |
| the other sixteen | 0.013 % – 1.543 %, and 0.013 % – 0.245 % on fifteen of them | 0.21 – 9.67 | 0.08 – 10.49 |

The two median columns are two different populations and are given separately because they were
measured separately: the split is the same in both, so it is the channel and not the immune
selection. 6551 is the one donor of the sixteen that reaches 1.543 %; the other fifteen sit at
0.013 – 0.245 %, and 6551 is still two orders of magnitude below the four.

Those four are **exactly `acq_run` 2026-06-16 and 2026-06-18** and nothing else — 6566 and 6591
share stain group S7 with 6450 and 6534 and read CD57 at background, so it is not the stain group
and it is not the panel variant. eta² by `acq_run` **0.998**, permutation p **0.0291** over 10,000
shuffles against a null median of 0.666; by `stain_group` eta² 0.691, p 0.065. The null median
matters: with 13 acquisition runs over 20 donors a shuffled label already reaches 0.67, so the
bare effect size would have meant nothing on its own. No correction is applied and none is
proposed here; the Fusion version is never the variable.

The consequence, written into the artifact: `Imm-NK-cell` on those four donors carries the status
`measured_not_cross_donor_comparable` with those numbers in its `reason`. Their counts
(6505 19,811 · 6534 9,387 · 6623 6,298 · 6450 815 against 203 – 1,985 on the sixteen) are a
measurement of that channel on that run and must not be pooled with the rest.

Four further terms cross the same descriptive criterion — `Mac-M1-like` / `Mac-M2-like` via
CD11b, CD66, MPO and HLA_DR, on one or two donors each and only just over it (0.504 – 0.665) — but
their permutation p runs 0.054 – 0.67 on `acq_run`, so nothing distinguishes them from chance
labelling. They are flagged with their
numbers and that is all the evidence supports. The full table of 152 (rule, marker) pairs is
`term_comparability.csv`.

**Reading across from the pLN lane:** it measured that CD8 never carries the NK gate (AUC
0.107 – 0.428 on all 32 clusters) and that CD56 is the term that works. This lane's rule table
already agrees without being changed — `Imm-NK-cell` requires CD56 and names CD8 nowhere at all.

## 5. T-regulatory: what the label names

FOXP3 is in `lineages.SUPPORTING_ONLY` — the maintainer's statement that it may not carry a call
alone — and `lineages.py` records that Nimbus scores it above its own boundary in 90.4 % of cells. The rule table licenses it
only in company: `T-regulatory` is `CD3e + CD4 + FOXP3` with `(CD39 or ICOS)`.

**So `T-regulatory` names the CD3e-positive, CD4-positive cells that also clear the level on CD39
or ICOS *and* on FOXP3. It is not a Treg count, and this panel carries no independent marker that
would make one.** The sentence is generated from the rule itself (`supporting_only_note`), so it
cannot drift from it, and it is written into the `supporting_only_note` column of
`per_donor_counts.csv` beside every T-regulatory row, into `coverage.json`, and onto stdout.

Two more of the same kind: `Imm-NK-cell` rests on CD56, also SUPPORTING_ONLY, and carries the
equivalent sentence. `T-CD4-helper` *excludes* cells above the level on FOXP3, and an exclusion
inherits the same uncertainty a requirement on it would.

`T-regulatory` is the smallest population in the build (2,178 cells over 20 donors, 5 – 454 per
donor). Nothing here is a reason to prefer a number that would make it larger.

## 6. The per-cell offset, carried over from the pLN lane

The pLN lane's measured point is that a dim cell reads low on every channel, so "cleared no rule"
and "is not that cell type" are different statements. `panel_offset` is written per cell: the
median, over that donor's acquired panel, of the cell's within-donor rank on each marker, taken
over the analysis-eligible cells only — a cell with no nucleus is not in `typing/cells/` and does
not enter the ranks. It is invariant to multiplying any channel by any positive constant.

**It is a diagnostic, not a gate.** No call anywhere in this build depends on it. What it shows is
that the uncalled population is not simply the dim one: the median offset of named cells runs
0.426 – 0.733 and of uncalled cells 0.472 – 0.758, and on 13 of 20 donors the *uncalled* cells sit
higher. On 6539 — whose islets are mostly one-cell classifier fragments — named 0.610 against
uncalled 0.646 on a 26.9 % uncalled rate, the highest in the cohort. Brightness does not explain
these cells.

## 7. What this does not unblock, stated plainly

**No app artifact reads this layer.** `build_workbook.py`, `aggregate_islets.py`,
`build_parquet_views.py` and `build_islet_marker_means.py` contain no reference to
`typing/immune`; their donor coverage comes from their own `COHORT.resolve` at build time, not
from a read of this layer. Extending immune subtyping to twenty does not by itself move a card.

The one place a reader is waiting is `neighborhood.py`, which says so in its own header: the named
immune sub-metrics `cd8_to_macro_ratio` and their siblings "are not written … because immune
subtyping is a later phase", and writing them would have put all-NaN columns behind live-looking
controls. That phase has now happened on all twenty donors, so those columns are derivable. That
file is outside this change.

## 8. Rules held

* Nothing was deleted under any `data/` directory. The previous ten partitions were overwritten in
  place by their twenty-donor successors; a copy of the superseded generation was taken to this
  session's scratchpad first, and it is not a durable archive.
* `object_id` is preserved verbatim and every join is `validate="one_to_one"`; 0 immune cells on
  any donor failed to match an aligned row.
* No donor id, cohort size or typing method appears as a string literal in the producer — asserted
  by `test_no_donor_id_is_written_into_this_module_as_a_string`, which reads the syntax tree.
* Donor exclusions are derived, never restated: 6457 stays out through
  `phenocycler.cohort.DONOR_EXCLUSIONS`; 6579 has no partition in this lane and is not claimed.
* `cell_class` is never modified and no cell outside `Immune` is read, scored or written.
* No count is a verification target, and no rule was adjusted toward a number.
