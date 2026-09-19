# The islet trajectory, rebuilt within donor — 2026-09-17

The maintainer decided on 2026-09-17 to refit the single-cell trajectory **within donor** for every
donor the inputs support, accepting that the ten donors already carrying a pseudotime would get
different values. This is what was built, what it cost, and what it is now safe to say.

Everything below was produced by running something. Where a number in the brief could not be
reproduced against the artifacts on disk, the measured value is given and the difference is named.

> **Read §6 before using the ordering.** The rebuild did what was asked — one fit per donor, twenty
> donors, all three Aab+ cases, the cross-donor relation reported separately under a permutation
> null. It also turned up something the pooled build had all along and nothing had measured: at
> default parameters on this data, **destiny → slingshot returns a different ordering on a different
> seed from identical input**, and only 1 of 20 donors reproduces its own axis. The counts, the
> coverage, the Aab+ answer and the batch relation stand. Per-islet *positions* on the axis do not
> yet, and that is measured rather than suspected.

---

## 1. What was there, and why it covered ten donors

| artifact | rows | donors | written |
|---|---|---|---|
| `typing/trajectory/singlecell_features.parquet` | 800,562 | 10 | 2026-09-14 10:43 |
| `typing/trajectory/singlecell_pseudotime.parquet` | 49,999 | 10 | 2026-09-14 11:13 |
| `typing/trajectory/singlecell_islet_summary.parquet` | 7,179 | 10 | 2026-09-14 11:14 |
| `typing/app/adata_ins_root.h5ad` | 7,302 islets | 10 | 2026-09-14 17:12 |

**The ten are not a choice anyone made about the trajectory.** `trajectory_singlecell_features.py`
resolved its donors through `COHORT.resolve(None)`, which returns `cohort.CLEAN` and nothing else.
On 2026-09-14 `CLEAN` was ten donors: the seven re-exported batch-2 donors had not been rebuilt yet
and so were not in it, and the three `MARGINAL` donors need `--include-marginal`. The seven joined
`CLEAN` on 2026-09-15, the day after the fit. The coverage is the default of a helper, a day stale.

**What it was fit on.** The producer already reads the current typing generation —
`data/development/typing/cells/donor_id=<d>/data.parquet`, not
`typing/subpanel/subpanel_calls.parquet`. So the generation was right; the *snapshot* was not.
`typing/cells/` was itself rewritten on 2026-09-15 13:12–13:16, after the features were built.
Measured against the current content of that tree, for the ten donors:

* **3,706 cells (0.46 %)** in the 2026-09-14 feature build are not eligible under the current
  `typing/cells/` — not 32,661 / 4.1 %.
* **36,367 cells** are eligible now and were not in that build. The larger number is the one that
  mattered: 6476 alone gained 18,987 eligible cells.

**How the 49,999 were sampled.** `N_FIT <- 100000` by default but the deployed run passed `50000`;
the R script took proportional per-donor quotas, `round(50000 * table(donor) / 800562)`, under seed
20260914, which rounds to 49,999. Per donor that ranged from 1,576 (6539) to 7,679 (115). Cells
outside the sample carried no pseudotime, and nothing was projected onto the fit. **6.2 % of the
population was ordered**, and that belongs in the caption, which it now is — the layer carries
`n_ordered_cells` and `n_population_cells` and the card prints them.

### The η² = 0.988 figure

It does not describe either artifact on disk, and it should stop being quoted as though it does.
Measured today:

| ordering | η² donor | η² stain_group | η² acq_run | η² donor_status |
|---|---|---|---|---|
| deployed single-cell slingshot (`singlecell_pseudotime.parquet`, 49,999 cells) | **0.0093** | 0.0064 | 0.0082 | 0.0064 |
| deployed islet DPT (`adata_ins_root.h5ad`, 4,371 ordered islets) | **0.5289** | 0.5201 | 0.5273 | 0.3563 |
| legacy-domain islet DPT (`data/app_data/adata_ins_root.h5ad`) | 0.0913 | — | — | 0.0445 |

0.988 is from the 2026-07 audit recorded in the header of
`scripts/senior/compute_trajectory_senior.py`: it described the **legacy CODEX-domain** scVI
trajectory *before* that script's Harmony fix, and the artifact it described no longer exists. The
current legacy-domain file measures 0.091 on donor.

**None of that weakens the decision, and the decision did not need it.** The alias it was invoked
for is present and measurable in this cohort, on the twenty donors that now enter:

| factor | Cramér's V against `donor_status` | levels |
|---|---|---|
| coarse batch (batch 1 / batch 2) | **0.604** | 2 |
| `stain_group` | **0.716** | 9 |
| `acq_run` | **0.882** | 13 |

Every ND donor is in batch 1. Batch 2 contains no ND donor at all — 2 Aab+ and 5 T1D. A pooled
embedding over this cohort has acquisition and disease status pointing the same way, and the islet
DPT at η² 0.529 on donor / 0.520 on stain_group is what that looks like in practice.

---

## 2. Coverage after the extension

`trajectory_singlecell_features.py` now defaults to **every eligible donor** — `CLEAN` plus
`MARGINAL`, with the centrally excluded donors removed by `phenocycler.cohort.filter_eligible_donors`
rather than by any list restated in this lane. `--clean-only` reproduces the old default.

**20 donors, 1,639,174 islet-core and peri-islet cells**, against 10 donors and 800,562 before.
Donor 6457 is excluded centrally and does not appear. Donor 6579 is reinstated as of 2026-09-17 but
has no partition under `typing/cells/`, so it is not in this build either — see
`docs/DONOR_6579_REENTRY_2026-09-17.md`.

| donor | status | stain | population cells | eligible | dropped | ordered | % | islets | islets ordered | islets with a core value | with a paired difference |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 115 | ND | S2 | 128,187 | 128,187 | 0 | 25,000 | 19.5 | 1,304 | 1,304 | 441 | 436 |
| 6374 | ND | S3 | 120,759 | 120,759 | 0 | 25,000 | 20.7 | 1,002 | 1,000 | 519 | 516 |
| 6442 | ND | S4 | 51,528 | 51,528 | 0 | 25,000 | 48.5 | 490 | 491 | 274 | 273 |
| 6461 | ND | S4 | 87,918 | 87,918 | 0 | 25,000 | 28.4 | 933 | 931 | 424 | 422 |
| 6476 | ND | S0 | 122,068 | 122,068 | 0 | 25,000 | 20.5 | 1,275 | 1,274 | 518 | 513 |
| 6482 | ND | S3 | 63,299 | 63,299 | 0 | 25,000 | 39.5 | 489 | 490 | 393 | 393 |
| 6493 | ND | S4 | 84,096 | 84,096 | 0 | 25,000 | 29.7 | 696 | 695 | 349 | 347 |
| 6539 | ND | S4 | 25,140 | 25,140 | 0 | 25,000 | 99.4 | 571 | 572 | 170 | 169 |
| **6450** | **Aab+** | S7 | 34,584 | 34,313 | 271 | 25,000 | 72.9 | 448 | 449 | 213 | 213 |
| **6505** | **Aab+** | S6 | 188,997 | 184,822 | 4,175 | 25,000 | 13.5 | 3,117 | 2,993 | 521 | 474 |
| **6521** | **Aab+** | S0 | 47,591 | 47,591 | 0 | 25,000 | 52.5 | 401 | 401 | 234 | 232 |
| 6380 | T1D | S1 | 5,444 | 5,444 | 0 | 5,444 | 100.0 | 67 | 68 | 56 | 56 |
| 6436 | T1D | S5 | 88,564 | 88,564 | 0 | 25,000 | 28.2 | 600 | 599 | 292 | 291 |
| 6523 | T1D | S8 | 35,086 | 34,196 | 890 | 25,000 | 73.1 | 364 | 365 | 223 | 222 |
| 6534 | T1D | S7 | 162,201 | 161,004 | 1,197 | 25,000 | 15.5 | 2,025 | 1,994 | 422 | 406 |
| 6551 | T1D | S1 | 47,643 | 47,643 | 0 | 25,000 | 52.5 | 427 | 428 | 246 | 244 |
| 6566 | T1D | S7 | 63,772 | 59,258 | 4,514 | 25,000 | 42.2 | 881 | 882 | 314 | 312 |
| 6591 | T1D | S7 | 132,097 | 132,051 | 46 | 25,000 | 18.9 | 1,172 | 1,169 | 414 | 411 |
| 6593 | T1D | S4 | 61,664 | 61,664 | 0 | 25,000 | 40.5 | 420 | 421 | 253 | 253 |
| 6623 | T1D | S6 | 88,536 | 88,389 | 147 | 25,000 | 28.3 | 1,702 | 1,702 | 352 | 348 |

**Every one of the twenty was ordered. No donor carries a "cannot be ordered" reason.** The floor
that would produce one — `MIN_CELLS_FIT = 500` eligible cells — is not reached by any donor; the
smallest, 6380, has 5,444 and was ordered whole.

**11,240 cells were dropped with a reason**, all of them on the seven re-exported batch-2 donors:
`n_dropped_no_aligned_marker_row`. Their rebuilt `object_id` has no row in the aligned values that
the rebuild reuses by join, so their 38 marker values are **unavailable** — not zero, not imputed —
and a cell with no values cannot enter a kernel that is a distance over those values. The count is
recorded per donor in `singlecell_trajectory.json`. 6566 is the worst affected at 4,514 of 63,772
(7.1 %), then 6505 at 4,175 of 188,997 (2.2 %).

### The sampling, stated in the artifact

`singlecell_trajectory.json` records `n_population` 1,639,174, `n_ordered` 480,444 and the
per-donor `n_ordered` / `n_eligible` / `sampled_fraction`. The layer carries `n_ordered_cells` and
`n_population_cells` and the card prints "Ordering was fit on 480,444 of 1,639,174 cells."

**Fitting within donor made the sample larger, not smaller.** The cost is quadratic in one donor's
cells, so twenty fits of N cost 20·N² where one pooled fit of 20·N costs 400·N². The build orders
**480,444 of 1,639,174 cells — 29.3 %**, against the pooled build's 49,999 of 800,562 — **6.2 %**,
in comparable wall time (≈3 min per donor, ≈55 min total).

---

## 3. Aab+ is reachable, and here is what enters

This was the sharpest gap in the deployed build: `donor_status` in `adata_ins_root.h5ad` has two
categories, ND and T1D. The tab could not show the middle of the progression the study question is
about, because **none of the ten donors is an Aab+ case**.

All three Aab+ donors in the cohort now enter, and none of them is marginal on cells:

| donor | stain | acq_run | coarse batch | ordered cells | islets | islets with a core value | why it was missing before |
|---|---|---|---|---|---|---|---|
| 6450 | S7 | 2026-06-18 | batch 2 | 25,000 | 448 | 213 | re-exported batch-2 donor; joined `CLEAN` on 2026-09-15, a day after the fit |
| 6505 | S6 | 2026-06-16 | batch 2 | 25,000 | 3,117 | 521 | same |
| 6521 | S0 | 2025-09-10 | batch 1 | 25,000 | 401 | 234 | `MARGINAL` (1–2 missed islets of 395); needs `--include-marginal`, which the default did not pass |

**Read the batch column before reading anything else into an Aab+ comparison.** Two of the three
Aab+ donors are in batch 2, which contains no ND donor; the third, 6521, is the only Aab+ donor in
batch 1 and is there at a stated cost. An ND↔Aab+↔T1D comparison drawn across these three is a
comparison across acquisition as much as across status — 6521 is the only donor that can carry the
Aab+ contrast inside batch 1, and one donor is not a variance decomposition. The cohort supports
*showing* Aab+ islets; it does not yet support an ND < Aab+ < T1D claim.

---

## 4. The fit, and the two things it costs

One `destiny::DiffusionMap` and one `slingshot` curve **per donor**, default parameters, no start
or leaf nodes, no cluster labels, 20 diffusion components — the published algorithm unchanged
(Steenbuck / Damond / Bodenmiller, bioRxiv 2025.03.05.641526). destiny 3.22.0, slingshot 2.16.0,
R 4.5.1. **No donor's cells enter another donor's kernel, and there is no pooled embedding
anywhere in the build.** Every donor returned exactly one lineage.

**Each donor's axis has its own scale.** slingshot pseudotime is arc length along that donor's own
curve; the units are not shared. So the column every consumer reads, `pseudotime`, is the
**within-donor quantile** — rank/(n+1) in (0,1). The raw arc length is kept as `pseudotime_raw`
and nothing downstream reads it.

**Each donor's axis has its own direction**, because slingshot given no start node picks an endpoint
arbitrarily. The direction is fixed by a declared convention applied identically to every donor:
flip if the donor's pooled core cells sit above its pooled peri cells, so low pseudotime is the core
end. **The cost is that the cohort-level *sign* of peri-minus-core is a convention, not a
measurement.** Ten of twenty donors needed the flip; the other ten already ran that way, and
`orientation_flipped` records which. What the convention does *not* force is the per-islet sign: it
is set on pooled cell medians, and islets inside a donor still go both ways — 6521 has 53.0 % of its
paired islets with a negative difference and a median of −0.024, 6591 49.6 %, 6461 46.4 %, while
6505, 6523, 6551 and 6593 have none at all.

### Consequence for the layer's batch diagnostic

`build_trajectory_layer.py` computed a within-status stain_group η² on the sampled-cell pseudotime.
Under a within-donor fit that number is **0 by arithmetic** — every donor is rank-uniform on (0,1),
so pooling donors inside a status gives every stain group the same mean. A zero that arithmetic
forces is not evidence that batch was removed, and printing it beside the word "batch" would read
as exactly that. The builder now reports it `unavailable` with the reason
`within_donor_fit_makes_pooled_batch_eta_sq_arithmetic(each_donor_is_rank_uniform_on_its_own_axis)`
whenever the source metadata records `fit_scope: within_donor`, and the real question is answered
in §5 instead. The 23 tests in `core/phenocycler/tests/current/test_trajectory_layer.py` still pass:
the synthetic fixtures carry no `fit_scope`, so their numbers are untouched.

---

## 5. Relating the donors to each other — a separate step, with a null

`core/phenocycler/scripts/typing/trajectory_crossdonor_relation.py`. The question is whether the
twenty orderings are the *same* ordering, described in a vocabulary the donors share:

> `profile[donor, marker] = Spearman(within-donor pseudotime, marker)` over that donor's ordered
> cells — 20 donors × 38 markers, no marker missing on any donor.

Pairwise profile correlation between donors: **median +0.128, range −0.901 to +0.952.** PC1 of the
profile matrix explains 0.621 of the variance, PC2 0.147.

**Every cross-donor number is reported against a permutation null over DONORS**, 10,000 draws,
because the donor is the unit that carries a batch label and because twenty donors spread over
thirteen acquisition runs give a large η² *on random labels*. p is the add-one estimate.

| quantity | factor | η² | null median | null p95 | p |
|---|---|---|---|---|---|
| profile PC1 | `acq_run` | 0.5365 | 0.6362 | 0.8592 | 0.742 |
| profile PC1 | `stain_group` | 0.3781 | 0.4171 | 0.6800 | 0.600 |
| profile PC1 | coarse batch | 0.0262 | 0.0263 | 0.2022 | 0.501 |
| profile PC1 | `donor_status` | 0.0913 | 0.0774 | 0.2923 | 0.439 |
| profile PC2 | `acq_run` | 0.4836 | 0.6432 | 0.8562 | 0.817 |
| profile PC2 | `stain_group` | 0.2926 | 0.4126 | 0.6879 | 0.769 |
| profile PC2 | coarse batch | 0.1920 | 0.0263 | 0.1960 | **0.052** |
| profile PC2 | `donor_status` | 0.1663 | 0.0790 | 0.2893 | 0.215 |
| median peri-minus-core | `acq_run` | 0.8128 | 0.6304 | 0.8603 | 0.103 |
| median peri-minus-core | `stain_group` | 0.4485 | 0.4145 | 0.6660 | 0.403 |
| median peri-minus-core | coarse batch | 0.1088 | 0.0243 | 0.1992 | 0.152 |
| median peri-minus-core | `donor_status` | 0.1235 | 0.0806 | 0.2962 | 0.324 |

**The cross-donor relation is not detectably a batch axis, and the null is the entire reason we can
say so.** The raw η² of 0.81 for the donor-median peri-minus-core on `acq_run` would have read as
an acquisition axis; the null median is 0.63, because thirteen groups over twenty donors explain
that much of anything. The one value near the line is profile PC2 on coarse batch, η² 0.192 against
a null p95 of 0.196, p = 0.052 — an excess that does not clear its own null and is reported as
such, not as a finding.

**This relation orders nothing and gates nothing.** It does not sort the donors, it is not offered
to the app, and no downstream quantity is conditioned on it.

---

## 6. What the rebuild cost the deployed ten

The maintainer accepted that the ten already-deployed donors would get different values. They got
**largely unrelated** ones, and that is the number rather than a characterisation of it.

`core/phenocycler/scripts/typing/trajectory_deployed_change.py`, comparing per islet — the level
the app keys on — for islets with a core value in both builds:

| donor | status | islets in both | Spearman (old vs new) | median &#124;Δ percentile&#124; | max &#124;Δ percentile&#124; | Spearman of peri-minus-core |
|---|---|---|---|---|---|---|
| 115 | ND | 164 | −0.351 | 0.364 | 0.950 | −0.226 |
| 6374 | ND | 183 | −0.295 | 0.297 | 0.975 | −0.018 |
| 6436 | T1D | 161 | −0.138 | 0.315 | 0.901 | −0.094 |
| 6442 | ND | 89 | −0.037 | 0.298 | 0.920 | −0.096 |
| 6461 | ND | 125 | −0.185 | 0.276 | 0.915 | +0.118 |
| 6476 | ND | 115 | −0.318 | 0.361 | 0.910 | −0.287 |
| 6482 | ND | 118 | −0.233 | 0.323 | 0.959 | +0.135 |
| 6493 | ND | 151 | +0.185 | 0.238 | 0.834 | +0.146 |
| 6539 | ND | 33 | +0.146 | 0.246 | 0.741 | +0.267 |
| 6593 | T1D | 125 | **+0.491** | **0.189** | 0.889 | +0.063 |

Both builds have an arbitrary direction, so the **magnitude** of ρ is the comparison and the sign is
the two conventions meeting. |ρ| at the islet level: **median 0.209, range 0.037–0.491**. At the
cell level, over the 1,298–1,615 cells each pair of samples happens to share: |ρ| median 0.137,
range 0.009–0.791.

**The yardstick is 1/3**: the expected absolute difference between two independent uniform
percentiles. The median |Δ percentile| across the ten is **0.298**, and nine of ten sit between
0.238 and 0.364. In plain terms, for most of the deployed ten **the new ordering is close to
unrelated to the old one**, not a shifted version of it. 6593 is the exception at ρ +0.491 and a
median shift of 0.189; 6493 and 6539 are weakly positive.

The superseded pooled artifacts are kept unmodified at
`data/development/typing/trajectory/superseded_pooled_2026-09-14/` so this comparison stays
re-runnable. Nothing reads them.

### How much of that change is the decision, and how much is the draw

This is the part the size of the change made necessary, and it is the one result here that was not
asked for. A replicate pass — the same fit, per donor, on an **independently drawn** subsample of
the same donor's cells — was run and compared, because "the ordering changed a lot" has two very
different causes and the artifact could not tell them apart.

`core/phenocycler/scripts/typing/trajectory_reproducibility.py`, against a replicate written by
`trajectory_singlecell.R 25000 99001` — the same fit, per donor, on an independently drawn
subsample. Three numbers per donor, because they fail separately: agreement over the cells the two
passes share, agreement of the two axes' 38-marker Spearman profiles, and agreement of the per-islet
core medians the app keys on.

| donor | status | overlap of the two draws | ρ shared cells | ρ marker profile | ρ islet median | sampled core cells / islet | verdict |
|---|---|---|---|---|---|---|---|
| 115 | ND | 0.198 | −0.387 | −0.366 | −0.185 | 3 | axis not reproducible |
| 6374 | ND | 0.206 | −0.042 | +0.137 | +0.279 | 6 | axis not reproducible |
| 6380 | T1D | **1.000** | −0.287 | −0.337 | +0.448 | 13 | **same cells, different axis** |
| 6436 | T1D | 0.283 | −0.168 | +0.081 | −0.354 | 7 | axis not reproducible |
| 6442 | ND | 0.483 | +0.025 | −0.011 | +0.018 | 7 | axis not reproducible |
| 6450 | Aab+ | 0.729 | +0.163 | +0.571 | −0.339 | 4 | axis not reproducible |
| 6461 | ND | 0.286 | +0.325 | +0.173 | +0.236 | 5 | axis not reproducible |
| 6476 | ND | 0.205 | +0.033 | +0.138 | +0.223 | 4 | axis not reproducible |
| 6482 | ND | 0.394 | +0.094 | +0.160 | +0.189 | 12 | axis not reproducible |
| 6493 | ND | 0.296 | −0.100 | +0.176 | −0.093 | 6 | axis not reproducible |
| 6505 | Aab+ | 0.136 | **+0.935** | **+0.990** | −0.034 | 2 | axis holds, islet median thin |
| 6521 | Aab+ | 0.526 | −0.064 | −0.468 | −0.553 | 11 | axis not reproducible |
| 6523 | T1D | 0.732 | −0.033 | −0.062 | +0.415 | 7 | axis not reproducible |
| 6534 | T1D | 0.156 | −0.063 | +0.085 | +0.334 | 2 | axis not reproducible |
| 6539 | ND | **0.994** | −0.371 | −0.731 | −0.002 | 2 | **same cells, different axis** |
| 6551 | T1D | 0.522 | −0.401 | −0.303 | −0.675 | 7 | axis not reproducible |
| 6566 | T1D | 0.426 | +0.124 | +0.344 | +0.280 | 3 | axis not reproducible |
| 6591 | T1D | 0.190 | +0.240 | +0.590 | +0.573 | 3 | axis not reproducible |
| 6593 | T1D | 0.408 | **+0.994** | **+1.000** | **+0.886** | 12 | reproducible |
| 6623 | T1D | 0.283 | +0.373 | +0.925 | +0.087 | 2 | axis holds, islet median thin |

**One donor of twenty reproduces.** Two more (6505, 6623) recover the same axis but cannot carry it
to the islet level because only 2 sampled core cells per islet survive the draw. The other
seventeen find a different axis.

### It is not the sampling, and it is not the within-donor decision

6380 was ordered whole — 5,444 of 5,444 eligible cells — so the two passes ran on **literally the
same matrix** and still disagreed at ρ −0.287. 6539 overlaps at 0.994 and disagrees at −0.371. That
rules the draw out as the cause, so it was tested directly on 6380:

```
same seed : rho(pseudotime) 1.0000   rho(DC1) +1.0000   eigenvalues identical
diff seed : rho(pseudotime) 0.0857   rho(DC1) −1.0000   eigenvalues identical
```

**The fit is exactly deterministic given a seed, and unstable across seeds on identical input.** The
diffusion map is not the unstable part: the eigenvalues agree to six decimals and DC1 comes back
exactly negated, which is a sign convention, not a different embedding. What changes is the
principal curve slingshot then fits through that 20-dimensional cloud — a different curve, giving a
pseudotime correlated at 0.086 with the first.

This is a property of destiny → slingshot at default parameters on this data. **The pooled build had
it too**; nothing measured it, and a single seeded run looks identical whether it is stable or not.

### What that does to §6

It very largely explains it. Comparing like with like on the same ten donors, at the islet level:

| comparison | median &#124;ρ&#124; | range |
|---|---|---|
| pooled (old) vs within-donor (new) | **0.209** | 0.037 – 0.491 |
| within-donor vs its own replicate, same method, different seed | **0.206** | 0.002 – 0.886 |

Paired Wilcoxon across the ten donors, p = 0.695. **The deployed ten changed no more than this
method changes when you reseed it.** The maintainer accepted a change; most of what they got was
not the within-donor decision but the seed.

### What this does and does not invalidate

*Unaffected.* The coverage (§2), the Aab+ answer (§3), the batch alias (§1) and the cross-donor
relation (§5) are counts and label statistics, not orderings. §5 in particular is measured on
marker-correlation profiles under a permutation null and its conclusion — no detectable batch axis —
does not turn on the axis being stable.

*Affected.* Any statement of the form "this islet is earlier than that one" — including
`within_donor_rank` and the sign and magnitude of `sc_pseudotime_peri_minus_core` for an individual
islet — is, for seventeen of twenty donors, a property of seed 20260917 + donor id rather than of
the donor. 6593 is the exception.

*Not fixed here, deliberately.* Stabilising it means changing the method — averaging over seeds,
fixing the eigenvector signs, supplying cluster labels or a start node, cutting `N_DCS` from 20 —
and every one of those is a method change made after seeing which way the results moved. The
instability is measured, recorded per donor in `reproducibility.json`, and handed over. It is the
first thing to settle before anything reads this ordering as a position in a process.


---

## 7. `sc_pseudotime_peri_minus_core`, and what may be said with it

This is the quantity the card leads on, and it is the right one: it is the difference between the
peri-islet compartment's position and the core's **within the same islet, on the same donor's own
axis**. It never needs a pooled embedding to mean something, and after this rebuild it never had
one.

**6,531 islets carry a measured paired difference**, over all 20 donors, against 1,020 over 10
before. `sc_n_core` and `sc_n_peri` travel beside every value; an islet with fewer than 5 sampled
cells in a region is Unavailable there with its reason, never zero.

Cohort-pooled per-islet numbers are wrong without stratifying, so both required strata were checked:

* **Islet size.** 4,273 of the 6,531 paired islets have ≥10 sampled core cells. Median difference
  +0.2270 above the floor against +0.2269 over all paired islets; the share of islets running
  negative is 21.3 % above the floor against 22.6 % overall. The difference does not depend on the
  size stratum here — which is itself worth knowing, since the endocrine *shape* signal does.
* **6539's fragments.** 6539 contributes 169 paired islets of 571 in its population; its islets are
  mostly one-cell classifier fragments and most of them never clear the 5-sampled-cell floor, which
  is the correct outcome rather than a gap.

Islet membership throughout is geometric — `cell_region_20um` from the signed distance to the islet
border — never QuPath's `Parent`, which mislabels 2.7 % of islet-core cells. The 20 µm peri-islet
expansion is derived from the signed distance, not from a drawn annotation. Cells without a nucleus
were excluded upstream in `build_cells.py` and never reach this lane. `object_id` is carried
verbatim from QuPath through every join.

**What must not be said.** Two limits, and they are different in kind.

The median is +0.227, but its sign is the orientation convention of §4, not a measurement. What that
convention leaves measured is the spread: how far apart the two compartments sit inside a given
islet, how much islets differ from each other inside a donor, and which islets run the other way.

The second limit is heavier. For seventeen of twenty donors the axis itself does not survive a
reseed (§6), so an **individual islet's** difference — its size, its sign, its rank against other
islets in the same donor — is not yet a stable quantity. The 6,531 paired islets are a real count of
where the difference is computable, and that is what this section establishes. Reading any one of
them as this islet's position in a process needs the instability settled first. 6593 is the one
donor where that reading is currently supported.

---

## 8. Files

Rebuilt or added, all under this lane:

| path | what |
|---|---|
| `core/phenocycler/scripts/typing/trajectory_singlecell_features.py` | defaults to every eligible donor; records per-donor coverage and source generation |
| `core/phenocycler/scripts/typing/trajectory_singlecell.R` | one diffusion map and one slingshot curve **per donor**; quantile output; declared orientation; replicate pass |
| `core/phenocycler/scripts/typing/trajectory_crossdonor_relation.py` | **new** — the cross-donor relation with a permutation null |
| `core/phenocycler/scripts/typing/trajectory_deployed_change.py` | **new** — what the rebuild cost the deployed ten |
| `core/phenocycler/scripts/typing/trajectory_reproducibility.py` | **new** — does the ordering reproduce itself, per donor |
| `scripts/senior/build_trajectory_layer.py` | refuses the pooled batch η² when the source fit is within-donor |
| `data/development/typing/trajectory/` | `singlecell_{features,pseudotime,islet_summary}.parquet`, `singlecell_pseudotime_replicate.parquet`, `singlecell_trajectory{,_replicate}.json`, `crossdonor_relation.json`, `crossdonor_axis_profiles.parquet`, `deployed_ten_change.json`, `reproducibility.json`, `superseded_pooled_2026-09-14/` |

To reproduce end to end:

```bash
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
.venv/bin/python core/phenocycler/scripts/typing/trajectory_singlecell_features.py
/opt/R/4.5.1/bin/Rscript core/phenocycler/scripts/typing/trajectory_singlecell.R 25000
.venv/bin/python core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py --no-merge
.venv/bin/python core/phenocycler/scripts/typing/trajectory_crossdonor_relation.py
.venv/bin/python core/phenocycler/scripts/typing/trajectory_deployed_change.py
/opt/R/4.5.1/bin/Rscript core/phenocycler/scripts/typing/trajectory_singlecell.R 25000 99001
.venv/bin/python core/phenocycler/scripts/typing/trajectory_reproducibility.py
ISLET_LAYERS_DIR=... .venv/bin/python scripts/senior/build_trajectory_layer.py
```

---

## 9. Left for the maintainer

> **CORRECTED 2026-09-18 — items 1, 2 and 3 were closed later the same day** by the app-side track,
> `docs/TRAJECTORY_APP_2026-09-17.md`. Read them as a record of what was outstanding at the time of
> writing, not as open work.
>
> * **1 (availability layer).** Rebuilt. `scripts/senior/build_availability.py` was run; `trajectory`
>   is `measured` for 20 of the 25 roster donors, and all three Aab+ donors are visible to the card.
> * **2 (`card_s2_trajectory.R`).** Edited. The false sentence — *"the fit was not recomputed per
>   donor"* — was removed, and the standing prose was cut to one measured line plus a methods modal.
> * **3 (`typing/app/adata_ins_root.h5ad`).** Rebuilt by `scripts/rebuild_adata_ins_root.py`:
>   17,786 islets, 20 donors, three ordered statuses, carrying the within-donor ordering. The
>   superseded ten-donor file is kept beside it as
>   `adata_ins_root.superseded-20260914-1712.h5ad`.
>
> **Items 4 and 5 are still open.** The ordering's seed-dependence (§6) is unsettled and is the
> first thing to settle before anything reads this ordering as a position in a process; and 6579
> still has no `typing/cells/` partition, with `data/donor_acquisition_runs.csv` still recording it
> `excluded=True` against `cohort.py` — verified on disk 2026-09-18.

Five things this session could not finish, most of them because the file belongs to a concurrent
lane this session must not write.

1. **The availability layer still hides the ten new donors, including all three Aab+ cases.**
   `/data/islet_layers/availability/data.parquet` records `trajectory: unavailable /
   no_islet_trajectory_exported` for 6380, 6450, 6505, 6521, 6523, 6534, 6551, 6566, 6591 and 6623,
   because it was built at 12:35 today against the ten-donor summary. The card checks that ledger
   before it draws, so **until it is rebuilt the Aab+ donors are in the data and not on the screen.**
   `scripts/senior/build_availability.py` reads
   `typing/trajectory/singlecell_islet_summary.parquet` directly and will flip all twenty to
   `measured`; it is outside this session's ownership, so it was not run. It is one command:

   ```bash
   .venv/bin/python scripts/senior/build_availability.py
   ```

   `/data/islet_layers/trajectory/data.parquet` **was** rebuilt — that one is this lane's own output,
   and leaving it stale left `test_built_artifacts_layer.py` failing. The previous version is in
   this session's scratchpad as `trajectory_layer_before_2026-09-17.parquet`.

2. **`apps/islet/R/card_s2_trajectory.R` now states something false.** Its caption ends with
   *"Ranks summarize the existing sampled-cell fit; the fit was not recomputed per donor."* The fit
   **was** recomputed per donor. The same caption prints the batch η², which the layer now returns
   as `unavailable` with a reason, so the card will take its `else` branch and print
   *"Within-status batch η² unavailable: within_donor_fit_makes_pooled_batch_eta_sq_arithmetic(…)"* —
   correct, but long, and it should probably say the shorter true thing and point at §5. The card is
   outside this session's ownership and was not edited.

3. **`typing/app/adata_ins_root.h5ad` is still the pooled ten-donor islet-composition DPT**, at η²
   0.529 on donor and 0.520 on `stain_group`, with `donor_status` still carrying only ND and T1D.
   It needs the same treatment the single-cell lane just got, and its producer
   (`trajectory_islet.py`) writes into `typing/app/`, which a concurrent lane owns. Until it is
   rebuilt, the Trajectory tab's H5AD-backed views and the single-cell card disagree about how many
   donors exist and whether Aab+ is in the data.

4. **The ordering's seed-dependence is the thing to settle next**, and settling it is a method
   change that should be made deliberately rather than by whoever next notices the numbers move.
   The candidates are visible from §6 — average the pseudotime over several seeds, pin the
   eigenvector signs before slingshot sees them, reduce `N_DCS` below 20, or supply the cluster
   labels the paper deliberately withheld. Each is a departure from the published algorithm and
   each will change every number in §7, so none was made here.

5. **Donor 6579** is reinstated but has no `typing/cells/` partition, so it is in neither build.
   `data/donor_acquisition_runs.csv` still records it `excluded=True`, which disagrees with
   `cohort.py`, and that file is read by this lane for `stain_group` and `acq_run`.

---

## 10. A note on source control

Commit `3373361` ("Every app artifact reaches all 20 donors…"), made by a concurrent lane at
18:27 while this work was in progress, swept four of this lane's in-progress files into itself:
`trajectory_singlecell_features.py`, `trajectory_singlecell.R` (at an intermediate state, before the
replicate pass), `trajectory_crossdonor_relation.py` and `trajectory_deployed_change.py`. It also
rewrote `trajectory_islet.py` (+244 lines) and `apps/islet/R/mod_trajectory_server.R`, both of which
were assigned to this track or flagged as coordinate-only. Nothing was lost — the working tree holds
the finished versions and this session committed nothing — but the history now attributes this
lane's files to that commit's message. Worth knowing before reading `git log` for provenance.
