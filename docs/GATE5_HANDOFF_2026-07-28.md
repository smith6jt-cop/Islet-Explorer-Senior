# Handoff — 2026-07-28: Gate 5 evidence built, review figures rebuilt three times, back at the manuscript

**Start here.** Supersedes the "open items" of
[`GATE4_HANDOFF_2026-07-27.md`](GATE4_HANDOFF_2026-07-27.md) and adds to
[`GATE5_ROUND_2026-07-27.md`](GATE5_ROUND_2026-07-27.md), which holds the measured findings in full.

**Nothing is committed.** Changes sit in the working tree of the parent repo and the
`Phenocycler_Analysis` submodule. 252 tests pass. No process is running.

---

## The shape of this round, stated honestly

We built the Gate-5 evidence package, then spent the second half finding that **the instruments were
wrong, not the data**. Three defects in the review figures, each found only by the maintainer opening
finished PDFs, none of which touched a divisor or a call. Then the fourth question — "control and
positive doesn't make sense" — sent us back to the manuscript, which resolved it in one line.

That is the circle. It is worth naming what it actually taught, because it is not "the manuscript
method doesn't work out of the box":

> **Every defect we hit came from a place where our implementation had drifted from the manuscript's
> model and nothing checked the drift.** The categories multiplied because implementation bookkeeping
> (NNMF component membership, fit exclusions, double-highs) was given the same visual status as the
> manuscript's actual groups. The Step-2 line was drawn at the manuscript maximum after we froze the
> policy to p99 and never asked what the figures drew. The "control and positive" contradiction is
> p99's 1% showing through a display that still assumed the manuscript's `max`.

The manuscript is not the problem. The problem is that we changed one thing (the divisor statistic)
and left every consumer of it — the figures, the categories, the legend — assuming the old one.

---

## Where things stand

**Data on disk is sound and unchanged.** Gate 4's output (`data/restore_gated_redsea/`,
`data/phenotype/broad/`, 25,151,028 cells) was never touched this round. Gate 6's accounting is
**proven exact**: every canonical cell assigned exactly once, `Other` is precisely the
all-gates-negative population, zero `Unresolved`.

**The screen is trustworthy.** `screen_v12_production` (820 rows, 41 pairs × 20 donors) reproduces
`apply_manifest.csv` **bit-identically on all 480 accepted donor-pairs** — divisor, full maximum,
separator, both input floors, state, feature compartments, both policy fields. Asserted by
`scripts/gate5/check_screen_matches_production.py`, not assumed.

**The review bundle is not yet reviewable.** `expanded_review_v12` has the three display fixes but
still uses the rejected 6-category scheme.

---

## THE ONE DECISION BLOCKING EVERYTHING

Implement the manuscript's category scheme and rebuild the bundle once. The proposal, from
`literature/RESTORE.pdf` Methods:

The manuscript defines **three** cell groups for a mutually exclusive pair (i, j): `Y₁` reference-
positive, `Y₂` target-positive, `Y₃` negative for both. Lemma 1 operates on exactly those. **There is
no double-high group** — a double-positive cell contradicts the mutual-exclusivity premise the method
rests on. Our six colours were mostly implementation bookkeeping.

**Four colours:**

1. Negative control (reference-positive) — `Y₁`, sets the divisor
2. Target-positive, **called** (> divisor)
3. Target-positive, **not called** (≤ divisor)
4. Negative for both — `Y₃`

**Double-high stops being a colour** and becomes a header number — under the manuscript it is a
violation of the pair's premise, i.e. grounds to reject the pair, not something to look at in the
cloud. The raw-floor guide lines added on 2026-07-28 already mark that region without a colour.

**"Control and positive" disappears as a concept.** Measured: `control_fraction_le_1` = **0.9900 on
all 479 VALID rows, zero manuscript-conformant**. The manuscript's normalization makes negative
controls ≤ 1 *by definition*, which holds only when the divisor is the **maximum**. At **p99**, exactly
1% of every control sits above it by construction. Those 1,458 cells (6591 CD11b←Pan-CK) are that 1%,
not a biological category.

**The policy consequence is a maintainer decision, not a code fix:** p99 is *what breaks* the
manuscript's ≤1 guarantee. It was frozen deliberately because `max` cannot run at all (halts on donor
115 `Vimentin ← E_cadherin`, MODEL_UNSTABLE). Defensible — but the panel should state it rather than
encode it as a contradictory category.

---

## What was fixed this round (all in the display path)

| defect | cause | fix |
|---|---|---|
| Step-2 line drawn at the wrong threshold, y-axis stretched to it | both used `target_full_stats["maximum"]`; production uses `divisor` (p99). Median 6×, up to 200× apart. Line sat at y≈0.95 on every panel by construction; cloud crushed on **283 of 479** VALID rows | both now use `["divisor"]` — `restore_validation.py` ~2561, `restore_pair_review.py` ~1214 |
| Marker crops rendered background as signal | black point was the 1st percentile of *positive* pixels (**below** background, so ~99% of pixels non-black) and gamma 0.65 **lifted** the dim end | black point = that donor-marker's **raw input floor**; upper = 99.8th pct of pixels *above* it; **linear, gamma 1.0**. No-signal crops render black. `_display_floors` **raises** on a missing floor |
| Figure legend clipped off both page edges | ~150-char labels, `loc="center"`, which matplotlib does not shrink | hard-wrapped + measured fit check |
| Control-arm contamination invisible | categories are exclusive by overwrite, so a control cell that is also double-high paints only as a control | raw input floors drawn as guide lines; contamination stated in the header |

**Guards added** (4 tests, 247 total) so none of these can return silently: the Step-2 line and axis must
use `divisor` not `maximum`; the black point must equal the measured floor with every sub-floor pixel
exactly black; a no-signal crop must render black; a missing floor must raise. The legacy window is
pinned in a test *as the defect*, including the gamma lift.

**`scripts/gate5/preview_review_row.py`** renders ONE row through the production path in ~1 minute.
Use it before any 34-pair rebuild. Both display defects shipped in full bundles because there was no
way to look at one figure cheaply.

---

## Scientific findings (full detail in `GATE5_ROUND_2026-07-27.md`)

1. **The islet-core Immune disease gradient is aliased with acquisition BATCH.** Kruskal p=0.0016 by
   batch vs 0.0094 by status. Batch 2 (6450 6505 6523 6534 6566 6591 6623) holds every high-immune
   donor and **zero** ND donors; within batch 1 the T1D donors sit inside the ND range. The Fusion
   version was a perfect proxy and is **not** the mechanism (maintainer ruled it out). A cohort-design
   constraint on the study's own question — no disease claim about immune infiltration is safe without
   a batch-aware analysis. `data/donor_acquisition_batch.csv` now carries a `batch` column with a
   cross-check that raises if it drifts from `panel_variant`.
2. **CD11b's positives are 81.1% epithelial-marker⁺ and 80.8% epithelial⁺ outside islets; 78.6% carry
   no other immune marker.** Immune gates first, so all become Immune. **My earlier "under-called ~4×"
   claim is RETRACTED** — there is no ground truth for the pilot's 7–10%, and raising CD11b would add
   more of the same cells.
3. **CD34 agrees with CD31 everywhere CD31 fires (85.8%) but only 39.7% of CD34⁺ is CD31⁺.** CD34⁺CD31⁻
   is **10.5% of all cells**, 44.5% Vimentin⁺, only 1.0% in islet core — it avoids islets, where vessels
   are densest.
4. **Counterfactuals** (one gate edit, nothing else): dropping CD34 takes Endothelial 13.5→7.4% and
   nearly doubles Mesenchymal 2.0→3.9%, **and the islet-enrichment control survives** (8.6% core vs
   7.4% overall). Dropping CD11b takes Immune 21.3→19.7%.
5. **SST's "disease" pattern is withdrawn** — bimodal `target_fold` with nothing between 1.81 and 2.54;
   failures spread across all groups, Fisher p=0.167.
6. **INS is inert on its 4 weak donors** — dropping it changes 0–32 cells each, **zero** in islet core.
   Closed; no image review needed on that account.
7. **30% of all cells had ≥2 compartments eligible**, resolved silently by `COMPARTMENT_ORDER`.
   Mesenchymal loses **90%** of its eligible cells to higher gates, so "Mesenchymal = 2%" describes the
   ordering, not fibroblast abundance.
8. **The gold-standard scorer had gone stale** — it defaulted to a dead 7-marker anatomy-check arm and
   crashed. Repointed at production and dry-run end-to-end on synthetic labels.

---

## The pair-validity criterion that was missing, and the re-run (2026-07-28)

**Maintainer:** *"The RESTORE method should handle all cases."* Correct — and testing that split the two
"gate decisions" I had wrongly lumped together.

The manuscript's Lemma 1 partitions a pair's cells into three groups: reference-positive,
target-positive, negative for both. **A cell positive for BOTH contradicts the premise the method rests
on.** `double_high_fraction` measured exactly that on every screen and **was never used as a selection
criterion**. Ranking the 24 accepted pairs by it:

| | double-positive, median |
|---|---|
| `CD11b ← EpCAM` | **17.3%** (max 41.1%) |
| `EpCAM ← CD31` | 1.02% |
| the other 22 | **0.12% median** |

So `CD11b ← EpCAM` — frozen by a maintainer override on 2026-07-24 — is the single worst of all 41
candidate pairs, and its "negative control" is 17% target-positive. That is why its divisor was the only
one that rose, and why 81% of its calls are epithelial. **It is a RESTORE-level failure, not something
to patch in the lineage gate.**

`CD34 ← E_cadherin` is the opposite: **0.09%** double-positive, target_fold 16.9 — a clean pair, premise
satisfied, RESTORE answering its question correctly. CD34 is simply expressed by fibroblasts *and*
endothelium, which no intensity method can resolve. That one genuinely is a gate decision, because the
lineage gate is the only layer that can express it.

**Implemented** — `PairValidationConfig.max_double_positive_fraction` (0.05; an order-of-magnitude
separator, not tuned) and, in `reference_selection`: `double_positive_median`, `double_positive_max`,
`double_positive_donors_over_limit`, `passes_exclusivity_rule`. Exclusivity is a **hard veto applied
before the override branch** — an override expresses a preference among admissible references and cannot
make a non-exclusive pair valid. 5 new tests (**252 total**), including that an override *is* still
honoured when its reference is exclusive.

**Re-run: `data/restore_pair_validation/reference_selection_v3/`** (from `screen_v12_production`).

- **CD11b: all four epithelial references rejected as non-exclusive** (7.5–17.3% median; 80–90% of
  donors over the limit). **Only `CD11b ← CD3e` is exclusive** — 0.1% median, **0% of donors** over the
  limit. The EpCAM override is VOIDED and CD11b is now `open`, because CD3e fails the *frequency* rule
  (ratio 0.04). **The two criteria disagree, so this is a maintainer decision, not a computation** —
  and it is the decision that replaces "drop CD11b from the Immune gate".
- **A median-only rule is not enough.** Three pairs pass on the median and violate the premise on many
  individual donors: **`SST ← Pan_CK` 45% of donors over the limit** (max 46.1%), **`INS ← Pan_CK` 35%**
  (max 56.6%), `B3TUBB ← EpCAM` 5% (max 26.5%). Reported, not yet gating — whether the rule should be
  per-donor rather than median is the second open question.

## Next steps, in order

1. **Implement the 4-category scheme** (the decision above), then rebuild `expanded_review_v12` once.
   Preview one row first. ~50 min for the bundle.
2. **Maintainer reviews the bundle** — priority CD34, CD11b, MPO (MPO uniquely claims 2.6% of all cells
   and has had no image review) — and fills `data/gate5/gate5_signoff.csv`.
3. **CD11b's reference — a RESTORE decision, not a gate edit.** `reference_selection_v3` voids the
   EpCAM override and leaves CD11b `open`: CD3e is the only exclusive reference but fails the frequency
   rule. Three ways out, all maintainer calls: (a) accept CD3e and record that exclusivity outranks
   control size; (b) defer CD11b entirely, as CD20 was; (c) state a scope for the frequency rule the way
   `E_cadherin ← CD31` already did. **Whatever is chosen replaces "drop CD11b from the Immune gate",
   which would have hidden a broken pair behind a lineage patch.**
4. **CD34 out of the Endothelial gate — this one IS a gate decision.** Its pair is clean (0.09%
   double-positive); the marker genuinely marks fibroblasts and endothelium. Counterfactual already
   measured: Endothelial 13.5→7.4%, Mesenchymal 2.0→3.9%, islet-enrichment control survives.
5. **Decide whether the exclusivity rule is per-donor or median.** `SST ← Pan_CK` violates the premise
   on 45% of donors while passing on the median.

   Any of 3–5 requires re-running `restore_apply` + `lineage`, **which invalidates the Gate-5
   dashboard** — so settle them together, before anything downstream.
6. **Gold standard** — built, blinded, scorer proven, needs two annotators. Still the only thing that
   can verify any of this.
7. **Gate 7** (unstarted, blocked on Gate-6 acceptance): four broken downstream scripts —
   `build_phenotyped_h5ad.py` reads `broad_lineage`/`epi_default` (now `compartment`/`cell_type`);
   its input `singlecell_protein.h5ad` was purged; `build_cell_distance_parquet.py` hardcodes
   `CD20_pos`/`CD163_pos` which no longer exist; `compute_islet_markers.py` +
   `build_composition_and_tissue.py` use `X > 0.5` on raw MFI. Plus `data/app_data/` does not exist and
   `/srv/shiny-server/islet-explorer` still points at the **Junior** repo.

## Disk

`data/restore_pair_validation/` pruned 2026-07-28, **1,211 MB reclaimed**, now 6.4 GB (5.3 GB of that
is `qupath_compartment_full_20donor`, which is the irreplaceable input). See its README for what was
deleted and why. **Recommended next**, not done because it is outside what was asked:

- `data/_arm/` — **8.1 GB**. The anatomy-check isolated arms, built with the OLD 7-marker gate set. The
  gold-standard scorer no longer points at them, and every number the anatomy check relies on is in its
  tracked markdown. Strong delete candidate.
- `data/cells_preqc/` (9.0 GB) and `data/redsea_scratch/` (73 GB) are much larger but are provenance
  for the cell-QC and REDSEA stages — do not delete without a decision.

`/home` is at 80% (698 GB free).
