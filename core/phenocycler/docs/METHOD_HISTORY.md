# Method history

> **Superseded current-state notice (2026-08-13):** this file records historical
> method decisions; it does not describe the active workflow. There are now no
> verified production target/reference pairs. Ordinary execution stops after
> expression, and pair-dependent controls, calibration, typing, and assignment
> export are blocked. Pair-independent state calling remains explicit-only.
> The only active pairwise review is the
> mutable two-comparison upstream QuPath diagnostic documented in
> `scripts/review/README.md`.

This file was `MIGRATION_RESTORE_TO_CALIBRATION.md`, describing a migration that
was later reversed in place — so the filename had come to contradict its own
contents. It is renamed and records both directions.

## RESTORE was replaced, then restored (historical)

At that point, the RESTORE-to-calibration migration was judged a mistake:
RESTORE had been replaced by donor-marker calibration even though it was
intended to remain the production phenotyping method. RESTORE was consequently
restored as the producer behind the immutable `calibrate` stage key. That
historical restoration is not an active workflow now.

The stage name was not changed then because existing content-addressed manifests,
donor receipts, adoption plans, and the optional probabilistic interface refer
to `calibrate`. The method version and run identity now make its actual meaning
unambiguous:

```text
calibrate stage key
  -> restore-marker-evidence-v4-selected-inputs
     + restore-reference-seeded-v2-selected-inputs
  -> marker_evidence/ (tri-state calls + RESTORE normalized values)
  -> audit/restore_pair_models.parquet
```

The donor-calibration modules remained in the repository as historical and
test code, but that production pipeline did not call them.

## Improvements preserved by that historical restoration

Restoring RESTORE did not discard the later provenance and safety work:

- immutable cohort/run/stage manifests and donor receipts remained authoritative;
- geometry-QC-ineligible cells stayed in the object universe but could not estimate
  a RESTORE model or receive a marker call;
- the marker registry declared a measurement source separately from spillover
  policy: Iba1, CD68, and M2Gal3 use exact signed QuPath `Mean (norm)` fields,
  while other active markers retain their declared REDSEA source;
- signed QuPath norm values remained unchanged through ingest/expression and were
  left-censored only at RESTORE's nonnegative boundary, with the transform and
  counts audited;
- missing QuPath norm fields failed closed; raw, REDSEA, and
  RESTORE-normalized values cannot substitute;
- the then-active pair list had to match its versioned registry;
- CD11b is assay-QC-only with no RESTORE pair or typing role;
- M2Gal3 was ingested/selected from QuPath norm but deferred until a reference
  is validated;
- technical pair failures were `unavailable`, never silently negative;
- simultaneous broad typing was the label-free production default; the ordered
  residual hierarchy was NOT restored (see the section below), and the
  parent-constrained probabilistic path was an explicit signed-bundle option;
- old valid upstream stages could be adopted without recomputation.

## What to compare

Run `c8412757d12de36ca056` is evidence of a failed phenotyping implementation,
not a source of labels. At the time this was written its upstream
image/geometry/REDSEA artifacts could be adopted into a new run without
recomputation; the `adopt` command that did that was removed on 2026-09-16
along with the run framework it served, so nothing today can carry those
artifacts forward automatically. Its marker evidence, assignments, and QuPath
export must be recomputed regardless.

The then-current [REFINEMENT_WORKFLOW.md](REFINEMENT_WORKFLOW.md) described the
error-resolution loop. Its essential distinction was that QuPath review
identified a failing marker or rule while RESTORE refit donor-local thresholds
automatically. That file now documents the replacement development-only review
boundary.

## The ordered residual lineage was NOT restored

Reinstating RESTORE's marker calling and reinstating its ordered residual
hierarchy were two separate decisions, and only the first survived review.

Measured on run `3bf4fc8e391d007766d1`'s own marker evidence, over its
25,150,537 analysis-eligible cells, the declared gate order alone moved the
cohort Immune fraction from 20.50% (Immune first, as shipped) to 6.75%
(Epithelial first) to 3.33% (Immune last). 29.3% of cells are positive for two
or more broad gates and only 16.2% of Immune-gate-positive cells are positive
for no other gate. Under simultaneous typing 84% of the ordered residual's
Immune calls and 75% of its Endothelial calls are `Ambiguous`, against 4% of its
Epithelial calls.

The study question is the immune environment around islets across disease
progression, so the ordering was setting the headline number rather than the
evidence. `simultaneous_typing.py` and `broad_gates.json` replaced
`ordered_residual_typing.py` and `ordered_residual_rules.json`; the history is
at commit `ec97cbc`. See `WORKFLOW.md` § "Stage 7".

## Later corrections

- `restore-marker-evidence-v4-selected-inputs` separates a marker's
  `measurement_source` from its spillover policy, routes Iba1/CD68/M2Gal3 to
  exact QuPath `Mean (norm)` measurements, and renames evidence to the
  source-generic `<marker>__selected_intensity`. Finite signed negatives are
  preserved upstream and audited as `qupath_norm_left_censor_zero` only at the
  RESTORE boundary; M2Gal3 was deferred.
- `restore-marker-evidence-v2`: a reviewed absence of target expression no
  longer writes `negative` over cells that were never measured.
- `restore-pair-validation-v11`: `divisor_reproducibility` is NaN, not 1.0, when
  fewer than two locked seeds produced a divisor, so a pair can no longer pass
  the stability gate on one seed. Adds `seed_divisor_n`.
- SST is out of the broad Epithelial gate. `marker_registry.json` declares it
  `uses: ["subtype"]` with no `broad_compartment` and
  `hierarchical_typing.SUBTYPE_ONLY_MARKERS` forbids it as broad evidence; the
  ordered policy validated only the marker set, so nothing caught it.

## Historical divisor-statistic question (measured 2026-08-10; superseded)

This section records the evidence that rejected a single fixed `p999` rule.
The evaluator at that historical point used the registry-allocated finite-sample
`empirical_alpha` tail of the seeded control; the statistics below remain
historical comparisons, not production choices.

At the time, no `[restore] divisor_statistic` value was settled. It was frozen
to `p999` on 2026-08-10 (`4619bdf`) on cohort-median evidence. A per-pair
measurement on donors 6436, 6450 and 6566 — 72 donor-marker pairs, run through
the then-production `evaluate_locked_pair` — showed that no fixed statistic
worked for every marker and donor.

Agreement with each pair's own NNMF+floor anchor, the independent estimate of
the same positive population:

| statistic | within 0.5–2x of the anchor | median calls/anchor |
|---|---|---|
| manuscript `max` | 13/72 | 0.12 |
| `p999` (then-production) | 49/72 | 1.09 |
| `p99` | 18/72 | 2.58 |
| `p95` | 9/72 | 5.64 |
| mean+3SD | 21/72 | 2.56 |

`p999` is the best of the five and is still outside 2x on 23 of 72 pairs. Which
statistic is closest is **marker-dependent**: SMA wants `max` on all three
donors, EpCAM and INS want `p95`, most others want `p999`. A cohort median hides
this; the per-pair spread is the finding.

One replacement was tried and **rejected on measurement**, not opinion: purify
the control by fitting 1- vs 2-component Gaussian mixtures to log1p of the
control and keep the manuscript maximum over the negative mode. It fails because
the control's distribution is not two Gaussians — the fit splits distribution
shape, assigning ~95% of the control to the "contaminated" component for 18 of
24 markers, and returns `components_not_separable` for 31 of 72 pairs. Anchor
agreement drops to 10/41. It was not tuned into working.

Evidence: `data/pilots/divisor_survey/` (per-pair control distributions and every
candidate divisor) and `data/pilots/figures/fig1_divisor_landscape_<donor>.png`,
which draws each candidate on the marker's own control histogram.

The anchor is not ground truth — it is a second estimate from the same fit, so
agreement is evidence, not proof. Any replacement must be judged against each
donor's own data, never against an expected frequency or a disease ordering.

## Ambiguous is concentrated in one marker (measured 2026-08-10)

Under simultaneous typing, `Ambiguous` is 5.2–8.6% of cells across the three
pilot donors, and it is not spread evenly across gates:

- `Endothelial|Mesenchymal` 28–33% and `Immune|Mesenchymal` 20–42% of all
  Ambiguous cells;
- **Vimentin is positive in 55% of every Ambiguous cell**, and in 94% of the
  `Immune|Mesenchymal` collisions; SMA in 22%.

Vimentin is expressed by endothelium and leukocytes as well as mesenchyme, so as
a broad Mesenchymal gate member it fails the same mutual-exclusivity premise
that removed CD11b from broad Immune. This is recorded as a measurement, not yet
a decision: removing it is a reviewer's call, and the counterfactual has not
been run.

`Unavailable` needs no investigation: across all three donors it is 100%
`geometry_qc_ineligible`. No cell is Unavailable because of a failed pair.

## The intensity scale has no defined zero or maximum (measured 2026-08-11)

This heading and measurement apply to the raw/REDSEA source used in that
experiment. They do not redefine QuPath's signed `Mean (norm)` source: for
Iba1, CD68, and M2Gal3, QuPath explicitly anchors 0 at tissue p50 and 1 at
tissue p99. Those signed values are pre-RESTORE inputs, not positivity calls or
RESTORE-normalized outputs.

Every threshold in the workflow is a point on a scale whose ends were never
defined. Measured across the cohort:

- **No saturation.** 0 of 63 million sampled level-0 tissue pixels reach the
  16-bit ceiling; the brightest observed anywhere is 65324. The top of the scale
  is intact.
- **But the unit is donor-specific.** Exposure time differs between donors for
  16 of the 22 broad-gate markers (up to 1.67x, CD11c), and filter responsivity
  was recalibrated three times over the cohort (up to 1.41x, Cy5). The three
  calibration epochs track the acquisition batches, so this is a candidate
  mechanism for the batch confound.
- **Zero is a clamp, applied twice.** The scanner applies OffsetAndShading and
  AutofluorescenceSubtracted then stores unsigned, so 1-79% of tissue pixels are
  exactly 0 (11x variation between donors for DAPI). REDSEA clamps again at
  `redsea.py:529`, making 0.5-19% of analysis-eligible *cells* exactly 0, with
  ~3x donor variation for the same marker.

Scope, stated precisely: a per-donor-per-marker scalar **cancels** in RESTORE's
divisor arithmetic, so exposure normalisation will not move a single positivity
call, and REDSEA never mixes channels (`F` is a cell-cell contact matrix), so
mixed exposures do not corrupt spillover correction. What the censoring does
reach is `raw_otsu_floor`, which is computed over a distribution carrying a
donor-varying atom at zero — and the floor sets arm membership, hence the NNMF
separator and the target anchor that every divisor rule is judged against.

The specification is `docs/INTENSITY_SCALE.md`. Its governing rule is
standardize the unit, never the distribution: anchors may be properties of the
instrument or of the measurement's own censoring, never of the cell population.
