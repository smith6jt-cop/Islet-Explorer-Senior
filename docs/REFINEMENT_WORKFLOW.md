# Refinement workflow

## Purpose

This is the recourse when a QuPath review shows unacceptable calls. It must not
force the project into hand-labelling thousands of cells or training a
QuPath-style object classifier.

The review is an **acceptance test for a versioned method change**, not the
computational engine and not a training set. RESTORE refits every donor
automatically after each versioned change, so the loop needs tens of error
exemplars, not comprehensive labels.

## What the pipeline does

The immutable stage key is still `calibrate` so completed upstream runs remain
adoptable. Its producer is
`restore-marker-evidence-v2+restore-pair-validation-v11`:

1. For every donor and every accepted target/reference pair in
   `marker_registry.json`, find the two mutually exclusive groups by NNMF.
2. Use the reference-positive/target-negative group to estimate a donor-local
   divisor. The production policy is `reference_only` with the `p999`
   statistic.
3. Apply the divisor to every analysis-eligible cell and emit `positive`,
   `negative` or `unavailable`, plus the normalized intensity and the full
   pair-fit audit.
4. Pass those states to **simultaneous broad typing** (`broad_gates.json`).
   Every gate is evaluated against the same evidence, with no precedence:
   one supported gate is an `anchor`, more than one is `Ambiguous`, none with
   everything explicitly negative is `Other`, and an unmeasured gate marker
   leaves the cell `Unavailable`.
5. Treat any technical donor-marker failure as `unavailable` (or halt, if that
   policy is selected), never as a biological negative.

`--pipelined` changes only scheduling. It overlaps donor-ready stages under the
same dependency graph and the same scientific functions.

## The correction loop

Do not prepare a labelling package and do not run `probabilistic fit`.

1. Keep a small error ledger of representative visible failures. It needs only
   `donor_id`, `object_id`, `error` (`false_positive` or `false_negative`), the
   implicated marker or reported class, and a note. These are regression
   exemplars, not training labels.
2. Join those stable IDs to `marker_evidence/` and
   `audit/restore_pair_models.parquet` in the candidate run, and find **where
   the first wrong decision happens**:

   | first wrong decision | correction surface |
   |---|---|
   | raw image, segmentation or spillover is wrong | fix the manifest, geometry policy or REDSEA and rerun downstream |
   | the cell should not have been analysed at all | geometry QC — check `analysis_reasons`, not the marker calls |
   | a RESTORE marker state is wrong | change the target/reference rule or a declared divisor policy, version it, rerun `calibrate` |
   | marker states are right but the class is wrong | change `broad_gates.json`, version it, rerun `type` |
   | the cell is `Ambiguous` | that is not automatically an error — check whether the competing gates are both genuinely positive (`broad_supported_gates`, `broad_support_margin`) before changing anything |
   | evidence cannot decide because a gate is unavailable | preserve `Unavailable`; do not force a class |

3. Rerun the whole affected dependency suffix. Never overwrite an individual
   cell assignment. A new configuration, registry, source or audit input creates
   a new run ID automatically.
4. Reopen the new QuPath export and check both the recorded failures and a small
   untouched set of cells. Also check the cohort guardrails: per-donor class
   fractions, impossible co-expression, islet/endocrine localization, vascular
   CD31, acinar keratin retention. A correction is accepted only if it fixes the
   reported pattern without creating a new one.

### Judge a correction against the donor's own data

Never against another donor, a disease label, or an expected frequency. The
per-donor histogram and the pair's own negative-control population are the
reference. `audit/restore_pair_models.parquet` carries everything needed:
`target_fold`, `reference_fold`, `double_high_fraction`,
`divisor_reproducibility`, `seed_divisor_n`, `control_fraction_le_1`,
`canonical_divisor` beside `manuscript_divisor`, and the call counts.

## Known open questions in the current method

These are recorded so a reviewer does not have to rediscover them. Each is
measured on run `3bf4fc8e391d007766d1` (20 donors, 25,150,537 analysis-eligible
cells).

**The divisor statistic was p99 and is now p999.** *(Resolved.)* A quantile
divisor calls the upper (1−q) tail of the pair's own negative control positive
by construction. Under `p99` that was 1% of every control, and 417 of 480
donor-pairs called more positives than their own NNMF+floor target anchor
(median 2.15×, max 51×), worst on the sparse immune markers the study depends
on. Measured on 20 donors × 24 pairs under the real `evaluate_locked_pair`:

| | p99 | **p999** | max (manuscript) |
|---|---|---|---|
| callable pairs | 480/480 | **480/480** | 448/480 |
| calls exceeding own anchor | 417/480 | 240/480 | 6/448 |
| median calls ÷ anchor | 2.15 | **1.00** | 0.09 |
| control normalising ≤ 1 | 0.990 | 0.999 | 1.000 |

`p999` keeps every pair callable, as `p99` does, and centres the call on the
fit. The manuscript maximum makes every control normalise ≤ 1 by construction
but loses 32 pairs to `MODEL_UNSTABLE` and under-calls against its own anchor.
Note that 240 of 480 pairs still exceed their anchor under `p999`: the
distribution is centred, not tight.

**The seed-stability gate is non-binding on half the cohort.** 253 of 480 pairs
carry `stability_informative = False` — an arm was saturated, so seed agreement
is arithmetically guaranteed and is not evidence of robustness. All 253 are
`VALID` with median reproducibility exactly 1.000000. The flag is reported and
gates nothing, deliberately: making it gate would fail more than half the
cohort. **Surface it in every per-donor review** so a reader can tell a
reproducible pair from one whose stability was never tested.

**CD11b was dropped from broad identity.** *(Resolved.)* `double_high_fraction`
median 0.173 and max 0.407 against ≤0.010 for the other 23 accepted pairs: its
negative control was 17% target-positive. Screening every candidate reference
showed the failure belongs to the marker — the only clean control, CD3e
(0.0015), is undersized on all 20 donors. CD11b remains an active, fitted
RESTORE target and subtype evidence; it no longer gates broad Immune. Cost:
37,680 Immune anchors, 0.150% of eligible cells, because 88.3% of the cells
whose only immune evidence was CD11b were already `Ambiguous`.
`max_double_positive_fraction` is now a reported per-pair flag
(`double_high_exceeds_threshold`) rather than dead code — it still gates
nothing, because a violation is a reviewer's decision, but it can no longer go
unnoticed.

**14 of 24 pairs fit on whole cell alone, and that is accepted.** The extracted
compartment sample covers 12 markers, so most pairs get `Cell` only while the
rest get `Membrane+Cell`. Six declared feature sets were compared on 20 donors
under the production divisor: 99.5–100% of pairs VALID in every set,
reproducibility 0.9994–0.9999, `anchor_recovery` 0.980–0.996, `target_fold`
14.8–15.1, `double_high_fraction` 0.0023–0.0025, and positives within 6% of
production. That matches the repo's own earlier measurement — 8 features versus
4 moved 12 of 14 donor-pairs by ≤0.2% (`a64b949`). The NNMF outcome is
insensitive to the feature set here, so re-extracting compartments for all 24
markers was **not** done: it is a uniformity improvement, not an accuracy one.
The source CSVs do carry `Membrane:` and `Cytoplasm:` means for every gate
marker should that change.

> When screening any of this with an ad-hoc script, note that
> `PairValidationConfig()` defaults `divisor_statistic` to `max`, while
> production passes `p999` through `RestoreStageConfig.pair_config()`. Validate
> the harness against a real run's `audit/restore_pair_models.parquet` before
> believing it.

## Reusing a previous run

`pipeline adopt plan --from-run <id>` reports, per stage, whether the previous
output can be reused. It is strict on purpose, and two things currently block it:

- **A run produced from a dirty tree can never be adopted.** Adoption requires
  the source code to be reconstructible from the commit the manifest records.
  Runs `a43be26e1d576d3fe495` and `3bf4fc8e391d007766d1` recorded `dirty: true`;
  their code now exists at commit `ec97cbc`, but their manifests point at
  `ab9a578`, so the reconstruction check fails. Only `ingest` is reusable.
- **A producer-file change needs a recorded review.** `adoption_compatibility.json`
  holds transitions that name the exact source and target file hashes, the exact
  changed paths and the sha256 of the exact diff, and adoption re-verifies all of
  them. A transition is a human sign-off that the diff cannot change that
  stage's output; do not add one without making that judgement.

Consequence: after a code change to `donor_pipeline.py` or any producer, either
record a reviewed transition or recompute. Recomputing geometry and REDSEA cost
about 4.5 hours for this cohort (REDSEA alone was 238 minutes in
`c8412757d12de36ca056`), and yields a clean run every later calibrate iteration
can adopt from in roughly half an hour.

## Which artifact is authoritative

- `marker_evidence/donor_id=*/…parquet` — per-cell RESTORE states, normalized
  values, divisors, references and method version.
- `audit/restore_pair_models.parquet` — one row per donor × marker fit: pair
  state, divisor, feature set, stability, control counts and call counts.
- `assignments/donor_id=*/…parquet` — simultaneous broad typing applied to that
  evidence.
- `qupath_export/cell_assignments_*.csv` — visual review surface only.

## The optional probabilistic lane

A probabilistic rescue for genuinely unresolved cells remains a possible future
method. It is **not** the remedy for a typing rule that is wrong, and marker
evidence from a rejected run must never be supplied to a fit.

If that lane is ever authorized, its guarded release interface is
[06_probabilistic_validation_and_promotion](../notebooks/06_probabilistic_validation_and_promotion.ipynb),
and its contract is in [`PROBABILISTIC_RESCUE.md`](PROBABILISTIC_RESCUE.md).

Notebooks 01–06 are operator interfaces, not alternate production
implementations. Notebooks 01–04 retain their existing guarded package-stage
entry points; notebook 05 is read-only. Promotion stays disabled until the
probability-sample review is complete and its software-enforced provenance is
valid. The detached challenge sample may reveal a catastrophic problem, but it
cannot supply headline metrics or software promotion pass/fail.

Note that notebook 05 currently requires a configured model bundle and
therefore does not run against the production default.
