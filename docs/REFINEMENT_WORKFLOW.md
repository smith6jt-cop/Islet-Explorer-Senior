# RESTORE phenotyping refinement workflow

## Purpose

This workflow supplies the recourse that was missing after run
`c8412757d12de36ca056`: a QuPath review may demonstrate unacceptable false
positives or false negatives, but it must not force the project into manually
labeling thousands of cells or training a QuPath-style object classifier.

RESTORE remains the production marker-calling method. The review is an
acceptance test for a versioned method change, not the computational engine and
not a training set.

## What the pipeline does

The immutable stage key is still named `calibrate` so completed upstream runs
can be adopted safely. Its producer is now
`restore-marker-evidence-v1+restore-pair-validation-v10`:

1. For every donor and every active target/reference pair in
   `marker_registry.json`, use raw QuPath compartment means to find the two
   mutually exclusive groups by NNMF.
2. Use the reference-positive/target-negative group to estimate a donor-local
   target divisor. The production policy is `reference_only` with the robust
   `p99` divisor.
3. Apply the divisor to every analysis-eligible cell and emit
   `positive`, `negative`, or `unavailable`, plus the normalized intensity and
   complete pair-fit audit.
4. Pass those RESTORE states into the versioned ordered residual hierarchy in
   `ordered_residual_rules.json`: Immune, Endothelial, Epithelial, Neural,
   Mesenchymal, then Other. No reviewed cell labels or probabilistic bundle is
   required.
5. Treat any technical donor-marker failure as `unavailable` (or halt if that
   policy is explicitly selected), never as a biological negative.

`--pipelined` changes only scheduling: it overlaps donor-ready stages while
respecting this same dependency graph and these same scientific functions.

## Simplest correction after a failed QuPath review

Do not prepare a new labeling package and do not run `probabilistic fit`.

1. Keep a small error ledger of representative visible failures. It needs only
   `donor_id`, `object_id`, `error` (`false_positive` or `false_negative`), the
   implicated marker or reported class, and a note. These are regression
   exemplars, not training labels.
2. Join those stable IDs to `marker_evidence/` and
   `audit/restore_pair_models.parquet` in the candidate run. Determine where
   the first wrong decision occurs:

   | first wrong decision | correction surface |
   |---|---|
   | raw image, segmentation, or spillover is wrong | fix the manifest, geometry policy, or REDSEA and rerun downstream |
   | RESTORE marker state is wrong | change the target/reference rule or a declared RESTORE divisor policy, version it, and rerun RESTORE |
   | marker states are right but the ordered class is wrong | change `ordered_residual_rules.json`, version it, and rerun typing |
   | evidence cannot decide because a required gate is unavailable | preserve `Unavailable`; do not force a class |

3. Rerun the whole affected dependency suffix. Never overwrite an individual
   cell assignment. The new configuration, registry, source, and audit inputs
   create a new run ID automatically.
4. Reopen the new QuPath export and check both the recorded failures and a
   small untouched set of cells. Also check cohort guardrails: per-donor class
   fractions, impossible co-expression, islet/endocrine localization, vascular
   CD31, acinar keratin retention, and disease trends. A correction is accepted
   only if it fixes the reported pattern without creating a new one.

This loop usually requires review of tens of error exemplars, not comprehensive
cell labeling. RESTORE refits thresholds automatically for every donor after
each versioned change.

## Reuse the valid part of `c8412757d12de36ca056`

That run's marker evidence, assignments, and QuPath exports are rejected. Its
ingest, geometry, REDSEA, selected expression, and independent state artifacts
remain valid and are eligible for provenance-checked adoption. Controls are
recomputed because the accepted CD11b reference changed; RESTORE (`calibrate`)
and `type` are then recomputed before regenerating the QuPath export. The
generated adoption plan records the exact reusable stages and byte count.

```bash
python -m phenocycler.pipeline adopt plan \
  --config config.ini \
  --from-run c8412757d12de36ca056 \
  --out restore-adoption-plan.json

# Read the printed stage decisions. It should adopt ingest, geometry, REDSEA,
# expression, and states; controls, calibrate, and type should recompute.
python -m phenocycler.pipeline adopt apply \
  --config config.ini \
  --plan restore-adoption-plan.json \
  --continue
```

The configured exhaustive compartment source is
`../data/restore_pair_validation/qupath_compartment_full_20donor`. It is
fingerprinted as an external RESTORE input. It covers the current analysis-QC
cell universe and prevents an unnecessary rescan of the 148 GB QuPath
measurement CSV.

The accepted donor 6566 INS image review is separately pinned in
`phenocycler/restore_expression_reviews.csv`. RESTORE records that donor-marker
as `NO_TARGET_EXPRESSION_CONFIRMED`; this explicit biological decision must not
be confused with a technical failure or inferred automatically for another
donor.

## Which artifact is authoritative

- `marker_evidence/donor_id=*/...parquet`: per-cell RESTORE states,
  normalized values, divisors, references, and method version.
- `audit/restore_pair_models.parquet`: one row per donor × marker fit, including
  pair state, divisor, feature set, stability, control counts, and call counts.
- `assignments/donor_id=*/...parquet`: ordered residual hierarchy applied to
  RESTORE evidence.
- `qupath_export/cell_assignments_*.csv`: visual review surface only.

The old `c8412757d12de36ca056/marker_evidence` must not be supplied to a
probabilistic fit. A probabilistic rescue remains an optional future method for
genuinely unresolved cells, but it is not the remedy for having accidentally
removed RESTORE.

If that optional lane is ever authorized after RESTORE has passed review, its
separate guarded release interface remains
[06_probabilistic_validation_and_promotion](../notebooks/06_probabilistic_validation_and_promotion.ipynb).

That optional lane retains its prior safety contract. Notebooks 01–06 are operator interfaces,
not alternate production implementations. Notebooks 01–04 retain their existing guarded package-stage
entry points; notebook 05 is read-only. Promotion remains disabled until the
probability-sample review is complete and its software-enforced provenance is
valid. The detached challenge sample may reveal a catastrophic problem, but
it cannot supply headline
   metrics or software promotion pass/fail.
