# Development review workflow

## Current status

This file previously described refinement of a production RESTORE and typing
path. That path is superseded. There are no verified production
target/reference pairs, and ordinary execution stops after expression.
Pair-dependent controls, calibration, typing, and assignment export are
blocked. Pair-independent state calling is outside this review and remains an
explicit-only stage.

Historical probabilistic, marker-debug, and Otsu refinement commands are not
current alternatives and must not be used to create labels.

## Purpose

Current review is limited to determining whether a small set of upstream
QuPath measurements is technically interpretable despite segmentation and
mixed-scale effects. It does not tune a positivity threshold or assign a cell
type.

The mutable shipped roster contains only:

1. raw CD3e versus raw Ker8_18; and
2. exact signed QuPath Iba1 `Mean (norm)` versus raw Ker8_18.

CD79a is excluded as a pair candidate because its expression is too sparse in
normal tissue. CD68, MPO, EpCAM, E-cadherin, CD31, and Pan-cytokeratin are also
absent from this reduced review roster.

## Review loop

1. Refresh the narrow projection directly from the declared full QuPath
   exports.
2. Confirm exact source fields, UUID uniqueness, finite coverage, and signed
   Iba1 behavior.
3. Review donor-local rank association and exact co-high-tail overlap.
4. Use boundary/contact crops only to understand segmentation contamination;
   do not turn reviewed cells into training labels.
5. Edit and rerun the mutable roster if a comparison is uninformative. A roster
   edit has no production effect.

Commands:

```bash
python scripts/review/qupath_fundamental_inputs.py \
  --donors 6436 6450 6566

python scripts/review/fundamental_pair_survey.py \
  --donors 6436 6450 6566

python scripts/review/reference_overlap_audit.py \
  --donors 6436 6450 6566
```

The mixed Iba1/Ker8_18 comparison is valid only as a rank comparison. QuPath
Iba1 normalization is signed and marker-specific; it is not numerically
commensurate with raw Ker8_18 and `1` is not a positivity threshold.

## Output boundary

Development outputs live under `../data/development` and are replaceable only
when the destination carries the tool's ownership marker. They contain no run
ID, stage manifest, call, assignment, or production approval.

The content-addressed run/manifests architecture remains in force for ingest,
geometry, REDSEA, and expression. Development review must not write inside a
run directory or alter those upstream artifacts.

## What would be required to resume production work

A future pair-calling workflow requires a new, independently reviewed marker
registry and implementation contract, explicit mixed-source tests, donor-level
acceptance criteria, and a fresh production authorization. Historical RESTORE
pair lists, divisors, model bundles, classifications, or outputs cannot satisfy
that requirement.
