# Method history

This file was `MIGRATION_RESTORE_TO_CALIBRATION.md`, describing a migration that
was later reversed in place — so the filename had come to contradict its own
contents. It is renamed and records both directions.

## RESTORE was replaced, then restored

The RESTORE-to-calibration migration was a mistake: RESTORE was replaced
by donor-marker calibration even though it was intended to remain the
production phenotyping method. RESTORE is now restored as the producer behind
the immutable `calibrate` stage key.

The stage name was not changed because existing content-addressed manifests,
donor receipts, adoption plans, and the optional probabilistic interface refer
to `calibrate`. The method version and run identity now make its actual meaning
unambiguous:

```text
calibrate stage key
  -> restore-marker-evidence-v2+restore-pair-validation-v11
  -> marker_evidence/ (tri-state calls + RESTORE normalized values)
  -> audit/restore_pair_models.parquet
```

The donor-calibration modules remain in the repository as historical and test
code, but the production pipeline does not call them.

## Preserved improvements

Restoring RESTORE does not discard later provenance and safety work:

- immutable cohort/run/stage manifests and donor receipts remain authoritative;
- geometry-QC-ineligible cells stay in the object universe but cannot estimate
  a RESTORE model or receive a marker call;
- REDSEA compartment outputs remain the corrected source for divisors/calls;
- the active pair list must exactly match the accepted RESTORE v10 pair web;
- the reviewed `CD11b <- EpCAM` pair is restored, replacing the later
  unvalidated `CD11b <- E_cadherin` change;
- technical pair failures are `unavailable`, never silently negative;
- the historical ordered residual broad hierarchy is again the label-free
  production default; the simultaneous probabilistic path remains only as an
  explicit signed-bundle option;
- old valid upstream stages can be adopted without recomputation.

## What to compare

Run `c8412757d12de36ca056` is evidence of a failed phenotyping implementation,
not a source of labels. Its upstream image/geometry/REDSEA artifacts can be
adopted, but its marker evidence, assignments, and QuPath export must be
recomputed.

Use [REFINEMENT_WORKFLOW.md](REFINEMENT_WORKFLOW.md) for the error-resolution
loop. The essential distinction is that QuPath review identifies a failing
marker or rule; RESTORE refits donor-local thresholds automatically. Review
examples are regression tests, not an object-classifier training package.

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
evidence. `simultaneous_typing.py` and `broad_gates.json` replace
`ordered_residual_typing.py` and `ordered_residual_rules.json`; the history is
at commit `ec97cbc`. See `WORKFLOW.md` § "Stage 7".

## Later corrections

- `restore-marker-evidence-v2`: a reviewed absence of target expression no
  longer writes `negative` over cells that were never measured.
- `restore-pair-validation-v11`: `divisor_reproducibility` is NaN, not 1.0, when
  fewer than two locked seeds produced a divisor, so a pair can no longer pass
  the stability gate on one seed. Adds `seed_divisor_n`.
- SST is out of the broad Epithelial gate. `marker_registry.json` declares it
  `uses: ["subtype"]` with no `broad_compartment` and
  `hierarchical_typing.SUBTYPE_ONLY_MARKERS` forbids it as broad evidence; the
  ordered policy validated only the marker set, so nothing caught it.
