# RESTORE restoration note

The migration described by this filename was a mistake: RESTORE was replaced
by donor-marker calibration even though it was intended to remain the
production phenotyping method. RESTORE is now restored as the producer behind
the immutable `calibrate` stage key.

The stage name was not changed because existing content-addressed manifests,
donor receipts, adoption plans, and the optional probabilistic interface refer
to `calibrate`. The method version and run identity now make its actual meaning
unambiguous:

```text
calibrate stage key
  -> restore-marker-evidence-v1+restore-pair-validation-v10
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
