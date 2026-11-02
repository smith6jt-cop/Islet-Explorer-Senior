# GEOM-2: same-slide candidate-link remediation

The shared candidate-link implementation now refuses arbitrary centroid ties,
retains independent assay calls, reports every cell's disposition and evaluates
sensitivity to an explicitly supplied registration uncertainty bound. The
Stage 07/08 matching cap remains **5 µm**. These changes do not establish an
accepted registration or same-cell truth for any cohort pair.

## Implementation

- [paired_typing.py](../../../core/multiplex5000/multiplex5000/integration/paired_typing.py)
  emits unique mutual-nearest links, all nearby candidate edges and one disposition
  per ingested cell. Exact equal-distance nearest alternatives remain ambiguous in
  either direction; sorting by native cell keys makes these outputs invariant to
  independent input row permutations. Counts distinguish exclusions, ambiguity,
  nonmutual alternatives, no candidate within the cap and candidates reachable
  only under the supplied uncertainty bound. Missing counts for excluded cells
  remain unavailable. Multiple candidates are not asserted segmentation splits
  or merges without overlap evidence.
- Every supplied protein and RNA call column is retained under its assay prefix.
  Legacy `gene:*` columns keep their existing names and unavailable values. Actual
  RNA identity columns can replace or accompany that gene evidence. Original
  labels, call states and confidence columns remain separate. A conflicting
  reserved output name is rejected instead of overwriting independent evidence.
- Optional `same_slide_hierarchy_map` project settings declare the exact tissue,
  comparison level, independent call columns and their label mappings. Optional
  state filters declare which states permit comparison. Outputs say `agree`,
  `disagree`, `unmapped` or `unavailable`; no combined identity is produced.
- [sameslide.py](../../../core/multiplex5000/multiplex5000/integration/sameslide.py)
  applies the same no-tie rule to the legacy centroid matcher. Unambiguous
  candidates and its existing three-column return interface remain compatible.
- [Stage 07](../../../core/multiplex5000/multiplex5000/stages/same_slide_registration.py)
  verifies that the exact current cells match the accepted Stage 06 input hashes.
  [Stage 08](../../../core/multiplex5000/multiplex5000/stages/same_slide_integration.py)
  verifies the accepted Stage 07 receipt, exact transform and displacement,
  source manifest, cells, region table, and Stage 06 receipt and acceptance. It
  independently verifies the exact Stage 06 call output hashes. Missing, stale
  or duplicate bindings fail closed. Optional cell/call frame, unit and
  segmentation fields remain visible to validation, including fields present
  only in the call table; mixed source frames are rejected.

Stage 08 writes `candidate_links.parquet`, `spatial_candidates.parquet`,
`cell_link_dispositions.parquet`, `linking_metadata.json`, `pairing_counts.json`
and the receipt-bound `resolved_parameters.json`. The last artifact records the
actual fixed cap, selected matching/comparison methods, supplied uncertainty,
producer source hashes and input receipt hashes through the existing recording
helper. Recording does not grant approval or enable automatic application.

## Explicit uncertainty contract

The stage reads an optional `registration_uncertainty` record from the accepted
transform metadata. The shared function also accepts the same record explicitly:

```json
{
  "bound_um": 0.2,
  "frame": "crop_local_um",
  "kind": "radial_displacement_bound",
  "evidence_ref": "identifier-of-the-actual-bound-evidence"
}
```

This example is synthetic, not an estimated bound for any acquisition. A bound
must be finite, nonnegative, expressed in the exact fixed frame and accompanied
by an evidence reference. With no supplied bound, uncertainty, distance intervals
and stability are unavailable; no zero error or candidate-distance estimate is
substituted. The metadata distinguishes declared-only coordinate frames from
explicit cell metadata and preserves missing segmentation provenance as null.

For a supplied bound `u`, all edges within `5 + u` µm are visible. Nominal links
still require unique mutual nearest neighbours within **5 µm**. Each candidate
distance `d` has the conservative interval `[max(0, d-u), d+u]`. A nominal link
is marked stable under the provided bound only if its upper distance remains
within the cap and both nearest-versus-second-nearest margins exceed `2u`.
This is a sufficient conservative test when each transformed RNA centroid may
move within `u`; it does not estimate registration error or model segmentation
error. Alternatives and boundary sensitivity remain visible without replacing
the nominal result or tuning for assay agreement.

## Verification and remaining scope

**78 tests passed**, covering both focused same-slide files and the existing
crossmodal integration regression file; Ruff and `git diff --check` passed for
the changed files. Tests exercise exact ties in both directions, symmetric ties,
row permutations, duplicate identifiers, all-ineligible RNA, stale receipts,
stale call tables, changed cells, mixed frames, missing displacement, source
metadata, absent uncertainty, zero and positive supplied bounds, cap-boundary
sensitivity, nearest-neighbour competition, actual RNA calls, unavailable calls
and explicit-level disagreement. Nullable missing source metadata is rejected
when it differs from ingest; pandas' default skip-null aggregation cannot admit
an unavailable typing-eligibility value. The unambiguous legacy pilot retains its two
original candidate links and native cell identifiers.

Evidence is under
[`outputs/blocker_remediation_20261102/geom2`](../../../outputs/blocker_remediation_20261102/geom2):
seven labelled examples in `synthetic_diagnostics.json` and the original
77-test implementation record. The final nullable-metadata repair and 78-test
result are in
[`geom2_revision1`](../../../outputs/blocker_remediation_20261102/geom2_revision1),
with `final_tests.log/xml`, `lint.log`, `completion_summary.json` and
`artifact_manifest.json`. The revision links the prior immutable package by hash.
Earlier attempts are retained separately; the failed interpreter invocation did
not run tests. Synthetic review acceptance occurs only in temporary test projects.

Full-cohort routing and scientific acceptance remain separate. The audited
same-slide ledger has four pancreas declarations and eight verified spleen
pairs; pLN physical-section relationships remain unresolved. The existing
bounded Stage 07/08 project selects one spleen region. No cohort transform,
uncertainty bound, hierarchy map, segmentation overlap evidence or accepted
registration has been fabricated, and no historical output or scientific gate
was changed. Real pair evaluation still requires its actual independent calls,
exact section evidence, accepted geometry and review criteria. This remediation
closes the candidate-link implementation defect; it does not turn missing cohort
inputs or registration review into completed validation.
