# Additional CODEX reference-context correction received

The immutable `CODEX-handoff-20261102-lymphatic-isolation` package was received
and reviewed during this evaluation. Its checksum-list hash is
`67a56038efea9e03fa10f7ead0e5f620bf5d18cdcd6a425232fbe1077dbb47dc`;
all 17 payload checks pass. The [receipt](codex_received_handoff.json) binds each
file, package source commit `29f45b35dc1566d4fce89d0e8813fa5b55743c4c`, and
isolation fix `a425c05167a00300c0a644a19f52751baa6ba901`.

This is newly received source evidence, separate from the pinned CODEX parity
comparison and the running PhenoCycler cohort. No package code, thresholds,
labels or biological selections were applied to those runs. The supplied audit
reports a passing all-20-image CODEX structure-only replay; this evaluation
verified the delivered artifact, rather than rerunning that separate replay.

The correction prevents a Lymphatic geometry revision from silently changing
shared calibration populations and unrelated identities. It separates a
complete, explicitly selected baseline reference context from the current
structure-membership context. Its audit requires unchanged measurement, QC,
calibration, marker and non-target identity evidence. A successful change-scope
audit is distinct from biological accuracy or acceptance of a mixed v19/v20
candidate. The package states that no such mixed candidate is accepted.

## Adaptation requirements

The shared package currently exposes a full typing generation; it has no
equivalent `structure_update` / `calibration_context` command. Its existing
`scripts/typing/context.py`, `subpanel.py` and `run_calibration_arm.py` must not be
presented as an isolated structure revision. A future explicit route needs:

1. A selected baseline bound by hashes to the acquisition roster, native cell
   IDs, coordinates, QC, segmentation, measurements, panel/registry, fitted
   calibration and complete reference context. Read-only baseline inputs and
   separate writable destinations must prevent historical-output replacement.
2. Separate live structure membership and fixed calibration-reference readers.
   A structure-only invocation must reject recalibration, rebinding, missing
   reference fields and reuse of an already populated identity destination.
3. Panel/profile role mapping through the receiving registry. CODEX marker
   names, numeric thresholds, exclusions, endocrine membership and labels are
   not transferable defaults for PhenoCycler, Xenium, pLN or spleen.
4. An automatically invoked whole-cohort invariant audit. Verify unchanged
   UUIDs, measurements, QC, all marker evidence, calibration, non-target context
   and identities; constrain permitted changes to the declared target geometry
   and identity scope. Unknown fields remain immutable. An unchanged-geometry
   control must reproduce its baseline. Failed checks must prevent acceptance.
5. A separately reviewed coherent biological baseline and downstream rebuild.
   Preserve assay-specific measurement contracts and unresolved calls. Do not
   copy labels to conceal inconsistent evidence or infer acceptance from a
   diagnostic replay.

This specification extends the existing external-port dependency. It does not
change the 201-acquisition roster or any completion, pairing, registration or
approval disposition in this remediation delivery.
