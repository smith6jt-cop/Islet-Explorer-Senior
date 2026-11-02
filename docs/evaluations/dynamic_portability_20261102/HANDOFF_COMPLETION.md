# Handoff completion and remaining scientific acceptance

The requested evaluation, project-only recorder, reproducible evidence bundle and
workflow documentation are delivered. This does **not** mean every scientific
route passed. The handoff explicitly limits scientific changes to evaluation and
remediation specifications; only KINTSUGI recording is added to production code.
The [requirement matrix](handoff_requirements.json) maps the handoff to its evidence
and distinguishes completed evaluation from unmet scientific acceptance.

The [workflow user guide](../../PORTABILITY_EVALUATION_USER_GUIDE.md) covers setup,
inventory and identities, every assay route, registration and matching, MALDI,
parameter recording, retries, diagnostics and reproducible verification. The
[report](README.md), [readiness matrix](readiness.json),
[acquisition ledger](acquisition_dispositions.json) and
[remediation specifications](remediation.md) provide the current conclusions.

## Meaning of blocked

`Blocked` means an identified requirement prevents a route from producing the
requested supported output. It is an evidence-backed disposition, not a running
job or an omitted acquisition. All 201 declared acquisitions remain accounted
for, including 20 post-Xenium IF acquisitions excluded from this task's scope.

- Protein refinement was executed across the configured cohort and stopped at
  0112 after all 25 primary tables were built. A present but unacquired Ker8_18
  column eliminates the measurable population. Fixing that scientific admission
  rule is specified as PROT-1; it is outside this implementation scope.
- PC33 lacks a supported lineage registry. The dedicated-pLN candidate route has
  bound upstream results for only two of 19 acquisitions. Spleen RNA has audited
  inputs but no independently specified complete lineage route.
- Same-slide pairs lack accepted registration/link receipts. Serial MALDI lacks
  accepted neighbour/frame/registration evidence and complete pLN/Xenium routes.
  Morphological similarity and matching donor IDs cannot supply those facts.
- CODEX's original fixture and exact CSV checks expose reproducibility issues;
  the recovered journal, original-source comparisons and precision diagnostics
  are separate evidence. v20 has 22 mismatched current-file bindings.
- Biological transfer criteria and some resolved runtime values are unavailable.
  Static source tracing explicitly identifies those limits rather than filling
  them with assumptions or converting candidate results into approvals.

## Final additions to the baseline evaluation

The CODEX follow-up recovers the original execution journal, checks the input and
output chain, and executes all 20 images' numerical assertions without editing
the original scientific implementation or weakening tolerances. Historical guard
failures remain visible. The parameter trace covers every catalog declaration,
registry row and remaining literal candidate, with the actual historical RNA,
pLN and spleen producer versions kept separate from current code.

All 32 MALDI sections were visually compared against the inventory source figure
in addition to the 32 value/spatial diagnostics. The supplier archive and workbook
were inspected; they do not supply MALDI cutting order or section distance. The
pLN object audit resolves 19 acquisitions to 47 annotation UUIDs, including three
non-area LineStrings. The source guide and dated MALDI plan now document the
10.81% maximum mask exclusion and its denominator effect.

The unchanged scientific modules remain the baseline. Inactive, unvalidated
scientific draft patches are retained only under ignored evaluation outputs;
they are neither installed nor part of the committed implementation. Their useful
requirements are incorporated into MALDI-2 and MALDI-3. No registry, transform,
cell identity, production release or biological result was approved.

## Verification and storage

The recorder has 19 passing tests; the existing route suite has 49 and the pinned
port unit suite has 37. Scientific parity and cohort reconciliation have their
own results in [validation.json](validation.json), so numerical failures are not
hidden in a unit-test total. Actual MCP write/retry/history/export checks and the
delivery verifier check persistent records, immutable revisions, evidence hashes,
global-history isolation and zero approval. Detailed logs and input receipts are
hashed by [artifact_manifest.json](artifact_manifest.json).

SQLite databases, journals, locks, retry outboxes, environments and run output
trees are ignored. Reviewed JSON exports, inventories, matrices and reports remain
tracked. Large parameter traces are stored as deterministic `.json.gz` artifacts;
inspect them with `gzip -dc` or Python's `gzip.open`. Raw data stay in their
original stores, referenced by identifiers, paths and hashes.
