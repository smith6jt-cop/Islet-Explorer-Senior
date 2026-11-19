# GEOM-1: batch-routing assessment

Assessment only; no new route or scientific gate is implemented by this document.

The eight spleen pairs have a usable explicit source-manifest schema. Batch
selection is already possible at RNA ingest, and protein measurement/normalization
admit all eight sections. Three later assumptions enforce a single paired route:

1. `stages/spleen_typing_pilot.py::_selected` requires exactly one pair, even
   though its input declarations and execution loop are otherwise per pair.
   Its RNA output is identity-gene evidence, not independent tissue lineage calls.
2. `stages/same_slide_registration.py::_selected` requires one RNA acquisition;
   registration inputs, copied transform, summary and figures use singleton names.
3. `stages/same_slide_integration.py` inherits that selector, uses singleton
   inputs/outputs and authenticates a single accepted Stage 07 receipt. The
   project controller requires acceptance of the preceding stage as a whole.

Removing the selectors alone would therefore be insufficient. One missing
transform would fail the entire run; treating successful rows as the whole stage
would silently shrink the cohort. A stage-wide acceptance must not accidentally
stand in for review of every pair.

A bounded batch coordinator is tractable without new RNA lineage biology:

- Accept an explicit list of source-manifest pair IDs and per-pair project/release
  roots, with no donor or filename inference. Account for every requested pair.
- Run each pair through the existing validated Stage 07/08 contracts in its own
  output root, with its own input hashes, receipt and review acceptance. Reuse
  artifacts only by their content-bound source and measurement contracts; never
  fabricate Stage 06 or Stage 07 acceptance when preparing a project.
- Write a separate batch ledger with completed, blocked, unavailable and excluded
  dispositions, reason codes and receipt/evidence references. Missing geometry,
  cells, calls or acceptance blocks that pair only. A partial ledger is not
  full-cohort success or an accepted scientific stage.
- Keep spatial candidate linking and identity comparison separate. Existing
  gene-only RNA tables may support an explicitly labelled spatial/gene-evidence
  route; they cannot satisfy a request for independent RNA hierarchy calls.
  Full mapped comparison remains unavailable until real calls and a declared
  hierarchy mapping are present.
- Test two ready pairs plus one missing-transform pair, stale per-pair receipts,
  duplicate pair configurations, pair-frame mismatch, distinct output paths,
  complete roster reconciliation and unchanged single-pair behavior.

An in-project batch implementation is a larger alternative: namespace all pair
inputs/outputs, define independently inspectable per-pair gate acceptance, and
decide explicitly how those gates relate to the global Stage 07/08 gate. That
gate-model change should not be hidden inside a loop refactor.

IO60 pancreas/pLN joint typing is not a drop-in source for this spleen contract.
`io60_mixed_protein_pilot` writes `mixed/donor_id=<donor>/data.parquet`, consumes
`aligned:<donor>` inputs and embeds tissue/object context, whereas the same-slide
route expects a recorded physical section's `protein/<entry>.parquet`, exact
section cell keys, and `protein_cells:<entry>` input hashes. A future adapter
must establish those keys, segmentation and coordinate frames from source
evidence. The joint route's separate RNA outputs do not establish that mapping.
Unresolved pLN section relationships cannot be repaired by routing code.

No cohort transforms, uncertainty bounds, comparison mappings or acceptance
records were created during this assessment.
