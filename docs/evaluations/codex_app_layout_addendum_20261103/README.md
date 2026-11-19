# CODEX final layout and documentation addendum

Received and read the frozen code at `5e6580fe68af12dd18bee430988bc1ed64061526` and final documentation/evidence at `99bb82c87f4449c034d3c24276c7ebb2eb343e9d`. This separate [receipt](receipt.json) supplements the [earlier app/security review](../codex_app_fixes_receipt_20261103/README.md). It does not amend that review or the [completed portability delivery](../blocker_remediation_20261102/README.md).

| Commit | Reviewed change |
|---|---|
| `926c7ae` | Horizontal guide titles clear their gradients. |
| `1573a33` | OCI metadata moves after dependency layers; metadata values and non-metadata instructions remain equivalent. |
| `5e6580f` | Tissue legend text wraps; output and cached-render heights grow with the legend. |
| `99bb82c` | Final documentation, workflow guidance and compact evidence receipts; build source remains `5e6580f`. |

Only two Spatial rendering modules and the Containerfile differ from the earlier `b5fae5c` pin. The hierarchy helper and security patch/installer/verifier files remain unchanged. There was no source port, scientific recomputation, image build or deployment in this recipient review.

## What the evidence supports

Independent review passes 26 pure R checks and 10 synthetic JavaScript checks. All 73 tracked source hashes match the frozen code. Five final documents were captured directly from the documentation commit; 22 receipt, log, report and screenshot references match their supplied hashes. The [source review](../../../outputs/codex_app_layout_addendum_20261103/source_review.json), [saved UI review](../../../outputs/codex_app_layout_addendum_20261103/ui_review.json), and [document bindings](../../../outputs/codex_app_layout_addendum_20261103/final_document_bindings.json) preserve the scope and evidence.

Eight bound screenshots were visually inspected. Default and full 27-entry tissue legends are complete at both 1200 and 1600 pixel viewport widths. Reported tissue map heights are 406.82–482.54 pixels. Categorical UMAP legends and continuous colorbars sit below the plots, with their titles clear of the gradients. Final tissue captures use the final launch; the UMAP, Plot and Trajectory checks are reused from earlier checks of the unchanged relevant behavior.

The sender's DOM checks support stable output identity and two download links below the output. Brush/reset was exercised for the all-27, 1600-pixel case only. Download contents were not tested. Saved images cannot independently establish those interactions; no browser was relaunched here. The supplied expanded smoke result remains 45/49, with four Viewer/Range failures preserved. Container/public-dev acceptance remains separate.

The supplied logs record 35 parsed R files, 588 passing R expectations with zero test-framework failures/warnings/skips, and three passing Python tests. Those complete suites were not rerun for this addendum. The R log also contains a package-build-version warning outside the test summary, and the Python log contains an anndata/zarr deprecation warning; the logs are not warning-free.

## Open findings and adaptation requirements

Both 1200-pixel raw tissue images clip the rightmost `5000` x-axis tick at the image edge. The legends are complete. Any subsequent layout adjustment should retain the legend/map geometry and check the complete axis in both widths. This is a cosmetic caveat to the supplied visual acceptance, not a new scientific result.

The earlier mixed raw-Macrophage/profile merge defect remains unresolved: the complete `00_globals.R` blob is unchanged. Before generalizing that helper, calculate merged values from immutable original member columns and retain all destination families, with exact mixed-input tests for proportions, counts, enrichment and distance. The earlier probe does not establish a failure in the supplied current cohort, which contains profile labels only.

Carry forward the earlier requirements to preserve original identities and denominators, distinguish parent totals from additive leaves, retain ambiguity/unavailable evidence, and record the feature set for new trajectory fitting. Wrapped labels are display text, not renamed source identities. A future image needs its own installed-runtime and deployment evidence; historical security-image evidence does not certify a new image. Cache behavior was checked with pure tests, but cross-session performance was not benchmarked.

## Preservation and inspection

All eight [portability delivery checks](../../../outputs/codex_app_layout_addendum_20261103/prior_delivery_verification.json) still pass, including 2,973 fresh artifact hashes. The older app and component-resolution receipt hashes remain unchanged. The 201-acquisition ledger remains 132 completed, 44 blocked, 20 excluded and five unavailable. This addendum grants no biological approval and changes no parameter history or scientific acceptance.

The [artifact manifest](../../../outputs/codex_app_layout_addendum_20261103/artifact_manifest.json) binds local snapshots, screenshots and check outputs. The [documentation review](../../../outputs/codex_app_layout_addendum_20261103/documentation_review.json) records the final upstream documentation diff. Prior receipts stay immutable; any later source correction or deployment needs a new receipt. Acknowledgement delivery to the source task is reported separately from receipt creation.
