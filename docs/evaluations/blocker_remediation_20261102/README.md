# Portability blocker remediation — 2026-11-02

This follow-up implements the code repairs authorized after the original
[portability evaluation](../dynamic_portability_20261102/README.md). It preserves
that evaluation's source pins, failures, scientific outputs and 201-acquisition
roster. The new [user guide](../../PORTABILITY_REMEDIATION_USER_GUIDE.md) describes
current commands; the dated evaluation remains the record of its earlier code.

**Current accounting:** 201 acquisitions: **132 completed** within their stated
candidate/evaluation scope, **44 blocked**,
**20 excluded** by existing rules, and **5 unavailable**. These counts describe
computational evidence. No biological approval, registry approval or production
release activation is added. Both the full 25 mixed-protein and full 19 dedicated
pLN candidates passed their required computational audits and section inspections.

The [acquisition ledger](acquisitions.json) preserves every original row under
`baseline` and gives its current disposition separately. The [readiness matrix](readiness.json)
separates code, inputs, scientific specification, routing, computational checks
and biological review. An unchanged historical result is labelled as historical;
a source-code test is never substituted for a completed cohort run.

## Results and scope

| Route | Completed evidence | Current limit |
|---|---|---|
| CODEX pancreas, 20 images | Repaired content-bound replay: **60/60 exact assertions**, 3,199,606 eligible cells, 3,525,977 measured objects. **24/24** audited v20 content bindings recovered and reverified. | The external CODEX port remains separate; full v20 execution-input closure and biological acceptance are not established by these checks. |
| IO60 pancreas/embedded pLN protein, 25 acquisitions | Corrected normal publication, independent identity/count/primary-parity checks and inspection of all 25 map pairs pass. Final source comparison passes for 49 files. | Newly fitted refinement candidate; retained unavailable states and control-review limits. Stage06 is unaccepted. |
| Dedicated pLN protein, all 19 acquisitions | Fresh Stage03/05/06 candidate, **10,505,964 cells**, **9,883,811** nucleus eligible; exact identity/nucleus/coordinate/count reconciliation. | Raw QuPath review contract; no canonical own-compartment REDSEA equivalence or final pLN classes. Stage06 is unaccepted. |
| IO60 RNA, 65 acquisitions | Frozen full-cohort evidence retains **60 candidates and 5 unavailable**. Current wrapper applies validated naming policy, emits actual gene coverage/settings/fits, and remains independent of protein. | Full 65 was not rerun solely for default-preserving configuration wiring. Missing worker fits and biological review remain explicit. |
| Spleen PC67/PC33 protein | Historical PC67 eight-section candidate evidence retained. The ten previously unverified large inputs now match their declared hash algorithms; current full-file hashes are also recorded. | Four PC33 acquisitions remain blocked by missing reviewed lineage/registry specification. Two historical image fingerprints were sampled, so they do not prove historical full-byte identity. |
| Same-slide protein/RNA | Candidate linking preserves ties, alternatives, both assay calls, explicit mapped comparisons and supplied uncertainty. Batch coordinator plus matching regressions: **98 tests pass**; **12/12 pairs accounted**. | Eight spleen pairs lack their configured exact upstream receipts/transforms; four pancreas pairs lack a section adapter. No new accepted real-cohort geometry; pLN physical-section evidence remains unresolved. |
| Serial MALDI, 32 sections/18 donors | Original **1,666,450 spots** and mask/denominator audit retained. Xenium morphology route and receipt-bound object route implemented; **33 object/Stage08 tests pass**. Exact pLN manifest binds **19 acquisitions, 47 objects, 44 area geometries**. | Exact inspected tissue/neighbour declarations and accepted real Stage07 transforms are absent. No real Stage08 object result or same-cell evidence is claimed. |

## Mixed protein: all 25 completed

The [independent verification](../../../outputs/blocker_remediation_20261102/protein/verification/remediation_checks.json)
and [full-table audit](../../../outputs/blocker_remediation_20261102/protein/verification/release_audit/per_acquisition.json)
reconcile **39,394,852 segmented cells**: **32,874,241** nucleus eligible and
**6,520,611** ineligible. Embedded pLN contributes **1,767,936** segmented cells,
including **1,488,546** eligible cells. These are distinct from the
**31,385,695** eligible pancreas cells entering primary typing and refinement.
All 25 exact acquisition IDs match the declared inventory, and every primary
pancreas call matches its preserved baseline exactly.

| Mixed-table disposition | Cells |
|---|---:|
| Pancreas hierarchy | 31,385,695 |
| pLN evidence unavailable | 1,116,736 |
| Nucleus/QC unavailable | 6,520,611 |
| pLN no match | 161,063 |
| pLN unique candidate | 177,326 |
| pLN ambiguous candidate | 33,421 |
| Total segmented universe | 39,394,852 |

Within the pancreas refinement universe, **26,913,153** cells were evaluated,
**4,390,744** retain unavailable refinement and **81,798** are not applicable.
The per-section report preserves the specific missing defining-marker,
measurement, fallback-marker, fragment-measurement and fitted-bound reasons.
Unavailable refinement retains the primary call. Embedded pLN candidate states
do not become final cell classes.

All **253 input and 200 output byte-hash checks** pass. The **15 optional absent
pLN block inputs** have a separate, verified no-file disposition and remain in
the raw unavailable-hash audit; they are not successful byte hashes. The final
[source-stability check](../../../outputs/blocker_remediation_20261102/protein/producer_stability.json)
matches all 49 protein closure files against the actual attempt2 snapshot and
published runner manifest. The [inspection ledger](../../../outputs/blocker_remediation_20261102/protein/verification/inspection_ledger.json)
records inspection of all **25 call/availability map pairs** and **189 native-cell
records** from every observed review state. The [completion summary](../../../outputs/blocker_remediation_20261102/protein/completion_summary.json)
binds these results. Stage06 remains unaccepted.

The raw nucleus-admission audit used size/mtime binding; the current run verifies
the exact aligned/primary inputs and does not claim a new full raw-image/export
SHA audit. Existing coordinate-conversion differences remain documented in the
baseline audit, with no coordinate repair or registration acceptance. Seven
plasma-cell donor/rule cases remain unavailable in the immune component.
Retained inputs still resolve into the v3 checkout; the published gate requires
copying and verifying those retained inputs before an active release.

## Dedicated pLN: all 19 completed

The [completion summary](../../../outputs/blocker_remediation_20261102/config/completion_summary.json)
and [per-section audit](../../../outputs/blocker_remediation_20261102/config/pln_full19_validation/per_section.json)
bind the full candidate to its actual raw exports. All 19 consumed the fresh
shared 55-marker normalization fit. The run produced **1,045 section/marker
fits** and **1,311 section/rule evaluations**. Every section's spatial diagnostic
was inspected; 1,900 sampled cells provide **104,500 independently recomputed
marker-state checks**, all passing. All **204/204** upstream input/output hashes
matched. The 24 relevant pLN/UniFORM producer sources remained unchanged during
execution.

| Candidate disposition | Cells |
|---|---:|
| Unique candidate | 5,403,186 |
| Ambiguous | 683,914 |
| No match | 3,796,690 |
| Unavailable | 622,174 |
| Total | 10,505,964 |

The unavailable total preserves **622,153 nucleus-ineligible** cells and 21
eligible cells lacking sufficient candidate evidence. Fourteen rule rows require
markers absent from these exports; absent markers were not replaced by measured
zeros. The existing nonblocking normalization check flags **17 markers** above
its recorded spread bar. Those flags remain review evidence; no fit was tuned to
expected frequencies, disease labels or assay agreement. Final cell classes: **0**.

Scientific outputs are under
`/data/io60-current/evaluations/blocker_remediation_20261102/pln-full19-candidate-20261102`.
The [inspection ledger](../../../outputs/blocker_remediation_20261102/config/pln_full19_validation/inspection_ledger.json),
[receipt checks](../../../outputs/blocker_remediation_20261102/config/pln_full19_validation/receipt_checks.csv)
and [sample fidelity checks](../../../outputs/blocker_remediation_20261102/config/pln_full19_validation/sampled_evidence_checks.json)
remain available beside the overlays. Existing Stage03/05 computational gates
were automatic; no Stage06 acceptance was created.

## Protein retry provenance

The first full 25 attempt completed computation but the publication guard
[rejected it](../../../outputs/blocker_remediation_20261102/protein/stage-run.log)
because `refine.py` was read without being declared as an input. The declaration
has been repaired. The [recovery audit](../../../outputs/blocker_remediation_20261102/final_guard_review/protein_failed_attempt1_no_recovery_v1.json)
cannot establish the complete input contract for the old directory: 127
required inputs lack saved execution-time or historical content hashes, including
25 fundamental tables, 75 structure exports, 25 islet exports, donor metadata
and normalization acceptance. The 15 legitimate optional pLN block absences
remain a separate no-file category.

A fresh normal full-cohort run uses a new
[producer snapshot](../../../outputs/blocker_remediation_20261102/protein/producer_snapshot_attempt2.json).
Its [prelaunch input record](../../../outputs/blocker_remediation_20261102/protein/attempt2_input_bindings.json)
binds all 268 declarations: 253 existing inputs and 15 optional no-file
declarations. The failed directory and first snapshot remain unchanged.
The [current source check](../../../outputs/blocker_remediation_20261102/protein/producer_stability_attempt2_running.json)
matches all 49 files in the conservative protein import closure to the new
snapshot; this is a source comparison, not a completed runtime audit.
No receipt is fabricated and no old output is imported as a completed candidate;
all 25 acquisitions remained pending until the new run, independent audit and
per-section inspection passed.

The corrected attempt has now [published normally](../../../outputs/blocker_remediation_20261102/protein/releases/protein-availability-candidate/06_typing/receipt.json),
including all 25 mixed tables. Independent verification and the final source-stability
check pass, and all 25 sections have been inspected. The ledger now records this
computational candidate as completed while preserving the failed first attempt.

The [failure-recording review](../../../outputs/blocker_remediation_20261102/final_guard_review/protein_publication_failure_recording_review_v1.json)
binds a separate `protein_stage_publication` observation to the exact 25
acquisitions and six evidence artifacts. Its status is `failed`, with quality
unavailable: this records a publication-provenance failure, not a finding that
the scientific parameter values failed. The failed attempt has no complete outer
receipt. Only three failed metadata JSONs are newly hashed for the delivery;
the large failed scientific outputs remain unaudited and unpublished.

## Code repairs and verification

- **PROV-1:** content-based CODEX replay, exact reference decoding and runtime/thread
  validation; 10 focused fixture tests. The [CODEX report](codex.md) distinguishes
  the successful repair from preserved historical failures.
  A subsequently [received CODEX reference-context handoff](codex_received_handoff.md)
  has 17 verified payload files and an explicit adaptation specification. It was
  reviewed and not applied; its package pins do not replace the running producer
  or the original-method replay.
  The separate [component/Ki67 handoff](codex_received_component_handoff.md) is
  also received and not applied: 17 payload files and 214 regular sources were
  verified locally. Its scientific counts were supplied, not rerun; reported
  release/object failures remain, and remote delivery was not independently verified.
- **PROT-1:** acquisition/measurement availability now determines carried markers;
  unavailable refinement preserves primary calls and reasons. The full 25
  independent verifier, final source check and per-section inspection now pass.
  Representation fitting now includes carried-marker evidence in supported
  dimensions from four incomplete-panel sections previously discarded by
  union-column complete-case filtering. Fixed neighbour-count, concordance,
  bound and seed criteria are unchanged. Exact primary-call parity is required;
  refinement is newly fitted, with no claim of identical refinement labels to
  an unexecuted baseline or biological approval.
  The [first-attempt independent source review](../../../outputs/blocker_remediation_20261102/final_guard_review/protein_source_review.json)
  found no concrete defect across nine production/support sources and five test
  files; all nine producer hashes match the launch snapshot. This source review
  is separate from runtime completion.
  The [published producer warning](../../../outputs/blocker_remediation_20261102/protein/releases/protein-availability-candidate/06_typing/pancreas/build_cells_with_refinement.log)
  names `top_identity_marker`: it reads `FOXP3`
  on **4,487,012 cells (14.30%)** under the cohort scale. This diagnostic-field
  limitation is outside the availability repair and provides no Treg identity
  evidence. No scientific release acceptance is added.
- **CORE-1/CORE-2:** the existing protein cutoff and RNA naming settings now reach
  their consumers, reject unsupported values and prevent stale-setting reuse.
  Resolved parameter artifacts record fixed/selected/fitted/unavailable context
  and hash-bound inputs/sources. **39 focused tests** passed. This is instrumentation
  of the changed execution boundaries, not a claim that all historical literals
  have become dynamic settings or every worker fit is now emitted.
- **PANEL-3:** the existing dedicated pLN candidate contract has now run all 19
  sections; the previous two-section availability limit is removed for that contract.
- **GEOM-2:** [candidate-link repair](geom2.md) includes per-cell dispositions,
  ambiguity, independent calls, explicit hierarchy mappings and evidence-bound
  supplied uncertainty. Missing uncertainty stays unavailable; the existing
  **5 µm** cap is unchanged and is not a registration acceptance threshold.
- **GEOM-1:** the [batch coordinator](geom1.md) preserves independent pair roots,
  receipts and gates. The real preflight accounts for eight blocked spleen pairs
  and four unavailable pancreas adapters; unresolved pLN relationships remain
  separate. Twenty coordinator tests plus the 78 matching/crossmodal tests pass.
  A missing pair does not prevent another ready pair from proceeding, and no
  automatic review acceptance is created.
- **MALDI-2/3:** exact-region DAPI/calibration/frame inputs, source observation
  receipts, invalid-area availability, nonexclusive overlap summaries and finite
  denominators are checked. The [Stage08 review](../../../outputs/blocker_remediation_20261102/maldi_stage08_review.json)
  binds source hashes, 33 tests and visual inspection of the corrected overlay.
  The [MALDI implementation report](maldi.md) details source geometry, uncertainty
  and object-summary contracts.

The [spleen large-input audit](../../../outputs/blocker_remediation_20261102/hpg_maldi_xenium/spleen_algorithm_verification/verification.json)
read **207,864,852,220 current bytes** and verified **10/10 declared bindings**:
eight historical full geometry hashes and two historical sampled image
fingerprints (three 1 MiB windows). New full-file hashes of the images are recorded
separately; sampled history cannot establish their historical full-byte identity.
The first attempt compared incompatible hash algorithms and is retained as an
inapplicable comparison, not evidence of content drift.

The final indexed HiPerGator morphology validation accounts for **all 65 IO60
Xenium acquisitions: 62 validated, 3 unavailable because tissue was recorded
as damaged**. Both RNA-only acquisitions are included. The earlier indexed
60/65 pass rejected two explicitly relocated directories; the corrected reader
requires the exact acquisition/output index and matching experiment hash, then
records that relocation. The intermediate result is retained separately.
The earlier workstation pass found 30 available locally; its 35 unavailable
entries describe that machine's source access, not 35 damaged tissues. See the
[final source report](../../../outputs/blocker_remediation_20261102/hpg_maldi_xenium/results_relocated/source_validation.json).
No source-mask pass selects a MALDI neighbour or accepts a registration. These
source-morphology counts differ from RNA typing's five unavailable acquisitions,
which comprise three damaged regions and two low-RNA-signal regions.

The user subsequently supplied this serial-order context: **"For the most part,
the sections go first=xenium, second=phenocycler, third=maldi"**. The
[separate relationship-context record](relationship_context.json) preserves the
exact assertion and its qualifier; exceptions and exact acquisition/section IDs
were not enumerated. Within a confirmed three-section series following that
order, PhenoCycler immediately precedes MALDI. This context creates no exact
neighbour assignment or registration acceptance and does not overwrite the
inventory's separately recorded same-slide relationships.

All 62 source-mask panels were inspected, with flags and unapproved mask status
retained in the [source inspection record](../../../outputs/blocker_remediation_20261102/maldi_xenium_review.json).
Thirty locally available images match the corresponding remote full-file hashes.

The requested [MALDI scaling-plot viewer](../../../outputs/blocker_remediation_20261102/maldi_scaling/index.html)
shows all 32 sections side by side as linear intensity, `log1p` intensity and the
existing all-zero exclusion mask. The [plot manifest](../../../outputs/blocker_remediation_20261102/maldi_scaling/manifest.json)
binds both source exports, native micrometre coordinates and all 1,666,450 spots.
These are display transformations, with no intensity normalization or registration.
The existing mask excludes 12,486 all-zero spots; those spots remain visible in
the all-raster comparison.
The [final plot inspection](../../../outputs/blocker_remediation_20261102/maldi_scaling/inspection.json)
binds all 32 inspected figures and their final hashes.

The combined targeted regression checkpoint reports **413 passed, 9 deselected**,
with one known scikit-image FutureWarning. The
[log](../../../outputs/blocker_remediation_20261102/final_regression_revision6.log),
[XML](../../../outputs/blocker_remediation_20261102/final_regression_revision6.xml)
and [command record](../../../outputs/blocker_remediation_20261102/final_regression_revision6_command.json)
identify the run and tested source hashes, including the recorder's
primary/immune component evidence checks, input-declaration repair, fresh capture
command, explicit attempt bindings, publication-failure observation and precise
failed-attempt manifest scope. The
[source binding check](../../../outputs/blocker_remediation_20261102/final_test_source_binding.json)
bound 53 Python files to that tested revision. The subsequent
[source-stability guard check](../../../outputs/blocker_remediation_20261102/final_guard_review/protein_stability_guard_supplemental_v1.json)
passes 120 affected tests, including 21 new guard cases; these overlap the 413
checkpoint and must not be added to it. The final protein guard requires the
published runner manifest, actual attempt snapshot and current producer hashes.
The recorder archives this source-stability evidence with each protein observation.
The [current combined source binding](../../../outputs/blocker_remediation_20261102/final_current_test_source_binding.json)
matches 53 regression-tested Python files to the checkpoint or focused supplement,
including the four files updated by the final guard.
The separate evaluation-only review-materials helper is outside that 53-file
regression inventory. Its [CLI and Ruff checks](../../../outputs/blocker_remediation_20261102/protein/review-helper-checks.json)
and [actual full-cohort execution](../../../outputs/blocker_remediation_20261102/protein/review-helper-execution.json)
pass: nine unresampled contact sheets, 25 call/availability map pairs and 189
native-cell records were produced and separately inspected.
The earlier 223-, 298-, 307-, 339-, 360- and 398-pass
checkpoints are retained separately. The 14 recorded read-only
[CLI help/import checks](../../../outputs/blocker_remediation_20261102/cli_help_checks.json)
also pass; the review-materials helper has the separate workflow validation above.
Focused counts above overlap these suites and must not be summed as unique tests.

## What remains blocked

A blocker means a named required input or decision is absent; it does not mean
that all work stopped. Code defects with executable checks were repaired. The
remaining blockers require evidence the software cannot invent:

| Requirement | Evidence needed to proceed |
|---|---|
| MALDI registration/object results | Inspected tissue and exact neighbour acquisition for each section; existing source-frame/mask contract; receipt-bound registration review. All 32 remain in the ledger. |
| Real same-slide cohort links/comparison | Exact section, segmentation and frame provenance; accepted per-pair transforms; any uncertainty bound with its evidence; independent tissue calls and explicit hierarchy map for identity comparison. |
| IO60 same-slide adapter | Explicit mapping from mixed donor tables to the source section's cell universe, segmentation and crop coordinates. A joint typing stage does not supply this mapping. |
| PC33 and broader RNA tissue support | Reviewed panel identity cores/lineage rules and actual coverage; no PC67 or pancreas defaults substituted. |
| Broader CODEX port tissue hooks | Separate reviewed incorporation of the external port and active hook enforcement. The current package has no production `positivity` package. |
| Scientific transfer/release | Measurement equivalence, tissue/panel applicability and missing review criteria. Computational passes do not provide these approvals. |

Address missing physical source declarations and receipt-bound geometry first,
then the section adapter and missing panel/tissue specifications. The
[remaining panel/hook requirements](../../../outputs/blocker_remediation_20261102/config/remaining_panel_hook_requirements.json)
identify the actual source modules and unsupported contracts. The original
[remediation specifications](../dynamic_portability_20261102/remediation.md)
retain their reproduction evidence; the current readiness matrix identifies
which implementation portions have now been repaired. No unresolved numerical
review criterion is filled with a new cutoff.

Serial MALDI remains spot/object evidence. Zero means not detected; NaN remains
missing. Existing signal-mask denominators and explicit all-raster comparisons
remain separate. Lines are unavailable area objects, and annotation regions do
not become follicles, islets or cell types. No new distance bands are introduced.

## Evidence and operation

Use the [current user guide](../../PORTABILITY_REMEDIATION_USER_GUIDE.md) for isolated
candidate creation, route commands, review gates, diagnostics and recorder
retry/export. The [populated parameter history](parameter_history.json) contains
**159 observations and 159 records**, with **0 approved records**. The
[recording summary](../../../outputs/blocker_remediation_20261102/parameter-history/recording_summary.json)
confirms every retry was idempotent, automatic application and recommendations
are disabled, and global history remained absent and unchanged. The project
archive retains 58 distinct evidence files. Stage `resolved_parameters.json`
artifacts also preserve the actual boundary settings independently of database
availability.

The fresh [MCP protocol check](../../../outputs/blocker_remediation_20261102/mcp-protocol/protocol_check.json)
passes actual stdio record/read/export, idempotent retry, revision history and
rejection of approval, with global history unchanged. The remediation bundle
recorder also archives evidence and input receipts by hash inside its project,
so revisions retain earlier artifact bytes. The separate
[scientific MCP verification](../../../outputs/blocker_remediation_20261102/parameter-history/scientific_mcp_check.json)
passes actual protein-observation read, duplicate retry and exact export of all
159 records. Its export SHA256 remains
`32d889eaf2e2197f9e89691ee6390ed61cebe3548f229b5fef7ab98b83f92c7c`;
global history remains absent and unchanged.

The original delivery can still be verified with
`.venv/bin/python scripts/evaluation/verify_portability_delivery.py` in the documented
scientific environment. Evolved tracked code/document references resolve against
the explicit baseline Git blob where needed; scientific data and frozen evidence
remain byte-strict. New evidence is versioned under this remediation root and is
never written into the historical evaluation directories.
