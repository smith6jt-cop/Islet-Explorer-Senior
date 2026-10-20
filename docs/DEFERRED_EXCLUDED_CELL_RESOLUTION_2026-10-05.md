# Deferred resolution of quarantined IO60 pancreas protein cells

Recorded 2026-10-05. Status: **deferred until an explicit future user instruction**.

All **11,484 cells** remain quarantined: **11,459 native Excluded cells and 25 EndocrineSingle exclusions**. The current user decision excludes every one of these cells from **all analyses, including numerators, denominators, sensitivity analyses, compositions, fitting, validation and graphs**. The three donor 6579 cells with reviewed broad immune evidence have no exception. Retain observed marker-pattern labels and image observations as diagnostic evidence, not confirmed identities or permission to admit a cell.

This record preserves evidence and prerequisites for possible future resolution. It authorizes no reclassification, large job, change to the active release or app development. **App development remains halted.** There is no scheduled resumption; wait for an explicit user instruction. A request to resume evidence review does not itself authorize admitting cells or changing the active release.

## Evidence and population at deferral

The [cell evidence table](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/complete_cell_typing_paper_style_20261005/source_tables/excluded_identity_evidence_20261005.parquet) and [review summary](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/complete_cell_typing_paper_style_20261005/source_tables/excluded_identity_review_summary.json) account for 11,484 rows across 25 donors. Direct table inspection found no missing or duplicate donor/cell keys, and `candidate_is_confirmed_identity=False` for every row. Native geometry QC is true for all review rows; this does not override quarantine or establish identity.

The native Excluded population comprises 8,783 cells in Islet objects and 2,676 in EndocrineCluster objects. The [native review description](/data/io60-calibration/excluded_cell_identity_review_20261005/README.txt) explains that these cells fell through the ordinary chain to Exocrine inside endocrine outlines and were excluded as annotation conflicts. The exclusion reason checked measurement presence, not positive epithelial identity. Refinement skipped these cells, so `refined_class=Excluded` is not independent identity evidence.

The mutually exclusive **diagnostic pattern summaries** below cover the 11,459 native Excluded cells only. They are not biological class counts.

| Observed summary | Cells |
|---|---:|
| Exocrine-like | 3,851 |
| Endothelial-like | 1,155 |
| Immune-like | 752 |
| Stromal-like | 410 |
| Neural-like | 25 |
| Endocrine-like | 82 |
| Epithelial, subtype unclear | 370 |
| Overlapping patterns | 323 |
| No defining pattern | 4,491 |
| Total native Excluded | 11,459 |

The separate nonexclusive candidate-family fields can overlap and must not be summed as exclusive identities. The summary reports 509 native cells with marker-reference corroboration, 44 with strong hormone support that still failed the native endocrine gate, and zero exact native immune leaves. None of these facts establishes a recovered identity.

The [Single exclusion evidence manifest](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/complete_cell_typing_paper_style_20261005/typing_audit/protein_analysis_exclusion_manifest.json) records the other 25 cells: 19 with all scaled hormones at floor zero, two with all hormones unavailable after fragment masking, and four with two hormones tied at the upper clipping limit. Preserve those distinct reasons. The manifest's older assay-denominator policy is superseded by the current decision above.

## Evidence locations and limits

| Evidence | Location and interpretation |
|---|---|
| Full native measurements and provenance | [excluded_cells_all.parquet](/data/io60-calibration/excluded_cell_identity_review_20261005/excluded_cells_all.parquet), donor partitions and [manifest.json](/data/io60-calibration/excluded_cell_identity_review_20261005/manifest.json). Retains native IDs, calls, object membership, raw/aligned/scaled measurements, hormone masks, overlap fractions and geometry QC. Raw pre-mask hormones remain distinct from unavailable typing measurements. |
| Exact rule terms | [existing_rule_dictionary.csv](/data/io60-calibration/excluded_cell_identity_review_20261005/existing_rule_dictionary.csv), [observed prerequisite support](/data/io60-calibration/excluded_cell_identity_review_20261005/existing_workflow_observed_prerequisite_term_support.csv), [exact rule hits](/data/io60-calibration/excluded_cell_identity_review_20261005/exact_existing_workflow_rule_hits.csv) and the versioned [phenotype workflow](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/phenotype_workflow.csv). Diagnostic rule hits do not replace native producer calls. |
| Marker similarity and fits | [similarity manifest](/data/io60-calibration/excluded_cell_identity_review_20261005/marker_similarity_manifest.json), [query results](/data/io60-calibration/excluded_cell_identity_review_20261005/excluded_marker_similarity.parquet), [donor validation](/data/io60-calibration/excluded_cell_identity_review_20261005/marker_similarity_donor_validation.csv), and each donor's `marker_similarity_models.json` / `marker_similarity_heldout_validation.parquet`. Reference agreement measures reproduction of producer labels. There is no stromal reference class; nearest-reference resemblance cannot supply one. |
| Independent consistency checks | [native evidence validation](/data/io60-calibration/excluded_cell_identity_review_20261005/review_evidence_validation.json) reports PASS over 17,393 checks. [Evidence release validation](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/complete_cell_typing_paper_style_20261005/typing_audit/excluded_identity_release_independent_validation.json) reports PASS over 468 checks and preservation of 361 source columns. These validate extraction, bookkeeping and the earlier interfaces; they do not establish biological truth or compliance with the newer universal exclusion decision. |
| Matched native images | [image_context](/data/io60-calibration/excluded_cell_identity_review_20261005/image_context), [image availability](/data/io60-calibration/excluded_cell_identity_review_20261005/image_context/image_review_feasibility.json) and [image export validation](/data/io60-calibration/excluded_cell_identity_review_20261005/image_context/native_image_export_validation.json). The saved audit reports 42 inspected examples: 32 native representatives, six broad immune candidates and four clipping ties, with native pixels and exact geometry comparisons passing. This is an existing sample review, not individual image validation of all 11,484 cells. |
| Geometry and localization | [representative visual notes](/data/io60-calibration/excluded_cell_identity_review_20261005/image_context/excluded_native_representatives/inspected_representative_visual_notes.json), [Single tie review](/data/io60-calibration/excluded_cell_identity_review_20261005/image_context/single_clipping_ties/single_tie_visual_review.json), [geometry reference distributions](/data/io60-calibration/excluded_cell_identity_review_20261005/geometry_DAPI_reference_distributions.csv) and [QC reason counts](/data/io60-calibration/excluded_cell_identity_review_20261005/exact_QC_estimation_reason_counts.csv). Review exact cell and nucleus shapes together with measured overlap and raw localized signal. |

Similarity fitting used training references only and held out complete object/spatial reference groups. Nevertheless, **7,916 query cells share an object or validation tile with training references**. The remaining 3,568 are not thereby biologically validated. The group holdout is not independent spatial validation of query identities, and the validation tiles are not biological distance bands. Nearest-reference similarity alone is insufficient for any recovery.

## The three donor 6579 cells

Retain these stable cell UUIDs, all from donor `6579`:

| Cell UUID | Recorded image evidence and uncertainty |
|---|---|
| `0d9016f9-861c-4761-83fc-d2fb1d08fa38` | Localized MPO-associated evidence compatible with a broad immune candidate. Neither that pattern nor nuclear morphology proves neutrophil identity. |
| `c7ee9a71-51f2-4bae-a0bc-ee040bc40115` | Localized CD11c/MPO-associated evidence compatible with a broad immune/myeloid candidate. No specific dendritic, macrophage or neutrophil identity or state is established. |
| `fe87e8ee-444e-49df-9730-425370463c35` | Localized immune-associated evidence compatible with a broad immune candidate. Dendritic versus macrophage identity is unresolved and no exact native leaf call exists. |

The [six-case broad immune visual review](/data/io60-calibration/excluded_cell_identity_review_20261005/image_context/broad_immune_similarity/broad_immune_visual_review.json) records these observations and competing explanations for other examples. Preserve the three cells' measured evidence without promoting it to a confirmed identity. All three remain excluded from every analysis and denominator, including sensitivity analyses.

## Historical fields and policy conflicts

At initial inspection, the evidence table retained `reviewed_identity=Immune` and `reviewed_immune_sensitivity_included=True` for those three cells. The inspected versions of the repository's [review note](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/CLAUDE.md), [atlas README](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/complete_cell_typing_paper_style_20261005/README.txt), [Single exclusion manifest](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/complete_cell_typing_paper_style_20261005/typing_audit/protein_analysis_exclusion_manifest.json) and earlier sensitivity validation described a broader assay denominator and/or inclusion of the three reviewed cells. **Those permissions are superseded by the current user decision.** Historical fields and documents are audit evidence only and must not be used to admit any of the quarantined cells.

During final verification, a concurrent update changed the evidence table and summary. The updated summary explicitly records `active_analysis_inclusions=0`, `reviewed_immune_sensitivity_inclusions=0`, the universal exclusion policy and this deferred chat ID, `01a10dd3-8d09-7fd0-90bc-c10290fe4dcf`. The updated table has the reviewed sensitivity flag false for every row. The evidence snapshot below reflects that update. This documentation task did not make those source changes or verify every downstream analysis output; future work must reconcile the exact versions it consumes.

The old total assay-QC population was 31,397,179; subtracting all 11,484 quarantined cells gives 31,385,695 for that same population before any additional analysis-specific filtering. This arithmetic is a reconciliation reference, not proof that any current output was rebuilt or that every analysis has that denominator. Future authorized release work must check every applicable numerator, denominator, fitted population and sensitivity table against the exact quarantined keys. No active-release reconciliation or regeneration is performed by this deferred record.

## Required audits before any recovery proposal

| Audit | Evidence required for a defensible resolution |
|---|---|
| Endocrine annotation conflicts | Reconcile the exact original cell and endocrine-object IDs, membership, original outline and failed routing reason. Establish whether the apparent conflict comes from the cell, neighboring tissue or annotation. A fallthrough label is not exocrine identity; a candidate marker family does not cure an annotation conflict. Preserve original records and describe any future proposed correction separately. |
| Observed and missing prerequisites | Record each required measured positive, a measured positive witness for each OR group, and every required measured negative. Missing OR alternatives do not defeat a valid observed witness; a missing negative does not establish absence. Preserve conflicting measurements and exact saved calibration/rule provenance. The existing rule threshold is 0.10; this record introduces no positivity threshold. The saved 0.8715 epithelial immune-veto level is not a general epithelial identity cutoff. |
| Fragment-masked hormones | Keep the mask and unavailable typing values. Raw or aligned pre-mask brightness may explain the conflict but cannot be substituted as usable hormone evidence. Resolution requires attributable cell-body evidence and any measurement correction to be independently justified; never fill unavailable values with zero. |
| Cell body versus passing process | Show the native image with exact cell/nucleus geometry and the relevant structural object. Preserve centroid membership, whole-cell area coverage and creating-marker evidence. Structural identity requires the established centroid, majority-cell-area and creating-marker conditions together; overlap alone, a crossing vessel/nerve, or a touching immune process is insufficient. Keep continuous overlap fractions; introduce no expansion ring or replacement distance band. |
| Floor-zero and upper-clipped hormones | Distinguish the 19 floor-zero cells from two masked cells and four upper-clipping ties. Audit raw, aligned and scaled values, feature floors, clipping and cell compartments. Do not turn an all-zero argmax or a tie into a subtype, choose a convenient preclip winner, or invent a cutoff. The four saved tie reviews do not resolve subtype; one has identical cell and nucleus geometry and no nonnuclear cell area. |
| Geometry and segmentation | Establish a genuine nucleated cell with trustworthy boundaries and signal attribution. Audit nucleus/cell area, estimation flags, neighboring bright structures, segmentation merges/fragments and exact image registration. A cell without a nucleus remains excluded from every calculation. Passing historical geometry QC alone is insufficient to settle identity. |
| Independent identity evidence | Require cell-specific independent image evidence or correctly registered orthogonal measurements. Record source, channels, transformation/registration quality and exact cell correspondence. Adjacent sections, unverified same-cell matches, display brightness or agreement with protein-trained reference labels do not establish identity. Conflicting evidence remains unresolved. |

Keep the established Cluster/Islet endocrine rules: donor/size-stratum CD99 and its floor guard; the permitted large-stratum nucleus-elongation route; and a measured hormone before naming a subtype. Do not substitute an ordinary-tissue reference for an unavailable stratum gate. The EndocrineSingle unique-highest-usable-hormone exception does not relax Cluster/Islet rules or rescue missing/tied Single measurements. Audit against the saved producer rules and [typing reconstruction](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/core/phenocycler/scripts/typing/RECONSTRUCTION.md), retaining their provenance.

Never use disease ordering, expected class frequencies, donor uniformity or historical annotations as identity targets. Do not infer unsupported acinar/ductal or fine immune subtypes from this panel. Leave overlapping patterns and remaining uncertainty explicit.

## Prerequisites for resuming deferred work

1. Obtain an explicit future user instruction naming the authorized scope. Record that instruction and date. No elapsed time, new similarity score or changed file automatically resumes work.
2. Verify the evidence snapshots below and resolve any version drift. Confirm availability of the exact native images, channel mappings, cell/nucleus geometry, original object membership, native inputs and frozen calibration/rule provenance. Saved availability reports are evidence from their creation time, not a substitute for a fresh check when work resumes.
3. Preserve stable `donor_id`, `cell_id`, QuPath `object_id` and `endocrine_object_id`. Join through recorded keys, never row order, centroid guessing or donor-name parsing. The evidence table has `donor_id=112` for 159 cells whose native ID is `0112`; its `donor_id_native` preserves the exact source spelling. Direct inspection confirms exact native donor/cell/object key-set equality using that field. Retain both representations and the established mapping, without renaming the source IDs.
4. Prepare an isolated, versioned case ledger linked to the unchanged evidence. Each case needs original and analysis donor IDs; cell and object UUIDs; quarantine origin; native calls; diagnostic patterns; measured/missing prerequisites; masks and clipping; geometry/registration evidence; source hashes; competing explanations; remaining uncertainty; reviewer/date; and proposed resolution with evidence links. Keep diagnostic pattern, proposed identity and analysis eligibility separate. All entries begin and remain excluded until explicitly resolved and authorized for admission.
5. Before a proposal to admit any cell, complete every applicable audit above and document why its original exclusion is resolved without relaxing the rules. Record identity specificity no finer than the independent evidence supports. Cases lacking adequate evidence remain quarantined; there is no requirement to recover a target number or fraction.
6. Present the exact proposed cell keys, supporting evidence, unresolved limitations and downstream count/denominator effects for explicit user authorization before admission or active-release changes. Only an authorized future implementation may update versioned eligibility/assignment outputs, reconcile all analyses and validate absence of every still-quarantined key. Original measurements, identifiers and exclusion history remain preserved. App work requires a separate explicit instruction to resume it.

## Evidence snapshots

SHA-256 values independently read at final verification of this record:

| File | SHA-256 |
|---|---|
| Cell evidence Parquet linked above | `dbae50f2c6f4fb57e492d8637ce7aa9ec5b7017b04e63655dd36d3b3cf18ded6` |
| Review summary JSON linked above | `e3007f2a284ea6c5b69d30e95d00a1c1964f7042ce97a74a6dcc09991c84ad01` |
| Full native extraction Parquet linked above | `a160b98fec28824f6c39fb930d3dbde7e121264539dca116b6cabb59e78a6794` |
| Native review manifest | `49a52aec811fd8fb66d8977635ff2b1d1f971308030d67c3971da7a0516c2525` |
| Similarity manifest | `ae28fd6993740604f148062452053fc4f604473f782141f7418d1d287ba55807` |
| Native evidence validation | `48930ebf3a1b54af6332527f8b293c60f9c2d4bb28ef815dd370226de818f1e6` |
| Evidence release validation | `6bcfd428392a4409fbf0ed646bcbad126db3e2cf9ff27d8c856287cc5186d4d7` |
| Image export validation | `1414cace65d49f51c5d87d39e8f4ff4a11bb08070801080412857b1b0f9a7d50` |
| Broad immune visual review | `573509aad135201907c1ae23afadeabe2ccb406e2de69db6954f2d249e3c1b4c` |

The initial evidence table hash was `52761156b1e39a91134ac118762d8d4649703780c585e097de8a4e5e3220f507`, and the initial summary hash was `4d0a4cc11e4f62ced3b96811fd5ef3cc4b3ba688be817c883c40bf16765abb6f`. These identify the earlier policy snapshot, not the updated files. The saved evidence-release validation above refers to that initial evidence-table hash and has not been rerun by this task.

The source manifests retain deeper native, calibration, geometry and image provenance. This record preserves the current snapshot without copying or rewriting the evidence. Historical PASS results describe their original checks, including now-superseded sensitivity policy, and must not be treated as authorization for recovery.
