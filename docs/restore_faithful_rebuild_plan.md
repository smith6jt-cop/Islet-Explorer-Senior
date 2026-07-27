# RESTORE-faithful broad-lineage rebuild

> **This is the live gate plan for the phenotyping rebuild** (Gates 1–7). Committed
> 2026-07-22 from the Copilot session that produced it, where it existed only as
> `~/.copilot/session-state/7cb68b00-bc52-4fd8-b1cb-c804049aeb2b/plan.md` — a disposable
> location. The implementing code is the `Phenocycler_Analysis` submodule
> (`phenocycler/restore_validation.py`, `restore_pair_review.py`,
> `restore_validation_figures.py`, `cohort.py`, commit `6d49ae7`); the review bundles and
> screens it references now live under gitignored `data/restore_pair_validation/`.
>
> **Status: Gate 2 of 7.** No marker pair is accepted, no divisor is emitted, and no
> normalized, phenotype, or app output exists. As of 2026-07-24 one reference is FROZEN per target
> under the maintainer's reference≥target frequency rule (see the dated log entry and
> `data/restore_pair_validation/reference_selection_v1/frozen_reference_matrix.csv`); the active
> screen is `nnmf_all20_expanded_v10/` (METHOD_VERSION v8). The remaining Gate-2 work is the
> per-image biological review of the accepted pairs (mutual exclusivity / spatial coherence) in
> `data/restore_pair_validation/expanded_review_v11/review_decisions.csv` before Gate 3.
>
> `## Current status` below is an accumulating dated log — later entries supersede earlier
> ones. For what is active RIGHT NOW, read its last entry and
> `data/restore_pair_validation/README.txt`.

## Current status

- Full raw/REDSEA-only purge approved by the user.
- Purge completed for 17 named RESTORE, phenotype, QuPath, islet, trajectory,
  clustering, parquet, and deployed-app paths (about 9 GB).
- `data/` now contains only `cells/`, `cells_redsea/`, `redsea_scratch/`,
  `_redsea_validate/`, and `singlecell_protein.h5ad`.
- The retained H5AD has raw marker values and no phenotype/cluster embeddings, but it
  contains 23,280,387 cells and predates the 22-donor parquet rebuild; treat it as a
  historical raw input, not the canonical current cohort.
- Raw/REDSEA key reconciliation is complete:
  - root cause: REDSEA must rasterize every GeoJSON object for contact correction,
    whereas canonical `data/cells` excludes degenerate objects below 5 um2
  - REDSEA now retains all GeoJSON objects during correction but persists only
    canonical raw `object_id` values
  - all 22 retained REDSEA parquets were atomically reconciled
  - raw and REDSEA now contain the same 34,224,525 unique, non-null IDs with no
    donor-level mismatches or duplicates
  - all nine markers needed by the proposed eight target/reference pairs are present
    in every donor
  - full schemas retain the known physical panel difference: batch 1 has `CD38`,
    batch 2 may have `b_Catenin1`; neither marker participates in first-pass lineage
  - the invariant and reconciliation utility are covered by the submodule test suite
    (63 tests passing)
- First accepted taxonomy: `Immune`, `Epithelial`, `Endothelial`, `Neural`,
  `Mesenchymal`, and `Other`. Endocrine/exocrine is deferred to a separately
  validated second pass.
- Canonical background divisor: maximum target intensity in the accepted
  target-negative controls, as used by the manuscript equations. Released-code
  `mean + 3*SD` is sensitivity analysis only.
- Invalid or absent target/reference pair: fail the marker/image and redesign the pair
  before continuing. Do not use cohort imputation, adaptive per-donor reference
  selection, or an automatic fallback.
- The QuPath project is available at `/home/smith6jt/IO60panc2nd/project.qpproj`.
  Its saved detections contain nucleus, cytoplasm, membrane, and whole-cell
  measurements for all channels (mean/min/max/median/SD). Whole-cell means match the
  canonical raw parquet exactly by `object_id`; QuPath-only detections are the already
  understood sub-5 um2 objects and must be excluded downstream.
- NNMF is now a primary method candidate rather than SSC being assumed:
  - the manuscript permits biologically known pairs or data-driven pair discovery and
    permits NNMF, SSC, GMM, or another constrained separator for cell-group inference
  - two-column NNMF is the manuscript's trivial `W=Y, H=I` case
  - compartment-rich means provide a genuine higher-dimensional input without
    inventing area-derived totals
  - donor 6539 demonstrates clean two-factor separation for E-cadherin/CD31,
    B3TUBB/EpCAM, and Vimentin/E-cadherin, but sparse immune pairs need further
    component-count, missing-data, and control-purity validation
- Three-donor NNMF pilot completed on donors 6539, 6457, and 6476:
  - deterministic 5% QuPath sample matched canonical raw IDs exactly and retained one
    image per donor; whole-cell means matched canonical raw parquet within 0.00005 MFI
  - 240 model runs covered all eight pairs, three NNMF feature sets, two/three
    components, plus raw/REDSEA GMM comparators; no RESTORE or phenotype outputs were
    produced
  - membrane+cell features had complete coverage; non-nuclear coverage was 87-91%;
    all-compartment coverage fell to 40% for donor 6457, so adding nuclear features is
    not an automatic improvement
  - two-component membrane+cell/non-nuclear NNMF was strong and converged for reciprocal
    E-cadherin/CD31 and Vimentin/E-cadherin in all donors
  - B3TUBB/EpCAM was directionally separated by means in all fits but failed the
    median-direction criterion for donor 6457; retain as promising, not accepted
  - the immune pairs are not yet valid under this NNMF formulation: CD11b failed donor
    6457 and was component-sensitive; CD3e had only ~1.08-fold target separation;
    CD20 largely failed; CD68 did not converge in any two-component membrane/non-nuclear
    fit
  - canonical pilot tables are under session
    `files/restore_pair_validation/nnmf_pilot_verified/`
- Manuscript-style pilot figures completed and visually inspected:
  - corrected an initial display-only transposition caused by mapping the internal
    `(target, reference)` tuple directly to `(x, y)`; all pair clouds now follow
    RESTORE Figure 7 exactly: reference marker \(i\) on x and target marker \(j\)
    on y, with thresholds and marginal densities on the corresponding axes
  - superseded the initial log-transformed, unfiltered plots with raw-linear figures
    using explicit input QC before NNMF: estimate a per-donor, per-marker raw-MFI
    Otsu floor (upper 0.5% clipped only while estimating the floor), then retain a
    cell if either marker exceeds its floor; this removes the double-low lower-left
    population while preserving both mutually exclusive axis arms
  - after identifying severe B3TUBB/EpCAM arm imbalance and CD20/E-cadherin
    range dominance:
    - NNMF training now uses deterministic equal-size provisional target-only and
      reference-only arm samples; double-high cells are excluded from fitting, but all
      retained cells are assigned and actual arm counts remain reported
    - each membrane/whole-cell feature is robust min-max scaled using the balanced
      training sample's 0.5th and 99.5th percentiles; values outside those bounds are
      clipped to [0,1]
    - plots show the corresponding min-max-scaled whole-cell values on linear axes
      without a log transform and use equal-N component point samples while annotating
      actual group sizes
    - background maxima now use 200 repeated equal-N arm subsamples (smaller arm size
      capped at 2,000); plots show the median maximum and 5th-95th percentile band,
      with full-group maxima retained in the metrics CSV
    - B3TUBB's inferred target arm is only 0-1% of retained cells in these donors;
      the wide resampled-maximum bands now expose that threshold uncertainty rather
      than presenting a size-biased point estimate
  - added a donor-local low/no-target-expression check grounded in the manuscript's
    high-SBR requirement: target-arm mean intensity must be at least 2-fold above the
    inferred target-negative reference arm; this is not a lineage-frequency rule
    - confirmed donor 6539 CD20 is flagged `LOW/NO` (1.32-fold), consistent with the
      maintainer's image review
    - donor 6476 CD20 is also low-support (1.51-fold); donor 6457 remains supported
      (3.10-fold)
    - CD3e is likewise low-support in donors 6539 and 6476 under this diagnostic
  - `figure1_pair_cloud_examples.{png,pdf}` compares representative pairs across all
    three donors
    - removed the gray double-low input-QC rectangles from all marker-cloud panels
      without changing the underlying double-low exclusion
    - removed repeated `Dmed` and `conv` fields from the example subpanel annotations
    - updated the legend to distinguish dashed median equal-N sampled maxima from pale
      5th-95th percentile bands
    - green uncertainty intervals are computed independently for x and y and are drawn
      only when they retain non-zero width after min-max clipping; a band may therefore
      appear on only one axis when the other arm is sampled in full on every repeat or
      its interval collapses at the 0/1 display boundary
    - explanatory footers are wrapped inside the lateral span of the subplots rather
      than extending across the full figure canvas
  - `figureS1_all_pair_clouds.pdf` contains eight detailed pair pages with marginal
    component densities, NNMF arms, and illustrative equal-N sampled-maximum intervals
    versus mean+3SD background lines
  - `figure2_input_qc_summary.{png,pdf}` summarizes recomputed directional splits,
    retained fractions, and double-high fractions after input QC
  - individual pair PNGs, panel metrics, and a provenance manifest are in
    `files/restore_pair_validation/nnmf_pilot_verified/figures/`
  - all 91 submodule tests pass; syntax/whitespace checks pass; figure generation did
    not recreate any RESTORE-normalized, phenotype, or app-data directories
  - min-max-scaled linear plots show clear L-shaped geometry for reciprocal
    E-cadherin/CD31, Vimentin/E-cadherin, and B3TUBB/EpCAM; CD20 is supported only in
    donor 6457, CD68 remains weak, and CD3e/CD11b are donor-dependent with substantial
    double-high contamination in some donors, so no immune pair is accepted yet
- Gate 1 executable specification is complete and Gate 2 remains in progress:
  - the locked all-donor expanded screen contains 660 evaluations: 22 donors x 30
    fixed candidate pairs
  - computational separation never implies biological acceptance; all otherwise
    successful pairs remain `PAIR_REVIEW_PENDING` until a cohort-level pair review is
    explicitly recorded
  - small all-NaN REDSEA subsets in eight donors are tracked as cell-level
    `Unavailable`, not converted to biological negatives or image-wide failures
  - no marker pair has been accepted and no normalized RESTORE output has been created
- A deterministic Gate 2 review bundle is complete under
  `files/restore_pair_validation/expanded_review_v1/`:
  - 13 predeclared shortlisted fixed pairs, 77 donor-pair stress cases, 13 spatial-review
    PDFs, and 77 capped QuPath review CSVs
  - queue composition is 42 `PAIR_REVIEW_PENDING`, 20 `LOW_SIGNAL_UNRESOLVED`,
    11 `NO_REFERENCE_POPULATION`, 3 `MODEL_UNSTABLE`, and 1 `INVALID_PAIR`
  - case selection includes donor 6476, weakest target/reference folds, weakest
    seed stability, highest double-high overlap, greatest sampled/full-maximum
    divergence, a median-fold control, donor 6539 for CD20, and the sparsest B3TUBB arm
  - selection is deterministic, including exact median ties, and is a review-sampling
    rule only; it is never a cell-frequency acceptance rule
  - repeated generation produced byte-identical CSVs and manifests
  - `scripts/groovy/import_restore_pair_review.groovy` imports one overlay into the
    current QuPath image by canonical `object_id`; it preserves existing PathClasses by
    default and only assigns temporary review classes when explicitly enabled
  - all 103 submodule tests pass, including deterministic selection, category
    precedence, shortlist completeness, and atomic cleanup on write failure
  - the retained canonical `data/` still contains only raw/REDSEA inputs and scratch
    material; no forbidden RESTORE, phenotype, clustering, trajectory, or app output
    was recreated
- Initial figure inspection confirms the intended stress behavior but is not pair
  acceptance:
  - reciprocal structural anchors retain strong L-shaped separation and stable fits
  - B3TUBB/EpCAM remains computationally clean but requires image-level exclusion of
    endocrine TUBB3 before it can represent Neural
  - immune candidates remain heterogeneous: the review queue deliberately exposes
    low-signal, missing-reference, unstable, and double-high cases rather than hiding
    them behind cohort summaries
  - donor 6539 CD20 remains an explicit confirmed-absence case and must produce no
    positive calls in any later pilot
- The version-1 evaluations and review PDFs are superseded and must not be used:
  - review showed that the plotted maxima were calculated from entire NNMF components,
    allowing double-high and low-confidence assignments outside the intended axis arms
    to set the green lines
  - this was a scientific control-definition defect, not merely a legend or display
    issue
- Corrected method version `restore-pair-validation-v2` is complete:
  - Step 1 intersects the NNMF target component with the raw target-high/reference-low
    arm; its full maximum REDSEA reference intensity is the dashed vertical baseline
  - Step 2 requires NNMF reference-component cells to be raw
    reference-high/target-low and strictly right of Step 1; only those controls set the
    solid horizontal target-divisor candidate
  - double-high, double-low, and all other assigned cells set neither maximum
  - all 660 donor-pair evaluations were rerun; every nonempty inferred target population
    is bounded by Step 1, every reference control is right of Step 1, and no control count
    exceeds its candidate count
  - compared with version 1, the median vertical boundary is 0.48x and the median
    horizontal candidate divisor is 0.35x, confirming that the correction is material
  - no pair is accepted and no canonical divisor or normalized call was emitted
- Replacement review bundle is complete under
  `files/restore_pair_validation/expanded_review_v2/`:
  - 13 valid PDFs and 73 donor-pair QuPath overlays
  - a ready-to-run `import_restore_pair_review.groovy` is now included at the bundle
    root, recorded in the manifest, and preconfigured for the version-2 CSV directory;
    future bundles render their own absolute review path automatically
  - after a cross-donor import attempt exposed a misleading hardcoded 6539 CSV, the
    importer now detects the open image and presents only matching pairs in a QuPath
    dialog; it rejects image/UUID mismatches before mutation and applies temporary
    classes only when it will not overwrite unrelated PathClasses
  - non-sampled unclassified detections now receive a neutral, temporary
    `RESTORE review - Not sampled` class so the background can be toggled off; later
    overlay imports replace prior RESTORE review classes
  - pre-existing QuPath classes no longer block review: before the first overlay, the
    importer stores each original PathClass in persistent per-object metadata; the Clear
    dialog action restores every original class exactly, including hierarchical classes
    and Unclassified
  - `START_HERE_RESTORE_REVIEW.txt` is now the single authoritative operating guide
  - `review_decisions.csv` is the only reviewer-edited surface: 73 rows in the exact
    order of the 13 PDFs, with only Step 1, Step 2, spatial, disposition, and notes fields
  - QuPath is now an escalation step for unclear PDF rows, not a mandatory parallel
    workflow; the same importer dialog also clears temporary review classes
  - fixed a QuPath/Groovy 5 runtime incompatibility in the dialog option preparation:
    removed unsupported Java `ArrayList <=> ArrayList` sorting and sort only the final
    string labels; the reported CD3e/Vimentin case was reproduced successfully
  - every PDF legend explicitly defines the dashed Step 1 vertical line and solid Step 2
    horizontal line
  - each spatial row now shows the whole tissue with a zoom box and a deterministic
    representative-area zoom
  - representative structural, neural, and immune PDFs were rendered and visually
    inspected
  - the reversible class round-trip was validated directly in QuPath's Java 25 runtime
  - all 108 submodule tests pass; the two existing multiprocessing fork warnings remain
- Version 3 supersedes version 2 after the maintainer approved retention of every
  reference-negative NNMF target-component cell:
  - the high-confidence raw target-high/reference-low subset remains only the Step 1
    anchor; its maximum REDSEA reference MFI sets the vertical separator
  - lineage evidence is `normalized_target > 1 OR (target component AND reference REDSEA
    <= Step 1 separator)`; reference-positive target-component cells are not rescued
  - donors 6457 and 6579 are excluded until further notice from every analysis:
    6457 because poor DAPI staining prevented accurate channel registration, and 6579
    by maintainer decision as a pancreatitis outlier
  - every artifact containing 6457 or 6579 is quarantined under
    `superseded_includes_6457_2026-07-20/` or
    `superseded_includes_6579_2026-07-21/`; only raw/REDSEA provenance may remain
  - `qupath_compartment_sample_20donor/` is the only active validation input and contains
    exactly the 20 eligible donors
  - the v3 630-row screen and 77-row review are no longer active: maintainer review found
    that the universal raw-linear Otsu arm floor substantially over-called CD3e and CD20,
    including the purported high-confidence anchors
  - no pair or production divisor is accepted, and invalid/review-pending evaluations
    emit no calls
  - the new metrics explicitly count target anchors, supported target-component cells,
    additional supported cells, and supported cells retained at/below the divisor
  - version 4 superseded the universal-Otsu version 3 method by locking marker-specific
    raw arm floors:
    - CD3e and CD20 use `max(raw-linear Otsu, raw-linear Triangle)`
    - all 12 screened markers have an explicit policy entry; there is no silent fallback
    - every evaluation records the selected method, Otsu component, Triangle component,
      and final floor for both target and reference markers
  - version 5 is active after visual review exposed an ordering defect in v4 target
    support: NNMF target assignment could overwrite input-QC double-low status
    - across the 72 v6 review rows, 1,479,719 of 1,490,683 modeled lower-confidence
      cells (99.26%) were actually below both raw floors
    - v5 establishes the explicit order: required-measurement validity, raw double-low
      QC, NNMF assignment, Step 1 anchor/separator, then target-component support only
      among input-QC-retained cells
    - the 600-row rerun removed 11,580,847 double-low rescues while every raw floor,
      target anchor, reference control, Step 1 separator, candidate divisor, model state,
      and queue selection remained exactly unchanged
  - `nnmf_all20_expanded_v5/` contains the then-active 600-row screen: 20 eligible donors x
    30 fixed candidate pairs — SUPERSEDED 2026-07-23 by `nnmf_all20_expanded_v8/` (400 rows,
    20 pairs) when CD20 was retired as a reference; see the final entry of this section
    - state counts are 543 `PAIR_REVIEW_PENDING`, 30 `LOW_SIGNAL_UNRESOLVED`,
      20 `MODEL_UNSTABLE`, 5 `NO_TARGET_EXPRESSION_CONFIRMED`, and 2 `INVALID_PAIR`
    - all ordered-control and cell-accounting invariants pass
    - all 240 rows that involve neither CD3e nor CD20 are numerically identical to v3
    - CD20's floor increased above Otsu in all 20 donors; CD3e's increased in 10 of 20
    - donor 6539 CD20 produces no calls: five pairs are confirmed absent and two are
      invalid because their NNMF target signal contradicts the confirmed absence
    - no canonical divisor or production call was emitted because no pair is accepted
  - `expanded_review_v7/` was then the active review bundle: 13 PDFs, 72 checklist rows, and
    72 QuPath overlays; the checklist has zero reviewer edits and neither excluded donor
    appears — SUPERSEDED 2026-07-23 by `expanded_review_v10/`; v7, v8 and v9 all
    carried zero reviewer edits, so nothing recorded was lost
    - donor 6579 is the pancreatitis donor and remains excluded by maintainer decision
    - donor 6476 is not a pancreatitis donor; its obsolete fixed stress-case selection
      was removed, while nine evidence-driven 6476 stress rows remain
    - every canonical cell inside every representative crop is projected through the
      locked sample-fitted model and labeled, with no display or category subsampling
    - the 72 CSVs contain exactly 1,688,866 crop cells (2,704-48,405 per crop); every
      CSV row count matches its manifest count and every modeling-sample overlap retains
      exactly the original category
    - zero double-low cells are target-component supported; v7 has 6,699 lower-confidence
      crop cells versus 2,445,457 in v6
    - representative crops use a 10%-per-axis tissue window instead of 22%, producing a
      2.2x tighter linear close-up while still labeling every canonical crop cell
    - the overview inset retains the whole-tissue crop rectangle but has no title text
  - `expanded_review_v6/` is retained only as the superseded wide-crop bundle with the
    incorrect double-low rescue; `expanded_review_v5/` is the prior sampled/capped
    marker-image layout and `expanded_review_v4/` is the earlier point-scatter layout
  - each PDF uses an expression cloud plus four aligned views of one deterministic 2.2:1
    QPTIFF close-up: raw target pixels, raw reference pixels, review overlay, and merged
    target/reference pixels with classification rings; the overlay embeds whole-tissue
    context and the crop rectangle
  - marker crops decode only overlapping QPTIFF tiles and record the source image,
    pyramid level, physical extent, channel names, and display windows per row; the
    positive-pixel 1st-to-99.8th percentile window and gamma 0.65 affect display only
  - legends distinguish the vermillion high-confidence Step 1 anchor (`#D55E00`) from
    purple lower-confidence reference-negative target-component cells (`#6A3D9A`) and
    state that Step 2 is a normalization divisor, not their exclusion boundary
  - the final PDFs use 20 pt main titles, 12-14 pt row/axis labels, 11 pt ticks and
    legends, and a dedicated 3.2-inch legend-only band below all data axes; representative
    epithelial, immune, and stromal PDFs were visually inspected with no overlaps
  - the donor-specific anchor/target count summary is in the title area above each
    expression scatter, never inside the data axis; a 2.25-inch page-title band separates
    the first row summary from the main heading
  - the QuPath importer uses `RESTORE review - High-confidence target` and
    `RESTORE review - Lower-confidence target`; it labels every crop cell and leaves
    detections outside the representative crop at their existing PathClasses; the bundle
    copy and primary runnable copy are byte-identical
  - central cohort guards filter both 6457 and 6579 from automatic discovery and reject
    explicit analysis requests containing either donor
  - `marker_floor_audit_v1/` compares candidate statistics numerically and
    `marker_floor_review_v1/` contains five-page donor-local Otsu/Triangle/Yen PDFs for
    CD3e and CD20; neither excluded donor appears
  - CD31, E-cadherin, and Vimentin retain raw-linear Otsu as their reviewer-validated
    candidate statistic; CD3e and CD20 use the approved conservative Otsu/Triangle maximum
  - all 121 package tests pass; the two existing multiprocessing fork warnings remain

- 2026-07-23 — CD20 retired as a REFERENCE; this entry supersedes every earlier "active"
  designation in this section:
  - a reference marker's positive cells ARE the target's negative control, so a reference must
    label a population large enough to estimate a background from. CD20 cannot: its exclusive
    arm ran 41-2,026 cells per donor against 2,206-50,000 for the structural markers, it was the
    limiting arm in every degenerate screen row, and it produced inverted fits in which an
    abundant target was thresholded against a handful of reference cells — `6591 CD11b <- CD20`
    used 50 reference cells to threshold 30,507 target cells (610:1)
  - measured on the v6 screen, the three CD20-referenced pairs held 60 rows, 15 of them inverted
    fits, 4 of the 22 degenerate-bound rows, and 4 `MODEL_UNSTABLE` rows
  - substituting CD3e COLLAPSES those pairs rather than replacing them: `CD68 <- CD3e` and
    `CD11b <- CD3e` already existed and `CD3e <- CD20` would become a self-pair. So
    `IMMUNE_REFERENCE_SCREEN_PAIRS` drops 10 -> 7 and `CANDIDATE_PAIRS` 30 -> 27, and the screen
    600 -> 540 rows. CD20 remains a TARGET against CD3e/CD68/CD163 and the epithelial references
  - `INVALID_REFERENCE_MARKERS` records the rule with its rationale and is enforced at import,
    together with a self-pair guard, so a future pair set cannot reintroduce either quietly.
    `SHORTLIST_PAIRS` is checked against `CANDIDATE_PAIRS` so the two lists cannot drift apart
  - all 540 surviving rows kept divisors identical to v6, so nothing else moved.
    `MODEL_UNSTABLE` 18 -> 14 and degenerate-bound rows 22 -> 18
  - `nnmf_all20_expanded_v7/` was then the active screen — SUPERSEDED later the same day by
    `nnmf_all20_expanded_v8/`, see the final entry of this section: 540 rows (20 eligible donors x 27 pairs),
    state counts 489 `PAIR_REVIEW_PENDING`, 30 `LOW_SIGNAL_UNRESOLVED`, 14 `MODEL_UNSTABLE`,
    5 `NO_TARGET_EXPRESSION_CONFIRMED`, 2 `INVALID_PAIR`; `method_version`
    `restore-pair-validation-v6`; no pair accepted and no divisor emitted
  - `expanded_review_v9/` was then the active review bundle — SUPERSEDED later the same day by
    `expanded_review_v10/`, see the final entry of this section: 13 shortlisted
    pairs, 69 checklist rows, 69 QuPath overlays, 1,649,671 crop cells, zero reviewer edits.
    The two CD20-referenced shortlist entries became `CD68 <- CD3e` and `CD11b <- CD3e`
  - separately, in the lineage sub-split B/Plasma now gate at the LOWEST immune precedence so
    DC / Macrophage / Neutrophil / NK / T claim their cells first (a CD20+CD68+ cell is now
    Macrophage), and CD20 was removed from `COMPARTMENT_GATES["Immune"]` — the gate is an
    any-positive union, so a lone CD20 call could otherwise pull an abundant cell into Immune
    before a more frequent compartment claimed it. B cells still reach Immune via CD79a. That
    code is dormant until Gate 4 and changes no current artifact
  - the legacy `MARKER_PAIRS` / `FUNCTIONAL_PAIRS` web and its `restore.run_restore` /
    `run_restore_functional` drivers are DELETED rather than maintained known-wrong until
    Gate 4: ten of those pairs referenced CD20 and there is no mechanical substitute (CD3e is
    invalid for CD4/CD8/FOXP3, which are CD3e+; epithelium expresses PD-L1 and IDO1). The
    reusable machinery stays — `idx_select`, `neg_stat`, `fit_clusters`, the threshold LUT and
    apply helpers, the QC plots, and the standalone proliferation pass. The pipeline's `restore`
    stage is kept in `ORDER` but RAISES, so a full run stops there instead of silently typing
    cells from gated parquet no accepted method produced
  - COMPARATOR ONLY, emits no divisor and no call: `immune_joint_comparator_20donor.csv`
    (480 rows, 20 donors) fits one pooled NNMF per donor over CD3e/CD20/CD68/CD11b x
    Membrane+Cell on the union of immune-positive cells. Anchor recovery, pairwise -> joint
    k=4: CD20 worst case 0.217 -> 0.968 (every donor >= 0.968), CD11b median 0.920 -> 0.696.
    Pooling rescues exactly the sparse marker and costs the abundant ones — this does NOT
    license a switch to joint fitting; it independently reproduces the conclusion that retired
    CD20 as a reference. k=5 is worse than k=4 throughout. Decide at the pair review
  - all 132 package tests pass; the retired names are asserted absent so they cannot be
    resurrected quietly

- 2026-07-23 (later the same day) — CD20 EXCLUDED FROM RESTORE ENTIRELY. This entry supersedes the
  one above, which retired CD20 only as a reference:
  - maintainer decision after reviewing the CD20 PDFs in `expanded_review_v9/`: the CD20 rows show
    overcalling throughout, and both qptiff pixel panels are overexposed with occasional block-like
    artifacts. CD20 is to take no part in RESTORE and be added only after the broad cell types are
    worked out
  - measured corroboration from the bundle's own per-row display windows (99.8th percentile of
    positive crop pixels, recorded in `qupath_exports[].marker_zoom`): CD20's ceiling runs 320-698
    (median 503) and CD163's 60-425 (median 181), against E-cadherin ~6,814, CD31 ~7,552 and CD11b
    ~10,225. A [1st, 99.8th] window plus gamma 0.65 stretches a near-background channel until
    tile-to-tile exposure differences render as blocks — which is the `CD20 <- CD163` panel exactly.
    The windows are NOT specific to the CD20 rows (donor 6380 CD68: 3,221 in the CD20 row vs 3,441
    and 3,642 in the CD3e/CD68 rows), so this is channel weakness, not a display quirk of those rows
  - as a TARGET the case is independent of the reference case: Step 1 anchor down to 14 cells,
    anchor recovery to 0.217, and 35 of CD20's 140 screen rows sat outside `PAIR_REVIEW_PENDING`
  - `INVALID_REFERENCE_MARKERS` is REPLACED by `RESTORE_EXCLUDED_MARKERS`, enforced at import
    against BOTH roles (plus the self-pair guard and the SHORTLIST-subset-of-CANDIDATE check).
    `BASELINE_PAIRS` 8 -> 7 (`CD20 <- E_cadherin` gone, so the set no longer equals the eight in
    `BroadLineageRevised.md`), `IMMUNE_REFERENCE_SCREEN_PAIRS` 7 -> 4,
    `EPITHELIAL_REFERENCE_COMPARATOR_PAIRS` 12 -> 9, `CANDIDATE_PAIRS` 27 -> 20,
    `SHORTLIST_PAIRS` 13 -> 11
  - CD20 is still EXTRACTED and still carries its reviewer-validated `max(Otsu, Triangle)` input
    floor: `DEFERRED_EXTRACTION_MARKERS` feeds `QUPATH_MARKERS`, so the deferred second pass needs no
    QuPath re-export and `qupath_compartment_sample_20donor/` + `qupath_compartment_full_20donor/`
    remain schema-identical and reusable
  - `nnmf_all20_expanded_v8/` was then the active screen (SUPERSEDED later the same day by
    `nnmf_all20_expanded_v9`, the divisor-gated rerun -- see the final entry): 400 rows (20 eligible
    donors x 20 pairs), state counts 384 `PAIR_REVIEW_PENDING`, 10 `MODEL_UNSTABLE`,
    6 `LOW_SIGNAL_UNRESOLVED`. The two states that disappear (`NO_TARGET_EXPRESSION_CONFIRMED` x5,
    `INVALID_PAIR` x2) were all CD20 rows
  - VERIFIED: all 400 surviving rows are numerically identical to their v7 counterparts on every
    shared column, confirming pairs are evaluated independently and that dropping CD20 perturbed
    nothing else
  - `expanded_review_v10/` was then the active review bundle (SUPERSEDED by `expanded_review_v11`);
    `expanded_review_v9/` and `nnmf_all20_expanded_v7/` are superseded and carried zero reviewer edits
  - donor 6539's `CONFIRMED_ABSENT` CD20 expression review is now inert rather than an error: it is
    looked up by `(donor, target)` only for screened pairs, and is retained in the manifest for
    provenance. The CD20-only "maintainer-confirmed target absence" stress-case selector is removed
  - lineage is unchanged but now abstains: `pos("CD20")` is all-False with no CD20 threshold, and
    CD79a is unscreened, so the B/Plasma branch yields nothing and those cells stay `Immune-other`.
    That is an abstention, not a negative call. The inert CD20 term is kept so the deferred pass
    needs no edit there
  - CD20 stays in `IMMUNE_JOINT_MARKERS` ONLY as a co-factor of the pooled comparator's fit —
    removing it would change the other three markers and break comparability with the run already on
    record. Its 0.217 -> 0.968 joint rescue does not license reintroduction
  - INVESTIGATED and RESOLVED (2026-07-23), replacing the "same failure mode as CD20" flag: the
    blocky `CD3e_from_CD163` qptiff panel is DISPLAY-ONLY, and the calculation is sound.
    - the display window (`[1st, 99.8th pct] + gamma 0.65`) is computed in
      `restore_pair_review._display_marker_channel` and used only for the RGB panels and the manifest
      metadata; it never re-enters a computed metric. CD163's crop median is ~7 counts, so that window
      stretches its noise floor into the mottled/block-like look; a linear render removes it, and the
      per-cell spatial field shows no dead tiles -- i.e. a weak-channel stretch, not acquisition damage
    - the calculation is unaffected: across all 20 donors CD3e's candidate divisor via CD163 matches
      its divisor via the five strong references within +/-20% (median 1.08x, range 0.80-1.28), and on
      the worst donor (6450, CD163 display ceiling 78) the NNMF control arm is a real CD163+/CD3e-low
      macrophage population of 1,448 crop cells -- nothing like CD20's 14-cell anchor
    - CD163 IS the weakest valid reference (`reference_fold` median 4.4x vs 7.6x CD68 and 19-25x the
      keratins, but above the 2x SBR floor) and carries a mild +8% divisor bias, so it is weak-but-
      valid, not the CD20 failure mode
    - CD163 is reference-only (never a screen target), so it produces no threshold of its own; whether
      CD3e should lean on the weakest reference when CD68 and three keratin references exist is a Gate-2
      reviewer PREFERENCE, not a correctness fix. Nothing needs regenerating
  - CAVEAT for the reviewer, from the same investigation: CD163 is used ONLY as a reference
    (`CD3e <- CD163`) and is never a target, and the same is true of CD79a/CD206/Iba1/CD11c/MPO. Those
    markers carry NO RESTORE `_pos` threshold under the current pair web, so the lineage sub-split's
    `pos("CD163")` etc. are inert and Macrophage currently rests on CD68 alone. This is a deliberate
    abstention (like CD20/CD79a for B cells), not a bug -- but any reasoning that assumes a CD163 (or
    CD206/Iba1/CD11c/MPO) positivity call is not yet satisfied by the screen
  - 133 package tests pass (1 new: the excluded marker cannot be reintroduced in either role)

- 2026-07-23 — CD11b overcall FIXED: the positive call is now divisor-gated only (METHOD_VERSION v7):
  - reviewer flagged CD11b as drastically overcalled across nearly all donors/pairs (only 6593 via
    PanCK looked right). Confirmed: the screen placed 94% of cells in the CD11b NNMF component (vs 26%
    CD3e, 53% CD68) and called ~32% of cells positive (up to 95% via CD3e), because `target_supported`
    called every target-component cell left of Step 1 EVEN BELOW the Step 2 divisor
  - root cause is NOT the reference and NOT display: RESTORE's reciprocal divisor (max CD11b among the
    PanCK+ target-negative controls) sets the correct threshold and is consistent across references
    (6521: 14% via PanCK, 18% via EpCAM; PanCK and EpCAM agree). The overcall was the below-divisor
    component RETENTION, which is negligible for clean markers (CD68/CD3e/B3TUBB add <=0.3 pp) but
    doubled CD11b (6521 CD11b<-PanCK 14% -> 28%). CD11b simply has a high channel background, which is
    exactly the background the divisor is designed to measure and threshold above -- so CD11b stays a
    valid target; the method, not the marker, was leaking
  - fix (user-directed "option 1"): `normalize_and_call` positivity is strictly `target_redsea >
    divisor`; the `| target_component_supported` rescue is removed from the CALL. Below-divisor
    component cells remain tracked (`target_supported_below_divisor`, new `threshold_positive_n`) and
    are shown in the bundle as a distinct muted "Target component at/below Step 2 (retained, not
    called)" category (code 5) so any marker the retention would help stays visible for a per-marker
    decision -- they are simply not counted. Bundle review classes renamed: code 2 "Called target+",
    code 5 "Retained below Step 2"
  - `nnmf_all20_expanded_v9/` is the ACTIVE screen (v7, divisor-gated); `expanded_review_v11/` the
    ACTIVE bundle. Surviving-row divisors are unchanged from v8 (the divisor math did not change, only
    what is counted positive)
  - 133 package tests pass (call semantics, category split, importer classes all updated)

- 2026-07-24 — CD11b<-CD3e RESOLVED and one reference frozen per target (Gate-2 reference selection):
  - the screen never combines references (each pair emits its own divisor; extra references per immune
    target are comparators), so the CD11b<-CD3e overcall was an ISOLATION artifact: CD3e is the
    minimum/outlier divisor of CD11b's references on 19/20 donors, and any single epithelial reference
    -- or a median/max across the five -- gives the QuPath-validated ~10%. No resampling can fix it:
    CD3e's control arm is <=1,706 cells on every donor (below the 2,000 equal-N cap) so its equal-N
    max equals its full max exactly, and no scheme can sample a CD11b value above the 9,923 that is
    the brightest CD3e control cell; the false-positive epithelial cells occupy CD11b 9,923-13,026, a
    background region only an epithelial reference contains
  - MAINTAINER DECISION: one predeclared reference per target; selection condition = the reference-
    positive negative control must be at least as large as the target-positive population, a PURE
    PROPORTION reference_n/target_n >= 1 with NO absolute cell-count floor (images vary widely in total
    and per-marker counts, so an absolute threshold is not comparable across donors)
  - METHOD_VERSION v7 -> v8: restore_validation.py adds reference_to_target_ratio and
    reference_undersized as reported-only diagnostics (config min_reference_to_target_ratio, default
    1.0). They never gate PairState and never alter a divisor. The 400-row screen was regenerated
    (nnmf_all20_expanded_v10) and VERIFIED byte-identical to v9 on every divisor/metric -- only
    method_version and the two new columns differ
  - new module phenocycler/reference_selection.py builds the Gate-2 dossier (reference_selection_v1/):
    reference_matrix.csv per candidate reference, frozen_reference_matrix.csv (accepted reference per
    target), 7 per-target figures, SUMMARY.md, manifest.json. It reads the screen and emits NO divisor,
    normalized value, or phenotype call
  - FROZEN reference matrix: CD31<-E_cadherin, Vimentin<-E_cadherin, B3TUBB<-EpCAM, CD3e<-E_cadherin,
    CD68<-E_cadherin, CD11b<-EpCAM, E_cadherin<-CD31. CD11b took a maintainer override to EpCAM (median
    ratio 3.66x, undersized 15%) over the marginal E_cadherin baseline (1.50x, undersized 30%);
    E_cadherin<-CD31 is accepted DESPITE failing the proportion (0.19x) because endothelium is the only
    biologically exclusive epithelial reference in the panel -- the frequency rule is advisory for a
    structural target with a large absolute control
  - REJECTED by the rule: CD11b<-CD3e (0.04x, undersized 100%), CD68<-CD3e (0.20x, 100%),
    CD11b<-Pan_Cytokeratin (0.78x, 60%). The divisor-consistency cross-check independently agrees --
    the undersized CD3e reference is also the low-divisor outlier (0.78x the epithelial consensus)
  - the v11 PDF review bundle was NOT rebuilt: divisors are byte-identical, so every PDF/overlay is
    unchanged and the new columns are not shown in them. reference selection is recorded in the dossier
    (a rule-based decision), leaving the maintainer's per-image review_decisions.csv pristine
  - 142 package tests pass (7 new: pure-proportion diagnostic, no-absolute-floor, config validation,
    reference-matrix rejection + divisor-outlier cross-check, baseline-vs-open recommendation,
    maintainer overrides incl. a failing reference, frozen-overrides constant)
  - production wiring (BASELINE_PAIRS / config.py / lineage.py) is UNCHANGED -- that is Gate 4; no
    RESTORE-normalized, phenotype, or app output was produced

- 2026-07-24 (later) — Gate 3 full-cell stress pilot RUN and PASSED; two high-background gate defects
  fixed (METHOD v8 -> v9):
  - isolated pilot `data/restore_pair_validation/gate3_pilot_v2/` (noncanonical): the 7 ACCEPTED_PAIRS
    across 7 donors (6539, 6476 + the weakest-case donor for each accepted pair: 6442, 6505, 6551,
    6566, 6591), marked ACCEPTED via `pair_reviews_accepted.csv` so computationally-clean donor-pairs
    reach VALID and emit a divisor + calls, while weak donors stay in their technical state and emit
    no call and no divisor
  - the pilot confirmed the divisor policy on full cells: every full-cell divisor is HIGHER than the
    5% screen's (the negative-control max is taken over the complete control population, not a 5%
    sample; median ~1.15x, up to 1.65x), so the emitted call fraction is correspondingly lower. The
    5% screen mildly over-called by under-estimating the divisor. CD11b <- EpCAM lands at ~7-10% on
    the clean donors (the QuPath-validated level), confirming the EpCAM override
  - the FIRST pilot run (`gate3_pilot_v1`, v8) produced two abstentions, both on HIGH-BACKGROUND
    markers (maintainer: 6476 has high CD3e background, 6591 high CD11b background -- NOT extreme true
    expression), where the provisional QC gates misfired even though the divisor was well-defined:
    - 6476 CD3e MODEL_UNSTABLE: NNMF components 99.8%+ identical across seeds, but near-separator
      control MEMBERSHIP dipped to 0.8846 on one seed, tripping `min_control_jaccard >= 0.90`
    - 6591 CD11b LOW_SIGNAL: blocked purely by `target_arm_fold 1.29 < 2.0`, though 51% of the target
      arm clears the divisor (~8% of assigned cells, the same the clean donors call)
  - fix (METHOD v9), since RESTORE depends only on the divisor:
    - seed-stability is gated on `min(seed divisor)/max(seed divisor) >= min_divisor_reproducibility`
      (0.90), NOT control-set membership. `control_jaccard` stays a reported diagnostic. Verified on
      6476 CD3e: `divisor_reproducibility = 1.0000` while `control_jaccard = 0.8846` -- the divisor is
      rock-solid though membership shuffles
    - a low target-arm MEAN SBR is now a reported flag (`target_low_sbr`), not an abstention (the
      manuscript prescribes flagging weak/near-background calls, not zeroing them; the gate plan
      explicitly deferred this 2-fold check for calibration). Genuine absence (`target_n < min_arm_n`)
      and confirmed absence still abstain
  - re-run `gate3_pilot_v2` under v9: ALL 49 donor-pairs VALID; divisors BYTE-IDENTICAL to v1 (the fix
    touched only the gates); every invariant holds (0 VALID without a divisor, 0 non-VALID emitting a
    call). 6476 CD3e and 6591 CD11b now emit their calls (0.44% / 8.06%)
  - 145 package tests pass. Submodule `8c77d05`, parent `1d5b80e` (pushed). NOTE: the 20-donor 5%
    screen (`nnmf_all20_expanded_v10`, v8) was NOT regenerated -- its divisors and frequency ratios
    (what the frozen reference matrix depends on) are v9-identical; only diagnostic state columns would
    shift, and the pilot already demonstrates the v9 state behavior. Reference selection stays frozen
  - NEXT: Gate 4 -- production wiring (config.py / lineage.py) with tri-state
    positive|negative|unavailable calls; ordered residual Immune -> Epithelial -> Endothelial ->
    Neural -> Mesenchymal -> Other; remove the silent missing-marker fallback in lineage.py

- 2026-07-24 (later) — Gate 4 production wiring CODE implemented; PILOT RAN BUT REJECTED, output reverted:
  - new `phenocycler/restore_apply.py`: the production adapter. Runs the frozen `evaluate_locked_pair`
    over the 7 `ACCEPTED_PAIRS` on all cells (modulus-1, unsampled) per donor and writes ONE tri-state gated
    parquet per donor (`{target}_pos` / `_norm` / `_state`, state in positive|negative|unavailable) plus a
    per-donor + cohort QC-disposition manifest. Pure, unit-tested mapper `cell_calls_from_evaluation`. It
    reuses the frozen math (no re-implementation): VALID -> divisor-gated positive/negative + off-`valid_idx`
    unavailable; NO_TARGET_EXPRESSION_CONFIRMED -> confident negative for every cell (maintainer E-1
    decision); any technical non-VALID on an accepted pair -> HALT (no automatic fallback). Atomic writes.
  - `lineage.py` rewritten to a TRI-STATE ordered residual: per gate `positive > unavailable > negative`;
    first positive gate claims; the first unavailable gate reached BLOCKS -> `Unresolved` (with
    `blocked_by`); valid-and-negative for every gate -> `Other`. The silent absent->negative fallback
    (`_pos_norm_getters`) is REMOVED (fail-loud `_state_getters` on a missing `{m}_state`). Subtypes
    (`_cell_types`: endocrine/exocrine, immune/mesenchymal fine types, NK) DELETED -> deferred to level 2;
    `cell_type == compartment`. Provenance columns: `assignment_reason`, `blocked_by`, `{m}_state`
    passthrough. Atomic write.
  - `config.py`: `COMPARTMENT_ORDER` now Immune-first; `COMPARTMENT_GATES` reduced to EXACTLY the 7
    divisor-backed accepted targets (Immune = CD3e/CD68/CD11b; Mesenchymal = Vimentin; SMA + the 6
    unscreened immune markers + CD20 deferred); `UNRESOLVED_LABEL` added; new derived paths
    (`restore_allcells_sample_dir`, `restore_pair_reviews_csv`, `restore_apply_manifest_csv`).
    `pipeline.py` `_restore` now calls `restore_apply.run_apply` (was a `SystemExit` placeholder).
  - Maintainer decisions baked in: 5-compartment level-1 deliverable only; confirmed absence -> negative;
    a technical non-VALID accepted pair halts (immune anchors are universally expressed). The deferred
    level-2 hierarchy is recorded verbatim in `Phenocycler_Analysis/docs/level2_hierarchy.md` (NOT wired).
  - 170 package tests pass (was 145): new `tests/test_restore_apply.py` (pure mapper per state + driver
    orchestration/fail-loud/end-to-end apply->lineage), rewritten `tests/test_lineage.py` (tri-state,
    Immune-first, Unresolved, fail-loud missing column, provenance), updated `tests/test_config.py`
    (new order/gates + gates == accepted-targets closure invariant). The 64 `test_restore_validation.py`
    tests are untouched and green.
  - PILOT DRY RUN (6539, 6476, 6591) RAN and was REJECTED by the maintainer. The apply reproduced the
    accepted method's calls to the cell (byte-identical to the Gate-3 pilot — e.g. 6539 E-cadherin 93,749
    / CD11b 48,956), so this is NOT a Gate-4 wiring bug. But the resulting compartment composition is
    biologically implausible: ~15-19% Epithelial and ~65% `Other` across all three donors, when acinar
    epithelium should dominate the pancreas. The immediate driver is that the accepted method calls
    E-cadherin on only ~18% of cells (only ~18% exceed the E_cadherin<-CD31 divisor); the ROOT CAUSE is
    upstream of the Gate-4 wiring and is NOT yet resolved. The pairs are pre-decided/fixed, so it is not a
    pair-selection issue. The pilot output (`data/phenotype/`, `data/restore_gated_redsea/`) was DELETED to
    return `data/` to its pre-pilot state. **Gate 4 is therefore code-complete but NOT validated.** NEXT:
    diagnose and fix the upstream E-cadherin under-call, then re-run the pilot. No canonical 20-donor run
    and no downstream/app rebuild until it passes (Gates 5-7).

## Ground truth

The manuscript and supplement at
`/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/{RESTORE,RESTOREsup}.pdf`
override prior code, handoffs, notebooks, and documentation.

The paper's operative model is:

- Work locally per tissue sample/image on linear single-cell mean intensities.
- For target marker `i`, use a biologically validated mutually exclusive marker `j`
  to identify cells that are truly negative for `i`.
- Marker pairs may be fixed from biological knowledge or discovered from the data.
  Identify the cell groups with a constrained NNMF, SSC, GMM, or another justified
  method. The paper used two-group SSC in its reported implementation but explicitly
  permits more than two groups and other methods.
- Infer target background `B_i` from the target intensity in the reference-positive,
  target-negative group.
- Normalize every target measurement as `y_i / B_i`. The paper's equations use the
  maximum negative-control value, making negative controls `<= 1` and signal `> 1`;
  the released code instead uses negative-cluster `mean + 3*SD`, so this discrepancy
  requires an explicit decision before implementation.
- Do not borrow thresholds across donors, normalize class frequencies toward a cohort
  pattern, or manufacture a threshold when a valid positive reference population is
  absent.
- Treat high-SBR measurements as quantitatively reliable; explicitly flag weak,
  near-background calls because the paper states that normalization distorts that region.
- Validate mutual exclusivity visually and biologically for each sample. More than one
  candidate reference may be inspected, but the paper does not define the current
  median-of-reference-thresholds procedure.

Current data have one image per donor, so per-image and per-donor normalization units
coincide for all 22 donors. This must be asserted in code, not assumed.

## Phase 0 - Quarantine stale outputs and protect valid inputs

1. Treat these as valid inputs:
   - raw `data/cells/`
   - raw-source `data/singlecell_protein.h5ad` only after provenance confirmation
   - REDSEA `data/cells_redsea/`
   - REDSEA scratch/intermediates and REDSEA-only validation artifacts
2. Build a deletion manifest with path, size, mtime, producer, and dependency class.
3. Delete confirmed RESTORE outputs and direct descendants:
   `restore_redsea/`, `restore_gated_redsea/`, RESTORE threshold/pickle files,
   `phenotype/`, QuPath phenotype exports, and RESTORE/phenotype figures.
4. Older app H5AD/parquet/islet artifacts are included in the purge. The deployed app
   intentionally has no usable data until the accepted rebuild.
5. Verify retained raw/REDSEA files after deletion: donor set, schemas, row counts,
   unique `object_id`, and one-to-one key coverage. The first sampled donor already
   shows differing raw and REDSEA row counts, so key reconciliation is a hard blocker.

## Readiness remediation plan

This plan supersedes the earlier phase list. No RESTORE-normalized, phenotype, or
downstream app output may be generated outside an isolated pilot until Gate 5 is
accepted.

### Gate 1 - Lock the executable method specification

1. Separate the three data roles:
   - raw QuPath compartment means for target/reference group inference
   - canonical REDSEA whole-cell target values for background estimation and calls
   - raw whole-cell values only for matching, diagnostics, and raw-to-REDSEA comparison
2. Assert exact `object_id` lineage and one-image-per-donor before any model runs.
3. Define input QC as a fitting operation only:
   - remove the double-low population from separator fitting
   - robustly min-max scale each compartment feature for fitting only
   - train on equal-size provisional target-only and reference-only arms
   - assign every canonical cell; never discard low cells from the final dataset
4. Keep actual cell frequencies untouched and report them separately from balanced
   model/plot samples.
5. Resolve the threshold-estimation policy explicitly:
   - the manuscript full negative-control maximum remains the default canonical divisor
   - compare it with deterministic fixed-N sampled maxima over a predeclared sample-size
     grid and repeated seeds to quantify order-statistic inflation
   - use sampled intervals as QC, not as the production divisor, unless the user
     explicitly approves a documented pancreas-panel modification
   - keep `mean + 3*SD` as labeled released-code sensitivity only
6. Replace the provisional binary PASS/FAIL model with an auditable state machine:
   - `VALID`
   - `NO_TARGET_EXPRESSION_CONFIRMED`
   - `LOW_SIGNAL_UNRESOLVED`
   - `NO_REFERENCE_POPULATION`
   - `INVALID_PAIR`
   - `MODEL_UNSTABLE`
   - `TECHNICAL_FAILURE`
7. Calibrate the current provisional 2-fold target-arm support check against
   donor-local images reviewed by the user. It is not production-ready merely because
   it correctly flags donor 6539 CD20.
8. Add deterministic tests for every state, scaling boundary, balanced sample,
   canonical divisor, strict `>1` call, and prohibited fallback.

Exit criterion: one versioned method specification and test suite that states exactly
which values are used for fitting, divisor estimation, normalization, and calling.

### Gate 2 - Redesign and validate marker pairs

1. Run a two-pass screen:
   - deterministic sample of all 20 eligible donors to rank biologically plausible candidates
   - full-cell validation for every finalist
2. For CD3e, CD20, CD68, and CD11b, screen alternative counterparts from this panel
   using biological incompatibility first, then data-driven mutual-exclusivity ranking.
   Do not select a convenient counterpart independently for each donor.
3. Revalidate E-cadherin/CD31 reciprocals, Vimentin/E-cadherin, and B3TUBB/EpCAM with
   the locked method. Specifically test:
   - CD31-positive leukocyte contamination
   - genuine epithelial Vimentin expression
   - B3TUBB-positive endocrine cells in the EpCAM reference population
   - sparse-arm stability without assuming the sparse population is erroneous
4. For every donor/pair, report reference presence, target-expression support,
   double-high overlap, arm sizes, repeated-seed/sample stability, compartment
   completeness, raw-to-REDSEA changes, spatial coherence, and threshold uncertainty.
5. Review all `NO_TARGET_EXPRESSION` candidates and all failures in QuPath:
   - confirmed biological absence produces no positive calls for that donor-marker
   - weak staining or ambiguous absence remains technical/unresolved
6. Predeclare one fixed counterpart, or one fixed multi-reference rule, per target
   before production. If no rule works across the cohort, remove or replace that
   lineage anchor rather than using an automatic per-donor fallback.

Exit criterion: a user-approved donor-by-pair matrix with a disposition for every
required target and a frozen pair/marker set.

### Gate 3 - Run a full-cell RESTORE stress pilot

1. Use full-cell data, not the 5% evaluation sample, for donors 6539 and 6476,
   plus any eligible donor identified as the weakest case for an accepted pair. Donors
   6457 and 6579 remain excluded.
2. Infer controls from raw compartment features, calculate candidate divisors on the
   matched REDSEA target values, and call signal strictly above the accepted divisor.
3. Produce an isolated pilot bundle:
   - threshold and normalized-value tables
   - pair clouds and uncertainty diagnostics
   - QuPath import files and spatial overlays
   - complete provenance and state-machine dispositions
4. Confirm that `NO_TARGET_EXPRESSION_CONFIRMED` yields no positive calls without
   fabricating a divisor, while technical failure remains unresolved.
5. Stop and return to Gate 2 if any required anchor lacks a defensible control group.

Exit criterion: user approval of every required pilot gate and the sampling/divisor
policy. Pilot output is still noncanonical.

### Gate 4 - Replace obsolete production wiring

1. Update `config.py` to contain only the approved broad-lineage markers and validated
   target/reference rules.
2. Change the residual priority to:
   `Immune -> Epithelial -> Endothelial -> Neural -> Mesenchymal -> Other`.
3. Remove the current silent missing-marker behavior in `lineage.py` that converts an
   absent column to all-false/zero.
4. Implement tri-state marker calls (`positive`, `negative`, `unavailable`) plus the
   donor-marker QC disposition and provenance columns.
5. Implement ordered residual assignment so:
   - a valid earlier positive gate may claim a cell
   - an unavailable earlier gate blocks a later claim when it could change the class
   - `Other` is assigned only when every required gate is valid and negative
   - blocked cells are `Unresolved`, never biological `Other`
6. Stop the first accepted deliverable at five broad compartments plus `Other`;
   remove/defer endocrine/exocrine and other subtype logic from this run. The deferred
   level-2 subtype design (endocrine via CD56-in-E-cadherin+, T/B/myeloid/NK subsets,
   endothelial/mesenchymal subtypes, state markers, and their counterpart rules) is
   recorded verbatim in `Phenocycler_Analysis/docs/level2_hierarchy.md` — NOT wired; each
   level-2 pair needs its own Gate-1->3 validation before production.
7. Write outputs atomically and fail loudly on missing donors, columns, manifests,
   invalid pairs, or key mismatches.
8. Add unit and integration tests for ordering, multi-positive cells, confirmed
   marker absence, technical failure, `Other`, and `Unresolved`.

Exit criterion: a clean dry run on synthetic fixtures and the accepted pilot donors,
with no silent defaults and no downstream rebuild.

### Gate 5 - Validate all 20 eligible donors with immutable parameters

1. Run the production candidate once across all 20 eligible donors with frozen parameters,
   seeds, pairs, sample sizes, and software versions.
2. Generate a compact donor-by-marker dashboard and per-donor drill-down bundle.
3. Validate each donor against its own distributions, negative controls, image
   localization, and spatial coherence. Do not normalize toward cohort frequencies.
   Donor 6476 is not a pancreatitis donor; donor 6579 is the pancreatitis donor and
   remains excluded by maintainer decision.
4. Resolve every donor-marker state. No required marker may remain accidentally
   missing or silently negative.
5. Require user sign-off on the 20-donor dashboard and representative QuPath overlays.

Exit criterion: immutable all-donor threshold/call manifests and explicit user
acceptance. Only then may broad-lineage assignment become canonical.

### Gate 6 - Assign and validate broad lineage

1. Apply the accepted tri-state residual tree to all canonical cells.
2. Audit marker logic, multi-positive tie resolution, `Other`, `Unresolved`, and key
   preservation for every donor.
3. Review spatial overlays for each broad class, emphasizing immune/endothelial
   boundaries, endocrine B3TUBB exclusion, and Vimentin-positive residual cells.
4. Report donor variation descriptively without using class-frequency conformity as a
   QC target.

Exit criterion: accepted broad-lineage parquets plus QuPath review artifacts, with all
31,454,041 eligible canonical cells accounted for exactly once.

### Gate 7 - Rebuild descendants

1. Rebuild single-cell H5AD, islet aggregations, trajectory, app H5AD/parquet, and
   deployed app assets only from the accepted lineage version.
2. Update documentation to cite the manuscript and remove obsolete pair, threshold,
   ordering, subtype, and forced-assignment claims.
3. Run package tests, pipeline dry runs, key-integrity checks, app smoke tests, and
   provenance checks.

Exit criterion: downstream artifacts identify the accepted lineage version and pass
all package, pipeline, data-integrity, and app checks.

## Required user approval checkpoints

1. Final divisor policy after comparing full-control and fixed-N sampled maxima.
2. Final fixed marker/reference rules, especially the immune anchors and B3TUBB.
3. Donor-marker dispositions for low/no expression versus technical failure.
4. Full-cell pilot gates.
5. All-22-donor dashboard before canonical lineage assignment.

## Resolved user decisions

1. Purge all app/downstream artifacts as well as RESTORE and phenotype outputs.
2. Deliver five compartments plus `Other`; defer endocrine/exocrine.
3. Use the manuscript maximum negative-control divisor; retain `mean + 3*SD` only as
   sensitivity analysis.
4. Fail and redesign any invalid marker/image pair before continuing; no automatic
   fallback.
