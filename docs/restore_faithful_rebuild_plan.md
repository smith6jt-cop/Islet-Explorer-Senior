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
> normalized, phenotype, or app output exists. The immediate next action is recording the
> human pair review in `data/restore_pair_validation/expanded_review_v9/review_decisions.csv`
> (69 rows, against the `nnmf_all20_expanded_v7/` screen).
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
    30 fixed candidate pairs — SUPERSEDED 2026-07-23 by `nnmf_all20_expanded_v7/` (540 rows,
    27 pairs) when CD20 was retired as a reference; see the final entry of this section
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
    appears — SUPERSEDED 2026-07-23 by `expanded_review_v9/` (69 rows); v7 and v8 both
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
  - `nnmf_all20_expanded_v7/` is the ACTIVE screen: 540 rows (20 eligible donors x 27 pairs),
    state counts 489 `PAIR_REVIEW_PENDING`, 30 `LOW_SIGNAL_UNRESOLVED`, 14 `MODEL_UNSTABLE`,
    5 `NO_TARGET_EXPRESSION_CONFIRMED`, 2 `INVALID_PAIR`; `method_version`
    `restore-pair-validation-v6`; no pair accepted and no divisor emitted
  - `expanded_review_v9/` is the ACTIVE review bundle, built from that screen: 13 shortlisted
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
   remove/defer endocrine/exocrine and other subtype logic from this run.
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
