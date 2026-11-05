# Sequential-section feature support

The implementation is `multiplex5000.integration.serial_features`. It measures
the support for comparing features between independently sampled area objects.
It does not transfer cell identities, fit a registration, normalize expression,
or change a release acceptance receipt.

The motivating source is `MALDI/T1D_Nglycan_CRO_Report.pdf`, report date July 9,
2026, sections 9.4–9.5 on pages 19–20, with the masking and purity figures on
pages 8–10 and 29. The source's one-core masking scores are not thresholds or
expected results for another acquisition.

## What is applied

- **Area support:** `overlap_metrics` reports exact intersection/union area,
  source-only/target-only area, both directional coverage fractions, Dice and
  IoU. Native vector objects preserve holes and avoid introducing a raster
  approximation into the Xenium–PhenoCycler comparison. Neither section is
  designated biological ground truth. An observed disjoint pair has zero
  overlap; unavailable geometry has missing metrics.
- **Fractional assignment:** `fractional_links` records every positive-area
  object link, unassigned source area and area covered by competing targets.
  Competing area is excluded from each link's exclusive fraction. Fractions are
  not renormalized to make incomplete support sum to one. Splits and merges
  remain explicit. Inputs require the same declared calibrated `:um` frame and
  stable native object IDs; invalid area objects are refused without repair.
  If some target geometry is missing, coverage describes only known targets;
  exclusive, ambiguous and unassigned fractions remain unavailable.
- **Sampling fidelity:** `footprint_coverages` intersects sampling footprints
  with the union of reference objects, so overlapping objects do not count
  twice. `coverage_fidelity` implements the report's area precision, recall,
  Dice and IoU at coverage thresholds 0.3, 0.5 and 0.7 (other thresholds can be
  supplied). These measure representation loss, not registration accuracy.
  Sampling footprints must be non-overlapping for aggregate area metrics.
- **Feature application:** `paired_feature_summary` computes descriptive
  comparisons of the two native regional measurements on their common finite
  support. A geometric fraction may weight a comparison; it is not a calibrated
  probability. Measured zeros remain values. Missing measurements and zero
  support contribute no weight. The effective weight count is not an independent
  biological sample size. No features are imputed or copied into another assay.
- **Transform sensitivity:** `shift_sensitivity` keeps pair IDs fixed while
  translating the source by 25, 50 and 100 µm in eight directions. Reported
  ranges include the stored transform. They are sensitivity ranges, not
  uncertainty intervals, probabilities or a reason to choose a better-fitting
  transform. Calibrated uncertainty requires independently measured residuals.

## Applying it to the October 6 run

The reproducible consumer is:

`/home/smith6jt/Documents/papers/consolidation/sequential_methods_20261105/run_feature_validation.py`

Run it using the repository environment, with a new output directory:

```sh
.venv/bin/python /home/smith6jt/Documents/papers/consolidation/sequential_methods_20261105/run_feature_validation.py --out /absolute/new/output
```

The runner pins the already-audited input files, reconstructs the exact eight
donors' geometry in the PhenoCycler micrometre frame and measures all three
historical matching settings. It retains the original 473 primary candidates
and their original PASS/WARN fields. It adds:

- Pair-specific area and competing-object support, and all fractional links.
- Purity–yield curves at 0.3, 0.5, 0.7 and 0.9 for all primary candidates and
  the earlier 297-pair geometry/cell-support sensitivity set.
- A separately labeled sensitivity set requiring that earlier flag plus at
  least 0.5 exclusive reciprocal area support. This threshold is an explicit
  sensitivity choice, not a fitted accuracy threshold or a new primary result.
- Native regional feature comparisons at all four original RNA-depth settings,
  including continuous support weighting and thresholded sensitivities. Cell
  support is reassessed at each depth.
- Structured 500-µm displacement controls in 16 directions, preserving source
  shape and spacing. Because displacement can change tissue support, these
  controls produce descriptive overlap comparisons, not null p-values.
- Virtual 50-µm sampling at two grid phases, explicitly separated from native
  geometry and from cross-section registration accuracy. This is the direct
  analogue of the report's coarse-mask fidelity calculation, not a decision to
  downsample Xenium or PhenoCycler to MALDI resolution.

Low area overlap does not by itself prove a wrong anatomical pair: separate
sections can change shape, and segmentation can differ. It does expose weak
support for assigning/comparing a feature as though its entire region were
shared. Whole-section registration corroboration cannot substitute for this
pair-level measurement.

## Boundaries of the report adaptation

The report estimates scale from annotation point-cloud extent and fits tissue
footprint overlap. That does not replace instrument pixel calibration or
independent registration review here. Existing transforms stay fixed. Its
glycan-trained classifier, APPNP smoothing, disease-protected normalization,
peri-islet distance bins and disease results are not part of this implementation.
The report itself shows that random-split classification gains weaken under
spatial holdout and that the classifier overcalls islet regions. A future learned
assignment method must therefore evaluate held-out spatial blocks and donors;
there is no learned classifier to cross-validate in this geometry-only module.

No thresholds are selected to recover a biological result. Donor remains the
independent unit for disease inference. An area-support pass is not a production
gate acceptance, cell-level validation, or certification of the historical
typing labels.

## Validation

`tests/current/test_serial_features.py` uses analytically known geometry to
check fractional overlap, holes, three-way ambiguity, conservation of exclusive
plus ambiguous plus unassigned area, missing versus empty source sets, frame
errors, measured zeros, missing measurements and transform-shift sensitivity.
The real-data runner checks that exact virtual footprint coverage sums to the
native object-union area. The normal repository `make test-fast` includes the
module tests.
