# Sequential-section regional correspondence and diagnostics

## Correction: projected overlap is not required

Sequential sections sample different depths. Related anatomical regions can
change position, shape, extent and segmentation and can have zero projected
intersection after registration. Projected overlap is neither a prerequisite
nor sufficient evidence for relating their native regional features.

The earlier November 5 interpretation overextended the MALDI report's sampling
fidelity method into a criterion for serial-section correspondence. The 152-pair
high-overlap set and the 62-pair overlap-under-shifts set are geometric diagnostic
subsets, not accepted-pair counts. None of the 117 zero-overlap pairs is rejected
on overlap alone. All 473 historical candidates remain; the earlier 297-pair
stability/cell-support sensitivity set includes 24 zero-overlap pairs.

## Correspondence and feature comparison

The historical pairing uses independently measured regional objects and the
stored Xenium-to-PhenoCycler transform, with position, size, shape and competing
candidate costs. The audit also assesses the whole-section structural pattern.
All eight transforms were structurally corroborated, including six historically
marked WARN. This evidence supports qualified exploratory regional comparison;
it does not individually validate every pair.

The existing sensitivity set requires stable pair IDs across all three original
matching parameter settings, no low cost-margin flag, and at least 10 cells in
both assays. RNA cell support is reassessed at each original depth setting.
There is no projected-overlap cutoff, overlap-derived confidence or overlap
weight in the primary or correspondence-sensitivity comparison. Native features
are compared with equal pair weights and a common finite-measurement denominator;
measured zeros remain values and missing values remain unavailable. Donor remains
the independent unit for biological inference. Cell identities are not transferred.

Anatomical context, relative structural arrangement, independent registration
landmarks/residuals, competing candidates and pair-ID stability are the relevant
correspondence evidence. Local anatomical agreement has not been independently
reviewed for every historical pair. Relating tissue areas must also be kept
separate from asserting an identical islet or identical cells.

## What the MALDI methods measure

The source is `MALDI/T1D_Nglycan_CRO_Report.pdf`, 9 July 2026, sections 9.4–9.5,
pp. 19–20. Fractional spot coverage and precision/recall/Dice/IoU describe how a
sampling footprint represents an annotated compartment. They remain useful for
MALDI mixture and spatial-resolution diagnostics. They do not establish that
serial-section regions must intersect in projection to be related. The report's
one-core masking scores are not serial-section acceptance thresholds.

`multiplex5000.integration.serial_features` provides:

- `overlap_metrics`: exact projected intersection, union, directional coverage,
  Dice and IoU, preserving holes. Zero intersection is an observed geometric
  difference, not a failed anatomical relationship.
- `fractional_links`: an intersection graph, with uncovered and competing area.
  It is not a complete correspondence graph. Missing edges do not rule out
  relationships. Fractions are not renormalized; missing geometry remains
  unavailable and invalid geometry is rejected without repair.
- `footprint_coverages` and `coverage_fidelity`: native sampling/mixture and
  representation diagnostics. Footprints must not overlap when aggregating area.
- `paired_feature_summary`: equal weights by default. Optional weights require
  an explicitly labeled diagnostic; projected coverage must not silently filter
  serial-section comparisons or be called a correspondence probability.
- `shift_sensitivity`: change in projected overlap for fixed pair IDs under
  imposed translations. This does not measure stability of the pair identities.

The reported 30/50/70/90% coverage curves, 25/50/100 µm overlap perturbations,
500 µm displacement controls and virtual 50 µm masks are retained as diagnostics.
Their values do not determine eligibility, release status or biological agreement.
The report's classifier, normalization and expected biological findings are not
used to choose correspondences or thresholds.

## Reproduction and current assessment

The subsequent anatomical-neighborhood reassessment is documented in
`SEQUENTIAL_ANATOMICAL_CONTEXT.md` and implemented by
`multiplex5000.integration.serial_context`. It tests all six WARN registrations
using neighboring structures and morphology-weight sensitivities, with no
projected-overlap requirement.

The corrected consumer uses the verified historical pair measurements and native
own-section counts, preserving all existing transforms and geometry diagnostics:

```sh
.venv/bin/python /home/smith6jt/Documents/papers/consolidation/sequential_correspondence_20261105/correct_correspondence.py --out /absolute/new/output
```

The current report and machine-readable assessment are in
`/home/smith6jt/Documents/papers/consolidation/sequential_correspondence_20261105/`.
The earlier `sequential_methods_20261105` package is a historical diagnostic
snapshot; its overlap-based eligibility interpretation is superseded.

The consumer verifies that selection is invariant to removing or changing all
overlap measurements, retains the 24 zero-overlap pairs in the 297-pair set,
and reproduces the previously computed unweighted native feature comparisons.
The module tests cover exact geometry, conservation, missingness and the default
equal-weight comparison. Historical PASS/WARN fields, typing limitations and
source-provenance caveats remain explicit; no production receipt is changed.
