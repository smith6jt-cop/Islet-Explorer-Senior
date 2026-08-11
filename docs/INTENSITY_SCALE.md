# The intensity scale: what zero and maximum mean

Status: **measured, specified, not yet implemented.** Evidence is in
`data/pilots/acquisition_by_biomarker.parquet` and the pixel measurements below.

Every threshold in this workflow — the Otsu input floor, the NNMF arms, the
RESTORE divisor — is a point on an intensity scale whose two ends have never
been defined. This document defines them, and states precisely what may and may
not be standardized.

## The governing rule

**Standardize the unit, never the distribution.**

An anchor is legitimate when it is a property of the *instrument* or of the
*measurement's own censoring*, measured per donor and per marker from that
donor's own image. An anchor is illegitimate when it is a property of the
*cell population*, because that is biology and rescaling by it erases exactly
the differences the study is trying to measure.

This forbids, permanently:

- quantile normalization or z-scoring across donors;
- matching a marker's distribution between donors;
- any anchor taken as a quantile of the whole cell population — a quantile
  depends on how many cells are positive, so dividing by it would flatten real
  donor differences in abundance;
- any anchor chosen because it makes cell-type frequencies look expected.

## What was measured (2026-08-11)

### The ceiling is intact but the unit is not

Level-0 tile sampling, 63 million tissue pixels across donors 6436, 6450 and
6566 × 7 markers: **zero pixels at the 16-bit ceiling of 65535**, and single
digits at ≥64000. The brightest pixel observed anywhere was 65324. Saturation is
not currently corrupting this cohort.

The scale itself, however, is donor-specific. From the qptiff acquisition
metadata across all 22 slides:

- **Exposure time differs between donors for 16 of the 22 broad-gate markers**,
  by up to 1.67× (CD11c 60,000 vs 100,000; CD31 1.56×; CD206, GCG, CD34 1.50×).
- **Filter responsivity was recalibrated three times** over the cohort
  (2025-09-10, 2026-03-12, 2026-06-11), differing by up to 1.41× on Cy5.

So the same antigen density produces materially different counts on different
donors. The three calibration epochs also track the acquisition batches, which
[[gate5-batch-confound]] already shows are aliased with disease status — this is
a candidate mechanism for that confound and should be checked against it.

### Zero is a clamp, applied twice

The scanner reports `ImageProcessingApplied=OffsetAndShading` and
`AutofluorescenceSubtracted=True`. Both are subtractions, and the result is
stored as unsigned 16-bit, so every negative excursion is clamped to 0.

| | fraction of tissue pixels exactly 0 |
|---|---|
| DAPI | 0.93% (6566) – 10.50% (6450) |
| INS | 56% – 79% |
| E-cadherin | 30% – 50% |
| CD31 | 25% – 53% |

REDSEA then clamps a second time (`np.clip(corrected_sum, 0, None)`,
`redsea.py:529`), which is what produces exact zeros at the *cell* level:

| | fraction of analysis-eligible cells exactly 0 after REDSEA |
|---|---|
| CD31 | 6.6% (6566) – 19.0% (6436) |
| Vimentin | 6.1% (6450) – 14.1% (6566) |
| E-cadherin | 9.0% – 13.4% |

A value of exactly 0 is therefore **not a measurement of no signal**. It means
"≤ 0 after subtraction", and the amount of the distribution hidden under that
clamp varies by up to 11× between donors for the same marker.

## What this does and does not affect

Stated precisely, because the scope is narrower than it first appears.

**Not affected — RESTORE's positivity calls.** The divisor is estimated from the
same donor and marker as the values it divides, so a per-donor-per-marker scalar
cancels exactly. Exposure normalization will not move a single call. Anyone
expecting it to fix the divisor problem will be disappointed.

**Not affected — REDSEA spillover correction.** `F` is a cell–cell contact
matrix and `F @ edge` acts on each channel independently
(`redsea.py:11`, `:373`). Channels with different exposures are never mixed.

**Affected — the input floor, and therefore the calls.** `raw_otsu_floor` is
computed on a distribution carrying a donor-varying atom at 0. Otsu maximises
between-class variance, so a large atom at zero pulls the threshold. The floor
defines arm membership, which defines the NNMF separator and the target anchor —
the quantity every divisor rule is judged against. This is the path by which
censoring at the bottom becomes an error in the calls.

**Affected — every cross-donor intensity comparison.** The Shiny app's Card C
and the single-cell H5AD both carry REDSEA MFI directly, in units that differ
between donors by up to 1.67× for instrument reasons alone.

## The specification

Two anchors per donor per marker, both from that donor's own image.

### `z` — the zero anchor

Zero is a censoring boundary, not a value.

1. Estimate the background from the **strictly positive** part of the
   distribution only; the clamped mass at 0 must never be treated as observed
   background, because its size is an artefact of how far the offset and
   autofluorescence subtractions overshot.
2. Record `censored_below_fraction` per donor per marker at both pixel and cell
   level as a first-class audit quantity, not a diagnostic afterthought.
3. A cell whose value is exactly 0 is `censored_below`. It must not be reported
   as a confident negative: its true corrected value is unknown and ≤ 0.
   Under the tri-state contract that is `unavailable` with an explicit reason,
   never `negative`.

### `s` — the scale anchor

Divide counts by `exposure_time × responsivity`, both read from the qptiff's own
acquisition metadata for that donor and channel. This is an instrument
correction: it is a scalar per donor-marker, it is independent of every cell in
the image, and it restores a common physical unit without touching distribution
shape.

Join on the cohort manifest's `channel_map`, never on the qptiff `Biomarker`
string. The scanner writes `E-cadherin`, `Ker8-18`, `Keratin 5` and
`Pan-Cytokeratin`, while the registry uses underscores; `channel_map` already
carries `channel_index → channel_name → marker` per donor and is the reconciled
authority. Matching on the raw string silently drops four gate markers.

Add a ceiling guard even though nothing is saturated today: record
`saturated_pixel_fraction` per donor-marker, and treat a cell containing
saturated pixels as `censored_above` — safe to call positive, but its value is a
lower bound and it must not contribute to any statistic that summarises the top
of the distribution, including the manuscript maximum divisor. QuPath already
exports `Cell: <marker>: Max` for all 59 markers, so this is measurable per cell
from the existing measurement CSV without touching the images.

### The invariant that proves it is honest

The transform is affine and monotone within each donor-marker. Therefore
within-donor ranks, and every RESTORE call, are unchanged by it. If applying the
standard changes a donor's cell-type composition, the implementation is wrong —
it has stopped standardizing units and started reshaping distributions.

## Also found

`6523panc_Scan1.er.qptiff` carries the internal `SlideID` **`6623panc`**.
**Confirmed by the maintainer (2026-08-11): `6623` was typed in error when the
scan for donor 6523 was set up.** The image is donor 6523's tissue and the
filename is correct. The metadata is consistent with that — the two files are
distinct acquisitions with different `Identifier`, dimensions (30240×24000 vs
48960×28800) and channel counts (58 vs 59), so nothing was duplicated or
swapped.

Recorded because it is a standing reason never to key a donor on embedded image
metadata: the cohort manifest is the only donor authority, and any future code
that "corrects" the manifest from the qptiff's `SlideID` would silently reassign
a T1D donor's cells.

`6523panc`, `6566panc` and `6591panc` have 58 channels rather than 59, all
missing `b-Catenin1`. It is not a broad-gate marker, so broad typing is
unaffected.
