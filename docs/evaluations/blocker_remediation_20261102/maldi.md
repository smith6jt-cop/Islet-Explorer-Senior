# MALDI registration and object support remediation

The Xenium source reader and explicit object integration contracts are implemented.
All 65 inventory acquisitions were checked: 62 exact-region morphology sources
validated and the three recorded damaged regions remained unavailable. All 62
DAPI/mask overlays were inspected. This establishes source availability, not
physical adjacency, mask acceptance, or a valid MALDI transform.

The real MALDI cohort still contains 32 sections and 1,666,450 spots. Its original
source, identifier and denominator audits remain unchanged. Each actual
registration requires an inspected exact neighbour acquisition and receipt-bound
acceptance; those records have not been supplied or inferred.

## Changes and verified behavior

The registration reader binds the exact inventory acquisition, experiment hash,
declared DAPI channel, OME axes, physical pixel size and image content. It uses the
native region coordinate frame and a 50 µm grid with explicit pixel-centre origin.
An explicit source index may relocate a bundle to a differently named directory;
the acquisition/output mapping and experiment hash must still match. An undeclared
rename, wrong region, changed experiment, escaped image path, wrong channel or
incorrect calibration fails validation. Missing or empty source evidence produces
an explicit per-section unavailable reason without silently selecting another source.

The new `maldi_object_sources` contract supports independently identified area
objects in exact PhenoCycler or Xenium acquisition frames. Source files, calibration,
native object identifiers and producer or observation receipts are bound explicitly.
An observation receipt documents inspection of an existing export; it does not
claim to recover its historical producer. The prepared pLN source set contains
19 acquisitions and 47 native annotation objects. Only 44 are area geometries;
the three LineStrings from 6534 remain unavailable for area summaries. Bounding
rectangles remain annotation regions, without being relabelled as tissue outlines
or follicles.

Stage 08 checks the accepted Stage 04b/07 receipt and gate chain, every referenced
spot/feature/transform/displacement output, exact source and target frames, and
object-source bindings. It rejects duplicate identifiers and changed accepted
artifacts. Overlapping objects produce ambiguous spot membership rather than
first-object ownership. Per-object geometric summaries are explicitly nonexclusive.
An unregistered section remains in `registration_unavailable.csv` and section counts.

Measured zeros and unavailable values remain distinct. Summaries expose finite
denominators and detection fractions for the existing signal mask, all raster spots,
and uncertainty-eroded cores, with explicit detected counts for the signal-mask
and all-raster summaries. The existing outline RMSE + half-pitch erosion is a
geometric sensitivity/eligibility margin, not a statistically calibrated confidence
bound. The existing RMSE + 25 µm erosion is retained;
all-inside summaries are not silently relabelled as core summaries. Signed islet
distance remains continuous. Invalid or absent islet geometry makes that distance
unavailable. There is no cell-identity transfer, same-cell MALDI claim, or new
distance band.

## Source and image inspection

The local source check validated 30 regions. A corrected indexed remote check
validated 62; all 30 shared image files matched exactly by full SHA-256. The first
remote attempt lacked the relocation index. The next attempt rejected two
explicitly relocated directory names. Those attempts and logs remain separate;
the final run used the tested identity-preserving relocation rule.

Evidence is under
[`hpg_maldi_xenium/results_relocated`](../../../outputs/blocker_remediation_20261102/hpg_maldi_xenium/results_relocated/).
The report binds the validator, registration source, inventory, source index,
environment, image hashes, fitted Otsu values, masks and frame metadata.
The [inspection ledger](../../../outputs/blocker_remediation_20261102/maldi_xenium_review.json)
records all 65 dispositions and all 62 inspected panels. Specific concerns include
internal rectangular image gaps, small/coarse fields, separate fragments and
peripheral signal retained by the current mask. None of these masks was approved
as a tissue boundary merely because the reader completed.

## Display scaling and denominator interpretation

The [32-section plot index](../../../outputs/blocker_remediation_20261102/maldi_scaling/index.html)
shows linear and `log1p` displays of the same 48-feature exported intensity sum,
plus the existing mask's excluded zeros. Original instrument coordinates and
50 µm sampling are retained. Color limits are section-specific and shown on color
bars. These are display transforms, without scientific normalization or a new
intensity calibration. Source bindings and counts are in its
[manifest](../../../outputs/blocker_remediation_20261102/maldi_scaling/manifest.json).

The existing `sum > 0` mask excludes all 12,486 measured all-zero spots. The original
all-raster versus masked detection comparison remains available. New integration
columns make these denominator choices inspectable; an independently justified
tissue mask and any changed scientific denominator still require specification.

## Verification and remaining requirements

The user reports the typical cutting order **Xenium → PhenoCycler → MALDI**,
with the qualifier "For the most part". The [separate source-context record](relationship_context.json)
retains the exact statement. For a confirmed three-section series following this
order, PhenoCycler is the immediately preceding MALDI neighbour. Exceptions,
exact acquisition/section mappings and section spacing remain unspecified. This
context leaves the inventory's separate same-slide declarations, exact pair
assignments and registration acceptance unchanged.

Thirteen focused Xenium source tests cover identity, relocation, channel, axes,
calibration and unavailable-source behavior. Thirty-three object/unit/Stage 08
tests cover exact-frame source receipts, stale accepted artifacts, native IDs,
overlap, finite denominators, invalid geometry and missing registrations. These
tests use labelled synthetic fixtures; they are not real cohort registration
evidence. The broader regression results are recorded separately in this bundle.
The final receipt guard also checks unique experiment bindings for all 65 actual
inventory rows against the preserved source audit. It does not repeat the raw
image reads or relabel them as fresh checks.

To produce real object summaries, supply the exact inspected MALDI-to-neighbour
acquisition rows, review the available outlines/masks and transforms against the
existing registration contract, and bind any accepted result to those inputs.
Cutting order or slice spacing, if available, should be recorded as serial-section
provenance. Donor equality, filename similarity and typing agreement do not supply
that evidence. No new acceptance cutoff was invented in this remediation.

Runnable commands and source schemas are documented in the
[workflow guide](../../PORTABILITY_REMEDIATION_USER_GUIDE.md).
