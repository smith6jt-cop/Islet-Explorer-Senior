# Anatomical context in sequential-section correspondence

The November 5 audit is in
`/home/smith6jt/Documents/papers/consolidation/sequential_anatomy_20261105/`.
Run `run_anatomy_audit.py` with the repository `.venv/bin/python`, then
`verify_anatomy.py`. `analysis_plan.json` records the evidence rubric before
the new neighborhood results. Historical sources and transforms are read-only.

## What is assessed

`multiplex5000.integration.serial_context` keeps focal-object geometric
similarity separate from neighboring anatomical arrangement. It does not use
projected intersection or biological feature agreement to select candidates.

For each historical pair, exclude both focal objects. Take up to six nearest
independent islet objects in each section within 1,500 µm. Express their
positions relative to their own focal object, and compare those vectors using
a one-to-one assignment. An agreeing neighbor contributes
`max(0, 1 - vector_error_um / 200)`. The score is twice the sum divided by the
total number of observed neighbors on the two sides. It is unavailable when
either side has fewer than three neighbors. Missing/extra neighbors are allowed;
this score is a geometric diagnostic, not a correspondence probability.

The focal object's area, eccentricity, overlap and biological features are
excluded from this context score. It therefore supplies evidence not used in
the focal geometry assignment. It is not independent landmark truth: nearby
queries share neighbors, and the historical candidates were already selected
by position, size and eccentricity.

Wrong-location controls use actual target objects 300–2,000 µm from the source,
within the same 2.5-fold area-ratio bound and with a similar neighbor count.
Their local structure is retained. Up to 20 controls per query and 199 sampled
control-median replicates provide a conditional reference, not a calibrated
p-value or false-match rate. Results are examined in four spatial blocks and
at 750, 1,500 and 2,250 µm neighborhood radii. Absolute placement remains
supported separately by the prior image/structural registration evidence;
relative vectors alone cannot establish the correct absolute transform.

The independent structures available throughout the six WARN donors are islets.
Paired, independently reviewed duct/vessel annotations were not recovered for
those donors. Existing label-derived duct/vessel tables on the two PASS donors
are not treated as truth or used to fill that gap.

## How geometric similarity affects correspondence

Intrinsic descriptors are area, exact area-moment eccentricity, solidity
and circularity of the independently segmented masks in the fixed historical
physical coordinate frame. Shape descriptors are invariant to translation and
rotation; area remains a separate size measurement. The intrinsic cost is:

`0.70 * abs(log(area_ratio)) / log(2.5) + 0.20 * eccentricity_difference`
`+ 0.05 * solidity_difference + 0.05 * circularity_difference`.

The total candidate cost combines `distance / 200` with intrinsic cost. Fixed
position/geometry weights are 75/25, 50/50 and 25/75. Projected overlap has
weight zero in all cases. The distance limit remains 200 µm and the area-ratio
limit 2.5. One-to-one assignments permit unmatched objects and report cost
ambiguity; they do not force all objects to acquire a partner.

The three profiles are sensitivities, not fitted or calibrated model weights.
Their candidate assignments are evaluated by the neighbor score, which is
excluded from their cost. No profile is optimized on biological findings and
no alternate assignment silently replaces a historical pair. Native regional
feature comparisons continue to use equal pair weights; correspondence costs
must not be presented as expression weights or probabilities.

## Release scope and limitations

The audit supports all six WARN registrations for a qualified historical,
regional exploratory release. All four spatial blocks show positive context
advantage in each donor. Five WARN donors clear the control reference at all
three radii; 6551 clears it at 1,500 and 2,250 µm, with too few evaluable queries
at 750 µm. Donor 6374 has the weakest absolute neighborhood agreement and the
largest median neighbor-vector discrepancy, so its pair-level cautions matter.

Seven historical pairs across WARN donors 6374, 6482 and 6523 have a local
alternative with a higher neighbor score. These and other sparse or uncorroborated
local contexts total 51 WARN pairs enumerated for review. These flags do not
justify rejecting a whole registration or silently asserting a different islet
identity. The 473 historical
candidate pairs and all native measurements remain intact, including the 117
zero-overlap candidates. No projected overlap is required for correspondence.

Registration-level regional acceptance is not individual islet/cell validation,
typing acceptance, proof of biological concordance or a production gate change.
The earlier provenance and typing qualifications remain. Independent landmarks
and review of disputed local structures would support stronger accuracy claims.

## Tests

Analytic checks cover exact shape moments, coordinate invariance, focal exclusion,
perfect neighbor agreement despite zero focal overlap, missing evidence,
one-to-one neighbor links, morphology weighting, deterministic assignment and
unmatched objects. The real-data verifier preserves every historical measurement,
checks inputs and candidate gates, and independently reproduces selected neighbor
scores by exhaustive enumeration rather than calling the assignment implementation.
