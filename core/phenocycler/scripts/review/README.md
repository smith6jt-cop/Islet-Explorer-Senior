# Current review lane

The active review lane contains two mutable comparisons of upstream QuPath
measurements. It does not call markers or create cell types.

> There are no verified production target/reference pairs. Historical RESTORE
> pilots, probabilistic training/promotion, marker-debug/oracle, Otsu-floor,
> divisor, and NNMF pair-web scripts are superseded and are not current command
> paths.

Run commands from `core/phenocycler`. Development outputs live under the
repo-root `data/development`, outside the content-addressed run tree.

## Active roster

`phenocycler/fundamental_pair_diagnostics.json` ships with exactly two
comparisons:

| Target | Reference | Source semantics |
|---|---|---|
| CD3e | Ker8_18 | raw QuPath cell means |
| Iba1_norm | Ker8_18 | exact signed QuPath Iba1 `Mean (norm)` versus raw Ker8_18 |

The roster is editable during development, subject to exact source, no-self,
no-duplicate, and one-shared-reference checks. It is not a production registry.

CD79a is excluded as a marker-pair candidate because its expression is too sparse
in normal tissue. CD68, MPO, EpCAM, E-cadherin, CD31, and Pan-cytokeratin are
also deliberately absent.

## Refresh the narrow QuPath projection

```bash
python scripts/review/qupath_fundamental_inputs.py \
  --donors 6436 6450 6566
```

This scans the declared full QuPath measurement exports and atomically replaces
the tool-owned cache at:

```text
../data/development/qupath_fundamental_inputs/
```

The projection reads only exact fields needed by the diagnostics, including
`Cell: Iba1: Mean (norm)`. Raw Iba1, convenience `*_norm` aliases, and
`*__restore_normalized` outputs cannot substitute. Projection schema v3 also
binds every partition to its declared donor, image, row count, numeric storage,
and unique object IDs. Existing v1/v2 development caches must be refreshed
with the command above; they are not silently upgraded by consumers.

## Summarize inputs

```bash
python scripts/review/fundamental_pair_survey.py \
  --donors 6436 6450 6566
```

This reports coverage, missingness, signed/zero counts, native-scale quantiles,
and Pearson/Spearman association. It creates no positive column, threshold, or
assignment.

## Rank-only overlap

```bash
python scripts/review/reference_overlap_audit.py \
  --donors 6436 6450 6566
```

For each comparison and donor, the audit selects deterministic exact-rank upper
tails at 0.90, 0.95, 0.99, 0.995, and 0.999 and reports symmetric overlap,
Jaccard, and independence lift. UUID breaks score ties deterministically.
Spearman rank correlation is descriptive; native-scale Pearson is deliberately
absent from this rank-only artifact.

The Iba1/Ker8_18 comparison intentionally mixes scales. It is compatible
because each marker is ranked independently; a common numeric boundary is
never applied. QuPath's signed Iba1 normalization is
`(raw - tissue p50)/(tissue p99 - tissue p50)`, so `1` is the tissue p99 anchor,
not a biological positivity threshold.

## Boundary/contact crops

`boundary_contact_audit.py` may be used to inspect how segmentation contact and
pixel morphology affect the two comparisons. It selects Iba1 objects using
the exact QuPath-normalized measurement but displays the acquired raw Iba1
image channel. Its crops are qualitative development evidence, not labelled
training data.

## Render the current figures

After refreshing the projection, survey, rank overlap, and donor 6450 boundary
crops with their current commands, render the version-gated development figures:

```bash
python scripts/review/fundamental_pair_figures.py
```

The pinned review output uses q70. For an explicit sensitivity render only,
pass, for example:

```bash
python scripts/review/fundamental_pair_figures.py \
  --cluster-selection-quantile 0.60 \
  --label fundamental_pair_figures_q60
```

A non-default quantile is rejected under the default workspace label so it
cannot masquerade as the current q70 output.

The tool writes PNG and PDF versions below
`../data/development/fundamental_pair_figures/`:

- `fundamental_pair_joint_distributions`: Figure 7-style donor-by-pair native
  joint log-density panels with review-only two-cluster line marginals on
  both axes and no inferred baseline guides. Left and bottom marginals use
  identical fixed physical depth, padding, line normalization, and rotated
  styling around a physically square joint box;
- `rank_tail_overlap`: overlap and independence lift across the five exact-rank
  upper tails;
- `crop_montage_6450`: a compact RGB overview plus a two-page PDF retaining
  the full raw RGB/target/reference/DAPI composites.

This is mutable, development-only evidence. Joint panels use native-linear,
donor-specific display windows: raw measurements span 0 through q99.5 and
signed Iba1 spans q0.5 through q99.5. Cells beyond both-axis windows are filtered rather than
clipped, every panel reports the displayed fraction, and all annotated
statistics still use the full finite-pair population. These limits are
visualization aids, not inferred baselines, common numeric thresholds, or cell
calls. No q99 guide is drawn on the joint panels, and no smoothing or value
replacement is performed.

The marginal populations are a review-only, unverified approximation of the
RESTORE Figure 7 grouping, not production RESTORE output. Exact donor-local
q70 upper ranks select the union of cells in either marker's top 30%, with
`object_id` breaking ties. Q70 is pinned as the closest reciprocal longitudinal
analogue of the manuscript grouping; it was not optimized against these
panels. Marginal-line separation and overlap depend on the clusters fitted from
this selected union and can change with the selection quantile. The manifest
records the quantile and target, reference, intersection, fit-union, and
outside-fit predicted fractions for every panel. Each marker becomes a
continuous donor-local bright
percentile from 0 to 1 over every finite pair; the brightest value is 1. KMeans
therefore receives reciprocal, dimensionless rank features rather than mixed
native scales. Two-cluster KMeans fits every q70-union cell in canonical marker
column and ascending `object_id` order, then predicts every finite cell into
exactly one of the two fitted clusters. Thus cells below both q70 fit arms still
contribute to the red/blue native-axis marginal curves; no artificial q70
support edge is imposed on those curves. This follows the manuscript's
two-group construction, which assigns its double-low Y3 population to one of
the two groups rather than drawing a third marginal population. Red denotes the
higher target-percentile-center
cluster. Blue denotes the higher reference-percentile-center cluster only when
it is distinct; otherwise the panel is visibly marked `NON-OPPOSING`, red is
higher on both axes, and blue is the other cluster. Both all-finite predicted
clusters are projected without tail truncation onto native-value marginal axes.
This assignment is exploratory and does not
verify a biological population, threshold, baseline, or call. Crop examples
are deterministic rank extremes and deliberately nonrepresentative. For Iba1,
joint-distribution and crop selection use exact signed QuPath
`Cell: Iba1: Mean (norm)`, while the
crop display uses the acquired raw Iba1 channel. The figure tool has no
input-refresh option and refuses stale projection, survey, overlap, or crop
method versions rather than silently plotting them.

To keep the joint landscape compact, that figure does not carry an explanatory
paragraph along its bottom edge. Its legend is deliberately limited to
`red: y-axis marker` and `blue: x-axis marker`; the review-only clustering and
native-axis caveats remain in this section and in `figure_manifest.json`. Rank-tail
and crop figures retain their own short explanatory captions.

## Safety boundary

All active tools:

- read only upstream QuPath inputs and explicitly supplied geometry/scratch
  data;
- preserve signed QuPath-normalized values;
- require unique UUIDs and aligned arrays;
- replace only a tool-owned mutable workspace;
- never read RESTORE calls or a production assignment;
- never write a run/stage manifest, reference control, divisor, marker call,
  state call, cell type, or QuPath assignment export.

The production run/manifests architecture remains intact for ingest, geometry,
REDSEA, and expression. Review outputs must never be copied into `data/runs/`.

## Historical scripts

Retired pilot, marker-debug/oracle, divisor, Otsu-floor, preprocessing, and
NNMF pair-web executables have been removed from the supported tree. Dated
handoffs record their historical measurements but are not runnable workflow
instructions.
