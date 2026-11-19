# Xenium-derived islet boundary estimates

For the current all-size method, including endocrine clusters and singles, use
[Sequential endocrine objects](SEQUENTIAL_ENDOCRINE_OBJECTS.md). This document
describes the retained islet-only replay method.

`multiplex5000.integration.xenium_islets` supplies an executable RNA-derived
boundary path for sections without post-Xenium IF. Apply the RNA method across
the selected cohort to keep region definitions consistent. Available independent
IF outlines provide separate audit context; IF is not an input dependency.

The November 6 candidate run covers the 25 exact ordinary pancreas acquisitions
from the October 6 inventory, including all six WARN sections previously
unavailable solely because IF masks were missing. Its report and source pins are
under `/home/smith6jt/Documents/papers/consolidation/xenium_islet_boundaries_20261106/`.
This roster is not all Xenium images in the wider project: same-slide,
Xenium-only and other acquisitions need their own declared source records.

## Method and boundary meaning

The method reuses `structures.dbscan_labels`, with the existing primary defaults
`eps_um=50` and `min_samples=10`. A seed is an eligible RNA cell with an Endocrine
call and at least one detected naming gene for that call. Positive measured
nucleus count and positive RNA depth are also required. The adapter uses the
explicit `lineage_analysis_eligible` field; cell-core argmax is not used as an
identity call. It does not fall back to a broader seed set if evidence is absent.

The present candidate uses the pinned October 8 R3 review-generation calls.
Their review status remains attached to every object. This does not certify
those calls as current production typing or rebuild RNA expression/typing.
The library function `cluster_boundaries` also accepts explicit Boolean RNA
seed evidence from another declared producer, with the same nucleus/depth rules.

Each DBSCAN cluster becomes the convex hull of its seed centroids. This exports
the geometric support already summarized as convex-hull area in the existing
caller. The result is an **estimated boundary**, not a fluorescence segmentation
or a measured cell membrane. Convex hulls can bridge gaps, omit nonconvex edges,
and underrepresent an islet whose endocrine cells are not detected. DBSCAN's
radius is a local neighbor distance, not a maximum permitted islet diameter.
Degenerate point/line clusters remain unavailable; no artificial area is added.
Small/noise seed groups remain explicitly unassigned. No peri-islet expansion is
introduced.

Stable native cell IDs order DBSCAN traversal. Object IDs bind the acquisition
and sorted seed membership; seed memberships and the source file hashes preserve
the underlying segmentation evidence. Coordinates remain in the declared native
Xenium micrometre frame. An exported GeoJSON in micrometres is for geometric
analysis; it is not implicitly a Xenium Explorer image-pixel annotation.

The 40 and 60 micrometre radius profiles are sensitivity results. They are not
fitted against IF, disease groups, expression correlations or desired islet
counts. The historical IF masks had a 50 micrometre minimum equivalent diameter;
RNA clusters use a seed-count minimum, so their raw object counts are not directly
equivalent. The assessment reports RNA objects above that diameter separately.

## Execute and consume

From a checkout with the integration dependencies installed:

```bash
PYTHONPATH=core/multiplex5000 .venv/bin/python -m multiplex5000.integration.xenium_islets \
  --manifest inputs.json --output new_candidate_directory
```

The manifest schema is `xenium_islet_boundary_inputs_v1`. It declares `profiles`
(named `eps_um`/`min_samples` settings) and a `sections` list. Every section names
its donor, exact source serial, ROI, acquisition ID, micrometre coordinate units,
source release ID/status, vendor cells and RNA calls paths, and both SHA-256s.
The supplied review adapter expects the R3-style RNA evidence fields described
in `read_section`. Hashes, call-table acquisition fields, native ID uniqueness,
native-key joins and nucleus counts are checked before any clustering.

If upstream RNA QC omitted vendor cells, the manifest must explicitly declare
`calls_scope: declared_subset_of_vendor_cells`. Those raw cells stay in the input
accounting, with unavailable RNA call/depth evidence; they cannot become seeds.
Call IDs absent from the vendor cells always fail. Unknown cell labels are never
filled by matching coordinates or copying a neighboring cell's call.

For each profile/serial, outputs are:

- `objects.parquet`: WKB boundary estimates, source, frame, centroid, area,
  seed count and explicit status.
- `boundaries_um.geojson`: the same nondegenerate estimates in native micrometres.
- `seed_membership.parquet`: native seed IDs, coordinates, source and cluster/noise
  membership; no cross-assay identity transfer.
- `qc.json`: nucleus/depth admission, missing-call accounting and cluster counts.

The top-level summary and run manifest bind all output hashes. An existing run
directory is refused rather than overwritten. The November 6 assessment also
exports `registered_RNA_candidates.parquet` using the 14 existing non-FAIL
transforms, with native geometry, transform hashes and historical verdicts kept
beside the registered estimates. These are registration-review candidates;
six FAIL transforms and five absent transforms are not silently admitted.

## Validation and release limits

Analytic tests cover known hull area, native seed coverage, deterministic IDs,
input-order invariance, acquisition separation, nucleus/depth exclusion, empty
and degenerate cases, hash/serial checks, declared QC subsets and export safety.
The real-data assessment independently recomputes every nondegenerate hull area
using SciPy, verifies seed containment and all source/output hashes, reports
within-section hull intersections and compares cluster membership across radii.
No projected inter-section overlap threshold defines correspondence.

Where IF is available, the assessment compares fixed-method RNA estimates with
the independent outlines using position and geometry. Those candidate links and
figures are method evidence, not reviewed truth. IF is not used to tune seeds or
parameters. RNA-defined region membership is coupled to endocrine detection, so
endocrine enrichment inside those regions is not independent validation of the
boundaries or of cross-assay endocrine agreement.

The 473 historical IF-based pairs remain unchanged. Their release assessment
does not automatically transfer to RNA-derived pairs. A new boundary can remove
a missing-geometry prerequisite, but cannot establish a missing physical
relationship, repair a failed transform, certify individual islets or activate
a production release.

Primary implementation references:
[scikit-learn DBSCAN](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html)
and [Shapely convex hull](https://shapely.readthedocs.io/en/stable/reference/shapely.convex_hull.html).
