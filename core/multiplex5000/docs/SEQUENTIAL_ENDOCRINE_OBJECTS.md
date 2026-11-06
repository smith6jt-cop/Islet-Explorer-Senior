# Sequential islets, endocrine clusters and singles

The all-size method is `multiplex5000.integration.endocrine_workflow`. The earlier
`xenium_islets` command remains available to replay its historical, islet-only
centroid-hull generation. Do not use its ten-seed minimum for the all-size analysis.

## Live PubApp binding

On 2026-11-06 the running PubApp was verified as
`split-dev-20261105-plot-spatial-v2`, with data `split-dev-data-20261105`, app commit
`20baef5f4fceb193a46a8bff17fb426a263ef947`, and selected data producer
`378843e06c0cd66d8dfff5e3b60664f8a45569ee`. The scientific constructor is
`1b_islet_outlines_to_cells.groovy`, SHA-256
`533331b7ec1d39a6fd35062daad46152163e78829c322b0f1266b8b523aa483a`.
That hash matches the deployed data producer and the inspected CODEX source.
It incorporates the October 16 whole-cell pickup rule. The image-proposal
script also matches the live workflow: the CODEX port differs from IO60 only
in its header, platform pixel calibration and CD99-to-CHGA alias adaptation. The older IO60 October 2
export is not equivalent to this constructor.

The command executes those exact Groovy bytes in QuPath 0.7.0, with the pinned
common library, rather than approximating the ownership and splitting rules in
Python. Fresh scratch projects never read or save source qpdata. An image server
supplies construction calibration only. RNA geometry is in measured micrometres;
its scratch coordinate scale is one micrometre per geometry unit. Actual-image
review uses its own Xenium DAPI image and measured OME calibration.

The live PubApp has no Xenium proposal producer. The RNA proposal adapter below
is a declared extension feeding the same constructor. The PubApp's CODEX donor
cohort and cell calls are not substituted into the separate IO60 cohort. The
PubApp counts singles but omits them from aggregate islet rows. Native Single
construction is included here, but no Single is paired across serial sections.
The live app does not supply independent cell-to-cell correspondence evidence.

## Proposal and construction rules

Protein proposals are the original image-classifier outlines, recomputed against
native nucleated cell segmentations. Old exported Islet/Cluster/Single labels and
old `Whole cells` measurements are removed before construction.

RNA proposals use the pinned R3 calls joined by exact acquisition/native cell ID:
Endocrine lineage, lineage analysis eligibility, positive naming-gene evidence,
positive RNA depth and a measured nucleus. No post-Xenium IF is required. DBSCAN
uses 50 µm and `min_samples=1`, retaining isolated cells and sparse groups. The
proposal is the convex hull of the group's actual segmented cell areas; a
singleton uses its exact cell shape, preserving nonconvexity and holes. There is
no inferred circle fallback. All touching eligible nucleated cells supply
construction context; proximity does not assign endocrine cell identity.

The pinned constructor:

- picks a whole cell at at least 50% overlap with a proposal, resolving competing
  proposals by largest overlap; an otherwise empty proposal can pick its largest
  overlapping unclaimed cell;
- removes unowned cell areas and closes outlines around whole owned cells;
- applies its 0.63475 µm whisker opening with the original 90% retention safeguard;
- separates connected pieces while keeping parts of one segmented cell together;
- makes a one-cell object an `EndocrineSingle` with that cell's exact constructor
  geometry; with two or more whole cells, uses the mean of minimum enclosing
  circle diameter and minimum caliper width: below 35 µm is `EndocrineCluster`,
  at least 35 µm is `Islet`.

Construction counts include all owned cells, not only endocrine RNA seeds. An
RNA component without any seed is retained as a diagnostic and excluded from
correspondence. Stable structure keys derive from sorted native member IDs;
random QuPath object UUIDs are retained separately. RNA UUID translation is
recorded in the native cell sidecar and never replaces the original vendor ID.

All positive-area relationships to original proposals are retained after
splitting, including competing parents. Do not invent a unique parent where
more than one original region supports the final geometry.

## Independent checks and scientific limits

Independent area checks require unique cell ownership, no cut cells, agreement
with construction counts, and the same diameter category. Full-precision WKB
also checks that each Single equals the actual QuPath input cell exactly. The
raw vendor-to-QuPath discrepancy is reported separately: QuPath 0.7's native
GeometryTools uses a 0.01-unit precision model with float32 coordinate storage.
The coordinate tolerance is derived from that implementation, not fitted to
agreement in this cohort. It does not excuse a changed single-cell outline.

Before correspondence, render every section's native objects on actual DAPI,
with whole-tissue views and deterministic median-category, largest-object and
lowest-DAPI examples. Quantify image coverage and DAPI support per object. A
coarse-preview zero is a review flag; it is not proof of an absent nucleus.
Large RNA hulls can bridge internal structures, and eligibility/naming evidence
is not independent histological confirmation. Record seed membership changes at
40 and 60 µm; this is grouping sensitivity, not a full boundary sensitivity run.

## Correspondence and release scope

Retain all three categories in native outputs. **Exclude EndocrineSingle on both
sides before assignment**, controls and correspondence release statistics.
Neither anatomical context nor similar cell geometry establishes that isolated
cells on two sections are the same cell. Single-to-Islet and Single-to-Cluster
transitions are excluded too. Filtering completed assignments afterward is
incorrect: Singles would already have competed for multicell targets. The
`multicell_only_no_single_identity_v1` policy applies this restriction before
every primary and sensitivity assignment. Native Single counts remain available
per section; regional abundance comparisons require independently validated
tissue regions and comparable sampling denominators, not paired single-cell IDs.
Any future Single pairing requires a separate workflow with independent evidence
of cell identity, unavailable in this run.

Reuse the existing, checksum-pinned PASS/WARN
transforms; exclude failed or missing transformations. Do not refit a transform
to maximize endocrine agreement. The primary one-to-one candidate assignment
weights position 25% and geometry 75%, with 200 µm distance and 2.5 area-ratio
limits, an unmatched option, and no projected-overlap term. Position weights of
50% and 75% are sensitivity comparisons. Islet-to-Cluster and Cluster-to-Islet
transitions are allowed because adjacent sections can change a multicell
structure's observed size and cell count. These remain provisional tissue-area
relationships; small-cluster identity is not established by a good context score.

Anatomical evidence combines actual tissue context and focal-excluded
neighborhoods of up to six Islet-category anchors within 750, 1500 and 2250 µm.
Small objects are not validated by their own dense population. Compare relative
neighbor vectors with a 200 µm tolerance, requiring at least three anchors per
section. Controls use other targets 300–2000 µm away, the same target category,
area ratio within 2.5, and anchor counts within two. The 199 control-median draws
are descriptive negative controls, not a calibrated false-discovery estimate.
Use fixed median-x/y spatial blocks and report each RNA category separately.

The prespecified exploratory regional rubric requires prior tissue NCC z ≥5,
at least ten contextual queries, actual median above the control-median 95th
percentile at two of three radii, and three informative positive spatial blocks.
Release additionally requires support within the reported RNA object category;
a pooled signal cannot rescue failed category-stratified controls. This is a
conservative release safeguard, not a changed assignment parameter.
Pair-level ambiguity, image sampling flags, proposal sensitivity and category
changes remain visible. This rubric does not validate each pair, establish
same-cell identity, or activate a production release. Historical 473-pair evidence
remains a separate generation.

## Reproduction

Install the integration extra, including `ijson`, and use QuPath 0.7.0. The
manifest schema is `sequential_whole_endocrine_v1`: pinned constructor/common
sources, a live-workflow binding, and explicit sections with acquisition/native
source paths, SHA-256s, nucleus and cell geometry, calibration and scratch image
entry metadata. It refuses an existing output directory or an unreviewed
constructor generation.

```bash
python -m multiplex5000.integration.endocrine_workflow \
  --manifest /path/to/whole_endocrine_manifest.json \
  --output /new/output/directory \
  --qupath /path/to/QuPath
```

Construction writes native objects, native cell memberships, seed memberships,
40/60 µm grouping sensitivity, precise constructor and input-cell geometry, logs,
source pins and independent checks. Actual-image review and correspondence are
separate, required steps. The 2026-11-06 run bundle records the executed adapters,
reassessment code, images, source manifest and results under
`/home/smith6jt/Documents/papers/consolidation/sequential_all_endocrine_20261106`;
large native geometry lives under `/data/tmp/mplx5k/sequential_all_endocrine_20261106`.

## Completed 2026-11-06 assessment

All 39 native constructions completed: 27,333 RNA objects on 25 sections and
38,437 protein objects on the 14 nonfailed registrations. The 533,329 RNA seeds
were retained; 18 split RNA components without a seed remain diagnostic-only.
Every object passed whole-cell ownership/count, diameter-category and Single
geometry checks, and has positive-area original-proposal provenance. All 39
actual-image overview/category/outlier sheets were reviewed.

The first correspondence attempt incorrectly allowed Singles to compete in
object assignments. Its 6,644 candidates and 1,271 review candidates are
superseded. Calling Single pairs "regional" did not establish their identity.
That attempt contained 3,128 Single-involving pairs, including 305 in its review
subset; none is accepted for correspondence release.

The corrected run is
`/home/smith6jt/Documents/papers/consolidation/sequential_multicell_correspondence_20261106`.
It reuses unchanged native boundaries and transforms, excludes Singles before
matching on both sides, and recomputes assignment, sensitivity, controls and
release statistics. It contains 3,659 multicell candidates, 2,569 with zero
projected overlap, and 1,009 final candidates for exploratory regional review.
Compared with the former multicell subset, 3,440 pairs remain, 76 are replaced
and 219 are added. This is not a post-assignment deletion of Single rows.

All six original warning donors (6374, 6482, 6493, 6521, 6523, 6551) retain
regional support for both Islet and EndocrineCluster categories. Across the
expanded cohort, 12 of 14 donors qualify for exploratory regional review;
6461 and 6534 remain diagnostic-only and now both fail the pooled rubric.
No donor is accepted for Single correspondence. Native counts remain 10,907
RNA Singles across 25 sections and 16,452 protein Singles across 14 sections.
These are descriptive within-section counts, not validated matched-region
density estimates or cross-section cell identities.

The retained secondary local-alternative diagnostic checks other same-target-category objects
within 200 µm, without the area gate: a closer alternative within 0.05 neighbor
score of the assignment, or an alternative ahead by more than 0.05, triggers
review. Zero/0.05/0.10 score tolerances remain available as sensitivity columns.
This post-image-review flag changes no primary assignment or donor statistic.
These are review candidates, not individually validated matches.

A fresh portable replay of donor 6436 reproduced all 627 RNA and 3,245 protein
objects, categories, native cell memberships and geometries exactly (zero
symmetric-difference area). Other sections were executed and exhaustively
checked once in the final generation; they were not each independently rerun.
The prior historical geometry and 473-pair evidence retain their original hashes.

The audit command `scripts/audit_sequential_endocrine.py --run /path/to/run`
requires the completed native image-review receipt and the explicit
`multicell_only_no_single_identity_v1` policy in `analysis_plan.json`. It refuses
cached decisions from the previous all-category policy. Use a fresh directory
with `run_specs.json` pointing to the verified native outputs. The run bundle also contains
the image review, category-specific release filtering, report and packaging
programs. No production deployment or named cell-typing release was changed.

Run `scripts/check_sequential_local_alternatives.py --run /path/to/run` after
the primary audit and before final release-table assembly. The supplemental
receipt states that this diagnostic was added after visual review.
