# Dynamic typing and multimodal portability evaluation

Evaluation date: 2026-11-02. Work branch: `evaluation/dynamic-portability-20261102`.
Scientific baseline: `maldi/serial-platform@00ec0322`. No scientific algorithm,
registry approval, production deployment or release activation was changed.
The production implementation addition is the project-only KINTSUGI observation
adapter. Other new code is read-only evaluation and evidence assembly tooling.

**Portability is partial.** Shared numerical routines exist, but protein biology
and measurement contracts are not dynamically interchangeable. IO60 RNA-only
execution is supported across the frozen 65-acquisition roster. The later spleen
PC67 producer is computationally reconcilable across eight sections. Same-slide
cell links and serial MALDI registration/object integration cannot be validated
for the full requested cohorts with the available specifications and accepted
geometry. Those acquisitions remain explicitly blocked rather than omitted.

Use the [workflow user guide](../../PORTABILITY_EVALUATION_USER_GUIDE.md) for every
route and recording operation. The [handoff completion record](HANDOFF_COMPLETION.md)
and [requirement matrix](handoff_requirements.json) map this delivery to the original
task and separate finished evaluation work from remaining scientific acceptance.

## Results and scope

| Route | Cohort evidence | Conclusion |
|---|---|---|
| Original CODEX → typing-port parity | All 20 images; 60 exact original-source/port comparisons pass; 3,199,606 eligible identities reconcile | Computational port parity established. Original fixture still rejects historical mtimes; 40 exact CSV assertions fail. Producer-thread/round-trip parser replay reproduces all 20 tables exactly; no criterion changed. |
| CODEX v20 | Independent check of 24 release bindings | 22 mismatches: all 20 identity tables, levels and structure binding. Historical reconciliation says pass for 20 images/3,199,606 analysed cells, but is not current-file compatibility evidence. |
| IO60 PhenoCycler pancreas / embedded pLN | All 25 primary tables built; 31,385,695 eligible non-pLN IDs reconcile | Refinement fails at 0112: unavailable Ker8_18 is treated as carried, leaving zero measurable cells. Four acquisitions have that unavailable marker. No completed Stage 06 receipt or assembled mixed output. |
| Dedicated pLN protein | Two existing pilot sections, 832,062 cells, 786,255 eligible | Native/qualified IDs and eligibility reconcile. Original manifest recovered at its exact producer hash; current manifest differs. Candidate evidence only; remaining 17 of 19 acquisitions blocked. |
| Spleen PC67 protein | All eight sections, 4,346,039 eligible cells | Every eligible native ID reconciles. Later producer `a9549cc`, independently audited; ten large image/geometry input hashes unavailable. |
| Spleen PC33 protein | All four declared sections; 9,785,239 source cells, 8,337,560 QC eligible | Source IDs/coordinates audited and all four spatial panels inspected. Lineage roots/registry remain unsupported; no substituted PC67 typing. |
| IO60 Xenium-only | All 65 acquisitions, 33,864,036 cells | 60 candidates, three damaged and two low-signal regions unavailable; 510 receipt input/output hashes match; every identity/count and independent call-state reconciliation passes. |
| Same-slide protein/RNA | Four explicit pancreas declarations and eight source-verified spleen pairs | No accepted Stage 07/08 receipts in the inspected release roots. Full-cohort routing, independent spleen RNA lineage typing and link uncertainty/hierarchy comparison remain incomplete. pLN physical pair declarations unresolved. |
| Serial MALDI | All 32 IO60 sections, 18 donors, 1,666,450 spots | Source pairing, IDs, values, mask, denominators and source-image comparison complete. Accepted neighbour/frame/registration declarations remain unavailable; object summaries blocked. |

The [acquisition ledger](acquisition_dispositions.json) accounts for **201 source
acquisitions**: 129 in the IO60 baseline inventory, 32 IO60 MALDI sections, 20 CODEX
images, 12 spleen protein acquisitions and eight spleen RNA acquisitions. The
129 include 20 post-Xenium IF acquisitions explicitly excluded from the requested
typing routes. This is an acquisition count, not a donor, unique physical-section
or independently replicated sample count. Unknown section, segmentation, panel
or count fields remain null. Donor equality creates no relationship.
Final acquisition dispositions are 90 completed evaluations, five unavailable,
20 excluded and 86 blocked. CODEX's 20 completed evaluations retain the failed
historical artifact assertions as separate results; completed is not approval.

The existing donor 6457 exclusion applies only to PhenoCycler; its RNA and other
eligible modalities are not removed. Damaged RNA regions retain their declared
cell counts and unavailable outcomes. The two RNA-only acquisitions
`xenium:0041323:Region_1` and `xenium:0041326:Region_1` both complete reconciliation
without protein measurements, protein normalization or same-slide pairing.

## What is shared, configured or unresolved

[dependencies-v2.csv](dependencies-v2.csv) gives 33 evidence-backed grouped dependency
traces and remediation IDs. [parameter_catalog.json](parameter_catalog.json)
enumerates 572 module settings/function defaults and 13 registry tables with
source hashes, lines and references. The completed
[execution trace summary](parameter_trace_summary.json),
[full trace](parameter_traces.json.gz) and
[literal inventory](parameter_literal_candidates.json.gz) account for all catalog
items, 341 registry rows and 11,218 additional literals with explicit heuristic
dispositions. Exact source/report/receipt bindings supply 1,710 reported values
across eight declarations; 541 have static consumers without recorded runtime
values and 23 have no resolved read in the inspected sources. These are explicit
evidence limits, **not proof of exhaustive dynamic parameter coverage**. CORE-1
specifies future runtime emission; scientific instrumentation was not added. The
[readiness matrix](readiness.json) separates code, valid inputs, scientific
specification, executable routing, computational validation and biological review
by modality, tissue, panel and relationship.
The original `dependencies.csv` is retained unchanged because earlier immutable
observations reference its hash; v2 contains the final follow-up conclusions.

Shared components include fitting mathematics, percentile/mixture calculations,
identity/receipt checks, matrix operations and geometry primitives. Their existence
does not validate reference-population choices, normalization units or biological
resolution. Whole-cell positivity, membrane differentiation, coverage and
despeckling must follow their exact producer contracts.

The legacy pancreas vector uses one cohort scale fitted from QC-eligible cells
in the recorded pancreas alignment block, with at most 200,000 cells per donor
(seed 0), a 512-bin positive log-intensity mode and p99.9 bright end. A missing
background fit falls back to floor zero; an invalid span becomes one during
scaling. These execution defaults are recorded as unresolved transfer assumptions.

The first isolated attempt stopped on a stale acquisition-run table; its replacement
was recovered from the existing v3 checkout and verified against all 25 exact image
names. The corrected attempt completed primary fitting and assembly for all 25,
then stopped at refinement before any donor completed its refinement pass. Independent checks confirm each primary table contains
exactly the eligible non-pLN native IDs; differences from the whole aligned table
are expected tissue/QC filtering. This does not establish refined typing validity.
The independent nucleus audit covers 39,394,852 aligned cells, 32,874,241 eligible
cells and all 31,385,695 primary cells. None of the admitted or primary cells lacks
positive, finite source nucleus area, source nucleus/cell ratio or raster nucleus
ratio. All 25 geometry-QC hashes and published admission totals match. The existing
5.0 µm² ingestion minimum accounts for 16,983 filtered raw objects. Raw export
size/mtime bindings match, but their full file hashes were not recomputed in this
selected-column audit; that limitation remains explicit in `protein_nucleus.json`.

Stage 03 declares Ker8_18 unacquired and Stage 05 declares it unavailable in 0112,
6347, 6533 and 6580. Refinement instead treats a present all-NaN column as a carried
marker. Its finite-value predicate then leaves no measurable cells, and NNDescent
raises on the empty matrix. All 25 admission masks and primary spatial diagnostics
were inspected separately; no cell was imputed or relabelled to bypass the failure.
See `protein_primary/` and `protein_primary_expected_universe.json` in the evaluation
root. PROT-1 specifies availability-aware admission and an explicit empty/small-set
outcome; it is not implemented here.

Refinement actually fits its representation across all 25 configured donors,
despite the route settings label `single_donor_representation_fit`. Its reader uses
the geometry-QC-eligible inner join without an additional alignment-block filter;
the resulting population differs from the cohort-scale population. Cofactors,
centres and spreads were emitted and recorded. Per-class/stratum refinement bounds
were never produced. The fixed K=30, concordance 0.40, fifth-percentile bound,
minimum 50 clear members and representation fallbacks remain version-specific.

Native-ID agreement does not imply identical coordinate conventions. Across all
31,385,695 primary cells, primary coordinates exactly match the frozen geometry,
while aligned coordinates exactly match Stage 03 and the source context. The
islet adapter deliberately rescales QuPath's 0.5078 µm/pixel to the manifest's
0.5077663810243286 µm/pixel. The observed primary/aligned ratio matches that
declared conversion, with maximum floating arithmetic residual 7.28×10⁻¹² µm.
Unconverted coordinates differ by up to 2.194369 µm. All 25 source-context hashes
match. This documents two coordinate conventions, not a wrong-cell identity or
a registration acceptance result. Consumers must bind a common calibrated frame;
no tolerance or coordinate rewrite was introduced. Full evidence is in
`protein_primary_coordinates.json`.

At the pinned typing port, hormone-derived reference populations, endocrine/
exocrine context, lineage-specific canaries and pancreas defaults still enter
execution. Registry rows parameterize some marker choices, but the tissue-hook
API and registry checks do not govern all `levels`/`identity` paths. Missing context
can become implicit tissue/no-structure values. P1–P9 in the PC67 proposal table
are all proposed at this pin. The later HPG method and its changed registry are
separate evidence; they cannot retroactively approve the old port. PC33, dedicated
pLN and PhenoCycler pancreas lack complete registries for that per-image route.
Two additional configuration-control gaps are concrete: the historical pLN receipt's
`protein_cutoff_fraction=0.1` is not passed to its evaluator, which reads
`REQUIRE_MIN`; the RNA naming-policy receipt is descriptive while the actual worker
argument hardcodes `shared_er_exact_v1`. Equal metadata and code values do not
establish dynamic configuration control. Exact producer source snapshots are
retained; current wrappers are not attributed to historical runs.

The legacy dedicated-pLN candidate code accepts configured section rosters, but
only two sections have the required published measurement/normalization artifacts
in the inspected release. The other 17 are blocked on complete source bindings
and upstream artifacts/review, rather than on an intrinsic two-section array
limit. This does not supply a pLN registry for the separate per-image port.

The newer spleen producer uses lineage-based references and a 6551 spleen
hierarchy. Its Stage 05 prerequisite does not mean aligned values are consumed:
typing uses raw REDSEA measurements and section-fitted levels. No cross-assay
agreement was used to select or tune methods during this evaluation.

RNA uses tissue-specific workers and actual panel coverage. The detection target is
fixed at 10,000 counts/cell, with `log1p(raw * 10000 / depth) >= log1p(1)`.
The pLN clustering embedding separately calls Scanpy `normalize_total` with no
target argument, using the median depth after nucleus/cluster QC. Its original
reports omit that fitted target. Reconstruction from bound inputs gives 58–564
gene counts/cell across the 30 candidate pLN regions; values and the inspected
normalizer implementation are retained in `remote/rna/pln_embedding_normalization.json`.
Clustering and naming parameters are tissue-specific; fitted clusters, seeds,
depth, available genes and controls are retained per section. The region-level
20-gene median screen and worker QC/refinement thresholds remain scoped to their
recorded method, with no new transfer tolerance invented. Historical donor
lists, frozen partitions and old labels do not govern the exact Stage 06 route.
For example, `0041323/Region_1` reports 5,101 panel genes, 485 wanted genes present,
GNLY absent and INS/GCG/SST/PPY absent. Its surrogate endocrine core does not
support alpha/beta/delta claims. Missing genes remain distinct from measured zeros.

RNA state totals reconcile independently: 23,519,988 specific calls, 2,618,820
broad states, 1,303,679 unresolved, 3,773,554 uncalled and 2,647,995 unavailable.
The admitted nucleus-eligible universe is 32,561,109; 1,302,927 cells have no
nucleus. Unavailable-region counts also include other cells in those regions,
so these are different accounting partitions. No disease effect is inferred
from the two RNA-only donors or from expected lineage frequencies.

## Geometry and serial evidence

The [pair ledger](pair_dispositions.json) preserves four exact pancreas serial
declarations, eight verified spleen source pairs and unresolved pLN relationships.
Repeated donor-note rows are not additional pairs. The spleen source verification
checks the section table hash, export header, required bundle files and crop
offsets. A further input audit reconciles all eight RNA matrices to 4,167,720
native cells (4,159,120 with a nucleus; 8,600 without). All eight have the same
5,001-gene set; generic identity-core coverage is recorded by gene. This panel
is not assumed equivalent to the IO60 panel. This establishes source correspondence,
not registration acceptance or a validated spleen RNA lineage hierarchy.
Protein and RNA are independently segmented; a joint typing stage alone supplies
no same-cell evidence.

The current same-slide pilot selects one region and returns mutual-nearest
candidate links within 5 µm. It retains protein candidates and RNA gene evidence,
not a complete independently typed RNA hierarchy. There is no full ambiguity,
split/merge or empirical registration-uncertainty model. Fault probes demonstrate
an unflagged identical-centroid tie and changing link availability at offsets
4.999, 5.000 and 5.001 µm. Those are synthetic scenarios, not new acceptance
tolerances. Missing frames, duplicate IDs, stale transforms and changed manifests
are separately exercised. No real candidate links or concordance claims are
manufactured in the absence of accepted transforms.

Both MALDI exports contain 2,288,540 rows with identical spot keys; the IO60 subset
has 1,666,450 spots and 48 features (38 glycans, ten glycogen). There are no missing,
negative or duplicate-coordinate values in this subset. Observed grid residuals
reach 0.00390625 µm; this is descriptive precision evidence, not an invented
acceptance threshold. Source intensity calibration/units beyond exported values
remain undeclared.

The mask uses summed finite intensity greater than zero. **All 12,486 all-zero
spots are excluded from the mask.** Finite zeros otherwise remain measurements,
distinct from missing. Both all-raster and masked detection fractions are saved
for every section and feature. The largest observed difference is 0.104622 in
`6461_02_PLN`, feature `Hex5dHex1HexNAc4..Na`: 10,954 detected spots out of 12,691
raster spots (0.863131), versus 11,319 masked spots (0.967753). The ID token is not
used to infer inspected tissue identity.

The follow-up [source audit](maldi_source_evidence.json) checks the supplier ZIP,
canonical CSV source, inventory and workbook cells. The ZIP contains only the two
CSVs; no cutting order or MALDI section-distance manifest was found. The dated
plan permits visual comparison with same-series reference slides. All 32 sections
were inspected in the comparison figure; [per-section observations](maldi_source_review.json)
retain component differences and ambiguous region correspondence. Repeated Xenium
whole-slide thumbnails do not establish region-specific frames, and independent
thumbnail sizing does not establish physical scale. These observations do not
populate production acceptance fields.

The [pLN object-source audit](pln_object_availability.json) exactly joins all 19
acquisitions to QuPath entries and 47 native annotation UUIDs: 44 polygons and
three LineStrings in 6534. The latter cannot provide area summaries. Many polygons
are bounding rectangles, not inferred tissue outlines or follicles. Calibration,
frame and source contracts for future consumers are specified in MALDI-3; no
MALDI object association is inferred from their availability.

The MALDI registration code searches orientation/mirror and constrains physical
scale at 50 µm sampling. Dedicated pLN outlines have a DAPI-mask route; dedicated
pLN object summaries do not have a complete contract. Xenium-neighbour
registration is explicitly unbuilt. Pancreas object summaries depend on existing
islet/structure exports. Eligibility uses erosion by registration RMSE + 25 µm,
while values summarize all inside spots; overlapping objects use first-object
precedence. These choices need explicit review. Distance remains continuous and
signed. There are no MALDI-derived cell identities, same-cell claims or new bands.

## Evidence and review limits

All execution/evaluation outputs are isolated under
`outputs/dynamic_portability_20261102/`; the full RNA and spleen read-only audits
ran on HPG compute nodes and their evidence was copied locally. Receipt files
retain original producer versions and paths. [pins.json](pins.json) separately
pins main, remote main, the working MALDI branch, typing port, original CODEX
method, v20 producer, later RNA/spleen producers and KINTSUGI. Producer and
evaluation environments are distinct and retained in receipts/manifests.
The CODEX roster is frozen from its source-project manifest, with the original
cohort CSV recovered at `fd3c4de` to match the manifest hash exactly. Its 20 image
IDs match the tested table roster; this inventory recovery does not resolve the
separate v20 binding failures. A further [CODEX recovery](codex_completion_summary.json)
verifies 95 historical output hashes and the journal's source/input chain. All
60 direct unchanged source-to-port computations pass exactly, and every original
identity assertion passes (3,199,606 eligible cells from 3,525,977 measured objects).
The original mtime fixture remains unchanged and still rejects output produced
before its source commit, despite matching source bytes. All 40 levels/canary
CSV assertions fail exact equality under their default parser/thread settings.
The recovered producer uses one BLAS thread; that environment with round-trip
CSV parsing reproduces all 20 positivity/differentiation/canary tables exactly.
This diagnostic explains reproducibility limits without weakening a tolerance or
editing a scientific method. Historic endocrine inputs have filesystem bindings,
not original content hashes, and historical package versions were not locked.

[artifact_manifest.json](artifact_manifest.json) hashes the review bundle and
records scientific receipt checks, unavailable hashes and explicit scope limits.
It also retains the failed run’s primary-table, alignment and panel-availability
bindings. `audit_corrections.json` records evaluation-reader/configuration
corrections alongside their retained original logs.

[inspection_ledger.json](inspection_ledger.json) records per-section diagnostic
inspection and representative cells/spots. RNA, protein and MALDI scatterplots
are spatial diagnostics sampled deterministically (up to 12,000 points per
section); reconciliation always uses all cells/spots. The 144 panels were inspected
for every audited section with a spatial artifact, including all 25 primary-only
protein tables from the failed run. This is not raw-DAPI validation,
biological acceptance or exhaustive per-cell visual
review. Unavailable sections retain explicit diagnostic panels where cell data
exist. Unsupported acquisitions are accounted as blockers rather than plotted
with fabricated results. The separate 32-section source-comparison review above
is additional to these 144 panels; it is not registration acceptance.

Existing tests cover absent measurements/references, mixed panels, weak RNA,
identity errors, unmatched cells, source/transform changes, MALDI pairing,
zero/missing handling, orientation and object-size uncertainty. Additional
fault-probe results are saved under `fault_probes/`. Computational fidelity,
geometric validity and biological applicability are separate conclusions.
Frequency targets, disease ordering, historical labels and assay agreement were
not tuning objectives. [remediation.md](remediation.md) gives prioritized,
implementation-ready requirements for remaining gaps.

## Parameter recording and acceptance

The populated [parameter history](parameter_history.json) uses KINTSUGI's pinned
SQLite engine through a project-scoped CLI and actual stdio MCP adapter. Its
351 observations occupy 499 append-only history records, with zero approved
records. All requested acquisitions have observation context; excluded post-Xenium
IF acquisitions remain outside recording scope. It
records proposed/evaluated/failed/blocked observations, fixed settings, selected
methods, fitted values, review criteria and hashed evidence. Missing quality
is null. No fabricated marker is needed for geometry, pairs or feature sets.
See [recording.md](recording.md) for connection, inspection, export and retry.

Nineteen recorder tests and actual MCP protocol/scientific-project checks pass.
The existing route suite passes 49 tests and the pinned port unit suite passes
37. The follow-up executes the unchanged CODEX numerical assertions through a
journal-validated evaluation harness: 20 identity passes and 40 exact CSV failures.
Its independent source-to-port comparison passes all 60 exact checks.
Idempotent retry, append-only revisions, restart/read/export, project/global
isolation, approval denial and failure durability are verified. Original
observations are retained when attribution is corrected. Automatic application
and recommendations are disabled; `evaluated` does not mean user-approved.
The global database was absent before and after the real recording runs; tests
also confirm an existing global-history sentinel remains byte-for-byte unchanged.

Evaluation infrastructure is implemented. Full scientific acceptance is **not
met**: CODEX historical artifact replay/v20 binding, PhenoCycler refinement availability,
unsupported panel/tissue routes,
accepted same-slide/MALDI geometry, unrecorded dynamic parameter values and
biological review remain explicit blockers. No acceptance criteria were weakened
to convert those blockers into passes.
