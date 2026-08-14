# Handoff — RESTORE intensity inputs, CD11b, and NNMF diagnosis (2026-08-12)

> **Superseded snapshot (2026-08-13):** this handoff is historical evidence,
> not an executable plan. There are no verified production pairs; ordinary
> execution stops after expression, and pair-dependent controls, calibration,
> typing, and export are blocked; pair-independent state calling is
> explicit-only. The only active pairwise work
> is the two-comparison, rank-only upstream review in
> `scripts/review/README.md`.

This handoff records implementation through the three-donor preprocessing
comparison and the subsequent NNMF fit-design diagnosis. It is not a
production-release sign-off. The working tree is intentionally dirty and
contains unrelated marker-debug work that predated this effort; do not reset or
delete those files.

## 2026-08-13 development override

This section supersedes every later statement that calls the six-target,
16-arm reference web or the legacy RESTORE-normalized split the active
diagnostic. Those experiments and pilot artifacts remain historical evidence
only.

The active troubleshooting roster is two unique, symmetric comparisons of
upstream QuPath inputs:

- raw QuPath `Cell: CD3e: Mean` versus raw `Cell: Ker8-18: Mean`;
- QuPath `Cell: Iba1: Mean (norm)` versus raw `Cell: Ker8-18: Mean`.

CD79a was removed as a pair candidate because its expression is too sparse in
normal tissue. CD68, MPO, and EpCAM were also removed. Reciprocal rows were
removed because a rank-overlap comparison already measures both sides.
E-cadherin, CD31, and Pan-cytokeratin remain absent. The roster is schema v3 in
`phenocycler/fundamental_pair_diagnostics.json` and is editable during
development subject to its exact QuPath source contract plus basic
no-self/no-duplicate/one-reference checks.

`Iba1_norm` in the review scripts now means only QuPath's exact signed
`Cell: Iba1: Mean (norm)` measurement. QuPath created it before RESTORE for its
object classifier as `(raw - tissue p50)/(tissue p99 - tissue p50)`. Therefore
`1` is the tissue p99 anchor, not a positivity threshold, and 35--42% of the
three donors' Iba1 values are legitimately negative. The old `>1` cell calls
and fundamental assignments were invalid and have been removed. Missing
QuPath norm fails explicitly; raw Iba1, retired `restore_gated_redsea/Iba1_norm`,
and `Iba1__restore_normalized` are forbidden fallbacks.

Production ingest now preserves all four QuPath norm compartments for exactly
Iba1, CD68, and M2Gal3 under collision-proof
`QuPathNorm__<Compartment>__<Marker>` names. The registry records
`measurement_source=qupath_norm` independently of REDSEA spillover policy.
Signed values remain unchanged through ingest and selected expression; at the
nonnegative RESTORE boundary, finite values below zero are explicitly
left-censored to zero and audited as `qupath_norm_left_censor_zero`. M2Gal3 is
ingested and selected but remains deferred until it has a validated reference.
Marker evidence records the source-generic `<marker>__selected_intensity`;
`corrected_intensity` is a superseded REDSEA-specific evidence name. A finite
selected value at zero is source-generically `censored_below` and unavailable.

The active development tools no longer have immutable artifact behavior:

- `qupath_fundamental_inputs.py` scans the full current QuPath exports into a
  narrow mutable projection, and `fundamental_pair_survey.py`,
  `reference_overlap_audit.py`, and `boundary_contact_audit.py` consume it;
- they do not read run/stage manifests, receipts, adoption lineage, content
  IDs, or file hashes;
- `boundary_contact_audit.py` selects Iba1 objects from the exact QuPath norm
  measurement while displaying the acquired raw Iba1 pixel channel, accepts
  scratch object IDs plus a structurally valid CSR as its complete contact
  contract, and reports partial coverage instead of suppressing output; and
- the old `STALE_SCRATCH`, exact-area/replay, and fresh-label gates are gone.

Basic correctness checks remain: exact source headers, required columns,
unique IDs, aligned shapes, finite analysis rows, deterministic ties, and CSR
structure. Production run/stage immutability remains, but its ingest,
expression, and RESTORE input contracts were deliberately version-bumped for
the QuPath norm source correction.

### Reduced development artifacts must be regenerated

The existing reduced outputs below `../data/development` were generated from
legacy RESTORE-normalized columns and are superseded. Regenerate the shared
QuPath projection first, then rerun the rank-overlap, input-summary, and crop
tools. The corrected survey emits no per-cell calls or labels, and the old
reported `>1` counts/lifts must not be cited as QuPath-input results.

## Outcome at the stopping point

The original divisor landscape mixed coordinates from two different intensity
spaces: a raw QuPath Otsu floor was drawn over a REDSEA-corrected distribution.
That overlay was dimensionally invalid. Production RESTORE also reconstructed
reference populations from target-frequency-dependent raw floors even though a
versioned `reference_controls` artifact already existed.

The candidate implementation now:

1. selects the marker registry's explicit measurement source rather than
   treating every marker as REDSEA-derived;
2. uses exact QuPath Cell/Nucleus/Cytoplasm/Membrane `Mean (norm)` values for
   Iba1 and CD68 in active RESTORE fits/calls, retains the same inputs for
   deferred M2Gal3, and uses the declared REDSEA source for other active
   markers;
3. fits NNMF on a deterministic sample of all complete
   `estimation_eligible` cells, with no Otsu-selected or frequency-balanced
   arms;
4. intersects the NNMF reference component with the controls-stage,
   target-independent reference-positive mask;
5. estimates the boundary with the registry marker-alpha finite-sample
   empirical upper-tail order statistic, after a fixed contamination screen;
6. preserves signed QuPath norm values through ingest/expression, then
   left-censors finite negatives only at RESTORE's nonnegative boundary while
   retaining NaN; a finite selected value at zero is
   `censored_below`/unavailable and excluded from threshold estimation;
7. emits the measurement source, boundary transform, left-censor count, scale
   bounds, control source, marker alpha, contamination, censoring, and empirical
   tail probability as provenance, with marker evidence named
   `selected_intensity`;
8. fails closed when a pair cannot resolve the tail, is contaminated, is not
   directional, or is unstable; and
9. removes CD11b completely from RESTORE and phenotype rules. CD11b remains a
   REDSEA-corrected `assay_qc` measurement only.

The production-config candidate is `linear_upper_scale`: every NNMF feature is
divided by its all-estimation-eligible 99.5th-percentile upper bound, with a
fixed zero anchor and no lower-quantile subtraction. Fixed p95/p99/p999/max
divisors and shifted/nonlinear transforms are now review-only in the production
configuration path.

## Historical three-donor preprocessing comparison

This comparison predates the QuPath norm source correction and is retained as
method-development history. Its REDSEA-only input statements must not be
applied to Iba1, CD68, or M2Gal3 in the corrected workflow.

Donors: 6436, 6450, and 6566. Pair web: 23 targets after CD11b removal, hence
69 donor-pairs. Everything except the feature preprocessing transform was held
fixed: REDSEA values, estimation universe, reference masks, three NNMF seeds,
component rule, marker alpha, contamination rule, and call rule.

Artifact:

- `data/pilots/preprocessing_survey_linear/preprocessing_survey.parquet`
- `data/pilots/preprocessing_survey_linear/preprocessing_survey.json`

Comparison: `linear_robust_minmax` versus zero-preserving
`linear_upper_scale`.

| Result | Agreement |
|---|---:|
| Pair state | 69/69 |
| Reference-control membership fingerprint | 55/69 |
| Final call fingerprint | 66/69 |
| Valid pair states | 52/69 under either transform |
| Invalid pair states | 14/69 under either transform |
| Model-unstable states | 3/69 under either transform |

Only three call masks changed:

| Donor | Marker | Robust calls | Upper-scale calls | Change | Change as % of eligible cells |
|---|---|---:|---:|---:|---:|
| 6436 | B3TUBB | 37,101 | 37,196 | +95 (+0.26%) | +0.0060% |
| 6450 | EpCAM | 32,869 | 31,908 | -961 (-2.92%) | -0.2611% |
| 6566 | B3TUBB | 66,310 | 66,482 | +172 (+0.26%) | +0.0097% |

The selected upper-scale rule therefore preserves every accept/abstain state
and 95.7% of exact call masks while retaining a meaningful zero. The 14 control
membership changes show that subtracting the lower population quantile is not
algebraically inert, even when most final thresholds happen to agree.

### Nonlinear isolation result

A four-transform run on donor 6436 is retained at
`data/pilots/preprocessing_survey_6436/preprocessing_survey.parquet`.

Relative to `linear_robust_minmax`:

| Transform | Same states | Same call fingerprints |
|---|---:|---:|
| `linear_upper_scale` | 23/23 | 22/23 |
| `log1p_robust_minmax` | 19/23 | 6/23 |
| `asinh_robust_minmax` | 19/23 | 5/23 |

Thus a channelwise linear rescaling is nearly inert, while log/asinh reshaping
materially changes the clustering. Nonlinear transformation is not authorized
for production by this comparison.

### Reference-mask provenance check

The current registry fingerprint invalidates masks stamped by the source run.
The masks were recomputed for the same three donors from immutable expression
partitions and written under:

- `data/pilots/reference_controls_empirical_alpha/reference_controls/`
- `data/pilots/reference_controls_empirical_alpha/reference_control_models.parquet`

For all three donors and all four references, the selected 20,000 object IDs
and every seed rank exactly reproduced the source-run selection. Only registry
provenance changed.

## Important unresolved finding

The preprocessing choice is isolated, but the all-cell NNMF design is not yet
acceptable as a released replacement: only 52 of 69 donor-pairs are valid.
The same 17 failures occur under both linear transforms, so shifting/scaling is
not their cause.

- Donor 6436: invalid CD79a, CD206, CD4, SST; unstable CD3e, CD68, INS.
- Donor 6450: invalid CD79a, CD206, CD8, INS.
- Donor 6566: invalid CD3e, CD79a, CD206, CD4, CD8, GCG.

The 14 invalid pairs fail directional mean/median separation. The three
unstable pairs fail convergence/divisor reproducibility policy. In many sparse
markers, the all-cell two-component fit assigns an implausibly large bulk
component as the target arm. The fit-design experiments below diagnose this
failure without hiding it by weakening `PairState` gates or converting
unavailable calls to negatives.

## Resumed NNMF fit-design diagnosis

The 17 failures were then diagnosed with a frozen 69-pair review harness:
`scripts/review/nnmf_fit_design_survey.py`. Every experiment retained the
selected `linear_upper_scale` REDSEA inputs, estimation universe, reference
masks, empirical-alpha boundary, exact-zero censoring, final call rule, and
fail-closed state gates.

Twelve of the 14 invalid pairs pass directional separation by mean but not by
median. A rare bright right tail moves the mean while most cells assigned to
the nominal target component remain background-like. The other two invalid
pairs, plus 6436 CD4, can select the same winning target/reference component or
have the wrong directional fold. This is a component-partition problem, not a
slightly misplaced final threshold. Of the three unstable pairs, 6436 CD3e and
CD68 have reproducible divisors but hit the 2,000-iteration ceiling; 6436 INS
converges but varies with the sampled cells.

| Isolated design | Valid / 69 | Baseline-valid losses | Main finding |
|---|---:|---:|---|
| rank 2, Frobenius/CD, 2k (baseline) | 52 | — | 14 invalid, 3 unstable |
| iteration ceiling 10k | 54 | 0 | cleanly rescues CD3e and CD68; all 52 existing valid call masks are identical |
| fixed sample 0 | 53 | 0 | rescues INS; all 52 existing valid calls and all 69 controls are identical |
| fixed sample 0 + iteration ceiling 10k | 55 | 0 | diagnostic combination rescues all three unstable pairs; all 52 existing valid calls and controls are identical, but all 14 directional failures remain |
| Frobenius/MU | 55 | 0 | rescues the three unstable pairs, but changes 16/52 existing valid calls and 59/69 controls |
| KL/MU | 56 | 0 | additionally rescues one CD79a pair, but changes 19/52 existing valid calls and 61/69 controls |
| rank 3, 2k | 32 | 25 | widespread nonconvergence and complete control churn |
| rank 3, 10k | 50 | 7 | five gains do not offset seven losses; 41/45 retained-valid calls change |
| rank 5, 10k (stopped after donor 6436) | 21 / 23 | 2 / 16 | rescues that donor's seven failures but makes B3TUBB and MPO unstable; 0/14 retained-valid calls and 0/23 controls remain identical |
| nested fit cap 50k / 100k / 200k | 53 / 53 / 52 | 0 / 0 / 0 | larger uniform samples do not repair any directional failure |
| reference-stratified sample | 50 | 6 | four gains with six losses; rejected |

The solver/objective alternatives therefore cannot be promoted merely because
they suppress the three unstable states: they materially repartition controls
and move already-valid calls. Increasing the iteration ceiling is a clean
numerical correction, while the fixed-sample diagnostic proves the INS failure
is sampling rather than optimizer variability; fixing the sample is not itself
a biological solution or permission to stop checking sample stability.
Rank 5 was stopped after its first complete donor because it had already
violated both preservation safeguards; completing the other donors could not
make it a global candidate.

A final targeted rank-5 diagnostic assigned every component to target- or
reference-favoured roles from its target/reference feature loading energy and
then aggregated each role. This repaired the arbitrary single-reference-
component swaps (median minimum control Jaccard 0.999898 across the 17 failures)
but produced only 4/17 valid pairs; the other 13 still failed the unchanged
directional gates. Multi-component role aggregation is therefore useful causal
evidence, not the missing global separator. Its pair, seed, and component-
loading tables are retained under
`data/pilots/nnmf_fit_components_5_energy_aggregation/`.

All tables and three-seed diagnostics are retained under
`data/pilots/nnmf_fit_*/`. They are review artifacts and are ignored by Git.
No production config was changed as a result of these experiments.

## Historical pair-web diagnosis and superseded six-target diagnostic

The later review established that the 17 failures were not a representative
"sparse-marker NNMF" sample. The frozen 69-row artifact uses E-cadherin for
42 arms, CD31 for 15, Pan-cytokeratin for 9, and EpCAM for 3. Thus 66/69 arms
use one of the three structural references flagged by image review. Every
non-valid arm uses E-cadherin (13) or Pan-cytokeratin (4). After fixing the
sample and increasing the iteration ceiling, the three numerical instabilities
resolve, while all 14 directional failures remain on those two references.
This supersedes the earlier assumption that fit design alone was the primary
unresolved variable.

The then-active review diagnostic was locked in
`phenocycler/fundamental_pair_diagnostics.json`. It deliberately covers only a
fundamental epithelial-versus-immune split:

- immune targets CD3e, CD79a, CD68, and MPO are each checked independently
  against EpCAM and Ker8_18;
- epithelial targets EpCAM and Ker8_18 are each checked independently against
  CD3e, CD79a, CD68, and MPO;
- all 16 arms are required as registered and are never pooled or selected after
  looking at results; and
- E-cadherin, CD31, Pan-cytokeratin, endocrine/subtype/support markers, B3TUBB,
  Vimentin, and SMA do not enter this diagnostic.

This is review-only. It does not silently rewrite the production registry or
claim endothelial, neural, mesenchymal, or endocrine assignments. Those
lineages remain unavailable rather than becoming negative or `Other` until a
morphology-safe identity is validated.

### Historical three-donor six-target survey

`scripts/review/fundamental_pair_survey.py` builds one target-independent global
top-20,000 reference mask per donor/reference, then evaluates the fixed 16-arm
roster with rank-2 Frobenius/CD NNMF, linear upper scaling, fixed fit sample 0,
three model seeds, and a 10,000-iteration ceiling. Component roles use squared
loading energy after row-wise unit normalization. Mean and median direction,
reference fold, target fold at least 2, convergence, and divisor stability all
fail closed. A target is available only when every independent arm passes;
the retained divisor and empirical tail probability are the conservative
maxima across arms.

Artifact: `data/pilots/fundamental_pair_survey_v2_20260812/`.

| Result | Count |
|---|---:|
| Technically valid arms | 18/48 |
| Target-fold failures | 18/48 |
| Median-direction failures | 12/48 |
| Converged seed fits | 144/144 |
| Donor-target aggregates passing every arm | 2/18 |

Only donor 6436 MPO and donor 6450 CD68 pass every required arm. No target is
valid across all three donors, and all six targets are unavailable in donor
6566. Every failed aggregate has no divisor and no calls. Divisor
reproducibility and reference-control Jaccard are exactly 1.0 throughout, so
this is a stable-but-wrong all-cell partition, not an optimizer failure. The
strict result is intentionally retained; no gate was weakened.

### Historical outcome-blind overlap audit

`scripts/review/reference_overlap_audit.py` measures exact q90/q95/q99/q99.5/
q99.9 co-high tails, reference-seed target burden, spatially standardized and
spatial-plus-reference-rank-standardized complements, and Cell/Nucleus/
Cytoplasm/Membrane patterns. It reads neither calls, pair states, phenotypes,
disease, annotations, nor desired frequencies. Its primary selector is global
top-20,000 reference rank; 8x8 round-robin is a separately labeled OFAT
comparator.

Artifact: `data/pilots/reference_overlap_fundamental_20260812/`. All 480
pair-tail rows, 384 compartment rows, and 36 seed audits are valid and complete.
The dominant residual is donor 6436 CD79a with EpCAM/Ker8_18: q99 overlap is
21.1-22.0% (21.1-22.0x independence) and q99.9 overlap is 13.5-14.3%
(134.5-143.4x). Spatial plus reference-rank adjustment attenuates but does not
remove it. Donor 6450 has a weaker, mostly composition-driven version; donor
6566 is essentially clean. The 6436 signal persists across nucleus, cytoplasm,
and membrane proxies, so it is not solely a membrane-enrichment effect. Scalar
means cannot distinguish a whole-object segmentation/channel artifact from
true same-cell signal; raw-image crops remain the adjudicator.

Overlap is an important veto, not a complete explanation of the separator
failures. No technically valid arm has q99 independence lift above 1, but many
failed arms are also strongly co-high-depleted. The reduced roster therefore
removes the original structural-reference confound while independently showing
that the rank-2 all-cell partition still fails on sparse/bulk geometry.

### Boundary/contact audit

`scripts/review/boundary_contact_audit.py` validates a scratch archive against
the immutable REDSEA donor snapshot before releasing contact attribution. It
supports direct donor receipts and follows adopted runs only through sealed,
content-verified adoption receipts and fingerprinted source snapshots. It
never searches for an arbitrary run with matching outputs.

The prioritized donor-6450 artifact is
`data/pilots/boundary_contact_fundamental_6450_v2_20260812/`. The scratch is
`STALE_SCRATCH`: 54,574/441,534 immutable cells (12.36%) disagree in cell area,
and 54,270 disagree in derived cytoplasm area. Consequently no contact, replay,
or selected-edge metrics were released. The independent raw lane produced 36
target/reference/DAPI crops for the three retired structural sentinels. A
visual spot-check shows selected cells abutting or containing bright neighboring
membrane/process structures, including a slender CD31-positive process crossing
or adjoining an E-cadherin field. This supports the reported mechanism but is
morphology evidence, not quantitative contact proof. Donor 6566 scratch is also
known to lack the required compartment arrays; a fresh fingerprinted scratch
generation is required before a full contact audit.

The dominant replacement-roster outlier was then reviewed directly at
`data/pilots/boundary_contact_cd79a_6436_v2_20260812/`: 24 deterministic crops
cover CD79a:EpCAM and CD79a:Ker8_18. None of the eight dual-bright panels is a
clean isolated whole-cell double-positive example. They instead show diffuse or
cluster-wide CD79a signal, epithelial membrane/process signal at shared borders,
and small/high-nuclear-fraction segmentations; one selected cell has no measured
cytoplasm. The 6436 scratch also fails the immutable area check, so quantitative
contact metrics remain withheld. This image evidence supports treating the
CD79a/epithelial overlap as technical until a fresh run-bound scratch archive
can adjudicate it.

## CD11b decision now enforced

CD11b is no longer:

- an active RESTORE target;
- a reference;
- a broad or subtype marker;
- a member of the `broad_immune` error family; or
- an anchor, support, exclusion, or required marker in typing rules.

The active pair web now has 23 targets. The Immune family has nine members, so
CD3e and the other Immune markers receive `0.01 / 9`, rather than spending an
allocation on CD11b. Registry and typing-rule validators explicitly reject any
future attempt to reintroduce CD11b into phenotyping.

## Code surfaces changed

Primary implementation:

- `phenocycler/restore_seeded.py` — new reference-seeded production evaluator;
- `phenocycler/restore_stage.py` — control-mask loading, estimation eligibility,
  empirical evidence/censor columns, fail-closed production policy;
- `phenocycler/restore_validation.py` — empirical-alpha control statistic and
  23-pair review rosters;
- `phenocycler/config.py` and `config.ini` — explicit locked candidate policy;
- `phenocycler/marker_registry.json`, `marker_registry.py`,
  `typing_rules.json`, and `hierarchical_typing.py` — CD11b exclusion;
- `phenocycler/pipeline.py` and `donor_pipeline.py` — controls artifact wired
  into calibrate;
- `scripts/review/preprocessing_survey.py` — reproducible one-factor survey;
- `scripts/review/divisor_survey.py` and `make_figures.py` — raw floors and
  corrected divisor distributions no longer share an axis.

Focused regression coverage is in `tests/test_restore_seeded.py` plus updated
RESTORE, registry, typing, config, pipeline, control, and floor tests.

## Verification

Full suite after the original implementation:

```text
779 passed, 18 skipped, 2 warnings in 143.96s
```

The total includes four new fit-sampling/provenance tests. The skips are
optional integration dependencies/data (geopandas, anndata, SimpleITK, and the
XeniumPanelExplorer submodule), not RESTORE failures. The two warnings are
Python multiprocessing `fork()` deprecation warnings.

Full suite after the reduced diagnostic, overlap audit, boundary/contact audit,
and verified adoption-lineage resolver:

```text
810 passed, 18 skipped, 2 warnings in 116.78s
```

Ruff is clean for all new diagnostic/audit scripts and tests, the focused
diagnostic set is 31/31, and `git diff --check` is clean.

After the 2026-08-13 mutable-development revision (subsequently reduced to two
comparisons):

```text
819 passed, 18 skipped, 2 warnings in 116.32s
```

The three revised diagnostic test files pass 37/37, targeted Ruff is clean,
and `git diff --check` is clean. The skips and warnings have the same optional
dependency and multiprocessing causes described above.

## Deliberately not completed after the stop request

1. No full marker-evidence/typing pilot was written for the new method.
2. No cohort-wide 20-donor comparison was run.
3. The structural pair web and the subsequent six-target/16-arm survey are both
   retired from active diagnostics. The new two-comparison mutable diagnostic
   is running, but it is troubleshooting evidence and was not promoted to
   production.
4. The acquisition exposure/responsivity sidecar for cross-donor MFI was not
   implemented. Such a scalar should not move within-donor RESTORE calls, but
   it is still needed before comparing intensity magnitudes across donors.
5. This original documentation gap is superseded: `README.md`,
   `DATA_README.md`, `WORKFLOW.md`, and `METHOD_HISTORY.md` now state the
   selected-input/empirical-alpha contract and mark the old p999/Otsu path as
   historical. Separate refinement documents were not rewritten as production
   authority.
6. Historical figures were not regenerated from the corrected plotting code.
7. No commit or PR was created.

## Safe resume sequence

1. Read this file and inspect the two preprocessing survey artifacts.
2. Preserve the dirty marker-debug work; do not reset the tree.
3. Treat the old 23-pair fit-design table as historical evidence, not the active
   diagnostic. Do not promote rank 3/5, MU, KL, reference-stratified sampling,
   or a fixed sample as a fallback cascade.
4. Use only the two upstream-input comparisons declared in roster schema v3.
   Keep CD79a out of the candidate roster; its historical raw-image/contact
   evidence remains useful only as an explanation for its exclusion.
   Development outputs are mutable and may be refreshed in place; do not treat
   them as production evidence. Require a materially improved three-donor
   result before any 20-donor run.
5. Keep endothelial, neural, mesenchymal, endocrine, and subtype evidence
   unavailable until morphology-safe anchors are separately reviewed.
6. Only after that, write full marker evidence and simultaneous assignments,
   compare them with `baseline_p999`, finish the authoritative docs, and seek
   production sign-off.
