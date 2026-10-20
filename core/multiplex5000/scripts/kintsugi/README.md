# KINTSUGI analysis lane

Runs the PhenoCycler cohort through the analysis end of
[KINTSUGI](https://github.com/smith6jt-cop/KINTSUGI) — the maintainer's own published
pipeline (Zenodo `10.5281/zenodo.14984518`; Smith et al., *STAR Protocols* 6, 103976, 2025) —
as an **independent second method** beside `scripts/celltype/`, not as a replacement for it.
Both lanes answer against the same rule table, so where they disagree is measurable.

Only the analysis end is in scope. KINTSUGI's image processing — illumination correction,
stitching, `KDecon`, EDoF, registration, autofluorescence removal, InstanSeg segmentation —
is upstream of this repo and is not run here.

This lane derives no per-marker threshold from a frequency, needs no target/reference pair,
and does not depend on the blocked `controls` / `calibrate` / `type` stages. Nothing here
writes a `StageManifest` (no current stage writes one) or mutates the stable
`data/{cells,geometry_qc,redsea,expression}` roots or the historical `data/runs/` tree.
`pipeline.py` computes no `run_id` at all — the content-addressed run framework it fed was
removed on 2026-09-16 — so nothing under `scripts/**` can re-key a run that no longer exists.
All output goes to the mutable repo-root `data/development/kintsugi/`.

```bash
cd core/multiplex5000
python scripts/kintsugi/run_step0.py --smoke     # prove the SQL, write nothing
python scripts/kintsugi/run_step0.py             # both exports, ~3 min
python scripts/kintsugi/run_step0.py --validate-only
python -m pytest scripts/kintsugi/tests -q
```

## Which KINTSUGI, and which half of it

KINTSUGI's tree carries two analysis lineages and only one is KINTSUGI's own:

| lineage | where | method |
|---|---|---|
| **this lane ports it** | `notebooks/4_Segmentation_Analysis.ipynb` + `notebooks/Kanalysis.py` | feature block → **pyFlowSOM SOM** → scanpy → Leiden → merged clusters |
| separate `panc_CODEX` work | `CODEX_Panc_scimap_Analysis.ipynb`, `6551_Analysis.ipynb` | scimap `napariGater` → `rescale` → `phenotype_cells` |

Pinned to `external/KINTSUGI @ 716d1eb` (v1.2.3), a gitignored clone following the convention
`uniform_adapter.py` already uses. `som_adapter.py` asserts the commit at import, so a drifted
checkout stops the lane. The pin is on the *commit*, not the distribution version, because the
analysis end lives in `notebooks/` and is not part of the installed package.

**pyFlowSOM is not a declared KINTSUGI dependency.** It is a lazy import inside
`Kanalysis.visualize_som` and appears in no extra, in `requirements.txt`, in `constraints.txt`
or in any `envs/*.yml` at the pinned commit — so `kintsugi install analysis` does not supply
it and it must be installed explicitly.

## Environment

KINTSUGI documents its own, and this lane uses it rather than inventing one:

```bash
conda env create -f external/KINTSUGI/envs/env-linux.yml
conda activate KINTSUGI
kintsugi install gpu                 # first: provides CUDA torch + CuPy
kintsugi install analysis viz
pip install -c external/KINTSUGI/constraints.txt pyFlowSOM
kintsugi-check
```

Each step writes Parquet, so steps run in whichever environment they need.

**What actually runs this lane today.** The documented `conda env create` above was attempted
on this machine on 2026-08-27 and **did not complete** — it terminated during
`Solving environment` without writing an environment, after ~35 minutes. Nothing was
damaged: the pre-existing `KINTSUGI` conda env is untouched (and is itself dangling, editable
against `/home/smith6jt/KINTSUGI/src`, which no longer exists).

So the SOM steps run on a minimal virtualenv that honours KINTSUGI's own `constraints.txt`:

```bash
python3 -m venv external/.venv-som
external/.venv-som/bin/pip install -c external/KINTSUGI/constraints.txt \
    pyFlowSOM pandas scipy scikit-learn pyarrow shapely
```

That gives `numpy 1.26.4` — `pyFlowSOM` requires `numpy <2`, which is exactly why KINTSUGI
pins it, and why the repo `.venv` (numpy 2.5.1) cannot run the SOM. Steps 0, 4 and 5 run in
the repo `.venv`. The documented conda path stays the reference for a full KINTSUGI install
including the image-processing half; this lane needs only the numeric constraint.

## Step 0 — the measurement the run does not carry

`expression` carries `cell_region` (`core`/`peri`/`tissue`/`other`) and `islet_num`, both from
the QuPath `Parent` string. `cell_region == 'peri'` is the ring annotation, so it answers
"inside the ring or not" and cannot answer "how far". KINTSUGI's peri-islet analysis is built
on a continuous signed distance and its zone is `0 < d <= 50 µm`, which the ring cannot express.

That measurement already exists and was simply not selected into `cells/`: the cohort exports
carry `Signed distance to annotation Islet µm` and `Signed distance to annotation Tissue µm`
beside `Object ID`. Step 0 is therefore a **column-selective extraction, not a recomputation** —
1,228 columns in `Cellmeasurements.csv` and 992 in `CellmeasurementsBatch2.csv`, of which four
are read. 212 GB scanned in 178 s.

`data/Islets/*.geojson` could be used to recompute the distance, and deliberately is not: a
recomputation is a *second* measurement with its own boundary and rounding conventions.
QuPath's is the one the annotations were drawn in, and the one KINTSUGI's 6551 analysis used.
The polygons are used for something the distance column cannot tell you — see below.

### What the validation measured

`cell_region` is derived at ingest from a different column of the same export, so checking the
distance against it is a real cross-check rather than a restatement. Over all 20 donors,
31,454,041 cells:

| check | result |
|---|---|
| `object_id` resolving `expression` → distance | **0 unmatched, every donor** |
| `core` cells reading a positive distance | 0.33 %–2.72 % per donor |
| `peri` cells reading ≤ 0 | 0.65 %–2.08 %, uniform |
| `peri` cells beyond their donor's own rings | **0.0000 every donor** (6482: 0.0001) |

The two non-zero rates are geometry, not error, and they are reported separately because they
mean different things. A `core` cell reads positive when a cell straddling the islet boundary
has its **centroid** outside it. A `peri` cell reads ≤ 0 because the QuPath column is the
distance to the *nearest* islet annotation while `Parent` names one specific islet — so where
islets sit close together, a cell in one islet's ring can be inside another. That it is uniform
across every donor is what makes it geometric rather than a per-donor defect.

### Donor 6482's rings are 25 µm, and nothing in the data says so

The first validation pass put 6482's peri violations at **17.7 %**, ten times every other
donor. Relaxing the bar would have buried it. Measuring instead — from the `Islet` and
`IsletExpanded` polygon pairs in `data/Islets/*.geojson`, converted with each image's own
`pixel_size_um_x` from the QuPath manifest — gives the cause:

There is **one ring per islet** — 489 of them in 6482, 1,304 in donor 115 — so "the ring
width" is a median over rings, never a single number a donor possesses. Reported that way:
**19 of 20 run donors draw their rings at a median 20.1 µm; donor 6482 draws its at a median
25.1 µm.** The spread matters too and travels with it (`p5`, `p95`, `n_rings_measured`):
rings whose expansions merge with a neighbour's read narrower, which is why `p5` runs to 0 on
several donors while `p95` sits at the nominal. The annotations are still *named*
`Islet_N_exp20um`, so `cell_region` reads `peri` either way and `expression` cannot show the
difference. Against 6482's own measured median the violation rate falls to 0.0001.

Two consequences, and only the second is a repair:

* **Any statistic taken over `cell_region == 'peri'` inherits this**, including the existing
  cosine lane's islet refinement. A 25 µm annulus has ~25 % more area than a 20 µm one around
  the same islet, so 6482 contributes more cells at a lower density for a geometric reason.
  Recorded here so the number is not read as biology.
* **The continuous distance does not inherit it.** A zone defined as `0 < d <= 50 µm` is the
  same zone in every donor however the annotation was drawn. This is the strongest argument
  for Step 0: it does not merely add resolution, it removes a cohort defect the categorical
  label cannot express.

`ring_width.measure()` re-derives the median from the polygons on demand and returns the
ring counts beside it, and the test suite pins both the 25.1 µm and the 20.1 µm so a re-run
cannot quietly change either.

### Donor 6539 — the rings existed; the cells were not assigned to them (now fixed)

The inherited account, in `celltype/README.md` and commit `e8b40d8`, was that 6539's
"expansion had been run on 3 of 539" islets. **Measured, that is not what happened.** Its
GeoJSON carries **569 `IsletExpanded` polygons against 572 `Islet` polygons**, at a correct
**20.08 µm** median (`p5` 20.03, `p95` 20.09). The rings were drawn, for essentially every
islet, at the right width. What was missing was the **cell-to-ring parent assignment**: in the
run's `cells/`, 490,271 of 6539's 501,991 cells carried `Parent = Root object (Image)` and
only 167 cells fell in an `Islet_N_exp20um`, over 3 rings.

**Repaired upstream on 2026-08-27** and re-exported as
`/home/smith6jt/IO60panc2nd/results/Cellmeasurements_6539parentFixed.csv` (1,237 columns,
503,962 rows). The repair does exactly what the diagnosis predicted:

| | before | after |
|---|---|---|
| `Annotation (Tissue)` | 0 | 472,058 |
| `Root object (Image)` | 490,271 | 806 |
| `Islet_N` | 11,553 | 7,634 |
| `Islet_N_exp20um` | 167 | 23,464 |
| **rings carrying cells** | **3** | **569** |

569 matches the polygon count exactly.

**Object identity survived**: all 501,991 `object_id`s in the run's `cells/` are present in the
new export, none are missing, and the 1,971 extra rows are the ones ingest filtered. So the
join key is intact and nothing needs re-keying.

**And the distances did not move.** Over all 503,962 shared cells, `Signed distance to
annotation Islet µm` is bit-identical before and after (`max |diff| = 0.0`), and the tissue
distance is non-null and equal in both. The geometry was always right; only `Parent` was
wrong. **Step 0's extracted distances therefore stand as-is, 6539 included, and do not need
re-running.**

### A claim this corrected

An earlier version of this lane recorded 6539's `signed_dist_tissue_um` as *unavailable*,
citing `panel_availability.json`'s `has_tissue_annotation: false`. That field is inferred from
whether `Annotation (Tissue)` ever appears in the `Parent` column — a statement about
hierarchy resolution, not about whether the annotation exists. 6539's export carries a
non-null tissue distance for all 503,962 cells, identical before and after the repair, and
QuPath could not have computed it without the annotation. `NO_TISSUE_ANNOTATION` is now empty
and 6539's tissue distance is treated as the measurement it is.

**`cell_region` was re-ingested on 2026-08-28.** Run `30f6689e4d27d295f890` was built from the
fixed export (`cohort-spec.json` names `Cellmeasurements_6539parentFixed.csv`; `cells_parquet.py`
filters each CSV to the images it is contracted to), so 6539 reads 569 rings and 23,455 peri
cells in `cells/`, `expression/` and `celltype/uniform/`. Nothing in this lane patches around it.

## Step 1 — the feature block, and the scaling question

Markers are **derived, never hand-listed** — from the registry, `panel.py` and Step 2's own
`shifts.json`, so this lane and `celltype/README.md` cannot drift apart. Dropped: `IAPP`
(registry `failed_marker`), `DAPI` (registry `structural`, on every nucleus), `CD38` /
`b_Catenin1` (panel variants; a cohort-wide block has no honest value for a donor that never
acquired them), and `CD57` (`shifts.json`: `refused_reagent_change`, acquisition-group ratio
35.1 — a rigid shift models gain, not a different clone).

**CD11b stays**, matching `scripts/celltype/panel.py`: it was set aside for the retired RESTORE
lane, and the rule table's Myeloid gate is `CD68 | CD11c | CD11b`. Dropping it would remove a
third of that gate's evidence — and would mean the two lanes were not given the same evidence,
which would make the concordance measurement meaningless.

Arms: **A** 55 markers (over-bar markers as-is), **B** 55 (over-bar as within-donor rank),
**C** 43 (over-bar excluded). The 12 over-bar markers are read from `shifts.json`, not typed
out: `Beta_actin, CD11b, CD31, CD44, E_cadherin, HLA_A, INS, Ker8_18, PCNA, Pan_Cytokeratin,
TCF_1, Vimentin`. Six of them are the epithelial/endocrine/endothelial backbone of the rule
table, which is why arm C is not simply the safe default.

### `StandardScaler` is feature weighting, not donor correction

`Kanalysis.visualize_som` calls `StandardScaler().fit_transform` on the pooled block, and
`docs/INTENSITY_SCALE.md` permanently forbids "z-scoring across donors". These are not the same
operation: UniFORM has already done the cross-donor alignment, and the scaler that follows
equalises **marker** variance so a bright marker does not dominate Euclidean distance.

That distinction is checked rather than asserted. A per-marker affine map moves every donor's
median identically, so the spread of donor medians in units of the marker's pooled SD is
invariant; donor correction shrinks it. `features.scaler_is_feature_weighting` measures it and
**fails the arm** rather than warning. `tests/test_features.py` pins both halves: `StandardScaler`
passes at `< 1e-6`, and an actual per-donor centring is refused.

## Step 3 — two findings the pilot forced

### The floor is the whole gate

`REQUIRE_MIN = 0.10` is a position on the `[background mode, robust top]` scale, and an early
draft of `annotate.py` took the floor as `min()` instead of the background mode. `min()` is
the REDSEA clamp at zero, not the background, and the consequence is not subtle: **56.8 % of
all (cell, marker) values then sit above the gate**, every cluster shows 28–43 of 55 markers,
every reference's negative constraints are violated, and 94 % of cells come back `Ambiguous`.
With the aligned background mode — mirroring `celltype/run_step3.build_scale`, which is what
the cosine lane's own `REQUIRE_MIN` is measured from — that falls to **14.7 %**.

A known and inherited property is pinned beside it: the scale's top is `p99.9`, so for a
marker carrying *no* positive population the top is the background's own upper tail, the span
collapses, and a third of pure background clears the gate. `celltype/score.marker_scale` has
the same form. This is what Step 1's usability verdict exists to catch, and the test records
it so that changing the scale has to be deliberate.

### Nested references are not competing candidates

`rules.flatten` composes each phenotype with its ancestors, so a cluster satisfying
`Mac-M2-like` necessarily satisfies `Imm-macrophage` and `Immune` too. Scored against all 37
references at once, the matched null answered `cores_cleared_but_are_not_separable` for 10 of
25 clusters — **correctly**, since those cores share most of their markers by construction —
and its argmax-z awarded lineages the gate had refused (`Endo-delta` to an exocrine cluster,
`Neural` to a fibroblast one).

The hierarchy resolves that, not a statistic. `annotate.resolve` takes the maximal candidates;
if one, that is the call; if several siblings the evidence cannot separate, the call is their
lowest common ancestor within the candidate set — four T subsets under `Imm-T-cell` still
support "this is a T cell". Candidates spanning different top-level branches have no common
ancestor below `all` and stay `Ambiguous` with that reason. The null is then asked about the
*resolved* reference alone, which is the question that is actually open — is this cluster
enriched for **this** reference? — rather than a competition the hierarchy has already settled.

Measured on the 6551 pilot at 25 nodes, this took named clusters from 4 to 15 and `Ambiguous`
cells from 93.4 % to 45.8 %, and made the run faster (34.5 s → 9.5 s).

## Pilot result — donor 6551, arm A, 25 nodes

1,379,107 QC-eligible cells of 1,845,970 (466,862 QC-ineligible, 1 non-finite), 55 markers.

| call | share of cells |
|---|---|
| `Ambiguous` | 45.83 % |
| `Exo-acinar-ductal` | 16.14 % |
| `AMBIGUOUS-CD11b-pos-epithelial` | 10.65 % |
| `Immune` | 10.49 % |
| `Mes-CD34-fibroblast` | 6.57 % |
| `Endo-vascular` | 5.89 % |
| `Mac-M2-like` | 2.75 % |
| `Imm-T-cell` | 0.75 % |
| `Mes-smooth-muscle-pericyte` | 0.60 % |
| `Mes-fibroblast` | 0.31 % |

Six of the ten unnamed clusters fail with `no_reference_requirements_met` and they are among
the largest, which is the same signature the cosine lane records for its 169,734 no-rule cells:
pancreas is ~85 % acinar and this 59-plex carries no acinar enzyme marker. Consistent with
that lane rather than a separate defect — but **not yet demonstrated to be the same cells**;
that is what Step 4's concordance is for.

25 nodes is KINTSUGI's own value for a 26-marker lymph node and is coarse for 55 markers
against 37 references. The node count is chosen by the pre-registered sweep, never to make a
frequency come out.

## Step 4 — what the two lanes actually disagree about

Joined on `(donor_id, object_id)` over donor 6551's **48,326 islet-associated cells** (core +
peri, from `celltype/islet_residents_all.parquet` — the cosine lane's per-cell `calls/`
currently covers only donor 6450). At 25 nodes:

| | |
|---|---|
| cells called by both | 15,019 |
| exact-label agreement | **27.6 %** |
| top-level-branch agreement | **62.5 %** |
| refused by this lane | 31,135 |
| refused by the cosine lane | 6,808 |
| refused by both | 4,636 (Jaccard **0.139**) |

Two things fall straight out, and neither was visible from either lane alone.

**The refusals are largely not the same cells.** 68.1 % of the cosine lane's refusals are also
refused here, but only 14.9 % of this lane's are also refused there. Both lanes attribute
their `Ambiguous` fraction to the same cause — pancreas is ~85 % acinar and this 59-plex
carries no acinar enzyme marker — and agreeing on the *explanation* is not evidence of
refusing the *same cells*. At this resolution it mostly is not, so the acinar account covers
less than it appears to from one lane.

**The single largest disagreement is endocrine.** 18,668 cells — 39 % of the joined set — are
`Endo-alpha` to the cosine lane and `Ambiguous` here. That is what a 25-node SOM over 1.38 M
whole-tissue cells should be expected to do: islet cells are a few percent of the tissue, so
alpha/beta/delta cannot occupy their own nodes and the blended node median satisfies no
endocrine rule. It localises the disagreement to a resolution question rather than a
methodological one, and it is precisely what the node sweep is for.

The comparison is scoped to islet-associated cells of one donor and is a measurement of two
methods against each other — not a release gate for either. Neither lane has ground truth.

## The node sweep — reported in full, including what it does not settle

KINTSUGI's notebook uses a 5x5 grid (25 nodes) for a 26-marker lymph node. The count is swept
here rather than adopted, and **never chosen to make a cell-type frequency come out**. Donor
6551, arm A, 55 markers, 1,379,107 QC-eligible cells; every node was occupied at every grid.

| nodes | SOM | clusters named | label agr. | branch agr. | refused by this lane |
|---|---|---|---|---|---|
| 9 | 9.8 s | 4 / 9 | 0.291 | **0.846** | 38,520 |
| 16 | 13.2 s | 11 / 16 | 0.352 | 0.671 | 34,737 |
| 25 | 17.4 s | 15 / 25 | 0.276 | 0.625 | 31,135 |
| 36 | 22.0 s | 19 / 36 | 0.347 | 0.687 | 35,863 |
| 49 | 29.0 s | 25 / 49 | **0.382** | 0.756 | 35,580 |
| 64 | 39.5 s | — | — | — | — |
| 100 | 54.5 s | — | — | — | — |

Agreement with the cosine lane rises with resolution from 25 nodes up (label 0.276 → 0.382,
branch 0.625 → 0.756), and the fraction of nodes that earn a name holds near half throughout.
n9's branch agreement of 0.846 is not the best result on this table — it names 4 clusters and
refuses the most cells, so it is agreeing about very little.

**This sweep does not yet choose the node count**, and saying so is the point. The
pre-registered criterion is assignment stability across seeds, which needs repeat runs at each
grid that have not been done; agreement with the cosine lane is a *reported* measurement, not
the selection rule — choosing the grid that best matches the other lane would make the
concordance circular. Rows for 64 and 100 are the SOM timings only; their annotation was still
running when this was written.

## Deviations from the notebook, declared

* **PCA order.** Notebook cell 97 runs `sc.pp.pca` and cell 99 runs `normalize_total` / `log1p`
  / `scale` *after* it, so `sc.pp.neighbors` at cell 101 reuses an `X_pca` computed from
  unnormalised values and the neighbour graph, PAGA, UMAP and Leiden all descend from it.
  Nothing downstream reads a normalised PCA. That is an ordering mistake, not a design choice,
  and reproducing it would put a known defect under every number here. This lane runs
  normalise → log1p → scale → PCA → neighbours.
* **Cluster count.** Notebook cell 113 hardcodes `n_final_clusters = 23`. The mechanism
  (correlation distance, Ward linkage, `fcluster` at `maxclust`) is kept; the constant is not.
  A cluster count chosen to land near an expected number of cell types is the frequency
  reasoning `CLAUDE.md`'s first non-negotiable forbids. It comes from a pre-registered
  stability sweep instead, and the whole sweep is reported.
* **Progressive loop.** `visualize_som` refits the SOM from scratch for every `rlen_iter` to
  render a frame per iteration; only the last is the result and it equals a single fit at
  `rlen = num_passes`. `som_cluster` fits once.

## Files

| file | role |
|---|---|
| `distance.py` | the QuPath distance columns, the extraction SQL, and the validation |
| `ring_width.py` | per-donor ring widths (median over rings) from the annotation polygons |
| `run_step0.py` | Step 0 driver → `distance/`, `distance_validation.json` |
| `features.py` | marker derivation, the three arms, the scaler check, donor purity, block loading |
| `som_adapter.py` | the pinned KINTSUGI import, headless SOM, meta-clustering, provenance |
| `run_step2.py` | Step 2 driver → `clusters/` (single run or a node sweep) |
| `annotate.py` | the background-mode scale, the gate, hierarchy resolution, the matched null |
| `run_step3.py` | Step 3 driver → `calls/` |
| `concordance.py`, `run_step4.py` | the join against the cosine lane → `concordance/` |
| `tests/` | 28 tests: the claims above, pinned (22 collect in the repo `.venv`, where the pyFlowSOM test is skipped; all 28 in `external/.venv-som`) |

Reused unchanged: `scripts/celltype/panel.py` (cohort, markers, run pin),
`scripts/celltype/shifts.json` (marker classes), `data/qupath_manifest.json` (paths, pixel size).

## Status

Steps −1, 0, 2, 3 and 4 run end to end on the donor 6551 pilot. `python -m pytest
scripts/kintsugi/tests -q` — 28 passed in `external/.venv-som` (22 passed, 6 skipped in the repo `.venv`).

Not yet done: the three-arm scaling bake-off at cohort scale (arm A only so far, and on one
donor there is no cross-donor question for it to answer), the cohort run, and **Step 5** — the
three-level donor → islet → cell hierarchical bootstrap vendored from KINTSUGI's
`peri_islet_hierarchical_bootstrap.py`, which is what Step 0's continuous distance was
extracted for.

Two caveats travel with every downstream output: any ND → Aab+ → T1D reading is **aliased with
acquisition batch** (batch 2 is every high-immune donor and no ND donor), and **donor 6482's peri
rings are 25 µm, not 20 µm**, so any statistic over `cell_region == 'peri'` is comparing unequal
annuli. The former 6539 exclusion is lifted: the current run carries its repaired ring (§ Step 0).
