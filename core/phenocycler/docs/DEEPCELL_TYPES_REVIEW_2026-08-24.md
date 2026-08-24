# DeepCellTypes pilot — concordance gate and per-image review protocol

```yaml
status: pre-registered
written: 2026-08-24
supersedes: DEEPCELL_TYPES_REVIEW_2026-08-21.md
creates_cell_calls: true
creates_thresholds: false
production_outputs_written: false
owner: smith6jt
reviewer: TBD — record the name here before the first review field is opened
lane: mutable, data/development/. Writes no StageManifest and does not mutate data/runs/.
```

**Supersedes the 2026-08-21 protocol, which is left intact and unedited.** That document said a
change of mind must be a new dated protocol that says so; this is that document. No review field
was ever opened under it and its reviewer was never named, so nothing it governed has been spent.

## Why it is superseded

The 2026-08-21 protocol made a human per-image verdict the only release gate, and stratified the
reviewer's fields by the model's max-softmax confidence. On 2026-08-24 the feasibility check ran,
and two things it assumed turned out to be measurably false:

1. **Max-softmax does not predict agreement.** Under the 2026-06-15 checkpoint the *high*-confidence
   tercile has *worse* marker agreement than the low one on ten of the twenty-eight anchors below.
   A stratified review would have concentrated expert time on the stratum that is not better.
2. **The reviewed object was the wrong one.** The protocol reviewed the argmax `call`. On the ptft
   checkpoint the probability vector separates CD20-high from CD20-low cells at AUROC 0.951 while
   the argmax assigns `Bcell` to 0.0 % of them. A verdict on the calls alone cannot see that.

It also asked for adjudication of 51 specific classes, of which 16–20 are ever emitted and 25
cannot occur in this tissue or on this panel at all.

## What is being decided

Whether DeepCellTypes' calls on this cohort are right enough to build on, **at broad-compartment
granularity**. Islet location already comes from QuPath `Islet_N` geometry, so the model is never
required to call a beta cell; what it must deliver is the immune compartment and its separation
from the rest.

Nothing downstream happens inside this piece of work: no cohort run, no `apps/islet` change, no
promotion to a `vision_types` stage. A negative verdict is a real and expected outcome.

**Two gates, in order.** The quantitative concordance gate runs first and is cheap; the human
per-image review runs only on a configuration that passed it. Asking a person to adjudicate calls
whose probabilities move *against* this workflow's own measurements spends the scarcest resource in
the project on a question already answered.

---

## Gate 1 — Concordance with this workflow's own measurements

Implemented in `scripts/dct/concordance_gate.py`; thresholds in its `CRITERIA`. These are
acceptance bars on **rank agreement**, not marker thresholds. They derive no cutoff, need no
target/reference pair, and nothing in them may be tuned toward an expected cell-type frequency,
cross-donor uniformity or disease ordering.

### The anchors are derived, not chosen

Twenty-eight, produced by `derive_anchors()` and frozen in `scripts/dct/anchors.json`: every marker
in `phenocycler/marker_registry.json` with `kind == "type"` and a non-null `broad_compartment` that
resolves into the DeepCellTypes 278-marker vocabulary and is acquired on all 20 donors. Thirteen
Immune, nine Epithelial, three Endothelial, two Mesenchymal, one Neural. Caveolin and M2Gal3 drop
out with no vocabulary entry, recorded in the same file.

Hand-picking anchors after seeing a result would be the tuning this lane exists to avoid, so the
rule is code and the artifact is its output. `tests/test_concordance_gate.py` regenerates it and
fails if the two disagree.

### Both ends of the scale are reported

A correlation with neither end marked is unreadable, so every statistic is printed against:

- the **permutation null** — the same statistic on a probability matrix shuffled across cells, 200
  draws, 95th percentile. On donor 6539 that is |rho| ≤ 0.022 and |AUROC − 0.5| ≤ 0.030;
- the **measured ceiling** — `plane_identity_6539.txt`'s median self-rho of **0.920**, what
  agreement looks like when the two quantities really are the same thing measured twice, through
  this same `object_id` join.

### The label space is restricted first

Twenty-five of the model's 51 classes are removed before the argmax: 20 that cannot occur in human
pancreas and 5 with no supporting marker anywhere in this panel. The decisions and their written
reasons live in `scripts/dct/broad_class_crosswalk.json` (`in_label_space`, `exclusion_reason`), and
its `policy` block records that `expected_in_pancreas` was authored 2026-08-21 — before any
inference output existed — and is frozen.

Restriction is exactly a renormalization of the returned softmax over the kept columns, which is
algebraically identical to masking the excluded logits before the softmax. It is a re-read of the
model's own output, not a new decision, and `tests/test_label_space.py` asserts the identity.

Nothing is discarded: every cell keeps `p_impossible_mass`, `restricted_call` and
`unrestricted_call`, a cell whose removed mass exceeds its surviving best is `Ambiguous` with
reason `mass_on_impossible_classes`, and **every gate line and every figure reports the restricted
and unrestricted numbers side by side** — because the restriction moves some anchors up and others
down, and applying it silently would hide that.

Five exclusions carry a caveat rather than a certainty. Enterocyte, Goblet, Paneth, Foveolar and
ICC are excluded on the premise that the section holds no duodenum or stomach; nPOD pancreas-head
blocks can contain duodenum and this cohort already has a documented mixed-tissue case (donor 115,
pancreas + pLN). Each is verified per donor against the QuPath `Tissue` annotation and the answer
recorded, before the exclusion is trusted on that donor.

### The bar, stipulated in advance

Per donor, per checkpoint, at broad-compartment level, in the restricted label space. Counts are
fractions of the 28 anchors, because that set is derived and its size is a property of the panel.

| # | Criterion | Bar |
|---|---|---|
| 1 | **Sign** | Every anchor's rho against its own compartment > 0 |
| 2 | **Specificity** | ≥ 11/15 of anchors rank their compartment 1st of the compartments; none worse than 3rd |
| 3 | **Detectability** | AUROC ≥ 0.75 for every Immune anchor, ≥ 0.70 for the rest |
| 4 | **Argmax coherence** | ≥ 11/15 of anchors are *called* their own compartment more often in the marker's top decile than its bottom |
| 5 | **Impossible-class mass** | mean `p_impossible_mass` ≤ 0.10, and it beats the survivor for ≤ 10 % of cells |
| 6 | **Null separation** | Every anchor's AUROC lies outside the permutation null's 95th percentile |
| 7 | **Confidence informativeness** | ≥ 11/15 of anchors have top-tercile AUROC ≥ bottom-tercile. **Does not block release** — it decides whether Gate 2 may stratify at all |

Criteria 1 and 5 are the near-tautological ones. A negative sign means a compartment's probability
moves *against* the marker that defines it. A model spending a quarter of its mass on classes the
tissue cannot contain is not calling this tissue, whatever its anchor concordance says — and
criterion 5 uses no marker anchors, so it is genuinely orthogonal to 1–4. The rest are judgement,
and their one defensible property is that they were fixed before the confound factorial ran.

Absent, non-finite or under-powered evidence is a **failure**, never a pass — the rule
`scripts/dct/gates.py` already enforces, extended to this rectangular case by
`assess_anchor_concordance`. Below 50 usable observations an anchor is *unjudged*, and an unjudged
anchor fails.

**Criterion 4 replaced an invalid one, and the replacement is recorded rather than quietly made.**
The first draft capped each compartment's share of calls at 45 % and floored it at 0.5 %. Roughly
85 % of exocrine pancreas is acinar, so a *correct* model would have failed the cap; and requiring
≥ 0.5 % Neural calls is an assumption about how many nerve cells a pancreas holds. Both are the
expected-frequency reasoning `CLAUDE.md`'s first non-negotiable forbids. The replacement consults
no external frequency at all: if a compartment's probability separates a marker's top decile from
its bottom, the argmax must separate them the same way. The change was made before any verdict was
recorded, and it rescued nothing — every arm still fails.

### Aggregation, and what a pass is worth

Computed **per region**, pooled **per donor**, reported at both levels. A donor passes only if the
pooled statistic clears the bar *and* ≥ 6 of its 8 regions clear it individually, so one large
region cannot carry a donor. **Lane PASS requires all eight donors** — seven of eight is a fail with
a finding that names the donor, not a partial pass.

A pass is scoped to exactly what earned it: which compartments, which checkpoint, which window
size, which panel, which mask provenance, which donors. There is no aggregate pass across donors or
compartments and none may be computed later.

---

## The sampling frame — 8 donors × 8 regions

Fixed before any inference, from geometry and acquisition metadata only, all of which predate DCT.

The 2026-08-21 protocol drew 12 fields from each of 3 donors. Everything executed so far rests on
**one donor and one region**, and it is the wrong donor: peri-islet cells are **0.03 %** of donor
6539 — 167 cells on the whole slide — against a cohort median of 1.95 % and a maximum of 8.35 %. It
is a ~65× outlier at the bottom of the one quantity this study exists to measure. It was chosen for
cheapness and a validated archived mask, which were the right reasons for a first look and are the
wrong reasons for a verdict.

| donor | batch | panel | disease | mask | peri-islet % | why |
|---|---|---|---|---|---|---|
| 6539 | 1 | cd38-59ch | ND | archived | **0.03** | continuity with everything already run; the low extreme |
| 115 | 1 | cd38-59ch | ND | archived | **8.35** | the high extreme; 1,305 islets; the pancreas + pLN case |
| 6521 | 1 | cd38-59ch | Aab+ | archived | 1.95 | the only batch-1 Aab+; sits at the cohort median |
| 6593 | 1 | cd38-59ch | T1D | archived | 1.71 | batch-1 T1D |
| 6450 | 2 | b-catenin1-59ch | Aab+ | rasterize | 1.40 | smallest batch-2 slide (442 k cells) |
| 6534 | 2 | b-catenin1-59ch | T1D | rasterize | 3.91 | high-peri batch-2 |
| 6523 | 2 | neither-58ch | T1D | rasterize | 0.85 | the 58-channel panel, where every index from 16 shifts |
| 6591 | 2 | neither-58ch | T1D | rasterize | 2.36 | second 58-channel donor, different peri content |

**Disease coverage here is for tissue appearance, not contrast.** Batch 2 holds zero ND donors, so
every disease comparison available in this set is fully aliased with batch, acquisition software
and panel variant. This protocol answers technical validity only; no disease-ordered figure is
drawn from it, and none may be.

### Regions

Eight per donor, by a rule that cannot see the calls' plausibility:

1. Enumerate 1024-px cores holding ≥ 500 cells.
2. Order by `sha256(f"{donor}:{y0}:{x0}")` ascending; take the first **4**.
3. Take the first **4** in the same hash order among cores holding ≥ 1 cell whose QuPath parent is
   `Islet_N`. The study question is peri-islet; a blind draw over whole-pancreas area misses islets
   on the low-peri donors by construction.
4. All 64 origins are recorded before inference. A region that cannot be rendered is replaced by
   the next in hash order, and the substitution and its reason are recorded.

**Regions are anchor points, not tiles.** Each anchor gets concentric windows and only cells inside
the common 1024-px core are scored, so a window-size arm changes the normalization context and not
the cell population. This repairs the `--median-tile` arm below, which changed the tile *and* the
tissue and so could not attribute its own result.

### Masks

Thirteen donors have an archived `int32` mask; the seven batch-2 donors have none and a whole-slide
raster is ~2.4 GB. A windowed rasterizer streams each GeoJSON once for a bbox index, then draws
only the polygons intersecting a region, reusing `redsea.rasterize_mask`'s `label = feature_index+1`
invariant unchanged. It must be **bit-identical** to the archived mask over ≥ 20 windows on each of
the four archived-mask donors before it is trusted on the four that have none. A mismatch stops the
lane.

---

## Gate 2 — The per-image review, conditional and re-scoped

Runs only on a configuration that passed Gate 1.

- Fields are drawn by the same hash rule, over cells with a `restricted_call` — the 26-class space,
  never the `Ambiguous` ones.
- **6 fields × 8 donors = 48.** Eight donors at the old 12 fields each is 96 fields of expert
  adjudication, which will not happen; 48 fields over 26 classes is a denser per-class tally than
  36 over 51 was. If that is still too much, the reviewed **donor set** is cut and named — never
  the per-donor field count, so no donor is reviewed too thinly to judge.
- The reviewer adjudicates the **broad compartment**, with the specific class shown as context.
- **Stratification by max-softmax only if criterion 7 passed for that checkpoint.** Otherwise a
  single unstratified verdict, and the fact that confidence was uninformative is reported.

Carried forward unchanged from 2026-08-21, because they were sound: the four verdicts
(`Correct` / `Broad-only` / `Wrong` / `Cannot tell`), `Cannot tell` as a first-class answer excluded
from the denominator and reported separately, the per-class tally sheet with blank rows kept, the
partial-blinding statement, and the pass bar — `Wrong` ≤ 10 % of adjudicable cells, no single class
> 25 % of all `Wrong` calls, ≥ 50 % of reviewed cells adjudicable.

A **named reviewer before the first field opens.** `reviewer: TBD` is an unsatisfiable release gate.

---

## How to discard this cleanly

Unchanged from 2026-08-21 and still executable as one list. Nothing outside it is written, and in
particular nothing under `core/phenocycler/data/runs/`.

```text
data/development/deepcell_types_review/          # per-donor outputs, figures, review bundle
core/phenocycler/scripts/dct/                    # pilot code (git revert, do not rm -rf)
core/phenocycler/docs/DEEPCELL_TYPES_REVIEW_*.md # both protocols and this appendix
core/phenocycler/pyproject.toml                  # remove the [dct] optional-dependency group
~/.deepcell/models/                              # cached checkpoints, if no longer wanted
```

`git reset --hard` and `git clean -xfd` are **unsafe in this repository** — 157.9 GB of untracked,
ignored data lives under `core/phenocycler/data/`. Roll code back with `git revert`.

---

## Results appendix

Filled in below this line as results arrive. Everything here was produced before this protocol's
bar was written **only** in the sense that the inference ran first; the bar was fixed before the
confound factorial, which is what it governs.

### Executed 2026-08-24 — single-tile feasibility, donor 6539

Four arms, stock `predict()`, no halo, no hook. Densest 2048-px tile at origin **(16384, 18432)**,
10,794 of 10,800 cells joined to the run's expression parquet on `object_id`; the `--median-tile`
arm used a different tile with 7,401 cells. Run `3bf4fc8e391d007766d1`. Both checkpoints from
`~/.deepcell/models/`, `deepcell-types` at commit `59d74a9`.

Scored retrospectively by `concordance_gate.py`. Cells shown as **restricted (26 classes) / full (51)**.

| criterion | 2026-06-15 | 2026-06-15 +DAPI | 2026-06-15 median-tile | 2026-06-23-ptft |
|---|---|---|---|---|
| sign | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL |
| specificity | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL |
| detectability | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL |
| argmax_coherence | FAIL / FAIL | FAIL / FAIL | FAIL / PASS | PASS / PASS |
| impossible_mass | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL | PASS / PASS |
| null_separation | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL | FAIL / FAIL |
| confidence_informativeness | FAIL / FAIL | FAIL / FAIL | PASS / FAIL | PASS / PASS |
| **verdict** | **FAIL / FAIL** | **FAIL / FAIL** | **FAIL / FAIL** | **FAIL / FAIL** |

Per-arm detail, restricted label space:

**2026-06-15** (n=10,794) — null |rho| <= 0.022, |AUROC-0.5| <= 0.030
- `sign`: 18/28 positive; negative or unjudged: CD66, E_cadherin, EpCAM, GCG, INS, Ker8_18, Keratin_5, Pan_Cytokeratin, TP63, Vimentin
- `specificity`: 12/28 rank 1 (need 21); 11 worse than rank 3: CD34, CD66, E_cadherin, EpCAM, GCG, INS, Ker8_18, Keratin_5, Pan_Cytokeratin, TP63, Vimentin
- `detectability`: 13/28 clear their bar; below: CD34=0.648, Podoplanin=0.562, CD66=0.146, E_cadherin=0.283, EpCAM=0.322, GCG=0.179, INS=0.333, Ker8_18=0.404, Keratin_5=0.327, Pan_Cytokeratin=0.433, TP63=0.485, CD206=0.746, CD209=0.734, Iba1=0.654, Vimentin=0.474
- `argmax_coherence`: 18/28 anchors are called their own compartment more often in the marker's top decile than its bottom (need 21); not: Podoplanin(0.6% vs 5.0%), CD66(1.2% vs 16.0%), E_cadherin(0.4% vs 14.8%), EpCAM(1.1% vs 14.0%), GCG(0.4% vs 20.5%), INS(2.7% vs 15.9%), Ker8_18(0.8% vs 12.7%), Keratin_5(4.4% vs 14.9%), Pan_Cytokeratin(0.8% vs 12.2%), B3TUBB(0.0% vs 0.0%)
- `impossible_mass`: mean 0.240 (bar 0.1); outweighed 34.1% (bar 10%)
- `null_separation`: null |rho|<=0.022, |AUROC-0.5|<=0.030 at q0.95; inside the null: TP63, Vimentin
- `confidence_informativeness`: high >= low on 18/28 judged anchors (need 21); does NOT block release -- it licenses confidence stratification

**2026-06-15 +DAPI** (n=10,794) — null |rho| <= 0.022, |AUROC-0.5| <= 0.030
- `sign`: 18/28 positive; negative or unjudged: CD66, E_cadherin, EpCAM, GCG, INS, Ker8_18, Keratin_5, Pan_Cytokeratin, TP63, Vimentin
- `specificity`: 11/28 rank 1 (need 21); 11 worse than rank 3: CD34, CD66, E_cadherin, EpCAM, GCG, INS, Ker8_18, Keratin_5, Pan_Cytokeratin, TP63, Vimentin
- `detectability`: 12/28 clear their bar; below: CD34=0.646, Podoplanin=0.562, CD66=0.146, E_cadherin=0.281, EpCAM=0.325, GCG=0.180, INS=0.335, Ker8_18=0.399, Keratin_5=0.333, Pan_Cytokeratin=0.428, TP63=0.489, CD206=0.748, CD209=0.733, Iba1=0.653, MPO=0.750, Vimentin=0.483
- `argmax_coherence`: 18/28 anchors are called their own compartment more often in the marker's top decile than its bottom (need 21); not: Podoplanin(0.6% vs 4.7%), CD66(1.5% vs 16.7%), E_cadherin(0.4% vs 15.7%), EpCAM(1.7% vs 14.9%), GCG(0.4% vs 21.2%), INS(3.1% vs 16.8%), Ker8_18(1.1% vs 14.1%), Keratin_5(5.1% vs 15.7%), Pan_Cytokeratin(1.1% vs 13.6%), B3TUBB(0.0% vs 0.0%)
- `impossible_mass`: mean 0.243 (bar 0.1); outweighed 34.4% (bar 10%)
- `null_separation`: null |rho|<=0.022, |AUROC-0.5|<=0.030 at q0.95; inside the null: TP63, Vimentin
- `confidence_informativeness`: high >= low on 17/28 judged anchors (need 21); does NOT block release -- it licenses confidence stratification

**2026-06-15 median-tile** (n=7,401) — null |rho| <= 0.025, |AUROC-0.5| <= 0.034
- `sign`: 18/28 positive; negative or unjudged: Podoplanin, CD66, E_cadherin, EpCAM, GCG, INS, Ker8_18, Keratin_5, Pan_Cytokeratin, TP63
- `specificity`: 15/28 rank 1 (need 21); 10 worse than rank 3: CD34, CD66, E_cadherin, EpCAM, GCG, INS, Ker8_18, Keratin_5, Pan_Cytokeratin, TP63
- `detectability`: 14/28 clear their bar; below: CD34=0.606, Podoplanin=0.486, CD66=0.189, E_cadherin=0.293, EpCAM=0.428, GCG=0.292, INS=0.278, Ker8_18=0.381, Keratin_5=0.421, Pan_Cytokeratin=0.421, TP63=0.500, CD68=0.408, Vimentin=0.640, B3TUBB=0.627
- `argmax_coherence`: 17/28 anchors are called their own compartment more often in the marker's top decile than its bottom (need 21); not: Podoplanin(0.9% vs 4.0%), CD66(2.5% vs 13.6%), E_cadherin(0.0% vs 8.5%), EpCAM(0.4% vs 3.8%), GCG(0.0% vs 12.3%), INS(0.9% vs 11.9%), Ker8_18(0.0% vs 5.3%), Keratin_5(4.7% vs 5.5%), Pan_Cytokeratin(0.0% vs 4.9%), TP63(3.0% vs 3.0%), B3TUBB(0.0% vs 0.0%)
- `impossible_mass`: mean 0.253 (bar 0.1); outweighed 34.9% (bar 10%)
- `null_separation`: null |rho|<=0.025, |AUROC-0.5|<=0.034 at q0.95; inside the null: Podoplanin, TP63
- `confidence_informativeness`: high >= low on 21/28 judged anchors (need 21); does NOT block release -- it licenses confidence stratification

**2026-06-23-ptft** (n=10,794) — null |rho| <= 0.022, |AUROC-0.5| <= 0.027
- `sign`: 23/28 positive; negative or unjudged: Keratin_5, TP63, CD206, CD4, B3TUBB
- `specificity`: 16/28 rank 1 (need 21); 3 worse than rank 3: CD34, Keratin_5, B3TUBB
- `detectability`: 12/28 clear their bar; below: CD34=0.615, Podoplanin=0.613, CD66=0.491, INS=0.684, Keratin_5=0.304, TP63=0.624, CD11c=0.591, CD20=0.465, CD206=0.448, CD4=0.350, CD79a=0.538, CD8=0.702, FOXP3=0.732, MPO=0.610, SMA=0.582, B3TUBB=0.476
- `argmax_coherence`: 21/28 anchors are called their own compartment more often in the marker's top decile than its bottom (need 21); not: Podoplanin(1.9% vs 4.7%), CD66(12.4% vs 23.6%), Keratin_5(11.8% vs 35.6%), TP63(13.3% vs 16.3%), CD206(77.8% vs 81.4%), CD4(68.2% vs 81.8%), B3TUBB(0.0% vs 0.0%)
- `impossible_mass`: mean 0.087 (bar 0.1); outweighed 5.4% (bar 10%)
- `null_separation`: null |rho|<=0.022, |AUROC-0.5|<=0.027 at q0.95; inside the null: CD66, B3TUBB
- `confidence_informativeness`: high >= low on 24/28 judged anchors (need 21); does NOT block release -- it licenses confidence stratification

### What these four arms settle

**Every arm fails Gate 1, in both label spaces.** That is the honest state, and it is why the
confound factorial exists rather than a verdict.

- **The DAPI question is closed.** `DAPI → dsDNA` in versus out shifts every anchor by ≤ 0.013.
  Withholding the panel's only nuclear channel is not what is wrong. The marker-dropout arm did its
  job and needs no repeat.
- **Tile density is not a rescue.** The median-density tile makes the 2026-06-15 checkpoint *worse*
  on the anchors, not better. But that arm changed the tile *and* the tissue, so it cannot attribute
  its own result — which is why regions are anchor points with concentric windows from here on.
- **The two checkpoints fail differently, and the figure's naive reading inverts.** 2026-06-15 is
  the better *immune* separator (CD3e AUROC 0.891 against ptft's 0.830) and collapses on epithelium;
  ptft is the better *epithelial* separator (E-cadherin 0.926) and fails MPO, CD20, SMA and B3TUBB.
  Choosing between them by eye would have picked ptft for the wrong reason.
- **The argmax discards what the probability knows.** Under 2026-06-15, E-cadherin's *top* decile is
  called Epithelial **0.4 %** of the time and its *bottom* decile **14.8 %** — inverted. Under ptft,
  CD20-high cells separate at AUROC 0.951 while `Bcell` is assigned to 0.0 % of the tile.
- **Impossible-class mass is the cleanest single discriminator**, and it needs no anchors:
  **0.240** for 2026-06-15 against **0.087** for ptft, with the removed mass beating the survivor
  for **34.1 %** versus **5.4 %** of cells. Removing those classes does not repair 2026-06-15 —
  its `Tumor` calls go 51 % to Macrophage and 35 % to Neutrophil, sloshing into the largest immune
  sink rather than returning to the right class. ptft's `Enterocyte` goes 39 % to Epithelial, which
  is the correct redirection.
- **Confidence is not informative for 2026-06-15** (criterion 7 FAIL) and is for ptft (PASS). Under
  the superseded protocol the review would have been stratified for both.

Figure: `data/development/deepcell_types_review/6539/feasibility_gate_panels.png` — all four arms,
28 anchors, drawn against the permutation null and the 0.920 ceiling.

### Not yet run

The confound factorial and its confirmation set. Until they run, **no verdict about DeepCellTypes
is recorded here**: the panel handed to the model is 52 channels against a training median of 34,
20 of them activation/exhaustion state markers, with no Amylase, CD45, CHGA, CDX2, GP2 or MUC2 — so
acinar cells have no positive identifier and gut epithelium has no negative one. A failure under
that configuration is a statement about the configuration.

---
