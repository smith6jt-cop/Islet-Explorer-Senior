# Cell typing lane

Four steps plus an islet refinement. Step 1 decides which measurements the rest may lean on
and Step 2 puts them on a common scale — neither calls a cell type. Step 3 scores each cell
against the rule table by cosine, behind a requirement gate; Step 4 adds UTAG tissue domains
beside the per-cell calls; the islet lane uses the QuPath islet annotation as a prior to
improve endocrine assignment inside `Islet_N`, and lets every reference compete in the 20 µm
ring.

This lane derives no per-marker threshold, needs no target/reference pair, and does not
depend on the blocked `controls` / `calibrate` / `type` stages. Nothing here writes a
`StageManifest` (no current stage writes one) or mutates the stable
`data/{cells,geometry_qc,redsea,expression}` roots or the historical `data/runs/` tree.
`pipeline.py` computes no `run_id` at all — the content-addressed run framework it fed was
removed on 2026-09-16 — so nothing under `scripts/**` can re-key a run that no longer exists.
All output goes to the mutable repo-root `data/development/celltype/`.

```bash
cd core/phenocycler
python scripts/celltype/run_step1.py                  # whole cohort, ~10 min
python scripts/celltype/run_step1.py --donors 6450    # one donor, ~10 s
python scripts/celltype/figures.py                    # ~5 min
python scripts/celltype/run_step2.py --dry-run        # derive shifts only, ~1 min
python scripts/celltype/run_step2.py                  # + write normalised parquets
python scripts/celltype/verify_step2.py               # prove the shift is pure
python scripts/celltype/uniform_audit.py --help      # Step 8e: RESTORE residual-batch + marker-invariance audit
python -m pytest scripts/celltype/tests -q
```

Step 2 needs the UniFORM source, which carries no packaging metadata, so pip cannot pin it.
It is a gitignored checkout following this repo's existing convention (`.gitignore:74`), and
`uniform_adapter.py` asserts the commit at import so a drifted checkout stops the lane:

```bash
git clone https://github.com/smith6jt-cop/UniFORM.git external/UniFORM
git -C external/UniFORM checkout 7c33396d12e5a19dc0b36e46a94c3293de0cd814
pip install -e "core/phenocycler[uniform]"
```

## Which run this reads

`core/phenocycler/data/runs/30f6689e4d27d295f890/expression/` — the **only** run on disk, and
what `data/runs/LATEST` names.

An earlier version of this section argued the run choice was immaterial, because the two runs
then on disk hardlinked the same `expression` bytes. Both were deleted at the maintainer's
instruction (143 GB recovered), so there is nothing left to choose between, and the argument
describes nothing. What replaces it is the property that actually matters:

**This run is regenerable and the deleted one was not.** The registry used to ask for
`measurement_source = qupath_norm` on CD68, Iba1 and M2Gal3, which needs `QuPathNorm__<marker>`
columns that do not exist in `cells/` (0 present) — so re-running `expression` failed closed,
and the expression on disk was the only expression there was. Those three now take the `redsea`
default their `redsea_corrected` spillover policy already implied, so all 60 markers read one
source and `expression` regenerates.

**Measured: no value changes.** The expectation was that those three markers would move off
QuPath's signed scale onto the REDSEA whole-cell mean, so Step 1 and Step 2 would have to be
recomputed for them. They do not move. Re-running `expression` under the corrected registry
reproduced the old partition exactly **on donor 6380, the one donor compared** — same 62 numeric
columns, same `object_id` order, every column `allclose` (`panel.py:34-36` keeps that scope; the
pre-fix run has since been deleted, so the comparison cannot be widened).

The reason is that the requirement was never satisfied in the first place. `cells/` has never
carried a `QuPathNorm__` column, so the old run cannot have read one; the demand was added to
the registry *after* that run was built, and it is precisely what made the stage
un-regenerable. Removing it restores what the data always was rather than changing it. That is
also why `method_version` reads `selected-expression-v2-qupath-norm-inputs` while no number
moved — it is a provenance marker, not a value change.

Provenance: built 2026-08-28 from commit `0f5ea04` with a **dirty** worktree, so its source is
identified by `source_sha256 0542e887…` rather than by the commit alone. Its manifests record
a canonical `output.root`, so unlike the two deleted runs it was adoptable while the `adopt`
command still existed (it was removed with the rest of the content-addressed run framework on
2026-09-16; nothing here depends on that command running again).

`tests/` pins the content hash (`35fa90fa…`), the row count (31,454,041) and the agreement
between `panel.RUN_ID` and `LATEST`, so swapping the run or mutating the stage fails a test
rather than silently changing every number this lane produces. **Those two drifted apart once
already** — the pin outlived the run it named, and every output under
`data/development/celltype/` was unattributable until it was noticed. Repoint `RUN_ID`,
`EXPRESSION_CONTENT_SHA256` and the outputs together, never one without the others.

## Why these statistics and not the obvious ones

The measurement problem this lane exists for is that **raw intensities are not comparable
between donors** — donor-median fold-range is 4.2× across the panel, 22 of 55 markers span > 5×.
So every Step 1 statistic is computed *within* a donor and is invariant to any monotone
per-donor rescaling. That rules out means, variances, and anything compared to a fixed
number, and it rules in ranks, tails and spatial structure.

It also rules out the two things that look obvious and are not:

- **Whole-population correlation and PCA.** Background outnumbers signal by 100:1 or more,
  so those statistics describe the background. Average pairwise marker correlation is
  r̄ = 0.16–0.30, which forces PC1 to 21.8–36.3% of variance with 40–47 of 47 loadings
  sharing a sign. By the concentration statistic alone, CD20 (2.9% of signal in its
  brightest 1% of cells), CD79a (2.7%) and FOXP3 (1.9%) look as flat as DAPI (1.9%), which
  is on every nucleus. That reading is wrong: in the tail, CD20 & CD79a co-occur at **211×
  chance** and CD3e & FOXP3 at **41×**.
- **Decile contrasts as ground truth.** For a marker positive on ~0.1% of cells, the top
  decile is ≥ 99% not-that-cell. This is what made the DeepCellTypes verdict unattributable
  and it is why nothing here uses one.

## What Step 1 measures, per (donor, marker)

| family | statistic | reads as |
|---|---|---|
| shape | quantiles p10…p99.9, `zero_frac`, `top1_signal_share` | how much of the signal lives in the bright tail; DAPI ≈ 1.9% is the reference for "flat" |
| tail | `partner_best` vs `control_worst` — fold-enrichment over chance that a marker's brightest 0.1% are also another marker's brightest 0.1% | does this marker's bright tail land on the same cells as a marker of the same lineage, and *not* on the same cells as markers of other lineages |
| spatial | Moran's I on ranks and on the top-1% indicator, against a per-donor permutation null | is there any spatial structure at all — real staining is clustered because cells of a type sit together, noise is not |
| acquisition | between-group over worst within-group spread of donor medians, across calibration epoch × panel variant | does the level step with the acquisition group rather than with donors, i.e. a different antibody clone |

**Every tail enrichment is reported beside its negative controls.** A statistic with no null
is not a measurement, and this is the arm the DCT concordance gate never had. Controls are
*derived*, not hand-listed: for any marker they are the first two `CONTROL_REPRESENTATIVES`
that are neither the marker itself nor one of its declared partners. A hand-picked control
set can drift toward pairs that happen to behave; a rule cannot.

**A marker only has to show structure on one axis.** `morans_rank` and `morans_top1` answer
different questions, and a marker present on ~90% of cells has a scattered top 1% — measured:
E-cadherin `morans_top1` 0.009, `morans_rank` 0.148, and a 41× partner tail with Ker8-18.
That is the "flat because it is on everything" case, and condemning it would be wrong. So
`unusable` requires failing on *every* axis: no spatial structure above the null **and** no
partner tail.

## Verdict, and the two flags that are NOT part of it

The verdict answers exactly one question — **is this a valid measurement?** — and refuses to
answer two. It is fail-closed via `scripts/dct/gates.py` (imported, never copied), so absent,
non-finite or under-powered evidence is a FAILURE, never a pass.

- **unusable** — fewer than `MIN_CELLS` cells; more than half the cells reading exactly 0; or
  no spatial structure above this donor's own permutation null.
- **usable** — otherwise. **On this cohort that is all 1,120 (donor, marker) cells.**

Two further properties are **reported beside the verdict, never folded into it**, because
they constrain different later steps and neither makes a measurement invalid:

| column | means | constrains |
|---|---|---|
| `flag_acquisition_group` | the level steps with the acquisition group rather than with donors | Step 2 — the marker is fine but not *comparable* across those groups, and a rigid shift cannot repair a clone change |
| `flag_pair` | the declared partner tail does not clearly beat its non-partner controls | Step 3 — evidence about which cell types the panel can resolve |

An earlier draft folded `flag_pair` into the verdict and called 232 (donor, marker) cells
`degraded` for it. That was wrong, and the counter-examples are plain: CD11c is spatially
coherent (Moran's I 0.17), aligns cleanly under UniFORM (residual SD 0.042) and has a real
distribution — what fails is that **CD11c and CD209 do not mark the same cells**, which is a
fact about dendritic cells and this panel, not a defect in CD11c. Ki67 "failed" with a
partner tail of 276× because proliferating cells in this tissue are disproportionately
lymphocytes. CD79a "failed" with a partner tail of 170×.

## What changed after the first cohort run, and what did not

The first full run (2026-08-25) exposed three problems. **The acceptance bars were not
touched** — the failure that closed the DeepCellTypes lane was relaxing a pre-registered
criterion after seeing it fail, and that is not repeated here. What changed was declared
biology and one derivation rule, both checkable by someone who has never seen this cohort:

1. **Three declared partnerships were missing, so a partner was being used as a control.**
   `Granzyme_B ↔ CD3e` came out at **138× chance** as a "negative control" — Granzyme B is
   *in* CD3e+ cytotoxic T cells. Likewise `Keratin_5 ↔ Ker8-18` (12.2×) and
   `EpCAM ↔ Ker8-18` (8.9×) — all three are epithelial. Added to `LINEAGE_PARTNERS`.
2. **Controls could come from the marker's own compartment.** Cells of one compartment sit
   together and PhenoCycler spillover crosses touching cell boundaries, so a
   same-compartment marker is a poor negative control by construction. Measured on the first
   run: non-partner pairs where both markers are immune overlap at median **1.61×** (75th
   pct 5.82×) against median **0.00×** (75th 0.88×) for everything else. `control_map` now
   requires a different broad compartment, with the compartment inherited from a marker's
   partners where the registry declines to assign one (every state marker).
3. **`CD31 ↔ Podoplanin` stays a same-compartment control on purpose** — blood and lymphatic
   endothelium are different cells, and measuring them as *anti*-correlated (0.3×) is
   informative rather than a defect.
4. **The hormones were declared partners and are not.** `INS ↔ GCG` and `GCG ↔ SST` were
   listed with the note "islet endocrine, adjacent cells", which admits they fail the
   definition of a partner — *labelling the same cell*. Alpha, beta and delta are different
   cells, so their measured 1.2–1.4× is the right answer, and they are now declared
   controls. No hormone has a same-cell partner in this panel (no C-peptide), so they carry
   no `flag_pair`.
5. **The verdict was answering two questions and now answers one** (see above). The bars are
   unchanged; what moved is which column the answer is written to.

## Known marker-level findings this lane was built to record

- **CD57 is a reagent change, not a broken marker.** Its median steps with the *panel
  variant* and nothing else: 1.70–3.85 (2.3.1/CD38), 0.09–2.89 (2.4.0/CD38), 1.23–10.83
  (2.4.1/neither), **502–633 (2.4.1/b_Catenin1)** — same epoch, same exposure, same
  responsivity. Acquisition-group ratio 35.1; every other marker ≤ 1.67. It matters beyond
  CD57: **a rigid shift models a per-donor scale difference and does not model a different
  clone**, so a marker flagged here must not be normalised as though it were a scale
  problem, however good its shift statistic looks.
- **Eight markers are a different measurement, and not because REDSEA skipped them.**
  `DAPI, FOXP3, Ki67, PCNA, SOX2, TCF_1, TOX, TP63` are stored in `expression` as
  uncompensated **nucleus** means (`qupath:Nucleus__<m>`) rather than REDSEA-corrected
  whole-cell means. REDSEA *is* compartment-aware — it wrote `Nucleus__`, `Cytoplasm__`,
  `Membrane__` and `Cell__` for all 60 markers. What it did not do is **compensate** these:
  `marker_taxonomy.NO_SPILLOVER` lists them plus `Collagen_IV`, and `config.ini` sets
  `exclude_no_spillover = true`. Measured on donor 115, REDSEA's nucleus value differs from
  QuPath's by a median 2% — re-rasterisation, not correction — and the two carry the same
  tail evidence (partner medians 55.6 vs 53.6, controls 0.2 vs 0.2). So substituting it
  would change nothing measurable.
- **DAPI is usable as their denominator only nucleus-against-nucleus** — and on this rule
  table, that is nowhere. It removes the magnitude coupling (TOX ρ vs magnitude 0.632 →
  0.090) and destroys the identity signal (TOX & CD3e 115.7× → 0.0×), and it *manufactures*
  an association: the TOX & E-cadherin negative control goes 0.1× → **14.8×**. Regressing on
  `log(DAPI)` at its fitted slope (0.46–0.72) is better but still leaves that control at
  13.8×. DAPI carries nuclear size and chromatin density, which are cell-type traits, not
  nuisance. The one shape that survives is both sides carrying the same denominator —
  Ki67 & PCNA, 375.7× → 147.0×.
  So the ratio is decided **per reference**, not globally (`score.nuclear_only_references`):
  a reference qualifies only if every marker it constrains is nucleus-measured. **No
  reference in the current table qualifies**, because `rules.flatten` composes each phenotype
  with its ancestors, so even the proliferation split — which declares only Ki67 and PCNA —
  inherits Ker8-18/Pan-CK/E-cadherin and the immune negatives. The ratio is therefore inert
  here, by measurement rather than by switch, and turns on for a purely nuclear rule if one
  is ever added. `--no-nuclear-ratio` still disables it outright.
- **`CD11c & CD209` is the one declared lineage pair with no shared population** — 0.7× of
  chance in the top 0.1%, on every donor. That is a measurement about these two markers'
  bright tails, not a verdict on the panel: CD11c is shared with macrophages and CD209 is
  narrow, so a DC call resting on their co-occurrence has no support here. Whether dendritic
  cells are callable by some other combination is untested.

## What Step 2 produced

44 of 56 markers land at or under the 0.10 background-mode bar; median residual SD **0.048**.

| provenance | n | meaning |
|---|---|---|
| `automatic_registration` | 44 | shifted, under the bar |
| `automatic_registration_over_bar` | 12 | shifted and improved 1.2–6.0×, still over the bar, carried with the residual recorded |
| `refused_reagent_change` | 1 | CD57 — the shift is *not* applied, because a rigid shift cannot repair a clone change however good its statistic looks (716× → 2.79× would have looked fine) |

The twelve over the bar are `Vimentin 0.852, HLA_A 0.584, CD11b 0.454, Pan_Cytokeratin 0.405,
Beta_actin 0.339, E_cadherin 0.323, TCF_1 0.225, CD31 0.155, INS 0.145, PCNA 0.143,
CD44 0.118, Ker8_18 0.109`. Three attempted improvements were measured and all made things worse or no
better — see `uniform_adapter.py`, which records them so they are not retried.

**They fail for a reason no shift can fix: they differ between donors in distribution
*shape*, not in location.** After the shift each donor's histogram correlates with the
reference's at median 0.9704 for these against 0.9915 for the other 44 —
Spearman(residual SD, shape correlation) = **−0.607** over all shifted markers, Mann-Whitney
p = 1.6 × 10⁻⁵. Vimentin shows it plainly: bimodal, with the modes at different relative
positions in different donors, because it marks stroma, endothelium and immune cells in
donor-varying proportions. A rigid shift is a location model, and landmark refinement moves a
distribution rather than reshaping it, so neither applies.

The consequence for Step 3 is stated here rather than left implicit: **for these twelve
markers a single cohort-wide cutoff is not valid even after alignment.** Only within-donor
rank statistics are safe on them.

`verify_step2.py` proves the two claims the design rests on: the output is a pure
per-(donor, marker) multiplicative rescaling (max deviation from the recorded scale
2–9 × 10⁻¹⁶ over all 31.5M rows, `object_id` order identical), and therefore every Step 1
tail statistic is bit-identical before and after — CD20 & CD79a 169.9 → 169.9,
B3TUBB & CD56 483.0 → 483.0, and all four negative controls 0.0 → 0.0.

## What the islet lane measured, and the one donor it cannot measure

Cohort-wide over `islet_residents_all.parquet` (1,336,828 cells, 20 donors). Recorded here
because until 2026-08-27 these numbers existed only in commit messages.

**12.70% of islet-associated cells satisfy no rule** — `call = Ambiguous`,
`reason = no_reference_requirements_met`. All 169,734 of them are in `peri`; **zero** are in
`core`, which is what the endocrine prior is for. Per donor the peri share runs 0.159–0.339,
and it is **not** confounded with anything the cohort is unbalanced on: batch1 0.2547 vs
batch2 0.2494; ND 0.2551 / Aab+ 0.2614 / T1D 0.2480; acquisition epoch 2.3.1 0.2583 /
2.4.0 0.2525 / 2.4.1 0.2494.

Read `ambiguous_near_miss.parquet` carefully before calling these near misses. It reports
169,734 cells "one requirement short" for *eleven separate references* — the same 169,734
each time, i.e. every ambiguous cell. They are not close to a call; they are **below
`require_min` on every marker the table asks about**. That is a staining, segmentation or
panel-coverage question, not a rule-tuning one, and the leading candidate is acinar tissue,
which this 59-plex carries no enzyme marker for.

**Donor 6539's peri-islet ring was missing from the export, not from the tissue — and it has
been repaired.** The pre-repair export carried 541 islets with **124 peri cells** over 3 rings,
because `cell_region` is read off the QuPath `Parent` column (`cells_parquet.py`: `Islet_N` →
core, `Islet_N_exp20um` → peri) and 490,271 of 6539's cells were parented to `Root object
(Image)`. The rings themselves were drawn for essentially every islet at the right width (569
`IsletExpanded` polygons, median 20.08 µm — `scripts/kintsugi/README.md` § Step 0). The parent
assignment was re-resolved upstream and re-exported on 2026-08-27
(`Cellmeasurements_6539parentFixed.csv`, named in `cohort-spec.json`), and run
`30f6689e4d27d295f890` was built from it: 6539 now reads core 7,630 / peri 23,455 cells over 569
rings, and `islet_residents_all.parquet` carries 18,519 peri cells for it (the regeneration
section below). The former "excluded from peri-islet statistics" caveat is lifted by the data.

## What the 2026-08-28 regeneration found: the proliferating leaf is a PCNA background readout

The lane was re-run end to end against run `30f6689e4d27d295f890` — the first run of the
*current* code, since `run_islet.py` had not been re-run after the DAPI-ratio change. The
call distribution moved sharply, and the cause is measured, not inferred.

**What moved.** Over 1,351,845 islet-associated cells (up from 1,336,828 — the whole increase
is donor 6539's repaired ring, +18,395 peri and −3,378 core; **every other donor's cell count
is identical to the digit**):

| call | before | after |
|---|---:|---:|
| `Exo-acinar-ductal-proliferating` | 10,598 | **293,606** |
| `Epithelial` | 198,495 | 22,987 |
| `Epi-exocrine` | 37,385 | 0 |
| `Exo-acinar-ductal` | 33,996 | 0 |

**Why.** `Exo-acinar-ductal-proliferating`'s requirement gate is `any of [Ki67, PCNA]`, and
those two are the only thing separating it from `Exo-acinar-ductal` — the four epithelial
`anypos` markers are identical between them. So the leaf is decided by PCNA. Measured
fraction of islet-associated cells clearing `require_min = 0.10` on PCNA:

| donor | DAPI ratio inert (now) | DAPI ratio applied (before) |
|---|---:|---:|
| 115 | **22.26%** | 6.92% |
| 6374 | **21.84%** | 2.53% |
| 6476 | **25.52%** | 9.94% |

Ki67 barely moves either way (0.07–1.56%). 293,606 / 1,351,845 = **21.7%**, which is PCNA's
pass rate. Making the DAPI ratio inert — correct on its own terms, because the ratio
*manufactured* a TOX & E-cadherin association at 14.8× — raised PCNA's pass rate 3–9× and the
leaf absorbed the epithelial compartment.

**And PCNA cannot carry that call, by this lane's own Step 1 evidence.** Nothing here appeals
to how many cells "should" be proliferating; these are the pre-registered statistics, each
reported against its negative controls:

| | `top1_signal_share` | `control_worst` | `morans_top1` |
|---|---:|---:|---:|
| **PCNA** | **0.0295** | **15.1** (max 62.7) | 0.031 |
| Ki67 | 0.117 | 10.7 | 0.057 |
| Ker8_18 | 0.050 | 0.075 | 0.104 |
| CD20 | 0.030 | 4.2 | 0.062 |
| DAPI (the "flat" reference) | 0.0200 | 1.5 | 0.040 |

PCNA concentrates **2.95%** of its signal in its brightest 1% of cells — against DAPI's 1.9%,
and DAPI is on every nucleus. Its *negative controls* overlap its bright tail at a median
**15× chance**, two orders of magnitude worse than Ker8_18's 0.075×. Its top 1% is barely
spatially structured. On the concentration and control axes PCNA is the flattest non-structural
marker sampled, so a gate on it selects background, and a scale whose top is `p99.9` puts that
background above `require_min` — exactly the inherited property recorded in the Step 3 notes
above.

**Nothing was changed in response.** Not `require_min`, not the encoding, not the rule table.
Relaxing or tightening a pre-registered bar after seeing it fail is the failure that closed the
DeepCellTypes lane, and it is not repeated here. What follows is a question for the maintainer,
with three candidate answers that are all *declarations*, not tunings: drop `PCNA` from the
proliferation gate and let `Ki67` carry it alone; make the proliferation split a `state` score
beside the identity call rather than a competing leaf (decision #17's shape); or accept the
leaf and rename it for what it measures. **Until that is settled, no epithelial or exocrine
share from this lane means what its label says**, and the `Ambiguous` peri fraction — 24.38%,
170,405 of 698,974, still zero in `core` by construction — is the number to quote instead.

## Files

| file | role |
|---|---|
| `panel.py` | the cohort, its acquisition groups, the marker set, declared lineage partners and derived negative controls |
| `usability.py` | the Step 1 statistics and the fail-closed verdict |
| `run_step1.py` | Step 1 driver → `data/development/celltype/marker_donor_usability.parquet` |
| `figures.py` | one sheet per marker plus the cohort verdict sheet |
| `uniform_adapter.py` | the pinned UniFORM import, its histogramming, and the refuted-alternatives record |
| `run_step2.py` | Step 2 driver → `data/development/celltype/{shifts.json, uniform/donor_id=*/data.parquet}` |
| `verify_step2.py` | proves the shift is pure and rank statistics are invariant |
| `rules.py` | reads `phenotype_workflow.csv`, flattens the hierarchy onto each phenotype, and builds the reference matrix |
| `score.py` | the [0,1] scale, the cosine, the requirement gate, and the per-reference DAPI-ratio policy |
| `run_step3.py` | Step 3 driver → `data/development/celltype/calls/donor_id=*/data.parquet` |
| `utag_adapter.py`, `run_step4.py` | Step 4: per-donor UTAG domains → `data/development/celltype/domains/` |
| `run_islet.py` | the islet refinement → `data/development/celltype/islet_residents/` |
| `figures_ambiguous.py` | what the cells satisfying no rule actually are |
| `tests/` | the statistics on data whose answer is known, and the verdict failing closed |
