# Cell typing — preprocessing lane

Two steps, in order. Step 1 decides which measurements the rest may lean on; Step 2 puts
them on a common scale. Neither calls a cell type.

This lane derives no per-marker threshold, needs no target/reference pair, and does not
depend on the blocked `controls` / `calibrate` / `type` stages. Nothing here writes a
`StageManifest` or mutates `data/runs/`; `pipeline.py` hashes only `phenocycler/` and
`pyproject.toml` into `run_id`, so nothing under `scripts/**` can re-key a run. All output
goes to the mutable repo-root `data/development/celltype/`.

```bash
cd core/phenocycler
python scripts/celltype/run_step1.py                  # whole cohort, ~10 min
python scripts/celltype/run_step1.py --donors 6450    # one donor, ~10 s
python scripts/celltype/figures.py                    # ~5 min
python scripts/celltype/run_step2.py --dry-run        # derive shifts only, ~1 min
python scripts/celltype/run_step2.py                  # + write normalised parquets
python scripts/celltype/verify_step2.py               # prove the shift is pure
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

## Which run this reads, and why it does not matter

`core/phenocycler/data/runs/3bf4fc8e391d007766d1/expression/`. Three reasons, in order of
how much they matter:

1. **The choice is immaterial.** Both runs on disk carry the *same* expression stage — not
   equivalent, identical: same inode, `content_sha256 86e148ed…`, 31,454,041 rows, same
   `schema_sha256`. They diverge only at `controls` / `calibrate` / `type`, and this lane
   reads none of those. `METHOD_HISTORY.md`'s judgement that `c8412757d12de36ca056` is
   "evidence of a failed phenotyping implementation" is about its *assignments*; its
   expression is the same bytes.
2. `data/runs/LATEST` names it.
3. **It cannot be regenerated.** The current registry asks for `measurement_source =
   qupath_norm` on CD68, Iba1 and M2Gal3, which needs `QuPathNorm__<marker>` columns that do
   not exist in `cells/` (0 present) — so re-running `expression` today fails closed. The
   expression on disk is the only expression there is.

Provenance caveat worth carrying: it was built 2026-08-05 from a **dirty** worktree at commit
`ab9a578`, so its source is identified by `source_sha256 0c7489d2…` but is not reconstructible
from a commit alone. `tests/` pins the content hash and row count, so swapping the run or
mutating the stage fails a test rather than silently changing every number this lane produces.

## Why these statistics and not the obvious ones

The measurement problem this lane exists for is that **raw intensities are not comparable
between donors** — donor-median fold-range is 4.2× across 55 markers, 22 of 55 span > 5×.
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
- **Eight markers are a different measurement.** `DAPI, FOXP3, Ki67, PCNA, SOX2, TCF_1,
  TOX, TP63` are stored in `expression` as raw uncompensated QuPath **nucleus** means
  (`qupath:Nucleus__<m>`), not REDSEA-corrected whole-cell means.
- **DAPI is not usable as their denominator** — tested and rejected. It removes the
  magnitude coupling (TOX ρ vs magnitude 0.632 → 0.090) and destroys the identity signal
  (TOX & CD3e 115.7× → 0.0×), and it *manufactures* an association: the TOX & E-cadherin
  negative control goes 0.1× → **14.8×**. Regressing on `log(DAPI)` at its fitted slope
  (0.46–0.72) is better but still leaves that control at 13.8×. DAPI carries nuclear size
  and chromatin density, which are cell-type traits, not nuisance.
- **`CD11c & CD209` is the one declared lineage pair with no shared population** — 0.7× of
  chance in the top 0.1%, on every donor. Dendritic cells are not supported by this panel.

## What Step 2 produced

44 of 55 markers land at or under the 0.10 background-mode bar; median residual SD **0.048**.

| provenance | n | meaning |
|---|---|---|
| `automatic_registration` | 44 | shifted, under the bar |
| `automatic_registration_over_bar` | 11 | shifted and improved 1.2–6.0×, still over the bar, carried with the residual recorded |
| `refused_reagent_change` | 1 | CD57 — the shift is *not* applied, because a rigid shift cannot repair a clone change however good its statistic looks (716× → 2.79× would have looked fine) |

The eleven over the bar are `Vimentin 0.852, HLA_A 0.584, Pan_Cytokeratin 0.405, Beta_actin
0.339, E_cadherin 0.323, TCF_1 0.225, CD31 0.155, INS 0.145, PCNA 0.143, CD44 0.118,
Ker8_18 0.109`. Three attempted improvements were measured and all made things worse or no
better — see `uniform_adapter.py`, which records them so they are not retried.

**They fail for a reason no shift can fix: they differ between donors in distribution
*shape*, not in location.** After the shift each donor's histogram correlates with the
reference's at median 0.9704 for these eleven against 0.9915 for the other 44 —
Spearman(residual SD, shape correlation) = **−0.607** over all 55 markers, Mann-Whitney
p = 1.6 × 10⁻⁵. Vimentin shows it plainly: bimodal, with the modes at different relative
positions in different donors, because it marks stroma, endothelium and immune cells in
donor-varying proportions. A rigid shift is a location model, and landmark refinement moves a
distribution rather than reshaping it, so neither applies.

The consequence for Step 3 is stated here rather than left implicit: **for these eleven
markers a single cohort-wide cutoff is not valid even after alignment.** Only within-donor
rank statistics are safe on them.

`verify_step2.py` proves the two claims the design rests on: the output is a pure
per-(donor, marker) multiplicative rescaling (max deviation from the recorded scale
2–9 × 10⁻¹⁶ over all 31.5M rows, `object_id` order identical), and therefore every Step 1
tail statistic is bit-identical before and after — CD20 & CD79a 169.9 → 169.9,
B3TUBB & CD56 483.0 → 483.0, and all four negative controls 0.0 → 0.0.

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
| `tests/` | the statistics on data whose answer is known, and the verdict failing closed |
