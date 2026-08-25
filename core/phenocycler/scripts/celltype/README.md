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
python scripts/celltype/run_step1.py                  # whole cohort, ~12 min
python scripts/celltype/run_step1.py --donors 6450    # one donor, ~10 s
python -m pytest scripts/celltype/tests -q
```

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

## Verdicts

`usable` / `degraded` / `unusable`, each carrying its reason and the per-check detail in
`checks`. The verdict is fail-closed via `scripts/dct/gates.py` — imported, never copied —
so absent, non-finite or under-powered evidence is a FAILURE, never a pass.

- **unusable** — fewer than `MIN_CELLS` usable cells; more than half the cells reading
  exactly 0; or no evidence of signal on any axis.
- **degraded** — the level steps with the acquisition group; or the bright tail overlaps a
  non-partner more than its declared partner; or the partner tail does not clearly beat its
  controls.
- **usable** — otherwise.

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

## Files

| file | role |
|---|---|
| `panel.py` | the cohort, its acquisition groups, the marker set, declared lineage partners and derived negative controls |
| `usability.py` | the statistics and the fail-closed verdict |
| `run_step1.py` | Step 1 driver → `data/development/celltype/marker_donor_usability.parquet` |
| `tests/` | the statistics on data whose answer is known, and the verdict failing closed |
