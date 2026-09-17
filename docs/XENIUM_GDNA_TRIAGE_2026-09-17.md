# Xenium gDNA triage — the cohort, measured

**Increment 8 of the layer plan.** Supersedes "run the gate for the remaining 18 donors on the
joint arm", which was withdrawn on 2026-09-16 because its detection predicate spuriously
anchors immune lineages (see `docs/XENIUM_PIPELINE.md` and the plan's rewritten increment 8).

`metrics_summary.csv` was read by no code — trap 2 of `XENIUM_PIPELINE.md` — so the cohort's
genomic-DNA background had never been looked at. All 46 bundles were pulled from
`/orange/maigan/Xenium_T1D_PancPLN/` on 2026-09-16 (kilobytes each) and read here.

Reference: **0.00454**, the `adjusted_genomic_control_probe_rate` of sample 0041323, which is
what every detection constant in this repo was tuned against.

## Finding 1 — the cohort is not in the flagged pair's regime

| | rate | × tuned |
|---|---:|---:|
| cohort maximum (0041139 pLN, donor 6580) | 0.05522 | 12.2× |
| **0080870 panc, donor 6533** (same-slide pair) | **0.26100** | **57.5×** |
| **0081366 panc, donor 112** (same-slide pair) | **0.30600** | **67.4×** |

No cohort bundle reaches a fifth of the flagged pair. The gDNA correction specified for
0080870/0081366 is **not** automatically required cohort-wide, and nothing here justifies
applying it blind.

## Finding 2 — pLN sections carry more gDNA than pancreas sections, systematically

This was not expected and is the reason to read the file at all.

```
tissue   n   median rate     min        max
panc    23     0.002115   0.000598   0.013349
pln     23     0.010480   0.002888   0.055217
```

**18 serials carry both tissues. The pLN section is the noisier one in 17 of them**, median
**3.9×**, up to **39×** (serial 0041139: pLN 0.05522 against panc 0.00142).

The pairing controls for run, panel, reagent lot and operator — the two sections went through
the same Xenium run — so this is a property of the tissue or its handling, not of the batch.
It tracks median transcripts per cell inversely, which is what background does when real
signal falls: the three worst bundles have medians of 69, 49 and 217 transcripts against a
cohort panc median well above 400.

**Consequence for the plan.** The pLN arm of the study question — increments 9, 10 and 11, and
card S3's third compartment — rests on the noisiest Xenium data in the cohort. Any lineage
called in pLN on a raw `count >= 1` predicate inherits that. Per `XENIUM_PIPELINE.md`, the fix
is a change of scale rather than a higher floor, and the pLN sections are where it matters most.

## Finding 3 — three bundles above 10× the tuned rate, all pLN

| serial | tissue | donor | rate | × tuned | median tx |
|---|---|---:|---:|---:|---:|
| 0041139 | pLN | 6580 | 0.05522 | 12.2× | 69 |
| 0037637 | pLN | 6457 | 0.05475 | 12.1× | 49 |
| 0041140 | pLN | 6347 | 0.05094 | 11.2× | 217 |

Donor 6457 is centrally excluded (`core/phenocycler/phenocycler/cohort.py`), so two of the
three are live. Neither has been corrected.

## What this does and does not license

- It **does** partition the cohort: everything except those three sits under 5× the tuned rate.
- It **does not** license typing any donor. The detection predicate is still `raw count >= 1`
  in six places (`analysis/xenium/scripts/hierarchical_phenotype.py`), which anchors a
  129-gene Immune_T core by gDNA alone at P = 0.834 per cell against Exocrine's 0.016. That
  bias is independent of these rates and is unfixed.
- It **does not** replace the Y-probe validation. `XENIUM_PIPELINE.md` requires any predicate
  to put zero-template Y probes at or below ~0.3%; that check is unrun for every donor here.

## Reproducing

Per-bundle values are in `data/development/xenium_gdna/gdna_triage_2026-09-17.parquet`
(46 rows: serial, tissue, donor, rate, median_tx, n_cells). Sources are
`/data/xenium_pln/_metrics/*.metrics_summary.csv`, one per bundle, named for the bundle they
came from. Donor is resolved from the filename token where present and from
`core/phenocycler/data/integration/xenium_paths.csv` otherwise; all 46 resolved.
