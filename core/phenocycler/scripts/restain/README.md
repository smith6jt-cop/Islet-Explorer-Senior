# Post-Xenium IF restain lane

Registers the post-Xenium IF to its own Xenium section **on tissue outlines alone**, then
checks the result against a signal the registration never saw.

The post-Xenium IF images the *same physical section* Xenium measured. So the Xenium
segmentation **is** the segmentation — no QuPath run on the re-stains, no new cell objects,
and the correspondence is exact by construction rather than inferred. That is the maintainer's
route and it is what makes this lane possible at all.

```bash
cd core/phenocycler
../../.venv/bin/python scripts/restain/run_register.py              # every export found
../../.venv/bin/python scripts/restain/run_register.py --donors 6442
../../.venv/bin/python -m pytest -q scripts/restain/tests           # 6 passed
```

Inputs are the maintainer's QuPath exports under `~/IO60panc2nd/PostXeniumIHC/geojson/`
(two pixel classifiers, `Tissue.json` and `Islets.json`) and the Xenium bundle's
`cells.parquet`. Output goes to `data/development/restain/`. Nothing here writes a
`StageManifest` or touches `data/runs/`, and nothing under `scripts/**` is hashed into
`run_id`.

## No identity is used, on either side

The IF supplies a pixel-classifier tissue polygon; Xenium supplies a cell footprint rasterised
from centroids. Neither is a cell type, so this is not gated on `require_current_pheno` and it
could run on all 20 donors the moment their exports exist.

It also makes the islet check genuinely independent: **the registration never saw an islet and
never saw a gene.**

## Measured, 2026-08-28, on the two exports available

| donor | IF serial | Dice before → after | islets kept | scale fitted |
|---|---|---|---|---|
| 6436 | 0059907 | 0.753 → **0.890** | 125 of 4,671 | 0.960 |
| 6442 | 0059947 | 0.922 → **0.951** | 95 of 3,366 | 0.991 |

Both clear `qc_tissue_dice_min = 0.80`, in ~6 s each. The fitted transforms are near-identity
with sub-degree rotation and translations of ~1,000 µm and ~250 µm — a slide-placement offset,
which is what the same section imaged twice should look like. Non-rigid refinement is skipped
and says so: `SimpleITK` is not installed, so `register_image(nonrigid=False)` rather than an
affine result labelled as a B-spline.

**The 0.25 µm/px scale is confirmed, not assumed.** Donor 6442's IF tissue spans
9,190 × 15,020 µm against 9,436 × 15,100 µm of Xenium cell footprint — 0.5–2.6%. A wrong pixel
size would show up as a factor, and `geo.assert_scale_is_plausible` refuses one rather than
letting the registration fail obscurely.

## The `Islet` class needs a size filter, and it is not optional

The classifier emits thousands of sub-cell fragments beside the real islets:

| donor | objects | median area | ≥50 µm across | share of islet area |
|---|---:|---:|---:|---:|
| 6436 | 4,671 | 27 µm² (5.9 µm) | **125** | 84.6% |
| 6442 | 3,366 | 23 µm² (5.4 µm) | **95** | 81.1% |

The median object is one cell wide. Filtered, the counts land in the range the PhenoCycler
cohort shows (67–1,304 islets per donor) and give far more landmarks than a point-set
registration needs.

## The independent check: IF islets against Xenium endocrine RNA

Cells within 25 µm of a kept IF islet centroid, scored on the 18-gene
`01_pancreas_endocrine` identity core, against a placement null (the same number of arbitrary
on-tissue neighbourhoods, 20 draws):

| donor | inside | elsewhere | ratio | null (p95) |
|---|---:|---:|---:|---:|
| 6442 | 2.0340 | 0.1176 | **17.30×** | 0.116 (0.161) |
| 6436 | 0.3577 | 0.1265 | **2.83×** | 0.117 (0.167) |

Both clear the null. This is the smallest real cross-modality measurement the project can
make: no cell typing on either side, no verified marker pair, no hormone threshold. Islets
found by a protein classifier land, after a tissue-outline-only alignment, on cells whose
endocrine RNA reads 17× (6442) and 2.8× (6436) above their surroundings.

### Why 6436 is weaker is a geometry question, not yet a disease one

Split by the core's own `subtype_hint`:

| donor | alpha | beta | pan | delta |
|---|---:|---:|---:|---:|
| 6442 (ND) | 6.33× | **25.38×** | 19.79× | 0.56× |
| 6436 (T1D) | 2.83× | **2.89×** | 3.12× | 0.86× |

6442 shows beta dominant over alpha, which is what an intact islet should look like. **6436
shows alpha and beta enriched by the same factor** — and beta-cell loss depresses beta
*relative to* alpha. A uniform scale-down of every endocrine gene set is instead what a mask
landing partly off-target does, and 6436 is also the donor with the worse Dice (0.890 vs
0.951) and a 4% tissue-extent mismatch. So the current reading is registration quality, and it
is checkable: improve 6436's alignment and see whether beta separates from alpha.

**n = 1 per group. Nothing here may be ordered by disease state**, and this is a registration
diagnostic, not a result about T1D.

`HHEX` (delta) reads *below* background in both donors and is not interpretable — a one-gene
panel, which is the thin-panel trap decision #12 exists to prevent. It is reported rather than
hidden, and it should not be quoted.

## The hormone check: it does NOT deliver cell-level correspondence

With the channel map recorded (`Opal 690 = INS`, `Opal 520 = GCG`, `Opal 570` FAILED,
`Sample AF` disregarded), the same-cell question is now runnable:

```bash
../../.venv/bin/python scripts/restain/run_hormone.py
```

Neither side uses a threshold. The protein side is an AUROC. The RNA side's label is
`beta_surrogate > alpha_surrogate` **within the same cell** — a contrast, not a gate — over
cells already inside an islet, which is what makes it a test of *subtype* rather than of
whether the registration found islets.

| donor | informative cells | INS ranks beta-like | GCG ranks alpha-like | null |
|---|---:|---:|---:|---|
| 6436 | 6,034 (2,807 β / 3,227 α) | **0.533** | **0.517** | 0.500 [0.484, 0.513] |
| 6442 | 6,212 (5,194 β / 1,018 α) | **0.574** | **0.562** | 0.500 [0.479, 0.521] |

Both clear a label-permutation null, so the correspondence is real rather than chance. **But
0.53–0.57 is far below the 0.75 this was pre-registered at, and the bar is not moved.** On the
evidence below, the cell-level claim is not supported and should not be made.

### Where the signal actually lives

Spearman of the RNA contrast (β−α) against the protein contrast (INS−GCG), widening the
neighbourhood:

| radius | 0 (cell) | 10 µm | 20 µm | 30 µm | 50 µm | 80 µm | 120 µm |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6442 | +0.074 | +0.103 | +0.154 | +0.197 | +0.231 | +0.266 | **+0.287** |
| 6436 | +0.025 | +0.036 | +0.050 | +0.054 | +0.053 | +0.051 | +0.045 |

This is the diagnostic that says how to read the AUROCs, and it rules out the comfortable
explanation. A fixed registration **offset** would make agreement climb and then *plateau*
once the radius exceeded the offset. 6442 instead climbs monotonically and is **still climbing
at 120 µm, which is wider than an islet** — so what agrees is islet-to-islet, not cell-to-cell.
6436 is flat at every scale. Neither is an alignment artefact that better registration would
fix.

### Why each donor falls short — different reasons

| | 6436 | 6442 |
|---|---|---|
| islet-resident cells with ≥2 β genes detected | **5.7%** | 36.0% |
| …with ≥2 α genes | 1.7% | 2.2% |
| "informative" cells decided by a single transcript | **49.2%** | 26.4% |
| protein INS~GCG ρ (should be negative) | **−0.387** | **+0.085** |
| `Opal 570` (FAILED) ~ `Opal 690` (INS) ρ | +0.290 | **+0.712** |

**6436's ceiling is the RNA side.** Its protein is clean — INS and GCG anticorrelate at
−0.387, and GCG ranks INS-low cells at AUROC 0.687 — but only 5.7% of its islet cells detect
two beta surrogates, and half its "informative" labels rest on one transcript.

**6442's ceiling is the protein side.** Its RNA label is much deeper, but its two hormone
channels do *not* anticorrelate, and its INS channel shares a **0.712** rank correlation with
the channel recorded FAILED. Whatever is in a failed channel is background, so that is a
caution about what `Opal 690` is measuring in that section.

### Autofluorescence is reported, not applied

The instruction is to disregard `Sample AF`, so the headline numbers are raw. But it is not
inert and that is recorded rather than argued: `Sample AF ~ Opal 690` is +0.50 (6442) / +0.44
(6436), and regressing it out moves 6442's INS~GCG from **+0.085 to −0.006**. It does *not*
rescue the AUROC (0.574 → 0.570), so AF is not the cap.

`Sample AF` is **clipped at 255 in 18.1% (6442) and 28.8% (6436)** of islet-resident cells, so
a linear AF correction is not even well defined there — the covariate is censored exactly
where it is largest. Every AUROC ships with these numbers beside it.

## What this does not yet do

**Cell-level protein↔RNA typing.** The measurement above runs, clears its null, and does not
reach the bar; the scale profile says the agreement is structure-level. So this lane currently
supports islet-level cross-modality claims and **not** per-cell ones. What would change that
is on the two ceilings above — a deeper RNA endocrine label, and a hormone channel not sharing
0.7 with a failed one — not a better transform.

**pLN.** The exports are **pancreas only** (maintainer, 2026-08-29), which is what settles the
ROI column: each of serials 0059865/0059907/0059947 carries both a `__Panc__` and a `__pLN__`
Xenium bundle, so the ROI was never derivable from the filename. A lymph-node export would
need its own row in `restain_sections.csv`.

**The other three sections.** `0059865` (6539), `HDL115` (115) and `npod6591` (6591) have
qptiffs but no GeoJSON export yet. The registration half runs on them unchanged as soon as
they exist.

## Files

| file | role |
|---|---|
| `geo.py` | the exports in microns, the islet size filter, the scale guard, polygon rasterisation |
| `validate.py` | the endocrine gene-set score, the neighbourhood test and its placement null |
| `hormone.py` | the channel map, per-cell sampling, AUROC + label-permutation null, the AF and saturation diagnostics, and `agreement_by_scale` |
| `run_hormone.py` | driver → `data/development/restain/hormone_panc.json` |
| `run_register.py` | driver → `data/development/restain/registration_panc.json` |
| `tests/` | the machinery on geometry whose answer is arithmetic — no export or bundle needed |
