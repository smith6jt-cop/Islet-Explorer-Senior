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

## What this does not yet do

The hormone AUROC — *does protein INS, on the same cell, rank Xenium-called Beta cells above
non-Beta?* — needs the **Opal ↔ antibody assignment**, and that is in none of the five qptiffs:
they carry `DAPI, Opal 570, Opal 690, Opal 520, Sample AF` and **no `<Biomarker>` tag**. Three
stainable channels, unlabelled. Only the maintainer can supply that, and until then this lane
reads geometry, never intensity.

## Files

| file | role |
|---|---|
| `geo.py` | the exports in microns, the islet size filter, the scale guard, polygon rasterisation |
| `validate.py` | the endocrine gene-set score, the neighbourhood test and its placement null |
| `run_register.py` | driver → `data/development/restain/registration_panc.json` |
| `tests/` | the machinery on geometry whose answer is arithmetic — no export or bundle needed |
