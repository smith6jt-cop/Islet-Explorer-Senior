# Post-Xenium IF restain lane

Registers the post-Xenium IF to its own Xenium section **nucleus by nucleus**, then checks the
result against signals the registration never saw.

The post-Xenium IF images the *same physical section* Xenium measured. So the Xenium segmentation
**is** the segmentation — no QuPath run on the re-stains, no new cell objects — and the two DAPI
images show the same nuclei twice. That makes this a point-set problem rather than an intensity
one, which is the easier of the two: it needs no shared intensity model, and its residual is a
**per-cell distance in microns** rather than a per-window statistic.

```bash
cd core/phenocycler
../../.venv/bin/python scripts/restain/export_dapi.py   --donors 6436 6442   # the two DAPI planes
~/.venvs/valis/bin/python scripts/restain/valis_register.py --donors 6436 6442   # the seed
../../.venv/bin/python scripts/restain/nuclei.py        --donors 6436 6442   # the two point sets
../../.venv/bin/python scripts/restain/run_refine.py    --donors 6436 6442   # fit, score, keep
../../.venv/bin/python scripts/restain/run_align.py --transform refined      # the gate
../../.venv/bin/python scripts/restain/refine_figures.py --donors 6436 6442
../../.venv/bin/python scripts/restain/align_figures.py --transform refined
../../.venv/bin/python -m pytest -q scripts/restain/tests                    # 52 passed
```

Inputs are the maintainer's QuPath exports under `~/IO60panc2nd/PostXeniumIHC/`, the IF qptiffs,
and the Xenium bundle. Output goes to `data/development/restain/`. Nothing here writes a
`StageManifest` or touches `data/runs/`, and nothing under `scripts/**` is hashed into `run_id`.

## Where it got to

Per-cell: the distance from a Xenium nucleus to the IF nucleus it lands on, over every nucleus in
the section.

| donor | tissue outlines (ECC) | VALIS, as it stood | **nucleus-matched** | p90 | within 2 µm | null |
|---|---:|---:|---:|---:|---:|---:|
| 6436 | 42.13 µm¹ | 2.84 µm | **0.546 µm** | 1.571 | 53.1 % | 8.1 % |
| 6442 | 31.46 µm | 1.36 µm | **0.630 µm** | 1.653 | 57.1 % | 8.0 % |

¹ window residual; only 25 of 243 windows were even measurable at that alignment.

`within 2 µm` is the share of **all** Xenium nuclei with an IF nucleus that close. Its ceiling is
the IF detector's own recall — 56.7 % and 61.0 % — so 53.1 % and 57.1 % are 94 % of what was
available to match, and the placement null is 8 %.

The gate (`run_align.py`, below) returns **PASS on both donors**, at a median window residual of
**0.11 µm (6436)** and **0.13 µm (6442)** with every window measurable, against a bar of 1.80 and
1.88 µm.

## How, and what was wrong before

**The registration is `register.icp` on the two nucleus fields**, seeded by VALIS. Xenium's side is
`nucleus_boundaries.parquet` — already a segmentation of Xenium's own DAPI — and the IF side is
detected by `nuclei.py`: a band-pass at the ~6 µm nucleus scale, local maxima no closer than 4 µm,
and a parabolic sub-pixel fit. 395,048 and 584,723 IF nuclei inside the drawn tissue. No cell
objects and no cell types are created; these points exist only to be matched.

The seed only has to land inside half the spacing between nuclei (~4.5 µm) for the mutual
nearest-neighbour matching to be unambiguous. Three defects had been keeping VALIS outside that.
The third is pinned by tests; the first two are inside VALIS's own call and cannot be tested from
the pipeline venv, so what stands for them is VALIS's own reported error, recorded below and in
`valis_result.json`:

**1. VALIS's rigid transform had never been refined above a thumbnail.** `Valis.register` fits `M`
on the processed images — 512 px on a 20,880 px slide, about 40 µm/px. The only thing in
valis-wsi 1.2.0 that re-fits `M` is `MicroRigidRegistrar`, which `Valis.__init__` defaults to
`None` (registration.py:1790) and which `register()` runs between the rigid and non-rigid stages
(registration.py:4743-4745). Enabling it: `N low rez matches = 1372, N high rez matches = 31596,
low rez D = 1.190, high rez D = 0.231`, and VALIS's own error fell

    6436   rigid_D 15.28 -> 3.892 um      non_rigid_D 20.50 -> 3.855 um
    6442   rigid_D  6.89 -> 3.331 um      non_rigid_D 16.33 -> 3.414 um

Note the second column. The non-rigid stage used to make VALIS's own metric 1.3-2.4x **worse**
than the rigid fit it started from; it no longer does.

**2. The micro-registration pass was inert.** `register_micro` was called without
`reference_img_f`, so it forwarded `None` to `serial_non_rigid.register_images`, where
`warp_tools.get_ref_img_idx` returns index 0 for two images (warp_tools.py:57-59) — the **Xenium**.
That made the Xenium the fixed image, so it was handed a field of zeros
(serial_non_rigid.py:778-781) while the field computed for the IF was skipped by the update loop,
which skips the *Valis object's* own reference — the IF (registration.py:4941-4943). The two
references disagreed, the computed field was discarded, and the only surviving effect was that the
**old** field got resampled: measured, `non_rigid_D` went 20.497 → 21.331 on 6436, i.e. slightly
worse. Passing the argument makes the pass do what its name says. It would not have fixed the
residual on its own in any case: every path in `register_micro` ends at
`bk_dxdy = current + new` and none of them writes `M`, so it cannot remove a translation.

**3. Neither pyramid is `micron = pixel * um_per_px`, and they are wrong differently.** Measured on
donor 6436, three 384-px tiles, level 0 against level 2:

| | box mean | best decimation |
|---|---:|---:|
| Xenium `morphology_focus` | 0.9766 | **0.9936** at phase (0,0) |
| IF qptiff (DAPI page) | **0.99992** | 0.9891 |

So the Xenium pyramid *keeps* level-0 pixel `f·k` and the IF pyramid *averages* `f·k … f·k+f−1`.
Their pixel centres sit at `(f·k + 0.5)·base` and `(f·k + f/2)·base`, i.e. `origin_um` of 0.10625
and 0.5. Ignoring both left a fixed **+0.39 µm per axis** in the warped result. Separately,
`shapes[0][0] / shapes[lvl][0]` is not the downsample factor: on 6436 it is 68119/17029 = 4.000176
where the factor is exactly 4, a scale error worth 0.61 µm at the far corner. Both are now derived
from the integer factor with an assertion, and both are recorded in `export.json`.

## Every candidate, scored the same way

All of them start from the same seed and are scored by the same per-cell distance, and the losers
are reported rather than dropped (`refine_candidates.png`):

| candidate | package | 6436 | 6442 |
|---|---|---:|---:|
| `valis_affine` — the seed | valis-wsi 1.2.0 | 1.190 | 1.830 |
| `valis_field` — affine + its 250 µm binned residual | valis-wsi | 1.486 | 1.889 |
| `nuclei_translation` — seed + one constant | scipy | 0.801 | 1.809 |
| `nuclei_icp_similarity` | `register.icp` | 0.791 | 1.805 |
| **`nuclei_icp_affine`** | `register.icp` | **0.546** | **0.630** |
| `image_translation` — whole-section NCC | numpy FFT | 0.801 | 1.809 |
| `phase_translation` — phase correlation, upsample 100 | scikit-image | 2.259 | 1.891 |
| `nuclei_icp_affine+local` — plus a 250 µm field | scipy | 0.522 | 0.600 |

Four things this table says that a single number would not:

* **The affine matters and is not shear-hunting.** It cuts the residual 1.4x on 6436 and 2.9x on
  6442 against the similarity fit, and it does so by *reducing* anisotropy — 643 → 280 ppm on
  6436 and 1338 → 419 ppm on 6442. It is removing a stretch VALIS had, not adding one.
* **The local field does not earn its place.** It buys 0.02-0.03 µm, inside the declared 0.05 µm
  tolerance, so `choose()` keeps the simpler transform. More freedom cannot score worse on the
  data it was fitted to, which is exactly why the tolerance is declared rather than "lowest wins".
* **`valis_field` is now worse than `valis_affine` on both donors.** A 250 µm median-binned
  residual field with a third of its bins filled from a neighbour was compensating for a bad
  affine; against a good one it only adds noise.
* **Three mechanisms were asked for the same translation.** Nearest-neighbour nuclei and
  whole-section NCC agree to **0.046 µm (6436)** and **0.092 µm (6442)** — two different
  mechanisms, one answer. Phase correlation is the outlier, 2.15 and 0.53 µm from the other two,
  and it is reported rather than hidden: normalising the cross-power spectrum weights every
  spatial frequency equally, and a sparse point raster against a continuous DAPI image is exactly
  the case where that hurts.

SimpleITK 2.5.6 (in `~/.venvs/valis`) and itk-elastix 0.19.2 (in `~/miniconda3/envs/image2image`)
are installed and were **not** run. They would refine the same seed by the same intensity mechanism
as the two image candidates above, which the point-set fit beats by 1.5x (6436) and 2.9x (6442)
on the same score.
That is a judgement, and it is recorded here so it can be reversed.

## Is it a registration, or a fit?

Mutual nearest-neighbour matching can be self-consistently **wrong**: if the seed were out by more
than half the spacing between nuclei, every point would pair with a stranger and the fit would
converge confidently onto that, with nothing in the in-sample score showing it. So the section is
split in two and each fit is scored on the half it never saw:

| donor | fitted on | in-sample | **held out** |
|---|---|---:|---:|
| 6436 | left half (348,232 nuclei) | 0.542 | **0.614** |
| 6436 | right half | 0.545 | **0.585** |
| 6442 | left half (479,459 nuclei) | 0.657 | **0.605** |
| 6442 | right half | 0.603 | **0.671** |

Held-out is within 0.07 µm of in-sample, and on 6442 one fold's held-out score is *better* than its
in-sample one. The two halves' transforms disagree by 0.184 µm (6436) and 0.240 µm (6442) of
translation and by 19 and 2 ppm of scale. The transform transfers, so it is a property of the
sections rather than of the points it was fitted to.

## The gate, and what it is now worth

`run_align.py` warps Xenium nucleus **polygons** into the IF frame, rasterises them at 0.5 µm/px and
finds, per 500 µm tissue window, the shift that best lines that mask up with the IF DAPI by
normalised cross-correlation. Its thresholds are geometry, not preference: PASS at the displacement
where a disc of the donor's own median cell still covers 75 % of its own footprint.

| donor | windows | median | p90 | precision | DAPI contrast inside vs gaps | verdict |
|---|---|---:|---:|---:|---:|---|
| 6436 | 243/243 | **0.109 µm** | 0.233 | 0.015 | 2.04 → **5.14×** | PASS (bar 1.795) |
| 6442 | 365/365 | **0.134 µm** | 0.260 | 0.014 | 3.51 → **4.98×** | PASS (bar 1.880) |

**Read this as consistency, not independence.** For the tissue-outline and VALIS transforms the
gate was genuinely independent — neither had ever seen a nucleus. The nucleus-matched fit and the
gate now optimise closely related things, so the gate confirming at a tenth of a micron is expected.
The checks that are *not* circular are the hold-out above and the placement null.

**Which diameter anchors the bar is a decision, and all three are reported.** The gate takes its
diameter from `cell_area`, but it correlates *nucleus* polygons, and the hormone lane samples an
8 µm disc:

| | 6436 | 6442 |
|---|---:|---:|
| cell (9.08 / 9.51 µm) — the declared bar | 1.795 | 1.880 |
| nucleus (6.19 / 6.13 µm) — what the gate actually matches | 1.223 | 1.211 |
| the hormone sampler's 8 µm aperture | 1.581 | 1.581 |

The per-cell median (0.546 / 0.630 µm) is under all three, and so is the per-cell p90 for the cell
and aperture bars. The declared cell bar is what `run_align.py` reports, because it was declared
before the run and is not renegotiated in either direction.

**Tissue Dice went slightly down**: 0.8903 → 0.8813 (6436) and 0.9509 → 0.9497 (6442). That is not a
regression, it is the point. The ECC transform was *fitted to maximise* tissue-outline overlap, so
it wins on Dice by construction while sitting 42 µm out at cell scale. Dice is not the criterion.

## Where the remaining 0.55 µm is, and where it is not

* **Not a translation.** After the affine the median residual vector is (−0.006, +0.000) µm.
* **Not spatial structure.** A 250 µm field with full freedom to bend changes the median by 0.02 µm,
  and the residual quiver has no direction anywhere on the section (`refine_residual.png`).
* **Not peak localisation.** Detecting the IF nuclei at 0.5 µm/px instead of 1.0 µm/px moves the
  residual to 0.563 (6436) and 0.650 µm (6442) — very slightly *worse*, not better. Doubling the
  resolution buys nothing, so the limit is not how precisely a peak can be located.

What is left is the difference between where Xenium's nucleus-polygon centroid sits and where the
IF DAPI blob centre sits: two definitions of one nucleus, plus whatever genuinely differs between
two imagings of the same section. No transform of any class removes that.

**0.5 % of Xenium nuclei have no IF nucleus anywhere near them.** They are outside the drawn IF
`Tissue` polygon — on 6436 the drawn tissue stops 1.3 mm short of where the Xenium cells reach —
so there is nothing there to register to. Those are recorded as far distances rather than dropped,
and they are what the tail of `refine_pairs.png` is.

## The `Islet` class needs a size filter, and it is not optional

The classifier emits thousands of sub-cell fragments beside the real islets:

| donor | objects | median area | ≥50 µm across | share of islet area |
|---|---:|---:|---:|---:|
| 6436 | 4,671 | 27 µm² (5.9 µm) | **125** | 84.6 % |
| 6442 | 3,366 | 23 µm² (5.4 µm) | **95** | 81.1 % |

The median object is one cell wide. Filtered, the counts land in the range the PhenoCycler cohort
shows (67–1,304 islets per donor).

## The downstream checks have NOT been re-run

Everything below was measured on the **tissue-outline** transform, which this work now measures at
31–42 µm — three to four cell diameters. They are kept for the record and **must not be read as
current results**; re-running them on the new transform is deliberately out of scope for this pass.

*IF islets against Xenium endocrine RNA* (`validate.py`, `registration_panc.json`): 6442 17.30×
inside vs elsewhere, 6436 2.83×, against a placement null of ~0.117. 6436's alpha and beta moved
together (2.83× / 2.89×) where 6442's did not (6.33× / 25.38×), and `validate.py`'s own docstring
calls that the signature of a mask landing partly off-target.

*The hormone check* (`hormone.py`, `hormone_panc.json`): INS ranks beta-like at AUROC 0.533 (6436)
and 0.574 (6442) against a pre-registered bar of 0.75, and `agreement_by_scale` climbs monotonically
out to 120 µm on 6442. `hormone_panc.json` also records `um_per_px = 1.0018796` — the qptiff
resolution tag, not the recorded 0.25 — so it carries a second defect besides the transform and is
not usable as a "before" without regenerating it.

## Files

| file | role |
|---|---|
| `export_dapi.py` | the two DAPI planes, their pyramid factors and their pixel origins |
| `nuclei.py` | IF nucleus detection; Xenium nucleus centroids from the bundle |
| `refine.py` | the candidates, the per-cell score and its placement null |
| `run_refine.py` | driver, the parsimony rule and the spatial hold-out → `refine_panc.json` |
| `refine_figures.py` | candidate comparison, per-cell distribution, residual quiver |
| `valis_register.py` | the seed (VALIS venv) → `warped_if_um.npy` |
| `valis_transform.py` | VALIS's point map as a `Transform`, for the comparison candidate |
| `align.py` / `run_align.py` | the gate → `alignment_panc_refined.json` |
| `align_figures.py` | the gate's own QC sheets |
| `geo.py` | the exports in microns, the islet size filter, the scale guard |
| `tiles.py` / `rasterize_labels.py` / `display.py` | windowed reads, rasterisation, honest ranges |
| `validate.py` / `hormone.py` / `run_hormone.py` | the downstream checks (not re-run — see above) |
| `run_register.py` | the tissue-outline registration; also scores a supplied transform |
| `tests/` | 52 tests, all on synthetic geometry — no export or bundle needed |

`data/development/restain/valis_baseline/` holds the 2026-08-31 arm this replaced.
