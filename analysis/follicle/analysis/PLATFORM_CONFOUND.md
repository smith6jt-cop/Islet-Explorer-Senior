# The CODEX / Phenocycler Platform Confound

*Last updated 2026-08-11, after the region-layer rebuild.*

The fluorescent spleen cohort mixes two acquisition platforms whose imaged
fields differ by roughly an order of magnitude. This document is the single
reference for what that difference is, how it interacts with the rs3184504
genotype contrast, everything that has been done about it, and what remains
broken.

Scope: the CODEX vs Phenocycler comparison within the **fluorescent** cohort.
The separate question of how "tissue area" is defined in the H&E cohort versus
the fluorescent cohort is handled in `H19_meta_analysis.ipynb`, which
meta-analyses per-cohort effect sizes rather than pooling raw measurements, and
is not developed here.

---

# Part I — Summary for reporting

Spleen sections were imaged on two multiplex platforms. CODEX acquisitions
cover a fixed field of approximately 13.3 mm² at 0.374 µm/pixel, whereas
Phenocycler acquisitions are whole-section scans covering 189–433 mm² at
0.508 µm/pixel. After segmentation this yields a median of 12.7 mm² of
annotated tissue per CODEX donor against 123.4 mm² per Phenocycler donor, a
roughly tenfold difference, together with a comparable difference in the number
of follicles available per donor.

This difference is not a simple matter of statistical precision. It acts
through four distinct channels: unequal tissue sampled per donor; non-random
selection of which region of the spleen a small CODEX field captures;
a difference in optical resolution that changes what the segmentation models
resolve; and a difference in marker detection sensitivity, with CD20 detection
sensitivity of 0.42–0.54 in CODEX donors against 0.87–0.92 in Phenocycler
donors. The first is correctable by area normalisation; the remaining three are
not.

Critically, platform is not balanced across genotypes. In the annotation-level
cohort (14 donors) CODEX contributes 1 C/C, 2 C/T and 2 T/T donors while
Phenocycler contributes 3, 5 and 1 respectively, so 2 of the 3 T/T donors are
CODEX. Platform is nonetheless estimable there: corr(platform, dosage) = 0.288,
VIF 1.09, design matrix full rank.

The single-cell notebooks are the weaker case. They read `Cells.csv` directly
and so cover 11 donors, in which **both** CODEX donors are T/T — platform and
the T/T contrast are collinear and no platform-adjusted genotype effect is
identifiable. That restriction is now an artefact rather than a data limit: the
three remaining CODEX donors do have single-cell data, and the only thing that
kept them out was that their CSV lacks the QuPath Follicle/PALS/RedPulp
distance columns those notebooks used to assign regions. Region assignment is
point-in-polygon as of the 2026-08-11 rebuild, so the columns are no longer
required. See §2.

A range of strategies was applied. Every test was run both on the full cohort
and on Phenocycler-only samples, with agreement between the two required before
a result was treated as real. Region boundaries were re-derived from single-cell
marker expression rather than from pixel classifiers trained on Phenocycler
images, and vessel burden was expressed as a dimensionless area fraction; after
this change no platform difference remained detectable across nine metric ×
region combinations. Images were cropped to a common spatial extent and the
primary analyses repeated. Five explicit deconfounding approaches were
evaluated, including regression residualisation on a platform indicator and
platform-invariant feature selection; none survived permutation testing. The
current preferred treatment models counts directly with a tissue-area offset
rather than discarding tissue by cropping.

The conclusion carried forward is that platform is now *estimable* but not
*resolved*. The annotation cohort separates the two effects statistically;
mechanisms 2 to 4 remain uncorrected in both cohorts, and at 4 / 7 / 3 donors
only very large effects are detectable, so a null result is weak evidence of
absence rather than evidence of no effect. Analyses restricted to Phenocycler
donors are internally consistent but leave only one T/T donor. Answering the
question properly still requires a platform-balanced cohort with adequate
representation of each genotype, on the order of ten donors per group.

---

# Part II — Internal record

## 1. The four mechanisms

They are frequently collapsed into "CODEX images are smaller". They behave
differently and are mitigated differently.

| # | Mechanism | Magnitude | Correctable by area normalisation? |
|---|---|---|---|
| 1 | **Imaged field** | 13.3 mm² vs 189–433 mm² (~20×) | Yes — this is unequal effort |
| 2 | **Field selection** | CODEX is a fixed 3.76 × 3.58 mm tile | **No** — a choice of where to look |
| 3 | **Optical resolution** | 0.374 vs 0.5078 µm/px | **No** — changes what is resolved |
| 4 | **Detection sensitivity** | CD20 sens. 0.42/0.54 vs 0.87–0.92 | **No** — biases counts downward |

Two structural facts matter more than the fold-changes:

- **CODEX extent carries no between-sample variance.** The three `ForSH2B3`
  images span 3.76 × 3.58 mm to three decimal places; `1901HBMP004` and
  `HDL073` span 4.73 × 2.79 mm. These are fixed acquisition tiles, not
  tissue-driven extents. `residualize_platform()` fits a single additive offset
  on what is effectively a constant, from n=5.
- **Mechanisms 2 and 4 push in opposite directions.** If CODEX fields were
  placed on follicle-rich tissue, densities are inflated in the (mostly T/T)
  CODEX donors. If CODEX CD20 sensitivity is low, B-cells are undercounted in
  the same donors. Either way the pooled T/T estimate is not a clean genotype
  effect, and the sign of the residual bias is unknown.

## 2. Two cohorts, two different problems

This distinction is the single most important thing in this document, and
existing notebook prose conflates it. `H20`'s statement that "both CODEX donors
are T/T" is true of its cohort and false of `H1`/`H7`'s.

**Annotation cohort** — `H1`–`H13`, n=14, source `AnnotationsFinal.csv` +
`ForSH2B3.csv`:

| | C/C | C/T | T/T |
|---|---|---|---|
| CODEX | 1 | 2 | 2 |
| Phenocycler | 3 | 5 | 1 |

CODEX spans all three genotypes. Platform is strongly associated with T/T
(2 of 3 T/T donors) but is **statistically separable** —
`region_layer_identifiability.csv` records corr(platform, dosage) = **0.288**,
VIF **1.09**, rank 3 of 3.

**Single-cell cohort** — `H14`/`H15`/`H16`/`H17`/`H20`, n=11, because those
notebooks read `Cells.csv` directly. CODEX is `1901HBMP004` and `HDL073` only,
and **both are T/T**. Platform is collinear with the T/T contrast; the
platform-adjusted genotype effect is **not identifiable at n=11**. Confirmed by
`n_CODEX=2, n_PC=9` in `H15_platform_diagnostic.csv` and
`H16_platform_diagnostic.csv`, and by the empty CODEX-only stratum in
`H20_platform_stratified.csv`.

**This restriction is now removable, and that is the single most actionable
item in this document.** The three `ForSH2B3` CODEX donors (`HBMP006`,
`HBMP007`, `HDL098`) *do* have single-cell data — 289,500 / 284,836 / 205,792
cells in `Measurements/Spleen_Upload2_2024Sept17_v2_distances.csv`, which
`data_utils.load_single_cells` already merges and which the region-layer build
consumes for all 14 donors.

What kept them out was mechanical, not biological: that CSV carries
`Signed distance to annotation Sinusoid/Vessel µm` and **not** the
Follicle/PALS/RedPulp distance columns those notebooks took an argmin over to
assign regions. Since 2026-08-11 region assignment is point-in-polygon against
`analysis/regions/`, which covers all 14 donors, so the missing columns no
longer matter. Moving H15/H16/H17 from `pd.read_csv(CELLS_CSV, …)` to
`load_single_cells()` would take the single-cell cohort to n=14 with CODEX
spanning all three genotypes — the same identifiable design the annotation
cohort already has. That work has not been done; until it is, every
single-cell-cohort statement below stands at n=11.

## 3. Quantified disparity

**Acquisition** — from `data/*/server.json`:

| | CODEX | Phenocycler |
|---|---|---|
| Field of view | 13.3 mm² (all 8 acquisitions) | 188.7–433.0 mm² |
| Megapixels | ~95 | 732–1,680 |
| Pixel size | 0.374 µm | 0.5078 µm |
| Channels | 25–29 | 25–33 |
| Raw file | 1.1–3.5 GB | 24.6–44.5 GB |

CODEX is ~11.7× fewer pixels covering ~20× less tissue at 1.36× finer linear
resolution.

**Annotated tissue, post-fix** (see §6 — these numbers changed on 2026-08-05):

| Platform | n | Range (mm²) | Median (mm²) |
|---|---|---|---|
| CODEX | 5 | 12.1 – 14.2 | 12.7 |
| Phenocycler | 9 | 50.9 – 222.6 | 123.4 |

Median fold difference **9.7×**. Note this is *larger* than the 6.8× reported
before the Sinusoid fix, because CODEX areas were previously inflated by
double-counting.

## 4. Evidence the confound bites

| Source | Finding |
|---|---|
| `H20_nb_rate_ratios.csv` | `Platform_CODEX` RR = **0.485, p = 3.4e-06** in the primary per-follicle fit (≥8400 µm²) — the largest platform main effect in the project. **Caveat:** the all-follicle sensitivity fit reverses sign (RR = 2.79, p = 0.067), so the magnitude is not stable to the area filter; only the presence of a strong platform term is |
| `H20_cd20_calibration.csv` | CD20 sensitivity 0.421 (`1901HBMP004`), 0.538 (`HDL073`) vs 0.873–0.918 for all Phenocycler donors |
| `H10_statistical_tests.csv` | PC1 of the per-donor feature matrix vs platform, ρ = 0.592, **p = 0.033** |
| `H11_baseline_diagnostic.csv` | R²(platform) = 0.286; dosage ρ = 0.478 collapses to 0.130 (p = 0.67) under leave-one-out |
| `H11_approach_comparison.csv` | permutation p = 0.89 / 0.87 / 0.79 / 0.57 / 0.36 — **nothing survives** |
| `H15_platform_diagnostic.csv` | M2 Macrophage and "Other" separate the platforms perfectly (rank-biserial ±1.0, p = 0.036) at n=2 vs n=9 |
| `H20_verdict.csv` | verbatim: *"Confound controlled? — No — platform collinear with T/T"*. True of its n=11 cohort. Re-derive after H20 moves to `load_single_cells` (§2) — the collinearity is a property of the 11-donor restriction, not of the data |

Counter-evidence, i.e. where the confound does **not** show up:

| Source | Finding |
|---|---|
| `H1_statistical_tests.csv` | LMM platform term **p = 0.808** (n = 2169 annotations, 14 samples) |
| `H7` cell 12 | aggregate density CODEX vs Phenocycler, all five metrics p = 0.364–0.699 |
| `P1_perifollicular_platform_check.csv` | after the QuPath-sourced region layer + area fraction, the perifollicular platform check is clean. Supersedes `H14_platform_check.csv`, which came from the abandoned cell-derived KDE boundaries |
| `H16_platform_diagnostic.csv` | after RESTORE normalisation only Endothelial survives (p = 0.036), down from four cell types in H15 |

The pattern: platform is undetectable in simple per-region density aggregates,
and severe in multivariate axes, single-cell composition, and count-rate models.

## 5. Mitigation history

| # | Where | Approach | Verdict |
|---|---|---|---|
| 1 | H1 | Per-sample pixel-size audit; LMM platform term; Phenocycler-only arm | Platform term n.s. (p = 0.808) |
| 2 | H7 | Every test run All / Phenocycler-only / concordance | Standing practice; 0 SNPs survive FDR either way |
| 3 | H10 | Platform as a PCA input; platform-sensitivity PCA | Diagnosed, not fixed (ρ = 0.59) |
| 4 | H11 | Five deconfounding strategies incl. `residualize_platform()` | **All rejected** on permutation |
| 5 | H12 | Fixed-radius circular probes (r = 100/200/500 µm) + leave-one-image-out CV | Size-invariant sampling unit; **cropping removed 2026-08-11** |
| 6 | H13 | Crop to ~5000 µm, re-run H1/H2/H5 | Replication; stayed null. Discarded **75%** of annotations (120,110 of 483,145). **Cropping removed 2026-08-11** — H13 no longer differs from H1/H2/H5 except by re-including HDL018/HDL021 |
| 7 | H14 | Cell-marker KDE boundaries; 5000 × 3000 µm crop matched to CODEX; percentile rather than Otsu thresholds; dimensionless area fraction | Discarded **88%**. Both the boundaries and the crop are retired — see §7(f) |
| 8 | H15/H16 | Harmony (per-Image, deliberately gentler than per-Platform); RESTORE normalisation | RESTORE clears most of the cell-type confound |
| 9 | H20 | **Stop cropping.** Negative-binomial count-rate with area offset, donor cluster-robust SEs, platform-stratified inference, matched-window Monte Carlo | Concluded **non-identifiable** at n=11 |
| 10 | region layer | **Current position.** One geometry for every object: QuPath classifier polygons, scaled per donor, area-filtered, tissue-clipped. Region by containment, not by `Parent` or by signed distance | Removes the hierarchy defect and the CODEX area double-count, and unblocks the n=14 single-cell cohort (§2) |

The trajectory is worth noting: early work treated this as a nuisance to adjust
away (H7, H10, H11), the middle period treated it as a sampling problem to
equalise by cropping (H12, H13, H14, H15), and the current position is that
cropping throws away 75–88% of the data to solve only mechanism 1 while leaving
mechanisms 2–4 untouched, and that the honest answer is the design is
underpowered for the question.

**As of 2026-08-11 no analysis crops.** H20 had already stopped; H12, H13, H14,
H15, H16 and H17 followed. Every Phenocycler section is analysed at full extent.
The crop thresholds are set to `np.inf` rather than deleted, so the window
machinery still runs and returns "no crop" for every sample — verified by
replicating both window rules against `AnnotationsFinal.csv`, where the largest
extent in the cohort is 23,490 µm.

## 6. Defects found and fixed on 2026-08-05

**(a) CODEX RedPulp area was double-counted.** `data_utils.py` mapped both
`Red_Pulp` and `Sinusoid` to `RedPulp`, and `compute_density()` sums areas per
class. In the CODEX export these are *sibling* annotations (both parented to the
image root) covering the same territory — near-identical centroids, with
`Sinusoid` carrying a ~20× larger perimeter for a smaller area, i.e. the
sinusoidal network traced through the red pulp. Summed, they put total annotated
tissue at **132–149% of the image field of view**, which is impossible.

*Fix:* `_CODEX_CLASS_MAP` no longer maps `Sinusoid`, so its area is excluded;
`_CODEX_REGION_MAP` still does, so vessels parented to `Annotation (Sinusoid)`
are still counted under RedPulp. Numerator unchanged, denominator corrected.

*Effect* — confined to the three `ForSH2B3` images (verified: 0 other rows changed):

| Sample | RedPulp area | RedPulp density | Other regions' `Density_Normalized` |
|---|---|---|---|
| HBMP006 | 16.255 → 10.613 mm² | 460.0 → 704.6 | ×0.65 |
| HBMP007 | 16.684 → 12.370 mm² | 153.5 → 207.0 | ×0.74 |
| HDL098 | 16.721 → 10.864 mm² | 291.7 → 448.9 | ×0.65 |

Totals now sit at 91–107% of field of view.

**(b) `load_all_data()` mislabelled platform.** It blanket-assigned
`Platform = "Phenocycler"` to everything in `AnnotationsFinal.csv`, including
`1901HBMP004` and `HDL073`, which `CODEX_SAMPLES` declares CODEX. Impact was
latent — every downstream consumer re-derived platform from `CODEX_SAMPLES` —
but the trap was one line away. Now derived from `CODEX_SAMPLES` at source.
`H1` additionally carried a **stale local copy** of `load_all_data` that shadowed
the library version; removed.

**(c) Stale statistic in `H14_cell_boundaries.md`.** Cited "H1 found borderline
significant platform differences for Follicle density, p=0.055". No such value
exists in H1 or its outputs. Corrected to the actual p = 0.808, with a note.

**(d) Stale hardcoded sample counts in H7.** Seven `n=15` / `n=10` literals in
figure titles and printed headers, against actual n=14 / n=9. Replaced with
expressions derived from the data.

## 7. Defects — status after the 2026-08-11 rebuild

**(a) `1901HBMP004` has a broken annotation hierarchy. — RESOLVED 2026-08-11,
by routing around it rather than repairing it.**

In `AnnotationsFinal.csv`, only its `RedPulp` annotation received children. Of
4,710 SmallVessels, 3,696 are parented to RedPulp and **1,014 are stranded at
`Root object (Image)`** (21.5%; every other image strands 0.0–0.3%).
`Follicle`, `PALS`, `Trabeculae` **and** `LargeVessel` all have **zero** vessel
children, despite 35 / 5 / 44 / 18 annotations existing respectively.

That hierarchy is still broken in the CSV and there is still no plan to re-run
QuPath's parenting pass. It no longer matters, because nothing reads it: region
assignment is point-in-polygon throughout, via
`data_utils.assign_regions_by_image`. `1901HBMP004` now reports **811**
Follicle-resident SmallVessels where the `Parent` path reported 0.

The general point outlives this one image. `Parent` did not fail randomly — it
failed for one donor, by 21.5% against a cohort norm under 0.3%, and that donor
is one of only three T/T. A per-image failure of that size in a 14-donor design
is a genotype effect in everything but name. The `Parent` route is retained only
behind `get_vessels(..., from_parent=True)`, for reproducing old numbers.

Its stale GeoJSON export is also gone. It carried 54 Follicle and 28 PALS
objects against the 35 and 5 the CSV records, and reproduced only 0.82× of the
CSV Follicle area and 1.59× of PALS; re-exported from QuPath 0.7.0 it matches
the CSV exactly and reproduces 0.97× / 0.99×.

So three of the four main regions evaluate to a structural density of **0.0** —
not a measured value — for one of only three T/T donors. `H1` now warns rather
than silently reporting 0. Fixing it requires re-running the QuPath parenting
pass; it is not a Python change.

### Why this was hard to see

Two reasons, both worth remembering as failure modes rather than as history:

1. **Most committed outputs predated the defect.** Outputs lacking `HDL079`
   were generated when `1901HBMP004` was still intact and showed healthy values
   — committed `H10_feature_matrix.csv` had it at Follicle 354.7 / PALS 37.3 /
   Trabeculae 1196.7 where the broken data yielded 0 / 0 / 0. Presence of
   `HDL079` is a reliable staleness test for fluorescent-cohort tables (not for
   `HE_*` tables, where HDL079 is an unrelated H&E donor).
2. **Pooled per-annotation plots hid it even when current.** In H1 the T/T group
   pooled 510 annotations, of which only 35 were `1901HBMP004`'s zeros, so the
   group median was unaffected. The zero was visible only in per-sample
   aggregates. A structural zero that survives pooling is the kind of defect
   that reaches a figure.

The exposure map that used to sit here — which notebooks resolved
vessel→region through `Parent`, and which assigned spatially — no longer
describes anything. Every notebook assigns spatially. `data_utils.get_vessels`,
`compute_density` and `assign_vessels_to_follicles` all take the polygon route,
and `compute_density` no longer sums annotation areas either, which removes the
CODEX Red_Pulp/Sinusoid double-count at the same time.

**(b) The committed H1 tables were stale.** `analysis/tables/H1_*.csv` had been
generated against an earlier `AnnotationsFinal.csv`: they omitted `HDL079`
entirely and reported `1901HBMP004` with 54 follicle annotations and 801 vessels
against the current 35 and 0. Re-running moved n from 13 to 14 samples and 1903
to 2169 annotations. **Most of the change in H1's p-values on 2026-08-05 comes
from this data refresh, not from the Sinusoid fix.** Assume other committed
tables may be similarly stale; `analysis/tables/` is gitignored apart from
`H1_statistical_tests.csv`, so staleness is invisible to version control.

**(c) `Peripheral_White_Pulp` handling is questionable. — RESOLVED
2026-08-11.** It was mapped to `PALS`, but peripheral white pulp is closer to
marginal zone than to periarteriolar lymphoid sheath, and some PWP annotations
nest inside Follicle/PALS/Sinusoid annotations and were double-counted into
PALS, adding 0.37–0.63 mm² per CODEX image.

The region layer maps it to **Perifollicular** instead
(`build_region_layer_qupath.CLASS_TO_REGION`), which is what it anatomically is:
the zone around the follicle, where PALS sheaths the central arteriole. Folding
it into PALS had given the CODEX donors a collar hugging every follicle that
then consumed the perifollicular band — HBMP006 had a band on 2 of 7 follicles,
HDL098 on 1 of 8 — and doubled their PALS area (HBMP006: 0.90 mm² reported as
1.76). `data_utils._CODEX_REGION_MAP` still folds it into PALS for the
annotation-CSV path; the difference is deliberate and `tests/test_region_layer.py`
asserts both.

**(d) Two Phenocycler RedPulp annotations are mis-parented.** In `HDL055` and
`HDL086` the single RedPulp annotation (106 and 112 mm²) is parented to a
`PALS` / `Follicle` annotation far smaller than itself. Harmless to current
sums, but it means **a "use only root-parented annotations" rule would silently
delete ~90% of those two images' tissue** — do not implement one.

**(e) Corrupt `server.json` entries.** `data/10`–`data/13` and `data/17`
re-register `1901HBMP004_PC29` / `HDL073_PC19` with `pixelWidth = 1.0 µm` and
channels loaded as timepoints. They carry no annotations, but any code walking
`data/*/server.json` for calibration picks up a 2.67×-wrong pixel size. Also
`data/15` (HDL070) has 1.65M cells but only 5 annotations; `data/19` is the
annotated copy.

## 8. What would actually resolve it

**Available in this dataset, and not yet done:** move H15/H16/H17 (and H20) to
`load_single_cells()` so the single-cell cohort becomes n=14 with CODEX
spanning all three genotypes. This is now only an ingest change — see §2. It
converts the single-cell arm from non-identifiable to identifiable, and it
costs nothing but the edit.

**Not available in this dataset:** enough donors. Even at n=14 the annotation
cohort supports platform adjustment with 1 Phenocycler T/T donor, and 4 / 7 / 3
detects only very large effects. Mechanisms 2 to 4 — field selection, optical
resolution, detection sensitivity — are not correctable by any amount of
normalisation and remain uncorrected. The requirement is a **platform-balanced
cohort with on the order of 10 donors per genotype**, or re-imaging the existing
CODEX donors on Phenocycler so that platform and genotype are crossed rather
than nested.

Defensible claims today: Phenocycler-only analyses (accepting n_T/T = 1); area
fractions and densities within the region layer, whose geometry is now common to
both platforms; and H13's crop-matched replications. Still indefensible: any
platform-adjusted genotype effect from the **n=11** single-cell cohort, and any
reading of a null as evidence of no effect.

## 9. Reproducing the numbers

Use `/blue/maigan/smith6jt/miniforge3/envs/xenium_analysis/bin/python`
(Python 3.12.13, numpy 2.0.2). This is the **only** environment with full
coverage; KINTSUGI, which earlier revisions of this document named, lacks
scanpy/anndata/harmonypy/leidenalg and cannot run H15/H16/H17 at all. Notebooks
must run with the working directory set to `analysis/`.

The region layer must exist before any of this runs — `analysis/regions/` is
gitignored, so a fresh clone has to build it:

```bash
QP=/blue/maigan/ebutterworth/software/QuPath/bin/QuPath   # 0.7.0, headless
"$QP" script -p "$PWD/project.qpproj" "$PWD/scripts/export_all_regions_geojson.groovy"
python scripts/radial_follicle_profiles.py     # -> radial_pfz_width.csv
python scripts/build_region_layer_qupath.py    # -> analysis/regions/
python scripts/validate_region_layer.py
```

```python
import sys; sys.path.insert(0, '.')
import data_utils as du, pandas as pd
den = du.compute_density(du.load_all_data())
t = den.groupby('Sample').agg(tissue_mm2=('Region_Area_mm2', 'sum'),
                              Genotype=('Genotype', 'first')).reset_index()
t['Platform'] = t['Sample'].apply(
    lambda s: 'CODEX' if s in du.CODEX_SAMPLES else 'Phenocycler')
pd.crosstab(t['Platform'], t['Genotype'])   # -> 1/2/2 and 3/5/1
```

Spot-checks that should hold: every CODEX acquisition is 13.3 mm² of field of
view; CODEX annotated tissue totals 91–107% of that; `2006_CC2B` RedPulp is
10.613 mm² (not 16.255).
