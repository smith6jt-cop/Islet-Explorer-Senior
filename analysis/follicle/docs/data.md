# Data formats and conventions

> This is the detailed source-era inventory. After consolidation, all large
> paths are relative to `$FOLLICLE_DATA_ROOT` and generated paths are relative
> to `$FOLLICLE_OUTPUT_ROOT`; see `../README.md`. Sample identity is governed by
> `../config/sample_manifest.csv`.

This document describes every file under `Images/`, `Measurements/`,
`Groups.xlsx`, and the supporting metadata that the analysis notebooks
read. It also defines the sample-ID conventions used throughout.

## Images

### Phenocycler / CODEX OME-TIFF (`Images/`)

The fluorescent cohort comprises **9 unique Phenocycler donors**
(HDL011, HDL043, HDL052, HDL053, HDL055, HDL063, HDL070, HDL086,
1901HBMP004) plus **5 CODEX donors** (HBMP006, HBMP007, HDL073, HDL098,
1901HBMP004 — the last is the same donor as the Phenocycler entry,
imaged on both platforms). HDL073 additionally appears as two
Phenocycler OME-TIFF variants (`HDL073_PC19` pre-signal-isolation,
`HDL073_PC29` post-signal-isolation) — both files, one donor. Common
specs:

- **Format:** pyramidal OME-TIFF, 6 levels (1×, 2×, 4×, 8×, 16×, 32×),
  512×512 tiles, ADOBE_DEFLATE compression (code 8), uint16,
  big-endian.
- **Pixel size:** varies by donor, and **is not readable from the OME-XML**.
  Always use `data_utils.get_pixel_size`:
  0.5077663810243286 µm for the Phenocycler donors;
  **0.374** for `1901HBMP004` and `HDL073`, whose headers wrongly declare
  0.5078 (a stale value written by `convert_to_ome_tiff.py`);
  0.37740007578193524 for `HBMP006`, `HBMP007`, `HDL098`.
  Getting this wrong scales every polygon area by 1.8×.
- **Axes:** CYX — one IFD page per channel.
- **Channel count:** 19–33 per slide.
- **Channel order:** DAPI first, then alphabetical.
- **OME-XML:** Manually constructed with `BigEndian="true"`,
  `Interleaved="false"`, `<LightPath/>` per channel, per-channel
  `<TiffData>` blocks. µm is encoded as raw UTF-8 (`\xc2\xb5`), not
  the XML entity `&#181;m`. OME-XML is passed to tifffile as `bytes`
  to bypass the ASCII check.

Platform for analysis is decided by `data_utils.CODEX_SAMPLES` membership, not
by the file the data came from. The analysis cohort is 9 Phenocycler + 5 CODEX
donors = 14:

| Sample | Channels | Analysis platform |
|---|---|---|
| HDL011, HDL043, HDL052, HDL053, HDL055, HDL063, HDL070, HDL079, HDL086 | varies | Phenocycler |
| HBMP006, HBMP007, HDL098 | varies | CODEX |
| 1901HBMP004 | 29 | CODEX |
| HDL073 (`_PC19` pre-, `_PC29` post-signal-isolation — two files, one donor) | 19 / 29 | CODEX |

### HiperGator single-channel TIFs (`FromHipergator/`)

Raw outputs from the HiperGator processing pipeline. Three directories:

| Directory | ALT ID | HANDEL ID | Output OME-TIFF | Channels | Status |
|---|---|---|---|---|---|
| `19-001_SP_CC2-A28` | 1901 | — | `1901HBMP004_PC29.ome.tiff` | 29 | converted (1.74 GB) |
| `19-002_spleen_CC2-A_D200210` | 1902 | HDL073 | `HDL073_PC19.ome.tiff` | 19 | converted (2.44 GB) |
| `HDL73_SPL_Processed/ImageJ` | 1902 | HDL073 | `HDL073_PC29.ome.tiff` | 29 | signal-isolated (3.30 GB) |
| `19-002_spleen_CC3-C` | — | — | *skipped* | — | — |

Conversion is done by `scripts/convert_to_ome_tiff.py`.

### HDL73 signal isolation

`scripts/process_hdl73_channels.py` performs autofluorescence subtraction
using matched blank pairs via the KINTSUGI `kintsugi.signal` module.
Pass `--force` to re-process already-isolated channels.

Blank position mapping (from
`FromHipergator/HDL73_SPL_meta/channelnames.txt`):

| Position | Blank pair | Markers |
|---|---|---|
| a | Blank1a + Blank13a | CD20, CD31, CD34, CD35, Lyve1, PanCK, SMActin |
| b | Blank1b + Blank13b | CD8, CD15, CD21, CD44, CD45RO, CD5, CollagenIV, ECAD, FoxP3, Ki67, Podoplanin |
| c | Blank1c + Blank13c | CD3e, CD4, CD11c, CD107a, CD163, CD1c, CD45, CD68, HLADR, Vimentin |

- **Failed markers:** CD1c, CD5 (processed but may have poor signal).
- **Low signal-preservation warnings:** PanCK (0.14), Podoplanin (0.19),
  SMActin (0.18).
- **DAPI** is copied directly (no subtraction).

### InstanSeg cell-detection channels (26)

DAPI, CD45RO, Ki67, FOXP3, CD38, CD20, CD4, CD44, CD31, CD11c, CD34,
CD107a, PDL1, CD163, HLA-DR, CD68, CD8, CD21, CD66, CD141, CD57, CD3e,
HLA-A, PD-1, CD45, Podoplanin.

## Measurements

### `Measurements/AnnotationsFinal.csv`

The current primary analysis input — ~499K rows, 13 images, 17 columns.
Supersedes `AllAnnotations.csv` (which only covered 9 images).

| Column | Type | Notes |
|---|---|---|
| `Image` | str | OME-TIFF filename, e.g. `HDL011_PC33.ome.tiff` |
| `Object ID` | UUID | unique per annotation |
| `Object type` | str | always `Annotation` here |
| `Name` | str | optional human label |
| `Classification` | str | Follicle / PALS / RedPulp / Trabeculae / LargeVessel / SmallVessel |
| `Parent` | str | encodes the annotation hierarchy, e.g. `Annotation (Follicle)` |
| `ROI` | str | `Polygon`, `MultiPolygon`, … |
| `Centroid X µm`, `Centroid Y µm` | float | in µm |
| `Area µm^2` | float | in µm² |
| `Perimeter µm` | float | in µm |
| `Num Detections` | int | detections inside the annotation |
| `Length µm` | float | for linear ROIs |
| `Circularity` | float | 4π·Area / Perimeter² |
| `Solidity` | float | Area / ConvexHullArea |
| `Max diameter µm`, `Min diameter µm` | float | for ellipse-fit features |

The `Parent` column is the key for hierarchical aggregation: SmallVessel
annotations are children of region annotations, so vessel density per
region is `count(SmallVessel where Parent == region) / region area`.

### `Measurements/Cells.csv`

Single-cell table: 22.6M cells, 13 images, 71 columns, 11.8 GB on disk.

- 42 marker intensities (`Cell: {marker}: Mean`)
- 6 signed distances to region annotations (negative = inside)
- Per-cell morphology (Area, Circularity, …)
- Per-nucleus morphology

Read by H11, H14, H15.

### `Measurements/ForSH2B3.csv`

CODEX cohort annotations: 17.5K rows, 3 images
(`2006_CC2B`, `2007_CC1A`, `2008_CC3A`). Uses different region names
than the Phenocycler cohort (e.g. `Red_Pulp`, `Sinusoid`, `Trabecula`,
`Vessel`, `Peripheral_White_Pulp`). Harmonized into the Phenocycler
vocabulary by `_CODEX_REGION_MAP` / `_CODEX_CLASS_MAP` in
`analysis/data_utils.py`.

### `Measurements/HE_follicles_unet.csv`

Output of `scripts/predict_HE_unet.py`. Schema matches
`27_donor_spleen_measurements1.csv` exactly so it can be swapped into
H18 with a one-line change.

### `27_donor_spleen_measurements1.csv` (repo root)

H&E cohort, 27 donors balanced 9 / 9 / 9 by rs3184504. Per-follicle
records exported from a separate QuPath project (the polygons were
drawn by hand). Columns:

```
Image, Object ID, Object type, Name, Classification, Parent, ROI,
Centroid X µm, Centroid Y µm, Area µm^2, Perimeter µm
```

Four donors overlap the fluorescent cohort: HDL055, HDL063, HDL073,
HDL086.

## H&E U-Net training data

`scripts/export_HE_training_pairs.groovy` writes one training triple per
H&E QuPath entry (those with ≥ 1 `Follicle` annotation) into `data/`, all
resampled to a common **2 µm/px** grid:

- `data/images/<donor>.png` — downsampled RGB H&E.
- `data/masks/<donor>.png` — single-channel uint8 binary Follicle mask
  (255 = follicle, 0 = background; antialiasing off).
- `data/meta/<donor>.json` — calibration sidecar: `native_px_size_um`,
  `target_px_size_um` (2.0), `downsample`, native/export dimensions, and
  `n_follicle_annotations`.

These are gitignored (PNGs can be hundreds of MB). `data/` therefore holds
two unrelated things: the numbered `data/[1-9]/` QuPath segmentation output
and these `images/masks/meta/` U-Net subdirs.
`scripts/make_synthetic_he_data.py` writes the same triple layout for the
test suite and the CI smoke-train step.

### `analysis/HE_unet_models/`

Training output from `scripts/train_HE_unet.py`:

- `fold_{i}.pt` — best checkpoint per fold (state_dict + encoder name +
  best epoch + per-donor val IoU).
- `cv_summary.json` — per fold: best IoU, the **full per-epoch history**
  (loss / val IoU / lr), early-stop flag, and best-checkpoint follicle
  count + mean area vs. ground truth.
- `pooled_metrics.json` — per-donor predicted-vs-GT count/area plus their
  pooled Spearman ρ (written only with ≥ 3 validation donors, i.e. the
  `--folds 5` run). ρ > 0.9 against the manual set is the deployment
  success bar.

> The `fold_0.pt` / `cv_summary.json` currently tracked in git are a stale
> *undertrained* smoke artifact (a 5-epoch run). Regenerate with a real
> `--folds 1 --epochs 60 --patience 12` (to a plateau) or `--folds 5` run.

## DAPI / Phenocycler U-Net training data

`scripts/export_DAPI_training_pairs.groovy` writes one training tuple per
fluorescent QuPath entry (those with ≥ 1 `Follicle` annotation) into
`data_dapi/`, on the **same 2 µm/px grid** as the H&E set:

- `data_dapi/images/<donor>.png` — single-channel uint8 DAPI, normalised
  per slide and **inverted** so high nuclear density is dark, matching the
  polarity of the H&E hematoxylin images.
- `data_dapi/tissue/<donor>.png` — tissue mask (255 = tissue) from
  3-class Otsu on the raw DAPI histogram, lower threshold.
- `data_dapi/masks/<donor>.png` — binary Follicle mask, **derived from the
  region map** so the two cannot disagree.
- `data_dapi/regions/<donor>.png` — uint8 region label index:
  `0 Background, 1 RedPulp, 2 Trabeculae, 3 LargeVessel, 4 PALS, 5 Follicle`
  (painted in that order, later winning, so Follicle is painted last).
- `data_dapi/meta/<donor>.json` — the H&E sidecar keys plus
  `dapi_transform` (`otsu_bin`, `otsu_upper_bin`, `lo`, `hi`,
  `tissue_frac`, `polarity`), `dapi_channel_index`, `region_classes`, and
  per-class annotated area in µm².

`images/ + masks/ + tissue/ + meta/` is exactly the contract
`HE_unet_utils.discover_samples()` consumes, so the two cohorts can be
pooled with `--extra-data-dir`; `regions/` is the DAPI-only extra.
Gitignored, like `data/`. `scripts/make_synthetic_dapi_data.py` writes the
same layout (plus synthetic multiplex OME-TIFFs) for the tests and CI.

**Labels are machine-generated.** They come from the `SpleenRegions*` pixel
classifier, which used 25–33 marker channels; the model sees DAPI only. That
is deliberate (knowledge distillation), and it means Follicle vs PALS — which
the teacher separated with CD20 vs CD3e — is the irreducible error. See
`scripts/DAPI_unet_README.md`.

### `analysis/DAPI_unet_models/`

Output of `scripts/train_DAPI_unet.py`; same file names as the H&E models.
`cv_summary.json` additionally carries a `region_metrics` block (per-class
IoU + the Follicle↔PALS confusion) and `pooled_metrics.json` carries
`region_per_class_iou`. Checkpoints record `n_region_classes` /
`region_classes` / `region_class_weights` so inference knows how to slice
the `1 + K` channel head.

`analysis/DAPI_unet_predictions/` holds per-slide GeoJSON in **pixel**
coordinates (so they re-import into QuPath), unlike
`analysis/HE_unet_predictions/` which are in µm.

## Genotype metadata

### `Groups.xlsx`

Fluorescent-cohort metadata. Key columns the loaders read:

- `HANDEL ID` — canonical sample ID (e.g. `HDL011`, `1901HBMP004`).
- `ALT ID` — 4- or 5-digit short ID (e.g. `1901`, `19001`).
- `rs3184504` (or `rs3184504 (SH2B3)`) — the primary genotype column.
- `rs*` and `rs*_risk_het_prot` — 39 SNP genotypes + risk
  annotations for H7 and H9.
- `Age (yrs)`, `Gender`, `Ethnicity`, `C-pep (ng/ml)`, `HbA1c`, `BMI`
  — clinical covariates for H8.

### `Spleen_rs3184504_Genotypes.xlsx` (repo root)

H&E cohort genotypes. Key columns:

- `Sample ID` — canonical donor identifier. H&E image aliases are mapped to it
  explicitly in `config/sample_manifest.csv`.
- Some column containing `rs3184504` in its name — picked up by
  `next(c for c in cols if "rs3184504" in c)`.

## Sample ID conventions

The fluorescent cohort uses several aliases for the same donor:

| Form | Example | Where seen |
|---|---|---|
| HANDEL ID | `HDL073` | preferred canonical form, used by `data_utils` |
| ALT ID (4-digit) | `1902` | HiperGator filename prefixes |
| ALT ID (5-digit) | `19002` | `Groups.xlsx ALT ID` column |
| HBMP form | `1901HBMP004` | some HiperGator outputs |

`config/sample_manifest.csv` records every accepted alias and its evidence.
`data_utils.extract_sample_id(image_name)` performs an exact lookup through
`identity.resolve_sample_id`; an unknown, repeated, or conflicting alias is an
error. Filename parsing and ALT-ID digit slicing are not part of the workflow.

## Tissue-region classes

### QuPath classes

Defined in `classifiers/classes.json`:

- **Follicle** — B-cell zone of the white pulp.
- **PALS** — periarteriolar lymphoid sheath (T-cell zone).
- **RedPulp** — splenic red pulp.
- **Trabeculae** — fibrous trabeculae.
- **LargeVessel** — central arteries and large veins.
- **SmallVessel** — small vessels (capillaries, venules, sinuses).
- **Ignore*** — regions to exclude from analysis.

CODEX-cohort region names are harmonized into this vocabulary. Note the two
mappings differ deliberately — the annotation-CSV path
(`data_utils._CODEX_REGION_MAP`) and the region-layer path
(`build_region_layer_qupath.CLASS_TO_REGION`) disagree on
`Peripheral_White_Pulp`, and `tests/test_region_layer.py` asserts both:

| CODEX class | Annotation CSV | Region layer |
|---|---|---|
| `Red_Pulp` | `RedPulp` | (unused — RedPulp is the tissue remainder) |
| `Sinusoid` | `RedPulp` for vessel parenting; **area excluded** (it overlaps `Red_Pulp`, and summing both put annotated tissue at 132–149% of the field of view) | unused |
| `Trabecula` | `Trabeculae` | `Trabeculae` |
| `Peripheral_White_Pulp` | `PALS` | **`Perifollicular`** |
| `Vessel` | `SmallVessel` | — |
| `Vein`, `Arteriole` | `LargeVessel` | — |

`Peripheral_White_Pulp` is peri-follicular, not periarteriolar. Folding it into
PALS gave the CODEX donors a collar hugging every follicle that then consumed
the perifollicular band, and doubled their PALS area.

### Analysis regions — the region layer

`data_utils.REGION_ORDER` = **Follicle**, **PALS**, **Perifollicular**,
**RedPulp**, **Trabeculae**. This is what every notebook uses, and it is not
the same as the QuPath class list:

- **Perifollicular** has no QuPath class. It is a measured band around each
  follicle: `scripts/radial_follicle_profiles.py` walks outward in 30 µm bands
  from the follicle *boundary* (negative inside) and reads the width off the
  SMA fold-change peak, which sits at +15 µm in 12 of 14 donors. Per-donor
  width 45–105 µm, median 75. Unioned with `Peripheral_White_Pulp` where it
  exists; PALS wins where the two overlap.
- **RedPulp** is tissue minus everything else, where tissue comes from cell
  density (`data_utils.tissue_extent`) rather than from the RedPulp annotation.
  It reproduces QuPath's total annotated area at 0.92–1.12×, median 1.00.
- **SmallVessel/LargeVessel** are absent by design. Vessels come from
  `AnnotationsFinal.csv`, which carries their centroid and area.
- Objects outside every polygon are `"Unassigned"` — 4–9% of SmallVessels,
  sitting just outside the KDE tissue contour at section edges.

## Exclusion rules

Hard-coded in `data_utils.EXCLUDE_SAMPLES`:

- **HDL018, HDL021** — image-quality issues.
- **HDL172** — no genotype available.

Every loader applies these exclusions before returning. Rows without
a mappable genotype are also dropped.

## GeoJSON exports and the region layer

`scripts/export_all_regions_geojson.groovy` writes one GeoJSON
FeatureCollection per image into `analysis/geojson_regions/`, in **pixel**
coordinates, containing Follicle, PALS, RedPulp, Trabeculae (both spellings)
and Peripheral_White_Pulp. Read-only; give `--project` an absolute path or
QuPath fails with `dirBase is null`.

It supersedes `scripts/export_regions_geojson.groovy`, which filtered to
Follicle + PALS. That filter is why no Phenocycler donor had Trabeculae
geometry, and their trabecular area sat inside RedPulp instead.

`scripts/build_region_layer_qupath.py` turns those exports into
`analysis/regions/<image>.geojson` — **µm** coordinates, one feature per
polygon with `properties.region`. It reads `geojson_regions/` first and falls
back per file to `geojson/`, which is how the three ForSH2B3 CODEX sections
(absent from `project.qpproj`) still get geometry. Pixel size comes from
`data_utils.get_pixel_size`, never from the OME-XML: `1901HBMP004` and `HDL073`
declare 0.5078 in their headers and are really 0.374.

QC goes to `analysis/tables/region_layer_qc.csv`. Every donor's Follicle and
PALS area is checked against the annotation CSV on each run; 13 of 14 reproduce
at 0.95–1.00×, and HDL073 at 0.86× is real — it has many small follicles, so
the fragment filter removes more of its area than the cohort average.

All three directories are gitignored.

`scripts/predict_HE_unet.py` also writes GeoJSON (in µm coordinates)
into `analysis/HE_unet_predictions/` — one per H&E WSI. These can be
dropped back into QuPath via *Objects → Import objects* for visual QC.
