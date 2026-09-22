# The viewer's region layer, 2026-09-17

The maintainer reported that TissUUmaps draws the wrong `IsletNonEndocrine` regions. It does.
This is what was measured, which of the three candidate causes it is, and what was rebuilt.

Every number below was produced by running something against the files on disk, and the file each
one came from is named beside it.

## The join, stated first

The comparison that had to happen was between a feature in
`data/development/typing/app/json/<donor>.geojson` and a row of
`data/development/typing/geometry/endocrine_union_<donor>.parquet`. It is not the GeoJSON `id`:
the `Islet` features carry a QuPath UUID there and the `IsletNonEndocrine` features carry
`null`. The join is on **`properties.name` with its suffix removed**:

| class | `properties.name` | `islet_key` |
|---|---|---|
| `Islet` | `Islet_17` | `Islet_17` |
| `IsletExpanded` | `Islet_17_exp20um` | `Islet_17` |
| `IsletNonEndocrine` | `Islet_17_nonendo` | `Islet_17` |

`name.removesuffix("_nonendo")`, not a fixed character count. `"_nonendo"` is eight characters,
and `name[:-9]` turns `Islet_17_nonendo` into `Islet_1` — a join that still matches, still
produces a populated frame, and compares the wrong islets. That is what an earlier attempt did:
run against 6539 it matched 29 islets instead of 221 and reported that no non-endocrine region
touched its own islet at all, which is an artefact of the join and not a fact about the data.
With the suffix removed properly the join is exact and total: **221 of 221 on 6539, 586 of 586 on
6436**, and the same row counts as `islet_geometry.parquet`'s `core` rows for those donors (572
and 600 islets, of which those with an endocrine union get a non-endocrine region).

`scripts/check_build_tissuumaps_regions.py` now pins the suffix rule.

## The three tests

Run on **6539** and **6436**, per islet, against the shipped `app/json/<donor>.geojson`,
`geometry/endocrine_union_<donor>.parquet` and `app/islet_geometry.parquet` (rebuilt 2026-09-17,
20 donors, 18,404 core islets). Pixel size 0.5077663810243286 µm/px from
`core/phenocycler/data/qupath_manifest.json`; areas in px² times that squared.

### (a) Is it inverted? No.

If the shipped region were the endocrine union itself, every area would be complemented.

| ratio of the shipped `IsletNonEndocrine` area to | 6539 (median) | 6436 (median) |
|---|---|---|
| the endocrine union | 16.612 | 6.783 |
| core − endocrine union (what correction E specifies) | 15.666 | 8.005 |
| **(core ∪ 20 µm expansion) − endocrine union** | **1.000** | **1.000** |
| the drawn `IsletExpanded` annotation | 1.053 | 1.308 |

Fraction of islets within 2% of the third row: **1.000 on 6539 (221/221) and 1.000 on 6436
(586/586)**. Nothing is inverted — the region is a different region.

### (b) Is it scaled wrong? No.

- `json Islet area × (0.5077663810243286)²` over `islet_geometry` `core.region_um2`:
  **median 1.0000, 100% of islets within 2%**, on both donors.
- Region bbox against the image (`/data/tissuumaps/channels/<donor>/00_DAPI.tif`):
  6539 is 24960×23760 px and the regions span **94.2% × 85.5%** of it; 6436 is 29760×47520 px and
  the regions span **94.4% × 93.0%**. At the builder's nominal ×2 they would span 188% × 171% and
  189% × 186% — outside the image, which is why `infer_region_scale` falsifies the nominal and
  settles on ×1. Confirmed on the shipped copy: `/data/tissuumaps/regions/6539.geojson` has the
  identical bbox to its source, `1084.8, 1903.3, 23559.4, 20356.2`.

Worth noting for whoever touches the builder: `DEFAULT_UM_PER_PX = 0.2539` in
`scripts/senior/build_tissuumaps_project.py` is not this cohort's slide pixel size — level 0 of
these OME-TIFFs is 0.5078 µm/px, the same as the GeoJSON — so the nominal ×2 is wrong here and is
caught only because the geometry falsifies it. A slide whose tissue covered less than half the
frame would not falsify it, and the regions would come out at twice size. It is the constant that
is wrong, not the fallback.

### (c) Is it the wrong source? **Yes — this is the cause.**

The shipped region is the **pre-correction-E** non-endocrine region:
`(core ∪ 20 µm expansion) − endocrine union`, agreeing with that formula at ratio 1.000 for 100%
of islets on both donors (table above).

`core/phenocycler/scripts/typing/build_islet_geometry.py` corrected the definition on 2026-09-14
to `core.difference(endocrine)` — the part of the islet **core** that is not endocrine cells — and
says so in the file, with the same symptom recorded there that the measurements above reproduce.
The corrected outlines are in `data/development/typing/geometry/nonendocrine_<donor>.geojson`,
written 2026-09-17 13:23–13:48 for all 20 donors. They are right:

| check | 6539 | 6436 |
|---|---|---|
| `nonendocrine_<donor>.geojson` area over core − endocrine union | 1.0000 (221/221 within 2%) | 1.0000 (586/586) |
| the same over `islet_geometry` `non_endocrine.region_um2` | 1.0000 (221/221) | 1.0000 (586/586) |
| islets where the region is larger than its own core | 0 of 221 | 0 of 586 |

Against the shipped region: larger than its own core on **97.7% of islets (6539)** and **92.8%
(6436)**, and larger than the whole drawn `IsletExpanded` annotation on **67.9%** and **99.7%**.
The shipped file was written 2026-09-14 11:11 — after the correction landed in the code, before
the rebuild that produced its inputs. It is stale in *definition*, which is why the feature counts
looked fine: 600 / 598 / 586 on 6436 against exactly 600 core islets and exactly 586
`endocrine_union` rows. The counts were never the defect.

## A second defect in the same layer: `IsletExpanded`

The shipped `IsletExpanded` is the drawn QuPath `_exp20um` annotation, carried through unchanged.
It is not one shape, and it is not 20 µm everywhere. Measured over 40–60 sampled islets per donor
on the islet export each donor's geometry was built from:

- **Shape.** On 6539 it is a filled dilated islet — it contains the core on 57 of 60 sampled
  islets and none of the 60 has an interior ring. On 6436 and 6482 it is an annulus with the core
  punched out — 58 of 60 and 60 of 60 carry an interior ring, and the core's area lies inside it
  for only 7.1% and 3.4%. Two different regions under one class name.
- **Width.** The 95th-percentile distance from the drawn boundary to the core is 19.92–20.51 µm on
  19 donors and **27.72 µm on 6482**, and every one of them is named `_exp20um`. That is a
  different estimator from the 18.8–19.3 µm / 23.1 µm figures recorded in
  `build_islet_geometry.py` and `card_s3_compartments.R`, so the numbers are not comparable to
  those; what both find is a single donor drawn wider than the rest.
- **Area.** Drawn area over `core.buffer(20 µm)` runs 0.485 (6593) to 0.859 (6539) — mostly the
  annulus/disc difference, with 6482 high for its width.

Every expansion is 20 µm and derives from the signed distance, not the drawn annotation. So the
layer now carries the recomputed band.

## What was rebuilt

`scripts/senior/build_tissuumaps_regions.py` (new) writes
`data/development/typing/app/json/<donor>.geojson` from the sources the areas are measured from:

- **`Islet`** — the QuPath outline, unchanged, but resolved through
  `build_islet_geometry.islet_export`: recursive (so a `2026SeptNew/` re-export is seen), newest
  mtime wins, and a detection export is refused with its reason. `build_app_data.py:write_outlines`
  uses a flat `glob("*.geojson")` on the same directory, which cannot see a re-export at all.
  The export's own QuPath annotation UUID is carried through as the feature `id`, not regenerated,
  so a region in the viewer still traces back to the object in the project it came from.
- **`IsletExpanded`** — `core.buffer(20 µm / pixel_size).difference(core)`, the same geometry
  `build_islet_geometry.py` measures as the `expansion` region. A ring, oriented exterior-CCW and
  interior-CW because TissUUmaps fills with the nonzero winding rule
  (`glUtils._regionFillRule`, default `"never"` — outline only — and `"nonzero"` when filling is
  switched on), so the core is not filled in.
- **`IsletNonEndocrine`** — the features of `geometry/nonendocrine_<donor>.geojson` verbatim.
  Absent for a donor, the class is simply absent: that region is then Unavailable, which is what
  it was before correction E had a source, and is not a licence to draw something else.
- Anything else in the export is carried through untouched and counted in the log.

Coordinates stay in full-resolution pixels. Colours match the app's drill-down
(`segmentation_helpers.R`) so the same region is not two colours in two views.

**The build refuses to write a donor whose derived region disagrees with the area
`islet_geometry.parquet` records for it**, per islet, at a relative tolerance of 1e-6 — they are
the same computation on the same inputs, so anything larger means they have diverged. That check
is the thing the previous generation did not have. Run against the superseded 6539 file it
reports **0 of 221 islets agreeing, worst relative difference 175×**, and refuses to write it.

The ten superseded files are not deleted; they are renamed
`<donor>.geojson.superseded-2026-09-14` beside the new ones, an extension no glob in
`resolve_regions` matches.

## Per-donor counts

Built 2026-09-17 over the 20 donors that have an `endocrine_union_<donor>.parquet`. `agree` is the
number of islets whose derived region's area matches `islet_geometry.parquet` to within 1e-6,
per class.

| donor | Islet | IsletExpanded | IsletNonEndocrine | file | agree core / exp / non-endo |
|---|---|---|---|---|---|
| 115 | 1,304 | 1,304 | 1,293 | 191.6 MB | 1,304 / 1,304 / 1,293 |
| 6374 | 1,003 | 1,003 | 1,000 | 137.8 MB | 1,003 / 1,003 / 1,000 |
| 6380 | 67 | 67 | 66 | 10.9 MB | 67 / 67 / 66 |
| 6436 | 600 | 600 | 586 | 140.8 MB | 600 / 600 / 586 |
| 6442 | 490 | 490 | 486 | 81.5 MB | 490 / 490 / 486 |
| 6450 | 453 | 453 | 438 | 41.0 MB | 453 / 453 / 438 |
| 6461 | 933 | 933 | 926 | 152.2 MB | 933 / 933 / 926 |
| 6476 | 1,275 | 1,275 | 1,210 | 201.0 MB | 1,275 / 1,275 / 1,210 |
| 6482 | 489 | 489 | 488 | 35.4 MB | 489 / 489 / 488 |
| 6493 | 696 | 696 | 685 | 133.6 MB | 696 / 696 / 685 |
| 6505 | 3,123 | 3,123 | 2,974 | 294.7 MB | 3,123 / 3,123 / 2,974 |
| 6521 | 403 | 403 | 377 | 62.5 MB | 403 / 403 / 377 |
| 6523 | 364 | 364 | 359 | 40.9 MB | 364 / 364 / 359 |
| 6534 | 2,030 | 2,030 | 1,983 | 211.7 MB | 2,030 / 2,030 / 1,983 |
| 6539 | 572 | 572 | 221 | 48.6 MB | 572 / 572 / 221 |
| 6551 | 427 | 427 | 413 | 67.0 MB | 427 / 427 / 413 |
| 6566 | 881 | 881 | 870 | 86.8 MB | 881 / 881 / 870 |
| 6591 | 1,172 | 1,172 | 1,132 | 131.1 MB | 1,172 / 1,172 / 1,132 |
| 6593 | 420 | 420 | 416 | 92.6 MB | 420 / 420 / 416 |
| 6623 | 1,702 | 1,702 | 1,630 | 139.8 MB | 1,702 / 1,702 / 1,630 |
| **total** | **18,404** | **18,404** | **17,553** | 2.3 GB | every islet, every class |

18,404 is exactly the core-islet count in the rebuilt `islet_geometry.parquet`, and every donor
now has a region layer — ten did not before. `--check-only` measures and writes nothing.

**Where the non-endocrine count is short of the islet count**, the islet has no endocrine union:
`build_nonendocrine_outline.py` writes a row only where at least one endocrine cell in that islet
has an outline in `data/redsea_scratch/geojson/`. Those islets get no non-endocrine region, and
`islet_geometry.parquet` carries `NaN` for their `non_endocrine` area — unavailable, not zero.
6539 is the extreme case at 221 of 572, which is consistent with its islets being mostly one-cell
fragments; do not read its 351 missing regions as islets without non-endocrine tissue.

Verified after the rebuild, on four donors, against `islet_geometry.parquet` and
`endocrine_union_<donor>.parquet`:

| donor | non-endo / (core − endocrine) | Islet / core_um2 | IsletExpanded / expansion_um2 | islets larger than their core | annuli |
|---|---|---|---|---|---|
| 6380 | 1.0000 (66/66 within 2%) | 1.0000 (67/67) | 1.0000 (67/67) | 0 | 67 of 67 |
| 6436 | 1.0000 (586/586) | 1.0000 (600/600) | 1.0000 (600/600) | 0 | 600 of 600 |
| 6482 | 1.0000 (488/488) | 1.0000 (489/489) | 1.0000 (489/489) | 0 | 489 of 489 |
| 6539 | 1.0000 (221/221) | 1.0000 (572/572) | 1.0000 (572/572) | 0 | 572 of 572 |

Before the rebuild the first column read 15.666 on 6539 and 8.005 on 6436, with 0 of 221 and 0 of
586 islets within 2%, and the expansion was an annulus on some donors and a filled disc on others.

## Still open, for whoever owns those files

1. `core/phenocycler/scripts/typing/build_app_data.py:write_outlines` is a second producer of
   `app/json/<donor>.geojson`. It already merges `geometry/nonendocrine_<donor>.geojson`, so it
   would no longer ship the superseded non-endocrine region — but it still carries the drawn
   `_exp20um` annotation through and still resolves the islet export with a flat glob. Running it
   would undo both fixes above. It should call this builder rather than reimplement it.
2. `app/gson/<donor>.geojson.gz` is the app's fallback when the `.geojson` is missing. Those
   copies are still 2026-09-14 and still carry the superseded region. They are not mine to write.
3. `/data/tissuumaps/regions/<donor>.geojson` is the copy the viewer actually serves, taken at
   project-build time. It was refreshed from the rebuilt json for the ten donors that have a
   `.tmap` project, using `build_tissuumaps_project.py`'s own `resolve_regions` / `_load_geojson` /
   `copy_regions` and its own `infer_region_scale` against each image's real dimensions — every
   one settled on ×1, with the nominal ×2 falsified by the geometry as before. The superseded
   copies are renamed `<donor>.geojson.superseded-2026-09-14` in place — **except 6436's**, where
   the refresh was run twice and the second run renamed the first run's corrected output over it.
   Nothing unique was lost: that file was only `app/json/6436.geojson` plus a stamped
   `collectionIndex`, and the 2026-09-14 original is still at
   `data/development/typing/app/json/6436.geojson.superseded-2026-09-14` (84,742,877 bytes).
   The other ten donors have
   no project yet; `build_all_tissuumaps_projects.py --regions data/development/typing/app/json`
   is what would give them one, and that is a channel build, not a region one.
4. `DEFAULT_UM_PER_PX = 0.2539` in `build_tissuumaps_project.py`, as above.
5. **The app's drill-down never draws the non-endocrine region at all**, whatever this layer
   contains. `apps/islet/R/segmentation_helpers.R` selects it
   (`.keep_selected(polygons$IsletNonEndocrine, paste0(islet_key, "_nonendo"))`), gives it a
   display name (`Non-endocrine`) and a colour (`#8C6BB1`) — and then builds the plot from
   `layer_order <- c("IsletExpanded", "Islet", "Lymphatic", "Capillary", "Nerve")`, which does not
   name it. It is never added to `poly_list`, so it is never in `combined_sf` and never reaches
   `geom_sf`. TissUUmaps draws it; the Shiny drill-down does not. Adding `"IsletNonEndocrine"` to
   that vector is the whole change, but that file is not mine to edit.

A sanity note on the reader for whoever makes that change: class matching goes through
`get_classification_name()`, which parses the `classification` JSON and takes `name`, so the
colours this builder writes do not affect which features land in which layer.

Nothing here is a release and the `NO_VERIFIED_TYPING` latch has not moved. The regions are a
display of the geometry; they are not a result.
