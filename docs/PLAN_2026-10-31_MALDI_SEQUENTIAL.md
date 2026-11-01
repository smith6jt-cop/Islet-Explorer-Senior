# Plan: add MALDI glycan and glycogen imaging to the IO60 sequential-slide analysis (MultiPlex5000)

Status: written 2026-10-31 from the maintainer's answers. **Implemented 2026-10-31** on branch
`maldi/serial-platform` (not yet committed). Waiting on the maintainer's inspection rows
(`maldi_tissue_inspection.csv`) before MALDI Ingest can run on the real export.

**Where the implementation departs from the text below, and why:**
- `maldi_acquisitions.csv` is a Source Inventory **output**, not a tracked generated file. A tracked
  copy of facts derived from the export would be a second copy that can go stale. The tracked
  record is only the maintainer's `maldi_tissue_inspection.csv`, pre-filled with the 32 IO60 Tissue_ids.
- **Tissue mask:** a spot is on the tissue if it is inside its raster and at least one ion was
  detected in it (`maldi.ON_TISSUE_METHOD`). The ion images showed each raster drawn to its tissue's
  outline, with 0–2.1% of spots per tissue holding no ion at all, so no per-image background fit
  was needed.
- **Images:** MALDI Ingest writes the summed-ion image and mask as `images/<Tissue_id>.npz` on the
  spot grid, which is what registration reads. An OME-Zarr export for the viewer belongs to the
  app step.
- **Registration to a Xenium neighbour is not built.** A section recorded next to a Xenium slide
  (expected for 6457 only) is reported unavailable with that reason.
- **Fixed image:** a PhenoCycler pancreas scan is registered by its `Tissue` or `pLN`
  annotation (whichever matches the inspected tissue); a dedicated pLN slide by its DAPI Otsu
  mask.
- The registration and integration routes run through the new routers
  `stages/registration_route.py` and `stages/integration_route.py`. A project naming no route keeps
  the same-slide pilot.
- **New site source:** `typing_export_root` (`results/v3/export`: `typing_islets/`, `structures/`),
  alongside `maldi_spots`.

## Context

Preliminary MALDI imaging mass spectrometry data arrived as two CSVs in the untracked repo folder `MALDI/`
(from `~/Downloads/wilkes_pancreas_glyc.zip`). The maintainer wants it added to MultiPlex5000's **sequential-slide**
analysis next to PhenoCycler and Xenium.

**What the files hold** (read-only profile, 2026-10-31):
- `wilkes_pancreas_glycan_pbp.csv` (433 MB) holds 38 sodiated N-glycan compositions (`Hex5HexNAc4NeuAc1..Na`, …).
- `wilkes_pancreas_glycogen_pbp.csv` (218 MB) holds 10 glycogen m/z features (`X689` … `X2147`).
- Both have the **same 2,288,540 spots**. `spotId`, `x`, `y`, `Tissue_id` and `raster` are identical, with 0 mismatches on a full join. So this is one acquisition with two feature sets.
- **Grid:** spots sit on a 50 µm grid. `x`/`y` are µm in one stage frame shared by all 60 tissues.
  - There are 82 rasters, and `spotId` is unique across the file.
  - `z` is constant and `Date` is NA throughout.
- **Values** are peak-picked: 0, or an intensity (mostly integers).
  - Example: Hex5HexNAc4NeuAc1 is 0 in 68% of spots.
  - 1.1% of spots are 0 on every glycan.
  - A spot has a median of 27 of the 38 glycans detected.
- **Tissues:** 18 IO60 donors, 6426 (nPOD ND, not IO60) and 27 HPAP sections.
  - The 18 IO60 donors are 8 ND, 1 Aab+, 9 T1D (`donor_metadata_panc.xlsx`). HDL115 is donor 115 by the cohort profile's alias.
  - 8 IO60 PhenoCycler donors have no MALDI: 0112, 6347, 6450, 6505, 6533, 6534, 6580, 6623.
- **Reference:** `literature/MultimodalImagingPipeline.pdf` (Pham et al., bioRxiv 2025.04.07.647675, posted 2026-02-14).
  - It describes MALDI IMS with CODEX on the **same** section, registered with image2image. That is method context only, since ours are serial sections.
  - The maintainer's fork `smith6jt-cop/image2image`, with `image2image-reg` and `imzy`, sits in the `image2image` conda env. Its editable install points to a removed folder, so nothing here depends on it.
- Nothing in the repo mentioned MALDI before this plan.

**Maintainer's decisions (2026-10-31):**
1. **The registration images are the ion images in these files.** Each feature rasterised on the 50 µm grid is an image of the section. No separate optical scan is assumed.
2. **The MALDI sections are sequential slides from the same series of slides** as the PhenoCycler and Xenium data.
   - Each one is a serial section: different cells, never a same-cell claim.
   - Which tissue each Tissue_id is (pancreas or pLN/lymphoid), and which slide of the series it neighbours, is recorded from what is seen, never parsed from the ID.
   - This follows the Xenium pLN/lymphoid inspection precedent ("settled only by looking at the slide").
3. **A 0 is "not detected", kept as a measured value.** It counts in denominators. Off-tissue spots are removed by a tissue mask, never by value.
4. **Scope: the 18 IO60 donors only.**
   - 6426 and HPAP are counted in the inventory as out of cohort and are not ingested.
   - 6457's exclusion covers its PhenoCycler data only (`cohort.EXCLUDED_PLATFORMS`), so its MALDI counts. Its PhenoCycler section is not available as a neighbour; its Xenium is.

**Rules carried through:**
- Serial sections never transfer identity or measurements onto cells (`docs/INTEGRATION.md:241-252`, `scripts/restain/README.md:460-463`).
- No expansion rings or distance bands; continuous signed distance only.
- Missing is unavailable, not zero.
- No result is a verification target.
- No analysis beyond what is asked: this plan stops at the integration outputs.
- Steps are named, never numbered, in prose.

**What the repo does not record:**
- The PhenoCycler block for any donor; the `block_id` column in `integration/manifest.py:78` is empty everywhere.
- The block of the dedicated pLN slides.
- Xenium blocks are only on the gitignored run sheet `Xenium_Visium request.xlsx`.
- These IDs cannot be resolved from the ID text:
  - 6476, 6566, 6593 and HDL115 carry no block;
  - 6380_01, 6442_03, 6482_02, 6521_02, 6523_01 and 6551_03 carry a bare block number with no `_pln`.

  The inspection below settles them.

---

## Source records (not in Git)

- Copy both CSVs, read-only, to `/data/io60-current/sources/maldi/wilkes_2026-10/` with a `SHA256SUMS` file.
- Declare them to the site as one source, `maldi_spots` (both paths and their sha256).
- Add `MALDI/` to `.gitignore`, so the 622 MB of raw data is never committed. The folder itself stays until the maintainer moves or deletes it.

## MALDI as a platform (branch `maldi/serial-platform`, PR)

**Platform registration.** There is no single registry, so every place that names platforms changes:
- `multiplex5000/cohort.py`: `PLATFORMS += ("maldi",)`.
  - `excluded_on()` then accepts the name and still refuses misspellings.
  - 6457 stays PhenoCycler-only.
- `profiles/platforms/maldi.yaml` (schema `profiles.py:59-66`):
  - `id: maldi`, `display_name: MALDI imaging (N-glycan + glycogen)`;
  - `pixel_size_um: 50`, the spot pitch;
  - `panel_names`, and `sources` (the acquiring lab and the `maldi_spots` source).
- `multiplex5000/panel_names.py`: a namer that builds the panel name from the feature columns, e.g. `38 N-glycans + 10 glycogen`. The maintainer confirms the wording at the Source Inventory gate.
- `stages/inventory.py`: entries in `PANEL_PROFILES`, `PLATFORM_OF`/`PLATFORM_NAMES`/`ACQUISITION_COLUMNS` and `FILLS`; and `_blocks`.
- `integration/contract.py` `MODALITIES`: change only if the integration contract is reused by the Integration step.

**Tracked records** (`core/phenocycler/data/integration/`):
- `maldi_acquisitions.csv`, generated from the source by `core/multiplex5000/scripts/maldi/inventory_maldi.py`. It holds literal facts only:
  - `tissue_id`;
  - `donor_id`, through the cohort aliases (HDL115 → 115);
  - `id_block_token`, the ID's text after the donor, kept as text and never interpreted;
  - `rasters`, `n_spots`, the bounding box, and the source sha.
- `maldi_tissue_inspection.csv`, filled by the maintainer from the inspection figure. Columns:
  - `tissue_id`;
  - `tissue` (pancreas or pLN/lymphoid);
  - `neighbour_slide`: the PhenoCycler pancreas, PhenoCycler pLN or Xenium slide it was cut next to;
  - `orientation_note`.

**Relationships.** `multiplex5000/relationships.py` reads the inspection file as a new source.
- Kind: `serial_section` (`:45-51`).
- Note: "sequential slides from the same series (maintainer, 2026-10-31)".
- `python -m multiplex5000.relationships` prints it. A disagreement with the donor notes blocks the gate, as with the other sources.

**Tests** (`tests/current/test_maldi_*.py`), on a small synthetic MALDI CSV:
- `maldi` is accepted and a misspelling is refused.
- 6457's MALDI counts, and 6426/HPAP are refused from the cohort.
- The panel name is built from the columns.
- A missing inspection row blocks the gate.

## Source Inventory: MALDI rows and the inspection gate

- **Acquisitions:** one per IO60 Tissue_id, with `acquisition_id = maldi:<tissue_id>`. The two feature files make one acquisition, so one dataset.
- **Out of cohort:** 6426 and the HPAP sections appear in the gate as counts, not as rows.
- **Coverage figure:** a **MALDI slide** group beside the PhenoCycler, Xenium and dedicated-pLN groups.
  - Each cell is filled by its inspected tissue.
  - Each carries a thumbnail: the summed-ion image on the 50 µm grid.
- **Inspection figure** (new, beside `_inspection_figure`). For each IO60 Tissue_id it shows:
  - the summed-ion image and two contrasting single-feature images (e.g. a high-mannose glycan and a glycogen peak);
  - beside them, the donor's PhenoCycler pancreas overview with its Tissue and pLN annotations, the dedicated pLN slide, and the Xenium overview.

  Every image is shown as scanned.
- **Blocking check:** the gate blocks until every IO60 MALDI tissue has a `maldi_tissue_inspection.csv` row.

## MALDI Ingest (new stage; code id `04b_maldi_ingest`, after RNA Ingest in `stages/__init__.py` PLAN/BUILT)

- **Module:** `stages/maldi_ingest.py`, patterned on `stages/rna_ingest.py`.
- **Inputs:** the `maldi_spots` source and the accepted Source Inventory.
- **Reading:** DuckDB reads both CSVs joined on `spotId`, keeping IO60 tissues only.
- **Keys:** `spot_uid = identity.cell_uid("maldi:<tissue_id>", "spot50", spotId)` (the existing function), with `check_cell_keys` for uniqueness.
- **Per-tissue output:** `data/spots/tissue_id=*/data.parquet` with `spot_uid`, `x_um`, `y_um`, `raster`, the 48 features as measured, and `on_tissue`.
- **Zeros:** 0 is kept as "not detected". A NaN would mean the value was missing from the file; none is expected, and any is counted.
- **Tissue mask (`on_tissue`):**
  - Computed per tissue on the summed-ion image, from that tissue's own background (the low peak of the summed-signal histogram, or a mixture gate), chosen as the typing lane chooses its level methods.
  - Spots at the raster edge with no tissue signal are off tissue.
  - The method is named in the receipt. It is a property of each image, never a value cutoff on a feature.
- **Images:** per-feature images on an `integration/rasterize.py:45` `GridSpec` (50 µm pixels, µm origin), written as OME-Zarr for the viewer and for registration.
- **Normalization:** none. Intensities stay as measured; MALDI intensity normalization (TIC or other) is a later, separate decision.
- **Gate:** automatic, as for the measurement stages. Its non-blocking indicators:
  - spots per tissue on and off the mask;
  - detection fraction per feature per tissue;
  - all-zero spots;
  - the mask figure per tissue.

## Serial Registration: MALDI onto the PhenoCycler frame (Registration stage, new serial MALDI route)

`stages/same_slide_registration.py` handles same-slide pairs only.
- **Route:** add one selected by `settings.registration_route: maldi_serial`, following the routes in `stages/typing.py`, in a new module `stages/maldi_serial_registration.py`.
- **Fixed image:** the inspected neighbour at 50 µm.
  - Normally the PhenoCycler image's tissue mask and DAPI, downsampled from the qptiff level the manifest names.
  - For 6457, and for any section recorded next to a Xenium slide, the Xenium DAPI. That section then reaches the PhenoCycler frame through the existing Xenium→PhenoCycler `Transform`, where one is recorded.
- **Method:** reuse `integration/register.py`:
  - `coarse_orientation`: a rotation and mirror search, since serial sections can be flipped;
  - `refine_ecc` on the summed-ion image and its tissue mask;
  - an affine only, with `refine_bspline` off unless a section shows deformation the affine cannot follow. If B-spline is used, it stays capped by `reg_max_disp_um`.
  - The result is an `integration/transform.py` `Transform`: moving = MALDI stage µm, fixed = PhenoCycler µm.
- **QC,** measured and reported, never tuned to a result:
  - tissue-mask Dice;
  - outline nearest-neighbour RMSE (µm);
  - scale and shear within physical bounds (the re-stain lane's "scale stays physical");
  - where the neighbour has QuPath objects, the overlap of each object with a feature that visibly follows it, reported as evidence only.
- **Gate: a review stop.**
  - Per tissue: an overlay of the MALDI ion images on the neighbour's DAPI and object outlines, with the QC table.
  - The maintainer accepts the transforms, because a registration cannot be verified by its own metrics.

## Serial Integration: spots to objects, never to cells (Integration stage, serial MALDI route)

- **Route:** `stages/maldi_serial_integration.py`, under the same `maldi_serial` route setting.
- Each on-tissue spot is mapped through the accepted transform into the PhenoCycler frame.

**Per spot:**
- `signed_dist_islet_um`: the continuous signed distance to the nearest islet object on the neighbour section. It is computed as the existing islet distances are, from the neighbour's islet objects. No bands.
- The containing PhenoCycler object, by spot centre: islet (`Islet_N`), Nerve, Blood vessel, Lymphatic, Tissue, or pLN annotation.
- The spot's own registration uncertainty, carried from the QC RMSE.

**Per object** (for each MALDI tissue and neighbour object):
- the number of spot centres inside;
- the detection fraction per feature (a 0, meaning not detected, counts);
- the median intensity among detected spots;
- `evaluable`: true only if the object holds at least one whole spot footprint **after** allowing for the measured registration RMSE.
  - An object that fails is `unavailable`, with the reason recorded, and never silently dropped.
  - At 50 µm many small islets will be unavailable. That is reported, not worked around.

**Outputs:**
- Keyed by `spot_uid` and by `(donor_id, neighbour acquisition, object id)`.
- Every table carries `relationship = serial_section`, and no column names a cell.

**Gate:** automatic, with integrity checks plus a figure per tissue (spots coloured by containing object, and the distance field). Biological review happens where the results are used.

## Next, outline only

- **App and viewer:** a MALDI layer in the run20 domain config (ion images and per-object tables), once the Integration outputs exist.
- **MALDI intensity normalization:** TIC or reference-based, decided separately.
- **Cohort Analysis:** comparisons by disease status read these outputs, but only when asked.

---

## Critical files

**MultiPlex5000** (`core/multiplex5000/`):
- `multiplex5000/cohort.py`, `profiles/platforms/maldi.yaml` (new), `multiplex5000/panel_names.py`;
- `stages/inventory.py`, `stages/__init__.py`;
- new stage modules `stages/maldi_ingest.py`, `stages/maldi_serial_registration.py`, `stages/maldi_serial_integration.py`;
- `multiplex5000/relationships.py`, new `scripts/maldi/inventory_maldi.py`;
- `projects/io60.yaml` (the MALDI stages and settings, as each route lands), the site files (`maldi_spots`);
- `core/phenocycler/data/integration/{maldi_acquisitions,maldi_tissue_inspection}.csv`;
- `docs/INTEGRATION.md` (a MALDI serial section), repo `CLAUDE.md` (one paragraph), `.gitignore`.

**Reuse:**
- `identity.cell_uid` / `check_cell_keys`;
- `integration/rasterize.py` `GridSpec` / `tissue_mask`;
- `integration/register.py` `coarse_orientation` / `refine_ecc` / `refine_bspline` / QC;
- `integration/transform.py` `Transform`;
- the existing islet signed-distance code;
- `tasks.run_stage`, `gates`, the inventory figure and thumbnail helpers.

## Verification

- `make test-fast` and `make test-lanes` pass. The new synthetic tests cover:
  - the platform (accepted, misspelling refused, 6457 counts);
  - ingest keys and zeros;
  - the tissue mask on a synthetic grid with glass;
  - registration recovering a known rotation, mirror and shift on a synthetic MALDI/DAPI pair;
  - integration on synthetic objects, including one below the evaluable size.
- `make test-workflow`: the graph and the notebooks give the same outputs for the MALDI stages.
- **Pilot:** two IO60 tissues whose inspection rows the maintainer has filled. They are picked by data completeness, never by disease status. Then all 18 donors.
- **Ingest:** spots in = CSV rows for IO60 tissues, and the two feature files join with zero mismatches (already 0 on the full files).
- **Registration:** the maintainer accepts the overlays, and Dice, RMSE and scale are reported per tissue.
- **Relationships:** `python -m multiplex5000.relationships` reports the MALDI serial relationships with no unexplained disagreement.

## Asked when they arise

- The panel-name wording, at the Source Inventory gate.
- Each tissue's inspection row: its tissue and neighbour slide.
- Registration acceptance.
- MALDI intensity normalization, later and separately.
