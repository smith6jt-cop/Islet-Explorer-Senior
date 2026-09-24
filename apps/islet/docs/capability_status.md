# Islet Explorer (Senior) — App Capability Status

> **HOW TO READ THIS FILE.** It has two sections, each a dated snapshot of a *running app*, and
> they describe two different datasets, selected by `ISLET_DOMAIN_CONFIG`. A tab verdict here is
> what the app did on the date above it, in a headless browser — it is not a standing claim, and
> nothing in it has been re-run by editing this document. Where a number has been re-measured, the
> date of that measurement is given inline. Where it has not, the original date stands and the
> sentence says so.
>
> - **Run cohort** (`config/run20_default.yml`) — tab verdicts assessed **2026-09-14**; the cohort
>   and data-on-disk figures re-measured **2026-09-18**.
> - **Legacy CODEX cohort** (`config/islet_default.yml`) — assessed **2026-06-29**, with two
>   corrections made inline on 2026-09-18 and dated there.
>
> The earlier version of this note flagged three passages as "stale" without correcting them. Each
> is now corrected in place instead.
>
> **What this file does NOT cover:** the typing itself. A tab verdict here says a card drew; it
> says nothing about how the cells under it were named. For that, read
> `core/phenocycler/scripts/typing/RECONSTRUCTION.md` and the header of
> `apps/islet/config/run20_default.yml`. The one typing fact a reader of this file needs is below,
> under "What the run cohort's cells were typed by".

## Run cohort (`config/run20_default.yml`) — assessed 2026-09-14

**Method:** `tests/tab_smoke.R` — Chrome via shinytest2, reporting per tab the R output errors
Shiny surfaced, the JS console errors, and how many outputs exist but drew nothing. Screenshots are
what the verdicts rest on, not a code read.

**The table below covers five tabs because the app had five tabs mounted on 2026-09-14.** It now
has nine `nav_panel_hidden` panels and `tab_smoke.R` drives eight of them (`tab_smoke.R:14`):
Neighborhoods, CrossModal and LymphNode are not in this table and have **not** been assessed here.
Their absence is a gap in this document, not a verdict on those tabs.

| Tab | R errors | JS errors | Verdict |
|-----|---------|-----------|---------|
| **Plot** | 0 | 0 | Works — islet-size scatter and donor-group violin render on ND/T1D; the domain's denominator note shows in the sidebar |
| **Trajectory** | 0 | 0 | Works — composition pseudotime, both UMAPs, feature selector carrying the eight class names |
| **Spatial** | 0 | 0 | Works — tissue scatter 122,951 foreground / 729,312 background for donor 115; neighbourhood Cards A–C live on 7,634 islets |
| **Statistics** | 0 | 0 | Works — narrative and controls render; result panels populate on demand after Run Statistics |
| **Viewer** | 0 | 1 | Wired correctly, image fetch fails **as expected**: Avivator mounts, the dropdown lists `115.ome.tiff`, the URL builds, and only the tile request fails. `runApp` cannot serve HTTP Range. Verify under the deployed server (`deploy/local/README.md`). |

**Data on disk for this cohort — re-measured 2026-09-18.** The cohort moved from 10 donors to 20,
and the trees did not all move with it. Stated separately, because the difference is the point:

| tree | donors | measured |
|---|---|---|
| `parquet/islets.parquet` | **20** | 18,384 islets, 104 columns; 8 ND / 3 Aab+ / 9 T1D |
| `cells/` (drill-down per-islet CSVs) | **10** | 7,780 files, from 115, 6374, 6436, 6442, 6461, 6476, 6482, 6493, 6539, 6593 |
| `donors/` | **10** | |
| `gson/` | **10** | |

The four trees at 10 donors are the pre-2026-09-15 cohort; `islets.parquet` was rebuilt on the 20.
A drill-down click on an islet belonging to one of the ten donors the `cells/` tree does not carry
therefore has nothing to open. That is a build gap, not a measurement.

The old line here read "7,780 islets, 13.37 M cells, 10 donors (8 ND / 2 T1D, no Aab+)". The islet
and donor counts are superseded above; the "no Aab+" is simply no longer true — 3,966 of the 18,384
islets are Aab+, across 6450, 6505 and 6521.

DuckDB views `islets`, `tissue`, `cell_distance`, `uns_targets`, `uns_composition`; `uns_markers`
deliberately absent because this dataset has no per-marker positivity.

**What the run cohort's cells were typed by — measured 2026-09-18.** The app data is built from
the **subpanel classification**, not from scVI: context from the QuPath objects, then one question
per context answered by the markers that can answer it, against reference levels drawn from the
ordinary tissue and fitted **per donor**. `config/run20_default.yml` declares
`classifier_label: subpanel classification` and the Plot tab's measure-group headings read it, so
no heading in this app says "(scVI)" over this dataset.

The source is `data/development/typing/cells/` — **20 donor partitions, 25,230,109 rows, rewritten
2026-09-19 with the islet chain** (measured 2026-09-24). Its `cell_class` column carries eight
classes plus `Excluded` and nothing else: Exocrine 20,099,247 · Immune 2,875,751 · Blood vessel
1,422,676 · Alpha 294,514 · Beta 235,052 · Excluded 165,431 · Delta 99,015 · Nerve 32,836 ·
Lymphatic 5,587. Three things follow, and all three are easy to get wrong from older text:

- **`typing/subpanel/subpanel_calls.parquet` is a PREVIOUS generation** — 25,150,537 rows, written
  Sep 11 — and nothing in this app reads it. That tree stacks generations side by side with
  nothing marking which is which; check mtimes before citing anything in it.
- **The 2026-09-18 islet chain is on disk, but not yet in everything this app draws.** Every
  `typing/cells` partition carries `islet_route` (cd99 408,950 · elongation 29,060 · hormone
  190,571). The `typing/app` artifacts were written at different times: islet geometry, islet
  composition, neighbourhood metrics, islet marker means and islet Leiden on 2026-09-21, after the
  rewrite; `parquet/islets.parquet` on 2026-09-17 and the per-cell tissue, drill-down, outline and
  distance artifacts on 2026-09-14, before it. `make check` names each stale one, and the DATA step
  rebuilds them.
- **The ninth vocabulary entry, `Unresolved`, is carried by no cell yet.** It is an *absence*, not a
  class — what the refinement pass (`core/phenocycler/scripts/typing/refine.py`, producer label
  `Unresolved (no lineage evidence)`) writes on a cell demoted out of Exocrine. The pass was ruled in
  on 2026-09-23 and is wired as the `build_cells.py --refine` join, which adds `refined_class`
  beside `cell_class`; the join has not run on disk, `build_cells.py` still sets `phenotype` from
  `cell_class`, and `parquet/islets.parquet` carries `count_*` / `prop_*` / `peri_*` /
  `enrich_z_*` / `min_dist_*` over the eight classes with **no `Unresolved` column of any kind**.
  Its 5th-percentile bound demotes about 5 % of every class by construction; whether any builder
  should serve `refined_class` instead of `cell_class` is open.

`Exocrine` is one label here and is not divisible: the 60-plex carries no acinar-specific and no
ductal-specific marker. `Immune` is one label for the same reason — the immune second pass
(`scripts/typing/immune_subtype.py`) never modifies `cell_class`, and **no `typing/app` builder
reads `data/development/typing/immune`** (the `grid_composition` layer does, for `fine.csv`);
coverage comes from each builder's own `COHORT.resolve`.

**OME-TIFFs, measured 2026-09-18:** 16 files in `/data/islet_ome_tiff` — 115, 6374, 6380, 6436,
6442, 6450, 6457, 6461, 6476, 6482, 6493, 6521, 6539, 6551, 6579, 6593. Against the 20-donor
cohort that is **14 of 20**: the six batch-2 slides 6505, 6523, 6534, 6566, 6591 and 6623 are not
converted. Two of the 16 (6457, 6579) are not in the cohort at all — 6457 is centrally excluded,
and the Viewer's dropdown does not filter on the cohort (see `deploy/local/README.md`).

**Eight cohort facts had been written into the UI as literals** and printed as measurements of
whichever dataset was loaded — a donor id (6521) that the run cohort did not contain when those
strings were written, "N=15 total", "N=5 per group", "all fifteen", a hard-coded
`c("ND","Aab+","T1D")` filter, and "(scVI)" on a dataset with no scVI. (6521 has been in the run
cohort since 2026-09-15, as one of the three MARGINAL donors, which is exactly why a cohort fact
belongs nowhere but the loaded data: the literal was wrong when written and would be wrong the
other way now.) All now derive from the loaded data or the DOMAIN. **None of them was visible to the testthat
suite**, because every path resolved and every column existed — the suite was 238 assertions when
that was written and is 737 over 69 `test_that` blocks as measured on 2026-09-18, and the size of
it was never what would have caught these.

---

## Legacy CODEX cohort (`config/islet_default.yml`)

**Assessed:** 2026-06-29 · **Method:** the dev server was booted
(`Rscript -e 'shiny::runApp(".", port=7777)'`, `RETICULATE_PYTHON` → `scvi-env`) and every tab was
loaded and screenshotted in a headless browser. Verdicts below are what the **running app actually
did**, not a code read. Screenshots are in `docs/img/capability/`.

## Headline

- The app **boots cleanly** and **all 5 data tabs are functional** on the data currently on disk
  (Viewer needs images — see below).
- **Trajectory works** — an earlier "broken" call was a **premature 8s screenshot**. The tab was
  merely **slow to first render** (~10–30s): `traj_data_clean` ran a 9,130-iteration
  `get_islet_annotations` loop, each scanning the full 9.6k-row spatial lookup (~88M comparisons per
  render). **Fixed 2026-06-29** by vectorizing that lookup into one `match()` (renders instantly).
  The "islet-key join" was a red herring — the keys (`Islet_755` both sides) match perfectly.
- **Viewer serves per-donor OME-TIFFs** (Jul 2026). **Corrected 2026-09-18: not all 22.** There are
  16 files in `/data/islet_ome_tiff`; the six batch-2 slides 6505, 6523, 6534, 6566, 6591 and 6623
  are not converted. The per-donor list is in the run-cohort section above. The `.qptiff`
  slides are converted to pyramidal, tiled OME-TIFF via `scripts/senior/convert_viewer_ometiff.sh`
  (bioformats2raw 0.12.1 → raw2ometiff 0.10.0, Baseline series only, LZW) into
  `/data/islet_ome_tiff/<case_id>.ome.tiff` (+ `<case_id>.offsets.json`), surfaced through the
  `www/local_images` symlink. Donor **6539 verified end-to-end** in Avivator (islets render;
  INS/GCG/SST/DAPI defaults correct; 2 mm scale bar; 5-level pyramid) —
  `/data/islet_ome_tiff/verification_6539.png`. All 22 convert via the resumable batch
  (`… all22`). The cohort **mixes panels** (batch-1 SST@35 / batch-2 SST@1, some 58-ch) so the
  channel_config was made **name-based** (`build_channel_config_b64` emits INS/GCG/SST/DAPI by name;
  Avivator resolves them against each image's own OME-XML) — correct defaults for every donor.
  **Avivator requires HTTP Range (206)** — production nginx/shiny-server provides it; the dev
  `runApp` (httpuv) does **not**, so images render only under the deployed server, not `shiny::runApp`.
- **Every tab still serves the OLD provisional fine phenotypes** (Beta cell, Alpha cell, M2
  macrophage, NK cell, Regulatory T cell, Stromal, …) — true of **this domain**, the CODEX one,
  which is the only dataset this 2026-06-29 assessment covers. **Corrected 2026-09-18:** the
  sentence that followed it, that the new REDSEA/RESTORE/broad-lineage outputs are "intentionally
  not wired into the app" and live only in
  `notebooks/senior_phenotyping_redsea_restore.ipynb`, no longer describes anything. RESTORE is
  retired for the source workflow (2026-08-14) and the flat-file apparatus with it; the run cohort
  is typed by the subpanel classification and is wired in as a second domain with its own config.
  Neither is a change to the CODEX domain, which still shows what it showed.

Startup log confirmed: DuckDB views `islets` + `tissue` registered; 4 mirai workers; 9,677 islet
spatial records. The channel sidecar is now the **59-channel Senior panel** (`apps/islet/Channel_names`,
regenerated from the qptiff biomarker order; the old 35-plex Junior file is backed up as
`Channel_names.junior35.bak`); `www/local_images` → `/data/islet_ome_tiff`.

## Per-tab matrix

| Tab | Verdict | Evidence (screenshot) | Notes |
|-----|---------|-----------------------|-------|
| **Plot** | ✅ Works | `plot.png` | Sidebar filters, "Insulin⁺ fraction vs Islet Size" scatter (clean ND↘ / Aab+↗ / T1D-low biology), donor-group violin, AI panel — all render. Drill-down controls present (9,706 per-islet CSVs + 15 GeoJSONs on disk). |
| **Trajectory** | ✅ Works (was slow) | `trajectory.png` (8s, mid-load) → after 30s the pseudotime scatter (INS vs pseudotime, ND/Aab+/T1D + trend lines) + both UMAPs render fully. Slow first render was vectorized away 2026-06-29. |
| **Viewer** | ✅ Works (images added Jul 2026) | `verification_6539.png` | Per-donor 59-plex OME-TIFFs from `convert_viewer_ometiff.sh`, served via `www/local_images` → `/data/islet_ome_tiff`. Donor 6539 verified rendering in Avivator with correct INS/GCG/SST/DAPI defaults + `6539.offsets.json` (Viv IFD index). Requires Range-serving (prod nginx/shiny-server), not `runApp`. |
| **Statistics** | ✅ Works | `statistics.png` | Full narrative + Configure Analysis controls + Run Statistics / Download CSV. Result panels (Hypothesis Testing, Effect-Size Forest, AUC) populate **on demand** after Run Statistics; the donor-aware (N=15) mixed-effects workflow is wired. |
| **Spatial** | ✅ Works | `spatial.png` | Tissue scatter renders **145,651 foreground cells** for donor 115 via the **DuckDB tissue backend** (no per-donor CSV fallback needed); Leiden UMAP (9,130 islets, clusters 0–14), Donor-Status UMAP, and "Mean Phenotype Composition by Cluster" all render. |
| **Drill-Down** *(embedded in Plot/Trajectory)* | ✅ Works via Plot | — | Per-islet segmentation + single-cell overlay is driven by a click on a scatter point. Functional from **Plot** (cell CSVs + GeoJSONs present); unusable from **Trajectory** only because that scatter is empty. |
| **AI Assistant** *(Plot side-panel)* | ✅ Renders | `plot.png` | "Panc Floyd" panel + model selector (Navigator `gpt-oss-20b`) shown and ready. Live chat needs UF Navigator API creds in `.Renviron` (not validated here). |

## The Trajectory "failure" — what it actually was (and the fix)

Not a bug in the data or the join — a **performance** problem plus a premature screenshot. Diagnosis
(2026-06-29):

- The trajectory data is fully valid: all **9,130** islets pass the filters for INS (`dpt_pseudotime`
  finite 0–1, INS expression finite, `islet_key`/`case_id` present). The server log even shows
  `[traj_data_clean] Successfully processed 9130 observations for INS`.
- The islet keys match perfectly on both sides (`Islet_755`), so there is **no join defect** —
  `prep_data()` already calls `add_islet_key()` on `comp`, and `traj()` synthesizes
  `combined_islet_id`. The "0 islet metrics" line in the first screenshot was the **initial frame**
  before `input$traj_feature` propagated.
- The real cost was `traj_data_clean` calling `get_islet_annotations()` **once per islet** (9,130×),
  each scanning the full 9.6k-row `islet_spatial_lookup` — ~88M comparisons **on every render**
  (feature change, filter change). First render took ~10–30s, so the 8s screenshot caught an empty
  tab.

**Fix (`mod_trajectory_server.R`):** replaced that per-row loop with a single vectorized `match()`
against the lookup (identical result, O(n)). The full render (pseudotime scatter + both UMAPs) was
verified working via screenshot; the `match()` returns the same `centroid_x/y_um`/`has_spatial`
columns as the loop (same first-match semantics, same lookup columns), so it only changes speed.

## Other observations

- ~~**No Excel fallback.**~~ **Wrong as written, corrected 2026-09-18.**
  `data/app_data/master_results.xlsx` is **present** — 11,254,291 bytes, dated 2025-09-17. So
  `load_master_auto()` does have its Excel fallback on this domain, and the "single point of
  failure" that followed from the claim does not hold. Whether that path was ever exercised is a
  separate question this document does not answer; the file's existence is all that was measured.
  The run cohort resolves `master` through its own config key and does not use this file.
- **DuckDB backend is doing real work** — the Spatial tissue scatter (145k+ cells) comes from
  `data/parquet/tissue` via partition pruning, confirming `ISLET_USE_DUCKDB=TRUE` is engaged.
- **OLD-vs-new phenotypes.** The Spatial phenotype list still shows the full provisional taxonomy
  (Acinar, Ductal, CD4/CD8 T cell, M2 macrophage, NK cell, Regulatory T cell, Stromal, …) — true of
  **this domain**, the CODEX one, whose upstream typing produced those labels. ~~The notebook's
  broad-lineage output (6 classes) is a separate, QC'd artifact and is not surfaced in the app by
  design.~~ **Corrected 2026-09-18:** that sentence no longer describes anything. The
  notebook's broad-lineage lane went with RESTORE (retired 2026-08-14) and the broad-typing
  producer comparison it fed was withdrawn 2026-09-08; neither is a fallback and neither is to be
  revived. The run cohort is typed by the subpanel classification instead — see "What the run
  cohort's cells were typed by" above — and `Acinar` / `Ductal` are CODEX-domain labels that the
  60-plex panel cannot support and that no run-cohort card may show.

## Reproduce

```bash
cd apps/islet
RETICULATE_PYTHON=/home/smith6jt/miniconda3/envs/scvi-env/bin/python \
  Rscript -e 'shiny::runApp(".", port = 7777, host = "127.0.0.1")'
# then load http://127.0.0.1:7777 and click through Plot / Trajectory / Viewer / Statistics / Spatial
```
