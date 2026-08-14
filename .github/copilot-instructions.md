## Islet-Explorer — coding agent quickstart

Purpose: Interactive R Shiny app to explore human pancreatic islet measurements by donor group and islet size, plus a lightweight embedded image viewer (Avivator/Viv) for local OME-TIFFs.

Big picture
- App entrypoint: `apps/islet/app.R` (single-file Shiny app). Data come from `data/master_results.xlsx` with sheets: `Islet_Markers`, `Islet_Targets`, `Islet_Composition`, and optional `LGALS3`.
- Core server data flow:
  1) `load_master()` reads sheets via readxl.
  2) `prep_data()` normalizes and joins into three tables: `targets_all`, `markers_all`, `comp`, and computes islet diameter (µm) from Islet_Targets CORE area: `islet_diam_um = 2*sqrt(core_region_um2/pi)`.
  3) App UI tabs use these tables for plots, stats, QA, and the embedded viewer.
- Donor groups and ordering are fixed: ND, Aab+, T1D. Colors: ND `#1f77b4`, Aab+ `#ff7f0e`, T1D `#d62728`.

Project-specific conventions you must preserve
- Keying: All joins/grouping use (`Case ID`, `Donor Status`, `islet_key`). `islet_key` is derived by `add_islet_key()` from `region` like `Islet_200_core`, or `name`/`islet_id` fallback.
- Region labels are normalized to exact tokens: `islet_core`, `islet_band`, `islet_union`. Do not introduce new spellings.
- AAb flags: columns `AAb_GADA`, `AAb_IA2A`, `AAb_ZnT8A`, `AAb_IAA`, `AAb_mIAA`. AAb filters apply only within the Aab+ group (Any/All logic); ND/T1D are unaffected.
- Synthesis rules (important):
  - Targets: If union rows are missing, synthesize counts only as `union = core + band`. Do NOT fabricate `region_um2`, `area_um2`, or `area_density` for union.
  - Markers: If union rows are missing, synthesize `n_cells` and `pos_count` as sums (compute `pos_frac` from them). If band missing but core+union exist, backfill band = union − core (clamp to non-negative). Reattach AAb flags after synthesis.
- Normalization options for plotting: `none`, `global z-score`, `robust per-donor` (median/MAD; MAD×1.4826). Apply consistently in both main plot and distribution plots.
- Stats tab: Global model is `value ~ donor_status + islet_diam_um`; pairwise t-tests on residuals from `value ~ islet_diam_um` with BH adjustment. Respect factor order ND < Aab+ < T1D.

Embedded viewer (Avivator/Viv)
- Static bundle lives under `apps/islet/www/avivator/`. Install via `scripts/install_avivator.sh` (mirrors public site; or builds if source present and Node+pnpm are available).
- Optional channel mapping from `apps/islet/Channel_names` or `Channel_names.txt` lines like `INS (C26)`. `build_channel_config_b64()` encodes a config sent to the viewer; it highlights INS (red), GCG (blue), SST (yellow), DAPI (grey) if present.
- Local OME-TIFFs can be served by setting `LOCAL_IMAGE_ROOT` or placing files in `local_images/` (resource path `/local_images`).

Developer workflows
- Install R deps: run `scripts/install_shiny_deps.R` (adds CRAN pkgs; Bioconductor extras are optional and not required for current app features).
- Run the app:
  - VS Code task: “Run Shiny app to reproduce issues”.
  - Or shell: `R -q -e "shiny::runApp('apps/islet', launch.browser=FALSE)"`.
- Build `data/master_results.xlsx` from TSVs: `python scripts/build_master_excel.py` (expects `data/results/*.tsv` and `CODEX_Pancreas_Donors.xlsx`; writes `data/master_results.xlsx`).
- Alternate donor compilation: `python scripts/compile_donors.py --input donors.zip --output out.xlsx [--donor-map-xlsx CODEX_Pancreas_Donors.xlsx]`.

Data expectations (examples)
- Islet_Targets: `region` like `Islet_64_core`, `class`, `region_um2`, `area_um2`, `area_density`, `count`.
- Islet_Markers/LGALS3: `region` or `region_type`, `marker`, `n_cells`, `pos_count`, `pos_frac`.
- Composition: `cells_total`, `Ins_any`, `Glu_any`, `Stt_any` (fractions computed as 100×value/cells_total in app).
- Donor metadata: `Case ID` (string; zero-padded), `Donor Status` in {ND, Aab+, T1D}, optional AAb columns.

Quality checks and diagnostics
- QA tab enforces the “3× regions” rule per donor: markers/targets should have 3× distinct region rows vs composition islets. Console emits NA audits via `audit_na()` on load.

When making changes
- Keep the three prepared tables’ schema and keys stable (`islet_key`, region tokens, donor ordering). If you add metrics, thread them through `raw_df_base()` → `plot_df()` → `summary_df()` and update labels/metric selectors accordingly.

Questions to confirm with maintainers
- Provide a sample `CODEX_Pancreas_Donors.xlsx` layout and minimal `data/results/*.tsv` examples if contributing data pipelines.