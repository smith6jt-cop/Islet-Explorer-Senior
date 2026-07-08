# Plan (ACTIVE): Remove DAPI, IAPP, Ker8_18, Keratin_5 from all figures & analysis

## Context

The Step-1 identity heatmap (`data/phenotype/celltype_marker_heatmap.png`) and the app's marker
selectors surface four markers to drop: **DAPI** (nuclear stain), **IAPP** (failed marker), and the
two redundant epithelial keratins **Ker8_18** and **Keratin_5**. **`Pan_Cytokeratin` STAYS** — it is
the RESTORE universal negative-control reference and the sole Epithelial gate; removing it would break
phenotyping. Goal: these four appear in no figure, plot, selector, heatmap, or analysis output.

**Functional safety (verified):** none of the four is used as a gate/threshold/phenotype rule.
`Ker8_18`/`Keratin_5` are TYPE-display markers only; `DAPI` appears only in a guarded REDSEA
alignment-QC probe + an image-viewer color; `IAPP` is unused. `data/phenotype_rules.csv` contains none
of them. So moving all four to the central `EXCLUDED` set is safe.

## Change set

### A. Pipeline identity figures — one central edit + regenerate
`scripts/senior/marker_taxonomy.py` is the single source of the TYPE/PROCESS/EXCLUDED split:
- **Line 36:** `EXCLUDED = ["DAPI", "IAPP", "Ker8_18", "Keratin_5"]` (add the two keratins; DAPI/IAPP
  already present).
- **Line 44:** remove `"Ker8_18", "Keratin_5"` from the `TYPE` epithelial block.
- Update the docstring (lines 6, 14): "epithelial = Pan_CK/Ker8_18" → "Pan_CK"; note the two keratins
  are now excluded.

This propagates automatically to `plot_celltype_markers.py` (uses `heatmap_markers()`) and the notebook
Step-F figure. Then **regenerate the stale figures**: the on-disk `celltype_marker_heatmap.png` is a
stale cached PNG (it still shows DAPI/IAPP *and* PROCESS markers the current code already excludes) —
delete `data/phenotype/celltype_marker_{heatmap,dotplot}.png` and run
`python scripts/senior/plot_celltype_markers.py`. **Note:** the refreshed heatmap will also drop the
PROCESS/state markers (Ki67, IDO1, CD107a, …) — that is the *existing* intended TYPE-only identity
design that the stale PNG predates, not a new change. The dotplot (10 RESTORE + B3TUBB/CD99/MPO) has
none of the four and is unaffected. No data rebuild needed for the figure.

### B. Shiny app (running deployment) — hide them from every selector now, no rebuild
The four surface at **4 UI enumeration points**, derived on-the-fly from `var_names`/CSV columns (no
central marker list today). Add one constant and filter the CHOICE vectors (not the underlying data):
- `app/shiny_app/R/00_globals.R`: add `EXCLUDED_MARKERS <- c("DAPI","IAPP","Ker8_18","Keratin_5")` and
  mirror it into `globalenv()` in the existing `local({...})` block (~:269-274) so auto-sourced modules
  see it.
- `R/mod_plot_server.R:1168` — append `EXCLUDED_MARKERS` to the drill-down `setdiff` blocklist.
- `R/mod_trajectory_server.R:870` — same append (drill-down duplicate).
- `R/mod_trajectory_server.R:176` — `available_features <- setdiff(tr$var_names, EXCLUDED_MARKERS)`
  (Select-Feature dropdown).
- `R/mod_trajectory_server.R:1347-1348` — filter the multi-heatmap `choices`/`selected` with
  `setdiff(tr$var_names, EXCLUDED_MARKERS)`.

**CRITICAL:** filter only the UI `choices`/`selected` — do **NOT** filter `tr$var_names` at its source;
`mod_trajectory_server.R:1359` uses it as a positional index into `adata$X`, so shortening it selects
the wrong expression columns. (Optional DRY: extract the duplicated drill-down blocklist into a
`drilldown_marker_choices()` helper in `R/drilldown_helpers.R`.)

### C. Rebuild-time data filter — drop them from app-facing data at the 22-donor rebuild
The app's marker table and drill-down CSVs PHYSICALLY carry the four today
(`data/parquet/uns_markers.parquet`, `data/app_data/cells/*.csv`), produced by schema-derived (not
taxonomy-driven) scripts. Apply the `marker_taxonomy.EXCLUDED` filter at the producers so the already-
planned rebuild drops them:
- `scripts/senior/compute_islet_markers.py:22` — filter `markers` against `marker_taxonomy.EXCLUDED`
  before `marker_pos_frac` (propagates to the `markers` sheet → `uns_markers.parquet`).
- `scripts/senior/extract_islet_cells.py:49,70` — filter `markers` for the drill-down CSVs.
- (`scripts/senior/build_app_h5ad.py:105` — optional backstop filter on `mk_names`.)

Base pipeline files (`data/cells/`, `data/cells_redsea/`, `singlecell_protein*.h5ad`) KEEP the physical
columns — REDSEA/RESTORE need all physical channels and downstream selects markers by name, so the
columns are inert there.

### Scoping decisions (stated for approval)
- **Raw image Viewer (Avivator) = OUT OF SCOPE.** It displays the physical OME-TIFF channels; the four
  are real stained channels and can't be removed without rebuilding the OME-TIFFs. (Optional one-liner:
  flip DAPI to default-off at `R/viewer_helpers.R:369-374` — say if you want it.)
- **Docs/comments** mentioning the four (`assign_broad_lineage.py`, `restore_normalize.py`,
  `scripts/CLAUDE.md`, `docs/*`) are cosmetic — optional cleanup, no behavior change.

## Critical files
- `scripts/senior/marker_taxonomy.py` (lines 36, 44, docstring) — central figure/taxonomy lever.
- `scripts/senior/plot_celltype_markers.py` — regenerate (no edit).
- `app/shiny_app/R/00_globals.R`, `R/mod_plot_server.R`, `R/mod_trajectory_server.R` — app selectors.
- `scripts/senior/compute_islet_markers.py`, `extract_islet_cells.py` (+ optional `build_app_h5ad.py`)
  — rebuild-time data filter.

## Verification
- Re-render `celltype_marker_heatmap.png`: no DAPI/IAPP/Ker8_18/Keratin_5 columns; `Pan_Cytokeratin`
  retained; panel is TYPE-only.
- Launch dev app (`cd app/shiny_app && Rscript -e 'shiny::runApp(".", port=7777)'`): the Trajectory
  "Select Feature", the multi-feature heatmap, and both drill-down "Color by" dropdowns no longer list
  the four; `Pan_Cytokeratin` still present; the trajectory heatmap still colors correctly (positional
  index intact).
- After a filtered producer run / rebuild: `uns_markers.parquet` `marker` column and a drill-down CSV
  header no longer contain the four.

---

> **STATUS: DEFERRED — keep on file, do NOT execute (per user, 2026-07-07).** The REDSEA reassessment
> below is preserved for when the pixel data is available; the ACTIVE task above is the marker removal.

# Plan (DEFERRED): Reassess & improve REDSEA correction — reduce keratin & Vimentin spillover into endocrine / other cells

## Context

The Step-1 QC figure `data/phenotype/celltype_marker_heatmap.png` shows the **Endocrine** lineage
lit up for the epithelial block — **Ker8_18, Pan_Cytokeratin, Keratin_5** (and EpCAM/E-cadherin) —
i.e. residual lateral **keratin spillover** sitting on cells that were *correctly* typed Endocrine.
The goal is to reduce keratin spillover into small endocrine clusters / single endocrine cells
outside large islets, and to reduce Vimentin spillover into non-mesenchymal cells — **without**
eroding real epithelial/ductal keratin or real vessel/fibroblast Vimentin, and **without** assuming
cohort-uniform frequencies (per the repo working rule: validate per-donor on that donor's own data).

### What the diagnostic already established (read-only, all 22 donors)

Mean **corrected** (α=1) marker intensity for cells *called* Endocrine, by islet context:

| Endocrine stratum | n cells | Pan_CK | Ker8_18 | Vimentin |
|---|---|---|---|---|
| core, large islet (>60 cells) | 680k | 180 | 315 | 357 |
| core, small islet (≤15) | 1.6k | 300 | 486 | 398 |
| peri-islet | 218k | 576 | 1029 | 540 |
| **lone / tissue** | **1.0M** | **685** | **1130** | **659** |

- Contamination is **spatial**: lone/scattered endocrine carry ~3.8× the keratin of large-islet-core
  endocrine, and are the **largest** stratum (~63% of all endocrine calls). Pooled endocrine keratin
  is ~0.85× the epithelial mean (per-donor 0.5–1.4). **Ker8_18 and Pan_Cytokeratin dominate;
  Keratin_5 is minor** (~30–50 units).
- These peri/scattered endocrine cells are **weak-hormone** (median hormone `_norm` ≈ 1.6× threshold
  vs ≈ 16× in islet cores) but morphologically normal single cells (solidity 0.91, not doublets) —
  a boundary population sitting right at the hormone threshold *and* soaked in acinar keratin/Vimentin.

## Why REDSEA under-serves these cells today

`scripts/senior/redsea_full.py::compensate` (line ~221):
`corrected = clip(data − α·(F@edge), 0) / size`. It is a **single global α=1**, **channel-diagonal**,
purely **geometric** (F = row-normalized contact fractions), **subtract-only**, **1-px band**. No
per-marker selectivity, no brightness/direction awareness. Small/peri/scattered endocrine cells have
high perimeter-to-area and sit in Pan_CK+ acinar/duct tissue, so they receive disproportionate keratin
spillover that α=1 + a 1-px band only partly removes — exactly the heatmap residual.

## Constraints & non-obvious facts that shape the fix

- **Only `Pan_Cytokeratin` gates cell *typing***; `Ker8_18`/`Keratin_5` are inert for the lineage
  call but appear in the heatmap → include all three for **marker-value** cleanliness.
- **RESTORE largely cancels a global per-channel α** for `_pos`/`_norm` (the threshold scales with the
  channel), so re-running RESTORE after any REDSEA change is **mandatory**; the surviving effect is the
  per-cell *differential*. The heatmap itself plots **corrected MFI** (not `_norm`), so the keratin α
  acts on it directly.
- **The parquet α-simulation** `max(0, corr − (α−1)·clip(raw−corr,0))` is **exact only at α=1**; for
  α>1 it is directionally right but the **guardrail cannot be certified** from parquets (QuPath's mean
  ≠ REDSEA's internal raw mean; r≈0.998 not 1.0; NaN edge-cells). Empirically it blew up the epithelial
  guardrail ratio → **definitive tuning requires the pixel re-run** (data being transferred).
- **Vimentin is the sole Fibroblast marker** (no redundancy) and genuinely co-expresses on vessels
  (CD31+SMA+Vimentin+) → a **global Vimentin α is unsafe**; use a **context-gated / directional** rule.
- Real CD31+ vessels are protected by the lineage hierarchy (CD31 override), but SMA+/Vim+ mural cells
  and Vimentin-homogeneous fibroblast bands are **not** — over-subtraction silently defaults them to
  Epithelial. Guard against this explicitly.

## Approach (staged)

### Phase 0 — Formalize the diagnostic as the acceptance yardstick (this machine, read-only)
Turn the prototypes in `scratchpad/diag_*.py` into a committed
`scripts/senior/redsea_reassess_diag.py` that, per donor and **stratified by islet context**
(core large/med/small · peri · lone/tissue, from `islet_num`/`cell_region`/`dist_islet`), reports:
- **Objective:** endocrine-typed mean corrected keratin (3) + Vimentin vs the Epithelial reference;
  impossible co-expression = fraction of confident-endocrine that are `Pan_Cytokeratin_pos`.
- **Guardrails (must NOT regress):** confident-acinar keratin retention (high-corrected-Pan_CK,
  `tissue`, large `dist_islet` stays Pan_CK+/Epithelial); confident-vessel Vimentin retention (CD31+
  /SMA+ stays Vimentin+/Endothelial); Fibroblast-compartment size; T1D hallmarks (β-loss INS⁺, T-gain
  CD3e⁺ across ND→T1D).
- **Real-vs-false endocrine:** hormone `_norm` strength + morphology (solidity, `cell_maxdiam`) to
  separate genuine scattered endocrine from spillover false-positives and segmentation doublets — so we
  never "fix" a segmentation artifact with α. "Confident endocrine" = donor-percentile hormone (not a
  cohort constant), so the metric stays frequency-free.

### Phase 1 — Make REDSEA correction selective (additive; default behavior unchanged)
Edit `scripts/senior/redsea_full.py::compensate` and CLI:
- **Per-channel α vector**: scalar `alpha` → length-59 vector broadcast over the columns of `(F@edge)`
  (one-line change at line ~221). New `--alpha-per-channel "Pan_Cytokeratin=2,Ker8_18=2,Keratin_5=1.5"`;
  **default keeps α=1 for every channel** so current output is bit-for-bit unchanged until opted in.
- **Directional (brightness-gated) subtraction** for the aggressor channels — the key lever: scale the
  subtraction for cell *i*, channel *c*, by a **receptivity weight** derived from *i*'s **interior
  density** `= (data−edge)/(size−bandcount)`. Cells bright-in-interior (real acinar/ductal keratin;
  real fibroblast/vessel Vimentin) are **protected**; dim victims (endocrine soaked in neighbor keratin;
  immune/endocrine soaked in Vimentin) are **cleaned**. This protects real biology *by construction* and
  is what makes the **Vimentin** reduction safe. Needs one extra `np.bincount(mask_band,…)` (per-cell
  band-pixel count) in `process_donor` + persisted in `save_intermediates`.
- Optionally gate Vimentin on **orthogonal context** (protect a Vimentin+ cell that is CD31+/SMA+ or in
  a stromal band) — non-circular because the gate uses *other* channels/space, not the corrected channel.

### Phase 2 — Execute on pixels, re-RESTORE, re-assign, re-plot (needs transferred qptiffs + GeoJSONs)
1. Regenerate intermediates + masks for **2–3 representative donors** spanning ND/Aab+/T1D
   (`redsea_full.py --donor <id> --save-intermediates --keep-mask`).
2. Sweep α / gating cheaply via the **extended validator** (below), pick the operating point on the
   Phase-0 objective.
3. Full re-run: `redsea_full.py` (new params) → `restore_normalize.py --cells-dir data/cells_redsea …
   --no-robust` (causal check; `--robust` only as production preview — it imputes toward the cohort
   median and would violate the per-donor rule) → `assign_broad_lineage.py` → `plot_celltype_markers.py`.
4. **Accept only if** endocrine keratin/Vimentin drop materially **AND** all Phase-0 guardrails hold.

### Phase 3 (OPTIONAL — only if you want the false-endocrine / small-endocrine-recovery angle)
Restrict the `INS/GCG/SST` RESTORE negative-control reference to Pan_CK+ acinar with **large
`dist_islet`** (uncontaminated background → correct threshold), so weak scattered endocrine are
correctly *kept vs dropped*. This is orthogonal to keratin and targets the weak-hormone boundary
population the diagnostic surfaced. `reference_exclusivity_qc()` already computes `refpos_dist_islet_med`.

### Extend the validator (`scripts/senior/redsea_full_validate.py`)
Add SST to the Endocrine compartment; add keratin-in-endocrine and Vimentin-in-non-vascular metrics;
**stratify by islet context**; and — critically — **re-RESTORE per α instead of holding gates fixed**
(the current gate-fixed harness *overstates* the effect). Bootstrap SSC thresholds over a few seeds so
the lever effect is distinguished from threshold noise.

## Critical files
- `scripts/senior/redsea_full.py` — `compensate` (line ~221: scalar α → per-channel/gated vector),
  argparse, `process_donor` (add band-pixel count), `save/load_intermediates`.
- `scripts/senior/redsea_full_validate.py` — extend objectives + **re-RESTORE**.
- `scripts/senior/restore_normalize.py` — re-run on the new corrected MFI; (Phase 3) the
  `dist_islet`-restricted hormone reference in `DEFAULT_MARKER_PAIRS` / `idx_select`.
- `scripts/senior/assign_broad_lineage.py` — **hard-coded** `GATED`/`EXTRA_GATED` paths (lines ~66–71)
  need CLI/env indirection so the reassessment loop can point at a sim/gated dir.
- `scripts/senior/plot_celltype_markers.py` — regenerate the heatmap to confirm the endocrine row.
- `scripts/senior/redsea_reassess_diag.py` — **new**, the Phase-0 diagnostic (from `scratchpad/diag_*.py`).
- Inputs already present: `data/cells/`, `data/cells_redsea/`, `data/restore_gated_redsea/` (22 donors).

## Verification
- Re-render `celltype_marker_heatmap.png`: endocrine keratin (Ker8_18/Pan_CK/Keratin_5) z-scores drop
  toward baseline while the **Epithelial** keratin row is essentially unchanged (guardrail).
- Phase-0 per-donor objective/guardrail table: improves on target strata, holds on all guardrails;
  T1D hallmarks (β-loss, T-gain) unchanged.
- QuPath UUID overlay spot-check on a handful of small/scattered endocrine cells
  (`scripts/groovy/import_broad_lineage.groovy`) — they read endocrine, not epithelial.

## Dependencies / what I need from the data transfer
- The **59-plex `.qptiff` slides** (were at `/home/smith6jt/IO60panc2nd/Images/`) and the **per-cell
  GeoJSONs** (`data/redsea_scratch/geojson/cells__<image>.geojson`, feature `id` == `object_id`) — the
  two inputs `redsea_full.py::resolve_paths` needs. I'll regenerate the pixel intermediates + masks
  (`redsea_scratch/` was deleted). `data/cells/` etc. are already here.
- Note: this is a **plan only** — no pipeline runs until you approve and the pixel data lands.
