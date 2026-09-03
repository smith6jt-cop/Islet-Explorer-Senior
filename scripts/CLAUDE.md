# Scripts guide

## Ownership boundary

`core/phenocycler/` owns the raw QuPath-to-expression workflow. Root
scripts will consume a declared identity release after a run-aware H5AD adapter
is implemented (`docs/PLAN_2026-09-02_AUXILIARY_PACKAGES.md`, Track D.3); current
app scripts still consume a legacy H5AD. Root scripts must not invent marker
thresholds, positivity frequencies, or alternate cell-typing rules.

Current source stages (`docs/WORKFLOW_2026-09-02_MAP.md` § 1):

```text
manifest -> ingest -> geometry -> redsea -> expression
```

`controls`, `calibrate`, `type` and the QuPath assignment export are blocked
(`NO_VERIFIED_TYPING`); `states` is explicit-only. Cell identities exist only in
the mutable development lanes under `core/phenocycler/scripts/` (map § 2).

Run and inspect them through the package:

```bash
python -m pip install -e core/phenocycler        # into the repo .venv
cd core/phenocycler
python -m phenocycler.pipeline run --config config.ini --pipelined --through expression
python -m phenocycler.pipeline status --config config.ini
```

`senior/run_pipeline.sh` is a convenience wrapper for `run`, `status`, and
`export`; it adds no workflow logic. Set `PHENOCYCLER_CONFIG` only when a
different configuration file is intentional.

See:

- `../core/phenocycler/README.md`
- `../core/phenocycler/docs/WORKFLOW.md`
- `../core/phenocycler/DATA_README.md`
- `../docs/WORKFLOW_2026-09-02_MAP.md`

## Downstream scripts and current boundary

The existing aggregation, trajectory, spatial, and application scripts operate
on a legacy single-cell H5AD. They remain useful for the deployed application,
but they are not yet a supported consumer of the content-addressed source
workflow and must not be presented as a current end-to-end rebuild.

Before those scripts are used with new assignments, implement one canonical
run-to-H5AD adapter that reads:

```text
<run>/cells
<run>/expression
<run>/marker_evidence
<run>/assignments
<run>/cell_states
```

The adapter must join on `(donor_id, object_id)`, preserve uncertainty, and
record the run ID and relevant manifest content IDs. Existing H5AD phenotype
columns are not valid calibration or typing inputs. `fix_tcell_subtypes.py` is
legacy H5AD maintenance only and must not modify current assignments.

## Environment separation

- `scvi-env` is the established downstream Scanpy/scVI/Harmony environment.
- The source workflow runs from the repo `.venv` today; `core/phenocycler/environment.yml`
  is the HiPerGator recipe (the `phenocycler-analysis` env does not exist on this workstation).
- Do not add embedding, trajectory, Shiny, or app-only dependencies to the
  focused source environment.

## Donor and identity contracts

- Central exclusions for `6457` and `6579` are enforced by
  `phenocycler.cohort`; scripts must not recreate or override the list.
- Join source artifacts by exact `object_id`.
- Do not use row order, rounded centroids, filenames, or provisional
  `classification` values as identity.
- Missing/invalid marker evidence stays unavailable or indeterminate.
- Disease labels and annotations are audit-only.

## Retired phenotyping code

RESTORE is retired for the source workflow (2026-08-14; root `CLAUDE.md`):
there are no verified target/reference pairs, `restore_stage.py` is reachable
only through the blocked `calibrate` stage, and no compatibility fallback may
revive it. The flat-file apparatus below is retired with it. Note the repo-root
`data/restore_gated_redsea/` and `data/phenotype/broad/` are still read by the
Shiny app's Card C — retired for the source workflow, not deleted for the app.

The following approaches are historical and must not be used as live inputs or
fallbacks:

- the standalone Gate 2-5 divisor/floor screens and their dossiers;
- the GMM/Otsu gate searches run outside the pair evaluation;
- ordered residual broad-lineage assignment;
- `data/cells_redsea*`, `data/restore*`, and `data/phenotype/broad*` as current
  source-stage roots;
- Gate 5 dashboards/counterfactuals and anatomy-check release arms;
- hard-coded QuPath import/export paths and the old broad-lineage importer.

The QuPath handoff is
`core/phenocycler/scripts/groovy/import_cell_assignments.groovy`, which
uses UUIDs and preserves uncertainty; it has no current input while assignment
export is blocked. Historical experiments remain
recoverable from Git history.

After importing, `groovy/view_phenotypes.groovy` can show or hide any exact
broad, specific, ambiguous, or unavailable QuPath class.

## Viewer and application assets

OME-TIFF conversion and browser-specific viewer behavior are downstream
application concerns. Keep channel order explicit and follow
`../apps/islet/CLAUDE.md` for viewer, plotting, reticulate, and deployment
details.

The Viewer runs on either backend (`ISLET_VIEWER_BACKEND=avivator|tissuumaps`):

- `senior/build_tissuumaps_project.py` — builds the per-donor `.tmap` project,
  naming channels from each slide's own OME-XML (the cohort mixes panels).
  `--regions` takes the QuPath export directory (`data/Islets`) and matches by
  case id, refusing to guess when the match is ambiguous; the coordinate scale
  into full-res pixels is derived from the two pixel sizes and then checked
  against the image bounds rather than assumed. `--structures` adds vessels,
  nerves and lymphatics as a **second region layer that does not auto-load** —
  those classes run to ~200k features per donor. `--markers` is falsified the
  same way regions are: cells that would land outside the slide are refused, not
  attached, and the run exits 3.
- `senior/build_all_tissuumaps_projects.py` — builds the whole cohort. Donors are
  parsed from the `MAP` in `convert_viewer_ometiff.sh` and the exclusions are
  imported from `phenocycler.cohort`; neither is restated, and the script exits 1
  rather than proceeding if it cannot import the exclusions. Ends on a per-donor
  coverage table naming any donor missing outlines or cells.
- `senior/export_viewer_cell_markers.py` — writes the whole-slide per-cell CSV the
  builder's `--markers` wants, from the Parquet tissue partitions (legacy
  `app_data/donors/*.csv` as fallback), into `data/tissuumaps_inputs/markers/`.
  **Not** into `data/app_data/cells/`: that is the per-islet drill-down namespace
  and `drilldown_available()` globs it for any `*.csv`.
- `senior/build_viewer_assets.py --structures` — writes
  `app_data/json/structures/<case>_structures.geojson` for the viewer only.
  `CLASS_MAP` stays islets-only; merging structures into the app's own GeoJSON is
  what the v1 deferral rules out (`sf::st_read` on 200k features).
- `test_build_tissuumaps_regions.py` — unit tests for that region handling
  (directory lookup, gzip, scale inference). No slide or server needed.
- `senior/tissuumaps_smoke_test.py` — boots a server against a built project and
  asserts every route the app depends on; `--browser` also checks rendering. It
  logs the version and install path of the TissUUmaps it exercises — that should
  be our fork (`smith6jt-cop/TissUUmaps`), which is what the deployment runs, not
  the PyPI release.
- [`../docs/tissuumaps_evaluation.md`](../docs/tissuumaps_evaluation.md) —
  findings, costs, and the deployment recipe.

HTTP Range (206) is an **Avivator/Viv** requirement: it streams pyramid tiles by
range, so production must serve it. TissUUmaps tiles images itself and needs no
Range support, which is why the Viewer also works under `shiny::runApp` on that
backend.

Generated data, figures, run directories, caches, and executed notebooks stay
out of source control.
