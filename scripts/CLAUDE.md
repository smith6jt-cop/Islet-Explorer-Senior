# Scripts guide

## Ownership boundary

`Phenocycler_Analysis/` owns the raw QuPath-to-cell-identity workflow. Root
scripts will consume completed assignments after a run-aware H5AD adapter is
implemented; current app scripts still consume a legacy H5AD. Root scripts
must not invent marker thresholds, positivity frequencies, or alternate
cell-typing rules.

Current source stages:

```text
manifest -> ingest -> geometry -> redsea -> expression -> controls
         -> calibrate -> type -> states -> QuPath export
```

Run and inspect them through the package:

```bash
conda activate scvi-env
python -m pip install -e Phenocycler_Analysis
python -m phenocycler.pipeline run --config config.ini
python -m phenocycler.pipeline status --config config.ini
```

`senior/run_pipeline.sh` is a convenience wrapper for `run`, `status`, and
`export`; it adds no workflow logic. Set `PHENOCYCLER_CONFIG` only when a
different configuration file is intentional.

See:

- `../Phenocycler_Analysis/README.md`
- `../Phenocycler_Analysis/docs/WORKFLOW.md`
- `../Phenocycler_Analysis/DATA_README.md`

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
- `Phenocycler_Analysis/environment.yml` is the focused reproducible source
  workflow environment.
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

The following approaches are historical and must not be used as live inputs or
fallbacks:

- RESTORE normalization and divisor/floor searches;
- NNMF/GMM/Otsu Gate 2–5 screens;
- ordered residual broad-lineage assignment;
- `data/cells_redsea*`, `data/restore*`, and `data/phenotype/broad*` as current
  source-stage roots;
- Gate 5 dashboards/counterfactuals and anatomy-check release arms;
- hard-coded QuPath import/export paths and the old broad-lineage importer.

The active QuPath handoff is
`Phenocycler_Analysis/scripts/groovy/import_cell_assignments.groovy`, which
uses UUIDs and preserves uncertainty. Historical experiments remain
recoverable from Git history.

After importing, `groovy/view_phenotypes.groovy` can show or hide any exact
broad, specific, ambiguous, or unavailable QuPath class.

## Viewer and application assets

OME-TIFF conversion and browser-specific viewer behavior are downstream
application concerns. Keep channel order explicit and follow
`../app/shiny_app/CLAUDE.md` for viewer, plotting, reticulate, and deployment
details.

The Viewer runs on either backend (`ISLET_VIEWER_BACKEND=avivator|tissuumaps`):

- `senior/build_tissuumaps_project.py` — builds the per-donor `.tmap` project,
  naming channels from each slide's own OME-XML (the cohort mixes panels).
- `senior/tissuumaps_smoke_test.py` — boots a server against a built project and
  asserts every route the app depends on; `--browser` also checks rendering.
- [`../docs/tissuumaps_evaluation.md`](../docs/tissuumaps_evaluation.md) —
  findings, costs, and the deployment recipe.

HTTP Range (206) is an **Avivator/Viv** requirement: it streams pyramid tiles by
range, so production must serve it. TissUUmaps tiles images itself and needs no
Range support, which is why the Viewer also works under `shiny::runApp` on that
backend.

Generated data, figures, run directories, caches, and executed notebooks stay
out of source control.
