# Scripts guide

## Ownership

`core/phenocycler/` owns the raw QuPath-to-expression workflow; its development lanes live under
`core/phenocycler/scripts/`. Root `scripts/` builds application data and viewer assets.

Root scripts must not invent marker thresholds, positivity frequencies, or alternate typing rules.
Where a level is needed, it comes from a reference population that a lane already fitted.

## The boundary that still binds

The app's aggregation, trajectory and spatial scripts were written against a legacy single-cell
H5AD. They are not yet a consumer of the content-addressed run, and must not be presented as an
end-to-end rebuild of it.

The adapter that closes the gap reads the run's own stages, joins on **`(donor_id, object_id)`**,
preserves uncertainty, and records the run id and the manifest content ids. Existing H5AD phenotype
columns are never calibration or typing inputs.

Join by exact `object_id`. Never row order, rounded centroids, filenames, or a provisional
`classification`. A missing or invalid measurement stays unavailable.

## Environments

- Repo `.venv` — the pipeline, the integration lane, and most development lanes. **Globus Connect
  Personal's daemon runs under this interpreter; never recreate it.**
- `scvi-env` — the Xenium typer and the downstream Scanpy/scVI/Harmony work.
- `~/.venvs/valis`, `external/.venv-som`, `external/.venv-nimbus` — single-purpose, incompatible pins.
- R 4.5.1 with `~/R/…/4.5` for every R lane.

Do not add embedding, trajectory or app-only dependencies to the focused source environment.

## What is retired

RESTORE is retired for the source workflow (2026-08-14), and the flat-file apparatus with it. The
following are historical and must never be live inputs or fallbacks: the Gate 2–5 divisor/floor
screens, the GMM/Otsu gate searches, ordered-residual broad-lineage assignment, `data/cells_redsea*`
and `data/restore*` as source roots, and the broad-typing producer comparison withdrawn 2026-09-08.
Git history has all of it.

`core/phenocycler/scripts/dct/` is a closed lane but exports two symbols other lanes import
(`gates.require_evidence`, `tiff_window.read_window`) — move those before deleting anything there.

The QuPath handoff is `core/phenocycler/scripts/groovy/import_cell_assignments.groovy`: UUID-keyed,
preserves uncertainty, verifies each CSV's size and SHA-256 against its manifest before changing a
classification, and records each detection's previous class.

## Cohort and batch tables

- `scripts/build_acquisition_batch_table.py` → `data/donor_acquisition_batch.csv` — panel variant and
  Fusion version, read from the qptiff headers.
- `scripts/build_acquisition_runs.py` → `data/donor_acquisition_runs.csv` — the **real** batch
  structure: acquisition run (scan day) and staining group (consecutive days), from the TIFF
  `DateTime` tag. Use these for any batch term; the Fusion version alone has no ND donor at 2.4.1.
- `scripts/build_donor_metadata.py` → the active domain's `donor_metadata.csv` — cohort sheet joined
  to the gitignored nPOD export. Fails loudly on an unresolved donor rather than dropping it.

## Viewer and application assets

The Viewer runs on either backend via `ISLET_VIEWER_BACKEND=avivator|tissuumaps`.

- `senior/build_tissuumaps_project.py` — per-donor `.tmap`, naming channels from each slide's own
  OME-XML. `--regions` takes a QuPath export directory and refuses to guess an ambiguous match;
  `--markers` refuses cells that would land outside the slide rather than attaching them.
- `senior/build_all_tissuumaps_projects.py` — the cohort, with exclusions imported from
  `phenocycler.cohort` rather than restated.
- `senior/export_viewer_cell_markers.py` — the whole-slide per-cell CSV, into
  `data/tissuumaps_inputs/markers/`. **Not** into the drill-down's `cells/` namespace, which is
  globbed for any `*.csv`.
- `senior/convert_viewer_ometiff.sh` — qptiff → pyramidal OME-TIFF on `/data`.

**The cohort mixes panels**, so default-on channels resolve **by marker name, never by index**.

HTTP Range (206) is an Avivator/Viv requirement and production must serve it; TissUUmaps tiles
images itself and needs none, which is why that backend also works under `runApp`.

Generated data, figures, run directories, caches and executed notebooks stay out of source control.
