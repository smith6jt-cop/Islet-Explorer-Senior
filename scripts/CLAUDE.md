# Scripts guide

## Ownership

`core/phenocycler/` owns the raw QuPath-to-expression workflow; its development lanes live under
`core/phenocycler/scripts/`. Root `scripts/` builds application data and viewer assets.

Root scripts must not invent marker thresholds, positivity frequencies, or alternate typing rules.
Where a level is needed, it comes from a reference population that a lane already fitted.

## The boundary that still binds

The app's aggregation, trajectory and spatial scripts were written against a legacy single-cell
H5AD. They are not an end-to-end rebuild of the source workflow and must not be presented as one.

The adapter that closes the gap reads the pipeline's own stages and joins on
**`(donor_id, object_id)`**, preserving uncertainty. Existing H5AD phenotype columns are never
calibration or typing inputs.

Join by exact `object_id`. Never row order, rounded centroids, filenames, or a provisional
`classification`. A missing or invalid measurement stays unavailable.

**There is no run id and no content id to record.** The content-addressed run framework was removed
on 2026-09-16. Stages write to stable roots — `data/cells`, `data/geometry_qc`, `data/redsea`,
`data/expression` — and a donor is done when its partition has content; a stage whose output is
already present is skipped, and `--force` reruns it. `core/phenocycler/data/qupath_manifest.json`
still records a `content_id` and per-file sha256, but nothing validates them: it is a description,
not a check, and it can be hand-edited when a source file moves. Runs made before 2026-09-16 remain
under `core/phenocycler/data/runs/<run_id>/`; that tree is historical and nothing writes to it. An
earlier version of this section told you to record the run id and the manifest content ids — do not
add either back, and do not write a script that reads one.

## What the typing lane hands these scripts

`data/development/typing/cells/` was rewritten on **2026-09-19** with the islet chain the
maintainer ruled on 2026-09-18 applied to all 20 donors — 25,230,109 cells, eight classes plus
`Excluded`. Measured 2026-09-24, every partition carries `islet_route` (cd99 / elongation / hormone)
and `islet_stratum`. There is one class per cell, reached in two steps (maintainer, 2026-10-19):
`cell_class` is the first step's column, and `refined_class` -- written by the refinement
(`core/phenocycler/scripts/typing/refine.py`) through the `build_cells.py --refine` join, which ran on
2026-09-24 -- holds the result. `Unresolved`, the ninth entry `run20_default.yml` declares, is a value
only the result reaches (1,063,097 cells). Every root script still reads `cell_class`, so none serves
the finished class yet; changing them is a code change not yet made.

`data/development/typing/subpanel/subpanel_calls.parquet` is the PREVIOUS generation (25,150,537
rows, written Sep 11). `senior/build_distance_profiles.py` refuses it by path with the counts in the
error. That tree stacks generations side by side, so check mtimes before reading anything else in
it.

**No `typing/app` artifact reads `data/development/typing/immune`.** That second pass
(`core/phenocycler/scripts/typing/immune_subtype.py`, rewritten 2026-09-17, 20 donors, 2,934,898
Immune cells) adds a finer subtype over the `Immune` class only and never modifies `cell_class`;
each app builder's coverage comes from its own `COHORT.resolve`. One layer does read it:
`senior/build_grid_composition.py` takes the fine subtype for `fine.csv` (:585), and
`build_layers.py` declares it as that node's input. It was last rebuilt 2026-09-26 over the
current Immune class (2,875,793 cells); `--check` reports its state under the table
(`SOURCE_CHECKS`) without turning the gate red.

## The layer builders

Nine builders write the **layer artifacts** the app's cards read: small, pre-aggregated files under
`$ISLET_LAYERS_DIR` (default `/data/islet_layers`). A card never scans a source tree at render
time; it reads one layer root and, if it is absent, renders a named note. The contract every one of
them is written against is [`docs/ISLET_LAYERS_SPEC.md`](../docs/ISLET_LAYERS_SPEC.md).

| layer root | builder |
|---|---|
| `availability/` | `senior/build_availability.py` |
| `sections/` | `senior/build_sections_layer.py` |
| `distance_profiles/` | `senior/build_distance_profiles.py` |
| `restain/` | `senior/build_restain_layer.py` |
| `islet_summary/` | `senior/build_islet_summary.py` |
| `trajectory/` | `senior/build_trajectory_layer.py` |
| `compartments/` | `senior/build_compartment_gradient.py` |
| `grid_composition/` | `senior/build_grid_composition.py` |
| `pln/` | `core/phenocycler/scripts/pln/build_pln_layer.py` — the one that does not live here |

`senior/derive_pln_node_membership.py` writes into `pln/` as well; it derives node membership
geometrically and is upstream of the layer, not a tenth layer.

**Do not invoke these by hand in an order you remember.** `senior/build_layers.py` declares the
graph — `(name, command, inputs, outputs)` per node, covering these nine and the `typing/app`
artifacts — and is the only thing that knows the order:

```bash
.venv/bin/python scripts/senior/build_layers.py --check   # stale/fresh table; non-zero if stale
.venv/bin/python scripts/senior/build_layers.py --plan    # the ordered commands, runs nothing
.venv/bin/python scripts/senior/build_layers.py --run     # executes in dependency order
```

`--check` compares each node's oldest output against its newest input, and **a builder's own source
file is one of its inputs**: a builder edited after its artifact was written means the artifact no
longer describes that builder. Five of the nine layers were stale against their own inputs on
2026-09-17 and nothing said so, which is what this driver exists to stop.

Each builder still runs standalone, and should, when you are working on one. What it must not do is
be the thing that decides the order.

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

`core/phenocycler/scripts/dct/` is a closed lane. One of its modules is still imported from outside
it: `regions.py` (the hash rule), by `broadtype/frame.py:34` — move it before deleting anything
there. The two symbols this line used to name, `gates.require_evidence` and
`tiff_window.read_window`, already live in `core/phenocycler/scripts/common/`.

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
