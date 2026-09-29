# Phenocycler analysis

This package currently turns an explicitly declared QuPath cohort into
UUID-keyed, source-audited per-cell expression data. There are **no verified
production target/reference pairs** in the current workflow, so ordinary
execution stops after expression and does not create marker calls, cell types,
or a QuPath assignment export.

The current workflow is:

```text
QuPath cohort manifest
  -> measurement ingest
  -> geometry QC
  -> REDSEA boundary-spillover correction   (a logged no-op: no marker reads it since 2026-09-29)
  -> authoritative per-marker expression
```

Stages write to stable roots under `data/` (`cells/`, `geometry_qc/`, `redsea/`,
`expression/`); a donor is done when its partition has content, a present
stage is skipped, and `--force` reruns it. There is no run manifest, no
run ID, and no fingerprint gate on these upstream artifacts — see
[Outputs](#outputs). Pair-dependent `controls`,
`calibrate`, `type`, and assignment `export` are blocked until
target/reference pairs are independently verified. Historical RESTORE,
probabilistic, marker-debug, and Otsu material is retained only as method
history; it is not a current command path or fallback.

## What the workflow guarantees

- Every source image is declared in the cohort manifest
  (`qupath_manifest.json`), which records a content ID and per-file SHA-256 as
  a description; nothing validates them, so the file can be hand-edited when a
  source moves. Donors are never inferred from filenames.
- The QuPath `object_id` UUID is the cell key throughout the upstream workflow.
- Geometry QC precedes REDSEA and uses geometry only. Analysis eligibility,
  model-estimation eligibility, and physical spillover context are separate
  decisions. Spillover-ineligible polygons are removed before REDSEA contact
  construction; analysis-ineligible cells remain in the UUID universe. An
  InstanSeg fragment without a usable nucleus is
  never eligible for analysis, estimation, or physical spillover context.
  Known duplicated islet detections use an explicit deterministic survivor
  policy; duplicate copies are excluded from all three roles.
- Each registered marker has one authoritative intensity source and
  compartment. An acquired marker with a missing required source is an error;
  a marker absent from a panel remains explicitly unavailable.
- Every marker reads one source: the raw QuPath mean in its declared
  compartment (`uncompensated`), or the declared REDSEA compartment for a
  marker routed to REDSEA. Since 2026-09-29 all 72 markers are
  `uncompensated`: a 13-donor audit found no marker whose REDSEA correction
  removed biologically impossible co-positives without also removing real
  signal, so the REDSEA stage is a logged no-op and expression reads only
  `cells` (the evidence is in [the data guide](DATA_README.md#spillover-policy)).
  Since `a30a517` (2026-08-27) no marker may select `qupath_norm` — the signed
  QuPath `Mean (norm)` scale was a RESTORE-era exception and `cells/` carries no
  `QuPathNorm__` column.
- The only active pairwise work is a mutable, rank-only development diagnostic:
  raw CD3e versus raw Ker8_18 and exact QuPath `Iba1_norm` versus raw Ker8_18.
  CD79a is excluded as a pair candidate because its expression is too sparse in
  normal tissue. These comparisons make no positivity call and write no
  production artifact.
- Each stage writes to a stable, non-content-addressed root. Existing output
  for a donor/stage is reused as-is; there is no fingerprint check against the
  source data, registry, or code, so a stale input is not detected — rerun
  with `--force` after any upstream change.

The detailed scientific and operational contract is in
[the workflow guide](docs/WORKFLOW.md). Input, output, and column contracts are
in [the data guide](DATA_README.md). The
[method history](docs/METHOD_HISTORY.md) records superseded pair-calling work
without authorizing it for current use.

## Install

Python 3.10 or newer is required. In an activated analysis environment:

```bash
python -m pip install -e .
```

For development:

```bash
python -m pip install -e '.[test]'
pytest -q
```

The reference `environment.yml` includes Jupyter and Matplotlib for the review
notebooks. In another environment, install them with
`python -m pip install -e '.[notebook]'` before opening a notebook.

## 1. Export the QuPath inputs

For every image, retain:

- the QuPath cell-measurement CSV;
- the matching `.qptiff`;
- full-resolution cell and nucleus geometry;
- the exact value in the CSV `Image` column;
- pixel calibration, panel identity, channel mapping, and segmentation version.

Use [`scripts/groovy/export_cells_geojson.groovy`](scripts/groovy/export_cells_geojson.groovy)
to export geometry with the QuPath detection UUID intact. Nucleus geometry may
be embedded as `nucleusGeometry` in the cell GeoJSON or supplied as a separate
UUID-matched GeoJSON. New exports must contain a usable nucleus for every cell.
Historical cell-only fragments may remain in a source export for provenance,
but geometry QC rejects them from analysis, model estimation, and REDSEA
spillover context. Always pass the intended output directory as the first
script argument (`--args "<output-dir>"`); the exporter refuses to run without
it.

## 2. Create the cohort manifest

Generate a versioned specification template:

```bash
python -m phenocycler.pipeline manifest template --out cohort-spec.json
```

Edit the specification, then fingerprint and validate its sources:

```bash
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json
```

Measurement CSVs receive a full SHA-256 fingerprint by default. Large qptiffs
use an explicitly recorded deterministic sampled fingerprint by default; use
`--qptiff-hash-mode full` when full-file hashing is practical.

The specification schema and a complete example are in
[DATA_README.md](DATA_README.md#cohort-specification).

## 3. Configure and run

At minimum, the configuration identifies the data root and cohort manifest:

```ini
[paths]
data_dir = data
qupath_manifest = data/qupath_manifest.json
```

The package default supplies the marker registry. A project can pin an
explicit path to a reviewed copy:

```ini
marker_registry = phenocycler/marker_registry.json
```

Run and inspect the current upstream workflow:

```bash
python -m phenocycler.pipeline run --config config.ini --pipelined
python -m phenocycler.pipeline status --config config.ini
```

`--pipelined` keeps ingest as a cohort barrier, then schedules independent
donor work as soon as its own prerequisites are complete. For example, REDSEA
for one donor can run while geometry rasterization continues for another.
Each stage writes directly to its stable output root
(`data/{cells,geometry_qc,redsea,expression}`); a donor-stage is considered
done once its partition has non-empty content, with no receipt or manifest
step in between.

The operational limits live in `[compute]`:

```ini
n_jobs = 4
geometry_workers = 2
redsea_workers = 1
downstream_workers = 4
```

`n_jobs` is the total number of simultaneously running tasks. The other values
are group-specific ceilings within that total. Geometry, REDSEA, and
downstream work use separate process pools, preventing a REDSEA GPU context
from migrating across workers. With `use_gpu = true`, REDSEA is always limited
to one task on the configured GPU.

The ordinary run currently ends at expression. The equivalent bounded command
is:

```bash
python -m phenocycler.pipeline run --config config.ini \
  --pipelined --through expression
```

Do not invoke `controls`, `calibrate`, `type`, or assignment export. They
require a verified marker-pair workflow, which does not currently exist. The
orthogonal `states` stage is pair-independent and may be requested explicitly
after expression with `--only states`; it is not part of the ordinary run.

There is no adoption, resume, or recovery command any more —
`adoption.py`, `resume.py`, and `recovery.py` were deleted on 2026-09-16 along
with the run framework they served (see [Outputs](#outputs)). Granularity
depends on which invocation you rerun: with `--pipelined` (the invocation
above), each donor-stage is checked separately, so rerunning picks up any
donor still missing output for that stage without redoing donors already
written. Without `--pipelined`, the check is whole-stage: `run_stage` skips a
stage the moment its output root has *any* content, so a partially-written
stage missing one donor is not topped up automatically — `--force` reruns
the whole stage (or donor-stage) whose output is already present.

The donor-local dependency graph is fixed:

```text
ingest (cohort barrier)
  └─ geometry[d] -> redsea[d] -> expression[d]
```

While no registry marker reads REDSEA (the case since 2026-09-29), `redsea[d]` is
a logged no-op (`[donor-noop] redsea/<d>`) that computes and writes nothing, so
`expression[d]` follows `geometry[d]` without a REDSEA partition.

`--only` does not infer or rebuild dependencies. This makes every requested
transition explicit and prevents an interactive notebook from quietly
substituting a different upstream result. Pipelined execution accepts
dependency prefixes through `--through`; it does not accept `--only`.

## Current development diagnostics

The active review lane works directly from upstream QuPath measurements and is
mutable by design. Its shipped roster contains only:

1. raw `Cell: CD3e: Mean` versus raw `Cell: Ker8-18: Mean`; and
2. exact signed `Cell: Iba1: Mean (norm)` versus raw
   `Cell: Ker8-18: Mean`.

CD79a is not a candidate pair marker because its expression is too sparse in
normal tissue.

The mixed Iba1/Ker8_18 comparison is rank-only; the two marker scales are not
numerically combined. The tools report coverage, association, and exact
co-high-tail overlap without creating a threshold, positive call, cell type,
or production artifact. See
[`scripts/review/README.md`](scripts/review/README.md) for commands.

The first two notebooks remain useful for upstream inspection. Calibration,
typing, refinement, probabilistic-promotion, marker-debug, and Otsu notebooks
or scripts are historical and do not define a current execution path.

## The typing development lane (outside this package's outputs)

The statement above — that the workflow creates no cell types — is about **this package's
pipeline**, which is latched. A separate development lane under
[`scripts/typing/`](scripts/typing/) does classify cells, and it writes **outside this tree**, to
the repo-root `data/development/typing/`. It is not a release, it does not move
`NO_VERIFIED_TYPING`, and nothing in `data/` here reads it.

- What it does and how it was recovered: [`scripts/typing/RECONSTRUCTION.md`](scripts/typing/RECONSTRUCTION.md).
  Every statement in that file is marked `[ON DISK]` or `[CHAIN NOW]`; read the mark before quoting
  a number.
- What is on disk: `data/development/typing/cells/` — rewritten on **2026-09-19**, 25,230,109 cells,
  20 donors, eight classes plus `Excluded`, with the islet chain the maintainer ruled on 2026-09-18
  (structure, then exocrine-inside-an-islet excluded, then CD99 at a level fitted per donor **and per
  size stratum**, with nucleus elongation as a second route in the large stratum only) applied to
  every donor. Measured 2026-09-24: every partition carries `islet_route` and `islet_stratum`.
- `data/development/typing/subpanel/subpanel_calls.parquet` is the PREVIOUS generation (25,150,537
  rows) and is not current. The refinement pass
  ([`scripts/typing/REFINEMENT.md`](scripts/typing/REFINEMENT.md), `refine.py`) was ruled in on
  2026-09-23 and is wired as the `build_cells.py --refine` join, `refined_class` beside `cell_class`;
  the join has not run on disk, and whether any builder should read `refined_class` is open.
- Reference levels are fitted **per donor**. Pooling was measured wrong on 2026-09-15; the four
  pooled levels survive only so `subpanel.py --gate` can reproduce the older table.

## Outputs

**The content-addressed run framework was removed on 2026-09-16.** Each
upstream stage now writes straight to a stable root under `data/`:

```text
data/cells/
data/geometry_qc/
data/redsea/
data/expression/
data/audit/expression_availability.parquet
```

A donor is done when its partition (`donor_id=<donor>/`) has content; a
present stage is skipped on the next run, and `--force` reruns it. There is no
`run_id`, no stage manifest, no donor receipt, and no fingerprint check
against the cohort, registry, configuration, or code — a stale or moved input
is not detected, so correct the declared input yourself and pass `--force`
where output already exists. `resume`, `recover`, and `adopt` are gone as
commands; rerunning `run` (optionally with `--force`) is the only re-entry
path now that there is no locked run state to re-enter.

Runs made before 2026-09-16 remain, untouched, under
`data/runs/<run_id>/manifests/` and the corresponding stage directories; that
tree is historical and nothing writes to it going forward, but the full
20-donor cohort output still lives only there. Donors 6533 and `0112` were added
on 2026-09-16 and exist only in the flat `data/<stage>/` tree above (the partitions are
`donor_id=6533` and `donor_id=0112`).

See [DATA_README.md](DATA_README.md#run-layout) for the full tree and schemas.

## PhenoCycler ↔ Xenium integration

`phenocycler/integration/` integrates this repository's PhenoCycler output with 10x Xenium
spatial transcriptomics, in two modes — `sequential` (serial sections from one block, where
the cells are *not* the same cells, so pairing is at structure / niche / donor level) and
`same_slide` (one physical section imaged twice). It is an optional extra with its own
environment:

```bash
# On the workstation the repo .venv is the environment; the conda recipe is HPG-only.
../../.venv/bin/pip install -e '.[integration]'
../../.venv/bin/python -m phenocycler.integration.manifest --report
```

Design and data caveats: **[docs/INTEGRATION.md](docs/INTEGRATION.md)**.
Running the whole thing on HiPerGator: **[docs/HIPERGATOR.md](docs/HIPERGATOR.md)**.

**Current status (verified 2026-09-03; run framework removed 2026-09-16; release
provenance changed to coverage-checked on `a1cd264`).** Both modalities can be
exported as geometry-profile cell tables. The PhenoCycler exporter reads the
`expression` stage output plus same-run cell area; the Xenium importer ignores
embedded label-shaped columns. Manifest, geometry export, Xenium import,
raster/image registration, geometry QC, and geometry figures are live.
Identity-dependent stages accept only exact native-key calls from configured
release directories; with those paths blank they fail before reading caches.
Releases are **not** identity-pinned any more — the old `run_id` +
`content_id` check compared against a stage manifest, and stage manifests
were removed on 2026-09-16. What runs now is a coverage check: does the
declared expression input resolve, and does it contain every donor the
release claims. It refuses only a release that claims a donor its own
declared input lacks; an input that does not resolve at all is recorded as
`unverifiable:<reason>` on `PhenotypeRelease.provenance` rather than refused,
and `integration --status` prints that value per modality. Nothing in the
current upstream or integration workflow reads legacy phenotype output.

## Historical pair-calling note

Prior RESTORE, probabilistic, marker-debug, and Otsu artifacts remain historical
evidence only. They are not copied into a current run and do not authorize a
target/reference pair, marker call, identity, or export.

## Scope

This repository begins with QuPath measurement, geometry, and qptiff outputs.
Image registration, illumination correction, segmentation-model development,
and downstream spatial-neighborhood or trajectory analyses are separate
workflows. The current goal is a reproducible upstream expression layer and a
small mutable diagnostic surface from which future pair hypotheses can be
evaluated without being mistaken for production calls.
