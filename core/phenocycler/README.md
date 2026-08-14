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
  -> REDSEA boundary-spillover correction
  -> authoritative per-marker expression
```

The run/stage manifest architecture remains authoritative for these upstream
artifacts. Pair-dependent `controls`, `calibrate`, `type`, and assignment
`export` are blocked until target/reference pairs are independently verified.
Historical RESTORE, probabilistic, marker-debug, and Otsu material is retained
only as method history; it is not a current command path or fallback.

## What the workflow guarantees

- Every source image is declared in an immutable cohort manifest. Donors are
  never inferred from filenames.
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
- Iba1, CD68, and M2Gal3 use QuPath's exact signed `Mean (norm)` measurements
  in Cell, Nucleus, Cytoplasm, and Membrane compartments. Those measurements
  were created upstream for the QuPath object classifier. Raw, REDSEA, and
  historical RESTORE-normalized values are forbidden substitutes.
  Their signed scale is `(raw - tissue p50)/(tissue p99 - tissue p50)`, so 1 is
  the tissue p99 anchor rather than a cell-positivity threshold.
- Signed QuPath norm values are preserved unchanged through ingest and
  expression, including finite negative values and NaN.
- The only active pairwise work is a mutable, rank-only development diagnostic:
  raw CD3e versus raw Ker8_18 and exact QuPath `Iba1_norm` versus raw Ker8_18.
  CD79a is excluded as a pair candidate because its expression is too sparse in
  normal tissue. These comparisons make no positivity call and write no
  production artifact.
- Every run and stage is content addressed. Partial, stale, or donor-incomplete
  output fails closed instead of being silently reused.

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
Every donor-stage output is atomically written and receives an immutable
receipt before it can unlock downstream work. The ordinary cohort stage
manifest is sealed after all receipts for that stage are current.

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

To carry forward the maximum compatible work from any prior run with completed
stage manifests, first build an adoption plan and review every stage decision.
Planning does not modify either run; `--out` writes only the requested immutable
plan file:

```bash
python -m phenocycler.pipeline adopt plan \
  --config config.ini \
  --from-run <completed_run_id> \
  --out adoption-plan.json

python -m phenocycler.pipeline adopt apply \
  --config config.ini \
  --plan adoption-plan.json

# Inspect manifests/adoption/, then launch only the missing work.
python -m phenocycler.pipeline run --config config.ini --pipelined
```

Adoption follows the stage dependency graph. It compares the cohort and direct
inputs, stage-specific scientific configuration, method version, schema,
donor/UUID universe, recorded file fingerprints, and producing-code contract.
Identical code is accepted directly; code drift is accepted only by an exact,
checked-in compatibility record—there is no `--force` switch. Accepted Parquet
datasets and scientific audit sidecars are verified hard links, then receive
new target manifests under the current compatibility contract. Immutable
receipts retain the original producer manifests and current dependency IDs.
An incompatible node and its DAG descendants are recomputed. Current adoption
is limited to the upstream stages through expression; historical pair-dependent
artifacts are not adopted into the current workflow. The generated plan is
authoritative about the exact reusable stages and byte count.
Run these commands from `Phenocycler_Analysis`, whose `config.ini` owns that
source run. To launch recomputation immediately instead of inspecting first,
add `--continue` (and optionally `--jobs N`) to `adopt apply`. If code,
configuration, or source content changes, write a newly named plan rather than
overwriting the immutable reviewed plan.
Use `--validation-mode content` on `adopt plan` when an archival audit should
rehash every source byte rather than use the normal fingerprint cache.

If geometry stopped at the first duplicated-islet donor after completing a
duplicate-free prefix under `duplicate_policy = fatal`, resume the corrected
`deterministic_keep` run without rebuilding ingest or rerasterizing that
prefix:

```bash
python -m phenocycler.pipeline resume \
  --config config.ini \
  --from-run <failed_run_id> \
  --from-donor <first_unfinished_donor>
```

The resume command verifies the source ingest manifest, code, cohort, filter,
bytes, donor/UUID universes, geometry policy transition, and absence of
duplicate groups in every reused geometry donor. It hard-links ingest files
to avoid another data copy. The completed geometry prefix is value-preserving
re-encoded only to give its optional duplicate columns the same Parquet types
as later duplicate-bearing donors. Resume receipts are included in the new
stage manifests. The reusable geometry prefix immediately becomes eligible
for REDSEA while geometry resumes at the unfinished donor. By default the
command continues through all later stages; add `--through geometry` to stop
after sealing the resumed geometry stage.

If a run started before donor receipts were available and was stopped only
after sealing geometry, recover its manifested geometry and any fully written
REDSEA donor files into the receipt queue with:

```bash
python -m phenocycler.pipeline recover \
  --config config.ini \
  --from-run <stopped_run_id>
```

Recovery does not accept arbitrary partial files. It validates the source
ingest and geometry manifests, scientific configuration, code lineage, donor
set, row/UUID universes, and geometry content. Every surviving REDSEA Parquet
must pass a complete all-column scan and match the current schema, donor,
ingest row count, and UUID universe. Accepted files are hard-linked into a new
content-addressed run and receive donor receipts; scratch files are never
adopted. Use `--prepare-only` to create and inspect the verified handoff
without starting the queue.

The donor-local dependency graph is fixed:

```text
ingest (cohort barrier)
  └─ geometry[d] -> redsea[d] -> expression[d]
```

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
or immutable run artifact. See
[`scripts/review/README.md`](scripts/review/README.md) for commands.

The first two notebooks remain useful for upstream inspection. Calibration,
typing, refinement, probabilistic-promotion, marker-debug, and Otsu notebooks
or scripts are historical and do not define a current execution path.

## Runs and failure behavior

The run identifier is derived from the cohort content ID, marker-registry
fingerprint, upstream configuration, and relevant source-code hash. Outputs
live below:

```text
data/runs/<run_id>/
```

Each upstream stage writes a manifest containing its exact inputs,
configuration, source hash, expected donors, output schema, row counts,
object-ID universe, and fingerprints for its audit/availability sidecars. The
manifest architecture is preserved even though the current workflow stops at
expression.

Existing valid content is reused. The pipeline stops when it finds any of the
following:

- source data, registry, rules, configuration, or code changed;
- a donor or expected partition is missing;
- an output schema or UUID universe changed;
- an audit or availability sidecar changed or disappeared;
- an acquired marker's required source column is absent;
- geometry is missing or object IDs are duplicated;
- a non-empty output exists without a valid manifest.

Do not repair these failures by deleting manifests or mixing files between run
directories. Correct the declared input or configuration and create the new
content-addressed run.

## Outputs

The principal artifacts under a run are:

- `cells/`: canonical UUID-keyed QuPath measurements;
- `geometry_qc/`: geometry-only analysis and estimation eligibility;
- `redsea/`: compartment-resolved spillover-corrected measurements;
- `expression/`: one authoritative intensity per registered marker;
- `audit/expression_availability.parquet`: source selection and availability;
- `manifests/`: immutable upstream stage contracts.

See [DATA_README.md](DATA_README.md#run-layout) for the full tree and schemas.

## PhenoCycler ↔ Xenium integration

`phenocycler/integration/` integrates this repository's PhenoCycler output with 10x Xenium
spatial transcriptomics, in two modes — `sequential` (serial sections from one block, where
the cells are *not* the same cells, so pairing is at structure / niche / donor level) and
`same_slide` (one physical section imaged twice). It is an optional extra with its own
environment:

```bash
conda env create -f environment-integration.yml
conda activate phenocycler_integration && pip install -e '.[integration]'
python -m phenocycler.integration.manifest --report
```

Design and data caveats: **[docs/INTEGRATION.md](docs/INTEGRATION.md)**.
Running the whole thing on HiPerGator: **[docs/HIPERGATOR.md](docs/HIPERGATOR.md)**.

**Current status.** The integration layer's Xenium-only import and manifest reporting remain
available. Every Pheno-dependent orchestrator and public stage entry point fails before
reading cached historical `phenotype/broad/`, `restore_gated_redsea/`, or `cells_pheno`
outputs; the retired vocabulary crosswalk is blocked too. Treat the PhenoCycler side as
pending a newly verified pair/typing workflow and an explicit rewire. Nothing in the current
upstream pipeline reads legacy phenotype output.

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
