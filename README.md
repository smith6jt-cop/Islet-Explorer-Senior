# Phenocycler analysis

This package turns an explicitly declared QuPath cohort into donor-local
RESTORE marker evidence, simultaneous broad cell identities, and UUID-keyed tables
for review in QuPath.

The production workflow is:

```text
QuPath cohort manifest
  -> measurement ingest
  -> geometry QC
  -> REDSEA boundary-spillover correction
  -> authoritative per-marker expression
  -> retained reference-control diagnostics
  -> RESTORE donor-local target/reference fits and tri-state marker calls
  -> simultaneous broad typing
  -> separate process-state annotation
  -> uncertainty-preserving QuPath export
```

The command-line pipeline, its production APIs, and its versioned registries
are authoritative. Notebooks are review-oriented views that may invoke those
same immutable stage APIs; they are not an alternative implementation and must
not introduce thresholds, donor exclusions, or typing rules of their own.

## What the workflow guarantees

- Every source image is declared in an immutable cohort manifest. Donors are
  never inferred from filenames.
- The QuPath `object_id` UUID is the cell key from ingest through export.
- Geometry QC precedes REDSEA and uses geometry only. Analysis eligibility,
  model-estimation eligibility, and physical spillover context are separate
  decisions. Spillover-ineligible polygons are removed before REDSEA contact
  construction; analysis-ineligible cells remain in the UUID universe but
  export as unavailable. An InstanSeg fragment without a usable nucleus is
  never eligible for analysis, estimation, or physical spillover context.
  Known duplicated islet detections use an explicit deterministic survivor
  policy; duplicate copies are excluded from all three roles.
- Each registered marker has one authoritative intensity source and
  compartment. An acquired marker with a missing required source is an error;
  a marker absent from a panel remains explicitly unavailable.
- RESTORE uses versioned mutually exclusive target/reference pairs and
  donor-local NNMF/control fits; it never uses reviewed cell labels.
- Marker evidence preserves corrected intensity, the donor-marker divisor,
  RESTORE-normalized intensity, a zero-centered log2 threshold ratio, pair-fit
  stability, and a tri-state call: `positive`, `negative`, or `unavailable`.
- Broad types are evaluated simultaneously, with no precedence between gates.
  Exactly one supported gate is an `anchor`; more than one is `Ambiguous` with
  the competing gates recorded; none, with everything explicitly negative, is
  `Other`; an unmeasured gate marker leaves the cell `Unavailable`.
- State markers such as Ki67 and PCNA are reported separately and cannot
  change cell identity.
- Every run and stage is content addressed. Partial, stale, or donor-incomplete
  output fails closed instead of being silently reused.

The detailed scientific and operational contract is in
[the workflow guide](docs/WORKFLOW.md). Input, output, and column contracts are
in [the data guide](DATA_README.md). The detached, image-reader-agnostic
[error-aware refinement workflow](docs/REFINEMENT_WORKFLOW.md) defines plotted
stage optimization, donor-held-out method comparison, and separate
false-positive/false-negative review without overwriting production labels.
The [probabilistic rescue guide](docs/PROBABILISTIC_RESCUE.md) documents the
implemented state-safe fitting, bundle, and opt-in inference contracts. No real
bundle is enabled or shipped until a new authorized donor-disjoint label set
passes that process.

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

The package defaults supply the marker registry and broad gate rules. A
project can pin explicit paths to reviewed copies:

```ini
marker_registry = phenocycler/marker_registry.json
broad_gate_rules = phenocycler/broad_gates.json
```

Run and inspect the complete workflow:

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

Useful bounded runs are:

```bash
# Run the dependency prefix through expression selection.
python -m phenocycler.pipeline run --config config.ini \
  --pipelined --through expression

# Run named stages only; their prerequisite manifests must already be current.
python -m phenocycler.pipeline run --config config.ini \
  --only controls calibrate type states

# Build all analytical stages without the QuPath handoff.
python -m phenocycler.pipeline run --config config.ini \
  --pipelined --no-export

# Export a current typing result later.
python -m phenocycler.pipeline export --config config.ini
```

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
An incompatible node and its DAG descendants are recomputed while independent
branches remain eligible. `type` and `states` are separate branches, so a
typing change does not unnecessarily invalidate `states`.

For source run `c8412757d12de36ca056`, the restoration policy adopts only
stages whose values are unchanged, then recomputes the RESTORE-dependent
suffix and its QuPath export. The generated plan is authoritative about the
exact reusable stages and byte count.
Run these commands from `Phenocycler_Analysis`, whose `config.ini` owns that
source run. To launch recomputation immediately instead of inspecting first,
add `--continue` (and optionally `--jobs N` or `--no-export`) to `adopt apply`.
After `--no-export`, run the normal `export` command separately. If code,
configuration, or source content changes, write a newly named plan rather than
overwriting the immutable reviewed plan. Until RESTORE typing finishes, `status`
correctly reports the adopted target as incomplete and `LATEST` is unchanged.
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

Run `status` immediately before a standalone export and proceed only when all
analytical stages are `CURRENT` (the QuPath line may still be `MISSING`). Run
it again after export and require success before importing the CSVs into
QuPath. This validates the typing artifact and fingerprinted handoff together.

`status` checks the QuPath handoff as well as analytical stages. After an
intentional `--no-export` run it reports that export as missing and returns
nonzero until `export` is run.

The donor-local dependency graph is fixed:

```text
ingest (cohort barrier)
  └─ geometry[d] -> redsea[d] -> expression[d] -> controls[d]
                                      │              └─ RESTORE[d] -> type[d]
                                      └─ states[d]
```

`--only` does not infer or rebuild dependencies. This makes every requested
transition explicit and prevents an interactive notebook from quietly
substituting a different upstream result. Pipelined execution accepts
dependency prefixes through `--through`; it does not accept `--only`.

## Review notebooks

The notebooks follow the same four conceptual checkpoints:

1. [`01_qupath_ingest_and_geometry_qc.ipynb`](notebooks/01_qupath_ingest_and_geometry_qc.ipynb)
   reviews the cohort contract, ingest, and geometry decisions.
2. [`02_redsea_spillover.ipynb`](notebooks/02_redsea_spillover.ipynb)
   reviews REDSEA alignment, compartments, and authoritative expression.
3. [`03_marker_calibration.ipynb`](notebooks/03_marker_calibration.ipynb)
   (legacy filename) reviews RESTORE pair states, divisors, feature sets,
   stability, call fractions, and selected-pair evidence disposition.
4. [`04_hierarchical_cell_typing.ipynb`](notebooks/04_hierarchical_cell_typing.ipynb)
   reviews hierarchy, state, and the optional QuPath handoff.

[`05_error_aware_refinement.ipynb`](notebooks/05_error_aware_refinement.ipynb)
is a read-only companion over an unpromoted validation
`CandidateEvaluationManifest`. It uses projected queries to build a weighted
probability sample and a separately capped full-validation-split challenge
sample, then
plots FP/FN, reliability, and the calibration-frozen coverage/error point only
when complete detached adjudication and valid score semantics are available.
It does not require QuPath and never edits a production assignment. The
checked-in notebook reviews the broad target; each fitted subtype requires its
own blinded sample, review ledger, and structured held-out report.

[`06_probabilistic_validation_and_promotion.ipynb`](notebooks/06_probabilistic_validation_and_promotion.ipynb)
is the guarded release interface. It gates the complete validation report
before an explicit one-time locked-test acknowledgement, plots weighted
coverage and FP/FN confidence bounds for every required target, replays the
calibration and held-out source lineage, and signs the separate promotion
manifest with protected external keys. Its checked-in flags keep every write
and locked-test read disabled.

Optional execution cells in the earlier stage notebooks call `run_stage` or
`export_qupath`, which apply the same manifest validation and immutability
checks as the CLI. The RESTORE review notebook is deliberately read-only: while
the donor queue is active, it validates completed donor receipt trees and
reviews their donor-local audits and evidence without waiting for the cohort
RESTORE stage manifest. Notebook cells must not write stage files directly.

## RESTORE marker evidence

RESTORE is performed independently for each donor and active marker. The
marker registry declares the target's mutually exclusive reference. Raw
REDSEA's Nucleus/Cytoplasm/Membrane compartment decomposition supplies the
NNMF feature basis; the reference-positive/target-negative population defines a
donor-local divisor on the REDSEA target measurement. The divisor statistic is
configurable (`[restore] divisor_statistic`) and is an OPEN question: measured
over 3 donors x 24 pairs, no single statistic agrees with the pairs' own NNMF
anchors across all markers and donors. See `docs/METHOD_HISTORY.md`.

The key outputs answer different questions:

| Value | Interpretation |
|---|---|
| authoritative intensity | selected corrected signal in the marker's declared compartment |
| `log2_threshold_ratio` | signed distance from the donor-local threshold; zero is the threshold |
| `restore_normalized` | REDSEA target divided by the donor-marker RESTORE divisor |
| pair audit | target/reference, NNMF stability, feature set, controls, divisor, and disposition |
| expression state | positive, negative, or unavailable |

This representation addresses donor-to-donor scale and background while
keeping the full fit auditable. Technical failure is unavailable, never a
negative fallback.

Annotations and existing QuPath classifications are audit context only. They
are forbidden as RESTORE inputs and cannot define marker calls or
identity.

## Identity and state

Broad identity evaluates every gate in `phenocycler/broad_gates.json`
**simultaneously**. No gate has precedence over another:

```text
Immune | Endothelial | Epithelial | Neural | Mesenchymal
```

A cell supported by exactly one gate takes that identity. A cell supported by
more than one is `Ambiguous` — it is never awarded to whichever gate is listed
first. A cell with no supporting gate is `Other` only when every gate marker was
explicitly negative; if any gate marker was unmeasured it is `Unavailable`.
Geometry-QC-ineligible cells are also `Unavailable`.

The production default is broad-only and label-free; compatibility columns
such as `specific_type` contain the accepted broad label. The parent-constrained
probabilistic subtype implementation is retained only as an explicit
signed-bundle option and is dormant while `[typing] model_bundle` is blank.

CD11b is an active RESTORE target and subtype evidence, but it does **not**
gate broad Immune: measured over 20 donors its EpCAM control is 17.3% median
target-positive, so the pair's mutual-exclusivity premise fails for broad
identity. See `docs/METHOD_HISTORY.md`.

Process-state calls are stored in a separate `cell_states` artifact. They do not
rewrite broad or specific identity.

## Runs and failure behavior

The run identifier is derived from the cohort content ID, marker-registry
fingerprint, broad-gate rules fingerprint, identity-affecting configuration, and
relevant source-code hash. Outputs live below:

```text
data/runs/<run_id>/
```

Each stage writes a manifest containing its exact inputs, configuration, source
hash, expected donors, output schema, row counts, object-ID universe, and
fingerprints for its audit/availability sidecars. A complete analytical run
writes `run.json`; `data/runs/LATEST` records its run ID. QuPath export has its
own manifest because it is a separate handoff.

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

- `expression/`: one authoritative intensity per registered marker;
- `reference_controls/`: donor/reference membership and deterministic ranks;
- `marker_evidence/`: per-cell RESTORE calls, normalized values, and fit provenance;
- `assignments/`: simultaneous broad identity, supporting gates, and explicit
  unavailability;
- `cell_states/`: process-state annotations kept separate from identity;
- `qupath_export/`: UUID-keyed identity CSVs that retain uncertainty;
- `audit/`: expression availability, control diagnostics, RESTORE pair models,
  and state models;
- `manifests/`: immutable stage and export contracts.

See [DATA_README.md](DATA_README.md#run-layout) for the full tree and schemas.
Use
[`scripts/groovy/import_cell_assignments.groovy`](scripts/groovy/import_cell_assignments.groovy)
for a validation-first UUID import that keeps broad, specific, and uncertainty
decisions separate. Run `status` first: the importer validates every CSV row
but intentionally leaves QuPath detections outside the canonical ingest
universe (for example, fragments below the structural area floor) unchanged.

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

**Current status.** The integration layer's Xenium half is unaffected. Its
PhenoCycler half still reads the legacy RESTORE + broad-lineage layout (`phenotype/broad/`,
`restore_gated_redsea/`) and has **not** yet been repointed at this workflow's `assignments/`
and `marker_evidence/` outputs, so `export_pheno` and the vocabulary crosswalk are stale
against the new compartments. The config surface those 22 modules read is preserved (see the
quarantined block in `phenocycler/config.py`) so the package imports and the suite collects,
but treat the PhenoCycler side of the integration as pending that rewire. Nothing in the
eight-stage core pipeline reads any of it.

## Restoration note

RESTORE is the production marker-calling path. Historical loose RESTORE files
are not copied into a content-addressed run; the current stage refits from
fingerprinted raw compartments and REDSEA. For the failed run
`c8412757d12de36ca056`, reuse valid upstream work and recompute RESTORE plus
typing as described in
[`docs/REFINEMENT_WORKFLOW.md`](docs/REFINEMENT_WORKFLOW.md).

## Scope

This repository begins with QuPath measurement, geometry, and qptiff outputs.
Image registration, illumination correction, segmentation-model development,
and downstream spatial-neighborhood or trajectory analyses are separate
workflows. The goal here is a reproducible per-cell evidence layer that supports
eventual assignment of every cell while preserving honest ambiguity whenever
the panel or RESTORE evidence cannot support a specific call.
