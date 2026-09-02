# Current QuPath-to-expression workflow

This document is the operational and scientific contract for the current
Phenocycler workflow. It produces an auditable, UUID-keyed expression table
from declared QuPath inputs.

> **Development-only boundary:** no target/reference marker pair is currently
> verified for production. Ordinary execution stops after `expression`.
> Pair-dependent `controls`, `calibrate`, `type`, and assignment `export` are
> blocked. The pair-independent `states` stage remains explicit-only.
> Historical RESTORE, probabilistic, marker-debug, and
> Otsu paths do not provide fallbacks.

Start with [the project README](../README.md) for commands and
[the data guide](../DATA_README.md) for schemas. For the state of every step
on 2026-09-02 — this workflow, the development lanes, the Xenium lane, the
integration layer and the app — see
[`docs/WORKFLOW_2026-09-02_MAP.md`](../../../docs/WORKFLOW_2026-09-02_MAP.md).

## Intended outcome

The current product is one authoritative upstream measurement per registered
marker and cell:

```text
declared QuPath cohort
  + geometry-only eligibility
  + compartment-resolved REDSEA measurements
  + exact source selection and availability audit
  -> donor-partitioned expression
```

It is not a marker-call, phenotype, or cell-type dataset. Missing measurements
remain missing, signed QuPath-normalized values remain signed, and no expected
prevalence or disease ordering is used.

## Authorities and boundaries

Three versioned inputs govern the current workflow:

1. The QuPath cohort manifest owns raw-data identity, donor membership, image
   mapping, geometry, panel, channel, pixel calibration, and segmentation
   version.
2. `marker_registry.json` owns marker names, compartments, measurement sources,
   and spillover policy. Historical interpretation fields do not authorize a
   current pair.
3. `config.ini` owns fixed geometry, REDSEA, and compute parameters.

The command-line pipeline and its stage APIs are the only writers below a
content-addressed run. Review notebooks and development scripts must not:

- infer donors or panels from filenames;
- replace a required source with a similarly named raw, REDSEA, or historical
  RESTORE-normalized column;
- replace unavailable values with zero;
- create marker calls, cell types, or production manifests;
- overwrite a run/stage manifest or mix outputs from different runs.

## Stage graph

```text
QuPath cohort manifest
  -> ingest
  -> geometry
  -> redsea
  -> expression
```

Expression depends on geometry and REDSEA because the registry contains raw,
REDSEA-routed, and exact QuPath-normalized measurements. The run/stage manifest
architecture remains in place for all four stages.

## Stage 0: lock the QuPath cohort

The cohort specification is the only place raw source paths are entered. The
manifest builder resolves and fingerprints:

- the exact donor and QuPath `Image` mapping;
- measurement CSVs and qptiffs;
- cell and nucleus geometry and UUID universes;
- panel and channel mapping;
- X/Y pixel calibration and segmentation version;
- the expected donor set.

Donors 6457 and 6579 remain centrally excluded. A manifest that contains either
fails rather than silently filtering it later.

## Stage 1: ingest the measurement union

Ingest scans the declared CSVs, filters by exact image allow-list, assigns donor
IDs from the manifest, preserves QuPath UUIDs, and unions panel-specific columns
by name. It does not retain QuPath's signed `Mean (norm)` classifier fields:
`cells/` carries no `QuPathNorm__` column (verified 2026-08-27, `a30a517`), and
the registry refuses `measurement_source = qupath_norm`.

## Stage 2: geometry-only QC

Geometry QC uses polygons, areas, solidity, centroids, and nucleus/cell
relationships. It does not use marker intensity, annotation, classification,
disease status, or expected cell-type frequency.

The stage keeps distinct flags for:

- analysis eligibility;
- estimation eligibility; and
- physical REDSEA spillover context.

Every source UUID remains auditable. Missing nucleus geometry, prohibited
duplicates, and invalid cell geometry fail or exclude according to the declared
geometry policy.

## Stage 3: REDSEA

REDSEA uses the qptiff channel map and eligible cell boundaries to correct
lateral spillover independently for each channel. It emits aligned
Nucleus/Cytoplasm/Membrane/Cell measurements and preserves the exact UUID
universe. REDSEA does not assign positivity or cell identity.

## Stage 4: authoritative expression

For each registered marker, expression selects exactly one declared source:

- the declared REDSEA compartment for REDSEA-routed markers (Iba1, CD68 and
  M2Gal3 included — `qupath_norm` is refused since `a30a517`); or
- the declared raw QuPath compartment for uncompensated markers.

An acquired marker with a missing required source is fatal. A marker not
acquired in a panel becomes null with an explicit availability reason. Exact
QuPath-normalized values remain signed; finite negative values are not censored
or converted into calls.

The outputs are:

```text
expression/donor_id=<id>/*.parquet
audit/expression_availability.parquet
manifests/expression.json
```

## Mutable development diagnostics

The only current pairwise roster is
`phenocycler/fundamental_pair_diagnostics.json`. Its shipped default contains
two upstream comparisons:

1. raw `Cell: CD3e: Mean` versus raw `Cell: Ker8-18: Mean`;
2. exact signed `Cell: Iba1: Mean (norm)` versus raw
   `Cell: Ker8-18: Mean`.

CD79a is excluded as a pair candidate because its expression is too sparse in
normal tissue.

The Iba1/Ker8_18 comparison mixes measurement scales deliberately and therefore
uses donor-local ranks. It never applies a common numeric boundary. The tools
report coverage, correlations, and deterministic exact-rank co-high-tail
overlap only.

Development outputs live under `../data/development`, are replaceable when
owned by the tool, and carry no run/stage manifest or production meaning. See
[`scripts/review/README.md`](../scripts/review/README.md).

## Content-addressed execution

Run the upstream workflow from the package directory:

```bash
python -m phenocycler.pipeline run --config config.ini --pipelined
python -m phenocycler.pipeline status --config config.ini
```

The explicit bounded form is:

```bash
python -m phenocycler.pipeline run --config config.ini \
  --pipelined --through expression
```

Do not request pair-dependent stages or assignment export. A future production
workflow requires a new reviewed registry and code contract; it cannot be
enabled by treating historical artifacts as current.

Each current stage manifest records its method, inputs, scientific
configuration, code identity, donor set, schema, row and UUID universe, output
fingerprints, and owned audit sidecars. Donor receipts permit pipelined work
without weakening cohort-level completion checks.

Adoption, resume, and recovery may reuse only verified-compatible upstream
artifacts. Historical reference controls, marker evidence, assignments, state
calls, and exports are not adopted into the current workflow.

## Failure contract

| Failure | Required response |
|---|---|
| source fingerprint changed | verify the intended source and recreate the cohort manifest |
| unmapped `Image` value | correct the exact image mapping or source export |
| CSV parser reject | inspect `ingest_rejects.parquet` and repair the source |
| missing/duplicate geometry UUID | re-export or correct segmentation/geometry |
| acquired source column missing | correct the registry, channel map, or upstream export; do not substitute another scale |
| REDSEA alignment failure | correct UUID/channel/raster contracts |
| request for controls/calls/types/export | stop: no verified production pairs exist |
| stale or partial manifest | repair the declared input and create/continue the correct content-addressed run |

Never tune a threshold, pair, donor exclusion, or cell-type frequency to make a
failure disappear.

## Historical material

`METHOD_HISTORY.md`, dated handoffs, and older notebooks preserve the reasoning
and measurements behind prior RESTORE, probabilistic, marker-debug, and Otsu
experiments. They are evidence about superseded development, not current
commands or scientific authority.
