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
  + raw QuPath Cell/Nucleus means (REDSEA only for a marker routed to it; none since 2026-09-29)
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

The command-line pipeline and its stage APIs are the only writers below the
stable stage output roots under `data/`. Review notebooks and development
scripts must not:

- infer donors or panels from filenames;
- replace a required source with a similarly named raw, REDSEA, or historical
  RESTORE-normalized column;
- replace unavailable values with zero;
- create marker calls, cell types, or production manifests;
- overwrite another stage's output directory or mix outputs from different
  configurations.

## Stage graph

```text
QuPath cohort manifest
  -> ingest
  -> geometry
  -> redsea
  -> expression
```

Expression depends on geometry and, for a marker the registry routes to it, on
REDSEA. Since 2026-09-29 no marker is routed to REDSEA, so the `redsea` stage is a
logged no-op and expression reads only `cells` and `geometry_qc` (see
[Stage 3](#stage-3-redsea)). **The content-addressed
run framework that formerly gated these four stages was removed on
2026-09-16.** Each stage now writes directly to a stable root
(`data/{cells,geometry_qc,redsea,expression}`), keyed only by donor partition —
see [Stable-root execution](#stable-root-execution) below.

## Stage 0: lock the QuPath cohort

The cohort specification is the only place raw source paths are entered. The
manifest builder resolves and fingerprints:

- the exact donor and QuPath `Image` mapping;
- measurement CSVs and qptiffs;
- cell and nucleus geometry and UUID universes;
- panel and channel mapping;
- X/Y pixel calibration and segmentation version;
- the expected donor set.

Donor 6457 remains centrally excluded (`phenocycler/cohort.py` `DONOR_EXCLUSIONS`). A
manifest that lists it fails rather than silently filtering it later:
`ensure_eligible_donors` refuses the run (`pipeline.py:327`). 6579's exclusion was lifted
on 2026-09-17 (`REINSTATED`), and the manifest has declared it since 2026-09-23.

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

**Since 2026-09-29 the stage is a logged no-op** (`[no-op] redsea`, or
`[donor-noop] redsea/<id>` under `--pipelined`; `pipeline status` prints `NO-OP`):
every registered marker is `uncompensated`, so nothing reads its output. An audit of
13 donors found no marker whose correction removed biologically impossible
co-positives without also removing real single-positive signal: a median 5% fewer
co-positives, 0.91 of single positives kept, and a selectivity of about 1.3, close
to what dimming the channel gives. The evidence and the rule are in
[the data guide](../DATA_README.md#spillover-policy). The stage keeps its place in
the order and its code, and runs again once the registry routes a marker to it.

## Stage 4: authoritative expression

For each registered marker, expression selects exactly one declared source:

- the declared raw QuPath compartment for uncompensated markers, which since
  2026-09-29 is every marker: the bare `<marker>` column of `cells` for `Cell`
  and `Nucleus__<marker>` for `Nucleus`, named by ingest's own sanitiser; or
- the declared REDSEA compartment for a marker routed to REDSEA (none today;
  `qupath_norm` is refused since `a30a517`).

An acquired marker with a missing required source is fatal. A marker not
acquired in a panel becomes null with an explicit availability reason. Exact
QuPath-normalized values remain signed; finite negative values are not censored
or converted into calls.

The outputs are:

```text
expression/donor_id=<id>/*.parquet
audit/expression_availability.parquet
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
owned by the tool, and carry no production meaning. See
[`scripts/review/README.md`](../scripts/review/README.md).

## Stable-root execution

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

**The content-addressed run framework was removed on 2026-09-16.** There is no
stage manifest, no run ID, no donor receipt, and no fingerprint check against
the cohort, registry, configuration, or code. Each stage writes to its stable
root under `data/`; a donor-stage is done once its partition has non-empty
content, a present donor-stage is skipped, and `--force` reruns it. Runs made
before 2026-09-16 remain, untouched, under `data/runs/<run_id>/`; that tree is
historical and the 20-donor cohort output still lives only there.

`adoption.py`, `resume.py`, and `recovery.py` are deleted; there is no adopt,
resume, or recover command. They existed only to re-enter runs the removed
framework had locked, and there is no locked state left to re-enter — rerun
`pipeline run`, optionally with `--force`.

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
| donor/stage output looks stale after an upstream change | correct the input, then rerun with `--force` (no manifest or fingerprint check will catch this for you) |

Never tune a threshold, pair, donor exclusion, or cell-type frequency to make a
failure disappear.

## Historical material

`METHOD_HISTORY.md`, dated handoffs, and older notebooks preserve the reasoning
and measurements behind prior RESTORE, probabilistic, marker-debug, and Otsu
experiments. They are evidence about superseded development, not current
commands or scientific authority.
