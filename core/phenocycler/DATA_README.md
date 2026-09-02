# Data contracts

The workflow treats raw QuPath exports, biological registries, and every
derived dataset as explicit contracts. A directory's existence is never proof
that a stage completed.

**Current boundary:** there are no verified production target/reference pairs.
The ordinary workflow stops after `expression`. Pair-dependent
`reference_controls`, `marker_evidence`, `assignments`, and QuPath assignment
export are not current outputs. The run/stage manifest architecture is retained
for ingest, geometry, REDSEA, expression, and an explicit pair-independent
`states` stage.

For the scientific rationale and stage algorithms, see
[the workflow guide](docs/WORKFLOW.md). For installation and commands, see
[the project README](README.md).

## Source files

The current cohort contract contains exactly one image record per donor.
Different image records may point to the same batched measurement CSV.

| Source | Required contract |
|---|---|
| QuPath measurement CSV(s) | Exact `Image` and `Object ID` values; centroids, parent, cell area/solidity, cell marker means, nucleus means for registered nuclear/state markers, and exact Cell/Cytoplasm/Membrane/Nucleus `Mean (norm)` fields for Iba1, CD68, and M2Gal3 |
| qptiff | Original channel stack used for REDSEA; channel order and names must agree with the declared map |
| cell GeoJSON | One feature per detection, with feature `id` equal to the QuPath UUID and a usable cell geometry |
| nucleus geometry | Either `nucleusGeometry` embedded in the cell GeoJSON or a separate UUID-matched nucleus GeoJSON; every QC-accepted cell must have one |
| image metadata | Positive X/Y pixel size, panel ID, channel map, and a non-empty segmentation version; the current REDSEA path requires isotropic X/Y pixels |

The canonical key is the QuPath `Object ID`, represented downstream as
`object_id`. Centroid matching and filename parsing are not substitutes for
this UUID.

The measurement CSV is expected to contain these QuPath fields:

```text
Image
Object ID
Centroid X µm
Centroid Y µm
Parent
Cell: Area µm^2
Cell: Solidity
Cell: <marker>: Mean
Nucleus: <marker>: Mean       # where required by the registry
<compartment>: <marker>: Mean (norm)
                               # compartment = Cell, Cytoplasm, Membrane, Nucleus
                               # marker = Iba1, CD68, M2Gal3; exact fields required
```

The three `Mean (norm)` marker families were written by the original QuPath
object-classifier workflow before RESTORE, using
`(raw - tissue p50)/(tissue p99 - tissue p50)`. Thus 0 is the tissue p50 anchor,
1 is the tissue p99 anchor rather than a positivity call, and negative values
are valid. Ingest retains them as
`QuPathNorm__<Compartment>__<Marker>`. They are not RESTORE-normalized outputs,
and raw marker means, REDSEA measurements, `*_norm` convenience columns, and
`*__restore_normalized` columns may not substitute for a missing exact field.

The ingest stage derives `cell_region` and `islet_num` from `Parent`, but these
context fields do not define calibration controls or cell identity. Other
QuPath morphology, classification, and distance columns may remain in the raw
CSV but are intentionally not copied into the canonical cell table. Nucleus
and cell raster areas used by geometry QC come from the contracted GeoJSON.

Use
[`scripts/groovy/export_cells_geojson.groovy`](scripts/groovy/export_cells_geojson.groovy)
for UUID-preserving geometry export. Always pass the intended output directory
as its first script argument (`--args "<output-dir>"`); the exporter refuses
to run without it.

## Cohort specification

Create a starter specification:

```bash
python -m phenocycler.pipeline manifest template --out cohort-spec.json
```

Paths are resolved relative to the specification file. A minimal specification
is:

```json
{
  "spec_version": 1,
  "cohort_name": "islet-phenocycler",
  "expected_donors": ["D001"],
  "defaults": {
    "measurement_csv": [
      "raw/CellMeasurementsBatch1.csv",
      "raw/CellMeasurementsBatch2.csv"
    ],
    "segmentation_version": "qupath-project-commit-or-model-version"
  },
  "images": [
    {
      "donor_id": "D001",
      "image_id": "exact QuPath Image column value",
      "measurement_csv_index": 0,
      "qptiff": "raw/D001.qptiff",
      "cell_geojson": "raw/D001_cells.geojson",
      "nucleus_geojson": null,
      "pixel_size_um_x": 0.325,
      "pixel_size_um_y": 0.325,
      "panel_id": "panel-v1",
      "channel_map": null
    }
  ]
}
```

`measurement_csv` may be one path or a list of batched exports. When it is a
list, every image declares the zero-based `measurement_csv_index` containing
that image. The manifest fingerprints each unique export once, and ingest
combines all declared exports with `UNION ALL BY NAME` so panel-specific
measurement columns remain available.

The cohort manifest's exact `image_id` values are an ingest allow-list. Rows
for other images in a shared export are ignored before donor assignment. The
CSV files must still be scanned to locate those rows because CSV has no index,
but only contracted cohort rows are written to the cell Parquet dataset.

`nucleus_geojson: null` declares that nucleus polygons are embedded as
`nucleusGeometry`; a path declares a separate nucleus FeatureCollection.
Historical combined exports may also contain InstanSeg fragments without the
embedded member. Their absence is fingerprinted in the manifest and those
objects receive `missing_nucleus_geometry`; they can never be analysis- or
estimation-eligible and their polygons are removed from REDSEA spillover
context.
`channel_map: null` asks the manifest builder to derive canonical marker names
from the qptiff channel metadata. An explicit map uses:

```json
[
  {
    "marker": "DAPI",
    "channel_name": "DAPI",
    "channel_index": 0
  },
  {
    "marker": "E_cadherin",
    "channel_name": "E-cadherin",
    "channel_index": 1
  }
]
```

The array must enumerate every qptiff channel in contiguous zero-based index
order. `channel_name` must exactly match qptiff metadata; `marker` is the
canonical registry name. Every marker, channel name, and channel index must be
unique within an image. Images with the same `panel_id` must have equivalent
mappings. Different panel IDs may intentionally differ.

Build the immutable cohort manifest:

```bash
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json
```

The builder validates geometry and panel consistency, captures source
fingerprints, and assigns a cohort `content_id`. The default hash modes are:

| Source | Default | Meaning |
|---|---|---|
| measurement CSV | `full` | SHA-256 over every byte |
| qptiff | `sampled` | deterministic windows plus size and sampling layout; explicitly weaker than a full hash |
| geometry | `full` | exact file fingerprint plus UUID/geometry summary |

Select stronger or weaker modes explicitly when creating the manifest:

```bash
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json \
  --measurement-hash-mode full \
  --qptiff-hash-mode full
```

Do not hand-edit the generated cohort manifest. Edit the small specification
and recreate it so its identity can be verified.

## Panel availability

Ingest forms the union, not the intersection, of marker columns across declared
measurement exports. A marker that was not acquired for a panel is retained as
null for that donor. `panel_availability.json` records each image, panel, and
observed marker set.

Availability has distinct meanings:

| Condition | Representation | Downstream consequence |
|---|---|---|
| Marker not acquired in donor's panel | null authoritative value; `marker_not_acquired_in_panel` audit reason | per-cell marker state `unavailable` |
| Marker acquired, required declared source present, cell value missing | null for that cell | that cell is `unavailable` for the marker |
| Marker acquired, required declared source column absent | fatal contract error | stage stops |
| Marker acquired and its exact declared source is present | signed authoritative value retained | available for upstream review; no current marker call |

Missing markers never become negatives, never make a cell `Other`, and never
transfer their fixed error allocation to another marker.

## Ingested cell schema

Ingest writes one donor-partitioned Parquet dataset and preserves the union of
panel marker columns. Core columns are:

| Column | Meaning |
|---|---|
| `donor_id` | exact donor ID from the cohort manifest |
| `image` | exact QuPath `Image` value |
| `object_id` | canonical QuPath UUID |
| `X_centroid`, `Y_centroid` | centroid in micrometres |
| `parent_raw` | QuPath parent context retained for audit only |
| `cell_region`, `islet_num` | context derived from `Parent` |
| `cell_area`, `cell_solidity` | structural measurements retained from the CSV |
| `<marker>` | sanitized whole-cell mean |
| `Nucleus__<marker>` | retained nucleus mean where available |

The ingest floor removes only degenerate objects below
`[cells] min_cell_area` (5 µm² by default). Biological and geometry decisions
are not made at ingest. CSV parser rejects are written to
`ingest_rejects.parquet`; any rejected row makes ingest fail.

## Geometry QC schema

Geometry QC is an immutable per-cell decision made before REDSEA. Its principal
fields are:

| Column | Meaning |
|---|---|
| `analysis_eligible` | cell may receive an accepted identity; ineligible cells remain present but type as unavailable |
| `estimation_eligible` | cell may help estimate donor-local control models |
| `spillover_context_eligible` | polygon may provide physical REDSEA adjacency context; false labels are removed before contact construction |
| reason fields | explicit failed rules and audit context |
| raster/cell area ratio | agreement between rasterized polygon and QuPath area |
| nucleus/cell area ratio | identifies cells with no meaningful cytoplasm |
| duplicate metadata | UUID/centroid/area duplicate checks |

Analysis eligibility and estimation eligibility are intentionally different. A
real cell with no measurable cytoplasm, for example, may remain in the analysis
universe but must not estimate a cytoplasmic background model. Missing geometry
never produces an accepted cell: a missing nucleus excludes that object, while
a missing canonical cell polygon or prohibited duplicate remains fatal.

## Authoritative expression

[`phenocycler/marker_registry.json`](phenocycler/marker_registry.json) declares
the only accepted source for each marker:

- biological kind: identity marker or process-state marker;
- current availability/development status;
- whole-cell or nucleus compartment;
- measurement source: REDSEA or raw QuPath (`qupath_norm` is refused since `a30a517`);
- REDSEA-corrected, REDSEA passthrough, or uncompensated policy;
- historical interpretation metadata, which does not authorize a production pair.

> **2026-08-27 (`a30a517`):** `qupath_norm` is no longer a permitted
> `measurement_source` — `validate_registry` refuses it and all 60 markers
> resolve to `redsea` or `qupath_raw`. The RESTORE-era text that follows is
> historical.

Measurement source and spillover policy are separate contracts. Iba1, CD68,
and M2Gal3 had `measurement_source=qupath_norm`; their exact signed QuPath
measurements win even though their spillover policy remains recorded. Missing
exact measurements fail closed rather than falling back to raw, REDSEA, or a
historical RESTORE output. None of these markers currently receives a call.

The `expression` dataset contains exactly one selected column named
`<marker>` for every registered marker, plus:

```text
donor_id
object_id
image
X_centroid
Y_centroid
cell_region
islet_num
qc_analysis_eligible
qc_estimation_eligible
```

Keeping one authoritative value avoids ambiguous mixtures of raw,
whole-cell, nuclear, and corrected columns downstream.
`audit/expression_availability.parquet` explains every donor-marker selection
or absence, including its `measurement_source`. QuPath norm values remain
signed and unchanged in `cells` and `expression`; finite negative values are
not censored by the current workflow.

## Historical reference-control schema (blocked)

Historical runs may contain donor-local masks in `reference_controls`. The
current workflow has no verified references and does not create this dataset.
Its former canonical long representation recorded:

```text
mask_schema_version
method_version
registry_version
registry_fingerprint
source
donor_id
reference
object_id
reference_positive
seed_rank
```

Controls are selected from finite, positive reference-marker values among
estimation-eligible cells with finite coordinates. Ranking uses only the
reference intensity, deterministic spatial round-robin selection, and UUID
tie-breaking. Target intensity, target prevalence, QuPath `Classification`,
phenotype labels, and annotation-derived groups are forbidden inputs.

`audit/reference_control_models.parquet` records candidate counts, selected
counts, status, and provenance.

## Historical marker-evidence schema (blocked)

Historical runs may contain a wide `marker_evidence` table. The current
workflow does not create or consume it. Its former columns included one row per
UUID and a repeated set for each then-active marker `<m>`.

| Column | Meaning |
|---|---|
| `qc_analysis_eligible` | frozen geometry-QC analysis decision propagated to typing |
| `<m>__selected_intensity` | registry-selected target measurement on RESTORE's nonnegative input scale |
| `<m>__state` | `positive`, `negative`, or `unavailable` |
| `<m>__restore_normalized` | selected intensity divided by the donor-marker RESTORE divisor |
| `<m>__restore_positive` | exact boolean `restore_normalized > 1` call |
| `<m>__log2_threshold_ratio` | `log2(restore_normalized)`; zero at the donor-local divisor |
| `<m>__empirical_tail_probability` | finite-sample upper-tail probability against the clean seeded control |
| `<m>__censored_below` | finite selected value is zero after the source-declared RESTORE boundary transform, so its call is unavailable |
| `<m>__model_valid` | valid RESTORE pair and positive/negative cell call |
| `<m>__calibration_status` | RESTORE pair disposition (historical column name retained for schema compatibility) |
| `<m>__reference` | mutually exclusive reference marker used by RESTORE |
| `<m>__threshold` | donor-marker RESTORE divisor on the selected-intensity scale |
| `<m>__n_controls` | reference-positive/target-negative control count |
| `<m>__method_version` | pipeline RESTORE method version |

RESTORE pair dispositions include:

```text
VALID
NO_TARGET_EXPRESSION_CONFIRMED
LOW_SIGNAL_UNRESOLVED
NO_REFERENCE_POPULATION
INVALID_PAIR
MODEL_UNSTABLE
PAIR_REVIEW_PENDING
TECHNICAL_FAILURE
```

An invalid donor-marker fit does not erase its measurements or fabricate a
negative. Calls become `unavailable`. Full NNMF, arm, stability, control,
divisor, feature-set, call-count, measurement-source, boundary-transform, and
left-censor-count details are written to
`audit/restore_pair_models.parquet`.

For `measurement_source=qupath_norm`, finite signed negatives are explicitly
left-censored to zero exactly when RESTORE builds its nonnegative input and the
audit transform is `qupath_norm_left_censor_zero`. NaN remains NaN. The same
source-generic `censored_below` state also covers a selected measurement already
at zero; zero is unavailable rather than a confident biological negative.

## Historical assignment schema (blocked)

Historical runs may contain an `assignments` table. There is no current
production typing decision or stable assignment handoff. Former columns were:

| Column | Meaning |
|---|---|
| `broad_type` | the single supported gate, `Ambiguous`, `Other`, or `Unavailable` |
| `specific_type` | compatibility alias of `broad_type` in the broad-only default |
| `assignment_status` | `anchor`, `ambiguous`, `other`, or `unavailable` |
| `confidence` | null — this method estimates no probability |
| `best_broad_type` | the anchored gate, or the competing/blocking gate set |
| `best_specific_type` | compatibility alias of the broad result |
| `assignment_reason` | `single_supported_gate`, `multiple_supported_gates`, `all_negative`, `unmeasured_gate_marker`, or `geometry_qc_ineligible` |
| `blocked_by` | gates that could not be evaluated, `\|`-joined |
| `broad_supported_gates` | every supported (or blocking) gate, `\|`-joined |
| `broad_supported_gate_n` | how many gates support the cell |
| `broad_support` | RESTORE-normalized support of the best gate (1.0 is the threshold) |
| `broad_support_margin` | best minus second-best support, for ambiguous cells |
| `typing_mode` | `restore_simultaneous` |

Rules, registry, method, and categorical-score provenance were repeated on
each historical row. Former statuses were:

```text
anchor
ambiguous
other
unavailable
```

`Other` means every gate marker was explicitly negative. It does not mean
“could not classify” — that is `Ambiguous` (contradictory evidence) or
`Unavailable` (missing evidence), and the three are deliberately distinct.

`broad_best_probability`, `broad_margin`, `broad_stability` and `confidence`
are `NaN`: they were never estimated, and the same rows carry
`typing_broad_bootstrap_replicates = 0`. Use `broad_support` /
`broad_support_margin` for evidence strength; they are reported, not
thresholded.

## Optional pair-independent state schema

An explicitly requested `states` stage writes `cell_states` keyed by donor and
UUID, separately from assignments. It is omitted from the ordinary expression
run but does not depend on a target/reference pair. Columns are:

| Column | Meaning |
|---|---|
| `qc_analysis_eligible` | frozen geometry-QC decision controlling whether a state may be called |
| `<m>__value` | authoritative state-marker value |
| `<m>__state` | `positive`, `negative`, or `unavailable` |
| `<m>__threshold` | donor-local state threshold |
| `<m>__model_status` | `valid`, `marker_unavailable`, or `insufficient_measurements` |
| `<m>__method` | versioned state method |

These calls describe proliferation state and never alter broad or specific
identity. The threshold is fitted only on estimation-eligible cells; a state is
called only for analysis-eligible cells, and geometry-ineligible cells remain
`unavailable`. State thresholds do not claim the same frequency-neutral
mutually-exclusive-control interpretation as identity-marker calibration.

## Historical QuPath assignment export (blocked)

There is no current assignment export. Historical export stages wrote one CSV
per donor/image:

```text
qupath_export/cell_assignments_<donor>.csv
```

Columns are:

```text
object_id
image
broad_type
specific_type
assignment_status
confidence
best_broad_type
best_specific_type
broad_assignment_status
specific_assignment_status
ranked_broad_probabilities
ranked_specific_probabilities
ranked_probabilities
compartment
cell_type
```

`compartment` and `cell_type` were compatibility aliases for `broad_type` and
`specific_type`; they did not replace the uncertainty fields.

The old importer is retained for historical artifacts only. Do not run it as a
current workflow step.

## Run layout

Given a base `data_dir = data`, the content-addressed layout is:

```text
data/
├── qupath_manifest.json
└── runs/
    └── <run_id>/
        ├── run.json
        ├── manifests/
        │   ├── ingest.json
        │   ├── ingest_reuse.json                  # only for verified reuse
        │   ├── geometry.json
        │   ├── geometry_resume.json               # only for verified resume
        │   ├── geometry_recovery.json             # only for verified recovery
        │   ├── redsea_recovery.json               # only for verified recovery
        │   ├── redsea.json
        │   ├── expression.json
        │   ├── donors/
        │   │   └── donor_id=<id>/
        │   │       ├── geometry.json
        │   │       ├── redsea.json
        │   │       └── expression.json
        ├── cells/donor_id=<id>/*.parquet
        ├── ingest_rejects.parquet
        ├── panel_availability.json
        ├── geometry_qc/donor_id=<id>/*.parquet
        ├── redsea/donor_id=<id>/*.parquet
        ├── expression/donor_id=<id>/*.parquet
        ├── audit/
        │   ├── donors/<audit_name>/donor_id=<id>/data.parquet
        │   └── expression_availability.parquet
        └── scratch/redsea/
```

`run.json` records completion of the current upstream stage set. Historical
downstream directories may exist in older runs but are not current outputs and
must not be mixed into this layout.

## Manifest contents and immutability

Each stage manifest records:

- stage and method version;
- cohort, registry, and upstream configuration identities;
- relevant source-code hash and repository metadata;
- exact upstream artifact identifiers;
- exact expected donor set;
- output Parquet schema, row counts, file fingerprints, and UUID-universe
  hashes;
- exact fingerprints for stage-owned availability and audit sidecars.

With donor-pipelined execution, each post-ingest partition first receives a
donor-stage receipt under `manifests/donors/`. That receipt carries the same
config/code/input/output contracts at donor scope and owns its donor-local
audit file. Once every donor receipt for a stage is current, the pipeline
concatenates the donor audits, seals the ordinary cohort stage manifest, and
includes every donor receipt among that manifest's fingerprinted sidecars.
Thus donor receipts permit cross-stage scheduling without weakening the
cohort-level archival contract.

A donor partition is immutable within a run. In pipelined mode, a non-empty
donor partition without its valid donor receipt is partial and fails; current
receipts remain valid while the cohort stage manifest is still pending. In
stage mode, non-empty output without the stage manifest remains partial.
Changed inputs, configuration, code, schemas, donors, output bytes, receipt
bytes, or sidecar bytes make the owning contract stale and also fail. Correct
the source contract and run again; the changed content resolves to a new
`run_id`. Cross-run reuse is limited to the specialized `pipeline resume` and
`pipeline recover` contracts: they verify their documented transitions,
record reuse receipts, and seal new stage manifests under the new run ID.

## Storage, privacy, and retention

Multiplexed images and per-cell artifacts can be large. Keep `data/` outside
version control, retain the raw source files and generated cohort manifest
together, and archive the registries/rules used by published analyses. For
human tissue, follow the governing IRB, de-identification, access-control, and
retention requirements.
