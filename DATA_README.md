# Data contracts

The workflow treats raw QuPath exports, biological registries, and every
derived dataset as explicit contracts. A directory's existence is never proof
that a stage completed.

For the scientific rationale and stage algorithms, see
[the workflow guide](docs/WORKFLOW.md). For installation and commands, see
[the project README](README.md).

## Source files

The current cohort contract contains exactly one image record per donor.
Different image records may point to the same batched measurement CSV.

| Source | Required contract |
|---|---|
| QuPath measurement CSV(s) | Exact `Image` and `Object ID` values; centroids, parent, cell area/solidity, cell marker means, and nucleus means for registered nuclear/state markers |
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
```

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
| Marker measurable, donor-marker RESTORE fit technically invalid | measurement retained | state `unavailable`, never negative |
| Marker measurable, valid RESTORE fit | normalized value and divisor retained | positive iff normalized value is greater than 1 |

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
- active, deferred, or standalone use;
- whole-cell or nucleus compartment;
- REDSEA-corrected, REDSEA passthrough, or uncompensated policy;
- incompatible reference, marker family, and broad/subtype use.

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
or absence.

## Reference-control schema

`reference_controls` stores donor-local masks in a wide donor partition. The
canonical long representation records:

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

## Marker-evidence schema

`marker_evidence` is wide: one row per UUID and a repeated set of columns for
each active marker `<m>`.

| Column | Meaning |
|---|---|
| `qc_analysis_eligible` | frozen geometry-QC analysis decision propagated to typing |
| `<m>__corrected_intensity` | REDSEA whole-cell target measurement used for the call |
| `<m>__state` | `positive`, `negative`, or `unavailable` |
| `<m>__restore_normalized` | corrected intensity divided by the donor-marker RESTORE divisor |
| `<m>__restore_positive` | exact boolean `restore_normalized > 1` call |
| `<m>__log2_threshold_ratio` | `log2(restore_normalized)`; zero at the donor-local divisor |
| `<m>__model_valid` | valid RESTORE pair and positive/negative cell call |
| `<m>__calibration_status` | RESTORE pair disposition (historical column name retained for schema compatibility) |
| `<m>__reference` | mutually exclusive reference marker used by RESTORE |
| `<m>__threshold` | donor-marker RESTORE divisor on the REDSEA intensity scale |
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
divisor, feature-set, and call-count details are written to
`audit/restore_pair_models.parquet`.

## Assignment schema

The production-default `assignments` table records the simultaneous broad
typing decision. Stable handoff columns are:

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

Rules, registry, method, and categorical-score provenance are repeated on each
row. Probability/ranking and subtype columns remain in the schema so the
signed optional probabilistic lane can use the same handoff contract; they are
empty or not applicable in the production default. Default statuses are:

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

## State schema

`cell_states` is keyed by donor and UUID but remains separate from
`assignments`. The current state stage reports Ki67 and PCNA:

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

## QuPath export schema

The export stage writes one immutable CSV per donor/image:

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

`compartment` and `cell_type` are compatibility aliases for `broad_type` and
`specific_type`; they do not replace the uncertainty fields. Unavailable cells
must not be imported as confident named classes. The current
QuPath handoff exports identity; process-state annotations remain in the
separate `cell_states` artifact.

Use
[`scripts/groovy/import_cell_assignments.groovy`](scripts/groovy/import_cell_assignments.groovy).
First require a successful pipeline `status`, which verifies the fingerprinted
export. The script then validates the exact image and every row in the
canonical assignment CSV before changing a QuPath class. Detections outside
the canonical ingest universe, including fragments below the ingest area
floor, are intentionally left unchanged. At the selected broad or specific
level, accepted statuses (`anchor` and `other`) receive the corresponding
label; `unavailable` receives an explicit uncertainty class rather than the
best candidate. The importer also understands the additional statuses used by
an explicitly configured probabilistic lane.

## Run layout

Given a base `data_dir = data`, the content-addressed layout is:

```text
data/
├── qupath_manifest.json
└── runs/
    ├── LATEST
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
        │   ├── controls.json
        │   ├── calibrate.json
        │   ├── type.json
        │   ├── states.json
        │   ├── donors/
        │   │   └── donor_id=<id>/
        │   │       ├── geometry.json
        │   │       ├── redsea.json
        │   │       ├── expression.json
        │   │       ├── controls.json
        │   │       ├── calibrate.json
        │   │       ├── type.json
        │   │       └── states.json
        │   └── qupath_export.json                  # when exported
        ├── cells/donor_id=<id>/*.parquet
        ├── ingest_rejects.parquet
        ├── panel_availability.json
        ├── geometry_qc/donor_id=<id>/*.parquet
        ├── redsea/donor_id=<id>/*.parquet
        ├── expression/donor_id=<id>/*.parquet
        ├── reference_controls/donor_id=<id>/*.parquet
        ├── marker_evidence/donor_id=<id>/*.parquet
        ├── assignments/donor_id=<id>/*.parquet
        ├── cell_states/donor_id=<id>/*.parquet
        ├── qupath_export/cell_assignments_<id>.csv # when exported
        ├── audit/
        │   ├── donors/<audit_name>/donor_id=<id>/data.parquet
        │   ├── expression_availability.parquet
        │   ├── reference_control_models.parquet
        │   ├── restore_pair_models.parquet
        │   └── state_models.parquet
        └── scratch/redsea/
```

`run.json` exists only after all analytical stage manifests are current.
`runs/LATEST` contains the completed run ID. The QuPath export manifest is
separate and may be created after the analytical run.

## Manifest contents and immutability

Each stage manifest records:

- stage and method version;
- cohort, registry, rules, and configuration identities;
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
