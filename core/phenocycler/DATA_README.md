# Data contracts

The workflow treats raw QuPath exports, biological registries, and every
derived dataset as explicit contracts. A directory's existence is never proof
that a stage completed.

**Current boundary:** there are no verified production target/reference pairs.
The ordinary workflow stops after `expression`. Pair-dependent
`reference_controls`, `marker_evidence`, `assignments`, and QuPath assignment
export are not current outputs. Ingest, geometry, REDSEA, expression, and the
explicit pair-independent `states` stage each write to a stable root under
`data/`; see [Run layout](#run-layout) — **the content-addressed run framework
that used to gate these stages was removed on 2026-09-16.** REDSEA writes nothing
while no marker reads it, which has been the case since 2026-09-29 (see
[Spillover policy](#spillover-policy)).

**Adoption audit, 2026-09-02 (historical — the `adopt` command it exercised was
deleted on 2026-09-16).** A dry run of
`phenocycler adopt plan --from-run 30f6689e4d27d295f890` proposed zero adoptable
stages: ingest, geometry, REDSEA, expression, and states required recomputation;
controls, calibrate, and type remained blocked by `no_verified_typing`. No
`adopt apply` was run. At the time, lanes pinned inputs by run/content hash (the
source expression content ID was
`35fa90fa66887329a3fb90ce23c0c944ede2fa73594fb1caac744f1bc8931ea0`); that
content-ID pinning no longer exists. `pipeline status` now reports what each
stage's stable output root holds; a present, non-empty partition is the only
evidence a stage completed.

For the scientific rationale and stage algorithms, see
[the workflow guide](docs/WORKFLOW.md). For installation and commands, see
[the project README](README.md).

## Source files

The current cohort contract contains exactly one image record per donor.
Different image records may point to the same batched measurement CSV.

| Source | Required contract |
|---|---|
| QuPath measurement export(s) | v2: batch CSV/TSV. v3: per image, a cells Parquet plus a context Parquet (see [v3 Parquet pairs](#v3-parquet-pairs)). Exact `Image` and `Object ID` values; centroids, parent, cell area/solidity, cell marker means, and nucleus means for registered nuclear/state markers. No `Mean (norm)` field is required (see below) |
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

No `Mean (norm)` field is required. Ingest requires and keeps one only for a marker whose
registry `measurement_source` is `qupath_norm` (`cells_parquet.QUPATH_NORMALIZED_MARKERS`,
`cells_parquet.py:57-61`), and since 2026-08-27 (a30a517) no marker is: Iba1, CD68 and M2Gal3
took the `redsea` default like every other marker then, and since 2026-09-29 read the raw
QuPath mean like every other marker (`marker_registry.py:372-390`; see
[Spillover policy](#spillover-policy)).

The `Mean (norm)` families the batch exports still carry (24 columns: Iba1, CD68,
M2Gal3, B3TUBB, CD34 and Podoplanin across four compartments; the 0112 and 6533
revised exports carry none) were written by the original
QuPath object-classifier workflow before RESTORE, using
`(raw - tissue p50)/(tissue p99 - tissue p50)`. Thus 0 is the tissue p50 anchor,
1 is the tissue p99 anchor rather than a positivity call, and negative values
are valid. A marker contracted to `qupath_norm` would be kept as
`QuPathNorm__<Compartment>__<Marker>`, with no substitution for a missing exact
field; no marker is contracted to it today, so ingest keeps none.

The ingest stage derives `cell_region`, `islet_num`, `cluster_num` and
`single_num` from `Parent`, but these context fields do not define calibration
controls or cell identity. `cell_region` is `core` (`Islet_<n>`), `peri`
(`Islet_<id>_exp20um`), `cluster` (`Cluster_<n>`), `single` (`Single_<n>`),
`pLN` (`Annotation (pLN)`), `artifact` (`Annotation (Artifact)`), `tissue`
(`Annotation (Tissue)` or `Root object …`), else `other`. The v2 exports carry
only the core / peri / tissue parents, so their mapping is unchanged except that
an `Annotation (pLN)` or `Annotation (Artifact)` parent, previously `other`, now
has its own region. Other
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

### v3 Parquet pairs

The v3 re-export writes two Parquet files per image. `results/v3/cells/<tag>.parquet`
carries every detection measurement under QuPath's own column names; it is written
before islets are parented, so its `Parent` is not authoritative, and it has no
`(norm)` columns. `results/v3/context/<tag>.parquet` carries `Object ID`, `Parent`,
both centroids and every `Signed distance to annotation with <Class> µm` column
(QuPath 0.7 spells it with "with"; the v2 exports did not, and readers accept both).

Declare the pair on the image entry: `measurement_csv` names the cells Parquet and
`measurement_context` the context Parquet. Both are fingerprinted; a v2 entry omits
`measurement_context` and serialises exactly as before.

```json
{
  "donor_id": "6374",
  "image_id": "exact QuPath Image column value",
  "measurement_csv": "results/v3/cells/6374.parquet",
  "measurement_context": "results/v3/context/6374.parquet",
  "...": "qptiff, geojson, calibration and panel as above"
}
```

Ingest joins the pair on `Object ID` in DuckDB, taking `Parent`, centroids and signed
distances from the context file. It refuses the pair before writing anything if the two
`Object ID` sets differ (reporting both counts and a sample of each difference) or if
either file repeats an `Object ID`.

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

Build the cohort manifest:

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

`content_id` and the per-file SHA-256 values are a description, not an
enforced contract: nothing downstream validates them against the manifest
recreated for a hand edit, so `qupath_manifest.json` can be hand-edited when a
source file moves. Prefer editing the small specification and recreating the
manifest so the recorded fingerprints stay accurate, but a direct edit is not
rejected.

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
| `cell_region`, `islet_num`, `cluster_num`, `single_num` | context derived from `Parent` (`''` where the Parent is not that kind) |
| `cell_area`, `cell_solidity` | structural measurements retained from the CSV |
| `<marker>` | sanitized whole-cell mean (`Cell: <marker>: Mean`); the raw source `expression` reads for a `Cell` marker |
| `Nucleus__<marker>` | retained nucleus mean where available; the raw source for a `Nucleus` marker |

`<marker>` is `cells_parquet.sanitize` of the QuPath name: every run of characters
other than letters and digits becomes one `_`. Registry names follow the qptiff
channel spelling instead, which agrees for every marker but one:
`Histone_H3_(p_Ser28)` is stored here as `Histone_H3_p_Ser28` and
`Nucleus__Histone_H3_p_Ser28`, and `expression` looks it up under that name.

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
- REDSEA-corrected, REDSEA passthrough, or uncompensated policy — every marker
  is `uncompensated` since 2026-09-29 (see [Spillover policy](#spillover-policy));
- historical interpretation metadata, which does not authorize a production pair.

> **2026-08-27 (`a30a517`):** `qupath_norm` is no longer a permitted
> `measurement_source` — `validate_registry` refuses it and every marker
> resolves to `redsea` or `qupath_raw` (all 72 to `qupath_raw` since
> 2026-09-29). The RESTORE-era text that follows is historical.

Measurement source and spillover policy are separate contracts, which
`validate_registry` ties together: `uncompensated` selects `qupath_raw`, and the
two REDSEA policies select `redsea`. Iba1, CD68, and M2Gal3 had
`measurement_source=qupath_norm` until 2026-08-27 (a30a517), then took `redsea`
under their `redsea_corrected` policy; since 2026-09-29 they read the raw QuPath
mean like every other marker. None of these markers currently receives a call.

### Spillover policy

**Since 2026-09-29 all 72 registered markers are `uncompensated`**: each reads
QuPath's raw mean in its declared compartment from `cells`, and nothing in the
ordinary workflow reads REDSEA. Before that, 63 markers were `redsea_corrected`
(62 whole-cell markers and the nuclear `Histone_H3_(p_Ser28)`), `Collagen_IV` was
`redsea_passthrough`, and the eight nuclear markers (`DAPI`, `FOXP3`, `Ki67`, `PCNA`,
`SOX2`, `TCF_1`, `TOX`, `TP63`) were already `uncompensated`.

The approved rule keeps REDSEA for a marker only where it reduces biologically
impossible co-positives (cells positive for two markers of mutually exclusive
lineages, such as INS/GCG, CD3e/CD20, CD4/CD8 inside CD3e+ cells, or CD3e/E-cadherin)
without destroying the marker's real single-positive signal. **Evidence.** An audit
of 13 donors scored every corrected marker against each exclusive partner, correcting
the marker alone with the partner held raw, under the typing lane's own positivity
rule (Otsu where that rule's threshold was unstable between raw and corrected values).
Over 748 donor × marker × partner comparisons, correction lowered the impossible
co-positive rate by a median 5% (IQR 1–11%), kept a median 0.91 (IQR 0.85–0.93) of
the marker's raw single positives, and had a median selectivity of 1.3 (IQR 0.9–1.9),
where selectivity is co-positive cells cleared per single-positive cell lost and a
uniform dimming of the channel scores about 1. Inside the co-positive cells REDSEA
left the marker at a median 1.05 of its raw value, so it barely treats them as
spillover recipients: most impossible co-positives come from something its one-pixel
lateral model does not see. Keeping REDSEA needed a drop of at least 20% with
selectivity of at least 2 and retention of at least 0.90, in at least two-thirds of
the evaluable donors. Of the 748 comparisons, 9 met it (569 showed no benefit, 85
lost more than 20% of the marker's single positives); no marker met it in two-thirds
of donors for any partner, five markers were donor-dependent, and a marker with no
exclusive partner cannot earn REDSEA under the rule. The 12 donors the audit could
not yet read cannot change this: the best case, INS against CD31 (2 of 13 donors),
would reach 14 of 25 donors even if all 12 benefited, against the 17 that two-thirds
requires.

What this means for each stage:

- **`expression`** reads a `Cell` marker's bare `<marker>` column and a `Nucleus`
  marker's `Nucleus__<marker>` column from `cells`, naming the marker with ingest's
  own sanitiser (`expression.raw_measurement_column`). The lookup is exact: no
  similarly named column substitutes for a missing one. It opens `data/redsea` only
  when some marker's `measurement_source` is `redsea`, so it runs on donors REDSEA
  never reached. The availability audit's `source` reads `qupath:<column>`.
- **`redsea`** keeps its place in the stage order (`ingest, geometry, redsea,
  expression`) but is a logged no-op while no marker reads it: `[no-op] redsea`, or
  `[donor-noop] redsea/<id>` under `--pipelined`, with nothing computed or written,
  and `pipeline status` prints `NO-OP`. REDSEA partitions written earlier are left
  where they are and not read.
- **`controls`, `type` and `states`** read `expression`, not REDSEA. **`calibrate`**
  (RESTORE) still reads REDSEA compartments for its NNMF features; it is blocked by
  `NO_VERIFIED_TYPING` and cannot run under this registry.
- **The vocabulary is unchanged.** `validate_registry` still accepts
  `redsea_corrected` and `redsea_passthrough`, so a marker that later earns REDSEA
  is one registry edit away; the stage then runs again, and `expression` has to be
  rerun after it. `uncompensated`, which until 2026-09-29 was allowed only for
  `Nucleus` markers, is allowed for `Cell` and `Nucleus`, the two compartments whose
  raw mean ingest keeps; a `Cytoplasm` or `Membrane` marker has no raw source.
- **Expression written before the change** holds REDSEA values for 64 markers.
  A present partition is skipped and `write_donor_partition` refuses to overwrite,
  so move `data/expression` and its availability audits
  (`audit/expression_availability.parquet`, `audit/donors/expression_availability/`)
  aside before rerunning `expression`.

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
or absence, including its `measurement_source`. No QuPath norm value is carried
into `cells` or `expression` today, because no marker's source is `qupath_norm`.

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

**The content-addressed run framework was removed on 2026-09-16.** Given a
base `data_dir = data`, the current stable layout is:

```text
data/
├── qupath_manifest.json
├── cells/donor_id=<id>/*.parquet
├── ingest_rejects.parquet
├── panel_availability.json
├── geometry_qc/donor_id=<id>/*.parquet
├── redsea/donor_id=<id>/*.parquet               # written only while a marker reads REDSEA
├── expression/donor_id=<id>/*.parquet
├── audit/
│   ├── donors/<audit_name>/donor_id=<id>/data.parquet
│   └── expression_availability.parquet
├── scratch/redsea/
└── runs/
    └── <run_id>/
        └── ...                                    # historical, read-only; see below
```

A donor's partition under a stage root (e.g. `expression/donor_id=<id>/`) is
the only signal that stage is done for that donor, but only the `--pipelined`
donor queue checks it per donor: with `--pipelined`, a present, non-empty
partition is skipped and a missing one is picked up on the next run, and
`--force` reruns it regardless of whether the source data, registry, or code
changed since. Without `--pipelined`, `run_stage` checks the whole stage
output root, not individual donor partitions — any content there at all
skips the entire stage, so a stage missing one donor's partition is not
topped up by rerunning; `--force` there reruns the whole stage. There is no
`run.json`, no per-donor receipt, and no manifest sealing a stage as complete.

Runs made before 2026-09-16 remain, untouched, under `data/runs/<run_id>/` in
their original content-addressed layout (`run.json`, `manifests/`, and
per-stage output nested under the run). That tree is historical: nothing
writes to it any more, and it is not equivalent to — or interchangeable with —
the flat roots above. The full 20-donor cohort output still lives only there;
donors `6533` and `0112` were added on 2026-09-16 and exist only in the flat tree, as
`donor_id=6533` and `donor_id=0112`.
Discovery reads `data/expression` first and falls back to the historical
run tree with a warning; do not hand-merge files between the two.

### What is not in this tree

Development typing artifacts are **not** under `data/` here. The subpanel classification and its
downstream app artifacts live a level above `core/`, under the repo-root
`data/development/typing/` (`cells/`, `subpanel/`, `immune/`, `pln/`, `islet_identity/`, `app/`).
They are mutable review artifacts, not stage outputs: no contract in this file describes them, and
nothing in the workflow above reads them. Their contract is
[`scripts/typing/RECONSTRUCTION.md`](scripts/typing/RECONSTRUCTION.md).

## Manifest contents (historical runs only)

Runs made before 2026-09-16, under `data/runs/<run_id>/`, still carry the
manifests described in this section; nothing is generated in this shape any
more for current work. Each stage manifest in that historical tree recorded:

- stage and method version;
- cohort, registry, and upstream configuration identities;
- relevant source-code hash and repository metadata;
- exact upstream artifact identifiers;
- exact expected donor set;
- output Parquet schema, row counts, file fingerprints, and UUID-universe
  hashes;
- exact fingerprints for stage-owned availability and audit sidecars.

With donor-pipelined execution, each post-ingest partition first received a
donor-stage receipt under `manifests/donors/`. That receipt carried the same
config/code/input/output contracts at donor scope. Once every donor receipt
for a stage was current, the pipeline concatenated the donor audits and sealed
the cohort stage manifest.

None of that validation runs today. A donor partition under the current flat
roots is not immutable and is not receipted: `--force` overwrites it in place,
and nothing checks whether inputs, configuration, code, schema, or donor set
changed since it was written. Correcting a stale output is the caller's
responsibility — there is no `run_id` to change to, and no `pipeline resume`
or `pipeline recover` command any more (`resume.py` and `recovery.py` were
deleted on 2026-09-16; they existed only to re-enter runs the removed
framework had locked). Rerun `pipeline run`, with `--force` where output
already exists.

## Storage, privacy, and retention

Multiplexed images and per-cell artifacts can be large. Keep `data/` outside
version control, retain the raw source files and generated cohort manifest
together, and archive the registries/rules used by published analyses. For
human tissue, follow the governing IRB, de-identification, access-control, and
retention requirements.
