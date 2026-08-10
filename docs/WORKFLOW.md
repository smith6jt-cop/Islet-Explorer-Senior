# Current QuPath-to-cell-type workflow

This document is the operational and scientific contract for the production
workflow. It explains how raw QuPath compartments and REDSEA measurements
become donor-local RESTORE marker evidence and then a broad identity.

Start with [the project README](../README.md) for commands and
[the data guide](../DATA_README.md) for schemas.

## Intended outcome

The immediate product is a comparable, auditable RESTORE record for every
donor, marker, and cell:

```text
authoritative corrected intensity
  + donor-local target/reference pair and divisor
  + RESTORE-normalized intensity
  + NNMF/control/stability audit
  + explicit positive/negative/unavailable state
```

Those records support the larger objective of assigning all cells to broad
cell types. The workflow preserves unavailability when the measurement, panel,
or control population cannot support a gate; it does not manufacture
completeness by converting missing evidence to negative or `Other`.

## Authorities and boundaries

Four versioned inputs always govern a run:

1. The QuPath cohort manifest owns raw-data identity, donor membership, image
   mapping, geometry, panel, channel, pixel calibration, and segmentation
   version.
2. [`marker_registry.json`](../phenocycler/marker_registry.json) owns marker
   meaning, compartment selection, spillover policy, RESTORE activation,
   incompatible reference, uses, and testing families.
3. [`broad_gates.json`](../phenocycler/broad_gates.json)
   owns the production broad gate membership. It carries no precedence:
   gates are evaluated simultaneously.
4. The configuration owns fixed geometry, REDSEA, RESTORE, and compute
   parameters.

Probabilistic production typing adds one inseparable authorization pair: an
immutable threshold-frozen model bundle and its separately signed promotion
manifest, plus the separate `typing_rules.json`. Both paths must be configured together; the manifest must authorize
that exact bundle, donor split, registries, calibration source, validation and
locked-test reports, replay code, and numerical runtime. Leaving both paths
blank selects the simultaneous broad typing default. See
[`PROBABILISTIC_RESCUE.md`](PROBABILISTIC_RESCUE.md) for the fitting and release
contract.

The command-line pipeline and its stage APIs are the production writers.
Notebooks may read and visualize immutable artifacts, help review audits, and
optionally invoke `run_stage` or `export_qupath`. Those calls follow the same
manifest contracts as the CLI. Notebooks must not:

- infer donors or panels from filenames;
- write files directly into `data/runs/<run_id>/`;
- create their own marker pairs, thresholds, or lineage ordering;
- replace unavailable values with zero;
- select controls from phenotype annotations or target-positive frequency;
- overwrite a stage manifest or reuse an output from another run.

## Stage graph

```text
                    ┌───────── immutable biological registries ─────────┐
                    │                                                   │
QuPath manifest ────┴─> ingest ─> geometry ─> REDSEA                    │
                          │          │           │                      │
                          └──────────┴───────────┴──> RESTORE ──> type ─┘
                                     │                (calibrate)   │
                                     └──> expression ──> controls   │
                                              │        (diagnostics)│
                                              └──> states           │
                                                                    │
type ──────────────────────────────────────────────────> QuPath export
```

The fixed execution order is unchanged:

```text
ingest -> geometry -> redsea -> expression -> controls
       -> calibrate -> type -> states
```

but note what each stage actually reads. **RESTORE reads `cells/`, `redsea/` and
`geometry_qc/` directly** — it does not consume `expression/` or
`reference_controls/`, even though `calibrate` still declares both as stage
dependencies. `reference_controls/` currently has no reader anywhere in the
package; it is retained as a diagnostic. `expression/` is a real input to
`states`, not to identity.

State annotation depends on authoritative expression, not on identity.
QuPath export depends on typing and has a separate manifest.

## Stage 0: lock the QuPath cohort

The small cohort specification is the only place raw source paths are entered.
The manifest builder resolves those paths and captures:

- one donor and exact QuPath `Image` value per image record;
- full measurement-CSV fingerprint;
- qptiff fingerprint and exact marker/channel mapping;
- cell and nucleus GeoJSON fingerprints, feature counts, and UUID universes;
- X/Y pixel calibration;
- panel ID and segmentation version;
- exact expected donor set and cohort content ID.

One image record per donor is the current contract. Measurement CSVs may be
shared across image records when QuPath exported a cohort batch.

The resolver also enforces the central donor exclusions in
[`cohort.py`](../phenocycler/cohort.py). The current exclusions are donor 6457
(failed registration associated with poor DAPI staining) and donor 6579
(maintainer-designated pancreatitis outlier). A manifest containing either
donor fails; exclusions are not silently applied after ingest.

The current REDSEA raster path requires isotropic pixels. A declared X/Y
calibration mismatch, or disagreement between the manifest and qptiff
calibration, fails before correction.

Create it with:

```bash
python -m phenocycler.pipeline manifest template --out cohort-spec.json
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json
```

The manifest performs mapping; the pipeline never extracts a donor ID from an
image or file name. A CSV row whose `Image` value is not declared is fatal.

## Stage 1: ingest the measurement union

Ingest reads every unique measurement CSV named by the cohort manifest using
DuckDB. It:

- unions columns by name across exports;
- retains batch- or panel-specific markers as null where unacquired;
- maps the exact QuPath `Image` value to the manifest donor;
- preserves the UUID as `object_id`;
- retains whole-cell means and required nucleus means;
- retains context and morphology;
- applies only a small structural area floor for degenerate fragments;
- partitions the result by donor.

Panel union matters. Taking the intersection would silently discard useful
markers from donors that acquired them; filling an absent marker with zero
would falsely turn “not measured” into “measured negative.”

The CSV parser stores rejects in `ingest_rejects.parquet`. Any rejected row
stops ingest after the audit is written. `panel_availability.json` records the
observed union, intersection, and per-panel acquisitions.

## Stage 2: geometry-only QC

Geometry QC is deliberately before REDSEA and does not read marker intensity.
For each UUID it compares QuPath measurements with the declared cell/nucleus
polygons and records separate flags:

- `analysis_eligible`: may receive accepted identity and state calls;
- `estimation_eligible`: may help estimate donor-marker background;
- `spillover_context_eligible`: may participate as physical adjacency context.

The checks include minimum and donor-relative maximum area, solidity,
raster-to-QuPath area agreement, nucleus/cell area ratio, geometry presence,
and spatial duplicate detection. Current defaults include a 20 µm² cell area
minimum, four times the donor median maximum, 0.70 solidity minimum, 0.5–2.0
raster-area ratio, and duplicate checks within 1 µm with 15% area tolerance.
The cohort configuration explicitly uses `duplicate_policy =
deterministic_keep` for the known duplicated islet detections: a QC-viable
nucleus-bearing copy is preferred, then the copy retaining the most raster
pixels, then the lexicographically first UUID. Every non-survivor is excluded
from analysis, estimation, and REDSEA physical context. All identity-affecting
values are captured in the run ID.

The three flags must not be collapsed. A real cell with essentially no
cytoplasm may remain analysis eligible for nuclear measurements while being
unfit to estimate a cytoplasmic background. Conversely, a polygon excluded
from cell-level analysis may still be required to describe a real shared
boundary during spillover correction.

A missing nucleus fails closed at the cell level: the object receives
`missing_nucleus_geometry` and is ineligible for analysis and model
estimation. Its cell polygon is also removed from physical spillover context
before REDSEA. Missing canonical cell geometry and prohibited duplicates
remain stage-fatal. Geometry QC does not delete the canonical cell universe.

## Stage 3: REDSEA

REDSEA corrects lateral signal transferred across adjacent cell boundaries in
the qptiff. It consumes the declared image/channel mapping, UUID-matched
geometry, and pixel calibration. The geometry stage is a required upstream
contract. Before boundary and contact construction, canonical polygons whose
`spillover_context_eligible` flag is false are removed from the cell and
nucleus masks. GeoJSON-only fragments outside the canonical ingest universe
may remain as physical boundary context. The current production method emits
compartment-resolved measurements and preserves exact object alignment.

REDSEA rechecks qptiff channel names and contiguous indices against the
manifest and uses the manifest's canonical marker names. A changed or
misdeclared channel map fails rather than shifting marker identity by position.
Before output, an alignment contract also requires complete canonical UUID
coverage and agreement between rerasterized and QuPath means on available
probe channels. At least two of DAPI, INS, Pan_Cytokeratin, CD3e, Vimentin, and
CD31 must be acquired; each available probe requires Pearson correlation at
least 0.80 and a median raster/QuPath ratio between 0.25 and 4.

Nucleus geometry may be embedded as `nucleusGeometry` in each cell feature or
supplied as a separate UUID-matched GeoJSON. Separate nucleus labels are
relabelled into the cell-mask label space and intersected with their own cell
before correction. An embedded feature without a nucleus is retained only for
audit; its polygon is removed before correction and it cannot receive an
accepted cell identity or state call.

The marker registry decides whether an authoritative marker later uses:

- a REDSEA-corrected compartment;
- a REDSEA passthrough for a declared no-spillover marker; or
- an uncompensated raw QuPath compartment.

REDSEA therefore corrects a physical imaging artifact; it does not decide
whether a marker is positive and does not discover cell types. Nuclear and
other no-spillover policies are explicit rather than being silently treated as
whole-cell corrected signal.

Current method and parameter identities, including compartment mode, alpha,
gap bridging, normalization form, and any no-spillover exclusions, are stored
in the stage and run manifests.

## Stage 4: choose one authoritative intensity

The expression stage joins ingest, geometry QC, and REDSEA by exact UUID and
requires the same object universe. For every registered marker it follows the
registry's compartment and spillover policy and emits one column named for the
marker.

This is the boundary that removes competing meanings of “expression.” A
downstream analysis does not choose opportunistically among raw whole-cell,
raw nucleus, corrected cell, or other compartment columns.

For an acquired marker, a missing required source column is a contract failure.
For a marker absent from the declared panel, the selected column is null and
the availability audit records `marker_not_acquired_in_panel`.

The result also carries UUID, donor, image, coordinates, region context, and
the separate analysis/estimation eligibility flags. Existing QuPath
`Classification` or annotation context is not part of authoritative expression
and cannot be used to estimate calibration.

## Stage 5: retained reference-control diagnostics

An active marker has a predeclared biologically incompatible reference. The
reference is used to identify cells expected to be negative for the target; it
is not used to algebraically normalize the target and it does not set a target
positive fraction.

This stage retains the deterministic reference-control audit introduced during
the donor-calibration experiment. It is useful diagnostic context and remains
in the immutable graph so earlier upstream artifacts can be adopted. It does
not define RESTORE groups, divisors, or calls.

For each donor and distinct reference marker, the diagnostic selection:

1. starts with `estimation_eligible` cells;
2. requires finite positive reference intensity and finite X/Y coordinates;
3. divides the donor image extent into an 8 by 8 equal-width grid;
4. ranks cells within each occupied tile by descending reference intensity,
   with a stable UUID-derived tie break;
5. takes the brightest cell from every occupied tile before the
   second-brightest from any tile, continuing in deterministic rounds;
6. stops at 20,000 controls;
7. fails to an all-false mask when fewer than 10,000 are available.

The selection reads no target intensity. It therefore cannot change when the
target-positive population is rare, common, absent, or biologically shifted.
The mask records method, registry fingerprint, source, UUID, membership, and
rank.

The reference mask is built directly from its authoritative intensity, not
from a previously calibrated reference call. Reciprocal declarations such as
E-cadherin/CD31 therefore do not create an estimation cycle.

External control masks must use the canonical schema and exact UUID universe.
They are validated for registry provenance, ranks, counts, and methods.
Columns containing annotation, classification, phenotype, or cell-type labels
are rejected. Annotations remain audit-only context.

## Stage 6: RESTORE donor-marker evidence

The historical stage key remains `calibrate` for immutable manifest and
adoption compatibility. Its production method is
`restore-marker-evidence-v1+restore-pair-validation-v10`.

For every donor and active registry target/reference pair, RESTORE:

1. restricts model inputs to analysis-QC-eligible cells with finite,
   nonnegative raw target, raw reference, REDSEA target, REDSEA reference, and
   declared raw compartment features;
2. derives marker-specific raw input floors used to identify provisional
   target-only and reference-only arms and discard double-low cells from the
   fit;
3. fits a deterministic two-component NNMF on balanced arms for seeds 0, 1,
   and 2, using at most 50,000 cells;
4. orders the components by target/reference direction;
5. defines the reference-positive/target-negative control population to the
   right of the target-derived reference separator;
6. takes the configured `p999` of that population's REDSEA target values as the
   donor-marker divisor (`control_definition=reference_only`);
7. calls a cell positive iff its REDSEA target value is greater than the
   divisor and stores `restore_normalized = value / divisor`.

The marker registry is required to match the complete accepted RESTORE v10
pair web. In particular, it restores the reviewed `CD11b <- EpCAM` pair rather
than the later unvalidated `CD11b <- E_cadherin` replacement. Structural pairs use raw `Membrane+Cell` features;
markers absent from the historical compartment extraction use the validated
raw `Cell` fallback. The effective feature set is recorded per pair.

### Tri-state safety

- `positive`: measurable, valid pair, normalized value greater than 1;
- `negative`: measurable, valid pair, normalized value at or below 1;
- `unavailable`: panel/matrix measurement missing, geometry-ineligible cell,
  or technical donor-marker failure.

A separately documented biological absence may be a confident negative. A
technical failure never is. With the default
`on_technical_failure=unavailable`, the complete failure and reason remain in
`audit/restore_pair_models.parquet`; `halt` is the stricter alternative.
The accepted donor 6566 INS image review is versioned in
`restore_expression_reviews.csv` as `CONFIRMED_ABSENT`; it produces
`NO_TARGET_EXPRESSION_CONFIRMED` and explicit negative calls, not a fitted
divisor. The review file is a fingerprinted calibrate-stage input.

`marker_evidence/` also provides `log2_threshold_ratio =
log2(restore_normalized)` so the downstream hierarchy and optional model
interfaces retain one signed threshold-relative coordinate.

## Why this supplies recourse

| Observed problem | Versioned correction |
|---|---|
| marker false positives/negatives | inspect the RESTORE pair audit; change the registered reference or declared divisor policy, then refit all affected donors |
| bad spillover or segmentation | correct REDSEA/geometry upstream, then rerun the dependency suffix |
| marker calls correct but class wrong | change `broad_gates.json`, version it, and rerun typing |
| evidence cannot resolve the cell | retain `Unavailable`; never fabricate a negative or class |

QuPath error examples are regression checks for those changes, not training
labels. See `REFINEMENT_WORKFLOW.md`.

## Stage 7: simultaneous broad typing

Typing consumes the categorical RESTORE marker states, not raw intensities and
not reviewed labels. The default `typing_mode=restore_simultaneous` evaluates
every gate against the same evidence. **No gate has precedence over another.**

The gate membership is versioned in `broad_gates.json`:

| Gate | Members |
|---|---|
| Immune | CD3e, CD68, CD11b, CD79a, CD11c, Iba1, CD206, MPO, CD4, CD8 |
| Endothelial | CD31, CD34, Podoplanin |
| Epithelial | E_cadherin, Pan_Cytokeratin, Ker8_18, EpCAM, Keratin_5, INS, GCG |
| Neural | B3TUBB |
| Mesenchymal | Vimentin, SMA |

For each analysis-eligible cell:

| disposition | condition |
|---|---|
| `anchor` | exactly one gate is supported; that gate is the label |
| `Ambiguous` | more than one gate is supported; the competing gates are recorded and **no winner is chosen** |
| `Other` | no gate is supported and every gate marker is explicitly negative |
| `Unavailable` | geometry-QC ineligible, or unsupported with an unmeasured gate marker |

A positive member dominates an unavailable sibling in the same gate.

### Why the order was removed

The predecessor visited gates in a declared priority order and gave a
multi-gate cell to whichever gate came first. Measured on run
`3bf4fc8e391d007766d1`'s own marker evidence, over its 25,150,537
analysis-eligible cells, that one choice moved the cohort Immune fraction from
**20.50%** (Immune first, as shipped) to **6.75%** (Epithelial first) to
**3.33%** (Immune last). 29.3% of cells are positive for two or more gates and
only 16.2% of Immune-gate-positive cells are positive for no other gate. Under
simultaneous typing, **84% of the ordered residual's Immune calls and 75% of
its Endothelial calls are `Ambiguous`**, against 4% of its Epithelial calls.
A headline number that swings six-fold on a list ordering was being set by the
ordering, not by the evidence.

Gate names are iterated in sorted order, not file order, so reordering
`broad_gates.json` cannot change an assignment even by accident.

### SST

SST is **not** broad Epithelial evidence. `marker_registry.json` declares it
`uses: ["subtype"]` with no `broad_compartment`, and
`hierarchical_typing.SUBTYPE_ONLY_MARKERS` forbids it pulling a cell into the
broad Epithelial compartment. The retired ordered rules placed it in that gate
anyway, because the policy validated only that the gate marker *set* equalled
the active marker set. Gate membership is now checked against each marker's
registered `broad_compartment`, `calibration_status` and `uses`. Removing SST
moved 0.06 pp of cells from Epithelial to Other.

### Scores

`broad_best_probability`, `broad_margin`, `broad_stability` and `confidence`
are `NaN`. They are not estimated by this method, and the same rows record zero
bootstrap replicates; the predecessor stamped `1.0` into all three.
`broad_support` and `broad_support_margin` carry the RESTORE-normalized support
(each marker's intensity over its own donor-local divisor, so 1.0 is the
threshold for every marker) and are **reported only** — this stage does not
threshold them. A separately versioned rule may later use that margin to
resolve ambiguous cells on evidence.

This production default is deterministic and broad-only. The compatibility
columns `specific_type` and `cell_type` carry the accepted broad label; subtype
fields are empty because no subtype classifier is being fitted.

The probabilistic parent-constrained implementation remains available only when
an exact signed and promoted bundle is explicitly configured. It is a separate
optional method, not a fallback. The ordinary production type stage never fits
or retunes a model.

## Stage 8: keep process state orthogonal

Identity answers what a cell is; state answers what it is doing. Combining the
two creates contradictory labels and lets transient biology change lineage.

The current state stage reports Ki67 and PCNA proliferation separately. Each
available marker receives a donor-local `log10-mean-plus-k-SD` state threshold
with current defaults `k = 3` and at least 200 positive finite measurements.
Only estimation-eligible cells fit the threshold, and only analysis-eligible
cells receive positive/negative state calls; other cells are unavailable. This
method is explicitly not presented as the mutually-exclusive,
frequency-neutral identity calibration. It never changes `broad_type` or
`specific_type`.

State output lives in `cell_states`; model summaries live in
`audit/state_models.parquet`.

## Stage 9: uncertainty-preserving QuPath handoff

The export stage joins assignments to image context by UUID and writes:

```text
qupath_export/cell_assignments_<donor>.csv
```

The default handoff contains the ordered broad identity, status, decisive or
blocking gate, and compatibility aliases. Specific/probability columns remain
empty or repeat the broad identity unless the explicit probabilistic lane was
configured.

An importer or reviewer must inspect `assignment_status`; it must not turn
`Unavailable` into a confident named class. State is not folded
into this identity export and remains a separate artifact. Missing image
context, duplicate UUIDs, or an altered export file is fatal.

Use
[`import_cell_assignments.groovy`](../scripts/groovy/import_cell_assignments.groovy)
for the QuPath import. First require a successful pipeline `status` so the
export file fingerprint and typing stage are current. The script validates the
exact image, UUID uniqueness, and every CSV assignment in a first pass, then
imports either the broad or specific decision. QuPath detections outside the
canonical ingest universe are intentionally left unchanged. Unavailable
results receive an explicit uncertainty class rather than the blocking gate.

## Content-addressed execution

Before writing output, the pipeline resolves a run ID from:

- cohort content ID;
- marker-registry fingerprint;
- ordered-rules fingerprint;
- identity-affecting configuration;
- relevant source-code byte hash.

The ID is the first 20 hexadecimal characters of the resulting SHA-256. All
artifacts are written below:

```text
data/runs/<run_id>/
```

Every stage manifest records its exact upstream manifest IDs, method,
configuration, expected donors, relevant source hash, output schema, file
fingerprints, row counts, UUID universe, and fingerprints for its
availability/model-audit sidecars. Stage files are immutable.

If all eight analytical manifests are current, the pipeline writes `run.json`
and updates `data/runs/LATEST`. QuPath export has its own manifest and does not
change the analytical run identity.

This design gives three useful outcomes:

- an identical invocation validates and reports `[current]`;
- a legitimate scientific or source change resolves to another run directory;
- a partial or silently altered directory raises instead of being accepted.

Routine CLI status uses fast validation: unchanged path, size, and modification
time reuse the captured digest, while changed metadata triggers a rehash.
Archival audits can call manifest/cohort validation with `mode="content"` (or
`validation_mode="content"`) to recompute the recorded full or explicitly
sampled digest.

## Running, resuming, and checking

Run the complete workflow:

```bash
python -m phenocycler.pipeline run --config config.ini --pipelined
```

Run a dependency prefix:

```bash
python -m phenocycler.pipeline run --config config.ini \
  --pipelined --through calibrate
```

Run exact named stages whose dependencies already have current manifests:

```bash
python -m phenocycler.pipeline run --config config.ini \
  --only type states
```

### Adopt compatible completed stages from an older run

When workflow development changes the content-addressed run ID, do not repeat
every expensive stage merely because an unrelated downstream feature changed.
Create a deterministic plan against the finalized current checkout. Planning
does not modify either run; `--out` writes only this requested immutable file:

```bash
python -m phenocycler.pipeline adopt plan \
  --config config.ini \
  --from-run <completed_run_id> \
  --out adoption-plan.json
```

The planner evaluates the real DAG, not a linear prefix. Each source manifest,
dataset, required audit sidecar, donor set, schema, UUID universe, source input,
method, and stage-specific scientific configuration is validated without
mistaking historical code for current code. Producing-code hashes must either
match exactly or match a named, exact-hash transition in the checked-in
`adoption_compatibility.json`; an arbitrary runtime override is not available.
Review `action`, `reason_codes`, the reusable byte count, and the plan content ID.
The unit of reuse is a sealed stage: a historical root `run.json` is informative
but not required, so a prior run with individually complete manifests can still
be salvaged. Unmanifested donor output is never eligible.

Apply exactly that plan, then inspect `manifests/adoption/` before launching the
missing work:

```bash
python -m phenocycler.pipeline adopt apply \
  --config config.ini \
  --plan adoption-plan.json

python -m phenocycler.pipeline run --config config.ini --pipelined
```

Apply rebuilds the plan before writing, so source, policy, configuration, or
checkout drift causes a fail-closed stale-plan error. Accepted files are
hard-linked, never copied or symlinked. Each adopted stage gets a retained
source-manifest snapshot, a content-addressed receipt, and a new target manifest
under the current dependency, configuration, method, and code compatibility
contract; the receipt retains the original producing code. Old donor receipts,
assignments from an incompatible typing method, QuPath exports, `run.json`, and
`LATEST` are not transplanted. Add `--continue` to `adopt apply` to launch the
missing stages immediately. `--jobs N` and `--no-export` modify that continuation
only; after `--no-export`, run the normal `export` command separately.
Apply is idempotent and stage-incremental, not one filesystem-wide transaction:
after interruption, rerun the same immutable plan to validate existing links
and continue from the first unsealed stage.

The ordinary mode uses cached fingerprints and rehashes files whose metadata
changed. For a one-time archival read of every source byte, add
`--validation-mode content` when creating the plan. Cross-filesystem adoption
fails rather than silently allocating a second copy. Plans are immutable: after
source, configuration, policy, or code drift, write a newly named plan rather
than trying to replace the reviewed one.

For `c8412757d12de36ca056`, the reviewed plan adopts `ingest`, `geometry`,
`redsea`, `expression`, `controls`, `calibrate`, and `states`, preserving
99,139,122,797 bytes (92.33 GiB). `type` must rerun because its method changed
from `hierarchical-typing-stage-v1` to v3; the QuPath export then reruns because
it is bound to the new type-manifest ID. `states` remains reusable because it
depends on `expression`, not `type`. That source lives under this subrepository's
data root: from the parent workspace, first `cd Phenocycler_Analysis` (or pass
`--config Phenocycler_Analysis/config.ini` from the parent).

`--from-run` accepts a 20-character lowercase run ID, not an arbitrary path,
and resolves it only under the runs directory owned by the selected config.
The target records the reviewed plan at `manifests/adoption/plan.json`, retained
source manifests under `manifests/adoption/source/<source_run_id>/`, per-stage
receipts under `manifests/adoption/`, and the aggregate handoff at
`manifests/adoption.json`. `LATEST` is unchanged by adoption and advances only
after every analytical stage validates and the target `run.json` is sealed.

#### Reviewing future compatible code drift

An `unreviewed_code_change` decision is intentionally not overridable at the
command line. A maintainer must first decide whether the stage method should be
versioned and recomputed. Only a demonstrably non-scientific transition may be
added manually to `phenocycler/adoption_compatibility.json` in a reviewed code
change. Each record pins the stage, old/new method and aggregate code hashes,
source/target path contracts, exact changed-path set, old/new digest for every
changed file, canonical diff digest, accepted output schema, and the stage
configuration projection. Commits are recorded as provenance; compatibility
identity comes from the pinned code trees. Its rationale must explain why
producer values and parameter wiring are unchanged. The planner reconstructs
the historical aggregate code from the source manifest's commit and verifies
every pin; missing Git history or any mismatch leaves the stage at `recompute`.
Never add a policy record to avoid a legitimate method bump, schema migration,
or failed scientific comparison.

For the narrow case where `duplicate_policy = fatal` stopped geometry at the
first donor containing duplicated islet detections, reuse the verified ingest
and completed duplicate-free geometry prefix in the corrected
`deterministic_keep` run:

```bash
python -m phenocycler.pipeline resume \
  --config config.ini \
  --from-run <failed_run_id> \
  --from-donor <first_unfinished_donor>
```

This is not an unchecked partial-output override. It requires identical ingest
code, cohort input, ingest filter, donor set, and ingest content; requires the
only geometry configuration change to be `fatal -> deterministic_keep`;
verifies each prefix donor's exact row/UUID universe, nucleus-exclusion
semantics, and absence of duplicate groups; and records both reuse operations
in fingerprinted receipts. Ingest is hard-linked, not copied. Prefix geometry
values are preserved while its two all-null duplicate fields are re-encoded
with stable Parquet types. The default is to continue through the donor task
queue, so reusable-prefix REDSEA work can overlap geometry beginning at the
unfinished donor. Use `--through geometry` for a bounded resume.

For a stopped run that predates donor receipts but has a complete, current
geometry manifest, use the separate recovery contract:

```bash
python -m phenocycler.pipeline recover \
  --config config.ini \
  --from-run <stopped_run_id>
```

The source process must be fully stopped first. Recovery validates the source
ingest and geometry manifests, producing code, complete scientific
configuration, donor set, row counts, schemas, and UUID universes. It then
hard-links the complete geometry dataset and publishes one current geometry
receipt per donor. An unmanifested REDSEA donor is accepted only when its
single final `data_0.parquet` opens cleanly, passes a full-column scan, has the
exact current schema and donor identity, and matches that donor's ingest
row/UUID universe. Source scratch and incomplete donor directories are
ignored. The accepted REDSEA donors receive receipts and can enter expression
as soon as the queue starts; unfinished donors return to REDSEA.

Add `--prepare-only` to perform and record the handoff without launching work.
Recovery writes `geometry_recovery.json` and, when applicable,
`redsea_recovery.json`; both are fingerprinted by the donor receipts.

### Donor receipts and task resources

Ingest remains a cohort barrier because the shared measurement exports are
scanned together. Every later task is donor-local:

```text
geometry[d] -> redsea[d] -> expression[d] -> controls[d]
                                  │              └─ calibrate[d] -> type[d]
                                  └─ states[d]
```

A downstream task becomes ready only after each required upstream donor
receipt validates. A receipt records its donor, method, scientific
configuration, relevant code hash, upstream receipt or manifest identifiers,
output schema, row count, UUID-universe hash, file fingerprints, and any
donor-local audit sidecar. Writes use a temporary file and atomic rename; the
receipt is published last. Cohort stage manifests remain the archival
completion contract. They are sealed from the complete receipt set and
fingerprint every receipt. During an active queue, the calibration review
notebook may instead read a completed donor only after validating that donor's
receipt tree; it uses donor-local audit sidecars until the sealed cohort audit
is available. This read path does not authorize notebook-side stage execution
or make a partial cohort complete.

The queue has three dedicated process pools and enforces both a total task
limit and group limits:

```ini
[compute]
n_jobs = 4
geometry_workers = 2
redsea_workers = 1
downstream_workers = 4
use_gpu = false
gpu_device = 0
```

At most `n_jobs` tasks run at once. The group settings are ceilings within that
total, so this example permits two memory-heavy geometry tasks, one REDSEA
task, and any remaining slot for downstream work. When `use_gpu = true`, the
REDSEA pool is forced to one worker for the configured device even if
`redsea_workers` is larger. Separate pools keep the GPU context in the REDSEA
worker and keep geometry memory pressure bounded.

`status` reports `IN PROGRESS n/N donor receipts` until a cohort stage is
sealed. Re-running the same pipelined command validates and skips current
receipts. If all receipts exist but execution stopped before the aggregate
manifest was written, the next invocation seals that manifest without
recomputing donors.

Check every stage and the QuPath export:

```bash
python -m phenocycler.pipeline status --config config.ini
```

`status` returns nonzero if anything is missing or stale. Build the handoff
separately with:

```bash
python -m phenocycler.pipeline export --config config.ini
```

Run `status` immediately before a standalone export and proceed only if every
analytical stage is `CURRENT`; the QuPath line may be `MISSING`. Run it again
after export and require success before QuPath import. This validates the
typing stage and fingerprinted export as one handoff review.

This means an intentional `--no-export` analytical run is complete at
`run.json`, but `status` still reports the optional handoff missing and returns
nonzero until it is exported.

`--only` is intentionally literal; it does not guess or rerun prerequisites.
A missing required manifest stops the selected stage.

## Failure and recovery contract

| Failure | Meaning | Correct response |
|---|---|---|
| source fingerprint changed | raw export no longer matches the cohort manifest | verify the intended source and recreate the cohort manifest |
| cohort donor set differs | incomplete, extra, excluded, or duplicated donor record | correct the cohort specification |
| unmapped `Image` value | CSV contains an image not owned by the manifest | correct the exact `image_id` or source export |
| CSV parser rejects | one or more source rows could not be parsed | inspect `ingest_rejects.parquet` and fix/re-export the source |
| missing/duplicate geometry UUID | canonical geometry universe is invalid | re-export or correct segmentation/geometry |
| anisotropic or mismatched pixel calibration | current REDSEA raster units are not valid | correct the manifest/source calibration; do not rescale silently |
| acquired source column missing | panel says marker exists but its required compartment does not | correct channel mapping, registry, or upstream export |
| qptiff channel absent from marker registry | a channel has no declared biological/measurement policy | review and add a versioned registry entry or correct the channel map |
| REDSEA/QuPath alignment failure | UUID, channel, scale, or raster units do not agree | correct the manifest/export; do not accept the corrected values |
| insufficient reference controls | donor/reference cannot support the fixed estimator | retain indeterminate results; review panel/reference design |
| insufficient tail resolution | control count cannot resolve the fixed marker alpha | retain indeterminate results; do not loosen alpha post hoc |
| excess contamination | incompatible-reference controls contain too much target signal | audit biology/pairing and retain indeterminate results |
| marker panel absent | marker was not measured | retain unavailable state |
| required typing evidence invalid | broad absence cannot be established | retain unavailable/ambiguous; do not call `Other` |
| cohort output exists without manifest | interrupted or foreign stage-mode partial artifact | inspect it; never bless it by fabricating a manifest |
| donor output exists without donor receipt | interruption occurred after the atomic data write but before provenance publication | quarantine that exact donor partition and its donor audit, then rerun the pipelined command |
| older run has completed stages but current code resolves a new run ID | compatible work may be reusable on the true stage DAG | create and review an `adopt plan`, then apply that exact plan; incompatible nodes and descendants recompute |
| fatal duplicate policy stopped at the first duplicate-bearing donor | ingest and a duplicate-free geometry prefix may be reusable | use the verified `pipeline resume` contract with the source run and first unfinished donor |
| stopped pre-receipt run has a complete geometry manifest and some complete REDSEA donors | the old stage-barrier process cannot attach to the donor queue | stop every old worker, then use the verified `pipeline recover` contract |
| manifest or output stale | bytes, schema, UUID universe, config, or code differ | correct the cause and create the newly addressed run |
| audit/availability sidecar stale | a stage-owned diagnostic changed or disappeared | treat the stage as stale; do not regenerate the sidecar independently |

Never repair a contract failure by copying a donor partition from another run,
editing a generated manifest, reallocating alpha to available markers, or
adding a notebook-only cutoff. Cross-run reuse is permitted only through the
verified adoption, resume, or interrupted-run recovery contracts above. Each
creates a newly addressed run and fingerprints the exact provenance handoff.

For a verified transient interruption that left output before its manifest,
quarantine the exact unmanifested stage directory and its stage-owned
sidecars, then rerun the same stage. Do not remove current upstream manifests
or broadly clear the run root. Under pipelined execution, a current donor
receipt is resumable and must not be removed; only an exact orphan partition
without a receipt is quarantined. The aggregate sidecar ownership is:

| Stage | Sidecars |
|---|---|
| ingest | `ingest_rejects.parquet`, `panel_availability.json` |
| expression | `audit/expression_availability.parquet` |
| controls | `audit/reference_control_models.parquet` |
| calibrate | `audit/restore_pair_models.parquet` |
| states | `audit/state_models.parquet` |

## Registry and rule change control

A proposed marker or typing change should be made in the owning versioned file,
not scattered across notebooks or stage code.

For a marker-registry change, review:

- biological kind and broad/state use;
- authoritative compartment and spillover policy;
- active versus deferred RESTORE calling;
- incompatible reference and its feasibility in each panel;
- fixed family membership and downstream rule coverage;
- panel-absent behavior;
- accepted RESTORE pair-web invariants, including `CD11b <- EpCAM`.

For an ordered typing-rule change, review:

- gate priority;
- marker membership and whether any marker appears twice;
- required marker availability and unavailable-gate behavior;
- validation evidence and rules-version change.

Both file fingerprints enter the run ID, so a reviewed change cannot silently
reuse the old assignments.

## Review checklist

Before using a run for biological inference:

- `status` reports every required stage current.
- `run.json` and all eight analytical manifests exist.
- expected donor counts and UUID universes agree across stages.
- `panel_availability.json` matches the intended acquisitions.
- expression-availability audit has no unexplained acquired-source failures.
- geometry exclusions and estimation-only differences have been reviewed.
- reference-control counts and spatial coverage are plausible for every donor.
- RESTORE pair state, divisor, stability, controls, feature set, and call
  fractions are reviewed per donor and marker.
- unavailable fractions are reported, not discarded.
- ordered gate and blocked-gate fractions are reviewed.
- `Other` cells have all 24 gate markers explicitly negative.
- state results have not been merged into identity.
- the QuPath review uses UUID and preserves assignment status.

## Restoration and migration

RESTORE is the production marker-calling method. Historical loose
`restore_gated_redsea/` outputs are not copied into a content-addressed run;
the current stage refits RESTORE from fingerprinted raw compartments and
REDSEA, preserving full provenance. For the failed run
`c8412757d12de36ca056`, use the reviewed adoption contract to reuse valid
upstream stages and recompute only RESTORE plus typing. See
`REFINEMENT_WORKFLOW.md`.
