# Running and reviewing the portability evaluation

This guide explains how to inspect the completed 2026-11-02 evaluation, repeat its
supported computational checks, and prepare a new candidate evaluation. It covers
protein typing, RNA-only typing, same-slide comparison, serial-section MALDI and
project-only parameter recording. The [evaluation report](evaluations/dynamic_portability_20261102/README.md)
is the authority for the dated findings; this guide is an operating procedure.
The [handoff completion record](evaluations/dynamic_portability_20261102/HANDOFF_COMPLETION.md)
distinguishes delivery of the evaluation from scientific routes that remain blocked.

The evaluation changes no scientific method, approves no registry and activates
no release. Its implementation addition is the KINTSUGI observation adapter.
Scientific corrections are specified in the [remediation document](evaluations/dynamic_portability_20261102/remediation.md).
Files under `outputs/maldi_completion_20261102/proposals/` are inactive proposals,
not installed registration or object-integration capabilities.

## 1. Start with the evidence and the correct environment

From the repository root on the workstation:

```bash
export PORTABILITY_REPO="$PWD"
export PYTHONPATH="$PORTABILITY_REPO/core/multiplex5000"
export MPLX5K_SITE_FILE="$PORTABILITY_REPO/core/multiplex5000/sites/sn4622120742.json"
.venv/bin/python -m multiplex5000.cli status --project io60
.venv/bin/python -m multiplex5000.cli status --project io60_maldi_serial
```

Here `.venv/bin/python` is the existing scientific environment, backed by
`/data/envs/multiplex5000`. The site file names workstation `/data` and `/home`
sources. It is not a HiPerGator site file. Python version, dependency versions,
checkout commit and dirty diff belong in each new evaluation's provenance.

Use these entry points to inspect the delivered result:

| Question | Evidence |
|---|---|
| What was evaluated, and what completed? | [Report](evaluations/dynamic_portability_20261102/README.md), [acquisition ledger](evaluations/dynamic_portability_20261102/acquisition_dispositions.json), [cohort accounting](evaluations/dynamic_portability_20261102/cohort_accounting.json) |
| Which capabilities are available? | [Readiness matrix](evaluations/dynamic_portability_20261102/readiness.json), with separate code, input, specification, routing, computational and biological-review axes |
| Which settings depend on assay or tissue? | [Dependency traces](evaluations/dynamic_portability_20261102/dependencies-v2.csv), [parameter catalog](evaluations/dynamic_portability_20261102/parameter_catalog.json), [trace summary](evaluations/dynamic_portability_20261102/parameter_trace_summary.json) |
| Which physical pairs are established? | [Pair ledger](evaluations/dynamic_portability_20261102/pair_dispositions.json), frozen [inventories](evaluations/dynamic_portability_20261102/inventories/) |
| Where are diagnostics and their hashes? | [Artifact manifest](evaluations/dynamic_portability_20261102/artifact_manifest.json), [inspection ledger](evaluations/dynamic_portability_20261102/inspection_ledger.json) |
| What settings were actually recorded? | [Parameter history](evaluations/dynamic_portability_20261102/parameter_history.json), [recorder setup and contract](evaluations/dynamic_portability_20261102/recording.md) |
| What must be fixed before another route can run? | [Remediation specifications](evaluations/dynamic_portability_20261102/remediation.md), [fault coverage](evaluations/dynamic_portability_20261102/fault_coverage.csv) |

Tracked reports are portable. Large diagnostics, source files and SQLite history
are external or gitignored; the artifact manifest records their actual paths and
hashes. A missing local artifact is unavailable evidence, not a successful check.

The tracked [compressed parameter traces](evaluations/dynamic_portability_20261102/parameter_traces.json.gz)
and [literal candidates](evaluations/dynamic_portability_20261102/parameter_literal_candidates.json.gz)
preserve detailed source/consumer accounting without expanding large JSON in the
working tree. For example:

```bash
gzip -dc docs/evaluations/dynamic_portability_20261102/parameter_traces.json.gz | less
```

Static coverage accounts for catalog declarations, registry rows and remaining
literal candidates. It does not prove that every setting governed execution or
that every fitted runtime value was recorded; the trace summary explicitly
retains those limitations and the separately pinned historical producers.

Run the read-only consistency check after the dated delivery bundle is finalized:

```bash
.venv/bin/python scripts/evaluation/verify_portability_delivery.py
```

This checks the frozen roster, identities, observation revisions, evidence hashes,
outbox and artifact manifest. It does not establish numerical parity, biological
validity or registration acceptance. A failing hash means the referenced evidence
has changed or is inaccessible; investigate that specific reference.

### Keep producer versions separate

The full hashes and provenance are in [pins.json](evaluations/dynamic_portability_20261102/pins.json).
The following are distinct scientific/evaluation baselines:

| Component | Pinned commit |
|---|---|
| Local main baseline | `9b75b75ed812b68718704708b89309401edd4732` |
| Serial MALDI branch | `00ec0322d8b8da6529283a65db9105e62e0bdf85` |
| Per-image typing port | `df3fe1d6ecfe534bc1eaf9e1cf94b8ef6d3a7018` |
| Originating CODEX method | `578767eeaf0f95bb1df7272238f361888a0ae75d` |
| CODEX v20 producer | `9b1da13016bf9c9a1266c3ce212f087399a5e264` |
| Audited HPG RNA producer | `7fe0d32119687fe7f68792598e429678e7023284` |
| Audited later HPG spleen producer | `a9549cce6674d7a7536a325ae650c7210b1a9bf7` |
| KINTSUGI engine checkout | `716d1ebd87db916ab04910c6374e04d2d9ea508e` |

Auditing a later producer's saved output does not make its code part of the local
baseline. Original-method parity must be established independently of v20
compatibility. A historical pass does not validate files that now have different
hashes. Record the producer receipt, registry/panel hashes, source manifests and
environment alongside the code commit.

## 2. Freeze scope, preserve identities and choose a new output root

The dated ledger contains 201 acquisitions: 129 in the IO60 baseline, 32 MALDI,
20 CODEX, 12 spleen protein and eight spleen RNA. The IO60 baseline includes 20
post-Xenium IF acquisitions explicitly excluded from the requested routes. These
are acquisitions, not counts of donors, unique sections or biological replicates.

Use inventory identifiers exactly. Preserve native QuPath UUIDs/Xenium cell IDs,
qualified cell IDs, segmentation versions and MALDI spot IDs. Donor equality,
similar filenames and image proximity do not establish a physical relationship.
Frozen inventory CSV bytes are hash-bound; preserve their line endings.

Every declared acquisition receives one of these dispositions, with evidence:

| Disposition | Meaning |
|---|---|
| `completed` | The stated evaluation step completed; its scope still distinguishes computation, geometry and biological review. |
| `excluded` | An existing explicit rule excludes it from this route. Record the rule and its modality scope. |
| `unavailable` | Evidence cannot support a result, for example a damaged or insufficient-signal region. Preserve the acquisition and available counts. |
| `blocked` | A required source binding, method specification, implementation or accepted geometry is absent. Record what is missing and the concrete remediation. It does not mean the acquisition was discarded or that an unrun method failed numerically. |

The donor 6457 exclusion is PhenoCycler-specific. Nucleus admission still applies
to cells. Damaged RNA regions remain unavailable. Two selected RNA-only donors
do not establish population disease effects. Expected frequencies, previous
labels, disease ordering and cross-assay agreement are not tuning objectives.

For a fresh check, choose a new directory. The example below is intentionally
separate from the delivered evidence; choose another name if it already exists.

```bash
export PORTABILITY_CHECK="$PORTABILITY_REPO/outputs/portability-user-guide-recheck"
mkdir -p "$PORTABILITY_CHECK"
```

Read-only validators require a new `--out` directory outside the scientific
release. Do not overwrite evidence referenced by the parameter history. Several
dated assembly scripts have fixed output paths; inspect their constants before
reuse. `manifest_portability_artifacts.py` writes the delivery manifest and is a
finalization step, not a read-only check.

To prepare an isolated candidate, copy its project YAML into the new directory,
give it a unique `release`, and keep `id` equal to the YAML filename stem. Set an
absolute `releases_root`, or supply `--releases-root`. Review `sources`, exact
acquisition selection and `automatic_gates`. A new root alone does not import
upstream receipts or authorize reuse of another release's measurements.

For example, this prepares the existing RNA-only route with no automatic gates:

```bash
.venv/bin/python - <<'PY'
import os
from pathlib import Path
import yaml
base = Path(os.environ["PORTABILITY_REPO"])
out = Path(os.environ["PORTABILITY_CHECK"])
source = base / "core/multiplex5000/projects/io60_xenium_only_6512_6618.yaml"
project = yaml.safe_load(source.read_text())
project.update(id="rna_only_recheck", release="rna-only-recheck-1",
               releases_root=str(out / "releases"), automatic_gates=[])
with (out / "rna_only_recheck.yaml").open("x") as handle:
    yaml.safe_dump(project, handle, sort_keys=False)
PY
.venv/bin/python -m multiplex5000.cli graph-spec \
  --project "$PORTABILITY_CHECK/rna_only_recheck.yaml"
.venv/bin/python -m multiplex5000.cli stage run 00_inventory \
  --project "$PORTABILITY_CHECK/rna_only_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 00_inventory \
  --project "$PORTABILITY_CHECK/rna_only_recheck.yaml"
```

Review the emitted gate page and source hashes. Subsequent stages require the
previous gate's acceptance according to the project contract. Acceptance is an
explicit recorded review action; `gate show` and a passing integrity check do not
perform it. The actual syntax, after review, is `gate accept STAGE --project
PROJECT --by REVIEWER --note REVIEW_NOTE`. Do not manufacture a reviewer or mark a
missing inspection complete. The runner's `--force` moves an old stage aside;
fresh evaluation roots normally make it unnecessary.

## 3. Protein typing: preserve the measurement contract

Compare the source and consumer at each boundary: native cell identity and
nucleus admission; raw/REDSEA/normalized measurement; whole-cell positivity
versus membrane differentiation; coverage/despeckle policy; calibration and units;
reference populations; fitted levels; method selection; and review criteria.
An available marker column is not proof that the marker was acquired or measured.
A registry file's existence is not proof that every biological decision reads it.

### CODEX parity and v20 compatibility

The existing full-cohort parity runner uses the original port's own three tests
for every image. With the existing snapshot and baseline on this workstation:

```bash
.venv/bin/python scripts/evaluation/codex_port_parity.py \
  --snapshot "$PORTABILITY_REPO/outputs/dynamic_portability_20261102/typing_port" \
  --run /data/panc_CODEX/Islet-Explorer/data/development/uniform-20donor-20261029-v19-qc-work \
  --out "$PORTABILITY_CHECK/codex-parity"
```

Inspect each image's log/XML and `summary.json`. The initial dated run stopped in the
original provenance fixture for all 20 images before numerical comparisons.
Changing timestamps or weakening the fixture to obtain a pass is not parity.

A subsequent read-only recovery found the original execution journals. The
[journal proof](../outputs/handoff_completion_20261102/codex/journal_proof.json)
checks all 95 output hashes and the source bytes against the original CODEX pin;
the original pytest fixture remains unchanged and still rejects its mtime rule.
The final [CODEX completion summary](evaluations/dynamic_portability_20261102/codex_completion_summary.json)
records **60 of 60 exact direct source-to-port comparisons passing** across all
20 images: levels, canary and identity. The separate
[unchanged artifact assertions](../outputs/handoff_completion_20261102/codex/numerical_parity.json)
also pass all 20 identity comparisons, covering 3,199,606 eligible cells and
reconciling 3,525,977 measured objects, while retaining 40 exact levels/canary
CSV-comparison failures. Thus computation parity is established for these tested
source/port paths; the original artifact test suite still does not pass.

The [direct diagnostics](../outputs/handoff_completion_20261102/codex/source_to_port_diagnostics.json)
and final precision follow-up separate computation equality from CSV
parsing/serialization and runtime effects. Under the recovered driver's single
BLAS thread, original producer interpreter and round-trip CSV parsing, all 20
images have zero differences in the positivity, differentiation and canary
tables. This is diagnostic evidence: no original fixture, production method or
tolerance was changed. All 129 checked input/baseline files remained unchanged.
The proof also retains upstream binding limitations, including endocrine inputs
originally bound by filesystem metadata rather than producer SHA256.
The [recovery harness](../scripts/evaluation/recover_codex_parity.py) and
[diagnostic harness](../scripts/evaluation/diagnose_codex_numeric_parity.py) use
dated output roots; preserve those evidence files when preparing another run.

The v20 audit separately found changed identity/levels/structure bindings. Its
dated script, [audit_codex_v20.py](../scripts/evaluation/audit_codex_v20.py), uses
fixed paths; preserve the original result and choose new output paths before a
fresh source audit. Do not treat v20 as a replacement baseline for port parity.

### IO60 pancreas and embedded pLN

The [full candidate project](../core/multiplex5000/projects/io60_mixed_protein_full_candidate.yaml)
selects all 25 donors but needs exact upstream source bindings. The actual
corrected evaluated configuration is
[`protein_evaluation_v2.yaml`](../outputs/dynamic_portability_20261102/protein_evaluation_v2.yaml).
It binds the existing upstream release, geometry, registry log, shifts, exports
and recovered acquisition-run table. Copy this configuration to a new project
before another run; preserve its source hashes and change the project ID/release.

```bash
.venv/bin/python -m multiplex5000.cli status \
  --project outputs/dynamic_portability_20261102/protein_evaluation_v2.yaml
```

The existing attempt built all 25 primary tables and reconciled their eligible
non-pLN native IDs, then failed at refinement. Ker8_18 is unavailable in 0112,
6347, 6533 and 6580; refinement treats an all-NaN column as present and reaches an
empty fitting population. No donor completed refinement and there is no completed
mixed Stage 06 receipt. Repeating `stage run 06_typing` with unchanged code is a
failure reproduction, not a remedy. See PROT-1 in the remediation document.

Check the primary expected-ID universe, source nucleus fields, admission masks
and aligned/primary coordinate contracts separately. The primary adapter uses
`0.5077663810243286 µm/pixel`; aligned data use `0.5078`. The evaluation explains
the resulting scale difference without inventing a tolerance. Reusing one frame
as the other would require an explicit calibrated conversion.

The default `io60.yaml` runs inventory, protein measurement and normalization;
it does not implicitly run typing. Whole-cell and membrane signals, pancreas
hormone references, tissue context and refinement populations retain their
producer-specific meanings. See the [mixed pilot description](IO60_MIXED_TYPING_PILOT.md)
and the evaluated dependency traces before proposing a transfer.

### Dedicated pLN

Audit exact source exports first:

```bash
.venv/bin/python -m multiplex5000.cli sources pln
.venv/bin/python -m multiplex5000.cli sources pln \
  --image 6380pLN_Scan1.er.qptiff --image 6436pLN_Scan1.er.qptiff
.venv/bin/python -m multiplex5000.cli status --project pln_protein_pilot
.venv/bin/python scripts/evaluation/validate_typing_release.py \
  --root /data/io60-current/releases/pln-protein-pilot-candidate-1 \
  --out "$PORTABILITY_CHECK/pln" --kind pln
```

The existing [pilot](../core/multiplex5000/projects/pln_protein_pilot.yaml) has two
sections with bound Stage 03/05/06 artifacts. Its roster is configurable, but the
other 17 of 19 acquisitions do not thereby acquire those upstream artifacts or
review. The current and historical normalization manifests differ; use the
recovered producer binding for historical reconciliation. Pilot candidate calls
are not a completed pLN registry for the separate per-image typing port.

### Spleen PC67 and PC33

The local [spleen primary project](../core/multiplex5000/projects/spleen_primary.yaml)
runs Stage 03/05. The audited HPG producer `a9549cce` separately supplies PC67
typing. Its consumer reads raw REDSEA values and section-fitted levels; the
Stage 05 prerequisite does not imply that it consumes aligned values.

All eight PC67 sections reconcile computationally at that later producer. Ten
large image/geometry input hashes remain unavailable. P1–P9 in the original
port remain proposals; later results do not approve them retroactively. All four
PC33 sources were audited, but their lineage roots/registry are missing. Do not
apply PC67 registry settings to PC33 or count source-ID validation as typing.

The HPG command in section 4 can be reused with
`--root /blue/maigan/smith6jt/multiplex5000-spleen-primary/releases/spleen-primary-candidate-1`
and `--kind spleen`. The validator includes PC33 source dispositions; inspect
each row's status instead of interpreting the total acquisition count as 12
completed typing results.

## 4. Xenium-only typing and full-cohort RNA checks

The [RNA-only project](../core/multiplex5000/projects/io60_xenium_only_6512_6618.yaml)
selects exactly `xenium:0041323:Region_1` and `xenium:0041326:Region_1`. Its stages
are 00, 04 and 06; it requires neither protein measurements nor protein
normalization nor same-slide pairing. After the new inventory's actual review,
the commands are:

```bash
.venv/bin/python -m multiplex5000.cli stage run 04_rna_ingest \
  --project "$PORTABILITY_CHECK/rna_only_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 04_rna_ingest \
  --project "$PORTABILITY_CHECK/rna_only_recheck.yaml"
# After actual Stage 04 review and acceptance:
.venv/bin/python -m multiplex5000.cli stage run 06_typing \
  --project "$PORTABILITY_CHECK/rna_only_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 06_typing \
  --project "$PORTABILITY_CHECK/rna_only_recheck.yaml"
```

Local source availability is incomplete for some acquisitions. The exact HPG
[source locator](../core/multiplex5000/sites/hipergator/xenium_sources_20261026.csv)
and [setup notes](../core/multiplex5000/sites/hipergator/README.md) provide the
remote bundle identities. Do not search for a similar directory and substitute
it. The [full-cohort project](../core/multiplex5000/projects/io60_xenium_full_same_release.yaml)
declares 65 acquisitions. The older two-region `io60_xenium_pilot` is not a
full-cohort substitute.

The dated HPG release is
`/blue/maigan/smith6jt/mplx5k-io60-rna-20261027/releases/io60-candidate-1`.
The existing `venv-current/bin/python` and `site-current.json` under that HPG
project are recorded runtime paths, not a promise of the original producer pin.
For a fresh read-only audit, copy the audit script from this checkout and run on
a compute node. This does not rerun typing or alter the release:

```bash
# Workstation: prepare a separate remote audit directory and copy the script.
ssh hpg 'mkdir -p /blue/maigan/smith6jt/portability-user-guide-recheck'
scp scripts/evaluation/validate_typing_release.py \
  hpg:/blue/maigan/smith6jt/portability-user-guide-recheck/
```

In an allocated HPG compute shell, use the existing scientific runtime:

```bash
export PORTABILITY_HPG=/blue/maigan/smith6jt/mplx5k-io60-rna-20261027
export PORTABILITY_ENV=/blue/maigan/smith6jt/miniforge3/envs/xenium_analysis
export LD_LIBRARY_PATH="$PORTABILITY_ENV/lib:${LD_LIBRARY_PATH:-}"
export MPLX5K_SITE_FILE="$PORTABILITY_HPG/site-current.json"
export TMPDIR="$PORTABILITY_HPG/tmp"
"$PORTABILITY_HPG/venv-current/bin/python" \
  /blue/maigan/smith6jt/portability-user-guide-recheck/validate_typing_release.py \
  --root "$PORTABILITY_HPG/releases/io60-candidate-1" \
  --out /blue/maigan/smith6jt/portability-user-guide-recheck/rna --kind rna
```

For a new scientific run, install the intended pinned checkout and use a
host-specific site file; do not assume this audit interpreter reproduces the
recorded producer. HPG's default login Python is insufficient for the documented
new-environment setup; the site notes describe loading Python 3.12 and keeping
large scratch/release files under writable `/blue` storage.

Read `summary.json`, `receipt_checks.csv`, per-acquisition reports and overlays.
The delivered 65-region audit accounts for 60 candidates, three damaged and two
low-signal unavailable regions. Missing genes remain distinct from measured
zeros. Inspect actual panel genes, desired versus present identity-core genes,
controls, nucleus admission, fitted clusters and unresolved/broad/uncalled states.
For example, the recorded panel for `0041323/Region_1` lacks INS/GCG/SST/PPY;
surrogate endocrine evidence cannot establish alpha/beta/delta subtypes.

Keep fixed detection settings separate from fitted embedding values. RNA
detection uses the recorded target of 10,000 counts per cell and its documented
log predicate; pLN clustering separately uses Scanpy's section-median depth when
no target is supplied. The reconstructed pLN targets are evidence derived from
bound inputs, not originally emitted fit records. The fixed 20-gene region screen
and tissue-worker thresholds retain their recorded scope. Historical donor lists,
partitions and labels must not silently select a new acquisition's worker or
identity resolution.

For an independent state reconciliation, copy and run
[reconcile_rna_states.py](../scripts/evaluation/reconcile_rna_states.py) on that host
with `--root` set to the RNA release and `--out` to a new JSON path. Specific,
broad, unresolved, uncalled and unavailable totals must sum to the segmented-cell
universe; nucleus eligibility is a separate partition. Preserve both.

## 5. Same-slide protein/RNA comparison

First establish independent assay calls. The explicit `io60_joint_typing` route
produces assay outputs; it does not establish same-cell links. The pair ledger
contains four explicit pancreas declarations and eight verified spleen pairs;
pLN physical-section relationships remain unresolved. Do not pair by donor.

The [spleen source manifest](../core/multiplex5000/sites/hipergator/spleen_same_slide_20261026.yaml)
binds eight PC67 sections, their export headers, source images, section-table crop
offsets and Xenium regions. The local `spleen_pilot` selects only one RNA region,
`xenium:0110618:090`. Its bounded Stage 07/08 route is not full eight-pair support.

For that existing pilot, a supplied registration belongs under the site's
`sources.same_slide_registrations/<protein_entry>/`. Its `transform.json` must
declare the exact protein entry and Xenium acquisition, `moving_frame:
region_local_um`, `fixed_frame: crop_local_um`, and `reviewed: true`; include
`transform_displacement.npz` when declared. That metadata declares prior review,
not permission to skip the current Stage 07 gate. See the [source setup notes](../core/multiplex5000/sites/hipergator/README.md).

Once a real transform and upstream accepted receipts exist, the actual commands
are `stage run 07_registration --project spleen_pilot`, `gate show
07_registration --project spleen_pilot`, and, after review, `stage run
08_integration --project spleen_pilot`, through the CLI prefix used above. Use an
isolated project/root for a new evaluation. No accepted Stage 07/08 receipts were
found in the audited release roots, so these are conditional workflow commands,
not a claim that the current full cohort can complete them.

Review units, crop offsets, orientation, scale, fixed/moving frames, segmentation
provenance and receipt-bound transform hashes before matching. The existing
mutual-nearest centroid candidate cap is 5 µm; it is not a registration acceptance
tolerance. Report unmatched cells, ties/ambiguities and split/merge limitations.
Full uncertainty propagation and mapped hierarchy comparisons remain specified
work. Preserve both original assay calls, unavailable evidence and disagreement;
concordance is comparison evidence, not ground truth or a tuning target.

## 6. Serial-section MALDI: inspect measurements before registration

The complete IO60 source audit is runnable without approving tissue or neighbour
assignments:

```bash
.venv/bin/python scripts/evaluation/audit_maldi.py \
  --source /data/io60-current/sources/maldi/wilkes_2026-10 \
  --inventory docs/evaluations/dynamic_portability_20261102/inventories/maldi_sections.csv \
  --out "$PORTABILITY_CHECK/maldi"
```

It verifies source SHA256SUMS, glycan/glycogen pairing, spot IDs/counts, coordinate
grid residues, missing/negative values, the existing tissue mask and both
detection denominators for all 32 IO60 sections. Review `summary.json`,
`sections.csv`, `detection_denominators.csv`, every `overlays/<Tissue_id>.png` and
the representative spot rows in `samples/`. The source CSVs also contain material
outside the IO60 roster; only the frozen roster defines this evaluation cohort.

The delivered [MALDI overlays](../outputs/dynamic_portability_20261102/maldi/overlays/)
show instrument coordinates in micrometres, with `log1p` summed intensity used
only for display and cyan marking excluded spots. Their colour mapping is not
MALDI intensity normalization or a fitted registration scale. Acquisition
intensities remain as exported. Physical units beyond exported intensity remain
unavailable where the source does not declare them.

The existing mask admits a spot when summed detected intensity is positive.
Thus all-zero spots are excluded; measured zero at an otherwise admitted spot
still contributes to a feature's finite-value denominator as not detected.
Missing values do not become zero. The audit found 12,486 all-zero excluded spots;
compare all-raster and masked detection fractions before interpreting summaries.

The [MALDI project](../core/multiplex5000/projects/io60_maldi_serial.yaml) runs
00 inventory, 04b ingest, 07 registration and 08 integration. Its existing
automatic gates are 04b and 08; inventory inspection and registration remain
manual review stops. The current source inspection table has no accepted tissue
and exact-neighbour assignment for these sections. The audit does not fill it.
The [sequential MALDI plan](PLAN_2026-10-31_MALDI_SEQUENTIAL.md) defines the review.

After real inspection and accepted inventory, `stage run 04b_maldi_ingest
--project io60_maldi_serial` builds the measurement outputs. `stage run
07_registration --project io60_maldi_serial` then requires those receipts and
inspected neighbours. Current routing is:

| Declared neighbour | Existing registration | Existing object evidence |
|---|---|---|
| PhenoCycler pancreas scan | Inspected Tissue/pLN annotation outline, coarse orientation including mirror, rigid/affine refinement with physical-scale checks; review required | Declared pancreas islet/structure exports; availability and frame must match the receipt |
| Dedicated PhenoCycler pLN | DAPI tissue-mask outline route exists; review required | Active serial object integration does not consume dedicated pLN node exports |
| Xenium | Registration route is unbuilt | No accepted serial integration route |

All are serial-section relationships. There is no cell-identity transfer or
inferred cell type. MALDI samples at 50 µm. Current pancreas object availability
uses a surviving spot centre after erosion by outline RMSE plus half a spot
pitch (25 µm). Signed distance is continuous and negative inside real islets;
no new distance bands are introduced. Missing real islets means distance is
unavailable. Registration quality metrics require review, and are not tuned to
typing agreement or expected biology.

The read-only [pLN source audit](../outputs/maldi_completion_20261102/pln_object_availability.json)
reconciles QuPath entries to exact image/inventory acquisitions: 19 exports with
47 declared annotations, 44 area polygons and three LineStrings. Entry ID 21 is
not a count of 21 entries. Bounding rectangles are declared ROIs, not tissue
silhouettes or inferred follicles/islets. The audit assigns no MALDI neighbour.
It also records the difference between project pixel calibration and the active
registration profile. A future object consumer must explicitly bind geometry
units, frame, source hashes, producer receipt, object class and availability.
Source availability alone does not implement or approve that consumer.

## 7. Record and inspect parameters with KINTSUGI

The recorder uses KINTSUGI's pinned SQLite engine and checks both its commit and
engine-file hash before storage. It accepts only `proposed`, `evaluated`, `failed`
and `blocked` observations. `evaluated` does not mean successful or approved.
Automatic application and recommendation-driven tuning are disabled.

Inspect the populated scientific history without creating another observation:

```bash
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_REPO/outputs/dynamic_portability_20261102/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" history
```

The tracked JSON history is readable without KINTSUGI. The database is under the
selected project's `.claude_workflow/parameter_learning.db`; global history is
not used. Do not edit the database or set an approval flag directly.

For new observations use a new absolute project directory and a JSON payload
conforming to [the record contract](evaluations/dynamic_portability_20261102/recording.md).
Each payload includes:

- Stable observation ID, revision, preceding revision, status, operation and outcome.
- Modality, assay, tissue, panel, relationship type, acquisition/section/segmentation
  IDs, units, producer version, registry/panel hashes and input receipt references.
- A typed subject: markers, genes, feature set, acquisition pair, acquisition or
  geometry. Marker-independent work needs no invented marker.
- Settings labelled `fixed`, `selected_method`, `fitted` or `review_criterion`,
  with their source and units. Missing fit/review/quality values stay null.
- Versioned evidence references with artifact ID, path and SHA256. Keep full
  tables and plots in those artifacts instead of treating diagnostic counts as
  quality scores.

For a runnable round-trip example, copy one existing first revision into a
separate validation project. This preserves its original scientific meaning and
provenance; it does not claim that a new experiment was performed:

```bash
.venv/bin/python - <<'PY'
import json
import os
from pathlib import Path
base = Path(os.environ["PORTABILITY_REPO"])
out = Path(os.environ["PORTABILITY_CHECK"])
history = json.loads((base / "docs/evaluations/dynamic_portability_20261102/parameter_history.json").read_text())
observation = next(row["observation"] for row in history["records"]
                   if row["observation"]["revision"] == 1)
with (out / "observation.json").open("x") as handle:
    json.dump(observation, handle, indent=2, allow_nan=False)
    handle.write("\n")
PY
```

Given that example, or a prepared new schema-valid `observation.json`, the actual
operations are:

```bash
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_CHECK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" \
  record --input "$PORTABILITY_CHECK/observation.json"
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_CHECK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" history
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_CHECK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" export
```

Export writes `parameter_observations.json` inside that project. Repeating an
identical observation ID/revision/content is idempotent. A changed result appends
the next revision with `supersedes` naming its predecessor; it does not overwrite
history. Revised evidence must remain at a new immutable artifact path.

For a retryable storage failure, retain the scientific outputs and submitted
payload. The atomic outbox is `parameter_observation_outbox/`; replay it with:

```bash
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_CHECK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" retry
```

CLI status 2 identifies retryable recording failure. Revisions replay in order;
a predecessor's failure keeps later revisions pending. Invalid schema or changed
content under an existing ID/revision is a conflict requiring correction, not a
transient failure. If even the outbox cannot be written, `pending_path` is null:
retain the original input file and repair storage before retrying. Recording is
separate from scientific production and cannot repair or invalidate its outputs.

### Optional MCP connection

The CLI needs no MCP dependency. The existing optional environment is
`external/.venv-parameter-recording`, with `mcp==1.30.0`. On a new host, use the
[setup instructions](evaluations/dynamic_portability_20261102/recording.md) to
create that environment and a separate pinned KINTSUGI checkout. Do not change an
existing checkout just to make its pin check pass.

Launch the project-scoped stdio server:

```bash
PYTHONPATH=core/multiplex5000 external/.venv-parameter-recording/bin/python \
  -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_CHECK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" serve
```

Configure an MCP client with that executable/argument list and the absolute
`PYTHONPATH`; the project root is fixed at launch. The exposed tools are
`record_parameter_observation`, `get_parameter_observation_history` and
`export_parameter_observations`. They expose no approval, recommendation or
application operation. The evaluation verified a real stdio handshake and
write/retry/read/export; it did not modify global MCP client configuration.
Future reuse requires reviewed compatibility of modality, tissue, panel, units,
method version and relationship type; unknown fields cannot establish a match.

## 8. Verify, diagnose and close the evaluation

Use the applicable existing tests; they complement full-cohort evidence rather
than replace it. Recorder tests exercise the actual pinned SQLite engine:

```bash
.venv/bin/python -m pytest -q -o addopts= -m 'data or not data' \
  core/multiplex5000/tests/current/test_parameter_observations.py
.venv/bin/python -m pytest -q -o addopts= \
  core/multiplex5000/tests/current/test_maldi.py \
  core/multiplex5000/tests/current/test_paired_typing.py \
  core/multiplex5000/tests/current/test_same_slide_pair_stages.py \
  core/multiplex5000/tests/current/test_xenium_only_inventory.py \
  core/multiplex5000/tests/current/test_io60_rna_typing_pilot.py
```

Use `validate_typing_release.py --kind mixed` only after a completed mixed
Stage 06 receipt exists. The current failed refinement needs its primary-table
diagnostics and failure evidence instead. For all completed routes, reconcile
source versus admitted versus called/unavailable counts and native/qualified IDs,
inspect every section's spatial panel, and retain representative cell/spot
details. Report full counts separately from sampled plot rows.

| Symptom | Inspect and record |
|---|---|
| Marker column exists but fitting has zero rows | Acquisition/panel declaration, Stage 03/05 availability, finite values and minimum reference population; preserve explicit missingness and the failure. |
| Missing RNA gene or unsupported class resolution | Actual panel feature names, measured zero versus absent gene, tissue identity cores, controls and unresolved/broad outcome. Do not import a protein label. |
| Missing hormone/lineage reference | Source reference-population construction and method-selection branch; document absence instead of transferring pancreas references to another tissue. |
| Duplicate/mismatched IDs or counts | Native UUID/cell ID, acquisition-qualified ID, segmentation version, source QC and join cardinality. Never repair by row order. |
| Cells without nuclei appear admitted | Original nucleus fields and eligibility mask; separate segmented and eligible universes and retain modality-specific exclusions. |
| Coordinates disagree | Units, pixel calibration, crop origin, orientation, fixed/moving frame and transform receipt/hash before any distance comparison. |
| Stale transform or absent outline | Exact source and transform bindings, available reviewed geometry and explicit unavailable reason. Recompute only under the existing specified method. |
| Ambiguous same-slide matches | Candidate distances, ties, unmatched cells and segmentation splits/merges; record unsupported uncertainty assessment. |
| MALDI section cannot advance | Inspected tissue and exact neighbour declaration, source/ingest receipts, implemented neighbour route and accepted registration; donor or Tissue_id text is insufficient. |
| MALDI object table is empty | Declared real area objects, correct frame, source receipt and erosion eligibility. A DAPI outline or rectangular node ROI is not an inferred islet. |
| Recorder rejects input or storage | Pin/hash mismatch, schema conflict versus transient write failure, outbox path and predecessor revision. Preserve science artifacts and retry only the recording. |

Close with a report, machine-readable dependency/readiness matrices, a complete
acquisition/pair ledger, versioned diagnostics and hashes, exported parameter
history, and prioritized implementation-ready remediation. Keep computational
fidelity, geometric validity and biological applicability as separate conclusions.
A documented unsupported route can have a complete evaluation disposition while
remaining scientifically blocked; it must not be reported as a completed typing
or integration result.
