# Running and reviewing the remediated portability routes

Use this guide for the implementation after the user's blocker-fix authorization.
The [current report](evaluations/blocker_remediation_20261102/README.md),
[readiness matrix](evaluations/blocker_remediation_20261102/readiness.json) and
[201-acquisition ledger](evaluations/blocker_remediation_20261102/acquisitions.json)
separate completed candidates, active checks, missing evidence and scientific
approval. The [original guide](PORTABILITY_EVALUATION_USER_GUIDE.md) documents the
frozen pre-remediation evaluation; its statements about absent capabilities are
historical.

A passing computational gate does not validate biological transfer. Keep donor,
physical section, acquisition, segmentation, cell and MALDI spot identities
separate. No relationship is inferred from a donor name. No commands in this
guide activate a production release or approve a registry.

## Environment and read-only inspection

From the repository root on the workstation:

```bash
export PORTABILITY_REPO="$PWD"
export PYTHONPATH="$PORTABILITY_REPO/core/multiplex5000"
export MPLX5K_SITE_FILE="$PORTABILITY_REPO/core/multiplex5000/sites/sn4622120742.json"
.venv/bin/python -m multiplex5000.cli status --project io60
.venv/bin/python -m multiplex5000.cli status --project io60_maldi_serial
.venv/bin/python scripts/evaluation/verify_portability_delivery.py
```

The scientific interpreter is the existing `.venv/bin/python`; do not change
package versions merely to make a historical comparison run. The CODEX replay
has its own explicit Python/package/thread contract. HiPerGator uses its own
site file and source paths; workstation paths do not transfer automatically.
Record the commit, source diff/hashes, environment, input receipts and selected
roster for every new run. Missing artifacts are unavailable evidence.

Read the completed pLN and full 25 candidate status directly:

```bash
.venv/bin/python -m multiplex5000.cli status \
  --project outputs/blocker_remediation_20261102/config/pln_full19_candidate.yaml
.venv/bin/python -m multiplex5000.cli status \
  --project outputs/blocker_remediation_20261102/protein/protein_availability_remediation.yaml
```

Large outputs are external or gitignored. JSON reports reference their exact
paths and hashes. Stage receipts bind input/output content; an old filename or
existing folder alone does not establish compatible provenance.

## Create a fresh candidate

Do not force a rerun into the dated evidence roots. Make an isolated directory
and copy the appropriate project, changing its ID, release and output root.
This reusable example prepares the RNA-only project:

```bash
export PORTABILITY_WORK="$PORTABILITY_REPO/outputs/portability-remediation-recheck-1"
mkdir -p "$PORTABILITY_WORK"
.venv/bin/python - <<'PY'
import os
from pathlib import Path
import yaml
repo = Path(os.environ['PORTABILITY_REPO'])
work = Path(os.environ['PORTABILITY_WORK'])
source = repo / 'core/multiplex5000/projects/io60_xenium_only_6512_6618.yaml'
project = yaml.safe_load(source.read_text())
project.update(id='rna_only_recheck', release='rna-only-recheck-1',
               releases_root=str(work / 'releases'))
with (work / 'rna_only_recheck.yaml').open('x') as handle:
    yaml.safe_dump(project, handle, sort_keys=False)
PY
.venv/bin/python -m multiplex5000.cli graph-spec \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml"
```

The YAML ID must equal its filename stem. For another route, use its listed
source project and a distinct ID/release; retain the full explicit acquisition
roster and scientific settings. Preserve upstream reuse only when source and
measurement contracts match. Do not copy acceptance files to manufacture review.
`stage run` reuses an already completed stage; changing source code does not
silently rebuild it. Use a new candidate root for a fresh evaluation.

Capture the implementation version before execution, separately from the
[historical scientific producer pins](evaluations/dynamic_portability_20261102/pins.json):

```bash
git rev-parse HEAD > "$PORTABILITY_WORK/producer_commit.txt"
git status --short > "$PORTABILITY_WORK/producer_status.txt"
git diff --binary HEAD -- core/multiplex5000 scripts/evaluation \
  > "$PORTABILITY_WORK/producer_tracked_changes.patch"
git ls-files --others --exclude-standard -- core/multiplex5000 scripts/evaluation \
  > "$PORTABILITY_WORK/producer_untracked_files.txt"
sha256sum "$PORTABILITY_WORK/rna_only_recheck.yaml" \
  > "$PORTABILITY_WORK/project.sha256"
.venv/bin/python - <<'PY' > "$PORTABILITY_WORK/environment.json"
import importlib.metadata
import json
import platform
import sys
packages = {d.metadata['Name']: d.version for d in importlib.metadata.distributions()
            if d.metadata.get('Name')}
print(json.dumps({'python': platform.python_version(), 'executable': sys.executable,
                  'packages': packages}, indent=2, sort_keys=True))
PY
```

Use the actual selected project in the `sha256sum` command. Archive the exact
consumed source files too, including relevant untracked files: `git diff` does
not contain their bytes. The completed follow-up keeps producer snapshots and
source hashes beside its receipts. Never relabel an old output with the current
commit simply because it can be opened by newer code.

Commit reviewed implementation/tests and reviewed documentation/evidence in
separate, explicit groups. Inspect `git diff --cached --stat` and
`git diff --cached --check` before a commit; stage named owned files rather than
all untracked content. Generated outputs, source images, SQLite/runtime archives
and private user notes do not belong in those commits. Preserve each candidate's
actual producer version and dirty-source snapshot even after a later clean
commit. Committing code does not activate a release or rewrite its receipts.

Each stage's `receipt.json`, `gate/gate.json` and `gate/gate.html` explain what ran
and what passed. Use `gate show STAGE --project PROJECT` to inspect it. A review
acceptance is an action based on the actual evidence and existing criteria:

```bash
.venv/bin/python -m multiplex5000.cli gate show 00_inventory \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml"
# Only after the actual inventory review, with the real reviewer/evidence values:
.venv/bin/python -m multiplex5000.cli gate accept 00_inventory \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml" \
  --by "YOUR_REVIEWER_NAME" --note "REVIEW_EVIDENCE_ID"
```

Existing automatic computational gates are declared in the project YAML; they
are not biological approval. Missing acceptance blocks dependent stages. No new
scientific numerical criterion is supplied by this guide.

## Protein typing

### CODEX original-method replay and separate v20 bindings

The repaired validator uses the pinned external port and recovered producer
journal, preserves the producer's input parsing, and decodes expected CSV
references losslessly. It retains exact equality. Run into a new directory with
the verified runtime, before numerical imports:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/smith6jt/miniconda3/envs/scvi-env/bin/python \
  scripts/evaluation/codex_replay.py \
  --evidence outputs/handoff_completion_20261102/codex/journal_proof.json \
  --snapshot outputs/dynamic_portability_20261102/typing_port \
  --out "$PORTABILITY_WORK/codex-exact-replay"
.venv/bin/python scripts/evaluation/audit_codex_v20.py \
  --recovery-manifest outputs/blocker_remediation_20261102/codex/v20_recovered/recovery.json \
  --out "$PORTABILITY_WORK/codex-v20-bindings"
```

The explicit contract is `scripts/evaluation/codex_replay_contract.json`. A
runtime or content mismatch must be investigated; do not loosen equality or
edit the expected tables. The completed follow-up passed 60 exact assertions
across all 20 images. The 24 audited v20 bindings are independently verified;
that does not recover every v20 execution input or activate v20.

The later [received CODEX isolation handoff](evaluations/blocker_remediation_20261102/codex_received_handoff.md)
is separately verified source evidence: all 17 package-file checks pass, but its
code and thresholds were not applied to these runs. It proposes separate fixed
calibration-reference and live structure-membership contexts. The receiving
workflow still needs explicit baseline binding, a distinct structure-update route,
registry role mapping and whole-cohort invariants before that adaptation can run.
Keep its package/source pins separate from the actual producer; supplied audit
results are not a new replay executed by this project.

The separate [component/Ki67 handoff](evaluations/blocker_remediation_20261102/codex_received_component_handoff.md)
is also received and not applied. Its 17 payload files and 214 regular source
files were verified locally. The scientific counts are supplied results, not
reruns here; reported release/object failures remain, and HiPerGator delivery
was not independently verified. It changes no active producer or cohort scope.

### IO60 pancreas and embedded pLN, full 25

Use a fresh copy of
`outputs/blocker_remediation_20261102/protein/protein_availability_remediation.yaml`
(or the tracked `io60_mixed_protein_full_candidate.yaml` with the same explicit
source contract). This route consumes existing accepted measurement,
normalization, geometry, registry and acquisition-run evidence; it does not
replace canonical compartment definitions. The mixed candidate contains separate
pancreas and embedded pLN scopes. Its per-acquisition refinement availability
is derived from the actual source measurements, not merely a present column.

The completed dated candidate covers 39,394,852 segmented cells, including
31,385,695 eligible pancreas primary/refinement cells and a separate embedded
pLN scope. All 25 sections passed independent identity/count/primary-call checks
and inspection of their call/availability maps and 189 native-cell examples.
See the [full result and limits](evaluations/blocker_remediation_20261102/README.md#mixed-protein-all-25-completed).

For the fresh copy named `protein_recheck.yaml`, capture the current producer
bytes, configuration and complete declared input roster before running Stage06:

```bash
.venv/bin/python scripts/evaluation/capture_protein_launch.py \
  --project "$PORTABILITY_WORK/protein_recheck.yaml" \
  --base-snapshot outputs/blocker_remediation_20261102/protein/producer_snapshot_attempt2.json \
  --snapshot "$PORTABILITY_WORK/protein_producer_snapshot.json" \
  --archive "$PORTABILITY_WORK/protein_producer_sources" \
  --bindings "$PORTABILITY_WORK/protein_input_bindings.json" \
  --attempt 1
```

The base snapshot supplies the source-file roster; the capture hashes and archives
the current bytes rather than reusing its historical hashes. All three output
targets must be fresh. This evaluation helper neither executes a stage nor
accepts a gate. Keep the project and consumed sources/inputs unchanged between
capture and execution. If capture fails partway through, retain those artifacts
and use new output paths for the next attempt.

Then run the candidate:

```bash
.venv/bin/python -m multiplex5000.cli stage run 06_typing \
  --project "$PORTABILITY_WORK/protein_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 06_typing \
  --project "$PORTABILITY_WORK/protein_recheck.yaml"
```

Inspect `pancreas/refinement/refine_log.json`,
`pancreas/refinement/representation_provenance.json`, mixed tables and
`resolved_parameters.json`. Unavailable refinement preserves the primary call,
records its reason and does not demote a cell merely because evidence is absent.
Check the exact fit population, carried markers, source availability, marker
compartments, scales and numerical units. The full 25 representation fit is
recorded as full 25; it is not relabelled a single-donor fit.

The representation fit now uses carried-marker evidence in supported dimensions
from the four incomplete-panel sections that the previous union-column
complete-case filter discarded. The fixed neighbour count, concordance, bound
and seed criteria are unchanged. This is a newly fitted refinement candidate:
the required exact comparison is primary-call parity. It does not establish
identical refinement labels to an unexecuted baseline or biological approval.

The published producer also warns that `top_identity_marker` reads `FOXP3` on
4,487,012 cells (14.30%) under the cohort scale. This is a known diagnostic-field
limitation outside the availability repair; it is not evidence that those cells
are regulatory T cells. Preserve this warning when reviewing or exporting the
candidate, and do not use that field as a Treg identity call.

Seven plasma-cell donor/rule cases also remain unavailable in the immune
component. The raw nucleus-admission audit used size/mtime binding; current
aligned/primary checks do not constitute a new full raw-image/export hash audit.
Existing coordinate-conversion differences remain documented. Retained inputs
still resolve into the v3 checkout and must be copied and verified before an
active release, as required by the published gate. No Stage06 acceptance was
created by the computational evaluation.

The independent verifier checks a completed release against the preserved
primary tables. For the fresh copy with `release: protein-recheck-1`:

```bash
.venv/bin/python scripts/evaluation/verify_protein_remediation.py \
  --release "$PORTABILITY_WORK/releases/protein-recheck-1" \
  --baseline outputs/dynamic_portability_20261102/protein_primary_expected_universe.json \
  --producer-snapshot "$PORTABILITY_WORK/protein_producer_snapshot.json" \
  --input-bindings "$PORTABILITY_WORK/protein_input_bindings.json" \
  --out "$PORTABILITY_WORK/protein-verification"
.venv/bin/python scripts/evaluation/check_protein_producer_snapshot.py \
  --run "$PORTABILITY_WORK" \
  --snapshot "$PORTABILITY_WORK/protein_producer_snapshot.json" \
  --stage "$PORTABILITY_WORK/releases/protein-recheck-1/06_typing" \
  --out "$PORTABILITY_WORK/producer_stability.json"
```

The snapshot and input-binding files must belong to that exact attempt and must
be captured before execution. The independent verifier compares the complete
launch input dictionaries and settings against the final receipt, and checks
the unchanged project and snapshot hashes. The source checker reports a
conservative static import closure; it does not claim to trace every runtime
read or continuously monitor files. The delivery guard also requires a
post-publication `producer_stability.json` bound to the actual published runner
manifest, attempt snapshot and current source hashes; the earlier running check
cannot substitute for this final artifact.

The corrected dated attempt uses
`outputs/blocker_remediation_20261102/protein/producer_snapshot_attempt2.json`
and `outputs/blocker_remediation_20261102/protein/attempt2_input_bindings.json`.
Its launch record contains 268 input declarations: 253 existing inputs and
15 optional no-file declarations. These records cannot be reused as the launch
evidence for a new project copy. Attempt2 used the original ad hoc capture recipe;
the new reusable CLI did not create or rewrite those dated records.

After the independent audit passes, prepare the review material:

```bash
.venv/bin/python scripts/evaluation/build_protein_review_materials.py \
  --verification "$PORTABILITY_WORK/protein-verification"
```

This requires a fresh `review_materials/` directory beneath the verification
folder. It places existing call/availability maps into contact sheets without
resampling and selects the first native cell in every review state from each
full mixed table. The maps retain their existing deterministic 12,000-cell
samples. Inspect every section's maps and native-cell examples, then record
actual findings and evidence hashes in the separate inspection ledger. The
generated manifest keeps `visual_review` pending; preparing figures does not
perform inspection, review raw marker images, approve biology or accept Stage06.

A complete receipt plus successful identity/count/source verification is required
before reporting the full candidate complete. Review the per-section diagnostic
and unavailable states; no disease ordering or expected class frequency is a
repair target.

The verifier reads reported audit failures as well as the process exit status.
An absent `pln_blocks:<donor>` input can be classified as an optional no-file
declaration only when the exact receipt declares it absent, the consumed producer
declares it `required=False`, the full identity audit finds zero embedded pLN
cells, and the path remains absent. Such rows stay in the raw unavailable-hash
audit and receive a separate optional-absence disposition; they are never counted
as successful file hashes. A missing required or unclassified input fails the
verification.

An unsuccessful publication has no completed candidate receipt. The first full
25 attempt was rejected for an undeclared `refine.py` read. Although that source
declaration is now fixed, missing saved input bindings prevent safe publication
of the old outputs. The corrected attempt must execute normally with its own
launch input hashes and producer snapshot; retain the failed attempt separately.
Do not manufacture a receipt from files merely present after a failure.

### Dedicated pLN, all 19

The completed full 19 project is
`outputs/blocker_remediation_20261102/config/pln_full19_candidate.yaml`. It lists
all 19 exact images and uses raw QuPath review-candidate measurements, a fresh
shared 55-marker UniFORM fit and section-local candidate evaluation. It produces
no final cell classes and is not equivalent to canonical own-compartment REDSEA.
Audit raw sources with the existing CLI:

```bash
.venv/bin/python -m multiplex5000.cli sources pln
```

For a fresh copy named `pln_recheck.yaml`, execute the full upstream chain:

```bash
.venv/bin/python -m multiplex5000.cli stage run 03_protein_measurement \
  --project "$PORTABILITY_WORK/pln_recheck.yaml"
.venv/bin/python -m multiplex5000.cli stage run 05_normalization \
  --project "$PORTABILITY_WORK/pln_recheck.yaml"
.venv/bin/python -m multiplex5000.cli stage run 06_typing \
  --project "$PORTABILITY_WORK/pln_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 06_typing \
  --project "$PORTABILITY_WORK/pln_recheck.yaml"
```

The source project retains its existing automatic Stage03/05 computational gates;
Stage06 remains a review stop. Do not reuse the old two-section fit as a full 19
fit. The default `protein_cutoff_fraction: 0.1` is now actually consumed. Its
validation accepts finite fractions in `[0,1]`; range validity is not scientific
approval to choose a different threshold. Unknown or missing evidence remains
unavailable. The completed run retained all 10,505,964 cells, including 622,153
nucleus-ineligible cells, and recorded 17 existing normalization review flags.

### Spleen PC67 and PC33

The current `spleen_primary.yaml` declares measurement/normalization only. The
historical eight-section PC67 typing evidence comes from separately pinned
producer `a9549cc`; the external port's P1–P9 proposals are not automatically
approved or installed by that history. The four PC33 acquisitions remain in the
ledger with missing reviewed lineage/registry requirements. There is no supported
command that makes a PC33 final typing result from a PC67 registry.

The bounded `spleen_pilot.yaml` supplies protein and RNA gene evidence for the
same-slide workflow. Preserve its independent assay evidence and its scope;
RNA gene states alone are not full spleen RNA lineage calls.

The ten previously unverified large inputs now pass their declared historical
hash modes: eight full geometry hashes and two sampled image fingerprints. The
audit also reads and records all current full-file bytes. Do not compare a
historical sampled fingerprint directly with a full-file SHA-256 or describe
that difference as content drift. The original two sampled image bindings cannot
prove historical full-byte identity. The
[declared-mode verification](../outputs/blocker_remediation_20261102/hpg_maldi_xenium/spleen_algorithm_verification/verification.json)
records the manifest, algorithms and current hashes. Repetition with
`scripts/evaluation/verify_recorded_inputs.py` requires the exact `--receipt`,
its hash-bound `--manifest`, and a fresh `--out` on a host where those paths exist.
For the evaluated HiPerGator source paths, run on HiPerGator:

```bash
export SPLEEN_VERIFY_ROOT=/blue/maigan/smith6jt/blocker-remediation-maldi-xenium-20261102
PYTHONPATH="$SPLEEN_VERIFY_ROOT/snapshot/core/multiplex5000" \
  /blue/maigan/smith6jt/mplx5k-io60-rna-20261027/venv-current/bin/python \
  "$SPLEEN_VERIFY_ROOT/verify_recorded_inputs.py" \
  --receipt "$SPLEEN_VERIFY_ROOT/spleen_receipt.json" \
  --manifest "$SPLEEN_VERIFY_ROOT/spleen_qupath_manifest.json" \
  --out "$SPLEEN_VERIFY_ROOT/spleen-recheck-1"
```

## Xenium-only and full-cohort RNA

The RNA-only project includes exact acquisitions `xenium:0041323:Region_1` and
`xenium:0041326:Region_1`, with no protein measurement, normalization or pair
requirement. The full roster is in `io60_xenium_full_same_release.yaml`; copy it
to a new release before execution. Its 65 acquisitions remain 65 even when five
are unavailable under existing damaged/low-signal rules.

For the fresh RNA-only copy prepared above:

```bash
.venv/bin/python -m multiplex5000.cli stage run 00_inventory \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml"
# Inspect and accept the actual inventory before continuing.
.venv/bin/python -m multiplex5000.cli stage run 04_rna_ingest \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 04_rna_ingest \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml"
# Inspect and accept the actual ingest before continuing.
.venv/bin/python -m multiplex5000.cli stage run 06_typing \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 06_typing \
  --project "$PORTABILITY_WORK/rna_only_recheck.yaml"
```

The setting `pancreas_naming_policy` accepts the existing `shared_er_exact_v1`
(default) or `native` policy; unsupported values are rejected before worker
execution. This naming setting governs the pancreas worker; pLN retains its
tissue-specific worker. The resolved artifact includes the actual ordered gene list/hash,
identity-core coverage, selected tissue worker, fixed controls and emitted fits.
Absent genes remain different from measured zeros. Unemitted pLN embedding fit
values remain unavailable; the fixed detection normalization target is not
substituted for an unreported embedding value.

Nucleus admission and damaged-region exclusions remain in force. The 6457
PhenoCycler exclusion does not exclude its Xenium or MALDI acquisitions. The
frozen full RNA evaluation completed 60 candidates with five unavailable
(three damaged regions and two below the recorded minimum median detected-gene
criterion); current
configuration wiring was tested independently. Two selected donors support
provisional research comparisons, not population-level disease conclusions.

## Same-slide registration and candidate correspondence

The source manifest must establish the exact physical section, acquisitions,
coordinate frames, crop offsets, orientation and segmentation. Type the assays
independently. A joint typing output or matching donor label is insufficient.

For an isolated copy of `spleen_pilot.yaml`, retain the exact selected pair and
its independent upstream receipts. Once Stage06 is actually reviewed, the
existing stage commands are:

```bash
.venv/bin/python -m multiplex5000.cli stage run 07_registration \
  --project "$PORTABILITY_WORK/spleen_pair_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 07_registration \
  --project "$PORTABILITY_WORK/spleen_pair_recheck.yaml"
# Inspect the exact geometry and create its genuine review acceptance before Stage08.
.venv/bin/python -m multiplex5000.cli stage run 08_integration \
  --project "$PORTABILITY_WORK/spleen_pair_recheck.yaml"
```

Stage08 now writes `candidate_links.parquet`, `spatial_candidates.parquet`,
`cell_link_dispositions.parquet`, `pairing_counts.json`, `linking_metadata.json`
and `resolved_parameters.json`. Each native cell has a disposition. Equal nearest
distance ties remain ambiguous; unmatched and nonmutual alternatives remain
visible. Both original assay calls/states/confidences are retained. Multiple
centroid candidates are not asserted segmentation splits/merges without overlap
evidence.

The existing matching cap is 5 µm. It is not a registration quality cutoff. Optional
`registration_uncertainty` transform metadata must supply a finite nonnegative
`bound_um`, exact fixed `frame`, `kind: radial_displacement_bound`, and actual
`evidence_ref`. Without that evidence, uncertainty/stability remains unavailable.
The [GEOM-2 report](evaluations/blocker_remediation_20261102/geom2.md) specifies the
conservative interval and stability calculation; it does not estimate a bound
from candidate agreement.

An optional `same_slide_hierarchy_map` declares exact tissue, comparison level,
independent call columns, mappings and any availability-state filters. Comparisons
remain `agree`, `disagree`, `unmapped` or `unavailable`; no fused identity is made.
Do not tune registration or calls to improve agreement.

The [batch coordinator](evaluations/blocker_remediation_20261102/geom1.md) runs
each explicitly configured pair independently and always accounts for the full
hash-bound roster. Its real preflight accounts for all 12 declarations: eight
blocked spleen pairs and four unavailable pancreas section adapters, with
unresolved pLN relationships retained separately. The combined coordinator and
matching regression suite passed 98 tests.

Run a read-only preflight against the supplied evaluation configuration:

```bash
.venv/bin/python -m multiplex5000.same_slide_batch preflight \
  --config core/multiplex5000/projects/same_slide_batch_evaluation.json \
  --out "$PORTABILITY_WORK/same-slide-preflight"
```

Exit 1 means the cohort is incomplete, with exact per-pair reasons in
`ledger.json`; exit 2 means invalid configuration. The config requires one entry
for every frozen `pair_id`, an isolated `workspace_root`, and distinct pair
project release roots. `project: null` retains an unsupported route. The default
`required_evidence: independent_calls` requires actual independent calls and an
explicit hierarchy map. Gene-only RNA cannot satisfy it. An explicitly requested
`spatial_gene_evidence` mode has narrower scope and does not establish independent
RNA identity. No mode accepts gates or invents relationships.

To execute new candidates, prepare a fresh batch config and distinct pair project
copies under the new workspace, retaining the exact roster/hash and adjusting
relative project paths. Supply each pair's actual accepted upstream evidence:

```bash
.venv/bin/python -m multiplex5000.same_slide_batch registration \
  --config "$PORTABILITY_WORK/same_slide_batch.json" \
  --out "$PORTABILITY_WORK/same-slide-registration-ledger"
# Review each actual Stage07 result separately before accepting its gate.
.venv/bin/python -m multiplex5000.same_slide_batch integration \
  --config "$PORTABILITY_WORK/same_slide_batch.json" \
  --out "$PORTABILITY_WORK/same-slide-integration-ledger"
```

Every ledger directory must be new and outside the pair release roots. Ready
pairs can progress while other entries remain blocked; a partial ledger is not
full-cohort completion. Reused outputs are checked for current hashes/settings;
the coordinator neither forces replacement nor copies acceptance. IO60 mixed
donor outputs require a source-bound section adapter and cannot be passed directly
to the spleen cell-key contract. pLN physical-section links remain unresolved.

## Serial-section MALDI

The user describes the usual cutting order as **Xenium → PhenoCycler → MALDI**,
qualified by "For the most part". The [context record](evaluations/blocker_remediation_20261102/relationship_context.json)
preserves the exact assertion. Within a confirmed three-section series following
that order, PhenoCycler immediately precedes MALDI. Exceptions, exact section and
acquisition mappings, and spacing remain unspecified. Keep this general serial
context separate from the inventory's exact same-slide relationships; use it as
context when inspecting a declared series, without filling an exact neighbour
field from donor or filename labels.

MALDI is serial-section evidence at 50 µm sampling. The 32 sections/18 donors and
1,666,450 spots retain exact paired glycan/glycogen keys. Zero is not detected;
NaN is unavailable. No intensity normalization or new distance bands are added.
The new [scaling-plot viewer](../outputs/blocker_remediation_20261102/maldi_scaling/index.html)
shows all 32 sections as linear intensity, `log1p` intensity and the existing
signal mask. The original [spatial plots](../outputs/dynamic_portability_20261102/maldi/overlays/)
remain frozen. Plot display scaling is distinct from physical registration scale
and from intensity normalization.

Regenerate the display comparison only into a fresh folder:

```bash
.venv/bin/python scripts/evaluation/plot_maldi_scaling.py \
  --source /data/io60-current/sources/maldi/wilkes_2026-10 \
  --inventory docs/evaluations/dynamic_portability_20261102/inventories/maldi_sections.csv \
  --out "$PORTABILITY_WORK/maldi-scaling"
```

Current registration supports PhenoCycler pancreas outlines, dedicated pLN DAPI
masks and exact-region Xenium morphology. Xenium uses instrument-declared DAPI,
OME axes, micrometre calibration, region identity and native crop/orientation.
Whole-slide repeated thumbnails cannot substitute for region-local geometry.
The final indexed HiPerGator validation accounted for all 65 Xenium acquisitions:
62 source frames validated and three recorded damaged regions unavailable,
including both RNA-only acquisitions among the validated sources. Two explicitly
relocated directories are bound through the exact acquisition/output index and
matching experiment hash. The intermediate 60/65 result remains historical
evidence of the repaired reader check. This selects no
MALDI neighbour and accepts no transform.

Read-only source validation into a new root:

```bash
.venv/bin/python scripts/evaluation/validate_maldi_xenium_sources.py \
  --inventory-release /data/io60-current/releases/io60-maldi-candidate-1 \
  --out "$PORTABILITY_WORK/maldi-xenium-sources"
.venv/bin/python scripts/evaluation/prepare_maldi_object_sources.py \
  --audit docs/evaluations/dynamic_portability_20261102/pln_object_availability.json \
  --out "$PORTABILITY_WORK/maldi-object-sources"
```

On HiPerGator supply `--source-index` with the exact declared source CSV and use
its inventory-release path/site environment. A failed local source lookup is not
proof that an acquisition is scientifically absent.

A fresh `io60_maldi_serial.yaml` copy needs the inspection input
`maldi_tissue_inspection.csv`. Prepare an unchanged copy of the existing
[32-row template](../core/phenocycler/data/integration/maldi_tissue_inspection.csv)
and bind its absolute path in the new project:

```bash
.venv/bin/python - <<'PY'
import os
from pathlib import Path
import yaml
repo = Path(os.environ['PORTABILITY_REPO'])
work = Path(os.environ['PORTABILITY_WORK'])
template = repo / 'core/phenocycler/data/integration/maldi_tissue_inspection.csv'
inspection = work / 'maldi_tissue_inspection.csv'
with inspection.open('xb') as handle:
    handle.write(template.read_bytes())
project = yaml.safe_load((repo / 'core/multiplex5000/projects/io60_maldi_serial.yaml').read_text())
project.update(id='maldi_recheck', release='maldi-recheck-1',
               releases_root=str(work / 'releases'))
project.setdefault('sources', {})['maldi_tissue_inspection'] = str(inspection)
with (work / 'maldi_recheck.yaml').open('x') as handle:
    yaml.safe_dump(project, handle, sort_keys=False)
PY
```

This command leaves every assignment blank. The project source key
`sources.maldi_tissue_inspection` overrides the site source of the same name;
when neither is supplied, the inventory reads the repository template. Keep the
five CSV columns and all 32 tissue IDs:

| Column | Required content |
|---|---|
| `tissue_id` | Exact prefilled MALDI export ID; do not parse it to infer tissue or section order. |
| `tissue` | For the IO60 pancreas/pLN cohort, `pancreas` or `pancreatic_lymph_node`; `pLN` is an accepted alias. The profile parser also recognizes `spleen`, which does not itself establish a supported MALDI tissue route. |
| `neighbour_acquisition` | Exact acquisition ID from the inventory's `acquisitions.csv`, as titled in `maldi_sections.png`; it must be an actual inspected neighbour from the same donor. A matching donor alone is insufficient. |
| `inspected_by` | The person who actually inspected the section and supporting records. |
| `note` | Evidence reference and observed orientation, pieces, damage or an exception to the usual cutting order. This is retained free text, not a computed approval. |

Compare the ion images with the candidate slides in
`/data/io60-current/releases/io60-maldi-candidate-1/00_inventory/maldi_sections.png`
and its adjacent `acquisitions.csv`, then record only verified assignments in the
fresh copy. An unknown relationship remains blank and blocked. The
[original MALDI plan](PLAN_2026-10-31_MALDI_SEQUENTIAL.md) records historical scope;
it supplies no missing acquisition assignments. Stage00 checks unknown or
duplicate IDs, unknown tissue, donor mismatch, missing inspector and incomplete
rows. Its review gate remains necessary after those integrity checks pass.

Set `sources.maldi_object_sources` to
the explicit generated object manifest for object support; set
`sources.xenium_source_index` when source relocation is necessary. This manifest
binds acquisition/frame, units/calibration, selected classes, native object IDs,
export hashes and source receipts. An observation receipt records an existing
export; it does not invent its missing historical producer.

After actual inspection, the scientific sequence is:

```bash
.venv/bin/python -m multiplex5000.cli stage run 00_inventory \
  --project "$PORTABILITY_WORK/maldi_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 00_inventory \
  --project "$PORTABILITY_WORK/maldi_recheck.yaml"
# Genuine inventory acceptance is required before ingest.
.venv/bin/python -m multiplex5000.cli stage run 04b_maldi_ingest \
  --project "$PORTABILITY_WORK/maldi_recheck.yaml"
.venv/bin/python -m multiplex5000.cli stage run 07_registration \
  --project "$PORTABILITY_WORK/maldi_recheck.yaml"
.venv/bin/python -m multiplex5000.cli gate show 07_registration \
  --project "$PORTABILITY_WORK/maldi_recheck.yaml"
# Review and accept the exact registration before object integration.
.venv/bin/python -m multiplex5000.cli stage run 08_integration \
  --project "$PORTABILITY_WORK/maldi_recheck.yaml"
```

All 32 currently lack new accepted real transforms/object results. Source masks
and synthetic tests do not fill the missing inspected-neighbour rows. Use the
existing registration metrics/review contract; no Dice or RMSE cutoff is invented.

Inspect `registration.csv`, each transform/frame, overlay and uncertainty. Stage08
binds current accepted receipts/gates and the exact transform/source files. It
retains unavailable registrations in its declared-section accounting. Per-spot
outputs carry coordinates, uncertainty, object membership and continuous signed
islet distance where real islet evidence exists. Overlapping memberships have
no unique spot object ID; object summaries are explicitly nonexclusive.

Object eligibility uses the existing erosion by outline RMSE plus half the spot
pitch. This is a geometric sensitivity/eligibility margin, not a statistically
calibrated confidence bound. It is distinct from the explicitly supplied
uncertainty bounds supported by the same-slide matcher. LineStrings and invalid
areas remain unavailable entries. The 19 pLN source
exports contain 47 objects, 44 areas; these are recorded annotations, not inferred
follicles/islets/cell types. Main detection fractions use finite signal-mask
spots, retain measured zeros and exclude missing values. `all_raster_*` fields
show the separate all-acquired-spot sensitivity comparison. An all-zero spot
excluded by the existing mask is therefore visible without silently redefining
the principal denominator.

## Parameters: inspect, record, retry and export

Changed stages write `resolved_parameters.json`, bound by their enclosing receipt.
It identifies fixed values, selected methods, fitted values, unavailable evidence,
units, acquisition context, producer sources and input hashes. Recording failure
returns a retryable status and retains scientific outputs. This boundary artifact
is not itself the observation adapter's input schema: convert it to a complete
observation with context and evidence references before calling `record`.

KINTSUGI remains pinned at `716d1ebd87db916ab04910c6374e04d2d9ea508e`. The adapter
validates its checkout and engine-file hash. The project-only database lives at
`PROJECT/.claude_workflow/parameter_learning.db`. It never uses global history.
Allowed observation statuses are `proposed`, `evaluated`, `failed` and `blocked`;
none becomes an approved recommendation. Missing quality measurements remain
null. Typed subjects support genes, marker sets, features, geometry and acquisition
pairs without inventing a marker for geometric work.

Record the scope of a failure precisely. The bundle uses a separate failed
`protein_stage_publication` observation for the first full protein attempt, with its exact cohort,
error and six evidence artifacts. Its quality measurements remain unavailable,
and it has no complete outer receipt. This is not evidence that the scientific
parameter values failed, and its large unreceipted outputs remain unpublished.

The [completed remediation history](evaluations/blocker_remediation_20261102/parameter_history.json)
has 159 observations/records and zero approved records. All retries were
idempotent; 58 distinct evidence files are preserved in the project archive.
The [actual scientific MCP check](../outputs/blocker_remediation_20261102/parameter-history/scientific_mcp_check.json)
passed protein-observation read, duplicate retry and exact 159-record export.
Global history remained absent and unchanged, with recommendations and automatic
application disabled.

Inspect that history directly:

```bash
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_REPO/outputs/blocker_remediation_20261102/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" history
```

For a new schema-valid observation in a fresh project:

```bash
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_WORK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" \
  record --input "$PORTABILITY_WORK/observation.json"
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_WORK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" retry
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_WORK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" history
.venv/bin/python -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_WORK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" export
```

Export writes `parameter_observations.json` inside that project. Stable observation
ID/revision/content retries are idempotent. A later revision names its predecessor;
changed content under the same revision is a conflict, not a retry. CLI exit 2
means retryable recording failure. The outbox replays predecessor revisions first;
if the outbox itself cannot be written, preserve the original payload and repair
storage before retry. A database failure must not delete or rerun the science.
See the [observation schema](evaluations/dynamic_portability_20261102/recording.md)
for required context, units, method version, registry/panel hashes and evidence
references.

For the completed remediation bundle, the dedicated recorder validates each
Stage06 receipt/stage/release/output binding before recording. Protein records
also require the strict completed-cohort check and retain the independent
verification report, raw receipt audit, final producer-stability evidence and
inspection ledger. It archives every
referenced evidence file and input receipt under
`PROJECT/evidence/<sha256>/<basename>`, preserving the original artifact ID and
hash. A later revision therefore retains the earlier bytes even if a working
report changes. These archives are project runtime evidence, covered by the
existing `outputs/` ignore rule. For separately prepared generic CLI observations,
provide equivalently versioned evidence references; the generic payload does not
gain this bundle-specific archival step by itself.

Protein observations also retain emitted primary islet levels, stratum
diagnostics, selected method and quantile, plus the immune component's fixed
settings, review criteria, term comparability, donor availability, counts,
profiles and supporting-only notes. These come from five uniquely receipt-bound
outputs: `pancreas/subpanel_audit.json` and the immune `coverage.json`,
`term_comparability.csv`, `unavailable_rules.json` and `per_donor_counts.csv`.
The recorder archives those outputs; it does not refit them or turn a
supporting-only diagnostic into an identity call.

Once the full candidate and current readiness dispositions are finalized, the
bundle command is idempotent and uses the separate remediation history:

```bash
.venv/bin/python scripts/evaluation/record_remediation_observations.py \
  --project "$PORTABILITY_REPO/outputs/blocker_remediation_20261102/parameter-history"
```

The optional project-scoped MCP server uses the existing isolated environment:

```bash
PYTHONPATH=core/multiplex5000 external/.venv-parameter-recording/bin/python \
  -m multiplex5000.parameter_observations \
  --project "$PORTABILITY_WORK/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" serve
```

It exposes record/history/export, with no approval, recommendation or application
operation. Future parameter reuse requires compatible modality, tissue, panel,
units, method version and relationship type. Unknown compatibility is not a match.

After the scientific history is populated, verify a real protein observation's
idempotent retry, read and exact export through stdio MCP using the environment
that contains the MCP client:

```bash
PYTHONPATH=core/multiplex5000 external/.venv-parameter-recording/bin/python \
  scripts/evaluation/verify_scientific_history_mcp.py \
  --project "$PORTABILITY_REPO/outputs/blocker_remediation_20261102/parameter-history" \
  --kintsugi "$PORTABILITY_REPO/external/KINTSUGI" \
  --operation protein_typing
```

This writes `scientific_mcp_check.json` and `mcp_verified_export.json` in that
project and checks that global history is unchanged. The ordinary scientific
environment need not contain the MCP client; use the interpreter shown above.

## Verification and common stops

Run focused checks appropriate to the changed route; these are computational
checks, not cohort scientific review:

```bash
.venv/bin/python -m pytest -q -o addopts= \
  core/multiplex5000/tests/current/test_protein_evidence.py \
  core/multiplex5000/tests/current/test_protein_section.py \
  core/multiplex5000/tests/current/test_pln_typing_candidate_stage.py \
  core/multiplex5000/tests/current/test_io60_rna_typing_pilot.py \
  core/multiplex5000/tests/current/test_resolved_parameters.py
.venv/bin/python -m pytest -q -o addopts= \
  core/multiplex5000/tests/current/test_maldi_objects.py \
  core/multiplex5000/tests/current/test_maldi_objects_stage.py
.venv/bin/python -m pytest -q -o addopts= \
  core/multiplex5000/tests/current/test_paired_typing.py \
  core/multiplex5000/tests/current/test_same_slide_pair_stages.py
```

| Stop | What to inspect |
|---|---|
| Receipt/hash/frame rejected | Exact changed source, accepted receipt/gate, transform/displacement, region, units and calibration. Create reviewed new evidence; do not edit the old receipt. |
| Missing marker/gene or low signal | Actual acquired panel and emitted availability. Preserve unavailable calls; zeros are measured evidence, not missingness. |
| Refinement cannot evaluate a cell | Carried-marker decisions, fit population and `refine_state/refine_note`; retain primary calls and IDs. |
| No real same-slide result | Independent source calls, cell universes, exact section/frame and accepted transform; joint output alone is insufficient. |
| MALDI remains blocked | Inspected exact neighbour/tissue and accepted physical transform, not merely available morphology masks. |
| Object result unavailable | Native area geometry, exact source frame/calibration, erosion eligibility, finite feature support and actual islet evidence for distance. |
| Changed setting reuses old output | Use a new candidate root; changed-setting checks intentionally refuse an incompatible completed receipt. |
| Recorder retryable/conflict | Storage/outbox versus invalid schema/revision/content; preserve science and retry only the record. |

After all cohort verification, inspection, documentation and history are final,
build the new typed artifact manifest and run the delivery guard again:

```bash
.venv/bin/python scripts/evaluation/verify_remediation_delivery.py --build-manifest
.venv/bin/python scripts/evaluation/verify_remediation_delivery.py
```

The default paths target this remediation bundle. `--base` explicitly selects
another checkout containing the same declared bundle; `--out` selects a separate
verification report. The guard rejects pending acquisition dispositions, missing
cohort receipts/inspections, stale evidence, incomplete pair accounting and
unapproved-record promotion. It also verifies the frozen baseline separately.
It checks receipt-bound prior scientific audits without claiming fresh reads of
all large images or recomputing scientific typing. The manifest/provenance files
are written under the new remediation documentation directory, never over the
historical baseline.

Before reporting completion, reconcile every declared acquisition and pair,
inspect every section's diagnostics, retain unavailable/excluded entries, verify
source/output hashes and export parameter history. Keep candidate computation,
geometric validity and biological applicability as separate conclusions.
