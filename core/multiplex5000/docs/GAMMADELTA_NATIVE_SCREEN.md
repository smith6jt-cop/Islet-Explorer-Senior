# Native γδ candidate screening

`mplx5k gammadelta` runs the approved raw-count screen, independent selection
validation, native-coordinate extraction and independent coordinate validation.
The maintained producer is `multiplex5000.gammadelta`; the dated November 5
scripts remain frozen experimental evidence. There is one maintained scientific
implementation and one deliberately independent implementation of the predicates
for validation.

This command produces a dataset-analysis result. It does not assign a definitive
cell identity, accept Stage 06, or activate a release. It does not require protein
typing or an existing `Immune_T` call.

## Scientific contract and authorization

The base rule is `(CD3E >= 2 OR CD3G >= 2) AND (TRDC >= 2 OR TRGC2 >= 2)`;
every threshold is per gene in the same native cell. Fifteen concurrent foreign
programs are recorded as flags. All base candidates, including flagged cells,
are retained and passed to extraction. `TRAC >= 2 AND TRBC2 >= 2` is an annotation.
Missing genes are unavailable and stop contract evaluation; they are not zeros.

The user subsequently instructed: “Proceed with the 15 criteria. Preserve all
flagged cells, keep nucleus and tissue eligibility separate, and rerun the
corrected γδ selection and coordinate extraction. Show counts and membership
changes at each step.” The receipt was recorded on November 5, 2026; this is not
an inferred timestamp of the message itself. The earlier per-gene threshold and
the exact reviewed criteria are pinned in the historical runtime contract.
`gammadelta/source_provenance.json` records source hashes and the authorization
receipt digest. The bundled criteria file retains its historical proposal status
verbatim so its content hash remains valid; the later authorization governs this
execution. The proposal text is not a new approval or a filtering instruction.

For this route the later explicit instruction supersedes older general
no-nucleus exclusion wording. Zero, one, multiple and unknown nucleus counts
remain separate annotations; no nucleus filter is applied. Tissue eligibility
remains `not_evaluated`, independently of the instrument's region label.
Scientific eligibility can be reviewed later without destroying these candidates.

## Inputs, units and command

Use an external, frozen input directory containing:

- `approved_runtime_config.json`, including the input-inventory digest;
- `native_sources.json`, with complete hashes and paths of native H5, metadata
  and transcript sources, and explicit slide/region/experiment metadata;
- `reviewed_criteria.json` with SHA-256
  `b060488f17175d346277f112683c3e48805a663d3663ad16cbc81ce9d56eb1c3`;
- `historical_and/<region-key>.csv` and `historical_or2_cells.csv`, with the
  comparison hashes declared in the runtime config.

The November 5 study contains 25 native regions on 25 slides from 13 experiments.
Its legacy `run_id` field is a **slide identifier**, unique per region only in
that particular scope. `--region` selects that inventory row; it does not count
instrument runs or acquisitions. The inventory retains the actual region UUID,
experiment UUID and slide ID. Multiple regions on a slide need distinct inventory
keys, not duplicate `run_id` values. Expansion beyond the approved scope requires
an explicit new input contract and source inventory, not filename inference.

On the workstation, use the maintained environment explicitly; it resolves the
live Senior package from a neutral working directory without a `PYTHONPATH`
override. The repository `.venv` is a different environment and is not an alias
for this interpreter. The `/path/to/...` values below are placeholders for the
chosen frozen input contract and a new output folder:

```bash
MPLX5K_PY=/data/envs/multiplex5000/bin/python
env -u PYTHONPATH "$MPLX5K_PY" -B -m multiplex5000.cli gammadelta \
  --inputs-root /path/to/frozen/inputs \
  --region 0037634 \
  --output-dir /path/to/new/0037634 \
  --dry-run
```

Dry-run reports the producer/runtime, native unit, input bindings, required source
availability/stat checks, planned commands and output paths without creating the
output or running any stage. `READY_FOR_EXECUTION` is a prerequisite report;
`BLOCKED` (exit code 2) describes missing/inaccessible sources or changed source metadata. Invalid contracts and unsafe output paths are rejected with an error (exit code 1). It does not
rehash all raw sources, scan transcripts or imply biological acceptance. Remove
`--dry-run` only when deliberately executing the inspected new result. Execution
still performs full source verification and independent technical validation.

Equivalent installed entrypoint: `mplx5k gammadelta ...`, when its environment's
`bin` directory is on PATH. This dataset-analysis route uses `--inputs-root`,
`--region` and `--output-dir`; release-stage `--project`/`--releases-root` arguments
do not select its data or output location.
The output folder must not exist. Partial results remain as evidence; use a new
folder for a retry. The command also refuses output inside an immutable run
package identified by `run_manifest.json`.

For scheduled execution, `--stage selection` runs the first two steps.
`--stage coordinates --selection-dir /path/to/validated/selector` consumes that
selection and runs the last two into a separate new output folder. Each command
records a `pipeline.json` with exact producer/input hashes, commands, runtime,
step exit status and the independent validation-receipt hashes. A failure is
recorded as `FAILED`; it never becomes technical completion by directory presence.

Native paths are site-specific. A replica must match the complete original hash
before rebinding its path, size/mtime and enclosing bundle in a new inventory.
Keep the original inventory and a dated rebinding receipt. Do not rewrite the
original scientific package to make workstation paths fit HiPerGator.

## Outputs and validation

Selection preserves exact native IDs, every T-gate cell, all raw predicate counts
and masks, base/flagged/unflagged views, every flag overlap, ordered accounting
counts, historical membership changes and the selected full native feature matrix.
The independent validator does not import the selector or its predicate evaluator.

Extraction retains every original transcript column and all QVs in the union of
closed centroid ±30 µm X/Y squares, plus every row assigned to any base candidate.
It adds only ROI, assignment and QV≥20 annotation flags. Coordinates stay in native
micrometers. The independent validator uses different geometry, scans every
source row, compares every native field and reconciles assigned QV≥20 rows to the
original H5 across all native features, including zeros and controls. Background
and neighboring transcripts in these exports are not all assigned γδ transcripts.

Technical equality establishes transfer fidelity, not biological acceptance.
The existing dated package remains the authoritative completed 25-region result;
integrating its producer does not imply a new full-cohort rerun.

## Existing result and review walkthrough

The completed November 5 result remains available on the workstation at
`/home/smith6jt/Documents/papers/output/xenium_gammadelta_corrected_20261105`.
Its `final_report/REPORT.md`, `final_report/step_counts.csv`,
`final_report/historical_membership_changes.csv` and `coordinate_index.csv`
provide the complete counts, membership changes and coordinate exports. The
frozen result contains 2,331 candidates including all 35 flagged cells; no new
analysis is needed to open or review it.

Use the [run/review walkthrough](/home/smith6jt/Documents/papers/consolidation/module_usability_implementation_20261106/gamma_docs/WALKTHROUGH.md)
for the verified local selection-only input profile, concrete commands, expected
receipts and source-availability limits. The original full inventory points to
HiPerGator native files. Its existence on the workstation does not make those
sources locally available; the region 0037634 local replay contract rebinds only
its verified H5 and metadata, not native transcripts.

The [native segmentation review](/home/smith6jt/Documents/papers/consolidation/scientific_review_20261106_segmentation/REVIEW.html)
shows all 152 focused cell cases with cell/nucleus outlines and DAPI. Choose
**Foreign-program cells** for the 35 flags; use **Nucleus cases** and
`nucleus_zero` for the seven zero-nucleus cases. Record identity/segmentation
observations and tissue eligibility separately. Use **Export review progress**
before closing and **Resume saved review** to continue; store exports outside the
frozen evidence package. Original-packet judgments are preserved in import
history, with previously reviewed γδ cases marked as drafts for the new boundary
evidence. Review notes do not alter membership or accept a release.

Installed-wheel, neutral-directory CSV/Parquet fixtures already passed in the
[dated integration packet](/home/smith6jt/Documents/papers/consolidation/gammadelta_integration_20261106/README.md).
Reuse that evidence and repeat only checks affected by new code. The subsequent HiPerGator deployment/runtime verification is complete for one
bounded native region, as described below. It remains distinct from the frozen
full-cohort result.

## Verified HiPerGator deployment (November 6, 2026)

Use the [HiPerGator run/review walkthrough](/home/smith6jt/Documents/papers/consolidation/gammadelta_hpg_deployment_20261106/HPG_WALKTHROUGH.md)
for the exact verified launcher and Slurm commands. The additive deployment is
`/blue/maigan/smith6jt/mplx5k-gammadelta-20261106-6e70d4c9`; its
`bin/mplx5k-gammadelta` invokes the preserved `reader_arrow23/bin/python` and
selects the γδ subcommand. It supplies the verified wheel and a separately
hashed pure-Python PyYAML copy without changing either existing environment.

Slurm job `50019341` passed all four steps for native region `0059865` in
5 minutes 47 seconds. All 80 comparison checks passed: 26 candidates including
the flagged two-nucleus candidate and all 480,284 exported transcript rows
matched the historical result. The preservation check found no changes in its
10,776 protected entries; its receipt states the exact metadata/full-hash scope.
A separate job verified the reusable `scripts/run_region.sbatch` dry-run with a
new output location and no analysis output creation.

This verifies the tested region and native Parquet format, not a new execution
of all 25 regions. The original 2,331-candidate result, scientific criteria,
review decisions and active release pointers are preserved. Submit computation
through Slurm using the linked walkthrough. The remote
`deployment_verified.json` records the package, runtime and bounded result;
`VERIFIED_COMMANDS.md` preserves the verified operational instructions.

## Runtime and retention dependencies

The historical HiPerGator production run used Python 3.13.9, PyArrow 23.0.0,
h5py 3.14.0, NumPy 2.2.6, pandas 2.3.3 and SciPy 1.16.2. Its interpreter remains
`/blue/maigan/smith6jt/xenium_transcript_coordinates_20261020/reader_arrow23/bin/python`.
Preserve that environment and historical directory until a replacement runtime
has been installed and verified against the native source format. A workstation
synthetic or selection-only pass does not replace that native-transcript test.

The package declares these scientific dependencies already. No package versions,
raw files, old scripts, historical outputs or release pointers are changed by
this integration. The separate dated integration packet records the tested
runtime versions, source migration, full CLI fixtures and native replay scope.

Contract/preflight rejection occurs before output creation; it does not create a `FAILED` result directory.
