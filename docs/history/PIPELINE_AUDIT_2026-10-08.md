# PhenoCycler and Xenium pipeline audit

Audit date: 2026-10-08. The current workflow contains recoverable numerical components and substantial later work, but its end-to-end execution and provenance contracts are unreliable. Confirmed defects affect stage completion, forced recomputation, package resolution, acquisition identity, release completeness, coordinate consistency and notebook execution. These observations establish software behavior; they do not establish which assistant or author introduced a defect or whether every existing biological result is wrong.

Read this together with the [unified workflow plan](UNIFIED_WORKFLOW_PLAN_2026-10-08.md) and [repository evaluation](REPOSITORY_EVALUATION_2026-10-08.md). At the initial audit checkpoint, processing code, source data, environments, annotations and existing results were preserved; the deliverables were evaluation documents and an inventory, with no pipeline repair or cohort rebuild.

Later work on October 8 completed corrected protein/RNA typing generations and verified downstream consumers in the sibling v3 workspace. The [current typing review](/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3/results/v3/islet_immune_reassessment_20261006/cell_typing_validation/current_cohesive_typing_review_20261008/CURRENT_CELL_TYPING_REVIEW.md) binds those selected generations, independent processing checks and explicit scientific limits. The findings below describe the initial baseline; the unified CLI/HPC/notebook migration remains proposed, and independently measured biological accuracy remains unestablished.

## Scope and source generations

The connected account inventory contains 93 `smith6jt-cop` repositories; the next pagination page was empty. Metadata, recursive source trees and introductory documentation were screened for 79 imaging, omics, HPC or adjacent candidates. Fourteen unrelated repositories were inventoried only. The repository matrix distinguishes architectural screening from deeper code inspection and runtime checks. It is not an exhaustive survey of all repositories on public GitHub or a runtime benchmark of all 79 candidates.

Senior's PhenoCycler, Xenium, integration, restain, release, Slurm and notebook paths received executable-code review and synthetic checks. The local v3 checkout and selected current artifact manifests received recovery review. KINTSUGI received selected CLI, workflow and stage-script inspection; MultiSpaceNet, MaxFuse, UniFORM and MCMICRO received selected implementation inspection. Official 10x, scverse and Snakemake documentation informed the proposed interfaces.

| Source | Inspected identity | Relevant condition |
|---|---|---|
| Current Senior checkout | `main`, `ac4b47539514eef64c99ef6c3e7bbf66e44c8fd0`, September 26 | Initially one untracked user notes file; no tracked modifications |
| Published Senior main | `e0374fe9f86f0ee261de60cf3b66c793e23931f6`, September 29 | Seven commits ahead of this local main; includes v3 Parquet/context ingest changes |
| Sibling Senior v3 | `v3/typing-inputs`, `2714261a9e21cf2d575f97968bfd0e628089c152`, October 2 | 26 modified tracked files and 61 untracked entries at inspection, including a substantial results tree |
| Preserved Phenocycler repository | `main`, `e9874a9ffa6384857fc22179db0984bdf570c4d4`, August 14 | Separate recovery source; three deleted integration CSVs in its working tree |
| Published Phenocycler and Xenium sources | `fc06c0f7a394` and `8062e55407e5` | Already absorbed sources, useful for history/recovery rather than independent future producers |
| PostXeniumIHC | Non-Git QuPath project | 20 image entries/qptiffs/qpdata files; geometry/tissue exports, no standalone processing scripts |
| Current KINTSUGI main | `8a48ab04b9eb5fbe065e55f93303159f409261f6` | Newer CLI/Snakemake/Slurm work makes the older Senior description incomplete |
| IO60panc2nd-v3 GitHub scripts | `bcacffaa2c9e` | Current project-level QuPath processing/export evidence |

The sibling's newer timestamp does not make it a validated replacement baseline. Recover its committed changes, dirty changes and output provenance separately. The existing `docs/EXTERNAL_PACKAGES.md`, consolidation plans and `CLAUDE.md` are useful historical evidence, but multiple dated corrections contradict earlier sections; executable behavior and exact artifact manifests take precedence when describing current state.

## Confirmed lifecycle and environment defects

| Priority | Finding and evidence | Effect |
|---|---|---|
| P0 | `core/phenocycler/phenocycler/pipeline.py:547` skips a stage if its output root has any child; `donor_pipeline.py:275` uses the same presence rule for donor partitions | A partial directory or lone crashed `.tmp` file counts as completion; missing/stale partitions can be skipped |
| P0 | Full-stage `--force` reaches existing outputs, but `expression.py:280` refuses overwrite; `cells_parquet.py:438` refuses nonempty ingest | Advertised forced recomputation fails against real writers; a temporary-directory reproduction raised `FileExistsError` |
| P0 | `donor_pipeline.py:556` and `:596` do not pass forced execution to the donor worker | A forced donor queue enqueues the task and then skips it as already present; reproduced with temporary outputs |
| P0 | `make env-check` imports `phenocycler` from `Islet-Explorer-Senior-v3/core/phenocycler` | Scripts from this checkout can execute sibling code. Shared editable installation plus v3's symlinked `.venv` defeats checkout identity |
| P1 | `_aggregate_donor_audit()` at `donor_pipeline.py:475` has no production caller | Pipelined execution omits the promised cohort expression availability audit, although direct-helper tests pass |
| P1 | `pipeline.py:319` loads cohort sources without fingerprint validation; portions of `artifacts.py:878`/`:950` retain descriptive identity fields without validating current sources | Declared digests do not bind current stage results to actual input content |
| P1 | Both PhenoCycler notebooks import deleted `StageManifest`; notebook 01 uses removed `RunContext` members and notebook 02 hardcodes an old run ID | Review notebooks cannot execute against the current API; persisted notebook outputs are historical evidence |

The stage framework deletion is traceable in Git to `7f868d8`, followed by `776fb83` addressing force handling. This makes recovery possible without attributing cause to a particular agent. Restore selective dependency fingerprints and real lifecycle tests, rather than restoring a global hash that invalidates every stage after unrelated source changes.

## Existing donor, slide and ROI records

The project already records donor, acquisition and ROI identities. The repo-local [IO60panc2nd Notes.md](../IO60panc2nd%20Notes.md) contains 26 donor entries with pancreatic anatomy, pLN presence, relative block origin, assay relationships and cautions such as the 6457 reprocessing requirement. Its opening donor record matches the parent project's `ProcessingNotes.txt`. The metadata gap is inconsistent consumption and reconciliation of these records, rather than an absence of donor/slide information.

| Existing source | Verified scope | Role and limits |
|---|---|---|
| `IO60panc2nd Notes.md:2–27` | 26 donors | Lab descriptions of tissue/block relationships and sequential/same-slide assays; preserve the stated ambiguity around lymphoid aggregates/TLS |
| `core/phenocycler/data/integration/xenium_paths.csv` | 58 rows | Tracked donor/ROI/batch/bundle mapping, including additional 6533 and 112 acquisitions; an unavailable path is explicitly `NA` |
| `core/phenocycler/data/integration/manifest.csv` | 57 rows, 26 known donors, 20 populated PhenoCycler image mappings and two donor-unknown rows | Generated, ignored acquisition snapshot; all `block_id` and `section_gap_um` values are blank, so those columns do not carry the origin already described elsewhere |
| `section_modes.csv`, `section_overrides.csv` and `donor_overrides.csv` in the same directory | Six mode declarations and explicit image/donor exceptions | Acquisition-specific relationship/alias evidence; current Senior modes include serial and same-slide 6533/112, while v3 adds 6347/6580 |
| `core/phenocycler/cohort-spec.json` and `data/qupath_manifest.json` | 23 donor/image entries each | Explicit PhenoCycler image/source identities; cohort spec is tracked and the QuPath manifest is generated/ignored. The local v3 cohort spec and generated manifest have 25 entries |
| `pln_node_annotations.csv`, `restain_registration_inputs.csv` and `restain_sections.csv` | 19 dedicated-pLN records, 20 IF registration mappings and five IF channel-role records | Exact node/image/donor/ROI or IF image/serial identities. Dedicated-pLN records must remain distinct from pLN annotations inside pancreas images |
| [Parent-project `results/v3/pln_blocks.csv`](/home/smith6jt/IO60panc2nd/results/v3/pln_blocks.csv) | 12 annotations across 10 donors: seven same-block and five different-block | Image/entry identity, annotation UUID, piece count, measured gap/overlap and separately recorded `block_from_notes`; used by v3's block-aware alignment/typing code |

The block ledger already distinguishes the two pLN annotations in each of 6347 and 6580: one same-block and one different-block. It records 6533, HDL112 and HDL115 pLN annotations as different-block, and the documented 6539/6551 pLN annotations as same-block. Preserve recorded origin separately from measured geometry; a gap or contact cannot establish block provenance. Likewise, notes describing newer same-slide assays must not replace the older serial acquisitions retained in the mapping tables. Explicitly reconcile `112`, `0112` and `HDL112` aliases against their source images and acquisitions.

Recovery must preserve the untracked notes, ignored generated manifests, untracked v3 cohort spec and external annotation/block ledger alongside tracked tables. A Git-only snapshot would omit part of the identity evidence.

## Confirmed acquisition and relationship hazards

`core/phenocycler/data/integration/section_modes.csv:19` onward records sequential and same-slide pancreas acquisitions for the same donors, distinguished by Xenium serial. `scripts/senior/same_section_linkage.py:103` collapses these to a donor-level dictionary, then reads donor/ROI cells and a donor-level warp. Its duplicate-call handling at `:306` also drops duplicates by native cell ID without acquisition identity. A same-slide acquisition can therefore lend its relationship to another section from that donor.

`analysis/xenium/scripts/panc_xenium_release.py:299` writes each section to `calls/donor_id=<donor>/roi=panc/data.parquet`. Two pancreas sections for one donor would overwrite that path. This is a confirmed code-path hazard, not evidence that current one-section-per-donor pancreas outputs were overwritten. Cross-section ID collisions are real in the October pLN manifest: `duplicate_cell_ids_within_donor` records 8 for 6534 and 36 for 6579. Donor plus native cell ID cannot be the general observation key.

`scripts/senior/build_sections_layer.py:170` infers `same_block` from pLN annotations in a pancreas image and propagates it by donor at `:198`; the same defect exists in v3. That consumer and its test fixtures lack explicit block IDs, but the lab notes and annotation ledger already record relative block origin. The builder ignores those distinctions, including mixed-origin annotations in 6347/6580. Co-imaging alone cannot establish block origin. Use the recorded annotation-specific origin and acquisition mapping. In this session the user also confirmed the intended attached pancreas/lymphatic/pLN case as the **same block and physical section**; preserve its interface geometry without extending that relationship to every tissue from the donor.

The v3 `build_current_input_manifest.py:306` substitutes a global sequential-mode fallback for unoverridden inputs. Its October 3 roster has 54 of 58 RNA rows with that fallback. Those rows need evidence-backed relationships or an unknown state. `integration/pipeline.py:175` also selects one global stage order despite per-section modes; it cannot serve as the mixed-design execution graph.

## Release publication and availability defects

The pancreas release driver reuses a section whenever `<section>.report.json` exists (`panc_xenium_release.py:224`), independent of input, code, panel and scoring parameters. Worker failures are logged at `:235`; assembly skips missing reports at `:244` and continues. A `--sections` subset can write to the normal cohort destination, and existing partitions survive smaller rebuilds. Assembly hashes current h5ad files after loading cached reports (`:249`), so the published digest need not identify the input used to compute a cached call.

`scripts/senior/modality_concordance.py:399`, `:468` and `:606` build lineage count rows from nonzero `value_counts`. `shares()` at `:623` leaves omitted lineage entries as NaN. A donor with a measured zero B/plasma count therefore gets an unavailable share and can disappear from summaries. An in-memory two-donor reproduction confirmed this; the current pancreas release also explicitly reports a zero B/plasma count for donor 115. The fix must distinguish an assayed zero from an unsupported lineage instead of blanket zero filling.

`scripts/senior/xenium_release.py:159` classifies `np.nan` and `pd.NA` as `called`. An unmatched label join can therefore become a non-Endocrine negative in agreement analysis (`same_section_linkage.py:307`). Null labels, incomplete joins and absent identity releases need explicit unavailable states.

## Coordinate and cell measurement consistency

`scripts/senior/build_restain_layer.py:199` reads `warped_if_um.npy`, which `core/phenocycler/scripts/restain/valis_register.py:139` created with the VALIS seed transform. IF quantification in `core/phenocycler/scripts/restain/protein_rna.py:740` uses the refined transform. The viewer and the measurement thus use different coordinates; the same-section linkage source itself calls the earlier warp superseded.

A read-only comparison used 10,000 regularly spaced centroid indices per donor, comparing refined-transformed native coordinates with the stored seed-warp array. These numbers are coordinate differences, not errors against ground truth:

| Donor | Median difference in µm | p90 in µm | Maximum in sampled points in µm |
|---|---:|---:|---:|
| 6442 | 1.783 | 3.349 | 19.580 |
| 6436 | 1.639 | 2.579 | 26.348 |

The later v3 `analysis/xenium/scripts/build_registered_postxenium_pairs.py:220` loads legacy intensity vectors and newly loaded Xenium IDs, checking vector length but lacking an ordered-ID hash. The producer at `core/phenocycler/scripts/restain/protein_rna.py:757` records cell count/bundle rather than the original ordered IDs. A same-sized reordering is undetectable. Prefer keyed Parquet measurements and exact transform/segmentation identities; regenerate measurements when historical row order cannot be established.

## HPC and notebook workflow defects

| Evidence | Consequence |
|---|---|
| `core/phenocycler/scripts/slurm/run_full_pipeline.sh:24` invokes `run --jobs` without `--pipelined` | The HPC entry uses whole-stage presence skipping and assumes the invoking Python resolves the correct package |
| `analysis/xenium/hpc/run_pipeline.sh:48` sources absent `scripts/env_rapids.sh`; the file is under `hpc/` | Relocated HPC driver fails before processing |
| `analysis/xenium/hpc/slurm/run_all_samples.sh:14` submits absent `scripts/slurm/*` | Submitted paths disagree with the checked-in `hpc/slurm/` layout |
| Xenium Slurm phenotype job runs only notebook 02; local driver additionally runs finalize, hierarchical typing and validation | Slurm and local execution can produce different analytical stages |
| Xenium driver/notebooks use different data roots and historical sample lists | One `XENIUM_DATA_ROOT` or sample selection does not configure the entire chain |
| `scripts/senior/slurm/run_step.sbatch:12` reserves B200/32 CPUs/128 GB for every step and changes to a fixed checkout at `:23` | CPU jobs waste GPU allocations and submitted source identity is unpinned |
| Slurm output files reference `logs/` before the job creates it; environment activation falls back silently at `:28` | Missing submit-time directories or differing runtimes can break or alter processing |
| `Makefile:29` embeds local runtime paths; `build_layers.py:18` uses mtime freshness | Portability and reproducible dependency tracking remain incomplete |

The two protein notebooks and seven Xenium notebooks do not form the requested review series. For example, immune proximity hardcodes samples 0041323/0041326; phenotype QC hardcodes count/gene floors; group comparison notebooks use historical individual-sample cell-level tests. Preserve useful plots and helper logic, but move processing to shared modules and use donor/section-aware inferential analysis.

`core/phenocycler/docs/HIPERGATOR.md` still promises manifests, receipts and recovery capabilities removed from executable code. Treat that as stale documentation. Runtime/environment differences explain some test failures and must be separated from analytical defects.

## Later artifacts and recoverable progress

The resolver currently selects the October 3 pancreas release at `/data/xenium_typing/panc_20261003_bundle/release` with 25 declared/observed donors, and the pLN release at `/data/xenium_typing/pln_20261003` with 24 donors across 27 sections. Existing pancreas concordance files under `data/development/integration/concordance/` contain only two donors. Their inventory is stale relative to the currently selected inputs.

Three of the eight recorded pancreas producer-script hashes differ from the current older Senior checkout. All eight match the inspected v3 sources. This establishes a recoverable producer code generation; it does not validate the labels, the full runtime environment, or every cached input. The release declares development status. The pLN manifest records a dirty checkout without the same full script-hash ledger.

The October 6 v3 `cell_typing_checkpoint_manifest.json` declares `NOT_PASSED` and `accuracy_validated: false`, naming protein Cell/Nucleus source reprocessing/parity and pancreas RNA B/plasma adjudication as active blockers. These are declarations in a recovered artifact, not independent biological conclusions from this audit. They preclude treating its downstream figures as validation evidence by themselves.

Retain these useful components:

- QuPath UUID/image identity and exact cell-universe checks, acquired/missing marker distinctions, geometry-first eligibility and compartment-resolved expression.
- Sparse CPU REDSEA, per-channel image access and process-pool/cancellation primitives, with measured GPU limitations retained.
- The restain registered-export adapter: exact donor/serial/ROI, transform path/hash, PASS alignment and withheld-fold validation, followed by transactional export and decoded-image QC.
- The Xenium whole-input preflight, declared panel taxonomy, deterministic section-level naming, pLN partial-scope guard and release availability vocabulary.
- v3 cells/context Parquet pairing checks, pLN inputs, explicit block assignment and reference-donor preservation for UniFORM.
- The shared R domain/viewer layers and existing layer dependency graph, after replacing mtime/input discovery with pinned release artifacts.

The preserved earlier repository supplies content/manifests/receipts and method history, but its phenotype exporter explicitly refuses to proceed without verified identity. Its notebooks describe that earlier API and are recovery material, not a currently working cohesive pipeline.

## External component findings

Current KINTSUGI supports a Click CLI, a Snakemake workflow and Slurm profiles. Its image stages and resource routing are useful foundations. Selected inspection also found contracts to fix before adopting it: `src/kintsugi/cli.py:179` imports registration code from notebook paths; `workflow/scripts/registration.py:143` considers existing registered filenames sufficient for reuse; `workflow/scripts/signal_isolation.py:104` ignores error-status markers in its existing-output check and `:213` publishes success/sentinel metadata after counting failed channels. These are static code observations; KINTSUGI was not installed or executed in this audit.

The inspected MultiSpaceNet loader (`src/multispacenet/data/loader.py:152`) intersects RNA and protein observation names before constructing paired features. Its paired model belongs only downstream of validated same-section correspondences. The custom MaxFuse implementation and Senior's MaxFuse-class linking represent inferred similarity, not measured cell identity. Do not use either to bypass experimental-design or registration evidence.

Selected upstream comparisons were checked through GitHub: VALIS was identical; TissUUmaps was 16 commits ahead; ASHLAR, InstanSeg, MCMICRO, Scanpy, Squidpy, SCIMAP, RAPIDS SingleCell and Nextflow forks were behind their upstream default branches. Full values and inspected commits are in the matrix/inventory. These comparisons establish source relationships, not upgrade recommendations or runtime correctness.

## Verification results and practical limits

| Check | Result | Interpretation |
|---|---|---|
| Root `make env-check` | Failed | Current checkout imports sibling v3 package; console scripts resolve the same interpreter |
| Explicitly current protein lifecycle/geometry/expression/relationship subset | 46 passed | Does not cover realistic writer/orchestration failures; separate temporary reproductions exposed them |
| Full `core/phenocycler/tests/current` default selection | 1001 passed, 2 failed, 114 deselected, 2 xfailed | Both failures are neutral-directory wrong-checkout import guards |
| Restain synthetic suite | 224 passed, 3 deselected | Strong reusable adapter tests; real registration accuracy not established by these tests |
| Xenium suite with addopts override | 117 passed, 1 failed | Failure is a pandas string-dtype expectation in a denominator test, not demonstrated wrong nucleus handling |
| Release/restain/crosswalk subset | 80 passed, 4 deselected | Existing local component behavior; acquisition/release lifecycle gaps remain |
| Selected v3 recovery tests | 55 passed, 3 data tests deselected | Supports recovery of ingest/block/reference plumbing, not biological accuracy |
| Temporary/in-memory reproductions | Confirmed | Stray-file completion, force/writer conflict, force-queue skip, deleted notebook import, zero-share and null-label defects |

Suites overlap; these counts should not be added as independent coverage. Checks disabled pytest cache/bytecode where invoked and did not alter production data. Exact main commands are recorded below. The Xenium and root consumer suites used the shared environment, whose import contamination is itself a finding; explicitly-current protein checks and the package-level suite's pythonpath identify the protein code being tested.

```bash
make env-check

PYTHONPATH="$PWD/core/phenocycler" .venv/bin/python -B -m pytest -q \
  -p no:cacheprovider \
  core/phenocycler/tests/current/test_pipeline.py \
  core/phenocycler/tests/current/test_donor_pipeline.py \
  core/phenocycler/tests/current/test_geometry_stage.py \
  core/phenocycler/tests/current/test_expression.py \
  core/phenocycler/tests/current/test_section_relationship.py

# From core/phenocycler
../../.venv/bin/python -B -m pytest -q -p no:cacheprovider -rfEs tests/current
../../.venv/bin/python -B -m pytest -q -p no:cacheprovider -rfEs scripts/restain/tests

# From the repository root
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider -q \
  analysis/xenium/tests --override-ini='addopts='
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider -q \
  core/phenocycler/tests/current/test_xenium_release.py \
  core/phenocycler/tests/current/test_restain_layer.py \
  core/phenocycler/tests/current/test_crosswalk_reader.py
```

No raw image/RNA cohort processing, HPC submission, deployment, source repair or independent biological accuracy validation occurred. The implementation plan begins with source recovery, isolated environments, identity and lifecycle repairs, then a representative small-cohort workflow before full-scale analysis.
