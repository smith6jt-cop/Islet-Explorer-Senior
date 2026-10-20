# Dispositions — every branch, PR, worktree and repository (2026-10-20)

Phase B of the MultiPlex5000 plan (`~/.claude/plans/evaluate-the-local-implementation-wondrous-seal.md`)
requires a recorded disposition for each of these before code moves. **Default for anything the
maintainer has not answered: preserved (tagged or bundled) and inactive** — the integrator makes no
scientific choice. Preservation already done: 10 verified bundles and 2,318 uncommitted source
files in checksummed tarballs under `/data/io60-archive/code_bundles_20261019/`.

Status words: **merge** (lands on main through a PR with green CI), **merged**, **record** (closed
or kept unmerged, preserved by an `archive/…` tag), **tag + delete** (an
`archive/branch/<name>-20261019` tag, then the branch is removed — phase C), **remove** (a worktree
or clone, after checking it holds no ignored data), **archive** (a GitHub repository made read-only
with a pointer README — phase G), **owned** (another session's; not touched here).

## Done on 2026-10-20 (phase C, the parts that do not wait on a merge)

- **Monorepo branches.** 22 branches tagged `archive/branch/<name>-20261019` (pushed): the 19
  merged ones listed below as tag + delete, plus `v3/ingest`, `v3/ci-repair` and `v3/uniform-pin`.
  Deleted from GitHub: `fix/serial-registration-tissue-outline`, `stage0/groundwork`, `v3/ingest`,
  `v3/ci-repair`, `v3/uniform-pin`, and the closed-PR branches `v3/uniform-audit`,
  `v3/sameslide-islet-metrics`, `v3/qc`, `v3/galectin`, `v3/raw-expression` (each kept by its
  `archive/pr-*-20261020` tag). Deleted locally: the 19 merged branches, `v3/uniform-pin`,
  `v3/raw-expression`. Every deletion was checked first: an ancestor of `main`, or covered by a tag.
- **Worktrees removed** (each clean, holding only caches): `.worktrees/are-you-able-to-see-the-plan-and`,
  `/scratch/io60/work/senior-main`, `senior-pr17`, `senior-raw`.
- **Remote `hpg` removed** (its only ref, 990c5b2, is recorded in the A1 bundles); `origin` is the
  one remote.
- **IO60panc2nd-v3:** PRs #2 and #3 closed and tagged (`archive/pr-2-driver-20261020`,
  `archive/pr-3-islet-method-20261020`) after re-checking that every shared file was changed later on
  the plan branch and that #3's one unique file is the detector e3a43f2 withdrew; `v3/plan-2026-10-01`
  pushed and opened as **PR #4** into main (a fast-forward, 50 ahead, 0 behind); its worktrees
  `/scratch/io60/work/wt-{benchmark,driver,islet}` removed, after confirming `wt-islet`'s five
  untracked step-1c files are byte-identical on `preserve/wt-islet-20261019`.

**Merged later on 2026-10-20** (at the maintainer's request): #25 (45b02c8), #27 (3bedc11), #28
(9bd008c), #26 (be97b1d, after all 17 checks passed) and IO60panc2nd-v3 #4 (6b1cf6c, a merge commit,
so 340bf0c keeps its SHA). Then tagged `archive/branch/<name>-20261019` and removed: the
`consolidate/b2-v3-salvage`, `b4-red-checks` and `b5-dispositions` branches and worktrees, and the
GitHub `analysis/cohesive-finish` branch — its local branch and worktree stay, because the Codex
papers session reads that worktree. IO60panc2nd-v3's `v3/integration`, `v3/driver`,
`v3/islet-method` and `v3/benchmark-tools` were removed (each covered by an `archive/` tag).
`v3/typing-inputs` stays: it is checked out in the v3 worktree, which is removed last.

Then **#29** (82ad1a6, D1 without moving files: site-set data roots, the new environment, the env
check, the profile schema) and **#30** (4f802d4, **B3**: the HPG session's same-slide lane, 37f9311,
plus one test-only commit, 11654d3, naming the lane's two sha256 helpers in the checksum guard;
the maintainer accepted its CLAUDE.md edits as they are). Both branches were tagged
`archive/branch/<name>-20261019` and removed. **Phase B is complete.** The D1 rename
(`core/phenocycler` → `core/multiplex5000`) follows in its own PR.

**Scope change, maintainer 2026-10-20: Panc Floyd is out of scope** — no push request, no frozen
CODEX release, no archiving of its copies. Its repo, ported scripts, data and session are left
alone, and its data never appears on a deletion list. The A1 bundle of its repo stays as a passive
copy.

## Gate B answers (maintainer, 2026-10-20)

**Merge none of the four draft PRs as written; tag and close them, and carry the useful pieces into
the MultiPlex5000 stages**, so that no standalone v3 tool and not the legacy app builder becomes
another maintained route. Closed with comments and tagged the same day:

| PR | tag | what carries forward, and where |
|---|---|---|
| #18 UniFORM audit | `archive/pr-18-uniform-audit-20261020` (0fa93b7) | its batch diagnostics → the normalization gate (stage 05). The pooled threshold is diagnostic only; the disease-label test does not run before processing approval |
| #19 same-slide reader + islet metrics | `archive/pr-19-sameslide-islet-metrics-20261020` (a37eacc) | the Parquet reader is ported **with the HiPerGator same-slide lane after input parity is checked** (that lane imports `read_pheno_tsv` in `prepare_protein.py` and `protein_redsea.py`); the islet metrics are reviewed against the final typing labels before adoption |
| #20 QC flags + neighbour features | `archive/pr-20-qc-20261020` (a6ee8e9) | suitable reason-coded flags and the neighbour/bleed-through features → measurement QC (stage 03); **the fixed top-percentile flags (DAPI, E-cadherin) are omitted** — they mark a fixed share of every donor by construction |
| #21 galectin from typed cells | `archive/pr-21-galectin-20261020` (4e9b935) | the typed-cell approach (the 2026-09-28 ruling) → release export (stage 11). As written it rewires the old `build_layers.py` default and reads `cell_class`; by the 2026-10-19 ruling releases serve the final step of the one hierarchical typing method |
| #24 raw means / REDSEA off | `archive/pr-24-raw-expression-20261020` (a84eeeb) | comparison evidence only (approved plan): the October 14 own-compartment chain is the canonical protein path |
| IO60panc2nd-v3 #1 benchmark scorer | `archive/pr-1-benchmark-tools-20261020` (adfc495) | may serve the segmentation gate (stage 02) once its object-size convention is compared with the selected QuPath rule; its batch-confound inference waits for processing approval |

**Preserved work**
- **112 / 6533 pLN relationship rows** (preserve branch): ported **only after the HiPerGator
  session verifies the acquisition IDs**. Main already records these donors' pancreas same-slide
  sections; the rows add the `pln_1` regions while keeping the safe `sequential` default for an
  unspecified section (`core/phenocycler/data/integration/section_modes.csv`).
- **`work/H-restain`** (217d29a, tag via `preserve/work-H-restain-20261019`): not merged. Its two
  findings become **new registration-gate checks** (stage 07), with the old commit as evidence:
  declared vs producing pixel size (6436 and 6442: 0.2505 vs 0.25 µm/px) and registration failure
  (115: 7.198 µm).
- **Restain rework** (preserve branch, 8 files): stays tagged as a record.
- **SpleenFollicleCounterQP**: untouched, as requested — no port and no archive for now.

**Xenium tables: no edit to the older lane now.**
- `CTRL` is already an acinar **identity** marker in the canonical panel roles
  (`core/panels/data/panel_gene_roles.csv`: pancreas, `02_pancreas_exocrine`, `identity`,
  Exocrine, acinar) — checked 2026-10-20.
- **F11, F13A1, PROS1** are recorded for the **tissue-panel review at Gate D**. In particular the
  "hepatocyte-specific" F13A1 exclusion conflicts with the panel metadata
  (`core/panels/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv` lists macrophage, elicited and
  alternatively activated macrophage, microglia, Hofbauer cell and others — checked 2026-10-20);
  that deserves an evidence-backed correction **in the panel authority**, not a change to an unused
  older lane.

## Carried into MultiPlex5000 stages

| stage | piece | condition before adoption |
|---|---|---|
| 02 segmentation / QuPath export | IO60panc2nd-v3 #1 islet benchmark scorer | object-size convention compared with the selected QuPath rule |
| 03 protein measurement | #20 reason-coded QC flags and neighbour/bleed-through features | fixed top-percentile flags omitted |
| 05 normalization | #18 batch diagnostics | pooled threshold diagnostic only; disease-label test not before processing approval |
| 07 registration | H-restain pixel-size and registration-failure checks | re-tested on the current restain outputs |
| 08 integration (same-slide) | #19 Parquet reader; 112 / 6533 `pln_1` rows | input parity checked with the HPG lane; acquisition IDs verified by the HPG session |
| 10 analysis | #19 islet immune metrics | reviewed against the final typing labels; analysis stays disabled until processing approval |
| 11 release export | #21 typed-cell galectin | built from the final hierarchical typing step |
| Gate D panel review | F11, F13A1, PROS1 exclusions | evidence-backed correction in `core/panels` |

## Monorepo branches (Islet-Explorer-Senior)

| branch | tip | vs main | PR | disposition |
|---|---|---|---|---|
| `main` | 4f802d4 | — | — | integration target |
| `analysis/cohesive-finish` | 026315e | merged | #25 | **merged** (B1); GitHub branch removed, local branch + worktree kept for the papers session |
| `consolidate/b2-v3-salvage` | 25da443 | merged | #26 | **merged** (B2); tagged and removed |
| `consolidate/b4-red-checks` | 2f568d6 | merged | #27 | **merged** (B4); tagged and removed |
| `consolidate/b5-dispositions` | 2dc9cef | merged | #28 | **merged** (B5, this ledger); tagged and removed |
| `consolidate/d1-site-env-profiles` | 9e03a09 | merged | #29 | **merged** (D1, no file moves); tagged and removed |
| `consolidate/d1-rename-multiplex5000` | — | +1 | its own PR | **merge** (D1, the rename to `core/multiplex5000`) |
| `v3/sameslide-four` | 11654d3 | merged | #30 | **merged** (B3, 4f802d4): the HPG session committed and pushed its same-slide work (37f9311), and one test-only commit (11654d3) named the lane's two sha256 helpers in the checksum guard; tagged and removed |
| `v3/typing-inputs` (local only) | 2714261 | merged | — | its 6 commits are on main via #25; the branch stays checked out in the v3 worktree until that worktree is removed (last) |
| `v3/uniform-audit`, `v3/sameslide-islet-metrics`, `v3/qc`, `v3/galectin`, `v3/raw-expression` | — | — | #18–#21, #24 closed | **record** (tags above); branches tag + delete in phase C; remove `/scratch/io60/work/senior-raw` |
| `v3/uniform-pin` | 60a4166 / GitHub bed51d0 | +2 | #22 closed | tag + delete (closed by ruling: the pin stays 7c33396) |
| `preserve/worktree-2026-10-14` (GitHub) | 9eb1b39 | +1, 91 behind | — | **record**, never merged whole; the pLN rows are ported at stage 08 after HPG verification |
| `work/H-restain` = `preserve/work-H-restain-20261019` | 217d29a | +1, 133 behind | — | **record**; its findings become stage-07 checks |
| `fix/serial-registration-tissue-outline`, `v3/ingest`, `v3/ci-repair` | — | merged | #16, #17, #23 | tag + delete |
| `merge/hpg-sameslide`, `stage0/groundwork`, `agents/are-you-able-to-see-the-plan-and` | — | merged | — | tag + delete |
| `work/{A,D,E,F,G,I,J,K,L,M,N,O,P,Q,S}-*` (15) | — | merged | — | tag + delete |
| remote `hpg` → `/data/hpg_peek/*` | 990c5b2 | merged | — | remove the remote; the peek clones are bundled (A1) |

## Monorepo worktrees

| worktree | disposition |
|---|---|
| main checkout | integrator; stays (`/srv/shiny-server/islet-explorer-senior` points here) |
| `Islet-Explorer-Senior-v3` | read-only; holds 105 GB of data and `results/` (65.6 GB); removed **last**, after its data moves (phase H) |
| `.worktrees/cohesive-finish`, `b2-v3-salvage`, `b4-red-checks` | remove after their PRs merge |
| `.worktrees/are-you-able-to-see-the-plan-and`, `/scratch/io60/work/senior-main`, `senior-pr17` | remove (merged; check for ignored data first) |
| `/scratch/io60/work/senior-raw` | remove (#24 closed) |

## Other repositories

| repository | disposition |
|---|---|
| **IO60panc2nd-v3** (QuPath workflow; `~/IO60panc2nd/scripts`; public) | local `v3/plan-2026-10-01` preserved as `preserve/plan-2026-10-01-20261019`; merge it into that repo's main (C); close PR #2 (29 of its 36 files are identical on the plan branch, 7 are older) and PR #3 (it carries the withdrawn detector `1_panc_objects_v3.groovy`, removed by e3a43f2); PR #1 closed as a record (above); tag + delete `v3/{integration,driver,islet-method,benchmark-tools}`; remove `/scratch/io60/work/wt-{benchmark,driver,islet}` (step-1c files preserved as `preserve/wt-islet-20261019`); imported into `core/multiplex5000` at D9, then **archive** |
| **Phenocycler_Analysis** (absorbed as `core/phenocycler`, renamed `core/multiplex5000` in D1) | **archive** (G); `Phenocycler_Analysis.preserved` and `.git/modules/Phenocycler_Analysis` deleted at G with approval (bundled) |
| **XeniumPanelExplorer** (absorbed as `core/panels`) | **archive** (G); its dependabot branch bumps workflows deliberately left out — delete |
| **Xenium_Analysis** (absorbed as `analysis/xenium`) | HPG clone: unpushed 5d2d4ee plus uncommitted edits → preserved by the HPG session (A4); no older-lane edit (Gate B); then **archive** |
| **SpleenFollicleCounterQP** (absorbed as `analysis/follicle`) | **untouched for now** (maintainer, 2026-10-20): its three branches and its HPG clone stay as they are |
| **LymphVizkit** (absorbed as `ieCore` + `apps/lymph`) | its only other branch is one empty "Initial plan" commit with no files — nothing to port; **archive** |
| **KINTSUGI** | separate package, kept and pinned once (D1); `external/KINTSUGI` 716d1eb is 50 behind; HPG clone's 2 stashes preserved by the HPG session (A4); `/data/KINTSUGI__` (109 dirty files, 2025-03) and `/data/KINTSUGI_alt` (3 dirty) are not in any bundle — preserve before any deletion; none is scheduled |
| **Islet-Explorer (Panc Floyd)** | **out of scope** (maintainer, 2026-10-20): its repo, ported copies, data and session are left alone; nothing in this plan reads, moves, freezes, archives or deletes them |
| HPG `Islet-Explorer-Senior-main` | owned by the HPG session; its same-slide work was pushed (37f9311) and merged as #30 (B3). The session switches the clone to `main` before its next commit and, once the D1 rename merges, reinstalls the package from `core/multiplex5000` |
| HPG `Xenium_Celltype/code`, `~/IO60panc2nd/scripts`, `Phenocycler_Analysis` clones | remove at G after content checks (all their commits are on GitHub or in A1 bundles) |

## Evidence behind the Gate B items (as presented 2026-10-20)

- **#18** (+1,167): `celltype/uniform_audit.py` — a residual-batch test using exclusive marker
  pairs and reference-based partner positivity, and a marker-invariance audit (`shift_log ~ batch +
  disease`, permutation within batch); 14 synthetic tests; nothing reads it.
- **#19** (+1,217/−59): `sameslide/read_pheno_tsv.py` reads the v3 Parquet export as well as TSV;
  `typing/islet_immune_metrics.py`. The live HPG lane imports the reader.
- **#20** (+966): `phenocycler/qc_flags.py` (Tukey fences on shape, top-quantile DAPI and
  E-cadherin, saturation, edge proximity, artifact region, pan-channel autofluorescence) and
  `neighbour_features.py`; not wired.
- **#21** (+618/−11): `build_galectin_layer.py` from typed cells (`cell_class == "Immune"` and
  `immune_definition`), and the `build_layers.py` default rewired to it.
- **Preserve-branch restain** (8 files, +607/−105): DAPI export crop, `registration_inputs.py`,
  `run_field_refine.py` (+456), smaller edits to nuclei, refine, seed and VALIS registration.
- **H-restain** (217d29a): refined Xenium→IF transform with displacement; declared vs producing
  calibration with a 0.00005 µm/px tolerance; 6436/6442 overlays rejected; 115's 7.198 µm failure;
  app cards x2/x3; conflicts with main in 3 files.
- **112 / 6533 rows**: `pln_1` serial pairs 0041396 / 0041393 and same-slide pairs 0081366 /
  0080870; without them both donors' pLN fall back to the global `sequential` default
  (`phenocycler/integration/modes.py`).
- **Spleen branches**: `copilot/add-model-accuracy-validation` (+410 lines, 3 ahead),
  `copilot/evaluate-windows-11-compatibility` (5 ahead, 48 behind), `fix/he-unet-smoke-empty-split`
  (1 ahead, 40 behind; 3 of its 23 added lines are on monorepo main).
