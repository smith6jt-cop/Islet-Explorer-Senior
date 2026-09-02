# Plan — the planned workflow, and where every auxiliary package sits in it

**Status:** approved 2026-09-02 (plan file `create-a-concise-description-sparkling-flute`). Originally
written against HEAD `1a23105`; branch-convergence and consolidation-completion amendment written
against `31b7380`. Companion documents: `WORKFLOW_2026-09-02_MAP.md` (what runs today),
`EXTERNAL_PACKAGES.md` (the package index), `CONSOLIDATION.md` § 2026-09-02 (consolidation status).
**Decision owner:** maintainer (smith6jt).

This plan positions every auxiliary package identified across sessions, commits, history and memory
inside the consolidation plan's remaining phases (3b, 4, 5, 1-rewritten, 6), ordered by the
maintainer's stated priorities of 2026-08-25: **crossed-batch de-confounding** and **protein + RNA
measured on the same cell**, in service of the immune environment around islets across ND → Aab+ →
T1D. Four forks in the evidence are resolved by stated assumptions (§ Assumptions) so the plan is
executable as written; each is reversible by a one-line edit. The 2026-09-02 continuation also makes
branch convergence and completion of the seven-repository consolidation explicit: `main` is the sole
persistent destination, and the original repositories are not considered absorbed until their unique
work is accounted for, their replacement paths pass, and the originals can be archived read-only.

## Rules in force throughout

No frequency, cross-donor-uniformity or disease-ordering tuning. No analysis that was not asked for.
Every result ships with the figure that shows it, windowed on the measured floor, beside its null.
`Ambiguous` / `Unavailable` carry a recorded reason; a missing measurement is unavailable, not zero.
`object_id` is the only join key. Nothing writes into `data/runs/` except the pipeline; `phenocycler/**`
— which includes `integration/` — re-keys `run_id` whenever it is edited. QuPath annotations are audit
context or registration geometry, never identity. Existing classifications are not a gate: defining a
true positive is the maintainer's call (`1a23105`). The two Shiny apps stay separate.

## Track 0 — Make the record true (first; no science)

1. Inventory and reconcile every branch under Track R before opening the pull request to `main`.
   Do not delete a branch, force-push, or rewrite history merely because its name looks stale.
2. Run `phenocycler adopt plan --from-run 30f6689e4d27d295f890` as a **dry run only** — it will list
   zero adoptable stages — and record in `core/phenocycler/DATA_README.md` that the lanes pin runs by
   content hash and `status` is informational until a run exists under current code. **No `adopt apply`
   and no pipeline run until the last `phenocycler/**` edit of Track C.1 has landed** (Track C.6).
3. Land the three documents, the `CONSOLIDATION.md` section, the root `CLAUDE.md` authority-list entry,
   the `WORKFLOW.md` link, the `.gitignore` pointer, and the corrections listed in `CONSOLIDATION.md`
   § 2026-09-02 "Stale statements corrected". *(Done in this change.)*
4. Environments: mark the `phenocycler-analysis` / `phenocycler_integration` conda recipes HPG-only
   (`HIPERGATOR.md` still activates the former); install `SimpleITK` and `spatialdata-io` into `.venv`
   **additively** — `pip install --dry-run` first, no `--upgrade`, refuse any resolution that touches
   numpy or torch, during a Globus-idle window — which removes two integration test skips and makes the
   non-rigid registration arm and `import_xenium`'s bundle path live.

## Track R — Reconcile repository history to one persistent branch

Current evidence: `main` is at `412603d`; `feat/dct-step0-pinning` is 62 commits ahead of it and ten
commits ahead of its remote; `chore/remove-gold-standard-lane` is a sibling of the feature branch with
three commits absent from it. The latter is not safe to delete by name: `11d80c3` removes the retired
gold-standard lane, `be9406c` records the 2026-08-16 audit, and `b9be3bb` implements the geometry /
identity integration-gate split that Track C.1 still calls planned.

| step | what | exit |
|---|---|---|
| **R.1** protect and ledger | Verify the existing seven-repository bundles, create a current bundle of this monorepo outside its working tree, and record local branches, remote branches, tags, tips, merge bases, `git cherry` results and dirty state. Create a disposition ledger with one row for every commit unique to a non-`main` branch: ported, already present under another commit, deliberately superseded with evidence, or still blocked | every unique commit has an owner and a disposition before history changes |
| **R.2** reconcile the sibling branch | Merge the sibling branch's history into the active line without rewriting it, resolving its three changes deliberately: apply `11d80c3` after confirming no live import depends on `scripts/gold_standard`; carry forward only the still-true facts from `be9406c`; adapt `b9be3bb` to current integration code and tests instead of accepting stale conflict sides wholesale. The geometry/identity split must satisfy Track C.1's current contract, including geometry-only tables with unavailable identity rather than invented labels | `git merge-base --is-ancestor chore/remove-gold-standard-lane feat/dct-step0-pinning` succeeds; the ledger records every conflict resolution; current, research and integration tests pass |
| **R.3** converge on `main` | Push the reconciled feature tip, review the complete `main...feature` diff, and land it on `main` without squashing away the scientific history. No force-push. Run all repository CI plus the explicit lane tests added by Track E | local and remote `main` contain the reconciled feature history and the documentation commit |
| **R.4** remove redundant branches | Before deletion, prove `git cherry main <branch>` is empty and `git merge-base --is-ancestor <branch> main` succeeds for both the feature and chore branches. Then delete their local and `origin` refs. `main` is the sole persistent monorepo branch; any later short-lived review branch is deleted immediately after merge | local branches = `main`; branches under `origin` = `main` (apart from the symbolic `origin/HEAD`) |
| **R.5** keep it converged | Perform all remaining tracks against current `main`; tag or bundle before destructive repository-layout work. A new long-lived topic branch is a regression in consolidation, not a place to park unfinished work | the disposition ledger and branch count are checked at every phase boundary |

The `panels/main` and `phenocycler/main` remote-tracking refs belong to original source repositories,
not to `origin`. They remain only until Track D.4 reconciles and archives those repositories; the
final consolidation exit removes those auxiliary remotes as well.

## Track A — Protein identity (PhenoCycler): finish the cosine lane, then decide the producer

**Declarations needed first** (none may be resolved by tuning):

- **The PCNA leaf** (`scripts/celltype/README.md` § 2026-08-28 regeneration). Options, all declarations:
  drop PCNA from the proliferation gate and let Ki67 carry it; make proliferation a *state* beside
  identity (decision #17's shape — the never-run `states` stage is exactly Ki67/PCNA); or rename the leaf
  for what it measures.
- **The peri definition.** Adopt Track B2 Step 0's continuous signed distance (`0 < d ≤ 50 µm`) for every
  peri-islet statistic, because `cell_region == 'peri'` compares donor 6482's 25 µm annulus with 20 µm
  everywhere else.
- Inherited gate questions stay recorded, not decided: Vimentin in the Mesenchymal gate, CD34 in the
  Endothelial gate, the CD11b alpha-family drift.

The 6539 caveat is **lifted by the data**: run `30f6689e…` and `islet_residents_all.parquet` carry its
repaired ring (this change corrects the two READMEs that still say otherwise).

| step | what | packages |
|---|---|---|
| **A.1** cosine lane cohort pass | Step 3 and Step 4 on all 20 donors (currently 6450, and older than the Step 2 regeneration); `figures_ambiguous` regenerated; `verify_step2` made to write an artifact | UniFORM, UTAG, pSAM scoring |
| **A.4a** release contract, development-stamped (right after A.1) | `data/development/celltype/release/{phenotype_contract.json, label_space.json, calls/donor_id=*/data.parquet}` keyed by `run_id 30f6689e…` and expression `content_id 35fa90fa…`; producer commit stamped; `common_lineage` declared on broad labels only; `prior_source = qupath:Islet_N` recorded on islet-refined calls. The only thing integration and the app bridge may consume. Does **not** lift the registry latch | — |
| **A.2** SOM lane completion (off the critical path) | arms B/C bake-off; node count by the pre-registered seed-stability rule, never by agreement with A.1; cohort run; Step 5 bootstrap vendored from `external/KINTSUGI`; concordance cohort-wide | KINTSUGI, pyFlowSOM, stackview |
| **A.3** Phase 4, executed | the cosine-vs-SOM disagreement figure (endocrine resolution; the no-rule cells), with SpleenFollicle H16/H20 as historical data points once absorbed; a written producer choice | — |
| **A.4b** final release after A.3 | the declared producer; the latch stays down (a production typing workflow needs a new reviewed registry and code contract, per root `CLAUDE.md`) | — |

DeepCellTypes stays closed. Nimbus, KRONOS, TACIT and CellSighter are **not** scheduled; any of them
would enter only as a pre-registered pilot with a `gates.py` verdict, as DeepCellTypes did.

## Track B — RNA identity (Xenium): give the established typer its input

| step | what | packages |
|---|---|---|
| **B.1** | Add `--paths` to `hierarchical_phenotype.py` (mirror `finalize_panel.py:197`), make a missing input an error instead of a silent skip, fix `insulitis_analysis.py`'s flat path. Then run `build_phenotyped_base.py --donors 6436 6442` → `finalize_panel.py --paths` → typer (`scvi-env`, `--dl-workers 0`). This **is** the approved re-stain plan's Phase B and validates the fix | scvi-tools / scANVI, scanpy, core/panels |
| **B.2** | Extend to the other 18 donors, per-donor, the mixed ND/non-ND Xenium instrument runs first — MM035 (115, 6591), MM041 (6461, 6566), MM042 (6493, 6521), MM043 (6482, 6551), MM044 (6436, 6539), MM045 (6442, 6593) — then MM039, MM046, MM047, MM038, MM040. Report `broad_lowsignal` and `broad_n_identity_genes` so anchored vs propagated calls are visible. Assumption Q1 | same |
| **B.3** | Add a dated standing note beside decisions #16/#17 recording the 6539 null-gate measurement (5 of 8 lineages inside a detection-matched null) and the 2026-08-25 KODAMA decision beside the 2026-08-28 typer route. The `04_immune_B_plasma` core split (XBP1 = 39 % of the core mean) is a **flagged maintainer decision, one line**: accepting it moves the `076d1192` pin and the `0348d149` snapshot stamp on every Xenium output; the typer's runtime specificity filter already drops pan-ER genes, so it does not block B.2 | core/panels |
| **B.4** | KODAMA + nullgate is **parked** (assumption Q1). If reaffirmed it re-enters as a pre-registered comparison: k chosen by seed stability, never by agreement with the typer; the nullgate's ~30 % family-wise false-positive rate on noise treated as a gate; `--max-cells` fixed; M-parallelism measured on an idle machine first; HiPerGator needs a key | KODAMA, Rnanoflann |

## Track C — Cross-modality: the two prioritised gains

| step | what | packages |
|---|---|---|
| **C.3a** crossed-batch figure, Xenium-only — **now** | Panel C: the PhenoCycler-batch × Xenium-run occupancy matrix with donor counts and status composition, from `core/phenocycler/data/integration/xenium_paths.csv` (run = third segment of `source_path`) × `data/donor_acquisition_batch.csv`. Panels A/B: the Xenium lines plus a 1,000-draw status-permutation null once B.2 lands. On the figure: Aab+ n = 3 (6450, 6505, 6521); one ND vs one non-ND donor per mixed run. Implementation: `donor.comparable(cfg)` replacing the module-level `COMPARABLE` (`donor.py:47`); new `figures.crossed_batch_trend_figure`. Needs no registration, no matching; survives QC FAIL | — |
| **C.4** re-stain lane — the remaining approved phases | The lane already has: registration at 0.55 µm; `channel_check.py` (registration-free per-channel islet enrichment, i.e. Phase A); the AUROC re-run at the refined transform (bar 0.75 not met); a label-free rank association against the 18-gene endocrine core (`1a23105`); Explorer bundles. Remaining: polygon-exact zonal extraction on the Xenium cell outlines with a cytoplasm-only arm and ≥ 50 % area-overlap islet membership (NaN, never 0, for empty cells); regenerate `validate.py`; retire the void artifacts; figures A (islets at 0.25 µm/px with warped outlines), B (cells filled by value), C (the association). The typer's `celltype_fine` is **reported beside** the association as audit context, not used as a positive class — a classification statistic and any protein call wait for a maintainer-supplied definition of a true positive (assumption Q4). Repoint `INTEGRATION.md` § re-stain and `postxen.py` at this lane. Extension beyond 6436/6442 needs, from the maintainer: QuPath Tissue/Islet GeoJSON exports (2 of 20 exist), `restain_sections.csv` rows for the 15 qptiffs transferred 2026-08-31, and `morphology_focus` transfers only if the image registration track is ever needed | VALIS, geopandas, shapely, tiff_window (from the dct lane, see E.3) |
| **C.3b** crossed-batch figure, PhenoCycler line | from the A.4a development-stamped release with its caveats printed on the figure (assumption Q2), replaced by A.4b when it exists | — |
| **C.1** integration geometry split | conftest factory; contract profiles (`geometry` / `identity`), nullable `evidence_positive`, schema version in parquet metadata; `sections.resolve` from `qupath_manifest.json` with a tracked `section_overrides.csv`; repoint `_discover_pheno_sections` and the new `export_geometry` at `runs/<id>/`; `phenospace.py`; `GEOMETRY_STAGES` / `IDENTITY_STAGES`. Exit: 20 `cells_pheno/donor_id=*/roi=panc` partitions; `register` and `qc` no longer BLOCKED; `hormone_seed_mask` never reaches a geometry-only table | spatialdata, geopandas, shapely, opencv |
| **C.2** registration on all 20 donors — after C.3/C.4 | point-set track first (QuPath islet geometry as registration geometry only; bundles present); the DAPI ↔ DAPI image track only if the point-set gate fails (2 of 20 bundles have `morphology_focus`; do not request 18 transfers up front). The Xenium islet producer is `structures.hormone_seed_mask` reading the typer's `celltype_fine` (its `Endocrine` fallback over-merges), so C.2 depends on B.2 and must state which branch ran | SimpleITK (optional), opencv, scikit-image |
| **C.6** re-key once | after the last `phenocycler/**` edit: run the pipeline (CPU, `phenocycler/**` frozen; `--through states` if the PCNA declaration chose "state") or `adopt apply`; repoint `RUN_ID`, `EXPRESSION_CONTENT_SHA256` and every `data/development` output together (a completed run rewrites `LATEST`, which the celltype tests assert against) | — |
| **C.5** `xencontext` | the aggregate-only Xenium → PhenoCycler artifact (boundary-relative annuli; permuted-pairing, pairing-swap and rotation nulls) after C.2 and C.3. `crossmodal` S9 stays off with its reason recorded | — |

## Track D — Finish combining the repositories (Phases 3b, 5, 1-rewritten, 6)

This executes the unfinished portion of `docs/CONSOLIDATION.md`, which is the tracked execution
record for the earlier seven-repository plan. Importing a tree is not sufficient: unique commits,
working-tree changes, runtime entry points, tests, documentation ownership and archive status all
have to be resolved before a source repository counts as absorbed.

| step | what |
|---|---|
| **D.1** Phase 3b | `SpleenFollicleCounterQP` → `analysis/follicle/`: commit or stash its 11 dirty entries on a branch, fetch (GitHub is ahead, pushed 2026-08-24), then tree-copy as Xenium was — **never `git reset --hard` or `git clean -xfd` there** (a 12 GB ignored `Measurements/Cells.csv`); replace `extract_sample_id` / `_ALT_TO_HANDEL` (`analysis/data_utils.py:62`, `:48`) with a manifest lookup; the H&E / DAPI U-Net lane untouched; H16/H20 labelled historical for A.3; write `analysis/follicle/CLAUDE.md` |
| **D.2** Phase 5 | Complete the R-lane ownership table rather than copying one app wholesale: extract `apps/R/ieCore` from Senior's shared files, add LymphVizkit's `00_domain_config.R` and `config/lymphoid_default.yml` as `DOMAIN` fields (`PIXEL_SIZE_UM`, phenotype vocabulary), and harvest the selected Junior/modernist fixes under assumption Q3. Create `apps/lymph/` as configuration plus data, not another fork. Exit requires tab-by-tab parity tests and one implementation of every shared helper |
| **D.3** Phase 1, rewritten | the adapter root `CLAUDE.md` specifies, in `core/bridge/` (outside the hashed package): A.4b release + `expression` + geometry → run-aware H5AD (`uns` carries `run_id` and content ids) → islet / spatial / trajectory aggregates → `data/app_data/`; replaces the input-starved `scripts/senior` chain; `Ambiguous` / `Unavailable` stay first-class. Gated on A.4b plus a written declaration in `CLAUDE.md` § Downstream application boundary |
| **D.4** reconcile and archive source repos | Diff the imported `core/phenocycler` subtree against `Phenocycler_Analysis` `fc06c0f` / PR #18 commit by commit; port each unique change or record why current monorepo code supersedes it. Verify `core/panels` against `XeniumPanelExplorer@076d1192`, and audit the absorbed Xenium tree against its recorded source. After D.1–D.3 pass, give each of `Islet-Explorer`, `Phenocycler_Analysis`, `Xenium_Analysis`, `XeniumPanelExplorer`, `SpleenFollicleCounterQP` and `LymphVizkit` a final pointer to the monorepo and archive it read-only. Remove the `panels` and `phenocycler` remotes only after recording their final URLs and tips. The production symlink remains untouched unless the maintainer separately authorizes cutover |
| **D.5** viewer | `.tmap` for the 20 eligible donors, the remaining OME-TIFFs, the Avivator-removal decision — a listed maintainer decision beside E.4, not scheduled |
| **D.6** target-layout closure | Audit the original target tree: `core/phenocycler`, `core/panels`, `core/bridge`, `analysis/xenium`, `analysis/follicle`, `apps/R/ieCore`, `apps/islet`, and `apps/lymph` must all exist with one documented owner. Remove superseded runnable copies only after their replacement tests pass; preserve history in Git rather than parallel live entry points. Reconcile the retained rollback assets and create the first `phenocycler-v*` tag after the subtree is final |
| **D.7** consolidation exit | From a fresh clone of `main`, provision dependencies, run the Python and R suites, build one donor from source manifests through the bridge into a rendered app view, and verify provenance in the UI. Final Git state: clean tree, one local branch, one `origin` branch, only the `origin` remote, no submodules, and all six source repositories archived read-only |

## Track E — Package hygiene

1. `core/phenocycler/setup.sh`: replace the no-op `git submodule update --init --recursive` with the
   three `external/` clone + checkout commands at their pins and the `external/.venv-som` recipe. Today
   the only provisioning is prose in two READMEs and a `pyproject.toml` comment. Note the two KODAMA
   installers: pinned (`smith6jt-cop/KODAMA@d2dbe30e`) in `scripts/xenium_celltype/hpg/setup_envs.sh`,
   unpinned (`tkcaccia/KODAMA` HEAD) in `analysis/xenium/hpc/setup_envs.sh` — HPG-only, untestable here.
2. CI: run `core/phenocycler/scripts/*/tests` (celltype and kintsugi tests provisioned from E.1, or
   skipping with a stated reason) and add an `analysis/xenium/**` job for its 77 tests. `core/panels`
   testthat stays a listed gap.
3. Move `gates.py`, `tiff_window.py`, `window_mask.py` (and the `label_space.py` idiom) out of the
   closed `scripts/dct/` into `core/phenocycler/scripts/common/` — never into `phenocycler/` — and drop
   the dead `require_evidence` import in `celltype/usability.py`; then the review document's discard
   instruction is executable.
4. Maintainer decisions, listed not made: the three modified-but-unreferenced forks (maxfuse +29,
   samapi +8, Celloscope +4); whether isletscope, spatialbarrier, smuLymphNet or HPCsizer enter as
   pre-registered pilots or stay prior art; which dormant conda environments to delete.

## Sequencing

| row | items | needs from the maintainer | note |
|---|---|---|---|
| 1 | Track 0 → R.1 → R.2 → R.3 → R.4 | — | branch convergence precedes further implementation; no branch is deleted before its unique commits are resolved |
| 2 | A declarations (PCNA leaf, peri) ‖ B.1 ‖ C.3a Panel C | the two declarations | Panel C needs no identities |
| 3 | B.1 typer on 6436/6442 → B.2 (mixed runs first) ‖ A.1 ‖ C.1 | Q1 reaffirmed | |
| 4 | C.3a Panels A/B + null ‖ C.4 remaining phases | GeoJSON / rows / morphology for the other 18 only | **the headline figure lands here** |
| 5 | A.4a → C.3b protein line ‖ A.2 → A.3 → A.4b | Q2 reaffirmed | |
| 6 | C.6 re-key → C.2 point-set → C.5 xencontext | HPG only if KODAMA is reaffirmed | |
| 7 | D.1 → D.2; Track E | Q3 | repo absorption can proceed while the identity producer finishes |
| 8 | D.3 → D.4 → D.6 → D.7 | the version number, project name, and any production-cutover authorization | D.3 waits on A.4b; archival waits on verified replacements |

## Assumptions taken where the evidence forks

Stated so the plan is executable; each reversible by editing the line it governs.

- **Q1 — Xenium producer.** The established two-tier scANVI typer, per-donor, on all 20 donors: it is
  the route `48c15ef` built the input for and the approved re-stain plan uses; the 2026-08-25 KODAMA
  route is gated on an unresolved cost model and a denied HiPerGator route. KODAMA stays parked as a
  pre-registered comparison (B.4).
- **Q2 — Protein line.** The crossed-batch figure may consume the development-stamped cosine-lane
  release with its caveats printed on the figure (`prior_source = qupath:Islet_N`; the PCNA leaf);
  denominator = islet-associated cells, which is what the study question is about.
- **Q3 — Junior and Phase 6.** Harvesting *from* Junior/modernist into the monorepo is allowed (the
  two-apps rule forbids changes *to* Junior); Phase 6 records the serving facts and repoints nothing.
- **Q4 — Protein call.** No cell is called INS+ or GCG+ until the maintainer defines a true positive
  (`1a23105`). Until then the lane reports the label-free association; the approved rank-contrast call
  and `9e88272`'s mixture branch both stay unrun, and whichever the maintainer picks retires the other
  in the same change (root `CLAUDE.md`: no contradictory runnable paths).

## Verification for this document's own claims

```bash
cd core/phenocycler
../../.venv/bin/python -m phenocycler.pipeline status --config config.ini
../../.venv/bin/python -m phenocycler.integration.pipeline --status
../../.venv/bin/python -m phenocycler.integration.manifest --verify-paths
../../.venv/bin/python -m phenocycler.pipeline adopt plan --config config.ini --from-run 30f6689e4d27d295f890 --out /tmp/adopt.json
grep -rln "Phenocycler_Analysis\b" ../../README.md ../../scripts/CLAUDE.md docs/HIPERGATOR.md docs/INTEGRATION.md scripts/review/README.md   # expect none

# Branch/repository convergence (run from the repo root at the applicable phase boundary)
git cherry main feat/dct-step0-pinning                    # expect no output before deletion
git cherry main chore/remove-gold-standard-lane           # expect no output before deletion
git branch --format='%(refname:short)'                    # final: main
git for-each-ref --format='%(refname:short)' refs/remotes/origin  # final: origin/main (+ origin/HEAD)
git remote                                                # final after D.4: origin
test -d core/bridge -a -d analysis/follicle -a -d apps/R/ieCore -a -d apps/lymph
```
