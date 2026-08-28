# Consolidation — seven repos into one project

Execution reference for the plan agreed 2026-08-12. Phases and rationale live in the plan;
this file holds the parts later phases must execute *against*: the target layout, the
per-file ownership decisions, the archive list, and the reconciled rule hierarchy.

> **Read `docs/HANDOFF_2026-08-14_CONSOLIDATION.md` first.** The upstream workflow was cut
> back to end at expression on 08-13/08-14, which voids the plan's Phase 1 (the source→app
> bridge) and requires Phase 4 to be rewritten. The layout, ownership table and archive list
> below are unaffected.
>
> **Status 2026-08-14:** the retirement is merged to `main`, and **Phase 2 is complete**
> (PR #14, branch `consolidate/phase-2`). The layout below is now real: there are no
> submodules and no `.gitmodules`. Read § "Corrections from the 2026-08-14 audit" and
> § "What Phase 2 actually found" at the end of this file before executing Phase 3.

Backups: `/home/smith6jt/repo-backups/2026-08-13/` — seven verified bundles plus working-tree
snapshots for the two dirty repos. See that directory's `README.md`.

## Target layout

`[done]` marks what Phase 2 actually created; the rest is still target.

```
islet-explorer/                  (this repo, reshaped)
  core/phenocycler/      [done] Python package incl. integration/ — the scientific authority
  core/panels/           [done] XeniumPanelExplorer: panel taxonomy data + its app
  core/bridge/           run -> H5AD -> app aggregates     (VOID until a verified run exists)
  analysis/xenium/       [done] Xenium 5K lane: the two-tier typer, the panel machinery,
                         the benchmark harness (dormant), and hpc/ (HiPerGator-only)
  analysis/follicle/     SpleenFollicleCounterQP
  apps/R/ieCore/         shared R package: domain config, theme, stats, helpers
  apps/islet/            [done] the unified Shiny app
  apps/lymph/            a domain config + data dir, not a fork
  deploy/ docs/ data/    [done]
```

`core/bridge/` is listed for completeness only. The handoff records Phase 1 as void: three of
its five inputs no longer exist, and `CLAUDE.md` forbids describing the legacy H5AD rebuild as
current. A consequence to keep in view — `data/app_data/` does not exist, so **the Shiny app
cannot currently load data**. That is a known, accepted state, not a regression from Phase 2.

## Ground rules

1. This repo becomes the monorepo. It keeps its history and its GitHub URL.
2. Squash-import every incoming tree **except** `Phenocycler_Analysis`, which comes in with
   full history via `git subtree` — it is 3.3 MB and it is the scientific record.
3. Archive means archive: originals pushed complete, then set read-only on GitHub.
4. Absorption is not merging. Each phase ends with the old lane deleted, not left runnable
   beside the new one.
5. Data never enters git. Both multi-GB artifacts (`Islet-Explorer/CellCompartments.tsv`
    11 GB, `SpleenFollicleCounterQP/Measurements/Cells.csv` 12 GB) are **already untracked and
   gitignored**, verified 2026-08-13. They stay where they are — `Cells.csv` is read by
   `analysis/` notebooks H14/H15/H16/H17/H20 and moving it would break them. `git subtree`
   imports history, not the working tree, so neither can be swept in.

## R lane — per-file ownership

Senior and Junior share root commit `b78f0fd0` and forked 2026-06-01 at `ef051525`;
Senior is +73 commits, Junior +34, in both directions. Neither is a superset. Modernist
(`claude/modernist-redesign`, tip `f4e5b22`) is a strict superset of Junior HEAD.

| Family | Owner | Why |
|---|---|---|
| `ai_helpers`, `mod_ai_assistant_server`, `utils_safe_join` | any | byte-identical 3-way |
| `mod_ai_assistant_ui`, `mod_viewer_ui` | Senior==Junior | modernist differs by ≤12 lines |
| `segmentation_helpers` (658 ln) | either | 10 diff-lines, all path constants + one size const |
| `viewer_helpers`, `mod_viewer_server` | **Senior** | 9 TissUUmaps/MinIO/OME-Zarr functions Junior lacks |
| `rdeck_helpers`, DuckDB layer, mirai pool, `compute_lamian_pseudotime` | **Senior** | single-owner, additive |
| `utils_stats` | **Junior** | `flag_outliers` zero-exclusion fix; `size_companion_models`, `age_group_companion_model`, `bin_coverage_ok`, `format_p_value` |
| `spatial_helpers` | **Junior** logic, **Senior** backend | same public API over CSV vs DuckDB; keep Junior's `canonical_phenotype_order` |
| `mod_statistics_server` (1,315 diff-lines) | **Junior** | 831 lines longer; 4 correctness commits Senior never received |
| `mod_trajectory_server` | **both** | Junior's rank-axis stability fix + Senior's vectorised `match()` are non-overlapping |
| `mod_spatial_server` | **both** | Junior's download handlers + Senior's rdeck path |
| `drilldown_helpers` | **both** | Senior's 25-phenotype vocabulary + Junior's font fix; 33 lines total |
| `mod_inflammation_*` | **modernist** | edgeR replaced by a CLR compositional test |
| `cohort_filter`, `mod_cohort`, `ui_rail`, `theme_modernist`, `ui_fields`, `ui_modernist_bits`, bslib `app.R` | **modernist** | the shell everything else lands in |
| `00_domain_config.R` + `config/*.yml` | **LymphVizkit** | the generalization layer, already proven on a second tissue |
| `app/shiny_redesign/` | drop | superseded except `www/wiring.css` |

Two divergences are **data, not drift**, and become `DOMAIN` fields rather than constants:

- `PIXEL_SIZE_UM` — 0.5078 (Senior, IO60panc2nd) vs 0.3774 (Junior).
- Phenotype vocabulary — 25 (Senior/scimap) vs 21 (Junior) vs 12 (modernist, post-retirement).

## Rule hierarchy after consolidation

Today: ~2,600 lines of `CLAUDE.md` across 9 files in 5 repos, several mutually contradictory.
Target:

| File | Holds |
|---|---|
| `CLAUDE.md` (root) | project scope, the non-negotiable analysis rules, source-control hygiene, this layout |
| `core/phenocycler/CLAUDE.md` | the workflow as it stands (ingest → geometry → REDSEA → expression), the retired-methods list |
| `core/panels/CLAUDE.md` | panel taxonomy contract only |
| `analysis/xenium/CLAUDE.md` | decisions #1-#17, carried across **intact** (conflict #4 below beats this row's earlier "two invariants" summary -- other docs cite them by number), plus dated standing notes |
| `analysis/follicle/CLAUDE.md` | inventory + the per-image statistical framework |
| `apps/CLAUDE.md` | Junior's 457-line per-tab reference, plus Senior's viewer-backend section |

### Conflicts that must be resolved in writing, before code moves

1. **Junior's `CLAUDE.md` documents methods the root rules retire** (ordered residual lineage,
   the Gate 2–5 apparatus). Junior's app-level doc survives; its method sections do not.
2. **`SpleenFollicleCounterQP/CLAUDE.md` teaches donor-ID-from-filename** ("image name → sample
   ID by HDL prefix", implemented at `analysis/data_utils.py:62` with a hardcoded
   `_ALT_TO_HANDEL` map at `:48`). The root rule forbids exactly this. The guidance is deleted
   and the code replaced by a manifest lookup in Phase 3.
3. **"There is exactly one `config.ini`, and it is the submodule's"** — **now false**, as of
   Phase 2. There is no submodule. `core/phenocycler/config.ini` is the pipeline authority and
   the root `CLAUDE.md` still needs updating to say so (its "Run it from
   `Phenocycler_Analysis/`" instructions are stale). When `analysis/xenium/` lands, its
   `config.ini` is a reference mirror whose live values come from `XENIUM_*` env vars, and it
   must be labelled as such or deleted.
4. **Xenium's numbered "do not regress" decisions** are referenced by number across its docs
   ("decisions #16/#17"). Renumbering breaks those references — carry the list across intact.

## Archive list

Set read-only on GitHub after Phase 6, each with a final README commit pointing here:
`Islet-Explorer`, `Phenocycler_Analysis`, `Xenium_Analysis`, `XeniumPanelExplorer`,
`SpleenFollicleCounterQP`, `LymphVizkit`.

## Phase 3 — status 2026-08-28

**`analysis/xenium/` has landed** (115 files, 3.5 MB). **`analysis/follicle/` has not** and is
the remaining Phase 3 work, including the donor-ID-from-filename violation that conflict #2
above names (`analysis/data_utils.py:48,62`, the hardcoded `_ALT_TO_HANDEL` map).

Four things worth recording, because each contradicts an assumption above:

1. **Imported by tree copy, not `git subtree add --squash`.** `Xenium_Analysis/.git` is
   **889 MB** for a **3.4 MB** tracked tree (figures and `data/output*` committed before they
   were gitignored). `subtree` fetches all of it before the squash discards the graph, and
   recovering that needs `gc --prune=now` -- which this repo's own rules call unsafe with
   157.9 GB of untracked ignored data present. A `git checkout` of the tracked paths into one
   commit *is* a squash import and satisfies ground rule 2's intent, and it lets the pruning
   (~20 dropped files) be part of the import instead of a second deletion commit.
2. **`core/panels` IS the identity-core source, verified.** Everything
   `fetch_subpanels.py` fetched over HTTP from `XeniumPanelExplorer@076d1192` is present in
   the subtree -- 4 phases, 129 files, **0 missing** -- and `panel_roles.csv`,
   `panel_gene_roles.csv` (1,613,132 bytes) and the identity cores are **byte-identical** to
   the upstream pin. So the script now reads locally: one authority, one regenerable
   gitignored cache stamped with `git rev-parse HEAD:core/panels/data` (`0348d149...`), and a
   test that fails when the two drift.
3. **The `xenium_analysis` conda env is broken here** (`ImportError: cannot import name
   '_spropack'`). The lane runs in `scvi-env`, which has scvi-tools 1.4.0.post1 and CUDA.
   `environment.yml` was dropped rather than carried across as a recipe for a broken env.
4. **Two pre-existing defects surfaced on absorption and were fixed**: `cm.get_cmap` in
   `utils/xenium_explorer_export.py` was **removed in matplotlib 3.9** and its `ImportError`
   guard could not catch the resulting `AttributeError`, so the helper crashed rather than
   degrading on both local environments; and `tests/test_repo_health.py` asserted an
   `environment.yml`, a `config.ini` and a SLURM conda env that no longer exist.

## Blockers carried into later phases

- ~~**`Phenocycler_Analysis` has 136 uncommitted files**~~ — **CLEARED 2026-08-14.** The real
  count was **172** (75 deleted, 49 modified, **48** untracked; the "12 untracked" were 9 files
  plus 3 directories). Landed as four themed commits on
  `Phenocycler_Analysis` branch `agent/retire-typing-lane` (`8337e3a..f6953cf`), pushed.
  **`git subtree` reads `main`, so this branch must be merged to `main` before Phase 2 step C.**
- ~~**`Xenium_Analysis` hardcodes `ROOT = /blue/maigan/smith6jt/Xenium_Analysis`** in 12
  files.~~ **CLEARED 2026-08-28.** The real count was **39 files** (19 Python `ROOT` /
  `PROCESSED_ROOT` constants, 2 R, the rest shell/SLURM/notebook). Fixed with
  `Path(__file__).resolve().parents[n]` as planned, but on **two axes, not one**: `ROOT`
  (code) and `DATA_ROOT` (h5ads, `XENIUM_DATA_ROOT`, default `/data/xenium`). A single root
  would have put `data/processed/` under `analysis/xenium/` -- 60+ GB inside the repo, in a
  directory the anchored `/data/` ignore does **not** cover. `_panel_roles.py` needed one
  extra fix its module-level `SUBPANEL_DIR` had no override for. Asserted by
  `analysis/xenium/tests/test_repo_health.py`.
- **`XeniumPanelExplorer`'s live `renv.lock` and `.Rprofile` are missing** — only
  `.r46-backup` copies remain, so `renv::restore()` fails as documented in its README.

## Corrections from the 2026-08-14 audit

Five statements elsewhere in the consolidation record are wrong. Correct versions:

1. **`SpleenFollicleCounterQP` was not "2 dirty binaries".** Its checkout had **68 worktree
   deletions of tracked source** — all of `analysis/` (the H1–H20 notebooks and
   `data_utils.py`), `README.md`, `CLAUDE.md`, and the genotype/measurement spreadsheets — which
   the 08-13 2-entry tarball did not capture. Nothing was lost (HEAD `2040bb1` == `origin/main`);
   the checkout was **restored from HEAD on 2026-08-14**. It lives at
   `/home/smith6jt/SpleenFollicleCounterQP`, *not* `/data/SpleenFollicleCounterQP`.
2. **There is no three-way conflict for `core/panels`.** Phenocycler's gitlink `076d1192` **is**
   `XeniumPanelExplorer`'s `origin/main` (ancestor *and* zero distance). The local
   `/data/XeniumPanelExplorer` @ `0a23564` is simply that same line 11 commits back. Import
   `origin/main`; it is both the pinned/tested input and the published tip.
3. **`app/panel_app/` did exist.** It was removed in `544dd10`. `panel-app-ci.yml` is still dead
   and still gets deleted, but do not write "never existed" — it is checkable and false.
4. **Phenocycler has no `v*` release tags** (`git describe` → bare sha; the one tag,
   `archive/superseded-multiref-median-8gate`, is unreachable from `main`). So there is no tag
   history to migrate — but after the subtree import setuptools-scm walks up to the **monorepo**
   `.git`, so the version silently starts counting monorepo commits. CI cannot catch this
   (`__version__ != '0.0.0+unknown'` keeps passing). Namespace it: `root = "../.."`,
   `--match phenocycler-v*`, and cut `phenocycler-v0.4.0` right after the import.
5. **Renaming Senior's `app/shiny_app` cannot affect production.**
   `/srv/shiny-server/islet-explorer` → `/home/smith6jt/panc_CODEX/Islet-Explorer/app/shiny_app`,
   which is **Junior**. The Phase 6 constraint is satisfied by doing nothing.

### Two traps that will destroy data or leak a secret if missed

- **Never `git submodule deinit -f`.** `-f` makes it `rm -rf` the submodule working tree,
  including ignored content — that is the **148 GB** of real data at
  `Phenocycler_Analysis/data/`. Instead `mv` the whole submodule out of the repo first (same
  filesystem, so `rename(2)` is O(1) — `/home` has only 185 GB free, so a copy is not an
  option), then `git rm --cached` + drop the `.gitmodules` entry. Restore the data **last**, so
  `git reset --hard` stays a safe rollback until then.
- **`app/shiny_app/.Renviron` holds a live UF Navigator key** (`deploy/DEPLOY.md:177`) and is
  ignored by an *exact-path* rule. The six `app/shiny_app/...` lines in `.gitignore` must be
  repointed to `apps/islet/...` **in the same commit as the `git mv`**, or the next `git add -A`
  commits the secret.
- The root `.gitignore` `data/` pattern (no leading slash) matches at **any depth**, so it would
  shadow `core/phenocycler/data/` and make that subtree's own `!/data/integration/` negation
  unreachable. Anchor it to `/data/`, and anchor `external/` and `HANDOFF*` too.

## What Phase 2 actually found (2026-08-14)

Executed as PR #14. Order run: preflight → anchor `.gitignore` → dissolve submodule →
subtree `core/phenocycler` (full history, 117 commits) → subtree `core/panels` (squashed) →
drop the nested gitlink → CI → `git mv` the app → restore the data last.

**The data survived intact.** 157,857,924,055 bytes, 48 top-level entries, top-level listing
identical to the pre-flight invariant. Both runs present (`3bf4fc8e…` 99G, `c8412757…` 121G,
all stages). `df` was unchanged across both moves, confirming `rename(2)` rather than copy.

Four things were true only once the work was underway:

1. **The subtree import fixed the private-submodule problem outright.** `core/panels` is now
   vendored, so the taxonomy is always present. Integration went from *214 passed / 18
   skipped* to *229 passed / 3 skipped*, and the `PANEL_EXPLORER_TOKEN` workaround plus the
   non-fatal fetch are no longer needed — the strict zero-skip contract is restored.
2. **setuptools-scm broke silently on the move, in the opposite way to what was predicted.**
   It does *not* walk up to find a git root and start counting monorepo commits; it finds no
   root at all and emits `fallback_version` (`0.0.0+unknown`) — which CI asserts against, so
   this was a red build waiting to happen, not a cosmetic version drift. Fixed with
   `root = "../.."` plus `phenocycler-v*` tag namespacing.
3. **The nested gitlink was worse than inert.** A mode-160000 entry at
   `core/phenocycler/external/XeniumPanelExplorer` made `git submodule status` fail outright
   ("no submodule mapping found"), because git reads `.gitmodules` only at the worktree root.
4. **`.Renviron` was already safe.** The exact-path ignore rule was not its only protection —
   a bare `.Renviron` pattern matches at any depth. The exact-path lines were repointed
   anyway, since silently-stale rules are not something to leave to chance.

Confirmed as designed: the three grep forms find genuinely different sets (the split-pathlib
site at `scripts/senior/tissuumaps_smoke_test.py:43` is invisible to a literal search); the 16
depth-2 R data paths survive the move because `apps/islet` is still exactly 2 levels deep, and
`islet-r.yml` now asserts that; `git subtree` cannot sweep in ignored data, because
`read-tree --prefix` reads a commit's tree and never consults the working directory.

### Left deliberately undone

- **No `phenocycler-v*` tag.** No `v*` tag has ever existed in that project, and a published
  version cannot be reused, so the first version number is a decision rather than a guess.
  The current version is `0.0.1.dev<N>+g<sha>`, which satisfies CI and still pins the exact
  commit; the dev distance counts monorepo commits until a tag exists.
- **`/home/smith6jt/IO60panc2nd/Phenocycler_Analysis.preserved` and
  `.git/modules/Phenocycler_Analysis` are retained** as rollback assets until Phase 2 sign-off.
- **Production untouched.** `/srv/shiny-server/islet-explorer` points at the Junior checkout,
  so this rename structurally cannot affect it. Cutover remains Phase 6.

> **After the data was restored, `git reset --hard` and `git clean -xfd` are no longer safe in
> this repository** — 157.9 GB of untracked, ignored data now sits under `core/phenocycler/data/`.
> Roll back with `git revert`.
