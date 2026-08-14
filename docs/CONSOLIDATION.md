# Consolidation — seven repos into one project

Execution reference for the plan agreed 2026-08-12. Phases and rationale live in the plan;
this file holds the parts later phases must execute *against*: the target layout, the
per-file ownership decisions, the archive list, and the reconciled rule hierarchy.

> **Read `docs/HANDOFF_2026-08-14_CONSOLIDATION.md` first.** The upstream workflow was cut
> back to end at expression on 08-13/08-14, which voids the plan's Phase 1 (the source→app
> bridge) and requires Phase 4 to be rewritten. The layout, ownership table and archive list
> below are unaffected.
>
> **Status 2026-08-14:** the retirement is **committed and pushed** (four commits on
> `Phenocycler_Analysis` branch `agent/retire-typing-lane`), so Phase 2 is unblocked once that
> branch reaches `main`. Read § "Corrections from the 2026-08-14 audit" at the end of this file
> before executing Phase 2 — it lists five wrong statements in the earlier record and three
> ways to lose data or leak a secret.

Backups: `/home/smith6jt/repo-backups/2026-08-13/` — seven verified bundles plus working-tree
snapshots for the two dirty repos. See that directory's `README.md`.

## Target layout

```
islet-explorer/                  (this repo, reshaped)
  core/phenocycler/      Python package incl. integration/ — the scientific authority
  core/panels/           XeniumPanelExplorer: panel taxonomy data + its app
  core/bridge/           run -> H5AD -> app aggregates
  analysis/xenium/       Xenium notebooks on the shared cell-table contract
  analysis/follicle/     SpleenFollicleCounterQP
  apps/R/ieCore/         shared R package: domain config, theme, stats, helpers
  apps/islet/            the unified Shiny app
  apps/lymph/            a domain config + data dir, not a fork
  deploy/ docs/ data/
```

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
| `analysis/xenium/CLAUDE.md` | the two do-not-regress invariants |
| `analysis/follicle/CLAUDE.md` | inventory + the per-image statistical framework |
| `apps/CLAUDE.md` | Junior's 457-line per-tab reference, plus Senior's viewer-backend section |

### Conflicts that must be resolved in writing, before code moves

1. **Junior's `CLAUDE.md` documents methods the root rules retire** (ordered residual lineage,
   the Gate 2–5 apparatus). Junior's app-level doc survives; its method sections do not.
2. **`SpleenFollicleCounterQP/CLAUDE.md` teaches donor-ID-from-filename** ("image name → sample
   ID by HDL prefix", implemented at `analysis/data_utils.py:62` with a hardcoded
   `_ALT_TO_HANDEL` map at `:48`). The root rule forbids exactly this. The guidance is deleted
   and the code replaced by a manifest lookup in Phase 3.
3. **"There is exactly one `config.ini`, and it is the submodule's"** becomes false the moment
   the submodule is dissolved. Restate as: `core/phenocycler/config.ini` is the pipeline
   authority; `analysis/xenium/config.ini` is a reference mirror whose live values come from
   `XENIUM_*` env vars, and it must be labelled as such or deleted.
4. **Xenium's numbered "do not regress" decisions** are referenced by number across its docs
   ("decisions #16/#17"). Renumbering breaks those references — carry the list across intact.

## Archive list

Set read-only on GitHub after Phase 6, each with a final README commit pointing here:
`Islet-Explorer`, `Phenocycler_Analysis`, `Xenium_Analysis`, `XeniumPanelExplorer`,
`SpleenFollicleCounterQP`, `LymphVizkit`.

## Blockers carried into later phases

- ~~**`Phenocycler_Analysis` has 136 uncommitted files**~~ — **CLEARED 2026-08-14.** The real
  count was **172** (75 deleted, 49 modified, **48** untracked; the "12 untracked" were 9 files
  plus 3 directories). Landed as four themed commits on
  `Phenocycler_Analysis` branch `agent/retire-typing-lane` (`8337e3a..f6953cf`), pushed.
  **`git subtree` reads `main`, so this branch must be merged to `main` before Phase 2 step C.**
- **`Xenium_Analysis` hardcodes `ROOT = /blue/maigan/smith6jt/Xenium_Analysis`** in 12 files.
  Fixed in Phase 3 using the `Path(__file__).resolve().parents[n]` pattern already used at
  `thymic_zonation_analysis.py:54`.
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
