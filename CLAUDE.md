# Islet Explorer Senior

## Project scope

Three related but separate systems in one repository:

1. **`core/phenocycler/`** — the PhenoCycler source workflow, from QuPath exports to one
   authoritative intensity per marker, plus the development lanes under `scripts/`.
2. **`analysis/xenium/`** — the Xenium 5K lane (absorbed 2026-08-28). Two-tier identity-core →
   scANVI typer, the panel machinery reading `core/panels`, a dormant spatial-domain benchmark, a
   HiPerGator-only `hpc/`. Runs in `scvi-env`, **not** `xenium_analysis`. Its numbered
   do-not-regress decisions #1–#17 are cited by number elsewhere and **must not be renumbered**.
3. **`apps/islet/`** — the Shiny app. Shared R lives once in `apps/R/ieCore`; display contracts live
   in the per-domain YAML under `apps/islet/config/`. `apps/lymph` is configuration plus a synthetic
   fixture, not another app.

The panel taxonomy has **one** authority: the `core/panels` subtree.
`analysis/xenium/scripts/fetch_subpanels.py` regenerates a gitignored snapshot from it. Never commit
a second copy of those CSVs; never hand-edit the snapshot.

## Non-negotiable analysis rules

- Never infer positivity from a desired frequency, cross-donor uniformity, or disease ordering.
- Existing QuPath classifications, annotations, phenotypes, clusters, disease labels and
  gold-standard reviews are **audit context only** — never training, tuning or release-gating inputs.
- **A cell with no nucleus is excluded from every calculation**; we cannot be sure it is a cell.
- Preserve the QuPath UUID as `object_id` throughout. Never substitute centroid matching, row order,
  filenames or donor-name parsing.
- A missing measurement is unavailable — not zero, not negative.
- Donors `6457` and `6579` are centrally excluded in `core/phenocycler/phenocycler/cohort.py`.
  Configuration must not bypass it, and no lane may restate the list.
- Acquisition batch is aliased with disease status at the coarse level; correct only on the finer
  variables (`stain_group`, `acq_run` — see `data/donor_acquisition_runs.csv`), never on the Fusion
  version, and stratify where correction is not warranted.
- No result is a verification target. A check that passes by reproducing a previous biological
  finding measures the expectation, not the build.

## Vocabulary

Say **object** (a QuPath pixel-classifier region), **reference level** (a value derived from a
reference population), **peri-islet expansion** — not polygon, bar or ring. Say *method* or *the
comparison*, never "arm" or "bake-off". You cannot indict a marker, only a method.

## The source workflow

```text
QuPath cohort manifest -> ingest -> geometry QC -> REDSEA -> expression
```

Run it **from `core/phenocycler/`**, which owns the only `config.ini`. A repo-root `config.ini` is a
bug — one existed until 2026-08-10 and silently supplied stale defaults.

```bash
python -m pip install -e core/phenocycler
cd core/phenocycler
python -m phenocycler.pipeline run --config config.ini --pipelined --through expression
```

Ordinary execution stops after `expression`. `controls`, `calibrate`, `type` and assignment export
are blocked by `NO_VERIFIED_TYPING`, and that latch has not moved.

Raw CSV, qptiff, geometry, calibration, channel mapping, panel identity and segmentation version
belong in `core/phenocycler/data/qupath_manifest.json`. Results live under
`core/phenocycler/data/runs/<run_id>/`; stages and runs are immutable and content addressed. Mutable
review artifacts live under the repo-root `data/development/` — a level above `core/`, never writing
a stage manifest and never mutating `data/runs/`.

## Cell typing, as it stands

The current typing is the **subpanel classification**: context from QuPath objects, then one
question per context answered by the markers that can answer it, against reference levels drawn from
the ordinary tissue. 25,150,537 cells, all 20 donors, in
`data/development/typing/subpanel/subpanel_calls.parquet`.

- **How it works:** `core/phenocycler/scripts/typing/RECONSTRUCTION.md` — recovered by reproducing
  the calls from their inputs and verified against them. Read it rather than re-deriving it.
- It is a **development lane**. Nothing here is a release, and the latch stays down.
- Retired, and not to be revived as fallbacks: RESTORE pair webs, hierarchical / simultaneous /
  probabilistic typing, marker-debug and Otsu lanes, the Gate 2–5 apparatus, the broad-typing
  producer comparison withdrawn 2026-09-08, and flat `cells_redsea*` / `restore*` roots. Git has them.
- **Never claim acinar or ductal.** The 60-plex carries no amylase, trypsin, CPA1 or elastase, and
  no CK19, CA19-9, SOX9, HNF1B, CFTR or MUC1. Exocrine is one label.

Authoritative documentation: `core/phenocycler/README.md`, `docs/WORKFLOW.md`, `DATA_README.md`,
`docs/METHOD_HISTORY.md`; and at repo level `docs/WORKFLOW_2026-09-02_MAP.md` (what runs today),
`docs/PLAN_2026-09-02_AUXILIARY_PACKAGES.md`, `docs/EXTERNAL_PACKAGES.md`.

## The app

Two datasets, deliberately separate, selected by **`ISLET_DOMAIN_CONFIG`**:

| domain | config | data root |
|---|---|---|
| legacy CODEX, 15 donors | `apps/islet/config/islet_default.yml` | `data/app_data` |
| PhenoCycler run cohort, 20 donors | `apps/islet/config/run20_default.yml` | `data/development/typing/app` |

They share five donors and no vocabulary, so they are **never merged** — a cross-donor comparison
spanning both would compare methods, not biology. Every data path resolves through
`ieCore::domain_data_root()`; a path that escapes the active domain's root serves the wrong cells for
a right-looking islet, and the testthat suite asserts against exactly that.

Architecture: entrypoint `apps/islet/app.R`; constants and state `R/00_globals.R`; loading
`R/data_loading.R`; tabs `R/mod_*_server.R`; helpers `R/{drilldown,segmentation,spatial,viewer}_helpers.R`.

Conventions: donor status orders `ND < Aab+ < T1D`; donor colours are centralised in `00_globals.R`;
say **Sex** in labels though the stored column is `gender`; use `scales::pseudo_log_trans(base = 10)`
where zeros must stay visible; core parents are `Islet_N` and the expansion is `Islet_N_exp20um`
(every expansion is 20 µm — derive it from the signed distance, not the annotation); guard
root-level outputs by active tab and keep module inputs namespaced; build a downloadable ggplot once
and share it between the renderer and its `downloadHandler`.

Development server, which must stay exactly two directories deep:

```bash
cd apps/islet && ISLET_DOMAIN_CONFIG=config/run20_default.yml \
  Rscript -e 'shiny::runApp(".", port = 7777)'
```

R **4.5.1** with the `~/R/…/4.5` library — `/usr/bin/Rscript` is 4.6.1 and silently lacks the
dependencies. Production is nginx + shiny-server; see `deploy/DEPLOY.md` and `apps/islet/CLAUDE.md`.
Never treat the development server as a deployment test: Avivator needs HTTP Range, which `runApp`
does not serve.

App-specific UI, viewer, deployment and reticulate guidance: `apps/islet/CLAUDE.md`. Script lanes:
`scripts/CLAUDE.md`.

## Source control

- No submodules. `Phenocycler_Analysis` and `XeniumPanelExplorer` became the `core/phenocycler` and
  `core/panels` subtrees; one commit covers a change and its tests.
- Never commit generated `data/`, run manifests, executed notebooks, caches, or local model paths.
- **`git reset --hard` and `git clean -xfd` are unsafe here** — hundreds of GB of untracked,
  gitignored data live under `data/` and `core/phenocycler/data/`, and none of it is in Git. Roll
  back with `git revert`.
- When replacing a workflow, remove the obsolete entrypoint rather than leaving two runnable ones.
