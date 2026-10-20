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
- Donor `6457` is centrally excluded in `core/phenocycler/phenocycler/cohort.py`. Configuration
  must not bypass it, and no lane may restate the list — derive from `DONOR_EXCLUSIONS`, never
  from a literal. **`6579`'s exclusion was lifted 2026-09-17** and it is recorded in `REINSTATED`
  alongside why; the pancreatitis it was excluded for is still a disease context, not a reason to
  drop it. Both cohort manifests have declared it since 2026-09-23 (with its revised export,
  `data/RevisedQuPathExports/6579.tsv`), but it is not yet a 21st donor in any stage root — no stage
  has run on it; see `docs/DONOR_6579_REENTRY_2026-09-17.md` for what re-entry takes. `pln_node_annotations.csv` now records it `excluded=False`;
  `data/donor_acquisition_runs.csv` still carries `excluded=True` and disagrees with `cohort.py`
  until its generator reruns.
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

**The content-addressed run framework was removed on 2026-09-16.** Stages write to stable roots —
`data/cells`, `data/geometry_qc`, `data/redsea`, `data/expression` — and a donor is done when its
partition has content. There is no `run_id`, no stage manifest, no donor receipt and no fingerprint
validation. A stage whose output is already present is skipped; `--force` reruns it.

Raw CSV, qptiff, geometry, calibration, channel mapping, panel identity and segmentation version
belong in `core/phenocycler/data/qupath_manifest.json`. That manifest still records a `content_id`
and per-file sha256, but **nothing validates them** — it is a description, so it can be hand-edited
when a source file moves. Results live under `core/phenocycler/data/<stage>/`. Runs made before
2026-09-16 remain under `core/phenocycler/data/runs/<run_id>/`; that tree is historical, nothing
writes to it, and the 20-donor cohort output still lives only there. Mutable review artifacts live
under the repo-root `data/development/` — a level above `core/`.

## One entry point: `make`

`make help` lists every target; these are the ones a session uses, run from the repo root.

| command | what it does |
|---|---|
| `make test-fast` | the default suite without the data-state guards (`-m "not data"`) |
| `make test-data` | only the data-state guards; red means an artifact is stale |
| `make test-lanes` | each development lane's own tests, lane by lane |
| `make check` | is every layer/app artifact newer than its inputs? Reads mtimes only, and names every stale input, a run that did not finish, and a stale `typing/immune` |
| `make plan` | the rebuild the driver would run, in order; runs nothing |
| `make rebuild` | refuses unless `GO=1`; it writes `data/development/typing/app` and `/data/islet_layers` |
| `make env-check` | the `.venv` resolves `phenocycler` to this checkout |
| `make producers` | prints `docs/PRODUCERS.md` |

## Cell typing, as it stands

The current typing is the **subpanel classification**: context from QuPath objects, then one
question per context answered by the markers that can answer it, against reference levels drawn from
the ordinary tissue.

A **reference level** is the p99 of that marker set's max over the reference population — that
donor's cells in no pixel-classified object and outside every islet — and it is **fitted per
donor**. Pooling was measured wrong on 2026-09-15: the hormone level spreads 19.5x across donors
and separates by acquisition batch (batch 1 median 0.0418, batch 2 0.3167). The four pooled levels
(hormone 0.1010, Nerve 0.2272, Lymphatic 0.1776, Blood vessel 0.3073) belong to the OLD table and
survive only so `subpanel.py --gate` can reproduce it.

**The islet chain**, for a cell inside an islet outline, after the maintainer's ruling of
2026-09-18:

1. **Structure** — four steps, in order, none sufficient alone. (a) the cell's centroid is exactly
   inside the object's bounds; (b) **the whole cell is then checked — the majority of the cell's
   area must lie within the object**, which is what makes it a candidate; (c) the marker that
   *created* the object is simply higher than all the rest (`Nerve → B3TUBB`, `Lymphatic →
   Podoplanin`, `Blood vessel → Caveolin, CD34` — three objects, four markers); (d) otherwise
   **passthrough, carrying the fraction of the cell's area the object covers**, not a boolean.
   A vessel crossing a cell does not make the cell a vessel. Step (b) is the test that separates
   the two and it was absent from the code until 2026-09-19. Do not re-derive this rule or restate
   it as a choice — it is settled, and `structure-call-rule` in memory has the maintainer's words.
2. **Exocrine inside an islet is an annotation error** — `Excluded`, with the reason on the row.
3. **CD99 is the key differentiator**, at a level fitted **per donor and per size stratum** (1 /
   2–3 / 4–10 / >10 admitted cells per islet) from that donor-and-stratum's hormone-confirmed
   cells. A **floor guard** refuses a level at or below the feature's own floor (0.0 for a scaled
   marker, 1.0 for an elongation that is a diameter ratio) and returns it UNAVAILABLE, because
   `apply_scale` clips at background and a low quantile of a population sitting there returns 0.0,
   which admits every measured cell while looking fitted. The donor's ordinary-tissue CD99 level is
   carried alongside as `donor_reference_level` and **reported, never substituted** for the
   per-stratum level.
4. **Elongation is a second route and can only add** — a cell CD99 does not reach can still be
   endocrine if its nucleus is as elongated as that stratum's confirmed endocrine cells. Large
   stratum only (strictly more than 10 admitted cells), barred below it. It is
   `nucleus_elongation`, not a cell-boundary feature: cell circularity runs 0.514–0.523 across every
   class, a tessellation artefact.
5. **CD99 decides *whether*, never *which*.** The name is the hormone argmax over INS/GCG/SST, and
   at least one hormone must have been READ before anything may name a cell — an argmax over three
   NaNs returns index 0, "Beta", every time.
6. **Not endocrine** takes the ordinary chain: immune definitions, the supporting tier, then
   Exocrine.

**Current Single exception (human decisions, October 2026):** EndocrineSingle
uses its unique highest usable hormone directly, without CD99. On 2026-10-05
the human superseded retention of missing/tied subtypes with exclusion: the 25
current cases are excluded from identity analyses, object compositions and
trajectory fitting/graphs. Original object classes/counts and native measurements
remain in the audit. `Excluded` is an analysis status, never a biological type.
See the versioned `protein_analysis_exclusion_manifest.json` in the current atlas
typing audit; its keyed overlay is mandatory when using historical native calls.

**What is on disk carries that chain.** `data/development/typing/cells/` was rewritten on 2026-09-19
with the chain applied to all 20 donors — **25,230,109 cells**, `cell_class` eight classes plus
`Excluded`, no `Unresolved`. Measured 2026-09-26: every partition carries `islet_route` (cd99
408,934 · elongation 29,047 · hormone 190,571) and `islet_stratum` (`1` / `2-3` / `4-10` / `>10`),
and the endocrine triple is Beta 235,046 / Alpha 294,496 / Delta 99,010. The row count equals the
2026-09-15 pre-chain build's, so a count cannot tell the two apart; `islet_route` can. On 2026-09-24
the `--refine` join rewrote every partition to add `refined_class` and the pass's columns beside
`cell_class`; the count, the triple and the routes were unchanged by it. On 2026-09-26 donor 6374
alone was retyped (crossed fraction → `subpanel.py` → build → refine → join) on its re-exported
Nerve / Blood vessel / Lymphatic objects: 2,477 of its 2,132,181 cells changed `cell_class`, the
other 19 donors are identical row for row, and the numbers above are that build's. `RECONSTRUCTION.md`
marks every statement `[ON DISK]` or `[CHAIN NOW]`; read the mark before quoting it.

`data/development/typing/subpanel/subpanel_calls.parquet` is the PREVIOUS generation and must not
be read as current. Measured 2026-09-17: 25,150,537 rows, written Sep 11, against the 25,230,109
rows `typing/cells/` then held (written Sep 15) — 79,572 cells apart, and its endocrine counts differ
materially (Beta 191,485 / Alpha 246,055 / Delta 62,797 against 226,183 / 279,686 / 90,658 then, and
235,052 / 294,514 / 99,015 after the 2026-09-19 rewrite, 235,046 / 294,496 / 99,010 since 6374's
2026-09-26 retyping). That tree stacks
generations side by side with nothing marking which is which, so check mtimes before citing anything
in it. `scripts/senior/build_distance_profiles.py` refuses the superseded file by path, with the
counts in the error.

- **How it works:** `core/phenocycler/scripts/typing/RECONSTRUCTION.md` — recovered by reproducing
  the calls from their inputs and verified against them. Read it rather than re-deriving it.
- **Immune subtyping is a second pass over the `Immune` class only**
  (`scripts/typing/immune_subtype.py`, rewritten 2026-09-17, rerun 2026-09-26 over the current
  Immune class — 2,875,793 cells, 20 donors — so `make check` reports `typing/immune` fresh). It
  never modifies `cell_class`. A rule that cannot be asked is EMPTY, not 0; "undecidable" is its own
  outcome, distinct from "no rule fired"; reduced rules are reported on every donor and never
  written into `immune_class`; per-cell `panel_offset` is a diagnostic that gates nothing. **No
  `typing/app` artifact reads `typing/immune`** — the four app builders contain no reference to it,
  and coverage comes from each builder's own `COHORT.resolve`. The `grid_composition` layer does, for
  `fine.csv` (`scripts/senior/build_grid_composition.py:585`).
- **The refinement pass is part of the typing, through the `--refine` join only (ruled
  2026-09-23).** `scripts/typing/refine.py` (concordance / demotion / immune leaf, documented in
  `REFINEMENT.md`) writes a keyed sidecar under `data/development/typing/refine`, and
  `build_cells.py --refine` joins `refined_class` and the pass's per-cell columns **beside**
  `cell_class`; the islet chain's `cell_class`, `is_endocrine` and `endocrine_route` are not
  touched. Its marker set is 54, not 21, and its 5th-percentile bound demotes about 5% of every
  class by construction — not revisited by the ruling. **[ON DISK] the join ran 2026-09-24** on
  all 20 partitions (rerun 2026-09-26 with 6374's retyping), and every layer and app artifact was
  rebuilt from the joined tree on 2026-09-26 (Step 7, `make check` 24 fresh) — but every builder
  still reads `cell_class`, so what they serve is the islet chain's class, not `refined_class`.
  Whether any builder should read `refined_class` instead touches the islet chain, so it is
  asked, never inferred.
- **The ninth class is declared; only the join produces it.** `Unresolved (no lineage evidence)`
  is in `apps/islet/config/run20_default.yml`; `refine.py` writes it into `refined_class` for a
  demoted Exocrine cell, never into `cell_class`. Since the join, 1,063,097 cells carry it in
  `refined_class` and none in `cell_class`.
- Still open, and not to be resolved in prose: which class the builders serve now that the join
  has run, elongation's weight, and the size stratification. `islet_identity_mixture.py` was run for the first time on 2026-09-18 (20 donors,
  outputs under `data/development/typing/islet_identity/`); in the large stratum 14 of 20 donors
  were testable, and CD99 is the strongest feature at 5 of 10 donors firing. That is evidence about
  the question, not an answer to it.
- It is a **development lane**. Nothing here is a release, and the latch stays down.
- Retired, and not to be revived as fallbacks: RESTORE pair webs, hierarchical / simultaneous /
  probabilistic typing, marker-debug and Otsu lanes, the Gate 2–5 apparatus, the broad-typing
  producer comparison withdrawn 2026-09-08, and flat `cells_redsea*` / `restore*` roots. Git has them.
- **Never claim acinar or ductal.** The 60-plex carries no amylase, trypsin, CPA1 or elastase, and
  no CK19, CA19-9, SOX9, HNF1B, CFTR or MUC1. Exocrine is one label.

Authoritative documentation, all relative to `core/phenocycler/`: `README.md`, `DATA_README.md`,
`docs/WORKFLOW.md`, `docs/METHOD_HISTORY.md`; and at repo level `docs/WORKFLOW_2026-09-02_MAP.md`
(what runs today),
`docs/PLAN_2026-09-02_AUXILIARY_PACKAGES.md`, `docs/EXTERNAL_PACKAGES.md`.

## One producer per tissue × modality

| tissue × modality | producer | writes |
|---|---|---|
| pancreas × PhenoCycler | `core/phenocycler/scripts/typing/subpanel.py` → `build_cells.py`; the refinement joins through `build_cells.py --refine` | `data/development/typing/cells` — `cell_class`, and `refined_class` beside it (joined 2026-09-24; builders read `cell_class`) |
| pancreas × Xenium | `analysis/xenium/scripts/panc_xenium_cluster_name.py`, run by `panc_xenium_release.py` | `/data/xenium_typing/panc_2026-09-17/release` — `lineage`, `common_lineage` |
| pLN × Xenium | `analysis/xenium/scripts/pln_immune_typing_v3.py` | `data/development/xenium_typing/release_pln_immune_v3` — `lineage`, `common_lineage` |
| pLN × PhenoCycler | `core/phenocycler/scripts/typing/build_pln_immune_pheno.py` | `data/development/typing/pln/pheno` — `immune_class`, `immune_subtype`, `t_nk_class` |

Everything else that writes a class column is retired, quarantined or a declared library;
`docs/PRODUCERS.md` has each one with its evidence, and
`core/phenocycler/tests/current/test_declared_producers.py` holds the graph to these four. The
same-slide lane (6533, 0112) is typed on its own and is not one of them.

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

### Excluded-cell evidence review, 2026-10-05

Latest human ruling: identify observed patterns, exclude ALL 11,484 reviewed
cells from ALL analyses, and defer identity resolution to a future task.
All 11,459 native
Excluded cells were ordinary-chain fallthroughs inside endocrine outlines; the
reason did not require epithelial positivity and refinement skipped them. Do not
interpret the Excluded label as an exocrine biological class. Versioned evidence
review and exact existing workflow prerequisites are at
`/data/io60-calibration/excluded_cell_identity_review_20261005`.

Current package adds a candidate-evidence table and separate reviewed_identity
for three donor6579 broad immune cases (exact identifiers are in
`source_tables/reviewed_broad_immune_assignments_20261005.csv`). Their observed
broad immune rule, donor-matched marker resemblance, localized CD11c/MPO images
and zero structural overlaps support broad identity only. No fine immune
subtype/state is assigned. Native primary/refined fields and strict analysis
exclusions are retained. All three reviewed cells remain excluded, including
from sensitivity analyses; the old three-cell inclusion sensitivity is superseded.
Every pancreas analytical denominator, comparator and cell-derived aggregate
uses exact-key analysis_eligible masks (31,385,695 retained cells); 31,397,179
assay-QC cells is provenance only. Original object counts/geometry remain metadata.
Future local task: 01a10dd3-8d09-7fd0-90bc-c10290fe4dcf, "Resolve excluded IO60 cell identities".
Candidate families
elsewhere are qualified evidence, not production overrides. Missing negatives,
fragment masks, CD99/shape and structural body rules remain authoritative.

Reference-label resemblance is not independent biological validation. The
reference group holdout does not make queried cells spatially independent.
No expansion rings, replacement distance-band features or new trajectory floor
are introduced.
