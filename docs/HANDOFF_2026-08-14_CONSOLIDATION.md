# Handoff — repo consolidation (2026-08-12 → 08-14)

Two threads ran across these sessions: locating the app redesign work, and planning the
seven-repo consolidation. Phase 0 of the consolidation is complete. **Phase 1 as approved is
void** — the upstream scope changed underneath it on 08-13/08-14. That is the first thing to
read below.

Plan: `~/.claude/plans/steady-nibbling-canyon.md`
Execution reference: `docs/CONSOLIDATION.md` (untracked; see § Stale artifacts)
Backups: `/home/smith6jt/repo-backups/{2026-08-13,2026-08-14}/`

---

## 1. The app redesign code — found

Not in this repo. It is `claude/modernist-redesign` in the **Junior** repo
(`~/panc_CODEX/Islet-Explorer`), checked out in a worktree at
`~/panc_CODEX/Islet-Explorer-modernist`, tip `f4e5b22`, 18 commits, 46 files, +5485/−1693.

It was disk-only when found; it is now pushed and matches
`Islet-Explorer/claude/modernist-redesign`. Its own `docs/HANDOFF.md` documents the port,
the open items, and the Viewer's HTTP-Range root cause.

The loose `app_redesign.R` + `Shiny app aesthetic redesign.zip` at the Junior repo root are
the raw claude.ai Design deliverable, preserved on that branch under `design/`. They are the
input, not the work.

---

## 2. Consolidation — Phase 0 complete

Seven repos, disjoint root commits except Senior/Junior (shared root `b78f0fd0`, forked
2026-06-01 at `ef051525`, since then +73 Senior / +34 Junior in **both** directions).

**Done:**

- All seven bundled and `git bundle verify`-clean in `/home/smith6jt/repo-backups/2026-08-13/`
  (1.2 GB). `LymphVizkit` had no local checkout and was mirror-cloned for this.
  `/data/XeniumPanelExplorer` was 11 commits behind its remote — fetched before bundling.
- Working-tree snapshots for the dirty repos (patch + binary-safe tarball).
- `docs/CONSOLIDATION.md` written: target layout, R-lane per-file ownership table, the
  `CLAUDE.md` reconciliation, archive list.
- **Re-snapshot 2026-08-14** after the Phenocycler tree changed substantially:
  `/home/smith6jt/repo-backups/2026-08-14/` — fresh bundle, a 1.9 MB / 124-file patch, and a
  tarball of every modified + untracked file.

**Two plan deviations, both recorded in `CONSOLIDATION.md`:**

- The multi-GB file move was dropped. `CellCompartments.tsv` (11 GB) and
  `Measurements/Cells.csv` (12 GB) are already untracked *and* gitignored, and `git subtree`
  imports history rather than the working tree, so neither can be swept in. Moving `Cells.csv`
  would have broken notebooks H14/H15/H16/H17/H20 to solve a problem that does not exist.
- `Phenocycler_Analysis` gets full history via `git subtree` (3.4 MB) rather than a squash —
  it is the scientific record. Everything else still squashes; Xenium alone would otherwise
  pin 777 MB of committed figures into the monorepo permanently.

---

## 3. What changed underneath the plan — read before doing anything

`CLAUDE.md` was rewritten on 08-13/08-14, and `Phenocycler_Analysis` now carries **136
uncommitted files** (75 deleted, 49 modified, 12 untracked). Together they retire the entire
marker-calling and typing lane:

**The workflow now stops after expression.**

```text
QuPath cohort manifest -> ingest -> geometry QC -> REDSEA -> one intensity per marker
```

`controls`, `calibrate`, `type`, state calling and assignment export are **blocked**. The
run command is now `--pipelined --through expression`.

**Deleted in the working tree** (11 package modules, 48 test files, 8 notebooks,
`refinement.py`, 5 review scripts):
`hierarchical_typing.py`, `marker_calibration.py`, `probabilistic_training.py`,
`threshold_selection.py`, `typing_model_bundle.py`, `typing_stage.py`, `typing_rules.json`,
`candidate_evaluation.py`, `evaluation_reporting.py`, `refinement_contracts.py`,
`restore_expression_reviews.csv`. Notebooks 03–09 are gone; only
`01_qupath_ingest_and_geometry_qc.ipynb` and `02_redsea_spillover.ipynb` remain.

**New, untracked:** `empirical_tail.py`, `restore_seeded.py`,
`fundamental_pair_diagnostics.json` (now named in `CLAUDE.md` as an authority), the
`fundamental_pair_*` / `reference_overlap_audit` / `boundary_contact_audit` review scripts,
and a restructured `tests/{current,integration,research}/`.

**Rule changes that matter most:**

- *"There are no verified production target/reference pairs. Do not infer one from historical
  registries, outputs, or expected cell frequencies."*
- RESTORE pair webs, hierarchical/simultaneous/probabilistic typing, marker-debug/oracle and
  Otsu lanes are all now **retired**, not production.
- The diagnostic roster is only raw CD3e/Ker8_18 and signed Iba1_norm/raw Ker8_18. CD79a is
  excluded — too sparse in normal tissue. These are rank/overlap comparisons, **not**
  production pairs or calls.
- `data/pilots/` is superseded by `Phenocycler_Analysis/../data/development/`.

### Consequence: Phase 1 of the approved plan is void

Phase 1 was the source→app bridge:

```text
<run>/cells + expression + marker_evidence + assignments + cell_states -> H5AD -> app aggregates
```

Three of those five inputs no longer exist. `CLAUDE.md` now states the adapter is for a
*future verified run* and that until a verified pair/typing workflow and its adapter exist,
the legacy H5AD rebuild must not be described as current and provisional labels must not be
fed to the app scripts.

**Do not build the bridge.** The consolidation phases that remain viable are 0 (done), 2
(found the monorepo shell), 3 (absorb analysis lanes), 5 (unify the R app lane) and 6
(cut over and archive). Phase 4 needs rewriting — see below.

---

## 4. Open decisions and blockers

1. ~~**The 136-file retirement is uncommitted.**~~ **RESOLVED 2026-08-14.** The real count was
   **172** (75 deleted / 49 modified / **48** untracked). Landed as four themed commits on
   `Phenocycler_Analysis` branch `agent/retire-typing-lane` (`8337e3a..f6953cf`) and pushed:
   retire the lane → land the review lane → split the suite → document the workflow. Verified
   before pushing: 287 current + 18 research + 214 integration passed (18 integration skips,
   all the uninitialized XeniumPanelExplorer submodule, which CI checks out recursively); the
   ruff E9/F63/F7/F82 gate clean; `DEFAULT_STAGE_ORDER == (ingest, geometry, redsea,
   expression)`; all 28 live modules import; `status` still prints `BLOCKED NO_VERIFIED_PAIRS`
   for controls/calibrate/type/qupath; `git diff --check` clean.
   **`git subtree` reads `main`, so this branch must be merged before Phase 2 step C.**
   Two things were fixed as a precondition: `notebooks/02_redsea_spillover.ipynb` was stripped
   of its 9 executed output cells (256 KB → 18 KB, incl. a 141 KB embedded PNG and a leaked
   local path) and its hardcoded absolute `repo_root` replaced with `base_cfg.runs_dir`, which
   would otherwise have broken on relocation.
2. **Phase 4 (resolve duplicated science) needs rewriting.** It was built on comparing
   `SpleenFollicleCounterQP`'s notebook-resident RESTORE reimplementation against production
   `restore_stage.py`. There is no production caller now. The *finding* stands and arguably
   sharpens: `analysis/follicle` H16 cell 5 runs a 168-line GMM RESTORE whose approach
   `METHOD_HISTORY.md` records as tested and rejected, and H20 calibrates CD20 thresholds
   against H16's labels. What it should be compared *to* is now an open question.
3. **Junior's fate** — still unanswered from the planning session. Production
   (`/srv/shiny-server/islet-explorer`) still symlinks Junior's `app/shiny_app`. The plan
   assumes Junior serves until Phase 6 and its `R/` work is harvested in Phase 5.
4. **Project name.** `islet-explorer` no longer describes a scope covering spleen, lymph node,
   thymus and Xenium.

---

## 5. Stale artifacts to fix

- **`docs/CONSOLIDATION.md`** (untracked) was written 08-13 against the pre-retirement world.
  Its rule-hierarchy row *"the 8-stage workflow, RESTORE's status"* and its blocker count of
  53 uncommitted files are both wrong. Corrected in place on 08-14; re-read before executing
  Phase 2.
- **`~/.claude/plans/steady-nibbling-canyon.md`** still contains the void Phase 1 and the
  RESTORE-based Phase 4. Kept as-is for the record; this handoff supersedes it on both points.
- The memory `phenocycler-method-current.md` described RESTORE as the `calibrate` producer and
  simultaneous 5-gate typing as current. Updated 08-14.

---

## 6. State on disk

*2026-09-02: paths have moved since — Junior is at `/data/panc_CODEX/Islet-Explorer`,
`SpleenFollicleCounterQP` at `/data/SpleenFollicleCounterQP`; see `docs/CONSOLIDATION.md` § Status 2026-09-02.*

| Repo | Branch | Dirty | Note |
|---|---|---|---|
| `Islet-Explorer-Senior` | `agent/land-retirement-and-consolidation-docs` | — | branched off `main` @ `e190d42` |
| `Phenocycler_Analysis` | `agent/retire-typing-lane` @ `f6953cf` | clean | retirement landed + pushed; **merge to `main` before Phase 2** |
| `Islet-Explorer` (Junior) | `claude/spatial-cardc…` @ `de81dce` | clean | serves production; remote is named `Islet-Explorer`, not `origin` |
| `…-modernist` worktree | `claude/modernist-redesign` @ `f4e5b22` | clean | pushed |
| `Xenium_Analysis` | `main` @ `8062e55` | clean | at `/data/Xenium_Analysis`; 12 files hardcode `/blue/maigan/...` |
| `SpleenFollicleCounterQP` | `main` @ `2040bb1` | 3 binaries + 1 untracked | at `/home/smith6jt/...`; 68 deleted source files **restored from HEAD 08-14** |
| `XeniumPanelExplorer` | `main` @ `0a23564` | clean | local checkout is **11 behind** `origin/main`; `origin/main` == Phenocycler's gitlink `076d1192` |
| `LymphVizkit` | — | — | GitHub only; holds `00_domain_config.R` |
