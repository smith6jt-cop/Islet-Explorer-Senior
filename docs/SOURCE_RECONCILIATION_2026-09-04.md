# Source repository reconciliation — 2026-09-04

This is the commit-by-commit disposition ledger for plan Track D.4 and consolidation
Phase 6. It records what was carried into this monorepo, what was deliberately adapted or
superseded, and what cannot happen until the final identity release exists. It does not make
a scientific release out of development artifacts.

## Result

- Source-code reconciliation is complete for all six repositories in the archive list.
- The two source-audit remotes were removed after their URLs and exact tips were recorded
  below. The monorepo now has one local branch (`main`), one origin branch (`origin/main`),
  and one remote (`origin`).
- Archiving the six GitHub repositories is **not** complete. Track D.4 explicitly put that
  action after D.1–D.3; D.3 remains gated on a declared A.4b identity release and the
  downstream application declaration. Junior is also the publication build. No repository,
  production symlink, or service was changed during this audit.

GitHub was queried read-only on 2026-09-04. All six repositories were still writable
(`archived=false`) and their default branch was `main`.

| source | recorded URL | published `main` | monorepo owner |
|---|---|---|---|
| Islet-Explorer | `https://github.com/smith6jt-cop/Islet-Explorer.git` | `f2a5c540ab95584687b19683f2c23b8bd952beaa` | selected shared fixes in `apps/islet`; publication app remains separate |
| Phenocycler_Analysis | `https://github.com/smith6jt-cop/Phenocycler_Analysis.git` | `fc06c0f7a394c4c527352f61dce885cdd33de786` | `core/phenocycler` |
| Xenium_Analysis | `https://github.com/smith6jt-cop/Xenium_Analysis.git` | `8062e55407e57146cc658db3a88aa9a3c7e923fc` | `analysis/xenium` |
| XeniumPanelExplorer | `https://github.com/smith6jt-cop/XeniumPanelExplorer.git` | `076d1192a25946ac3f9cf941c98202d1f224579a` | `core/panels` |
| SpleenFollicleCounterQP | `https://github.com/smith6jt-cop/SpleenFollicleCounterQP.git` | `6ae30ba0a528b35e8be1b4abd8edfee84025b2c4` | `analysis/follicle` |
| LymphVizkit | `https://github.com/smith6jt-cop/LymphVizkit.git` | `820bb382280feb9f4af0231379943afd3ca211d0` | `apps/R/ieCore` plus `apps/lymph` configuration/data |

## Phenocycler_Analysis

The subtree was imported from `e9874a9ffa6384857fc22179db0984bdf570c4d4`.
The published source tip has exactly two unique non-merge commits:

| source commit | disposition |
|---|---|
| `33e4da8` — section-keyed cohort | Ported and adapted. Artifact discovery now keys sections by `(donor_id, roi)`, exposes `section_keys`, retains `donors()` for donor-only callers, rejects duplicate section keys, and validates manifest sections before any QuPath work. A donor can represent `panc`, `pln_1`, and `pln_2` without collision. |
| `ffba0a8` — narrow parse guards/template | Ported. Configuration/cohort parsing catches `SectionError` rather than hiding unrelated failures, duplicate messages name the section key, and the three-section HiPerGator example is retained. |
| `da5f18c`, `fc06c0f` | Merge commits only; their two content commits are accounted for above. The other merge parent is the already-imported `e9874a9` line. |

The port changes the same seven source paths but is not copied blindly: current monorepo
configuration and documentation have moved since the source fork. The 53 focused section-key
tests pass, as do the complete current/integration/research and re-stain suites listed below.

## XeniumPanelExplorer

Blob comparison against `076d1192` found 285 source files and 283 monorepo files:

- 282 blobs are byte-identical.
- `.github/workflows/lint.yaml` and `.github/workflows/tests.yaml` are the only omitted source
  files. Nested workflows do not run on GitHub; the repository-root Phenocycler workflow owns
  the `core/panels/**` path trigger and taxonomy gate. The source R suite remains a listed CI
  gap and was run locally during this reconciliation.
- `README.md` is the only changed blob. Its added monorepo note accurately records that only
  the R 4.6 backup lock/profile exist and that the host R 4.5 library is the current runtime.
- There are no monorepo-only code or data blobs. In particular, panel roles, gene roles,
  identity cores, and tissue subpanels remain byte-identical to the published source.

The local source checkout at `0a23564` is merely 11 commits behind its published `main`; it
was not treated as authority over `076d1192`.

## Xenium_Analysis

The local source checkout is clean at published `main` `8062e554`. A blob audit found 108
source files versus 90 files under `analysis/xenium`: 42 identical, 34 intentionally omitted,
16 monorepo path/provenance additions, and 32 adapted files.

The detailed per-path rationale already lives in `analysis/xenium/README.md` under “What was
dropped on absorption, and why.” The 34 omissions are source-repository scaffolding or
superseded material: nested CI/setup files, local Claude state, the contradictory second
`config.ini`, host-specific paths/environment recipes, retired embedding-gene panels,
`program_panels.py`, a demo R script, an unrun integration notebook, and external data files.
Fourteen run/environment scripts moved under `analysis/xenium/hpc`; `HANDOFF.md` became the
dated `docs/HANDOFF_2026-07-13.md`; `scripts/workflow_paths.py` and its test enforce the
monorepo layout. The 32 changed files are path adapters, documentation, or post-absorption
workflow fixes already owned and tested here. No unique source functionality remains stranded.

## SpleenFollicleCounterQP and LymphVizkit

Spleen source history is preserved without pretending the local work was on published
`main`: upstream is `6ae30ba`, the 15-entry worktree archive is `f6cbcd4`, and the clean
integration source tree is `1ed38d2`. Its 110 tracked files, strict 139-alias identity
manifest, adaptations, and test evidence are recorded in
`analysis/follicle/IMPORT_PROVENANCE.md`. Missing external measurement CSVs remain a data
preflight failure, not missing source code.

LymphVizkit `820bb38` was verified from the durable backup bundle. Its domain YAML and the
imported YAML are byte-identical; the shared loader became `apps/R/ieCore`, and
`apps/lymph` contains configuration plus a small synthetic schema fixture rather than a
second Shiny app. Exact hashes are in `apps/R/ieCore/IMPORT_PROVENANCE.md`.

## Islet-Explorer (Junior)

Junior publication `main` remains `f2a5c54`. The clean local review branch has one later
commit, `e589961` (“Harden UFRC reviewer deployment workflow”). Its portable changes are now
adapted into Senior:

- fixed-schema, allowlisted validation events that never log a rejected browser value;
- exact choice, bounded-number, asset-basename, and path-segment validation at browser-to-
  filesystem/query boundaries;
- production refusal of debug credential/viewer modes and of network AI credential use;
- generic user-facing failures instead of paths, endpoints, provider errors, or identifiers;
- sanitized Shiny errors, production launch flags, skip/focus semantics, landmarks, and
  higher-contrast sidebar buttons;
- a Senior-specific security boundary test suite.

Junior's OCI reviewer-release pipeline, immutable container evidence scripts, handoff text,
and request-form binary were not copied. They depend on Junior's release/data-backend and
container deployment architecture, while Senior is currently launched directly by systemd
and has no D.3 bridge release. Copying that pipeline would create a second, nonfunctional
deployment authority. Its provenance/release intent belongs in D.3/D.7 after the identity
release exists. Junior-only Barlow/PAGA helpers are absent here; equivalent redaction and
validation were applied to Senior's current trajectory implementation instead.

## Validation evidence

- Fatal Ruff (`E9,F63,F7,F82`): pass.
- Phenocycler: 322 current, 303 integration, and 173 combined research/re-stain tests pass.
- Xenium analysis: all 85 tests pass. The zarr export tail deadlocks at `zarr.open_group`
  under the filesystem sandbox but completes normally outside it; the bounded diagnostic
  located the wait in zarr's sync wrapper rather than repository code.
- Islet R lane: all 188 assertions pass; all package/app sources parse; a production-flagged
  app source smoke returns a server function.
- `core/panels`: the full local testthat suite passes. Two explicit environment skips remain:
  the Louvain fallback is not applicable because `leidenAlg` is installed, and Chrome cannot
  start in this sandbox for the `shinytest2` browser smoke.

The requested donor checks used the retained, immutable development artifacts rather than
rerunning two expensive scANVI models:

| donor | RNA base / sidecar | PhenoCycler / Xenium geometry | registration |
|---|---|---|---|
| 6436 / `0059907` | 641,487 rows in each; 5,101 genes / 0 variables; recorded SHA-256 pair matches | 1,579,440 / 695,113 cells | Dice 0.5855966, correctly `FAIL` at 0.80 |
| 6442 / `0059947` | 906,211 rows in each; 5,101 genes / 0 variables; recorded SHA-256 pair matches | 1,419,669 / 958,744 cells | Dice 0.7676936, correctly `FAIL` at 0.80 |

Both sidecars retain a 10-dimensional `X_scvi_broad` embedding. The build summary hash also
matches its recorded value. These are development typing sidecars, not declared integration
identity releases. The integration status therefore continues to block identity-consuming
stages even though geometry, imports, registration, and QC exist.

## Remaining gate

The final archive/pointer operation, `core/bridge`, first `phenocycler-v*` tag, fresh-clone
donor-to-app proof, and any production change remain out of scope until all of the following
exist:

1. A.3 selects the PhenoCycler producer and A.4b supplies a reviewed identity release.
2. The maintainer writes the downstream application declaration required by root `CLAUDE.md`.
3. The maintainer chooses the first package version/project name and separately authorizes
   any production cutover.

Until then, `Ambiguous`/`Unavailable` remain first-class, cached identities are ignored, all
six source repositories remain unarchived, and the publication service remains untouched.
