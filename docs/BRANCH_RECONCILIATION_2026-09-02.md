# Branch reconciliation ledger — 2026-09-02

> **DATED LEDGER — two facts it records have moved, 2026-09-18.** The branch dispositions and the
> safety record are provenance and stand as written. The workflow facts quoted in them do not:
> the **content-addressed run framework was removed on 2026-09-16** (no `run_id`, no stage
> manifest, no fingerprint validation; stages write to stable roots), so "the current immutable-run
> workflow" under § Conflict decisions and "immutable run `30f6689e4d27d295f890`" under § Merged-tree
> validation describe a framework that no longer exists. And the **25,150,537 analysis-eligible
> cells** recorded there is that run's geometry export, which is not the current cell-type
> population: `data/development/typing/cells/` holds **25,230,109 cells over 20 donors**, written
> 2026-09-15. The same 25,150,537 is also the row count of
> `typing/subpanel/subpanel_calls.parquet` — the **superseded** previous generation of the calls,
> 79,572 cells short of the current tree. Neither is the number to quote for the typing.


This ledger records the evidence and decisions used to converge the repository onto one
persistent branch. It is deliberately separate from the workflow map: branch history is
provenance, not workflow authority.

**Status: complete 2026-09-03.** The sibling merge is `0fc21e7`; the final merge on
`main` is `6cf301d6e8ee48470a50d34a73ce30c257e7d945`. Both retired local/origin refs were
deleted only after ancestry and empty-cherry proofs. `main`/`origin/main` are the sole
active monorepo branches. The temporary `panels` and `phenocycler` source-audit remotes were
retained through Track D.4a, then removed after their URLs, tips, unique commits, and blob
dispositions were recorded in `SOURCE_RECONCILIATION_2026-09-04.md`.

## Safety record

- The seven source-repository bundles under
  `/home/smith6jt/repo-backups/2026-08-13/` were verified with
  `git bundle verify` before reconciliation.
- A fresh bundle of this combined repository was created before the merge at
  `/tmp/Islet-Explorer-Senior-pre-reconcile-224ae9d.bundle`.
- The fresh bundle verified as complete and contains every local/remote branch, tag, and
  session checkpoint present at that point.
- SHA-256:
  `7af03237596fc1a62f72b82f9cce01b3bb4783ddfb92a2106c9625df3fa89d42`.

The `/tmp` bundle is a short-lived safety copy. The dated source-repository bundles remain
the durable pre-consolidation archives.

## Starting topology

| Ref | Starting tip | Disposition |
|---|---|---|
| `main` | `412603d` | Final persistent branch after validation and merge |
| `feat/dct-step0-pinning` | `224ae9d` | Working branch; includes the current workflow and continuation plan |
| `chore/remove-gold-standard-lane` | `b9be3bb` | Merged into the working branch, then retired after validation |
| `origin/feat/dct-step0-pinning` | 10 commits behind the feature tip before the plan commit | Retired after `main` is published |

At the start of reconciliation, the feature branch was 63 commits ahead of `main`; the
chore branch was a sibling containing three commits not present on the feature branch.

## Unique chore-branch commits

| Commit | Decision | Evidence |
|---|---|---|
| `11d80c3` — Remove gold-standard lane | Merge | No live code imported or invoked `scripts/gold_standard/`; the six retired scripts were deleted. |
| `be9406c` — dated 08-16 audit docs | Merge as dated history | Retained the fundamental-pair review, handoff, app-test guidance, and corrected Groovy path. Current 2026-09-02 workflow facts win wherever dated text conflicts. |
| `b9be3bb` — split integration geometry and identity gates | Merge and validate | Geometry export no longer requires an identity assignment; typing-dependent stages remain blocked by `NO_VERIFIED_TYPING`. Targeted integration tests passed: 53 tests. |

## Conflict decisions

- `README.md`: retain the current immutable-run workflow and combined-repository layout;
  retain the editable-install and legacy-app boundary clarifications from the chore branch.
- `scripts/CLAUDE.md`: retain the current `.venv`/HPG execution contract,
  `NO_VERIFIED_TYPING`, and the retirement of RESTORE-derived assignment export.
- `deploy/DEPLOY.md`: the current 2026-09-02 observation that the publication symlink is
  dangling overrides the dated 2026-08-20 deployment snapshot.
- `docs/CONSOLIDATION.md`: the 2026-09-02 correction and continuation plan are current
  authority; the older handoff remains explicitly labeled as a historical snapshot.

## Merged-tree validation — 2026-09-03

- The adapted C.1 boundary passed 635 current, research, and integration tests with no skips
  after additive installation of SimpleITK 2.5.6 and spatialdata-io 0.7.1. Blocking Ruff rules
  (`E9`, `F63`, `F7`, `F82`), `pip check`, and `git diff --check` pass. The repository's full
  advisory Ruff profile still has pre-existing debt and is not an R.2 gate.
- PhenoCycler geometry was regenerated from immutable run `30f6689e4d27d295f890`: 20 pancreas
  partitions, 25,150,537 analysis-eligible cells, schema v1 geometry profile, no excluded donors,
  no hormone identity features, and no duplicate UUIDs. Donors 6436 and 6442 contain 1,579,440 and
  1,419,669 cells, respectively.
- Xenium geometry imports for 6436 and 6442 contain 695,113 and 958,744 cells, respectively,
  with duplicate-free native `cell_id`, blank identity, and unavailable (nullable) evidence.
- Raster/image affine registration produced tissue Dice 0.586 for 6436 and 0.768 for 6442;
  both correctly grade `FAIL` at 0.80. A real 6436 B-spline probe was stopped after extended
  non-completion and is recorded only as performance evidence. The last completed affine
  artifacts remain authoritative.
- Direct selection of `structures` for donor 6436 fails in preflight because no declared
  PhenoCycler identity release is configured, even though cached structure outputs exist.
- The Xenium B.1 seam is closed and exercised end to end. Donor 6436 (`0059907`) produced a
  641,487-row slim sidecar from 416,962 clean anchors plus 224,525 propagated calls; donor 6442
  (`0059947`) produced 906,211 rows from 755,067 anchors plus 151,144 propagated calls. Required
  label/evidence fields have no nulls, IDs are unique and exactly ordered with their base files,
  and low-confidence calls remain visible (25.47% and 22.08% below 0.50). The full 85-test Xenium
  suite passes. No integration release was declared from these development artifacts.

## Exit conditions — met 2026-09-03

The completed reconciliation satisfies all five conditions:

1. the merged tree passes the relevant test suites and donor checks for 6436 and 6442;
2. `main` contains both former branch tips;
3. the reconciled `main` is published successfully;
4. the obsolete local and remote feature/chore branches are deleted; and
5. `git branch -a` shows one persistent development branch (`main`); after D.4a removes the
   source-audit remotes, the only remaining remote branch is `origin/main`.
