# Import provenance

The tracked tree in this directory came from
`https://github.com/smith6jt-cop/SpleenFollicleCounterQP.git`.

| Item | Value |
|---|---|
| upstream fetched tip | `6ae30ba0a528b35e8be1b4abd8edfee84025b2c4` |
| preserved dirty-work commit | `f6cbcd41157c57e718332ee4c262d56551bd07c6` |
| source integration/import commit | `1ed38d2e3b45e59606f219a11dc3b8d3155908c8` |
| source integration branch | `consolidation/import-2026-09-03` |
| preservation branch | `archive/pre-consolidation-worktree-2026-09-03` |
| imported tracked files | 110 |
| imported tracked bytes | 70,963,706 |
| source archive SHA-256 | `51f733f449c076dee554b6275919d112ac1def795180e84867648d445dbe9320` |

The source checkout was first protected without rewriting history: its 15
tracked/untracked worktree entries were committed on the preservation branch.
The upstream tip was then fetched, and that archive commit was cherry-picked
onto a new branch based on upstream `main`. The resulting tree was exported
with `git archive` and placed at `analysis/follicle/`.

The source repository remains intact and clean on its integration branch. Its
large ignored data were never copied into Git. The live source root has the
image and QuPath object trees but no `Measurements/` directory; this monorepo
therefore fails its strict real-data preflight until those exports are restored
through an explicit data-management decision.

Post-import monorepo changes are intentionally not mirrored back to the source
checkout. They establish explicit code/data/output roots, exact manifest-based
identity, root-level CI, and consolidation documentation.
