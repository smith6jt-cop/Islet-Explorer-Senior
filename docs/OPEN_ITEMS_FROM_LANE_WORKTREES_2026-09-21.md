# Open items recovered from the lane worktrees, 2026-09-21

Seventeen `wt-<letter>-<topic>` worktrees were removed on 2026-09-21 (Step 2 of the
cleanup plan). Their branches are kept. Before removal, the 26 `BRIEF.md` / `FIX.md` /
`NEEDS.md` files they carried were archived to
`/home/smith6jt/repo-backups/worktree-notes-2026-09-21/`.

Five of them were `NEEDS.md` — a worker's record of what it could **not** finish and why.
Those records were the only place this information existed. They are distilled here so the
open items survive the folders that held them.

Nothing below has been acted on. Each needs an owner decision, not a code change first.

## 1. A card exists that the app cannot display

`s1_size_composition` is listed in neither `IE_CROSSMODAL_CARDS`
(`apps/islet/R/mod_crossmodal_ui.R`) nor `IE_LYMPHNODE_CARDS`
(`apps/islet/R/mod_lymphnode_ui.R`), and no other module calls its UI or server. The card
is written and tested; it has no mount.

It now defaults to `tab_key = "Plot"`, which is a **proposed** destination, not an existing
one. Mounting means calling `card_s1_size_composition_ui(ns)` and
`card_s1_size_composition_server(...)` together on the chosen tab, forwarding
`active_tab = reactive(input$tabs)`. If another tab is preferred, its exact key is passed as
`tab_key`.

Reported independently by two workers (`wt-I-islet`, `wt-M-cards`), both of which recorded
that registry and app wiring were outside their file ownership.

## 2. The app's own startup assertion blocks its smoke test

`apps/islet/app.R:304` asserts `IE_TAB_KEYS` equals the old tab list — without `CrossModal`
and `LymphNode`. Those keys already exist in `apps/islet/R/00_globals.R`. The browser smoke
test therefore fails **during startup**, before any card renders:

```sh
cd apps/islet && NOT_CRAN=true ISLET_DOMAIN_CONFIG=config/run20_default.yml \
  /opt/R/4.5.1/bin/Rscript --vanilla tests/tab_smoke.R /tmp/tab-smoke
```

Separately, `s3_compartments` is omitted from both card registries though its UI and server
exist in `apps/islet/R/card_s3_compartments.R`.

Recorded by `wt-G-gradient`.

## 3. A standing-rule question: Xenium records carry no `object_id`

`wt-H-restain` stopped rather than guess, and the question is still open.

The repository rule is that the QuPath UUID `object_id` is preserved throughout and never
substituted by centroid matching, row order, filenames or donor-name parsing. **The Xenium
sources do not contain that key at all.** Measured in that worktree:

| donor | Xenium cell rows | warp rows | recorded nucleus IDs present |
|---|---:|---:|---:|
| 6436 | 695,113 | 695,113 | 668,284 |
| 6442 | 958,744 | 958,744 | 907,352 |

Both `cells_xen` tables carry native Xenium `cell_id` and `feat__if_INS` / `feat__if_GCG`
but no `object_id`; the vendor bundles' `cells.parquet` carry none either; and the warp
sidecars are `xenium_cell_ids.npy`, native Xenium ids rather than QuPath UUIDs.

Two ways forward, and the choice is the maintainer's:

1. supply a producer-recorded UUID crosswalk for these exact cells and warp records; or
2. explicitly permit aligning these Xenium-only records by native `cell_id`, with
   `object_id` recorded as `unavailable(no_qupath_object_id)`.

No centroid matching, fabricated UUID or positional fallback was added, and the builder's
acceptance test fails at this contract before any artifact is written.

## 4. Derived inputs that have not been produced for part of the cohort

From `wt-S-artifacts`. All five builders succeeded; these are upstream gaps, recorded as
*unavailable* rather than as zero.

- **112 and 6533** — raw PhenoCycler islet exports exist, but `typing/cells`,
  `typing/islet_membership` and `typing/app/donor_metadata.csv` were never extended to them.
- **UUID membership parquets never built** for seven donors, with the measured counts:

  | donor | membership UUIDs unavailable / eligible cells |
  |---|---|
  | 6450 | 364,379 of 364,379 |
  | 6505 | 2,145,798 of 2,145,798 |
  | 6523 | 809,328 of 809,328 |
  | 6534 | 1,115,452 of 1,115,452 |
  | 6566 | 1,779,170 of 1,779,170 |
  | 6591 | 1,099,734 of 1,099,734 |
  | 6623 | 1,189,718 of 1,189,718 |

- **Islet summary and trajectory sources** cover only 115, 6374, 6436, 6442, 6461, 6476,
  6482, 6493, 6539, 6593. The twelve uncovered cohort donors are 112, 6380, 6450, 6505,
  6521, 6523, 6533, 6534, 6551, 6566, 6591, 6623 — all have raw islet annotations; the
  derived aggregates, features and order need upstream production. For 6380, 6521 and 6551
  `cohort.py` additionally records marginal missed-islet deficits, and the layer builder
  must not override that provenance.

## 5. One commit that never landed

`wt-G-gradient` could not commit: `git add` failed with `Read-only file system` creating
`.git/worktrees/wt-G-gradient/index.lock`. Its work reached `work/G-gradient` by another
route (that branch is fully merged into `main`), but the intended message is recorded here
because the measurement in it is not otherwise written down:
*"Measure 1,808,128 immune cells across 13 donor distance profiles"* — 16,532,662 eligible
cells, 1,808,128 immune cells, 13 donors; seven donors had no matching membership UUIDs and
were explicitly unavailable.
