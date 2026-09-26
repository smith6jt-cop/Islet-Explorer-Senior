# PhenoCycler cell typing, and the app dataset built from it

Two things live here. The **typing** is the subpanel classification — context from QuPath objects,
then one question per context against a reference level, never a constant. **A structure's level is
the p99 of that donor's ordinary tissue — cells in no classified object and outside every islet —
fitted PER DONOR.** Inside an islet the question is the same statement upside down, so **CD99's
endocrine level is the 1st percentile of that donor's AND that size stratum's own hormone-confirmed
cells**, with a floor guard that returns the level UNAVAILABLE rather than accept one sitting on the
scale floor. The donor's ordinary-tissue CD99 level is carried and reported beside it as
`donor_reference_level`, never substituted for it. The shape criterion is fitted the same way and
only in the `>10` stratum, where alone it separates — and its weight is an open question (see
RECONSTRUCTION.md, "Elongation is not a level-able feature"). The **app build** turns those calls
into the dataset `apps/islet` serves under `ISLET_DOMAIN_CONFIG=config/run20_default.yml`.

Development lane: writes only `data/development/typing/`, never the stable
`data/{cells,geometry_qc,redsea,expression}` roots or the historical
`data/runs/` tree, and creates no manifest of any kind (there is no stage
manifest left to write); does not move the `NO_VERIFIED_TYPING` latch.

## The typing IS re-derived here

This section said the opposite until 2026-09-18, and it was true when it was written: the 25.15M-row
table had no producer, so nothing could re-decide `cell_class`. Two modules changed that, and both
are development lanes with `NO_VERIFIED_TYPING` still down.

```
celltype/uniform/  ->  subpanel.py     the classification  ->  typing/subpanel_rebuilt/
                   ->  build_cells.py  carries + derives   ->  typing/cells/
                   ->  refine.py       the refinement pass ->  a SIDECAR at --out, keyed on object_id
```

`subpanel.py` is the producer the table never had. `build_cells.py` reads
`typing/subpanel_rebuilt/` (`subpanel_io.TYPING`) and **carries** — every column the producer
writes reaches the output, `is_endocrine` is read off `resolved_by` rather than re-derived from the
class name, `Excluded` keeps its reason, and `--audit` attaches the per-stratum level each cell was
measured against. It decides nothing, but what it carries is whatever the last run of the producer
wrote, so **the producer runs before the app build whenever the typing changes**. `refine.py` does
not rewrite `typing/cells/`; it writes a sidecar beside it, joined only under
`build_cells.py --refine`.

### The 2026-09-18 islet chain is on disk, for every donor

Measured 2026-09-26. `typing/subpanel_rebuilt/` and `typing/cells/` both hold 25,230,109 rows over
20 donors and were rewritten on **2026-09-19** (09:41–09:46 and 09:55–09:59) with the islet chain,
and 6374 alone was retyped 2026-09-26 on its re-exported structure objects:
every partition carries `islet_route` (cd99 408,934 · elongation 29,047 · hormone 190,571) and
`islet_stratum` (`1` / `2-3` / `4-10` / `>10`). `cell_class` has nine values — the eight classes
plus `Excluded` — and **no `Unresolved`**: the refinement pass was ruled in on 2026-09-23 and is
wired as the `build_cells.py --refine` join, which has not run, so `typing/cells` carries no
`refined_class` yet. The row count equals the 2026-09-15 pre-chain build's; read `islet_route`,
not a count, to tell the two apart.

Read, rather than re-deriving them:

* [`RECONSTRUCTION.md`](RECONSTRUCTION.md) — the classification: the input chain and context
  reproduction (recovered by reproducing the calls from their inputs and checked against them), the
  islet chain ruled on 2026-09-18, the measured evidence for the size strata, and three defects the
  build found in the calls on disk.
* [`REFINEMENT.md`](REFINEMENT.md) — the refinement pass, recovered from the 2026-09-11 session:
  concordance, the per-donor bound, `immune_leaf`, and why it is applied once. **Wired, not run:
  the `build_cells.py --refine` join adds `refined_class` beside `cell_class`; no builder reads it
  and no cell carries `Unresolved` today.**

The two modules carrying the parts verified exactly are `subpanel_io.py` (input chain, scale,
reference levels) and `context.py` (centroid-in-object assignment, reproduced at 100.0000%).

**Two generations of calls sit in this tree and nothing on disk marks which is which.**
`typing/cells/` is current (25,230,109 rows); `typing/subpanel/subpanel_calls.parquet` is the
previous generation (25,150,537 rows) and its endocrine counts differ materially. Check mtimes
before citing anything here.

## The app build, in order

Run from the repo root with the repo `.venv`. **Seventeen** clean donors by default, twenty with
`--include-marginal` (`cohort.resolve`, measured 2026-09-18). It was ten batch-1 donors until the
2026-09-15 re-export: the seven re-exported batch-2 donors joined `CLEAN` then, and `resolve`
refuses a batch-2 donor only when no rebuild exists for it under `typing/revised_cells/`.
The three MARGINAL donors (6521, 6551, 6380) carry a stated islet deficit and are opt-in.

```bash
.venv/bin/python core/phenocycler/scripts/typing/build_cells.py                 # per-cell table (first build after new calls)
.venv/bin/python core/phenocycler/scripts/typing/refine.py --out data/development/typing/refine   # the refinement, ~10 min
.venv/bin/python core/phenocycler/scripts/typing/build_cells.py --refine data/development/typing/refine   # joined beside cell_class
.venv/bin/python core/phenocycler/scripts/typing/build_nonendocrine_outline.py  # ~30 s/donor
.venv/bin/python core/phenocycler/scripts/typing/build_islet_geometry.py        # four regions
.venv/bin/python core/phenocycler/scripts/typing/aggregate_islets.py            # per-islet counts
.venv/bin/python core/phenocycler/scripts/typing/build_galectin_layer.py        # correction J
.venv/bin/python core/phenocycler/scripts/typing/neighborhood.py                # Spatial A/B/C
.venv/bin/python core/phenocycler/scripts/typing/build_workbook.py              # the two sheets
.venv/bin/python core/phenocycler/scripts/typing/build_app_data.py              # tissue, cells, outlines, distance
.venv/bin/python core/phenocycler/scripts/typing/build_crossed_fraction.py      # correction D, ~2 min/donor
.venv/bin/python core/phenocycler/scripts/typing/build_phenotype_rules.py       # the rules modal
.venv/bin/python core/phenocycler/scripts/typing/build_parquet_views.py         # the DuckDB views
.venv/bin/python core/phenocycler/scripts/typing/summarise_crossed_fraction.py  # correction D -> islets
```

The refinement is part of the typing since the maintainer's ruling of 2026-09-23, through the
`--refine` join only: `refine.py` reads `typing/cells`, so the plain `build_cells.py` run is
needed first only when `subpanel.py` has written new calls, and the join ADDS `refined_class`
and the pass's per-cell columns beside `cell_class` without changing the islet chain's
`cell_class`, `is_endocrine` or `endocrine_route`. Pass `--include-marginal` to both
`build_cells.py` runs for the twenty-donor app cohort, in one invocation each, because the
build ledger is rewritten whole by every run.

Order matters in three places: the outlines must exist before the geometry (the non-endocrine
region has no shape without them), the geometry before the workbook, and
`build_parquet_views.py` before the two `summarise_*` scripts, which merge into the
`islets.parquet` it writes. Re-running `build_parquet_views.py` drops their columns, so the
summarisers run after it, never before.

### Viewer

`.tmap` region layers must be built from **`data/development/typing/app/json`**, not the
builder's default `data/Islets`. The QuPath export carries only `Islet` and `IsletExpanded`; the
merged files add `IsletNonEndocrine`, and a project built from the default ships without it:

```bash
.venv/bin/python scripts/senior/build_all_tissuumaps_projects.py \
  --cases 115 6374 6436 6442 6461 6476 6482 6493 6539 6593 \
  --regions data/development/typing/app/json
```

Channels are cached per donor under `/data/tissuumaps/channels/<case>`, so a restart costs
minutes rather than hours -- check the region source before letting a ten-donor run finish.

**Do not run this while QuPath is open.** One channel of a 128 GB slide held 175 GB RSS; an
interactive QuPath session was holding 332 GB, and the box has 503 GB. Check `free -g` (the
*available* column) first. The build is safe to interrupt -- it resumes from the channel cache.

### Trajectories

```bash
.venv/bin/python core/phenocycler/scripts/typing/trajectory_islet.py            # islet-level, ~3 min
.venv/bin/python core/phenocycler/scripts/typing/trajectory_singlecell_features.py
Rscript        core/phenocycler/scripts/typing/trajectory_singlecell.R 50000    # ~23 min
.venv/bin/python core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py
```

`trajectory_islet.py` gates on **method** checks and refuses to overwrite the deployed file when
one fails, writing `adata_ins_root.rejected.h5ad` instead. It does **not** gate on a biological
result: the legacy lane required pseudotime to anti-correlate with INS and to put T1D after ND, and
on the first run here those would all have passed while the ordering reversed under a change of
neighbourhood size. What gates is finiteness, a real ND root above the cell floor, stability under
neighbourhood size and re-rooting, and that the axis is not mostly donor. The disease ordering is
computed and reported, never required.

`trajectory_singlecell.R` runs the published algorithm — destiny diffusion map into slingshot,
default parameters, no start/leaf nodes, no cluster labels — over a proportional per-donor
subsample. Cells outside the subsample carry **no** pseudotime; none is imputed.

## Things this lane has learned the hard way

- **`pgrep -f` needs a bracket** (`positivit[y].py`), or the pattern matches the shell running it.
  A `pkill` written without one killed three unrelated background jobs in one go.
- **Long stages must be detached** (`nohup … &`). Streaming an 83 GB qptiff leaves little memory
  *free* while leaving ~480 GB *available*, and supervisors watching the wrong number kill the run.
- **anndata 0.13 writes pandas-3 strings as `nullable-string-array`**, which the app's R `anndata`
  reader does not understand and falls back from **without an error** — feature names arrive as
  `1, 2, 3 …` and string columns arrive empty. Force object dtype before writing, and check the
  HDF5 encodings rather than that the file opens.
- **The per-cell polygons exist in exactly one place**, `data/redsea_scratch/geojson/` — no
  PhenoCycler stage output (current flat roots or historical run tree) carries them. They are in
  pixels; convert with each image's own `pixel_size_um_x`.
- **R integer arithmetic overflows at 2^31.** `50000L * 122951L` silently becomes NA.
- **A proportion over an empty region is not zero.** Three islets have no peri cell at all; writing
  `0.0` there states "0 % of this islet's peri cells are Alpha" about an islet with no peri cells.

## Retired, and not to be revived as fallbacks

The Nimbus/Leiden chain this README used to describe (`substrate.py` → `positivity.py` → `call.py`
→ `validate.py` → `export_qupath.py`), the RESTORE pair webs, the hierarchical / simultaneous /
probabilistic typers, the marker-debug and Otsu lanes, the Gate 2–5 apparatus, and the broad-typing
producer withdrawn 2026-09-08. Git has them. A README describing a method that did not produce the
data on disk is how two runs were lost to the exocrine question.

**Never claim acinar or ductal.** The 60-plex carries no amylase, trypsin, CPA1 or elastase, and no
CK19, CA19-9, SOX9, HNF1B, CFTR or MUC1. Exocrine is one label.
