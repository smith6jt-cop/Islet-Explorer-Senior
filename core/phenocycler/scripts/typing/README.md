# PhenoCycler cell typing, and the app dataset built from it

Two things live here. The **typing** is the subpanel classification — context from QuPath objects,
then one question per context against reference levels drawn from the ordinary tissue. The **app
build** turns those calls into the dataset `apps/islet` serves under
`ISLET_DOMAIN_CONFIG=config/run20_default.yml`.

Development lane: writes only `data/development/typing/`, never `data/runs/`, never a stage
manifest, does not move the `NO_VERIFIED_TYPING` latch.

## The typing itself is not re-derived here

`subpanel_calls.parquet` (25,150,537 rows, 20 donors) already exists and the donors in it are
typed. **Nothing below re-decides `cell_class`.** How it was produced is recorded in
[`RECONSTRUCTION.md`](RECONSTRUCTION.md), recovered by reproducing the calls from their inputs and
verified against them — read that rather than re-deriving it. The two modules carrying the parts
verified exactly are `subpanel_io.py` (input chain, scale, reference levels) and `context.py`
(centroid-in-object assignment, reproduced at 100.0000%).

## The app build, in order

Run from the repo root with the repo `.venv`. Ten clean batch-1 donors by default
(`cohort.py`; batch 2 is refused until its export is corrected — see the root `CLAUDE.md`).

```bash
.venv/bin/python core/phenocycler/scripts/typing/build_cells.py                 # per-cell table
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
- **The per-cell polygons exist in exactly one place**, `data/redsea_scratch/geojson/` — the
  content-addressed run carries none. They are in pixels; convert with each image's own
  `pixel_size_um_x`.
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
