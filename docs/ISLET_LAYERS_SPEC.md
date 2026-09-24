# Islet Explorer layer artifacts — the contract

One spec, written so a worker can build against it without asking anyone. Every number in it was
measured on 2026-09-16, and the command that measured it is given so you can re-measure rather
than trust this file.

The app reads **layer artifacts**: small, pre-aggregated parquet written by a builder from the big
source trees. Cards never scan the source trees at render time.

- Builders write to `$ISLET_LAYERS_DIR`, default `/data/islet_layers`.
- **All nine roots now exist** (measured 2026-09-18: 430 MB — `availability`, `compartments`,
  `distance_profiles`, `grid_composition`, `islet_summary`, `pln`, `restain`, `sections`,
  `trajectory`). The line that stood here said nothing existed yet; that was true when this file
  was written on 2026-09-16 and stopped being true the same week. Sections 1.8 and 1.9 were added
  on 2026-09-18 for the two roots the original seven did not cover.
- A card that cannot find its artifact renders a **named note**, never a blank panel and never an
  empty plot. `ie_require_measured()` in `apps/R/ieCore/R/ie_availability.R` returns that note.

---

## 0. Rules that apply to every builder

### 0.1 Reading the source trees: `partitioning=None`, always

`pyarrow.dataset(root, partitioning="hive")` **fails** on four of the source roots. This is not a
corrupt file. The hive partition key inferred from the directory name collides with the file's own
`donor_id` column, which is `large_string`:

```
partitioning="hive"              FAILS  ArrowTypeError: Field donor_id has incompatible types
partitioning=None                OK     1,785,997 rows, donors ['0112', '6533']
explicit string partition key    FAILS  the file column is large_string, not string
```

So:

```python
import pyarrow.dataset as ds
d = ds.dataset(root, format="parquet")        # partitioning=None is the default
```

Affected: `core/phenocycler/data/expression`, `core/phenocycler/data/integration/cells_xen`,
`core/phenocycler/data/integration/cells_pheno`, `data/development/typing/immune`.

### 0.2 …except `typing/cells`, which has the opposite problem

`data/development/typing/cells` has **24 columns and no `donor_id` among them**. Its donor exists
only as the partition directory. Reading it with `partitioning=None` silently gives you 25 million
rows with no way to tell the donors apart.

| tree | has `donor_id` column? | how to get the donor |
|---|---|---|
| `typing/cells` | **no** | parse the `donor_id=` directory |
| `expression`, `cells_xen`, `cells_pheno`, `typing/immune` | yes (`large_string`) | the column; `partitioning=None` |

Read file-by-file and attach the donor yourself:

```python
import glob, pyarrow.parquet as pq
for path in sorted(glob.glob(f"{root}/donor_id=*/*.parquet")):
    donor = path.split("donor_id=")[1].split("/")[0]
    table = pq.ParquetFile(path).read()
```

### 0.3 `donor_id` is a string, and it is not normalised

`0112` and `112` are the same donor. The pLNpanc runs pad to four characters
(`donor_id=0112`); the cohort does not (`donor_id=115`, `donor_id=6374`). A glob for
`donor_id=112` finds nothing.

**Every builder normalises on read** with one shared helper, and writes the *unpadded* form:

```python
def norm_donor(value: str) -> str:
    """'0112' -> '112'. Never int(): 6533 must not become 6533.0 through a float column."""
    return str(value).strip().lstrip("0") or "0"
```

Two builders normalising differently is the exact failure the merge step looks for. Do not
write your own.

### 0.4 Excluded donors are filtered at the build side

**Only `6457`** is centrally excluded in `core/phenocycler/phenocycler/cohort.py`. **6579's
exclusion was LIFTED 2026-09-17** and it is recorded in `REINSTATED` there — the pancreatitis it
was excluded for is still a disease context, so stratify or report it separately, but do not
drop it. (Corrected 2026-09-18: this line named 6579 as excluded, which line 384 of this same
file already contradicted.) Derive from `DONOR_EXCLUSIONS`, never a literal, and filter in the
builder, not in the card. A card must never be the thing that remembers an exclusion.

### 0.5 Every coordinate column carries a frame

Frames are identities, not scales — `apps/R/ieCore/R/ie_frames.R`. A builder that writes
coordinates writes the frame beside them, as a column:

```
frame  = "xenium_um:6436:panc"      # kind:donor:roi
```

Valid kinds: `xenium_um`, `if_um`, `phenocycler_um`, `hex_qr`, `signed_distance_um`, `none`.
`hex_qr`, `signed_distance_um` and `none` are scale-free and take no donor or roi.

Two frames are the same frame only when all three parts match. **Two tissues on one slide are two
frames.** Sections mounted together came from different blocks.

### 0.6 Absence has five states, and none of them is zero

`measured | blank | unavailable | withheld | failed` — `apps/R/ieCore/R/ie_availability.R`.
A builder that cannot produce a value writes the state and the measured reason, e.g.
`unavailable(no_islet_outline_exported)`. It never writes `0`, `NA` or an empty frame and lets the
card guess.

---

## 1. The artifacts

Each builder writes exactly one root. Paths are relative to `$ISLET_LAYERS_DIR`.

### 1.1 `availability/` — worker F

`build_availability.py`. Donor × capability, one row per pair.

| column | type | meaning |
|---|---|---|
| `donor_id` | string | normalised, unpadded |
| `capability` | string | `phenocycler_panc`, `restain_if`, `islet_outline`, `trajectory` |
| `state` | string | one of the five |
| `reason` | string | measured, e.g. `no_islet_outline_exported`; empty only when `state == "measured"` |
| `n` | int64 | the count the state is based on; **null, never 0, when not measured** |
| `source` | string | the path the state was read from |

This is the table every other card consults before claiming a donor has something.

Each capability is read by a card: `phenocycler_panc` and `islet_outline`
(`card_x1_gradient.R:284`), `restain_if` (`card_x2_overlay.R:47`), `trajectory`
(`card_s2_trajectory.R:56`). Until 2026-09-23 it also carried `phenocycler_pln`, `xenium_panc`,
`xenium_pln` and `pln_nodes`; no card read any of them, so they are no longer written.

### 1.2 `sections/` — worker F

`build_sections_layer.py`. One row per physical section.

Source of truth is `core/phenocycler/data/integration/restain_sections.csv` (5 rows, hand-recorded)
and `core/phenocycler/data/integration/xenium_paths.csv` (56 rows, 26 donors).

Columns: `donor_id`, `roi`, `modality`, `image_file`, `xenium_serial`, `pixel_size_um`, `frame`.

**`pixel_size_um` comes from the CSV and from nowhere else.** Do not read it from the qptiff. The
maintainer declares all post-Xenium IF to be `0.2505`; two of those files record `0.50139` and the
CSV header documents that disagreement. The declared value governs.

**`roi` DOES select a Xenium bundle. An earlier issue of this spec said the opposite and was
wrong.** The `panc` and `pln_1` rows share a *serial* — one run — but name **different output
bundles**, distinguished by a region token:

```
115  panc   MM035_Xenium_nPODBatch1/output-XETG00298__0041117__PanT__20240919__203848
115  pln_1  MM035_Xenium_nPODBatch1/output-XETG00298__0041117__pLN__20240919__203848
```

They hold different cells. Reading the `panc` bundle for `roi=pln_1` gives the wrong tissue with no
error, and treating the split as geometric-within-one-bundle means the pLN cells are never loaded at
all. Resolve the bundle from the row's own `source_path`, never from the serial.

`6591,pln_1` is `NA`: recorded unavailable, and it stays that way.

### 1.3 `distance_profiles/` — worker G

`build_distance_profiles.py`. **Two arms that are not the same kind of measurement**, and the
artifact keeps them apart by `roi`, `tissue` and `modality`.

**The pancreas arm** (`roi=panc`, `modality=phenocycler`) is spatial: immune fraction against signed
distance to the islet edge, `frame=signed_distance_um`.

Reads `data/development/typing/cells` (20 files, **25,230,109** rows, 24 cols) and
`data/development/typing/islet_membership` (**20,509,916** rows, 7 cols: `object_id`, `in_islet`,
`in_expansion_outline`, `islet_num`, `islet_num_parent`, `signed_dist_islet_um`, `donor_id`).

Join on `object_id`. That UUID is the only join key; never centroid matching, never donor-name
parsing.

**The pLN arm** (`roi=pln_*`) is **compositional only**. `frame=none`, `distance_bin_um` is null and
no `compartment` is assigned. The pancreas and pLN sections mounted on one slide are cut from
**different blocks**: they share a pixel grid and nothing anatomically, and there is no islet in a
lymph node to measure a distance from. It joins to the pancreas arm on `donor_id` and `donor_status`
and on nothing else. `refuse_cross_tissue_distance` and `_assert_no_cross_tissue_distance` make that
structural: a pLN row carrying a bin, a compartment or a signed-distance frame raises instead of
being written. A pLN row carrying only COUNTS is the arm working, and is not refused.

Both pLN modalities are read, neither is waited on, and **they are never pooled into one fraction**:

| modality | source | when absent |
|---|---|---|
| `xenium` | the pLN release under `data/development/xenium_typing/` (found by its declared `tissue`, not its directory name) | falls back to the `cells_xen` ingest, which carries no call on any pLN row |
| `phenocycler` | `data/development/typing/pln/pheno/` at `donor_id=*/roi=*/*.parquet` | `unavailable(phenocycler_pln_calls_absent…)`, one row per donor the pLN **ingest** shows a section for |

Hive partitions are read with an explicit string schema. Every pLN donor id is all digits, so
inferred partitioning turns `0112` into the integer `112` and the id stops joining anywhere else.

Columns: `donor_id`, `donor_status`, `roi`, `tissue`, `modality`, `distance_bin_um`, `n_cells`,
`n_immune`, `immune_fraction`, `immune_subtype`, `resolved_by`, `frame`, `compartment`, `state`,
`reason`.

Five things this builder must do:

1. **Refuse `typing/subpanel/subpanel_calls.parquet` by path**, with the superseded counts in the
   error. That file is a previous generation — Beta 191,485 / Alpha 246,055 / Delta 62,797 against
   the current 226,183 / 279,686 / 90,658.
2. Carry `resolved_by` through, so the card can state provenance rather than imply it.
3. Use `islet_membership`, not QuPath's `Parent`. `Parent` mislabels 2.7% of islet-core cells,
   which fall through to Exocrine.
4. **Read which lineages are immune; never infer it from a label's spelling.** `is_immune` and
   `immune_tissue_scope` come from `core/phenocycler/data/integration/vocab_crosswalk.csv`, which is
   written by `phenocycler.integration.vocab` and by nothing else. Maintainer ruling 2026-09-17:
   *all immune cell calls from pancreas can be applied to pLN* — that licenses the **identity**, not
   the geometry, and it is not a licence to copy label rows, because the panc and pLN bundles are
   segmented separately and their cell id namespaces are disjoint. A lineage the crosswalk does not
   classify is withheld **by name**, never defaulted to non-immune.
5. **Count only what was answered.** A cell with no nucleus is excluded from every calculation. The
   pLN denominator is the cells stage 1 resolved (`Immune` or `Non_immune`); `Ambiguous` goes in
   neither the numerator nor the denominator and is emitted on its own `blank` row
   (`resolved_by=stage1_unresolved`) so it stays visible. The numerator is stage 1's `Immune`, so the
   stage-2 lineage rows **do not sum to it** — an immune cell stage 2 could not resolve is still
   immune.

Every source path the artifact records is validated by the app-side reader, not just
`source_domain_root`: `source_cells` and `source_membership` must lie in the active domain's tree,
and `source_xen_cells` / `source_immune_vocabulary` must be recorded once something depends on them.
`value$n` describes the rows a card draws, not the artifact's row count.

### 1.4 `restain/` — worker H

`build_restain_layer.py`. The post-Xenium IF re-stain overlay, **in the IF frame**.

Assert `n_xenium_warped_to_if == n_xenium_cells` before writing anything. `warped_if_um.npy` is the
Xenium centroids warped *into* the IF frame; the reverse direction is the mistake an independent
reader is most likely to catch, so state the direction in the docstring.

Columns: `donor_id`, `cell_id`, `object_id`, `x_um`, `y_um`, `frame` (`if_um:<donor>:panc`),
`feat__if_INS`, `feat__if_GCG`, `hormone_positivity_*`.

Only 6436 and 6442 carry `feat__if_*` today. Every other donor is
`unavailable(no_restain_pairing)`.

### 1.5 `islet_summary/` — worker I

`build_islet_summary.py`. Islet size against composition.

Reads, all of which exist (**re-measured 2026-09-18**; the figures that stood here — 31,128 /
111,823 / 7,311 — were the 2026-09-14 ten-donor build and were superseded by the 20-donor rebuild
of 2026-09-17, `docs/APP_COVERAGE_2026-09-17.md`):

| file | rows × cols |
|---|---|
| `data/development/typing/app/islet_geometry.parquet` | 73,616 × 13, 20 donors, 18,404 core islets |
| `data/development/typing/app/islet_composition_long.parquet` | 255,753 × 8, 20 donors |
| `data/development/typing/geometry` | 20 files, 17,553 rows × 5 cols in total, carries the `wkb` outline |

`islet_geometry` columns: `donor_id`, `islet_key`, `combined_islet_id`, `region`, `region_um2`,
`n_cells`, `cells_um2`, `islet_diam_um`, `centroid_x_um`, `centroid_y_um`, `area_um2_*`.
`islet_composition_long`: `combined_islet_id`, `islet_key`, `phenotype`, `n_cells`, `region`,
`region_total`, `prop`, `donor_id`.

**Two stratifications are mandatory, not stylistic**, and the builder emits the columns that make
them possible:

- `n_cells` per islet, because the endocrine shape signal is real only above ~10 cells per islet;
- a `donor_id` column kept unpooled, because 6539's islets are mostly one-cell classifier
  fragments and a cohort-pooled per-islet number is wrong without saying so.

The card shows the size distribution *before* any composition claim.

Do not use `typing/app`'s islet-level hormone composition as positivity: it is degenerate
(`Ins_any == Ins_single`, `Multi_Pos = Triple_Neg = 0` across all 18,063 rows — an argmax returns
one winner; re-measured 2026-09-18, the row count was 7,634 on the ten-donor build). The genuine per-donor thresholded positivity is
`data/app_data/parquet/uns_composition.parquet`. They are not interchangeable.

### 1.6 `trajectory/` — worker J

`build_trajectory_layer.py`. Within-donor islet ordering.

Reads `data/development/typing/trajectory/singlecell_islet_summary.parquet` (**18,208 × 9, 20
donors**: `case_id`, `islet_key`, `sc_pseudotime_core`, `sc_pseudotime_mean_core`, `sc_n_core`,
`sc_pseudotime_peri`, `sc_pseudotime_mean_peri`, `sc_n_peri`, …) and
`singlecell_pseudotime.parquet` (**480,444 × 33, 20 donors**: `donor_id`, `object_id`, `islet_key`,
`region`, `phenotype`, `donor_status`, `stain_group`, `X_centroid`, `pseudotime`, `DC1..DC20`, …).

*(Corrected 2026-09-18. The 7,179 × 9 and 49,999 × 32 figures were the 2026-09-14 pooled ten-donor
build. The ordering was refit **within donor** on 2026-09-17 — one destiny → slingshot fit per
donor, no pooled embedding — see `docs/TRAJECTORY_WITHIN_DONOR_2026-09-17.md`.)*

**Order islets within a donor. Never pool.** *(Corrected 2026-09-18: the η² = 0.988 that stood here
describes **neither** artifact on disk. It comes from the 2026-07 audit in the header of
`scripts/senior/compute_trajectory_senior.py`, of a legacy CODEX-domain trajectory that no longer
exists; the current legacy-domain file measures 0.091 on donor, and the deployed single-cell
ordering measures 0.0093. Measured in `docs/TRAJECTORY_WITHIN_DONOR_2026-09-17.md` §1. The rule
stands on the batch/status alias itself — Cramér's V against `donor_status` is 0.604 on coarse
batch, 0.716 on `stain_group`, 0.882 on `acq_run` — not on that figure.)* The builder emits
`batch_eta_sq_within_status`; under a within-donor fit that quantity is `unavailable` with a
reason, because every donor is rank-uniform on its own axis, and the card prints the reason
instead.

`sc_pseudotime_peri_minus_core` is the quantity that speaks to the study question directly: how far
the peri-islet compartment has moved relative to the core, in the same islet.

**State the sampling on the card.** The ordering is fit on **480,444 of 1,639,174 cells (29.3%)**
— corrected 2026-09-18; 49,999 of 800,562 was the superseded pooled build. That belongs in the
caption, not a footnote.

### 1.7 `pln/` — increment 9

`build_pln_layer.py`. The lymph node arm of the compartment comparison.

The structural unit is recorded in `core/phenocycler/data/integration/pln_node_annotations.csv`
(tracked, hand-recorded): 19 images, 18 eligible, 8,961,168 eligible cells, 44 nodes. Donor 6380's
single `Tissue` annotation **is** the node — it is the only entry whose whole annotation set is
`{"Tissue": 1}`, and every other entry carries `Tissue` alongside its `pLN_*` outlines, so reading
`Tissue` as the node anywhere else double-counts.

Node membership is derived **geometrically**, polygon against `Centroid X/Y µm`. Never from the
TSV's `Parent`: the export carries `Parent = "Annotation (Tissue)"` on 100% of rows.

Donor 6579 carries 1,231,818 cells and is **not** excluded: `pln_node_annotations.csv`
has recorded `excluded=False` for it since 2026-09-17 14:16, following `cohort.py`'s REINSTATED
entry of the same day. (Corrected 2026-09-18 — this line instructed an exclusion that had lapsed.)

### 1.8 `compartments/` — added 2026-09-18

`scripts/senior/build_compartment_gradient.py`. Immune composition by compartment, and the
pancreas/pLN pairing. Four files plus a manifest; CSV, not parquet.

`compartment_composition.csv` (465 × 16) is the table the card reads: `donor_id`, `tissue`,
`modality`, `compartment`, `stratum`, `lineage`, `n_cells`, `n_lineage`, `n_immune`,
`fraction_of_compartment`, `fraction_of_immune`, `state`, `reason`, `donor_status`, `stain_group`,
`acq_run`.

**The compartment is derived from the signed distance, never from the drawn annotation**, and the
manifest records why: one donor's drawn expansion measures 23.1 µm against 18.8–19.3 µm for every
other, and all of them are named `_exp20um`.

| compartment | definition |
|---|---|
| `core` | `signed_dist_islet_um <= 0` |
| `peri` | `0 < signed_dist_islet_um <= 20.0` |
| `exocrine` | `signed_dist_islet_um > 20.0` |
| `pln` | compositional only — no distance exists in a lymph node |

`stratum` carries the islet-size stratification: `min_islet_cells = 10`, because the endocrine
shape difference is real only above about ten cells per islet and 6539's islets are mostly
one-cell classifier fragments. `islet_size_summary.csv` (20 × 6) is the per-donor evidence for
that — `n_islets`, `n_islets_ge_10`, `median_core_cells_per_islet`, `fraction_islets_one_cell`,
`flag` — so a card can say which donors the stratum excludes rather than implying the cohort is
uniform.

`batch_eta_squared.csv` (105 × 14) reports η² on `stain_group` and `acq_run`, each against a
permutation null. **Nothing is gated on these numbers.** Coarse batch is aliased with disease
status, so it is not corrected on; the finer variables are what carry a correction.

The manifest's `pancreas_diagnostics` carries, per donor, `n_rows`, `n_eligible`, what was dropped
(`exocrine_inside_islet_border`), `cells_without_signed_distance`, `immune_calls_not_in_crosswalk`,
and a `membership_audit` comparing core membership by distance against core membership by outline
(`n_core_by_distance`, `n_core_by_outline`, `n_disagree`). That disagreement is a measurement to
report, not a defect to silence.

`pending_maintainer_classification` names the one thing that would license an islet-to-lymphoid
distance: an `attached_to_pancreas_section` column, arriving **as data**. Until it lands, every
pancreas/pLN pair joins on `donor_id` and `donor_status` only.

### 1.9 `grid_composition/` — added 2026-09-18

`scripts/senior/build_grid_composition.py`. The tissue × modality grid: pancreas/pLN against
protein/RNA, one common lineage vocabulary through
`core/phenocycler/data/integration/vocab_crosswalk.csv` (generated by `run_vocab` in
`phenocycler/integration/vocab.py`; read, never hand-edited).

| file | shape | what it is |
|---|---|---|
| `composition.parquet` / `.csv` | 825 × 19 | donor × tissue × modality × lineage counts and fractions |
| `presence.parquet` / `.csv` | 100 × 13 | donor × tissue × modality: state, reason, cells assessed, cells excluded for no nucleus |
| `fine.parquet` | 705 × 9 | the finer lineage each lane actually resolved, beside its common name |
| `crossmodal_tnk.csv` | 45 × 11 | the T/NK fraction in each modality, with the rank disagreement stated |
| `batch_eta_squared.csv` | 75 × 12 | η² per lineage on the batch variables, with a permutation null |

Measured in the manifest written 2026-09-17 21:19 UTC — the four grid cells and the sources they
were built from:

| tissue | modality | donors measured | cells assessed | no nucleus, excluded | source |
|---|---|---|---|---|---|
| pancreas | protein | 20 | 25,230,109 | 4,748,402 | `data/development/typing/cells` |
| pancreas | RNA | 18 | 13,725,222 | 775,234 | `/data/xenium_typing/panc_2026-09-17/release` |
| pln | protein | 14 | 5,122,269 | 0 | `data/development/typing/pln/pheno` |
| pln | RNA | 20 | 5,927,317 | 80,472 | `data/development/xenium_typing/release_pln_immune_v3` |

The pancreas RNA cell also records `n_nucleus_not_measured = 2,989,849` — cells whose nucleus
status is **unavailable**, kept distinct from the 775,234 measured as having none. They are not
added together and neither is zero.

Three rules the manifest states and the builder enforces:

- `nucleus_rule` — a cell with no nucleus is excluded from every calculation.
- `missing_rule` — a lineage a lane did not resolve carries `n_cells` NULL and a reason, never 0.
- `cross_tissue_rule` — compositional only, joined on `donor_id` and `disease_status`; no distance,
  adjacency, neighbourhood or KDE crosses the pancreas/pLN boundary.

`coherence` is the builder's own arithmetic check: every bucket of a grid cell sums to that cell's
assessed cells. 72 grid cells checked, 0 mismatched. It checks the build, not a result.

`n_donors_in_universe` is 25, wider than any one cell, and the manifest carries both
`cohort_exclusions` (6457) and `cohort_reinstated` (6579, lifted 2026-09-17 — the disease context
is not). A card reads the per-cell `donors_measured`, never the universe.

---

## 2. What a card may and may not do

- **Never co-plot two frames.** `ie_assert_same_frame()` inside every plot reactive; it `stop()`s on
  a mismatch, on a missing attribute, and on a frameless input.
- **No transform is applied in the app.** Reprojecting 6442's `x_um` through the stored affine and
  differencing against `x_um_reg` gives max 207 µm / mean 56.8 µm.
- **`PIXEL_SIZE_UM` is forbidden in new files**, by test. `segmentation_helpers.R` is an explicit
  pre-existing exemption and is not to be reused as precedent.
- **No expected-fraction gates.** Report raw numbers. Categorise each cell by what it expresses.
- Compartment comparisons across tissues join on `donor_id` and status only. No distance,
  adjacency, neighbourhood or KDE crosses a tissue boundary.
- **A layer artifact records the domain root it was built from, and a card checks it.** The two
  domains share five donors and an islet key exists on both sides, so a correct-looking artifact
  from the other domain is exactly the failure that looks like a result. Measured 2026-09-18:
  eight of the nine `card_*.R` read a layer directory and **five** of them carry the check —
  `card_s1_size_composition`, `card_g1_grid`, `card_s3_compartments`, `card_n1_sections`,
  `card_x1_gradient` — each refusing under its own named reason rather than a shared one. The
  other three do not check. That is a gap in the cards, not in this spec, and it is recorded here
  so a new card is written with the check rather than without it. `card_x1_gradient` goes
  further and checks every recorded source path, not only the one naming the domain: a builder
  can read the wrong cells and still stamp the right domain root.

---

## 3. Re-measuring this file

```bash
cd /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior
.venv/bin/python - <<'PY'
import glob, pyarrow.parquet as pq
for root in ("data/development/typing/cells",
             "data/development/typing/islet_membership",
             "core/phenocycler/data/integration/cells_xen"):
    fs = sorted(glob.glob(root + "/**/*.parquet", recursive=True))
    n = sum(pq.ParquetFile(f).metadata.num_rows for f in fs)
    s = pq.ParquetFile(fs[0]).schema_arrow
    print(f"{root}: {len(fs)} files, {n:,} rows, {len(s.names)} cols")
PY
```

**One correction to the approved plan.** It states `cells_xen` holds 15,981,833 cells. Measured
across its 20 files on 2026-09-16 the total is **14,980,833** — exactly 1,001,000 fewer. The plan's
figure has not been reproduced; use the measured one, or re-measure before citing either.
