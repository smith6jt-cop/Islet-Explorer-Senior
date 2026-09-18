# The app side of the within-donor trajectory — 2026-09-17

The build track refit the single-cell trajectory **within donor** for all twenty donors. This is
what it took to put that on the screen, what moved, and the four things it could not reach because
the file belongs to a concurrent lane. Every number below was produced by running something in this
checkout; the commands are named where they are not obvious.

---

## 1. Which script owns `adata_ins_root.h5ad`: neither of the two that look like it

The task asked which of `scripts/reaggregate_islets.py` and `scripts/recompute_viz_umap.py` owns
the file. **Neither does, and neither can run.** Both write `<repo>/data/adata_ins_root.h5ad` —
not the run-cohort domain root — and every input they name is absent from this checkout:

```
data/adata_ins_root.h5ad                                        MISSING
single_cell_analysis/CODEX_scvi_BioCov_phenotyped_newDuctal.h5ad MISSING
islet_analysis/                                                 MISSING
data/app_data/adata_ins_root.h5ad                               EXISTS   <- the CODEX domain's copy
```

They are the legacy CODEX lane's entrypoints. They were **left in place rather than removed**,
against the repo rule about not leaving two runnable entrypoints, for one reason:
`docs/TRAJECTORY_PARAMETER_AUDIT_2026-09-15.md` cites `reaggregate_islets.py` **by line number**
sixteen times as the provenance of parameters that are still live in `trajectory_islet.py`
(`n_neighbors=30` at :148, `K_ROOT=50` at :170, `n_comps=15` at :448, the 2026-04-25 alpha sweep at
:223–234). Deleting the file, or editing it and shifting its lines, silently voids those citations.
They are dead code that fails loudly on a missing input rather than quietly on a wrong one, which
is the tolerable failure of the two. **The maintainer's call, not this session's.**

The run cohort's producer is `core/phenocycler/scripts/typing/trajectory_islet.py`, which belongs
to the islet-level UMAP lane — and on 2026-09-17 it **refused to write the file**. Its pooled
twenty-donor islet-composition DPT failed its own `pseudotime_is_not_mostly_donor` check (donor
η² 0.444 within T1D, 0.314 within Aab+, against a 0.15 ceiling), so it wrote
`adata_ins_root.rejected.h5ad` and left the deployed ten-donor file alone. That is the behaviour it
was built to have, and it is why the deployed file was still 10 donors, 2 statuses and a superseded
typing generation at 20:50 today.

### The new producer

`scripts/rebuild_adata_ins_root.py`. It takes the twenty-donor **frame** out of the rejected file —
composition, ALR, the islet-composition UMAP, the batch variables — and replaces the one quantity
that was rejected, the pooled ordering, with the within-donor one summarised per islet.

```bash
.venv/bin/python scripts/rebuild_adata_ins_root.py
```

It refuses to run if `singlecell_trajectory.json` does not record `fit_scope: within_donor`, so it
cannot be pointed at a pooled fit by accident. It copies the previous output aside before writing
(`adata_ins_root.superseded-<mtime>.h5ad`) and deletes nothing.

---

## 2. Before and after

| | before (2026-09-14 17:12) | after |
|---|---|---|
| islets | 7,302 | **17,786** |
| `obs["imageid"]` | 10: 115, 6374, 6436, 6442, 6461, 6476, 6482, 6493, 6539, 6593 | **20**: the ten above **plus 6380, 6450, 6505, 6521, 6523, 6534, 6551, 6566, 6591, 6623** |
| `obs["donor_status"]` | `['ND', 'T1D']` | **`['ND', 'Aab+', 'T1D']`, ordered** |
| islets with a core ordering | 4,371 | 6,625 |
| islets with a paired peri-minus-core difference | — | 6,528 |
| ordering | pooled islet-composition DPT, 10 donors | **within donor, one destiny→slingshot fit per donor, 20 donors** |
| cells behind it | 49,999 of 800,562 (6.2%) | **480,444 of 1,639,174 (29.3%)** |
| typing generation | superseded (Sep 11 lane) | `typing/cells/` |

`donor_status` is written as an **ordered** categorical with categories in cohort order — ND, Aab+,
T1D — so every downstream `levels()` and every legend inherits it rather than restating it. Donor
colours come from the domain's declared palettes through `00_globals.R`; no palette was added.
Verified in the rendered tab: the legend draws ND `#3366CC`, Aab+ `#FF6600`, T1D `#109618`, in that
order, which are the domain's "Bright" swatches.

### Per donor

| donor | status | stain group | islets | with a core ordering | with a paired difference |
|---|---|---|---:|---:|---:|
| 115 | ND | S2 | 1,298 | 441 | 436 |
| 6374 | ND | S3 | 1,001 | 519 | 516 |
| 6380 | T1D | S1 | 66 | 56 | 56 |
| 6436 | T1D | S5 | 576 | 291 | 290 |
| 6442 | ND | S4 | 486 | 274 | 273 |
| 6450 | Aab+ | S7 | 438 | 213 | 213 |
| 6461 | ND | S4 | 926 | 424 | 422 |
| 6476 | ND | S0 | 1,253 | 518 | 513 |
| 6482 | ND | S3 | 487 | 393 | 393 |
| 6493 | ND | S4 | 686 | 349 | 347 |
| 6505 | Aab+ | S6 | 2,958 | 521 | 474 |
| 6521 | Aab+ | S0 | 377 | 234 | 232 |
| 6523 | T1D | S8 | 359 | 223 | 222 |
| 6534 | T1D | S7 | 1,981 | 422 | 406 |
| 6539 | ND | S4 | 442 | 168 | 167 |
| 6551 | T1D | S1 | 413 | 246 | 244 |
| 6566 | T1D | S7 | 865 | 314 | 312 |
| 6591 | T1D | S7 | 1,131 | 414 | 411 |
| 6593 | T1D | S4 | 417 | 253 | 253 |
| 6623 | T1D | S6 | 1,626 | 352 | 348 |

**6539 is the stratification case, visible here.** 442 islets, 168 with a core ordering — it has
the lowest ordered fraction of any donor, which is what one-cell classifier fragments look like
once a five-ordered-cell floor is applied. Nothing in this file pools per-islet numbers across
donors, so it does not have to be corrected for; it has to be *seen*, and it is.

**The Aab+ caveat is a cohort fact and is written into the file.** Two of the three Aab+ donors
(6450 S7, 6505 S6) are in batch 2, which contains no ND donor; 6521 is the only Aab+ donor inside
batch 1. `uns['trajectory_validation']['cohort_limit']` says so. The cohort supports **showing**
Aab+ islets. It does not support an ND < Aab+ < T1D claim, and this figure is not one.

---

## 3. What the file now carries, and what is deliberately not in it

`dpt_pseudotime` and `dpt_pseudotime_combined` keep their names because the app reads them by name
and its server module belongs to another lane. **What is in them changed completely**, and the file
says so in `uns['trajectory_construction']`: within-donor quantiles, one fit per donor, no pooled
embedding, and an explicit note that two donors' values are not on one axis. The honest names are
carried beside them in `obs` — `sc_pseudotime_core`, `sc_pseudotime_peri`,
`sc_pseudotime_peri_minus_core`, with `sc_n_core` / `sc_n_peri` beside each.

### `obsm['X_umap']` is pooled and the ordering is not — both, in one file

The islet-composition UMAP is carried over unchanged from the islet-level lane, which measured and
decided it separately (no cluster at any resolution more than 80% one donor). It is a description of
composition. It is **not** the ordering, it shares no construction with it, and reading the two as
one picture is the mistake the file's `umap_is_not_the_ordering` key exists to head off. The tab
now states it too.

### "Core + Peri" needed a second condition

The combined value is a median over the islet's ordered core **and** peri cells. A cell-count floor
alone is not enough: measured on this build, the five-cell rule by itself admits 15,265 islets, and
**1,840 of those have no ordered core cell at all**. A median of twelve peri cells labelled
"Core + Peri" tells the reader it contains a compartment it does not. The rule is therefore five
ordered cells *and* at least one in each region: 13,423 islets, 13,325 of them inside this frame.

### Unorderable islets carry a named reason

Not a blank, not a zero:

| reason | islets |
|---|---:|
| *(ordered)* | 6,625 |
| `below_5_sampled_core_cell_floor` | 7,602 |
| `no_core_cell_in_this_donor_ordered_subsample` | 3,414 |
| `no_ordered_cell_in_either_region_for_this_islet` | 145 |

And one more, recorded because it would otherwise vanish in a join:
`islets_ordered_but_absent_from_this_frame = 567`. Those islets have a within-donor ordering and no
row in the islet-composition frame, so the tab cannot draw them. They are measured; they are not
drawable here.

### No batch η² is printed on this file, and that is the deliberate part

The app's caption prints `uns['trajectory_validation']['within_status_donor_eta2']` under the words
"within-status donor η²", which on a **pooled** ordering is the batch check with a 0.15 ceiling. On
a within-donor fit it is not that check. Every donor is rank-uniform on its own axis, so what the
pooled number measures is where each donor's core cells sit inside that donor's own spread — a fact
about the summarising. **Printed beside a green "gates pass" badge it would read as a ceiling that
was cleared.** So the key is written `None` with a reason, and the number is carried in full under
its own name:

```
within_status_donor_eta2_of_islet_medians = {ND: 0.0801, Aab+: 0.5164, T1D: 0.2715}
within_status_donor_eta2                  = None
within_status_donor_eta2_reason           = not_applicable_to_a_within_donor_fit(...)
```

This is the same move `build_trajectory_layer.py` made for the same reason at the layer, and it is
why the card's caption now prints no η² either. The cross-donor question is asked separately and
properly — on the *relation* between the twenty orderings, against a permutation null over donors —
and that result is carried in `uns['crossdonor_relation']`, where it gates nothing and sorts
nothing.

### The gates

Four, all structural, all measured on the object being written, all passing:
`gate_every_donor_has_its_own_fit`, `gate_no_pooled_embedding`,
`gate_ordered_islets_all_finite`, `gate_every_missing_value_carries_a_reason`. `hard_pass = True`
is the conjunction of exactly those four and claims nothing else. No gate was invented to be
failed and none to be passed.

---

## 4. The card: five standing sentences down to one line

`apps/islet/R/card_s2_trajectory.R` printed five sentences above the figure on every render, one of
which — *"Ranks summarize the existing sampled-cell fit; the fit was not recomputed per donor"* —
became **false** on 2026-09-17. The fit is recomputed per donor. That is the rebuild.

Following the Lymph Node column, which was measured at ~470 words of standing prose and trimmed to
a single line with a modal behind it: standing text is stated once, in
`.s2_trajectory_methods_modal()`; **measured** state stays inline, because measured state is a fact
about what the reader is looking at and standing text is not. Rendered against the real layer:

```
fit within donor on 480,444 of 1,639,174 cells · 441 of 1,303 islets ordered ·
436 with a paired peri-minus-core difference · how this is measured →
```

and, beneath it, the islets that could not be ordered — which used to leave the figure silently,
because the plot filters to measured pairs and nothing said how many rows that dropped or why:

```
unavailable(867 of 1,303 islets have no paired difference): 546 — below the 5-ordered-cell
floor in the core; 134 — no ordered cell in the islet core; 113 — below the 5-ordered-cell
floor in both regions; 69 — no ordered cell in either region; 5 — below the 5-ordered-cell
floor in the peri-islet expansion. They are not zero and not absent data.
```

The downloaded PNG keeps a fuller caption (`.s2_trajectory_plot_caption`), because a figure that
travels alone has no modal to click.

`apps/islet/R/mod_trajectory_ui.R` gets the same one-line treatment on the tab itself, as a
`<details>` rather than a modal — a modal needs a server observer and that module is another lane's.
It carries **no cohort fact**: no donor count, no cell count, no donor id, nothing the loaded data
would have to agree with. Those live in `traj_status` and the figure captions.

---

## 5. The availability ledger, which was hiding all three Aab+ donors

`/data/islet_layers/availability/data.parquet` was built at 12:35 against the ten-donor summary and
recorded `trajectory: unavailable / no_islet_trajectory_exported` for 6380, 6450, 6505, 6521, 6523,
6534, 6551, 6566, 6591 and 6623. The card checks that ledger **before it draws**, so the three Aab+
donors were in the data and not on the screen. Rebuilt:

```bash
.venv/bin/python scripts/senior/build_availability.py
# {"rows": 200, "donors": 25, "states": {"measured": 100, "unavailable": 66, ...}}
```

`trajectory` is now `measured` for **20 of 25** roster donors. The five that are not are named and
correct: 112, 6347, 6533, 6580 have no PhenoCycler pancreas trajectory, and **6579 now appears in
the roster** — reinstated 2026-09-17, no `typing/cells/` partition yet, so `unavailable` with a
reason rather than absent.

---

## 6. Test results

### Tab smoke, per tab

```bash
cd apps/islet && NOT_CRAN=true ISLET_DOMAIN_CONFIG=config/run20_default.yml \
  /opt/R/4.5.1/bin/Rscript tests/tab_smoke.R /tmp/tabshots
```

| tab | R output errors | JS console errors | outputs that exist but drew nothing |
|---|---:|---:|---|
| Plot | 0 | 0 | 2 plotly, both drawn; 0 empty |
| **Trajectory** | **0** | **0** | 3 plots drawn / 3, 2 plotly drawn / 2; 2 empty html — `segmentation_viewer_panel` (no islet selected) and `traj_outlier_info` (table hidden by default) |
| Spatial | 0 | 0 | 2 plots, 2 plotly, all drawn |
| Neighborhoods | 0 | 0 | 3 plotly, all drawn |
| Statistics | 0 | 0 | 2 tables exist and drew nothing — downstream of a Plot-sidebar selection the driver never makes; 10 empty html |
| CrossModal | 0 | 0 | 6 plots drawn / 6 — **FAIL on a registry check, pre-existing, see below** |
| LymphNode | 0 | 0 | 1 plot drawn / 1 |
| Viewer | 0 | **194** | pre-existing: Avivator needs HTTP Range and `shiny::runApp` returns 200 with the whole body |

The **`s2_trajectory` card drew**: `card s2_trajectory ok | plots 1 | plotly 0 | tables 0`.

The **CrossModal FAIL is not from this work.** The check reads
`query(pane, '.ie-card-stack > .ie-card')` and expects seven named cards; it saw eleven, the four
extra having an empty `data-card-id`. They come from `apps/islet/R/card_g1_grid.R:273`, which
renders a **nested** `.ie-card-stack` of `.ie-card` panels with no card id, and a descendant
selector finds them. That file is committed (`3373361`), unmodified in the working tree, and owned
by another lane. Its two plots are counted twice, which is why the four phantom cards carry plots
0/1/0/1.

### R suite

```bash
/opt/R/4.5.1/bin/Rscript -e 'testthat::test_dir("apps/islet/tests/testthat")'
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 772 ]
```

---

## 7. Four things this session could not reach

All four are in `apps/islet/R/mod_trajectory_server.R`, which the running UMAP track owns. **Nothing
in that file was edited.** Each is stated with the measurement behind it so the owning lane can
decide.

### 7a. The coverage line renders empty — the one place the tab states what it drew

`mod_trajectory_server.R:541`:

```r
ms    <- traj_validate()
nm    <- if (is.null(ms)) 0 else nrow(ms)
...
sprintf("Loaded AnnData (%d obs). Joined %d rows (%d unique islets) of %d islet metrics. Source: %s",
        nrow(tr$obs), nj, ndist, nm, basename(src))
```

`traj_validate()` returns `list(feature = ..., status = "ready")` — a **list**, not a data frame —
so `nrow(ms)` is `NULL`, `sprintf()` with a zero-length argument returns `character(0)`, and
`tags$div(character(0))` renders `<div></div>`. Confirmed in the captured DOM:

```html
<div id="traj-traj_status" ...><div></div>
<div style="font-size:12px; margin-top:4px"><span ...>gates pass</span>
  Pseudotime: ρ vs beta proportion +0.16</div></div>
```

The validity badge renders; the sentence that would say *17,786 obs* does not. That sentence is the
tab's only statement of coverage, and the islet-level lane's own decision note requires coverage to
be stated in the UI so a partial fit can never draw silently. It has been silently empty. The fix is
one line — `nm <- if (is.data.frame(ms)) nrow(ms) else 0L`, or drop the term.

### 7b. "Colour by Donor ID" is not legible at twenty donors

`mod_trajectory_server.R:827` and `:1277` both do `brewer.pal(≤12, "Paired")` under thirteen donors
and `rainbow(n, s = 0.8, v = 0.8)` above. Twenty donors takes the second branch. Minimum pairwise
CIE-Lab distance between swatches, measured:

```
Paired(10) = 24.2    rainbow(10) = 27.1    rainbow(12) = 17.2    rainbow(20) = 7.4
```

**7.4 is below the threshold at which two swatches can be told apart side by side**, let alone
matched to a point in a 650px scatter. `rainbow()` is also not a palette from `00_globals.R`; the
palettes there (`DONOR_PALETTES`, `DONOR_COLORS`) are the domain's **disease-status** palettes, and
this repo has no per-donor-id palette at all. So the tab both ships an illegible legend and invents
a palette to do it.

**Proposed, not applied** (it is another lane's file, and the choice is a judgement):
colour by donor status and distinguish donors by facet or by hover, rather than by twenty hues;
or, if per-donor colour must stay, declare a twenty-value qualitative palette in the DOMAIN file
beside the status palettes — where the app already looks — rather than generating one at render
time. Shipping `rainbow(20)` was the one option explicitly ruled out.

### 7c. `pseudotime_unavailable_reason` is in the file and no R file reads it

```
$ grep -rn pseudotime_unavailable_reason apps/islet/R/*.R
(no matches)
```

The column has existed since 2026-09-14. On the tab, an islet that cannot be ordered is filtered out
of the scatter and nothing names why — 11,161 of 17,786 islets, now with four distinct reasons
between them. The **card** names them (§4); the **tab** does not. Making it do so means rendering a
count-by-reason line in `traj_status`, which is that lane's file.

### 7d. The Leiden colour modes are still a ten-donor fit

`traj_color_by` offers four Leiden resolutions, sourced from `app/islet_leiden.parquet`, which is
still the 2026-09-14 ten-donor file. The islet-level lane has the twenty-donor replacement staged at
`typing/umap_rebuild_2026-09-17/islet_leiden_20donor.parquet` and has not deployed it. Until it
does, those four modes colour half the islets and leave the rest grey. Its lane, its call — noted
here only so the two halves of the tab are known to disagree.

---

## 8. Files

| file | what |
|---|---|
| `scripts/rebuild_adata_ins_root.py` | **new.** The producer described in §1. |
| `apps/islet/R/card_s2_trajectory.R` | one measured line + methods modal + named-reason note; the false sentence removed. |
| `apps/islet/R/mod_trajectory_ui.R` | the one-line method strip with the detail in a `<details>`. |
| `data/development/typing/app/adata_ins_root.h5ad` | rebuilt: 20 donors, three statuses, within-donor ordering. |
| `data/development/typing/app/adata_ins_root.superseded-20260914-1712.h5ad` | the deployed ten-donor file, copied aside, not deleted. |
| `/data/islet_layers/availability/data.parquet` | rebuilt with `scripts/senior/build_availability.py`; 20 of 25 donors `measured` for `trajectory`. |

`scripts/reaggregate_islets.py` and `scripts/recompute_viz_umap.py` were **not** changed — see §1
for why, and for the decision that is the maintainer's rather than this session's.
