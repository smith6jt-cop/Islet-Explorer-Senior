# The refinement pass, and what measuring it at cohort scale changed

> **STATUS, 2026-09-23: RULED IN, WIRED THROUGH `--refine`, NOT YET RUN ON DISK.** The
> maintainer ruled that the pass "needs to be included in the cell typing for phenocycler data"
> and scoped it the same day to the join this file describes: `refine.py --out
> data/development/typing/refine`, then `build_cells.py --include-marginal --refine
> data/development/typing/refine`, which adds `refined_class` and the pass's per-cell columns
> **beside** `cell_class`. The islet chain — `cell_class`, `is_endocrine`, `endocrine_route`,
> the CD99 / elongation / hormone routes — is not touched. **[ON DISK]** the join has not run
> until the DATA step, so no cell in `typing/cells/` carries `Unresolved (no lineage evidence)`
> or a filled `refined_*` column yet, and every builder still reads `cell_class`. The
> 5th-percentile floor of section 6, which demotes about 5 % of every class by construction, was
> not revisited by the ruling and is unchanged.

`RECONSTRUCTION.md` documents the subpanel classification and, until 2026-09-18, **none of this**,
which is how the refinement pass got lost: it ran on 2026-09-11, its outcome was written to
`data/development/typing/umap/refinement_log.json`, and no producer for that file existed in the
repository. `refine.py` beside this file is that pass, recovered from the session transcript and
made maintainable. `RECONSTRUCTION.md` now carries a short section pointing here; this file stays
the full record of what the pass does, what was measured while making it run on 25,230,109 cells
instead of 320,000, and what was rejected.

Nothing here is a verification target. Every measurement below is reported because it qualifies a
decision in the module, not because the pass is expected to reproduce it.

## 1. What the pass does

    concordance(cell) = share of its K=30 nearest neighbours IN MARKER SPACE with the same label

`pynndescent.NNDescent(Z, n_neighbors=31, metric="euclidean", random_state=0, n_jobs=-1)`, self
dropped, `LOW = 0.40`. `Z` is the scaled marker matrix. From the recovered script's own docstring:

> Scatter is DETECTED in the embedding (that is where it is visible) but every cell is
> CHARACTERISED and acted on in marker space, so a rule is a statement about measurements and the
> tidier picture is a consequence, never the target.

`refine.py` builds no embedding at all — nothing downstream of it needs one.

    bound_X = 5th percentile of max(OWN markers of X) over the cells of X with concordance >= 0.40
    a cell of X whose max(OWN markers of X) falls below bound_X is DEMOTED

| class | own markers |
|---|---|
| Nerve | B3TUBB, CD56 |
| Lymphatic | Podoplanin |
| Blood vessel | CD31, Caveolin, CD34 |
| Beta / Alpha / Delta | INS / GCG / SST |
| Immune | CD3e, CD68, Iba1, CD20, CD79a, MPO, CD66, CD11c, CD209 |
| Exocrine | Ker8_18, E_cadherin, Pan_Cytokeratin |

The principle, verbatim: membership of a classified object — a centroid inside a nerve, lymphatic
or vessel object — and the islet prior are evidence of WHERE a cell is; they are not evidence that
the cell IS that thing when its own marker is absent.

**Where a demoted cell goes.** Back into the complement it came out of. A demoted Exocrine cell has
no complement beneath it and becomes `Unresolved (no lineage evidence)`. Anything else loses only
the override and is then `Immune` if it satisfies an immune leaf combination, else `Exocrine`.
Concordance is recomputed on the new labels, over the same neighbour index — the representation
does not change when a label does.

**`immune_leaf`** is what makes that routing correct. `rules.requirements(raw, with_negatives=True)`
over the `Imm-` / `T-` / `Mac-` leaves at `REQUIRE_MIN = 0.10`, per donor, over
`qc_analysis_eligible` cells, joined on `object_id`. Without it every demoted region cell went to
Exocrine and none could fall back to Immune; on the recovered subsample adding it moved 5,131 cells.
`REQUIRE_MIN` is imported from `run_step3`, not restated, and is not a tunable.

## 2. The marker set is 54, not 21

`REFINEMENT_RECOVERED.md` heads its section "The marker set — 21, derived not listed". **The
recovered run's own console output says 54**, twice: `saved 54 marker names`, `step 0: 320,000
cells, 54 markers`, and `used in the UMAP: 54` against a 60-marker registry. Re-deriving it on
2026-09-18 returns 54, with the same three drops:

    CD38, b_Catenin1        panel variants -- not acquired on all 20 donors
    CD57                    UniFORM refused the shift (reagent change) -- DERIVED from
                            shifts.json provenance, never named in the module

Six registry markers are absent: those three plus Collagen_IV, DAPI and IAPP, which carry no
non-zero entry in the rule table and so never enter the matrix. The "21" is an error in the
recovered prose, not in the recovered code. `derive_markers()` derives the list and
`test_refine.py` asserts CD57 leaves through provenance rather than by name.

## 3. Reference levels: per donor AND per size stratum

`bound_X` is a reference level, so the 2026-09-18 ruling applies to it: **per donor, and per islet
size stratum.** The strata are `subpanel_io.SIZE_STRATA` — `1`, `2-3`, `4-10`, `>10` cells per islet
— read through `subpanel_io.size_stratum`, so there is one authority for them. Islet size is counted
exactly as `islet_identity.load` counts it: over the cells the analysis admits, islet core,
non-structural, not `Excluded`. A cell outside an islet core carries `""`, which is the absence of a
stratum and not a fifth one: "the size of the islet this cell is not in" is not a number.

**Measured, three donors, per-stratum against per-donor-only:**

| donor | measured cells | demoted per stratum | demoted per donor only | labels differing | of those, islet-core |
|---|---|---|---|---|---|
| 6450 | 363,762 | 24,901 | 24,990 | 289 (0.079 %) | 216 |
| 6539 | 363,177 | 25,254 | 25,284 | 483 (0.133 %) | 289 |
| 6505 | 2,136,234 | 147,301 | 147,488 | 1,116 (0.052 %) | 746 |

So stratifying changes the fate of 0.05–0.13 % of cells, and 60–75 % of those cells are islet-core
although islet-core is a few per cent of the donor. The effect is small in total and concentrated
exactly where the ruling is about. Median |bound shift| is 0.039–0.079 SD; the largest shifts are
all in the smallest groups — Exocrine `4-10` on 6450 moves 1.268 SD on 58 clear members, Immune
`>10` on 6539 moves 1.193 SD on 110. Two shifts are large and not small-group noise: Exocrine's
bound rises monotonically with islet size on 6505 (−1.300 in tissue, −1.068, −0.916, −0.548 in
`>10`), and Delta's moves 0.41 SD between `4-10` and `>10`.

**An under-powered group gets no bound.** Fewer than `MIN_CLEAR = 50` concordant members and the
bound is UNAVAILABLE: not borrowed from another stratum, not replaced by the donor's, not defaulted.
Those cells keep their class and the reason is recorded per cell. Demotion requires positive
evidence that a cell falls below its own class's clear members, and an unfittable group supplies
none. On 6539 that is 23 of 40 groups holding 220 cells in total; on 6505, 13 groups holding 629.
`MIN_CLEAR` is the recovered run's own guard and is deliberately not `subpanel_io`'s 500 — the log
also reports `n_groups_under_500` so a group that clears 50 but not 500 is visible rather than
assumed adequate (6450: 10 of 18 fitted groups).

## 4. Neighbours are drawn within one donor

The recovered run pooled 20 donors into one 320,000-cell graph. `refine.py` builds one graph per
donor. **Measured on a 16,000-per-donor pooled graph (320,000 cells, 54 markers), which is the
recovered run's own geometry:**

* **35.2 %** of a cell's 30 neighbours come from its own donor, against **5.0 %** at chance. Marker
  space is strongly donor-structured after alignment — a 7x enrichment.
* A neighbour from another donor agrees with the cell's label **less often**, in every class but
  one:

| class | own-donor neighbours | same label, own donor | same label, other donor |
|---|---|---|---|
| Exocrine | 34.2 % | 96.3 % | 95.6 % |
| Immune | 35.5 % | 75.2 % | 64.8 % |
| Blood vessel | 45.5 % | 82.1 % | 71.1 % |
| Alpha | 47.1 % | 75.5 % | 67.0 % |
| Beta | 46.2 % | 89.2 % | 79.2 % |
| Lymphatic | 28.1 % | 24.8 % | 29.9 % |
| Delta | 51.9 % | 58.9 % | 49.8 % |
| Nerve | 22.1 % | 44.5 % | 38.9 % |

So in a pooled graph a cell's concordance is decided in part by another donor's staining, and
acquisition batch is aliased with disease status at the coarse level. Within-donor removes that
without a correction. It also removes an imputation: availability is a per-donor fact, and a
within-donor graph uses exactly the markers that donor carries. The recovered pooled matrix ended
with `np.nan_to_num(Z, nan=0.0)` — a marker a donor never acquired entered the distance at the
cohort mean. A missing measurement is unavailable, not a number; `to_Z` refuses a non-finite value
and the caller drops the cell and reports how many.

**The representation stays cohort-wide** — one arcsinh cofactor (each marker's own measured
background / 5) and one standardisation per marker, from a seeded equal-per-donor sample. That is
the opposite decision to the bound and deliberately so, for the reason `run_step3.build_scale`
already gives: a per-donor representation would re-impose the distribution matching UniFORM's step 2
avoided. The representation says what a marker's units are; the bound says how much of it this
donor's clear members show.

**Difference this makes to the numbers.** Within-donor concordance is higher than pooled, as the
table above predicts: scatter before the pass is 4.05 % (6450), 4.90 % (6539), 5.05 % (6505)
within-donor, against 5.25 % on the pooled 320,000-cell graph and 11.51 % in the recovered run. The
recovered run's 11.51 % is not comparable to either: its labels were built in that session from
object membership, the immune leaves and a hormone argmax over islet-core cells, four days before
`typing/cells/` was written. No number in this file is a reproduction of one in the recovered log,
and none should be.

## 5. Cost, measured

| step | scope | wall clock |
|---|---|---|
| representation parameters | 20 donors, 200,000 cells each, 54 markers | 14 s |
| one donor, end to end | 6539, 363,177 measured | 22.4 s |
| one donor, end to end | 6450, 363,762 measured | 23.7 s |
| one donor, end to end | 6505, 2,136,234 measured — the largest | 55.4 s |
| neighbour graph alone | 6539 | 17.9 s |

Twenty donors at that rate is roughly ten minutes serial, and the donors are independent. The
graph is reproducible: `NNDescent(random_state=0, n_jobs=-1)` returned a bit-identical index on two
runs over 200,000 cells, and `test_refine.py` asserts it on a small matrix.

## 6. A 5th-percentile bound cuts 5 % whether or not anything is wrong

**This is the most important measurement here and the recovered run did not make it.** `bound_X` is
the 5th percentile of the class's own concordant members, so a class with no low tail at all still
loses about 5 % of them. The floor is `0.05 x n_clear / n`; the excess above it is what the
scattered members contribute. Measured per class:

| class | 6539 floor / measured / excess | 6450 | 6505 |
|---|---|---|---|
| Exocrine | 4.95 / 5.16 / **+0.21 pp** | 4.96 / 5.32 / **+0.36** | 4.91 / 5.41 / **+0.50** |
| Blood vessel | 4.47 / 10.41 / +5.94 | 4.22 / 12.14 / +7.92 | 4.84 / 5.71 / **+0.87** |
| Beta | 4.74 / 6.95 / +2.21 | 4.61 / 7.60 / +2.99 | 4.51 / 6.18 / +1.67 |
| Alpha | 4.30 / 12.91 / +8.61 | 4.43 / 9.41 / +4.98 | 4.73 / 7.26 / +2.53 |
| Nerve | 4.02 / 14.72 / +10.70 | 2.44 / 27.38 / +24.94 | 3.54 / 15.48 / +11.94 |
| Immune | 3.50 / 18.53 / +15.03 | 3.70 / 16.91 / +13.21 | 4.42 / 10.91 / +6.49 |
| Delta | 3.88 / 12.90 / +9.02 | 3.84 / 22.47 / +18.63 | 3.54 / 20.47 / +16.93 |
| Lymphatic | 3.47 / 18.63 / +15.16 | 3.40 / 23.63 / +20.23 | 1.23 / 24.39 / +23.16 |

Read it this way: Lymphatic, Nerve, Immune and Delta are demoted far above the floor, and that is
the rule finding cells without their own evidence. **Exocrine is not.** Its excess is 0.2–0.5
percentage points, so essentially the whole Exocrine demotion — 15,265 cells on 6539, 16,247 on
6450, 82,522 on 6505 — is the percentile relabelling a fixed share of the class
`Unresolved (no lineage evidence)` by construction. Extrapolated, that is of the order of a million
cells cohort-wide whose new label carries no evidence that their keratin is absent, only that they
are in the bottom 5 % of a class that is 80 % of the tissue. Blood vessel on 6505 (+0.87 pp) is in
the same position; on 6450 and 6539 it is not.

The recovered run had exactly this property — Exocrine 6.0 % demoted against a ~4.9 % floor — and
nobody noticed, because the floor was never computed. It is a decision for the maintainer, not one
this module makes: the module implements the recovered rule and reports the floor in every log
(`groups[].n`, `groups[].n_clear`, `groups[].n_demoted` are all recorded so the floor can be
recomputed from any log). Two obvious alternatives — an absolute level rather than a percentile, or
a percentile only where the class's scattered share is above the floor — are both redesigns and
neither is adopted here.

## 7. Two more properties of the recovered rule, measured

**The routing is a no-op for some demoted Immune cells.** `route` sends a demoted non-Exocrine cell
to `Immune` when `immune_leaf` is true, with no exclusion for cells that were already Immune. So an
Immune cell demoted for carrying none of its own nine markers above the bound is relabelled
`Immune`. Measured: 267 cells on 6539, 2,930 on 6450, 34,163 on 6505 — 23 % of everything demoted on
6505. They are recorded with `demoted = True` and an unchanged `refined_class`, so they are
countable, but the demotion has no effect on them. This is the recovered rule verbatim; whether an
already-Immune cell should instead fall to `Unresolved` is an open decision.

**`immune_leaf` prevalence varies 3x across donors** — 6.4 % of measured cells on 6450, 7.0 % on
6539, 20.9 % on 6505 — and so does the share of demoted region cells it rescues: 5.8 %, 37.9 %,
57.1 %. The flag is computed per donor against the frozen cohort scale
(`broadtype/scratch/cosine_scale.parquet`), which is the scale the classification itself was made
on; the recovered run recomputed the scale with `run_step3.build_scale`, whose seeded 200,000-cell
sample would drift from it. Reading the frozen scale is a deliberate departure and the reason
`RECONSTRUCTION.md` gives for reading it elsewhere.

**Scatter does not necessarily fall.** Per class on 6539 every class improves — Immune 29.9 → 21.3 %,
Lymphatic 31.0 → 21.5 %, Blood vessel 10.6 → 6.1 %, Alpha 14.8 → 9.4 % — but Exocrine worsens
(1.1 → 1.8 %) because demoted cells from other classes land in it, and `Unresolved` is 26.6 %
scattered by construction. Totals: 4.90 → 4.87 % on 6539, 4.05 → 4.68 % on 6450, 5.05 → 5.69 % on
6505. The recovered run's total fell (11.51 → 10.44 %). **A rising total is not a failure of the
rule**: a tidier picture was never the target, and `Unresolved` is by definition the cells with the
least evidence. It is recorded here so no one reads the recovered 11.51 → 10.44 as a property the
pass has to keep.

## 8. Step 2, tested and rejected — and why the module cannot loop

Recorded on the recovered run and not re-litigated here:

> Re-applying the step-1 rule to an already-cleaned population is a ratchet, not a new cause. The
> cells it removes carry 3+ SD of their own defining marker and lie above the Exocrine 95th
> percentile in 91-100 % of cases (Nerve +3.00 SD 100 %, Delta +3.52 100 %, Alpha +3.31 100 %,
> Beta +3.60 100 %, Blood vessel +1.37 91 %), for a 0.12 percentage-point gain in scatter
> (10.44 % → 10.32 %). Only Lymphatic's removed cells were genuinely marginal (median Podoplanin
> +0.86 SD, 0 % above Exocrine p95). Rolling back; the refinement converged at step 1.

Three hypotheses were falsified before it was rejected. H1, that classes defined by one marker stay
scattered, fails: Beta has one marker and 6.5 % scatter, Lymphatic has one and 81 %; it is
separation, not count. H2, that the residue is bright in the non-specific channels, holds only for
`Unresolved`, whose scattered members score +0.96 SD on TP63/FOXP3/PCNA/SOX2/CD107a/Bcl-2/b-actin/
CD44. H3, that Podoplanin is bimodal in Lymphatic, fails: it is a smooth gradient (p5 +0.79 to p95
+3.53) and the class sits +1.47 SD above Exocrine, so a second pass would cut real cells.

The floor measured in section 6 is the mechanism behind that record: a percentile bound removes
about 5 % of a class every time it is applied, so iterating it converges on nothing and keeps
cutting.

**The module encodes the single pass rather than documenting it.** There is no `step` argument
anywhere — `test_refine.py` asserts that no public function takes `step`, `passes`, `iterations` or
`prev`. `apply_once` and `refine_donor` raise `RefinementAlreadyApplied` when handed labels that
already carry `Unresolved (no lineage evidence)` or a frame that already carries `refined_class` /
`refine_pass`. A caller that wants a second pass has to delete the guard, and then has to answer
the record above.

## 9. What the pass writes, and what it does not

    python core/multiplex5000/scripts/typing/refine.py --out <dir> [--donors 6450] [--params <file>]

`--out` is required and has no default, so the pass cannot write to a stage root by being run
without arguments. Per donor it writes `donor_id=<d>/data.parquet`, a sidecar keyed on `object_id`
carrying every row of that donor's typing in the typing's own order:

    object_id  donor_id  cell_class  refined_class  size_stratum  islet_cells  immune_leaf
    concordance  concordance_before  demoted  refine_note  refine_pass

It **does not rewrite `typing/cells/`**. A cell the pass could not measure keeps its class and
carries the reason in `refine_note` — a table overwritten in place could no longer say which cells
those were. `refine_log.json` carries the bounds, the group sizes, the routing and the before/after
composition per donor.

Cells the pass does not act on, and why:

* `cell_class == "Excluded"` — the `exocrine_inside_islet_border` set. An excluded cell is excluded
  from every calculation, so it is neither demoted nor allowed to be another cell's neighbour.
* not `qc_analysis_eligible` in **either** the typing table or the aligned table. On the seven
  re-exported donors the two disagree: 3,718 cells on 6450 (1.0 %) are eligible in the revised
  typing and ineligible in `celltype/uniform/`, which was aligned before the re-export. 6539, which
  was not re-exported, has no disagreement. The pass takes the intersection and reports the count
  rather than choosing a winner.

## 10. Open decisions this hands to the maintainer

Whether the pass runs at all is no longer among them: it was ruled in on 2026-09-23, through the
`--refine` join beside `cell_class`. What stays open is how the refined verdict is used. Items 1-3
below matter to any builder that reads `refined_class` rather than `cell_class`, and none does
today; pointing one at it needs decisions that touch the islet chain (a demoted endocrine cell's
`is_endocrine`; a demoted cell inside an islet outline routed to Exocrine), which are asked
specifically, never settled as part of the wiring (`docs/PRODUCERS.md`, "Not decided here").

1. **The Exocrine floor** (section 6). ~5 % of Exocrine becomes `Unresolved` by construction.
2. **A demoted islet-core endocrine cell becomes Exocrine**, and `aggregate_islets.py` drops the
   exocrine-inside-an-islet-border set as `Excluded` — so it leaves the denominator rather than
   becoming `Unresolved`. Measured: 1,299 cells on 6450, 650 on 6539, 5,368 on 6505, recorded in
   every log as `demoted_islet_core_endocrine_to_exocrine`. This is coherent with the
   classification's own chain, which also terminates at Exocrine inside an islet, and with the
   maintainer's rule that exocrine markers only reach islet annotations by mistake — the cell is
   dropped, not claimed as exocrine biology. It should be explicit downstream rather than incidental.
3. **The no-op Immune demotion** (section 7).
4. **`Unresolved (no lineage evidence)` is a new class value.** Two of the three consumers have
   since been told, on 2026-09-18: `build_cells.py`'s `DISPLAY` map carries it (→ `Unresolved`)
   and the app's `run20_default.yml` enumerates it as a ninth entry that is an absence of a call,
   not a ninth cell type. **`trajectory_singlecell_features.py` HAS now been told, 2026-09-19** —
   it loads `refine.OWN_MARKERS` by file path and no longer reads the 2026-09-11
   `umap/refinement_log.json` at all; that file now refuses to be read. Verified at runtime that
   the two names are the same object, so the table cannot drift into a copy.

   What is NOT closed, and is reported rather than fixed: 13 of that module's 38 markers still
   come only from `subpanel/subpanel_log.json`, which has no producer either — 688 bytes, written
   2026-09-11 14:10, 32 minutes after the orphan just closed. None of those files belongs to
   this module.

   Being told is not being wired: no cell carries the value today, because the pass has not been
   joined into a build. See the status banner at the top of this file.

## October 8 implementation of the CD8 refinement ruling

The October 3 `results/v3/release_20261003/production_calibration_ruling.json`
specifies shifted CD8 for primary pancreas typing and unshifted CD8 for refinement.
The historical selected arm did not implement that distinction. `refine.py` now
accepts `--cd8-alignment unshifted --shifts <matching-arm-shifts.json>`; the
calibration runner exposes `--refinement-cd8-alignment unshifted --shifts ...`.
Comparison runs retain their declared aligned-input behavior unless this option is
selected. This pancreas policy does not alter the dedicated pLN calibration.

The shared refinement reader reverses only CD8's recorded donor multiplier. It
preserves exact UUIDs, missing measurements, zeros, other channels and the primary
aligned artifact. The restored values enter fresh representation fitting and the
neighbor calculation. Immune fallback reads the restored measurement and converts
the frozen CD8 floor/top/span by the same inverse factor: measurement and reference
must use the same units. This preserves the frozen criterion instead of introducing
an uncalibrated threshold after changing units.

Unshifted mode refuses `--params`, so an aligned-CD8 representation cannot be reused.
Every refinement run saves `representation_params.parquet` and a provenance record
with its parameter hash, fit donors, source script, shifts, scale and CD8 policy.
Use a fresh output directory with the intended full donor cohort. The arm runner
rejects existing cell partitions outside the declared donor set. Before unshifting,
verify recovered CD8 against the exact pre-alignment expression keys and block mask;
recording a shifts path alone cannot prove it belongs to the aligned input.

The October 8 verification recovered CD8 for 2,048 exact keys in each of all 25
pancreas donors, matching the pre-alignment expression and recorded block mask.
The receipt is
`results/v3/islet_immune_reassessment_20261006/cell_typing_validation/provenance/CD8_refinement_input_parity_20261008.json`.
Focused tests passed 39 cases, with one data test deselected. This verifies the
measurement contract and its safeguards; it does not establish biological typing
accuracy or validate previously generated refinement outputs.
