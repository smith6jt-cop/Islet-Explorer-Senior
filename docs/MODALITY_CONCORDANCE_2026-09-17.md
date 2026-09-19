# Protein vs RNA lineage concordance, measured 2026-09-17

Produced by `scripts/senior/modality_concordance.py`. Every number below came out of one run of
that script; tables are in `data/development/integration/concordance/`.

Development lane. Nothing here is a release, and `NO_VERIFIED_TYPING` stays down.

> ## CORRECTION BANNER — 2026-09-18: every pancreas number below is off a superseded source
>
> **The pLN half of this document is unaffected.** What follows corrects the pancreas half only.
>
> **1. The 21-donor pancreas release exists and this document never reads it.** It is
> `xenium-panc-cluster-name-20260917`, written **2026-09-17T20:36:54+00:00** to
> **`/data/xenium_typing/panc_2026-09-17/release`** — outside the repository, which is why a
> repo-relative search never saw it. Measured from its own `release.json`: **21 donors, 21
> sections, 17,490,305 cells, 15,318,868 called, 775,234 excluded for no nucleus**, with 21
> `donor_id=` partitions under `calls/`. Before this banner the string `panc_2026-09-17` appeared
> nowhere in this file.
>
> **2. The run behind this document predates it by 49 minutes.**
> `data/development/integration/concordance/manifest.json` is stamped
> `generated_at_utc = 2026-09-17T19:47:22Z`. The release was written at 20:36:54Z. Every pancreas
> row below — the two-donor table, the `pancreas_*.csv` files, the empty eta-squared table — came
> out of that earlier run and describes `xenium-typing-scanvi-per-donor-dev-20260907`, the
> **superseded** lane, which the manifest names in `pancreas.rna_notes`. Read the two-donor table
> as a record of what a superseded lane said, not as the state of the pancreas comparison.
>
> **3. `release_panc_v1` will never land, so "re-run once it lands" is not a plan.** Measured:
> `data/development/xenium_typing/release_panc_v1` contains exactly four files —
> `_work/6539.calls.parquet`, `_work/6539.clusters.parquet`, `_work/6539.gene_stats.parquet` and
> `_work/6539.report.json`. One donor's scratch. **No `calls/` tree and no `release.json`**, so it
> declares no release id, no tissue and no write time; it cannot be selected by anything that
> reads a release's own metadata, and nothing is writing to it. The two sentences below that wait
> on it are waiting on something that is not coming.
>
> **4. The script has already changed, so "the script needs no change" is also wrong.** The
> pancreas pin was removed in commit `b57b48e` (2026-09-18). `modality_concordance.py` no longer
> holds a literal fallback chain; it calls `scripts/senior/xenium_release.resolve(tissue)`, which
> searches both roots and selects on each release's **own declared metadata** — never a filesystem
> mtime — and honours a declared `supersedes`. Run today it resolves:
>
> | tissue | resolver returns |
> |---|---|
> | `panc` | `xenium-panc-cluster-name-20260917` at `/data/xenium_typing/panc_2026-09-17/release` |
> | `pln` | `xenium-pln-immune-tnk-split-dev-20260917` at `data/development/xenium_typing/release_pln_immune_v3` |
>
> Note the second row: the pLN source named in this document, `release_pln_immune_v2`, has **also**
> been superseded, by `release_pln_immune_v3`. The pLN numbers below were correct for v2 when they
> were taken; they are not a reading of the current pLN release either.
>
> **5. What a re-run would have to work with, measured — and what it is not.** The protein
> pancreas side (`data/development/typing/cells`) holds 20 donors; the RNA release holds 21;
> **16 are in both**: 115, 6374, 6380, 6436, 6442, 6461, 6476, 6482, 6493, 6521, 6523, 6539,
> 6551, 6566, 6591, 6593. Protein-only: 6450, 6505, 6534, 6623. RNA-only: 112, 6347, 6533, 6579,
> 6580. 16 is past the four-donor floor, so the "fewer than 4 donors → no verdict at all" rule
> would no longer fire on the pancreas. **That licenses nothing.** No pancreas quantity has been
> recomputed here, no threshold has been applied to the new release, and **no pancreas lineage may
> be read as carrying or as refused a cross-modality claim on the strength of this banner.** The
> pancreas verdict is *unmeasured*, which is not the same as the "undetermined on two donors"
> stated below and is certainly not concordance.
>
> One thing a re-run inherits, recorded here so it is not discovered late:
> `docs/XENIUM_GDNA_TRIAGE_2026-09-17.md` measures `0080870` panc (donor **6533**) at 0.26100,
> **57.5x** the tuned rate, and `0081366` panc (donor **112**) at 0.30600, **67.4x**; both are
> typed in the 21-donor release with no gDNA correction applied and no flag on the rows
> (`docs/PANC_XENIUM_TYPING_2026-09-17.md` §Depth). Neither donor is on the protein pancreas side,
> so neither would enter the 16-donor comparable set — but both are inside any cohort-pooled RNA
> number taken off that release.

## What this is and what it is not

PhenoCycler and Xenium are **different measurements of different sections of different blocks**.
They are not two views of one cell. The comparison is therefore **compositional only**, joined on
`donor_id` and nothing else.

**No cell identifier crosses.** A PhenoCycler `object_id` is a QuPath UUID; a Xenium `cell_id` is
10x segmentation and unique only within a bundle. Nothing here joins them, and nothing here computes
a distance, an adjacency, a neighbourhood or a shared scatter across the modalities.

Labels are harmonised through `core/phenocycler/data/integration/vocab_crosswalk.csv`, which is
generated by `run_vocab` in `phenocycler/integration/vocab.py` and was read, never edited.

Neither side was tuned to improve the other. Agreement is informative and disagreement is
informative; which lineages may carry a cross-modality claim is the deliverable.

## The comparable set, and what the crosswalk keeps out

A lineage is comparable only when the crosswalk names it on **both** sides and neither side's
resolution disqualifies it.

- **Granulocyte is out.** `xenium_resolution = fine_only`: it has no broad class on the RNA side, so
  a Granulocyte share is not the same quantity in the two modalities. It is counted and reported,
  never compared.
- **Neural is out on the pLN and on the pancreas protein side** (`pheno_resolution = no`).
- **Stromal is out on the pancreas**: the PhenoCycler pancreas typing produces no Stromal call at
  all, so the protein share is **unavailable, not zero**.
- Two asymmetries are carried, not hidden: PhenoCycler `Epithelial_nonendocrine` is a **default
  sink** mapped to Exocrine (`epi_default` fallback when all positive gates fail), and Xenium
  Endocrine is a **surrogate** — there is no 1:1 marker, it is resolved through a curated identity
  panel. Both are marked in `manifest.json` under `asymmetries`.

**pLN comparable set: T_NK, B_Plasma, Myeloid.**
**Pancreas comparable set: Endocrine, Exocrine, Vascular, T_NK, B_Plasma, Myeloid.**

Shares are computed within the comparable set in each modality separately, so the two denominators
are each modality's own comparable cells. Both denominators are printed per donor in
`*_per_donor_per_lineage.csv` so a reader can see they differ.

## The disagreement threshold — fixed before the numbers

These thresholds are constants at the top of `scripts/senior/modality_concordance.py`
(`THRESHOLD_*`), written before the script was first run and unchanged since.

A share is a part of a composition over K comparable lineages, so the mean share is 100/K
percentage points — 33.3 pp for the three-part pLN composition. Some disagreement *must* be
sectioning and not method: which part of a node a section cuts (cortex, paracortex, medulla,
capsule) moves a B-versus-T share, and no modality can be blamed for that. What sectioning does not
do is reorder the composition.

| criterion | threshold | why |
|---|---|---|
| `abs(median over donors of (protein share − RNA share))` ≤ **10 pp** | concordant | under a third of the mean share; inside the range regional sampling of one node plausibly produces |
| 10 pp < median ≤ **20 pp** | grey band | usable only with the offset stated alongside every number |
| median > **20 pp** | structural | 60% of the mean share; a shift that large changes which lineage dominates |
| IQR of the per-donor difference > **20 pp** | structural | the offset is donor-dependent, and a donor-dependent offset cannot be carried by a pooled claim even when it cancels in the median |
| the lineage sits at a different rank position in the two modalities in > **50%** of donors | structural | the compositions are ordered differently, whatever the per-lineage medians say |
| fewer than **4 donors** | no verdict at all | a median and an IQR over two or three donors is not a measurement of a cohort |

One implementation note, stated because it is a change made after seeing output: the rank criterion
as first written compared the *whole* ordering of the comparable set, which yields one number for
the whole composition and therefore tags every lineage identically. The magnitude thresholds were
not touched; the rank criterion is now read per lineage (does *this* lineage hold the same rank
position in both modalities), and **both** numbers are in the summary CSV —
`rank_disagreement_frac` (per lineage, used by the verdict) and
`whole_composition_rank_disagreement_frac` (context).

## Sources, and who was available when this ran

| side | tree | donors present |
|---|---|---|
| pLN protein | `data/development/typing/pln/pheno` | 14 typed |
| pLN RNA | `data/development/xenium_typing/release_pln_immune_v2` | 20 called |
| pancreas protein | `data/development/typing/cells` | 20 |
| pancreas RNA | `data/development/xenium_typing/release` — **superseded; see the banner** | **2** |
| pancreas RNA, current | `/data/xenium_typing/panc_2026-09-17/release` — **not read by this run** | 21 called |
| pLN RNA, current | `data/development/xenium_typing/release_pln_immune_v3` — **not read by this run** | see banner |

**pLN: 11 donors are in both** — 115, 6380, 6436, 6442, 6461, 6476, 6482, 6493, 6521, 6523, 6539.
Not 18. The protein side has 14 typed sections: 6551, 6566, 6593 and 6623 are **unavailable** (their
QuPath export carries no intensity measurement, so there is nothing to cluster — never zero immune
cells) and 6579's node is **withheld** by `pln_node_annotations.csv`. Of the 14 typed, 6450, 6505
and 6534 have no Xenium pLN bundle in the v2 tree. 6457 is centrally excluded, read from
`DONOR_EXCLUSIONS` in `core/phenocycler/phenocycler/cohort.py` rather than from a literal.

**Pancreas: 2 donors, 6436 and 6442.** The only pancreas Xenium calls this run read are
`xenium-typing-scanvi-per-donor-dev-20260907`, a **superseded lane**. The pancreas numbers below
are observations on two donors from an old lane and carry **no verdict**.

> **CORRECTED 2026-09-18.** Two sentences stood here: "`release_panc_v1` — the fill-panc-rna track
> — did not exist when this ran" and "The script probes `release_panc_v1` first and will use it
> without modification once it lands." Both are wrong now, in opposite directions.
> `release_panc_v1` **does** exist — as `_work/6539.*` and nothing else, no `calls/` and no
> `release.json` — and it will not fill. And the script no longer probes it: the literal pin was
> removed in `b57b48e` and replaced by `scripts/senior/xenium_release.resolve()`. What did exist
> when this ran, 49 minutes later on the same day, is the 21-donor release at
> `/data/xenium_typing/panc_2026-09-17/release`. See the banner at the top.

### Nucleus handling

A cell whose nucleus was **measured absent** is excluded from every calculation. Six of the twenty
Xenium pLN sections (6436, 6476, 6493, 6523, 6539, 6579) carry `nucleus_count = -1` for every cell —
the nucleus was never measured there. That is a **missing measurement, not a measured zero**, and
dropping those cells would have discarded six whole sections on the strength of an absent
measurement. They are retained and the section is flagged by `nucleus_measured_frac` in every
output. Where the measurement does exist, measured-absent nuclei are 0.4–1.4% of cells, so the
retained sections cannot be carrying an error of the order of the thresholds above.

The pancreas Xenium release has no nucleus column at all; the filter could not be applied there and
that is recorded per donor.

### Controls (`pln_rna_controls.csv`)

Y-probe detection, which is a control only on donors with no Y template: the six female donors run
**0.027%–0.096%** (median 0.040%). The nine male donors run 12.2%–50.7% and track depth, as a real
template should. Detection is `lognorm >= log1p(1)` throughout, applied upstream; no raw-count floor
was introduced here.

Transcript depth is the flag that matters, and it is a flag on depth, not on gDNA: 6566 at a median
of **78** transcripts per cell, 6380 at **97** and 6580 at **67** are weak substrate for a 5K-panel
call. 6566 and 6580 are not in the compared set; **6380 is**, and its RNA composition should be read
with that in mind.

## The pLN result — 11 donors

Per-donor shares within {T_NK, B_Plasma, Myeloid}, protein / RNA / difference in percentage points:

| donor | status | T_NK prot | T_NK RNA | Δ | B_Plasma prot | B_Plasma RNA | Δ | Myeloid prot | Myeloid RNA | Δ |
|---|---|---|---|---|---|---|---|---|---|---|
| 115  | ND   | 12.2 | 55.1 | −42.8 | 35.0 | 33.9 | +1.0 | 52.8 | 11.0 | +41.8 |
| 6380 | T1D  | 43.1 | 51.5 | −8.4  | 25.8 | 33.2 | −7.4 | 31.1 | 15.3 | +15.8 |
| 6436 | T1D  | 39.5 | 64.1 | −24.6 | 19.0 | 19.2 | −0.3 | 41.5 | 16.7 | +24.8 |
| 6442 | ND   | 31.1 | 50.5 | −19.4 | 30.0 | 32.8 | −2.8 | 38.9 | 16.7 | +22.2 |
| 6461 | ND   | 35.5 | 40.9 | −5.4  | 29.9 | 44.9 | −15.0 | 34.6 | 14.2 | +20.3 |
| 6476 | ND   | 19.9 | 39.6 | −19.7 | 36.8 | 52.7 | −15.9 | 43.3 | 7.7 | +35.6 |
| 6482 | ND   | 28.4 | 60.9 | −32.4 | 21.5 | 18.8 | +2.7 | 50.1 | 20.3 | +29.7 |
| 6493 | ND   | 33.3 | 56.0 | −22.7 | 27.9 | 35.9 | −7.9 | 38.7 | 8.1 | +30.6 |
| 6521 | Aab+ | 21.8 | 48.0 | −26.2 | 26.2 | 32.2 | −6.0 | 52.0 | 19.8 | +32.1 |
| 6523 | T1D  | 34.9 | 50.8 | −15.9 | 22.5 | 37.0 | −14.5 | 42.5 | 12.2 | +30.4 |
| 6539 | ND   | 24.7 | 49.5 | −24.7 | 22.4 | 32.0 | −9.6 | 52.8 | 18.5 | +34.3 |

Summary, against the thresholds fixed above:

| lineage | n | median prot | median RNA | median Δ | IQR Δ | range Δ | rank disagree | verdict |
|---|---|---|---|---|---|---|---|---|
| **T_NK** | 11 | 31.1 | 50.8 | **−22.7** | 7.8 | −42.8 … −5.4 | 91% | **MAY NOT carry a cross-modality claim** |
| **B_Plasma** | 11 | 26.2 | 33.2 | **−7.4** | 10.5 | −15.9 … +2.7 | 73% | **MAY NOT carry a cross-modality claim** |
| **Myeloid** | 11 | 42.5 | 15.3 | **+30.4** | 9.7 | +15.8 … +41.8 | 100% | **MAY NOT carry a cross-modality claim** |

**None of the three pLN lineages may carry a cross-modality claim.** Two fail on magnitude — T_NK at
−22.7 pp and Myeloid at +30.4 pp, both past the 20 pp structural threshold. B_Plasma's median offset
of −7.4 pp is inside the concordant band and its IQR of 10.5 pp is inside the dispersion band; it
fails on rank alone, sitting at a different position in the composition in 8 of 11 donors. That is a
weaker failure than the other two and it is worth naming as such: **B_Plasma's magnitude agrees and
its ordering does not.** A statement about the B-cell share of a node that does not depend on its
rank against T and myeloid is the one statement in this table that comes close to surviving, and
even that is not licensed here — it fails the criterion as written, and the criterion was written
first.

The sign is completely consistent: protein calls **more myeloid and less T** than RNA in **11 of 11
donors**, with no donor reversing either sign. That is not sectioning. Sectioning would scatter.

What it is instead, this measurement does not determine. One candidate, stated as a candidate and
not as a result: the two methods do not have symmetric routes to the two labels. On the protein
side, the pLN naming set is 5 markers for Myeloid (CD11c, CD163, CD206, CD68, Iba1) against 3 for
T_cell (CD3e, CD4, CD8) and 2 for B_Plasma. On the RNA side, the identity cores are 129 genes for
T_NK, 111 for Myeloid, 58 for B_Plasma. Testing that would mean measuring each method against
itself, not against the other, and it is not done here.

Two further asymmetries that belong beside the numbers rather than inside them. NK is reachable on
the RNA side (`Immune_NK_ILC`, folded to T_NK by the crosswalk) and **not** on the protein side: on
this 60-plex CD56 is the NK subpanel's only member and `SUPPORTING_ONLY` holds it back from naming a
cell alone, so a count of zero NK clusters on the protein side is a property of the panel and not a
measurement that the node holds no NK cells. That asymmetry pushes T_NK **up** on the RNA side,
which is the direction observed. And donor 6461 has 10,155 protein cells that are immune but
unnamed; they are outside the comparable denominator on that donor, on that side only.

## The pancreas result — 2 donors from a superseded lane, off a superseded source

> **SUPERSEDED 2026-09-18 — read the banner at the top before this section.** Every number in this
> section is off `xenium-typing-scanvi-per-donor-dev-20260907`. The current pancreas release,
> `xenium-panc-cluster-name-20260917`, covers **21 donors** and **16 of them are shared with the
> protein side**; it has not been compared here. These rows are a record of a superseded lane, not
> a statement about the pancreas.

These are `xenium-typing-scanvi-per-donor-dev-20260907`, a superseded lane, on 6436 and 6442 only.
**Observations, not a verdict** — the fewer-than-4-donor rule fires on every row.

| lineage | prot 6436 | RNA 6436 | Δ | prot 6442 | RNA 6442 | Δ |
|---|---|---|---|---|---|---|
| Endocrine | 2.73 | 5.79 | −3.06 | 1.40 | 2.13 | −0.74 |
| Exocrine | 81.53 | 71.68 | +9.84 | 79.53 | 80.50 | −0.97 |
| Vascular | 4.51 | 8.25 | −3.74 | 3.89 | 1.15 | +2.73 |
| T_NK | 0.63 | 7.40 | −6.78 | 1.34 | 1.35 | −0.01 |
| B_Plasma | 2.95 | 0.86 | +2.09 | 3.04 | 3.23 | −0.19 |
| Myeloid | 7.65 | 6.02 | +1.64 | 10.81 | 11.62 | −0.82 |

The differences are an order of magnitude smaller than on the pLN — every one inside 10 pp — but
that comparison is between two donors from a withdrawn lane and eleven donors from the pLN lane
current at the time, and it means nothing. **Corrected 2026-09-18:** this sentence used to end "it
means nothing until `release_panc_v1` fills". `release_panc_v1` is not going to fill (banner, point
3). What would make this table mean something is a re-run against
`/data/xenium_typing/panc_2026-09-17/release`, which is 21 donors and 16 shared with the protein
side, and which has not been run. Note also that 6442's protein T_NK and RNA T_NK
agree to 0.01 pp while 6436's differ by 6.78 pp; with two donors there is no way to tell which is
the ordinary case.

Comparable totals differ widely — 1,517,058 protein against 576,120 RNA on 6436 — because the two
sides call different fractions of their cells into the comparable set, and because the nucleus
filter could not be applied on the RNA side. Both denominators are in the CSV.

## Batch

Acquisition batch is aliased with disease status at the coarse level, so correction is only ever on
the finer variables. Eta-squared was computed on `acq_run` and `stain_group` from
`data/donor_acquisition_runs.csv`, each against a permutation null of 10,000 shuffles of the batch
labels over donors, for **every** cross-donor quantity: each lineage's protein share, its RNA share,
and the protein−RNA difference. Full table: `pln_batch_eta_squared.csv`.

**Read the null before the statistic.** With 11 donors spread over 9 `acq_run` levels, the
permutation null median eta-squared is **0.83–0.89**. At that many groups per donor, a random
partition explains almost all the variance, so eta-squared on `acq_run` cannot discriminate here and
a large value on it is not evidence of anything. `stain_group` has 7 levels and a null median of
0.59–0.64, which leaves some room.

Three quantities clear their null at p < 0.05:

- **T_NK protein share on `stain_group`**: eta-squared **0.921**, null median 0.626, **p = 0.018**.
- **Myeloid protein−RNA difference on `acq_run`**: eta-squared **0.985**, null median 0.841,
  **p = 0.035**.
- The same Myeloid difference on `stain_group` is 0.752 against a null of 0.617, p = 0.251 — it does
  not clear.

Everything else sits at or below its null: T_NK RNA share on `stain_group` 0.723 (p = 0.319),
B_Plasma difference on `acq_run` 0.825 (p = 0.535), Myeloid RNA share on `stain_group` 0.301 against
a null of 0.592 (p = 0.956).

**Nothing is dropped and nothing is gated on this.** That the T_NK protein share carries stain-group
structure at p = 0.018 is a fact to report alongside the T_NK disagreement, not a reason to discard
either. It does mean the T_NK protein share is not a clean cross-donor quantity, which is one more
reason the T_NK cross-modality claim is refused.

The pancreas eta-squared table is empty: two donors is below the three needed for the statistic to
be defined at all. That is unavailable, not zero.

## What may and may not carry a cross-modality claim

**pLN — none of the three.**

| lineage | may carry a cross-modality claim? | the binding reason |
|---|---|---|
| T_NK | **No** | −22.7 pp median offset, past the structural threshold; rank differs in 91% of donors; the protein share additionally carries stain-group structure at p = 0.018; NK is reachable on the RNA side and unreachable on the protein side |
| Myeloid | **No** | +30.4 pp median offset, the largest measured; rank differs in 100% of donors; the difference carries acq_run structure at p = 0.035 |
| B_Plasma | **No**, but on rank only | median offset −7.4 pp and IQR 10.5 pp are both inside the fixed bands; it fails because its position in the composition differs in 73% of donors |
| Granulocyte | **Not comparable** | `xenium_resolution = fine_only`; no broad RNA class, so it is not the same quantity |
| NK, alone | **Not comparable** | unreachable on the 60-plex protein side; a zero there is a panel property |

**Pancreas — unmeasured against the current release.** No lineage may carry a pancreas
cross-modality claim today. As this run stood, the reason was arithmetic: two donors cannot
establish that two compositions agree, whatever their differences look like. **As of 2026-09-18 the
reason is different and stronger — the comparison has not been run against the release that
exists.**

> **CORRECTED 2026-09-18.** This paragraph used to end "Re-run once `release_panc_v1` lands; the
> script needs no change." Neither half holds. `release_panc_v1` holds one donor's `_work/` scratch
> and no `calls/` or `release.json`; it is not going to land. And the script was changed on
> 2026-09-18 (`b57b48e`) — its three literal pins were replaced by
> `scripts/senior/xenium_release.resolve()`, which now names
> `/data/xenium_typing/panc_2026-09-17/release` for the pancreas and `release_pln_immune_v3` for
> the pLN. The re-run is a re-run of **both** halves against sources neither half of this document
> read.

**What this forbids downstream.** Any statement of the form "the myeloid fraction of the node is X"
that mixes or substitutes modalities is refused by this measurement. Any statement that a Xenium
T-cell share corroborates a PhenoCycler one, or the reverse, is refused. What remains permitted is a
**within-modality** cross-donor statement, and a statement of the disagreement itself — the protein
side calls more myeloid and less T than the RNA side in 11 of 11 donors, with no reversal, which is
a property of the two methods worth carrying forward.

## What is not here, deliberately

No distance, adjacency, neighbourhood, KDE or shared scatter crosses the pancreas/pLN boundary. A
lymph node has no islet, so a signed distance to one is undefined there, not zero. The maintainer's
annotation of pLN portions of pancreas slides — where some lymphoid structures are anatomically
attached and an islet-to-lymphoid distance *is* meaningful — has not landed and is not inferred
here. This script is tissue-keyed and takes a tissue argument at every read, so honouring that
dataset when it arrives is adding a source, not reworking the comparison.

## Files

- `scripts/senior/modality_concordance.py`
- `data/development/integration/concordance/manifest.json` — sources, comparable set, exclusions with
  reasons, asymmetries, the thresholds as they were fixed
- `.../pln_per_donor_per_lineage.csv`, `.../pln_lineage_summary.csv`, `.../pln_batch_eta_squared.csv`
- `.../pln_per_donor_counts_protein.csv`, `.../pln_per_donor_counts_rna.csv`, `.../pln_rna_controls.csv`
- `.../pancreas_*` — the same shape, on the two donors available
