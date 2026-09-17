# The four-cell grid: both tissues, both modalities, every donor

Built 2026-09-17 by `scripts/senior/build_grid_composition.py` into
`$ISLET_LAYERS_DIR/grid_composition` (default `/data/islet_layers/grid_composition`), read by
`apps/islet/R/card_g1_grid.R`, tested by `core/phenocycler/tests/current/test_grid_composition.py`
(30 tests, all passing).

## What was wrong

The app read one of four cells. Every app composition artifact under
`data/development/typing/app/` is built off `data/development/typing/immune/`, which covers ten
donors, so the app showed ten. Three other cells of the grid had calls and no reader at all.

## What the data actually resolves — the ten "missing" donors were not missing

Rebuilding the immune lane was not necessary for the broad grid. `data/development/typing/cells`
carries `immune_definition` on **all 20 donors** — the per-marker definition gate the crosswalk
describes as "recovered from per-marker gates, not the broad label". Measured over the whole
25,230,109-cell tree:

| `immune_definition` | cells | common lineage |
|---|---|---|
| Imm-macrophage | 1,574,036 | Myeloid |
| Imm-B-cell | 372,562 | B_Plasma |
| Imm-dendritic-cell | 256,946 | Myeloid |
| Imm-T-cell | 237,066 | T_NK |
| Imm-plasma-cell | 212,247 | B_Plasma |
| Imm-neutrophil | 87,713 | Granulocyte |
| Imm-NK-cell | 38,804 | T_NK |
| *(empty)* | 155,524 of the 2,934,898 Immune cells | `Immune_unresolved` |

That total, 2,779,374, is exactly the `resolved_by == "definition"` count, and the 155,524 residue
is exactly `resolved_by == "supporting"`. `typing/immune` is opened only for the **fine** subtype
(Mac-M2-like, T-CD4-helper, …), which is the ten-donor tree. The broad pancreas protein grid is
20 donors, not ten.

## The layer

`(tissue, modality, donor_id)` keyed, written where its reader looks.

| file | what it is |
|---|---|
| `composition.parquet` / `.csv` | `(tissue, modality, donor_id, bucket_kind, lineage_common)` → `n_cells`, `denom_assessed`, `frac_of_assessed`, `state`, `reason`, the crosswalk's `pheno_resolution` / `xenium_resolution` / `is_immune`, `resolution_conflict`, `disease_status` |
| `fine.parquet` | the fine call inside a common lineage, where a lane has one: `definition` and `subtype` (pancreas protein), `t_nk_class` (pLN protein), `lineage` (both Xenium lanes) |
| `presence.parquet` / `.csv` | one row per donor per grid cell — 25 × 4 = 100 rows, each with one of the five availability states and a reason |
| `crossmodal_tnk.csv` | per-donor T_NK share of resolved immune, both modalities, per tissue, with ranks among donors measured in both |
| `batch_eta_squared.csv` | eta-squared of every per-donor lineage fraction on `disease_status`, `stain_group` and `acq_run`, each against its own permutation null (2000 permutations) |
| `manifest.json` | sources, `source_domain_root`, builder, commit, totals, coherence result |

`bucket_kind` separates a lineage from the book-keeping that makes the parts sum to the whole:
`lineage`, `partial` (`Immune_unresolved` — a positive immune call that stops short of a lineage),
`negative` (`Non_immune` — a measured "not that", which only the pLN lanes make), `unassigned`
(looked at, no call, one row per distinct producer reason) and `withdrawn` (called, then taken
back — the pancreas protein `exocrine_inside_islet_border` exclusions).

## Coverage, measured

| cell | donors | cells assessed | no nucleus, excluded | nucleus never measured |
|---|---|---|---|---|
| pancreas / protein | 20 | 25,230,109 | 4,748,402 | 0 |
| pancreas / RNA | 18 | 13,725,222 | 775,234 | 2,989,849 |
| pLN / protein | 14 | 5,122,269 | 0 | 0 |
| pLN / RNA | 20 | 5,927,317 | 80,472 | 0 |

**A nucleus that was never measured is not a nucleus that is absent.** In the pancreas Xenium
release, donors 6347, 6579 and 6580 carry `nucleus_measured == False` on every cell — 1,393,487,
665,135 and 931,227 cells. Treating that as "no nucleus" would silently exclude three whole donors
under a rule that never applied to them, so those three cells of the grid are `unavailable` with
that reason and contribute nothing to any denominator. Every other donor's `nucleus_count == 0`
cells match its `Unavailable_no_nucleus` label exactly and are excluded as the rule requires.

## Per-donor presence — which of the four cells each donor has, and why not

`yes` = measured. Everything else is a state with a reason, never a blank panel and never a zero.

| donor | status | panc/protein | panc/RNA | pLN/protein | pLN/RNA | cells measured |
|---|---|---|---|---|---|---|
| 115 | ND | yes | yes | yes | yes | 2,205,227 |
| 6374 | ND | yes | yes | blank | yes | 3,242,360 |
| 6380 | T1D | yes | yes | yes | yes | 1,781,877 |
| 6436 | T1D | yes | yes | yes | yes | 3,084,666 |
| 6442 | ND | yes | yes | yes | yes | 2,815,864 |
| 6450 | Aab+ | yes | blank | yes | blank | 761,379 |
| 6461 | ND | yes | yes | yes | yes | 2,462,127 |
| 6476 | ND | yes | yes | yes | yes | 3,465,466 |
| 6482 | ND | yes | yes | yes | yes | 2,400,880 |
| 6493 | ND | yes | yes | yes | yes | 2,208,276 |
| 6505 | Aab+ | yes | blank | yes | blank | 2,432,518 |
| 6521 | Aab+ | yes | yes | yes | yes | 2,814,060 |
| 6523 | T1D | yes | yes | yes | yes | 2,835,294 |
| 6534 | T1D | yes | blank | yes | blank | 1,285,691 |
| 6539 | ND | yes | yes | yes | yes | 1,227,715 |
| 6551 | T1D | yes | yes | unavailable | yes | 2,261,307 |
| 6566 | T1D | yes | yes | unavailable | yes | 3,248,372 |
| 6591 | T1D | yes | yes | blank | blank | 2,181,965 |
| 6593 | T1D | yes | yes | unavailable | yes | 3,444,368 |
| 6623 | T1D | yes | blank | unavailable | blank | 1,216,131 |
| 112 | ND | blank | yes | blank | yes | 738,324 |
| 6347 | Aab+ | blank | unavailable | blank | yes | 78,912 |
| 6533 | T1D | blank | yes | blank | yes | 1,124,476 |
| 6579 | T1D | blank | unavailable | withheld | yes | 569,223 |
| 6580 | T1D | blank | unavailable | blank | yes | 118,439 |

The reasons, as the producers themselves recorded them:

* **pancreas / protein, blank** (112, 6347, 6533, 6579, 6580) — no PhenoCycler pancreas partition;
  these donors reach the cohort through the Xenium lane only. For 6579 the note adds that
  `cohort.py` records it as reinstated on 2026-09-17 and the call tree predates that notice, so the
  absence may be the lapsed exclusion rather than a missing acquisition. Nothing here re-types it.
* **pancreas / RNA, unavailable** (6347, 6579, 6580) — nucleus never measured, above.
* **pancreas / RNA, blank** (6450, 6505, 6534, 6623) — the lane was not run on a pancreas section
  of these donors.
* **pLN / protein, unavailable** (6551, 6566, 6593, 6623) — the producer's words: "the clustering
  track produced no cluster for this donor; its QuPath export carries no intensity measurement, so
  there is nothing to cluster and nothing to name — never 0 immune cells."
* **pLN / protein, blank** — 6374 and 6591 have a PhenoCycler pancreas section but no pLN image in
  `pln_node_annotations.csv`; 112, 6347, 6533 and 6580 have no PhenoCycler acquisition at all.
* **pLN / protein, withheld** (6579) — see the citation conflict below.
* **pLN / RNA, blank** (6450, 6505, 6534, 6591, 6623) — the lane was not run on a lymph node
  section of these donors.

### A producer's reason that outlived the fact it cites

`data/development/typing/pln/pheno/per_donor_status.csv` (mtime 2026-09-17 16:03) withholds donor
6579 with the reason "withheld by pln_node_annotations.csv (excluded=True)". That file (mtime
2026-09-17 14:16) now records `excluded=False` for 6579 and 1,231,818 cells, following `cohort.py`'s
REINSTATED entry of the same day. The layer carries the state as the producer recorded it and flags
the citation in `presence.citation_conflict`; the card prints the flag next to the note. This
builder does not own the pLN protein lane and does not re-type a donor. **This one is worth a
maintainer's decision: 6579's pLN protein section exists and is no longer excluded.**

## Harmonisation — nothing dropped

Labels go through `core/phenocycler/data/integration/vocab_crosswalk.csv`, which is GENERATED by
`run_vocab` in `phenocycler/integration/vocab.py`. It is read, never hand-edited, and never narrowed
to the lineages that happen to agree. Both Xenium releases already carry their own `common_lineage`
computed against it, so those are used directly; the two PhenoCycler lanes map their fine labels
onto crosswalk `common` names through a table the build validates against the crosswalk (a target
the crosswalk does not declare refuses the build).

Every crosswalk lineage gets a row in every grid cell. Where a lane emitted none, the row carries
`n_cells` NULL — never 0 — and a reason drawn from the crosswalk's own columns: `immune_tissue_scope`
for whether a lineage is in scope for the lymph node at all, `pheno_resolution` / `xenium_resolution`
plus the matching note for whether the modality can resolve it.

### Two measured disagreements with the crosswalk, named rather than resolved

* **Granulocyte**, `pheno_resolution = no` ("not resolvable"), yet the pancreas protein lane emits
  87,713 `Imm-neutrophil` cells and the pLN protein lane 316,305 `Granulocyte` cells. The crosswalk
  row describes the broad label; these come from a per-marker gate and a cluster class.
* **Neural**, `pheno_resolution = no`, yet the subpanel classification emits 50,261 `Nerve` cells
  from a QuPath structure object.
* **Stromal**, `pheno_resolution = yes`, yet the subpanel classification emits no Stromal class at
  all — those cells fall through into Exocrine, which the crosswalk itself marks `default_sink`.

All three are flagged in `composition.resolution_conflict` and shown in the card. The crosswalk is
generated; the fix, if any is wanted, belongs to its producer and not to this lane.

## Cross-modality, measured, not smoothed

Different sections of different blocks. Agreement is informative and so is disagreement; neither
side is tuned toward the other and neither is the target.

T_NK as a share of the immune cells each lane resolved, over the donors measured in both:

| tissue | donors | protein | RNA | median offset (RNA − protein) | rank disagreement |
|---|---|---|---|---|---|
| pancreas | 16 | 8.6 % | 10.3 % | +1.0 pp | 14 of 16 (88 %) |
| pLN | 11 | 27.0 % | 51.8 % | +29.9 pp | 10 of 11 (91 %) |

The verdict already on record for the pLN — protein 31.1 % against RNA 50.8 %, median offset
−22.7 pp, rank disagreement in 91 % of 11 donors — is reproduced here in direction, in donor count
and in rank disagreement (10 of 11 is that same 91 %). The two percentages differ because the
denominator here is the immune cells each lane **resolved to a lineage**, excluding the
`Immune_unresolved` bucket; that moves protein from 31.1 % to 27.0 % and RNA from 50.8 % to 51.8 %,
and widens the median offset from −22.7 pp to −29.9 pp. The larger figure is presented rather than
the smaller one, and the denominator difference is named rather than reconciled.

Per lineage, over the donors measured in both modalities (the card renders this table):

| tissue | lineage | donors | protein % | RNA % | median offset (pp) | how to read it |
|---|---|---|---|---|---|---|
| pancreas | Endocrine | 16 | 2.22 | 2.58 | +0.20 | orthogonal: one side has no 1:1 marker and uses a surrogate panel |
| pancreas | Exocrine | 16 | 80.71 | 63.04 | −17.01 | against a residual: one side assigns this by fallback |
| pancreas | Vascular | 16 | 5.58 | 5.83 | +0.23 | like for like |
| pancreas | Myeloid | 16 | 6.57 | 5.77 | −0.92 | subclass: recovered from per-marker gates, not the broad label |
| pancreas | B_Plasma | 16 | 2.54 | 1.64 | −1.64 | subclass |
| pancreas | T_NK | 16 | 0.89 | 1.03 | −0.00 | subclass |
| pancreas | Granulocyte | 16 | 0.36 | 1.61 | +0.70 | broad to fine |
| pancreas | Neural | 16 | 0.17 | 0.55 | +0.43 | mixed: no / yes |
| pLN | B_Plasma | 11 | 15.75 | 17.72 | +0.11 | subclass |
| pLN | T_NK | 11 | 17.97 | 27.42 | +11.90 | subclass |
| pLN | Myeloid | 11 | 26.50 | 7.29 | −18.20 | subclass |
| pLN | Granulocyte | 4 | 8.10 | 1.15 | −6.40 | broad to fine |

The resolution pair says **what kind of comparison a row is**, not whether to make it. The Endocrine
row is the one the pair marks `yes` / `surrogate`: a protein hormone gate against an RNA surrogate
identity panel with no 1:1 marker. That is orthogonal evidence rather than a check of either side,
and the two routes land 0.20 pp apart on the median donor. The Exocrine row is the mirror image —
the protein side's Exocrine is a `default_sink` that absorbs the stroma the lane cannot name, and it
runs 17 pp above the RNA side, which names Stromal separately at 9.59 % of assessed cells over those same donors.

The pancreas and the pLN are joined on `donor_id` and `disease_status` alone. No distance, adjacency,
neighbourhood or KDE appears anywhere in this layer, and the builder raises `GridError` on any
quantity handed two tissues.

## Acquisition batch

Eta-squared of every per-donor lineage fraction, each against its own permutation null (2000
permutations, 75 quantity × variable rows):

| source | variable | median eta² | median null | p < 0.05 |
|---|---|---|---|---|
| pancreas / protein | `acq_run` (13 levels, 20 donors) | 0.675 | 0.637 | — |
| pancreas / protein | `stain_group` (9 levels) | 0.498 | 0.412 | 0 of 8 |
| pancreas / protein | `disease_status` | 0.093 | 0.081 | — |
| pancreas / RNA | `acq_run` (10 levels, 18 donors) | 0.479 | 0.529 | — |
| pancreas / RNA | `disease_status` | 0.157 | 0.083 | — |
| pLN / protein | `disease_status` | 0.086 | 0.120 | — |
| pLN / RNA | `disease_status` | 0.059 | 0.078 | — |

Pooled over all measured rows: `acq_run` median eta² 0.504 against a null median of 0.538, with 1 of
17 below p = 0.05; `stain_group` 0.498 against 0.412 with 0 of 8; `disease_status` 0.077 against
0.082 with 2 of 25. **At 14–20 donors over 9–13 acquisition levels a large eta-squared is the
expected value under no effect at all**, which is what the null medians say. Nothing is corrected on
any of them, and the Fusion version is never read. Coarse batch is aliased with disease status; only
the finer variables appear here.

The two pLN lanes record no batch variable in this checkout. Those rows are `blank` with the reason
"no `stain_group`/`acq_run` is recorded for pln/protein in this checkout" rather than filled from
another source's.

## Internal coherence

The one quantitative check, and it is a statement about the build rather than about the biology:
the buckets of each grid cell sum to that cell's assessed cells. 72 grid cells checked, 0
mismatched, recorded in `manifest.json["coherence"]` and asserted by the test suite. No result is a
verification target and nothing here reproduces a previous biological finding as a gate.

## The card

`apps/islet/R/card_g1_grid.R`, registered as `g1_grid` on the Cross-Modal tab. It refuses the
artifact when `manifest.source_domain_root` does not match `domain_data_root()` — the two domains
share five donors and no vocabulary, so a matching donor id is not a licence to read. It renders one
panel per cell of the grid for the selected donor: a composition plot when the cell is measured, and
a NAMED NOTE carrying the state and the producer's own reason when it is not. Never a blank panel,
never a zero for a hole. Donors in the layer but outside this domain's cohort roster (112, 6347,
6579, 6580 — 1,504,898 cells) are not selectable and are named in a `withheld` note rather than
dropped in silence.

No donor id, cohort size, status level or typing method is written into any string in the card; the
test suite asserts that by reading the layer and grepping the card's code for every one of them.

## Not claimed

Acinar and ductal. The 60-plex carries no amylase, trypsin, CPA1 or elastase and no CK19, CA19-9,
SOX9, HNF1B, CFTR or MUC1. Exocrine is one label, and the test suite refuses any label matching
those terms.
