# The compartment gradient — islet core, peri-islet expansion, exocrine tissue, pLN

**2026-09-17.** Part 3 of the study question: the immune cells near and in islets against
those in the exocrine tissue and in the pancreatic lymph node.

Everything below was produced by running
`scripts/senior/build_compartment_gradient.py` in this checkout. The artifact is
`data/development/integration/compartments/`. The app card is
`apps/islet/R/card_s3_compartments.R`; the contract tests are
`core/phenocycler/tests/current/test_compartment_gradient.py` (15 passing).

## What the four compartments are, and which three are the same kind of thing

Three are spatial and one is not.

| compartment | definition | tissue |
|---|---|---|
| islet core | `signed_dist_islet_um <= 0` | pancreas |
| peri-islet expansion | `0 < signed_dist_islet_um <= 20` | pancreas |
| exocrine tissue | `signed_dist_islet_um > 20` | pancreas |
| pLN | the whole node section; **no geometry at all** | pancreatic lymph node |

The three pancreas compartments are cut out of one section by one signed distance. The pLN
compartment is **compositional only**. A lymph node is a different block from the pancreas
beside it — for most donors a different slide entirely — so there is no islet in it, and a
signed distance to one is undefined there rather than zero. The pLN column joins to the
pancreas columns on `donor_id` and `donor_status` and on nothing else. No distance,
adjacency, neighbourhood or KDE crosses the boundary, and nothing is co-plotted in space
across the two tissues. `refuse_cross_tissue()` makes that structural: asking the builder
for a pLN distance raises `CrossTissueDistance` rather than returning a number.

One exception is coming and is **not** inferred here. The maintainer is annotating and
exporting the pLN portions of pancreas slides as a distinct dataset, and some of those
lymphoid structures are anatomically attached to the pancreas section — for those, an
islet-to-lymphoid distance *is* meaningful. That classification arrives as data. The
builder names the column it will honour when it lands (`attached_to_pancreas_section`) and
synthesises nothing in the meantime. Until it lands every pancreas/pLN pair is separate.

**The expansion comes from the signed distance, never from the drawn annotation.** Every
expansion is 20 µm, but one donor's drawn expansion measures 23.1 µm against 18.8–19.3 µm
for every other, and all of them are named `_exp20um`, so the annotation cannot be trusted
to mean what it is called. QuPath's `Parent` column is not read here or anywhere in this
module — it mislabels 2.7% of islet-core cells, which fall through to Exocrine.

Where `data/development/typing/islet_membership` exists (13 of 20 donors) the geometric
outline membership is read back as an **audit**, not an override, and the disagreement is
recorded per donor in `manifest.json`. Measured: on 6442, 20,671 core cells by distance
against 20,441 by outline, 1,130 cells disagreeing; on 6539, 9,193 against 9,042, 375
disagreeing. The two are near but not identical, and neither is silently reconciled.

## Denominators

`n_cells` is the eligible, non-excluded cell count of the compartment, repeated on every
subtype row — **do not sum it across rows**. A cell with no nucleus never entered this
tree at all, and what remains but is ineligible is dropped by its recorded name
(`exocrine_inside_islet_border`: 185 cells on 6442, 1,024 on 6539). No cell in
`typing/cells/` carries `qc_analysis_eligible = False`, and no cell in the cohort lacks a
signed distance (`cells_without_signed_distance: 0` for all 20 donors).

Lineage rows sum to **at most** the immune count, not to it: an immune cell whose lineage
was not resolved is still immune. That internal coherence — subsets summing to their
parent, never past it — is asserted by the test suite.

## Labels

PhenoCycler pancreas and PhenoCycler pLN are named by different lanes, so their labels are
harmonised to the common lineages of
`core/phenocycler/data/integration/vocab_crosswalk.csv` (generated; never hand-edited):
`T_NK`, `B_Plasma`, `Myeloid`, `Granulocyte`. Nothing was left unmapped
(`immune_calls_not_in_crosswalk` is empty for all 20 donors). The pancreas fine calls
`Imm-NK-cell` and `Imm-T-cell` both land in `T_NK` because that is the common lineage the
crosswalk records; the fine call is kept in the source column, not discarded.

A Myeloid call from CD68 protein and one from CD68 RNA are the same claim through
different evidence. The Xenium pLN release is reported as its own modality row beside the
PhenoCycler one, never pooled with it and never tuned against it. Cell ids do not cross: a
QuPath UUID and a 10x `cell_id` share no namespace, and nothing here joins on a cell.

## Immune fraction per compartment per donor

Median over donors of the immune fraction of eligible cells, one point per donor,
`ND < Aab+ < T1D`:

| compartment | modality | ND (n) | Aab+ (n) | T1D (n) |
|---|---|---|---|---|
| islet core | phenocycler | 0.0096 (8) | 0.0583 (3) | 0.0199 (9) |
| peri-islet expansion | phenocycler | 0.0885 (8) | 0.1064 (3) | 0.1363 (9) |
| exocrine tissue | phenocycler | 0.1016 (8) | 0.0717 (3) | 0.1460 (9) |
| pLN | phenocycler | 0.6816 (7) | 0.6121 (3) | 0.6977 (4) |
| pLN | xenium | 0.5073 (8) | 0.4086 (1) | 0.5558 (6) |

Three donors carry the whole Aab+ column and one donor carries the Aab+ pLN Xenium cell.
A median of one is that donor, and is written here as such rather than presented as a group.

The pLN rows are a **composition, not a location**: they share no axis with the three
pancreas rows above them and they are stacked in one table only because they are indexed by
the same donors.

## Immune-subtype composition per compartment

Median over donors of each lineage's share of that compartment's immune cells:

| compartment | status | B_Plasma | Granulocyte | Myeloid | T_NK |
|---|---|---|---|---|---|
| core | ND | 0.275 | 0.003 | 0.609 | 0.017 |
| core | Aab+ | 0.051 | 0.003 | 0.566 | 0.093 |
| core | T1D | 0.143 | 0.003 | 0.538 | 0.143 |
| peri | ND | 0.287 | 0.005 | 0.549 | 0.036 |
| peri | Aab+ | 0.110 | 0.004 | 0.621 | 0.143 |
| peri | T1D | 0.203 | 0.007 | 0.614 | 0.148 |
| exocrine | ND | 0.272 | 0.014 | 0.551 | 0.045 |
| exocrine | Aab+ | 0.123 | 0.005 | 0.627 | 0.121 |
| exocrine | T1D | 0.177 | 0.013 | 0.664 | 0.117 |
| pLN (phenocycler) | ND | 0.233 | 0.116 | 0.389 | 0.274 |
| pLN (phenocycler) | Aab+ | 0.223 | 0.126 | 0.398 | 0.307 |
| pLN (phenocycler) | T1D | 0.229 | 0.058 | 0.401 | 0.353 |
| pLN (xenium) | ND | 0.334 | 0.000 | 0.122 | 0.528 |
| pLN (xenium) | Aab+ | 0.314 | 0.024 | 0.193 | 0.468 |
| pLN (xenium) | T1D | 0.300 | 0.000 | 0.131 | 0.575 |

Two modality rows of the same node disagree substantially — the protein lane puts Myeloid
around 0.39–0.40 of pLN immune cells where the RNA lane puts it at 0.12–0.19, and the RNA
lane puts T_NK higher throughout. **That disagreement is a result, not a defect to tune
away.** Neither lane may be adjusted to move toward the other; what the table records is
that two measurements of the same tissue, by different evidence, do not agree about
myeloid share, and the reason is not established here.

The lineage rows sum to less than 1 wherever an immune cell's lineage was not resolved.

## Both stratifications, applied and shown

**1. Islet size.** The endocrine shape signal is real only above roughly ten cells per
islet, so every core and peri quantity is emitted a second time over islets of at least 10
core cells (`islet_size_strata.csv`). Measured per donor
(`islet_size_summary.csv`): median core cells per islet runs from 2 (6539) to 29 (6482);
the number of islets at or above 10 core cells runs from 36 (6380, of 66) to 914 (6505, of
3,038).

The restriction moves the core immune fraction by at most −0.0050 (6505: 0.0639 → 0.0589)
and by −0.0046 on 6539 (0.0234 → 0.0188); every core donor moves down or not at all. In
the peri-islet expansion it moves both ways, from −0.0165 (6505) to +0.0160 (6380). So the
stratification is not cosmetic on the expansion, and the two strata are shown as two claims
rather than merged.

**2. Fragmented islets.** 6539's islets are mostly one-cell classifier fragments: 477
islets, median 2 core cells, 35.0% of them a single cell. The flag is the measured
condition — median core cells below the 10-cell threshold — and it catches seven donors,
not one: 6450 (median 6), 6505 (4), 6523 (9), 6534 (7), 6539 (2), 6566 (7), 6623 (5). Each
is reported and flagged `fragmented_islets` in the card and the artifact. **None is
dropped**: dropping a donor for being inconvenient would be a cohort chosen for its answer.
What the flag forbids is pooling a per-islet number from these donors with the rest
untouched.

## Batch

Acquisition batch is aliased with disease status at the coarse level, so correction is only
ever considered on the finer variables (`stain_group`, `acq_run`), never on the Fusion
version. Eta-squared, each against its own permutation null over 2,000 permutations:

| quantity | grouping | donors | levels | eta² | null median | p |
|---|---|---|---|---|---|---|
| core immune fraction | stain_group | 20 | 9 | 0.630 | 0.386 | 0.157 |
| core immune fraction | acq_run | 20 | 13 | 0.861 | 0.636 | 0.092 |
| core immune fraction | donor_status | 20 | 3 | 0.421 | 0.087 | **0.009** |
| peri immune fraction | stain_group | 20 | 9 | 0.210 | 0.403 | 0.925 |
| peri immune fraction | acq_run | 20 | 13 | 0.545 | 0.634 | 0.703 |
| peri immune fraction | donor_status | 20 | 3 | 0.169 | 0.080 | 0.206 |
| exocrine immune fraction | stain_group | 20 | 9 | 0.202 | 0.409 | 0.944 |
| exocrine immune fraction | acq_run | 20 | 13 | 0.291 | 0.634 | 0.990 |
| exocrine immune fraction | donor_status | 20 | 3 | 0.023 | 0.084 | 0.818 |
| pLN immune fraction (phenocycler) | stain_group | 14 | 9 | 0.764 | 0.625 | 0.178 |
| pLN immune fraction (phenocycler) | acq_run | 14 | 11 | 0.848 | 0.782 | 0.342 |
| pLN immune fraction (xenium) | stain_group | 15 | 8 | 0.469 | 0.490 | 0.554 |
| pLN immune fraction (xenium) | acq_run | 15 | 11 | 0.661 | 0.724 | 0.656 |
| core immune fraction, islets ≥ 10 cells | acq_run | 20 | 13 | 0.867 | 0.619 | 0.082 |
| core immune fraction, islets ≥ 10 cells | donor_status | 20 | 3 | 0.425 | 0.083 | **0.007** |

**The null median is the point.** With 20 donors spread over 13 acquisition runs, an
eta-squared of 0.64 on `acq_run` is what no effect at all produces; the observed 0.861 sits
above it at p = 0.092. Reading 0.861 as "this quantity is a batch axis" without the null
would be reading the design, not the data. By the same token the exocrine fraction's 0.291
on `acq_run` is *below* its null (p = 0.990), which is not evidence of batch cleanliness
either — it is the same small-n geometry seen from the other side.

The only quantity that separates from its own null at p < 0.05 is the **core** immune
fraction on `donor_status` (0.421 against a null median of 0.087, p = 0.009), and it
survives the ≥ 10-cell restriction essentially unchanged (0.425, p = 0.007). Disease status
is itself aliased with coarse batch, so this is reported and **not** corrected on, and it
gates nothing. No quantity was dropped for being confounded and none of these numbers
gates anything on the card.

## The card

`apps/islet/R/card_s3_compartments.R` renders the three pancreas compartments and, in its
own column, the pLN composition. Every path resolves through `ieCore::domain_data_root()`,
and the artifact records the domain it was built from: the two domains share five donors
and no vocabulary, so a matching donor id is not a licence to read, and a mismatch renders
`unavailable(compartment_gradient_domain_mismatch)` rather than the wrong cells.

The pLN column never goes blank. Verified by rendering:

* 6539 — the compositional table, headed by the standing claim that this column joins on
  donor and status alone.
* 6623 — `unavailable(the clustering track produced no cluster for this donor; its QuPath
  export carries no intensity measurement, so there is nothing to cluster and nothing to
  name -- never 0 immune cells)`. The measured reason travels from
  `typing/pln/pheno/per_donor_status.csv` to the panel intact. It is never rendered as zero
  immune cells.
* a donor outside the active domain's cohort — `withheld(...)`, naming the domain boundary.

The artifact refuses itself if a pLN row ever arrives carrying a pancreas compartment, a
stratum or a distance column: that would mean a meaningless number had already been
computed upstream, and displaying it would be worse than refusing it.

## Withdrawn, and not reintroduced here

No number derived from one tissue reaches another tissue's call. No `context`, `fit_levels`
or `reference_level` is applied to a section with no reference population — the pLN has
none, and none is synthesised. There is no absolute `IMMUNE_REQUIRE_MIN` floor and no
raw-count Xenium detection floor anywhere in this lane.

## What this does not say

No acinar or ductal claim is made or could be: the 60-plex carries no amylase, trypsin,
CPA1, elastase, CK19, CA19-9, SOX9, HNF1B, CFTR or MUC1. Exocrine is one label.

Donor exclusions are derived from `phenocycler.cohort.DONOR_EXCLUSIONS` and never restated:
6457 is excluded (technical — poor DAPI prevented accurate channel registration), and
6579's exclusion was lifted 2026-09-17. 6579 has no `typing/cells` partition, so it carries
no pancreas compartment here; its pLN section is `withheld` by
`pln_node_annotations.csv`. Anything pooling 6579 with the rest is pooling a different
disease process, and the larger n is not free.

Five Xenium pLN donors (112, 6347, 6533, 6580 and 6579) are measured in the release and sit
outside this cohort's pancreas partitions; they are named in the manifest as outside the
cohort rather than quietly dropped. They are not missing data.
