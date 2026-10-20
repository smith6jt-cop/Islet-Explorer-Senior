#!/usr/bin/env python3
"""Immune populations in the PhenoCycler pancreatic lymph nodes, named on the CLUSTERS.

THE METHOD, as the maintainer states it:

    "All immune cell calls from pancreas can be applied to pLN"   -- the METHOD transfers.
    "Separate out the immune cells with canonical markers/genes and then use subpanels."
    "Clustering works well in pLN because of the regional cell types."
    "tissue does not mean anything at all."

The last sentence decides the shape of this file, and it is literal.  `cells_pln` has NO
`context` column; `parent_raw` is the single constant string 'Annotation (Tissue)' on all
8,961,168 rows; `qupath_classification` is NULL on every row.  The pancreas's context machinery --
islet / blood vessel / lymphatic / nerve / "tissue" from QuPath pixel-classifier objects -- has
nothing to attach to in a node.  There is no ordinary tissue to draw a reference population from
because there is no partition of the section at all.

**CLUSTERING IS THE SUBSTITUTE FOR CONTEXT.**  In the pancreas a QuPath object hands you the
context for free and a reference level is fitted against it.  In a node the clustering DISCOVERS
the structure, because the regional cell types are distinct.  So:

    1. cluster within donor, per section               -- scripts/pln/cluster_pln.py, consumed here
    2. name each cluster by canonical enrichment RELATIVE TO THE OTHER CLUSTERS ON THAT SECTION
    3. resolve within immune with the five identity_core subpanels
    4. SPLIT THE COMBINED T_NK COMPARTMENT BY NEGATIVE GATING on CD3e -- stage 3 below
    5. cells inherit their cluster's call, and per-cluster heterogeneity is reported so a mixed
       cluster is visible rather than hidden behind its majority label

NEGATIVE EVIDENCE IS EVIDENCE, AND NK IS NAMED.  An earlier generation of this file reported
class_NK_ILC = 0 on all fourteen donors and argued that NK could not be named because no single
channel on this panel is NK-specific.  THAT REASONING IS WITHDRAWN.  No channel has to be: NK is
CD3e-NEGATIVE with CD8 or CD56, T cells are CD3e-POSITIVE, and the panel resolves the two by
COMBINATION.  `lineages.LINEAGE_MARKERS["T_NK"]` carries them as one class because that is the
COARSE label, not because the split is out of reach -- stage 3 makes the split from the
presence/absence pattern inside that class, with the absence carrying as much as the presence.
Nothing is called unavailable that a combination of the channels on hand resolves.

THREE THINGS ARE WITHDRAWN AND ARE NOT IN THIS FILE.

 (a) Any number derived from pancreas cells reaching a pLN call.  An earlier attempt carried each
     donor's pancreas p99 onto their pLN slide; the transferred CD66 level left 1.9 % of the slide
     beneath it and on 12 donors sat below that slide's own 5th percentile.  `pancreas_levels()`,
     `acquisition_factors()`, `transfer_check()` and their `_reference/*.parquet` outputs are gone
     and are asserted gone by the test suite.  Their numbers are a property of the mistake.
 (b) `subpanel.fit_levels` / `reference_level` / `context == 'tissue'`.  The generation that
     followed (a) fitted each marker's level as the p99 of "context == 'tissue'", which on this
     export is the whole slide -- so every marker called 1 % of the slide positive BY
     CONSTRUCTION.  There is no context column and no reference population on a node, and this
     file synthesises neither.  `subpanel_io` is not imported.
 (c) `IMMUNE_REQUIRE_MIN = 0.10` as an absolute floor on pLN values.  It is defined against a
     reference level that does not exist here, so it cannot be applied here.

WHAT "STANDS OUT" MEANS, and why this statistic.

  For cluster *c* and marker *m* on one section, over the analysis-eligible cells of that section
  and that section only:

      AUC[c, m] = P( a random cell in c reads higher on m than a random cell elsewhere on this
                     same section )

  the Mann-Whitney common-language effect size, computed from within-section ranks.  It is:

    * a comparison BETWEEN CLUSTERS ON ONE SECTION -- the "elsewhere" is the rest of that slide,
      so nothing is imported and nothing has to mean the same thing on two slides;
    * threshold-free and unit-free -- it reads ranks, so no intensity value, level, quantile or
      scale factor enters, and 0.5 is "indistinguishable from the rest of this section";
    * robust -- an intensity distribution with a heavy tail cannot drag it, which a mean or a
      p99 both can;
    * exactly testable against an internal null -- under random relabelling of the cells of this
      section the rank-sum has a known mean and variance (tie-corrected), so "stands out beyond
      chance for a cluster of this size on this section" is a calculation, not a judgement.

  ENRICHED means AUC > 0.5 and the tie-corrected rank-sum p below 0.05 Bonferroni-corrected over
  every (cluster x marker) test on that section.  Direction plus an internal significance; no
  intensity threshold anywhere.

WHAT THE RANK STATISTIC CANNOT DO, AND WHAT ARBITRATES INSTEAD.  Corrected 2026-09-18.

  Enrichment is a statement about ONE channel: does this cluster read higher on it than the rest
  of this section.  AUC answers that well and still does.  The immune-vs-structural arbitration
  is a different question -- of TWO CHANNELS of the same cluster, which one is saying more -- and
  AUC cannot answer it, because it is scale-free BY DESIGN: the ordering of pure background is
  still an ordering, so a channel with no target on the section competes on equal terms with a
  channel that marks a population.  `dynamic_range` below has always computed the magnitude AUC
  discards; the generation before this one reported that number beside every exclusion and
  declined to read it, on the ground that reading it would be an absolute threshold.  It is not:
  a comparison between two channels of the same cluster is the comparison this file already makes
  everywhere else, and nothing is compared against a constant.

  MEASURED ON ALL FOURTEEN SECTIONS, with the arbitration on the rank statistic:

    * 1,669,194 of 5,122,269 cells -- 32.6 % of a cohort of lymph nodes -- were excluded as
      non-immune, per donor 20.0 % to 45.8 %.
    * 429,278 of those exclusions (25.7 %) were decided by an EPITHELIAL channel -- EpCAM,
      Pan_Cytokeratin or E_cadherin -- in a tissue with no epithelium, and 155,721 more by
      B3TUBB.
    * On 115 the deciding channel for the largest cluster on the section (115_c12, 35,251 cells)
      was EpCAM at raw AUC 0.936 against CD3e at 0.906.  That cluster's magnitudes are CD8 2.42,
      CD3e 2.05, CD4 1.79 against EpCAM 1.30 -- below its own DAPI, 1.34.  EpCAM's span on that
      whole section is 1.30x, and its median span over the five sections it decided on is 1.60x,
      against 20.3x for Caveolin.  Four more sections repeat it: 6521_c06 (39,676 cells, CD8
      20.5), 6482_c11 (27,890, CD8 6.7), 6505_c05 (27,060, CD3e 6.0) and 6476_c14 (45,864, CD3e
      4.5), each losing to an epithelial or B3TUBB channel under 2.4.
    * The arbitration was in fact on RAW AUC and not on S: the panel offset is a per-cluster
      constant and cancels in a comparison between two markers of one cluster.  Checked on all
      293 clusters -- the markers chosen and the sign of the margin are identical on S and on raw
      AUC -- so the offset never protected the arbitration and was never asked to.

  TWO REFERENCES THAT WERE TESTED AND DO NOT EXPLAIN IT, recorded so they are not retried.
  Recomputing the statistic against the OTHER CLUSTERS of the section, one unit per cluster, the
  way the Xenium pLN lane does, moves 17 of 293 clusters and the cohort's non-immune share from
  32.6 % to 31.1 %; on 115 the non-immune share does not move at all and 115_c12 stays excluded,
  its CD3e AUC RISING from 0.906 to 0.917 under that reference.  Adopting that lane's
  stage-1 rule as well -- at least two of a set enriched, then the larger FRACTION of its own set
  wins -- moves the non-immune share the WRONG WAY, to 47.8 %.  The abundance structure of a node
  is not what breaks this; discarding the magnitude is.

  SO THE ARBITRATION READS THE MAGNITUDE, and the rank statistic keeps the job it can do:
  `enriched` still says which of a set's markers are saying anything about this cluster, and
  `magnitude` -- the cluster's median over the section's median on that channel -- decides which
  of them is saying the most.  It is a ratio inside one channel and one section, so it survives
  any change of that channel's scale; it is dimensionless, so two channels of one cluster can be
  compared; and it is unavailable, never zero, where the section median is not positive.

WHAT THIS FILE DOES NOT LICENSE.  No spatial claim across tissues: pancreas and pLN are different
blocks on one slide and nothing here joins them at all, not even by donor.  No non-immune lineage
beyond what stage 1 must exclude -- the structural markers arbitrate and never become a label.
Never acinar or ductal: the 60-plex carries no marker that could split exocrine, and exocrine has
no meaning in a node anyway.

    python core/multiplex5000/scripts/typing/build_pln_immune_pheno.py --stage standout
    python core/multiplex5000/scripts/typing/build_pln_immune_pheno.py --donors 6442
    python core/multiplex5000/scripts/typing/build_pln_immune_pheno.py
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
import warnings

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from scipy.stats import norm, rankdata

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO))

import lineages as LIN                                                        # noqa: E402
import subpanel as SP                                                         # noqa: E402
from core.multiplex5000.multiplex5000.integration.vocab import (                  # noqa: E402
    PROTEIN_GENE_CANDIDATES,
    SURROGATE_PANEL,
)

PLN_CELLS = REPO / "core/phenocycler/data/integration/cells_pln"
INGEST_REPORT = PLN_CELLS.parent / "cells_pln_ingest_report.json"
IDENTITY_CORE = REPO / "core/panels/data/identity_core"
CLUSTERS = REPO / "data/development/typing/pln/clusters/phenocycler"

OUT = REPO / "data/development/typing/pln/pheno"
REF = OUT / "_reference"

#: The pLN partitioning, spelled out.  Bare `partitioning='hive'` INFERS `donor_id` as int32 on
#: this root -- every pLN donor id happens to be all digits -- which silently turns '0112'-style
#: ids into integers and breaks any join against the string donor ids the rest of the repo uses.
PLN_PARTITIONING = ds.partitioning(
    pa.schema([("donor_id", pa.string()), ("roi", pa.string())]), flavor="hive")

# --------------------------------------------------------------------------------------------
# Stage 1 -- what "immune" means here, marker by marker.
#
# The primaries are `lineages.LINEAGE_MARKERS` for the four immune lineages, MINUS
# `lineages.SUPPORTING_ONLY`.  That module is the maintainer's statement of what each antibody
# means, "made before any result is seen", and it already names the markers that may not carry a
# call alone: CD11b, M2Gal3, FOXP3, CD56, TP63.  Nothing is added to it by hand.
#
# CD45/PTPRC IS NOT ON THIS PANEL.  There is no pan-leukocyte channel, so stage 1 cannot be one
# gate; it is a union of lineage canonicals, arbitrated against the structural markers.
# --------------------------------------------------------------------------------------------
IMMUNE_LINEAGES = ("T_NK", "B_Plasma", "Myeloid", "Granulocyte")
NON_IMMUNE_LINEAGES = ("Exocrine", "Vascular", "Stromal", "Neural")

PRIMARY_IMMUNE: dict[str, tuple[str, ...]] = {
    lin: tuple(m for m in LIN.LINEAGE_MARKERS[lin] if m not in LIN.SUPPORTING_ONLY)
    for lin in IMMUNE_LINEAGES
}
SUPPORTING_IMMUNE: tuple[str, ...] = tuple(sorted(
    {m for lin in IMMUNE_LINEAGES for m in LIN.LINEAGE_MARKERS[lin]} & LIN.SUPPORTING_ONLY))

#: Structural markers, from the same module's non-immune lineages, with two removed by name.
#:
#: Vimentin  -- `lineages` puts it in Stromal, which is a statement about the PANCREAS.  Every
#:              lymphocyte carries vimentin, so using it to exclude would exclude the node.
#: TP63      -- already `SUPPORTING_ONLY`: it may not carry a call alone, so it may not exclude
#:              a cluster on its own either.
STRUCTURAL_DROPPED = {"Vimentin": "carried by leukocytes; excluding on it would exclude the node",
                      "TP63": "lineages.SUPPORTING_ONLY -- may not carry a call alone"}
STRUCTURAL: tuple[str, ...] = tuple(sorted(
    {m for lin in NON_IMMUNE_LINEAGES for m in LIN.LINEAGE_MARKERS[lin]} - set(STRUCTURAL_DROPPED)))

#: Channels on the pLN export that carry NO call, with the reason -- which is about the TAXONOMY,
#: never about a threshold.  `lineages.LINEAGE_MARKERS` is the authority on which markers may name
#: a lineage; a channel it does not list cannot name one, however suggestive the antibody.  They
#: are still measured and written to the per-cluster standout table, where a reader uses them to
#: recognise the regional structures -- CD21 for the follicular dendritic network, Ki67/PCNA/Bcl_2
#: for the germinal centre, TCF_1/TOX/PD_1 for T-zone state -- but they decide nothing.
NOT_A_NAMING_MARKER = ("CD14", "CD21", "CD141", "GP100", "IFNG", "HLA_DR", "HLA_A", "CD57",
                       "CD39", "CD44", "CD107a", "IDO1", "iNOS", "SOX2", "PDL1", "Bcl_2",
                       "Beta_actin", "PCNA", "Ki67", "PD_1", "TCF_1", "TOX", "ICOS", "LAG3",
                       "VISTA", "DAPI")

#: Columns of the export that are not channels.  The panel is everything else, so the offset
#: below is estimated over the WHOLE panel and no hand-picked set of markers can steer it.
NOT_A_CHANNEL = ("X_centroid", "Y_centroid", "nucleus_area_um2", "cell_area", "cell_solidity",
                 "nucleus_cell_area_ratio")

#: A cluster's marker is ENRICHED when three things hold, and each has a null of its own rather
#: than a tuned value:
#:
#:   AUC > AUC_NULL     it reads higher than the rest of its own section.  0.5 is not a threshold;
#:                      it is the value AUC takes when a cluster is indistinguishable from the
#:                      rest of the slide.
#:   p < ALPHA / tests  beyond what random relabelling of that section would give, Bonferroni-
#:                      corrected over every cluster x marker test on that section.
#:   S > 0              above the cluster's OWN panel-wide level, where
#:                      S[c, m] = AUC[c, m] - median over every channel of AUC[c, .].
#:
#: THE THIRD CONJUNCT IS THERE BECAUSE THE OFFSET IS MEASURED.  On donor 115 the per-cluster
#: median AUC across the panel runs 0.107 to 0.794 -- a cluster-wide offset of 0.687 AUC -- and a
#: cluster's EpCAM AUC correlates 0.973 with that median (CD3e, 0.810).  Without S a marker counts
#: as enriched merely because its whole cluster reads high on everything, and with 200k cells a
#: rank-sum p is significant almost for free, so S is what makes "enriched" mean anything.
#:
#: WHAT S DOES NOT DO, stated because it would be easy to claim otherwise: it does not change the
#: immune-vs-structural arbitration or the subpanel argmax.  Both are comparisons between two
#: markers of the SAME cluster, and a per-cluster constant cancels in a comparison -- the calls and
#: the separation margins are identical with and without it, measured.  S changes which markers
#: are admitted as evidence, nothing else.
AUC_NULL = 0.5
ALPHA = 0.05

#: Markers whose enrichment is reported but which may not name a cluster on their own, carried
#: from `lineages.SUPPORTING_ONLY`.  The one exception is licensed by the lane's own rule table,
#: never by hand -- see NK_LICENCE.
SUPPORTING_ONLY_SUBPANEL = {"NK_ILC": "Imm-NK-cell"}

SUBPANEL_FILES = {
    "T_cell": "03_immune_T_cell.csv",
    "B_Plasma": "04_immune_B_plasma.csv",
    "NK_ILC": "05_immune_NK_ILC.csv",
    "Myeloid": "06_immune_myeloid.csv",
    "Granulocyte": "07_immune_granulocyte.csv",
}

#: The section-internal rank definitions used ONLY for the coherence check.  A cell is "positive"
#: for a marker at q when its value is in the top q of that section for that marker.  This decides
#: NOTHING -- no call reads it -- and it is reported across a grid precisely so that the coherence
#: statement does not rest on one choice.  CD4 % + CD8 % among the same T cells cannot pass 100 %
#: at ANY q unless the two are calling the same cells twice.
COHERENCE_Q = (0.25, 0.10, 0.05, 0.01)


# --------------------------------------------------------------------------------------------
# Reading.  ONE accessor for the cells, one for the clustering.
# --------------------------------------------------------------------------------------------
def pln_source(donor: str) -> pathlib.Path:
    """The one place a pLN donor's cells are read from.

    `subpanel_io.cell_source()` is the equivalent for the PANCREAS and is deliberately not called
    here -- it resolves into the pancreas revised/run trees and this lane reads no pancreas cell
    at all.  The discipline it exists for is kept: this function is the only path to pLN values in
    this file, so no raw glob of a retired tree can creep back in.
    """
    part = PLN_CELLS / f"donor_id={donor}/roi=pln_1"
    if not part.is_dir():
        raise FileNotFoundError(f"{donor}: no pLN partition at {part}")
    return part


def pln_columns() -> list[str]:
    """Marker-bearing columns of the pLN export: whole-cell means, not the nucleus copies."""
    donor = declared_written_donors()[0]
    file = next(iter(pln_source(donor).glob("*.parquet")))
    return [c for c in pq.ParquetFile(file).schema_arrow.names if not c.startswith("Nucleus__")]


def panel_channels() -> list[str]:
    """Every marker channel the pLN export carries, derived from its schema and not listed here.

    A channel is a floating-point column that is not geometry.  The whole panel is used to
    estimate each cluster's offset, so the offset cannot be steered by which markers this lane
    happens to care about.
    """
    donor = declared_written_donors()[0]
    file = next(iter(pln_source(donor).glob("*.parquet")))
    schema = pq.ParquetFile(file).schema_arrow
    return [n for n, t in zip(schema.names, schema.types)
            if not n.startswith("Nucleus__") and pa.types.is_floating(t)
            and n not in NOT_A_CHANNEL]


def read_pln(donor: str, columns: list[str]) -> pd.DataFrame:
    """One donor's pLN cells, read through the explicit partitioning schema."""
    pln_source(donor)
    dataset = ds.dataset(PLN_CELLS, format="parquet", partitioning=PLN_PARTITIONING)
    have = set(dataset.schema.names)
    wanted = [c for c in columns if c in have]
    table = dataset.to_table(columns=wanted, filter=ds.field("donor_id") == str(donor))
    frame = table.to_pandas()
    for missing in [c for c in columns if c not in have]:
        # A marker this export does not carry is UNAVAILABLE, never zero.
        frame[missing] = np.nan
    return frame


def declared_written_donors() -> list[str]:
    return sorted(p.name.split("=")[1] for p in PLN_CELLS.glob("donor_id=*"))


def declared_manifest_donors() -> list[str]:
    report = json.loads(INGEST_REPORT.read_text())
    return [str(d) for d in report["donor_mapping"]["donors_in_manifest"]]


def withheld_reason(donor: str) -> str:
    report = json.loads(INGEST_REPORT.read_text())
    return report["donor_mapping"]["donors_refused"].get(donor, "refused by the pLN ingest")


def clustered_donors() -> list[str]:
    """Donors the clustering track produced an assignment for."""
    if not (CLUSTERS / "assignments").is_dir():
        return []
    return sorted(p.name.split("=")[1] for p in (CLUSTERS / "assignments").glob("donor_id=*"))


def cluster_assignment(donor: str) -> pd.DataFrame:
    """object_id -> cluster_id for one section, from the clustering track.

    This lane CONSUMES the clustering and never re-runs it.  The join key is the QuPath UUID; no
    centroid matching, no row order, no filename parsing.
    """
    path = CLUSTERS / f"assignments/donor_id={donor}/clusters.parquet"
    frame = pq.read_table(path, columns=["object_id", "cluster_id"]).to_pandas()
    frame["object_id"] = frame["object_id"].astype(str)
    frame["cluster_id"] = frame["cluster_id"].astype(str)
    if not frame["object_id"].is_unique:
        raise ValueError(f"{donor}: the clustering assignment repeats an object_id")
    return frame


def cluster_profile(donor: str) -> pd.DataFrame:
    """The clustering track's own per-cluster profile, read to CHECK the join, not to call on."""
    path = CLUSTERS / f"profiles/{donor}.parquet"
    if not path.exists():
        return pd.DataFrame(columns=["cluster_id", "n_cells"])
    return pq.read_table(path, columns=["cluster_id", "n_cells", "frac_of_donor"]).to_pandas()


# --------------------------------------------------------------------------------------------
# The statistic: a cluster against the other cells of its own section.
# --------------------------------------------------------------------------------------------
def standout(values: pd.DataFrame, codes: np.ndarray, n_clusters: int, markers: list[str],
             offset_over: list[str] | None = None) -> dict[str, np.ndarray]:
    """AUC, tie-corrected rank-sum z and p for every (cluster, marker) on ONE section.

    `AUC[c, m]` is the probability that a random cell of cluster c reads higher on m than a random
    cell elsewhere on this same section, computed from the section's own ranks.  Its null is not
    assumed: under random relabelling of this section's cells the rank sum of a cluster of size n
    has mean n(N+1)/2 and a tie-corrected variance, so z says whether the separation exceeds what
    a cluster of that size gets by chance HERE.  Nothing outside the section is read.

    `S[c, m] = AUC[c, m] - median over every channel of AUC[c, .]` is the same statistic with the
    cluster's panel-wide offset taken out -- see the ENRICHED comment for the measurement that
    motivates it.  `offset_over` names the channels the median runs over; it is the whole panel.

    A marker with no finite measurement on a cell contributes nothing for that cell: it is ranked
    only among the cells that carry it, and `n_finite` records how many those were.  A missing
    measurement is unavailable, never zero and never the bottom rank.
    """
    n = len(values)
    auc = np.full((n_clusters, len(markers)), np.nan)
    z = np.full((n_clusters, len(markers)), np.nan)
    p = np.full((n_clusters, len(markers)), np.nan)
    n_fin = np.zeros((n_clusters, len(markers)), dtype=np.int64)
    unit: dict[str, np.ndarray] = {}
    unit_mag: dict[str, np.ndarray] = {}

    for j, marker in enumerate(markers):
        v = values[marker].to_numpy(dtype=float)
        finite = np.isfinite(v)
        total = int(finite.sum())
        u = np.full(n, np.nan, dtype=np.float32)
        if total < 2:
            unit[marker] = u
            unit_mag[marker] = np.full(n, np.nan, dtype=np.float32)
            continue
        ranks = rankdata(v[finite])                      # average ranks, 1..total
        u[finite] = (ranks / total).astype(np.float32)   # the cell's own place on this section
        unit[marker] = u
        # The cell's own magnitude on this channel, on the same footing as the cluster's: its
        # value over the section's median for that channel.  A cross-channel comparison inside
        # ONE cell reads this, never `unit[marker]`, for the reason stage 1 does -- a rank is
        # scale-free, so a cell's place in a channel that is all background is still a place, and
        # comparing it against its place in a channel that marks something compares nothing.
        med = float(np.median(v[finite]))
        um = np.full(n, np.nan, dtype=np.float32)
        if med > 0:
            um[finite] = (v[finite] / med).astype(np.float32)
        unit_mag[marker] = um

        c = codes[finite]
        rank_sum = np.bincount(c, weights=ranks, minlength=n_clusters)
        size = np.bincount(c, minlength=n_clusters).astype(float)
        rest = total - size
        # Mann-Whitney U of the cluster against the rest of the section.
        stat = rank_sum - size * (size + 1.0) / 2.0
        with np.errstate(invalid="ignore", divide="ignore"):
            auc[:, j] = np.where((size > 0) & (rest > 0), stat / (size * rest), np.nan)
        _, tie_counts = np.unique(v[finite], return_counts=True)
        tie_term = float(((tie_counts ** 3) - tie_counts).sum())
        var = size * rest / 12.0 * ((total + 1.0) - tie_term / (total * (total - 1.0)))
        with np.errstate(invalid="ignore", divide="ignore"):
            zz = (stat - size * rest / 2.0) / np.sqrt(var)
        zz = np.where((size > 0) & (rest > 0) & (var > 0), zz, np.nan)
        z[:, j] = zz
        p[:, j] = 2.0 * norm.sf(np.abs(zz))
        n_fin[:, j] = size.astype(np.int64)

    over = [m for m in (offset_over or markers) if m in markers]
    cols = [markers.index(m) for m in over]
    with warnings.catch_warnings():
        # A cluster every one of whose channels is unmeasured has no offset.  That is the
        # unavailable case and it is carried as NaN, not as the null.
        warnings.simplefilter("ignore", RuntimeWarning)
        offset = np.nanmedian(auc[:, cols], axis=1) if cols else np.full(n_clusters, AUC_NULL)
    rel = auc - offset[:, None]

    # THE PER-CELL ANALOGUE OF THE CLUSTER'S PANEL OFFSET.  `offset` above takes out a CLUSTER's
    # panel-wide level because without it a marker counts as enriched merely because its whole
    # cluster reads high on everything.  The same is true one level down and it is what blocks a
    # split: inside a cluster, a dim cell reads low on CD3e, CD8 and CD56 alike, so comparing the
    # raw ranks of a CD3e-low side against a CD3e-high side compares brightness, not identity --
    # measured on all 14 sections, the CD3e-low side reads LOWER on CD8 in every one of the 32
    # compartment clusters (AUC 0.193 to 0.513, never above 0.52).
    #
    # `cell_offset[i]` is cell i's own median place on this section over every channel, and
    # `unit_rank[m] - cell_offset` is that cell's place on m with its brightness taken out.  It is
    # still a rank difference: no intensity, level or quantile enters, and a constant per cell
    # cancels in any comparison WITHIN one cell -- which is why the subpanel argmax is unchanged
    # by it, and why only comparisons BETWEEN cells (the two sides of a split) move.
    cell_offset = np.zeros(n, dtype=np.float32)
    if over:
        step = 200_000
        for a in range(0, n, step):
            block = np.vstack([unit[m][a:a + step] for m in over])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                cell_offset[a:a + step] = np.nanmedian(block, axis=0)

    tested = int(np.isfinite(p).sum())
    alpha = ALPHA / max(tested, 1)
    enriched = (np.isfinite(auc) & (auc > AUC_NULL) & np.isfinite(p) & (p < alpha)
                & np.isfinite(rel) & (rel > 0.0))
    return {"auc": auc, "rel": rel, "offset": offset, "offset_over": over, "z": z, "p": p,
            "n_finite": n_fin, "enriched": enriched, "alpha": alpha, "n_tests": tested,
            "unit_rank": unit, "unit_magnitude": unit_mag, "cell_offset": cell_offset,
            "markers": markers}


# --------------------------------------------------------------------------------------------
# Stage 2 -- subpanel membership, from core/panels and nowhere else.
# --------------------------------------------------------------------------------------------
def magnitude(values: pd.DataFrame, codes: np.ndarray, n_clusters: int,
              markers: list[str]) -> np.ndarray:
    """Per (cluster, marker) on ONE section: the cluster's median over the section's median.

    THIS IS THE MAGNITUDE AUC DISCARDS, AND STAGE 1 ARBITRATES ON IT.  See the module docstring
    for the measurement that put it there.  The statistic is:

      * a ratio inside ONE channel and ONE section, so multiplying that channel by any constant
        leaves it unchanged -- two sections acquired at different exposures give the same number,
        which is the whole of what "no level that has to mean the same thing on two slides" asks;
      * dimensionless, so two channels of the same cluster can be compared -- which is the only
        comparison stage 1 makes;
      * unavailable, never zero, wherever the section median is not positive: a channel whose
        section median is 0 or absent (``Ker8_18``, ``Keratin_5`` on this export) has no
        magnitude here and therefore cannot arbitrate at all.

    It is NOT invariant to an additive background, and that is the point rather than a defect: a
    channel whose whole section is background lands every cluster near 1 and stops deciding, which
    is exactly what a rank statistic cannot see.
    """
    out = np.full((n_clusters, len(markers)), np.nan)
    for j, marker in enumerate(markers):
        v = values[marker].to_numpy(dtype=float)
        finite = np.isfinite(v)
        if finite.sum() < 2:
            continue
        section = float(np.median(v[finite]))
        if not (section > 0):
            continue
        per = pd.Series(v[finite]).groupby(codes[finite]).median()
        for code, med in per.items():
            if 0 <= int(code) < n_clusters:
                out[int(code), j] = float(med) / section
    return out


def dynamic_range(values: pd.DataFrame, codes: np.ndarray, n_clusters: int,
                  markers: list[str]) -> dict[str, float]:
    """Per channel on ONE section: the brightest cluster's median over the section's median.

    The per-channel maximum of :func:`magnitude`, kept because the structural-veto accounting
    reports the deciding channel's whole span on the section beside every exclusion it caused --
    a reader needs the channel's span as well as the one cluster's lift to judge it.  ONE
    definition backs both: this function computes nothing of its own.
    """
    mag = magnitude(values, codes, n_clusters, markers)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        top = np.nanmax(mag, axis=0)
    return {m: float(top[j]) for j, m in enumerate(markers)}


def subpanel_membership(panel_markers: list[str]) -> dict[str, dict[str, str]]:
    """subpanel -> {protein: subtype_hint}, for the proteins this panel actually carries.

    Membership is `core/panels/data/identity_core/*.csv`, reached through the verified
    protein<->gene crosswalk in `phenocycler.integration.vocab`.  A protein whose gene is not in a
    subpanel's CSV is not in that subpanel, however obvious it looks: CD209 is a myeloid canonical
    and is used in stage 1, but `06_immune_myeloid.csv` does not list CD209, so it carries no
    stage-2 evidence; CD66's genes are CEACAM1/5/6 and `07_immune_granulocyte.csv` lists CEACAM3,
    so the granulocyte subpanel's only naming marker here is MPO.  core/panels is the one
    authority for the taxonomy and neither case is patched by hand.
    """
    identity = set(LIN.identity_markers())
    out: dict[str, dict[str, str]] = {}
    for name, filename in SUBPANEL_FILES.items():
        table = pd.read_csv(IDENTITY_CORE / filename)
        hints = dict(zip(table["gene"], table["subtype_hint"].fillna("")))
        members = {}
        for protein in panel_markers:
            # A protein with no gene on the Xenium panel reaches its subpanel through
            # `vocab.SURROGATE_PANEL`, which names which panel it routes to.  That is how Iba1 and
            # M2Gal3 -- myeloid antibodies whose genes AIF1 and LGALS3 are not in
            # 06_immune_myeloid.csv -- are placed without inventing membership.
            if SURROGATE_PANEL.get(protein) == filename[:-4]:
                members[protein] = ""
                continue
            for gene in PROTEIN_GENE_CANDIDATES.get(protein, []):
                if gene in hints:
                    members[protein] = hints[gene]
                    break
        # State markers may not name a cell.  `lineages` holds back proliferation, exhaustion and
        # activation channels by its own design -- "they say what a cell is doing, not what it
        # is" -- so PD_1, TCF_1, ICOS, LAG3 and VISTA carry no stage-2 evidence even though the
        # identity-core CSVs list their genes.
        out[name] = {p: h for p, h in members.items()
                     if p in identity or p in LIN.SUPPORTING_ONLY}
    return out


def naming_markers(members: dict[str, dict[str, str]]) -> dict[str, list[str]]:
    """Per subpanel, the markers that may name it -- its members minus `SUPPORTING_ONLY`."""
    return {name: sorted(m for m in hints if m not in LIN.SUPPORTING_ONLY)
            for name, hints in members.items()}


#: Which phenotype of `subpanel.immune_rules()` each identity_core subpanel is the coarse form of.
#: The mapping is by NAME, one line, so that a reader can check it against the table in one look;
#: nothing about which markers a phenotype requires is written here.  `NK_ILC` is deliberately
#: absent: it has no naming marker on this panel and reaches stage 2 only through `nk_licence()`,
#: which reads the same table's whole rule rather than its `all_of` clause alone.
SUBPANEL_RULE = {
    "T_cell": "Imm-T-cell",
    "B_Plasma": "Imm-B-cell",
    "Myeloid": "Imm-macrophage",
    "Granulocyte": "Imm-neutrophil",
}


def subpanel_requirements() -> dict[str, tuple[str, ...]]:
    """Per subpanel, the markers its own rule-table phenotype REQUIRES -- the `all_of` clause.

    Read from `subpanel.immune_rules()` and never written by hand.  On this panel the table gives
    `Imm-T-cell` `['CD3e']`, `Imm-B-cell` `['CD20', 'CD79a']`, `Imm-neutrophil` `['CD66']` and
    `Imm-macrophage` `[]` -- macrophage requires nothing because the table carries its evidence in
    an `any_of` group instead, and an empty requirement is an empty requirement, not a licence to
    invent one.

    WHY THE `all_of` CLAUSE AND NOT THE WHOLE RULE.  The `none_of` clause of every immune
    phenotype in that table is the four epithelial channels, and it is evaluated there against a
    LEVEL (`subpanel.IMMUNE_REQUIRE_MIN` on scaled values) -- "this marker is actually present".
    This lane has no level, and "enriched relative to this cluster's own median channel" is not
    the same statement: a cluster can sit above its own median on a channel that marks nothing
    anywhere on the section.  Carrying the negatives across on enrichment is exactly the
    substitution that put 429,278 cells of a lymph node behind an epithelial veto, so the
    negatives are NOT carried; what is carried is the requirement, which is a statement about
    evidence that must be present and which enrichment can express.
    """
    table = SP.immune_rules()
    return {name: tuple(table[rule][0]) for name, rule in SUBPANEL_RULE.items() if rule in table}


def nk_licence() -> tuple[list[str], list[list[str]], list[str]]:
    """NK_ILC's licence, taken from `subpanel.immune_rules()` and not written by hand.

    CD56 is `SUPPORTING_ONLY` -- shared with the Neural lineage -- and it is the NK subpanel's only
    member on this panel, so NK_ILC has no naming marker here.  The lane's own rule table says
    what CD56 needs before it may speak: `Imm-NK-cell` is CD56 with (CD57 or Granzyme B), CD3e
    negative and the epithelial channels negative.  Both partners are on the pLN panel, so no fact
    about the pancreas is needed to run it.  It is evaluated on ENRICHMENT here, never on a level.

    WHAT THIS MEANS FOR THE COUNTS, stated because a zero must never be read as a measurement.
    CD56 is `SUPPORTING_ONLY`, so it is not among the stage-1 primaries, so a cluster whose ONLY
    standing-out canonical is CD56 does not enter the immune compartment at all -- `lineages` says
    so directly: "A lineage whose evidence rests only on these is Ambiguous, not named."  NK_ILC
    can therefore only be reached by a cluster admitted on some other lineage's canonical that
    then satisfies this licence and reads higher on CD56 than on that canonical.  A count of zero
    NK_ILC clusters is a property of THIS PANEL -- one NK member, and it shared with the Neural
    lineage -- and not a statement that the node holds no NK cells.
    """
    all_of, any_of, none_of = SP.immune_rules()[SUPPORTING_ONLY_SUBPANEL["NK_ILC"]]
    # The second any_of group is the lane's broad "is it an immune cell at all" clause, which
    # stage 1 has already answered for a cluster that reached stage 2; the discriminating group is
    # the one that does not contain CD56's competitors.
    groups = [list(g) for g in any_of if "CD56" not in g and len(g) <= 3]
    return list(all_of), groups, list(none_of)


# --------------------------------------------------------------------------------------------
# Stage 1 and stage 2, ON CLUSTERS.
# --------------------------------------------------------------------------------------------
def _best(score_row: np.ndarray, enr_row: np.ndarray, index: dict[str, int],
          markers: list[str]) -> tuple[float, str, bool]:
    """The strongest of a marker set for one cluster on `score_row`, preferring ENRICHED.

    Returns (score, marker, enriched).  When no member of the set is enriched, the strongest
    member is still reported -- so a reader sees how close the cluster came -- with `enriched`
    False, and it cannot carry a call.

    TWO STATISTICS, EACH DOING THE JOB IT CAN DO.  `enr_row` is the rank statistic's verdict: it
    says which of the set's markers are saying anything about this cluster at all, direction and
    internal null and the cluster's own panel level, exactly as before.  `score_row` is the
    MAGNITUDE, and it decides which of those markers is saying the most.  A rank cannot make that
    second comparison -- it is scale-free by construction, so it cannot tell a channel that marks
    a population from a channel that is background, and the measurement in the module docstring
    is what happens when it is asked to.
    """
    present = [m for m in markers if m in index and np.isfinite(score_row[index[m]])]
    if not present:
        return float("nan"), "", False
    enriched = [m for m in present if enr_row[index[m]]]
    pool = enriched or present
    top = max(pool, key=lambda m: score_row[index[m]])
    return float(score_row[index[top]]), top, bool(enriched)


def call_clusters(donor: str, cluster_ids: list[str], sizes: np.ndarray,
                  st: dict, members: dict) -> pd.DataFrame:
    """The two stages, per cluster, on one section.

    STAGE 1.  A cluster is immune when a canonical immune marker is ENRICHED against the other
    cells of this same section and carries MORE MAGNITUDE on this cluster than any enriched
    structural marker does.  It is non-immune when a structural marker is enriched and the immune
    evidence does not beat it on that magnitude.  A cluster where nothing is enriched is
    UNCALLED -- a real category.  In the pancreas such a cell falls through to Exocrine because
    exocrine is the background there; a node has no such default, and inventing one would
    manufacture a majority class out of nothing.

    WHY THE MAGNITUDE AND NOT THE RANK.  The arbitration is a comparison between TWO CHANNELS of
    the same cluster, and a within-section rank cannot make it: see the module docstring for the
    measurement.  The rank statistic still decides which markers are admitted as evidence
    (`enriched`); the magnitude decides which of the admitted markers wins.

    STAGE 2.  Within the immune clusters only, the five identity_core subpanels compete on the
    strongest of THEIR OWN naming markers, on the same magnitude.  The winner is then held to TWO
    requirements, and fails to `Uncalled` rather than to the runner-up:

      * it must have an ENRICHED canonical of its own lineage.  That is the failure the Xenium
        draft made -- 55.4 % of its Granulocyte calls detected none of the granulocyte
        canonicals -- and the requirement is written so it cannot recur here.
      * every marker in its `all_of` clause in `subpanel.immune_rules()` must be enriched.  That
        table is the maintainer's statement of what each phenotype REQUIRES, it is the same table
        `nk_licence()` already reads, and it is taken whole rather than marker by marker: it is
        why no cluster is named T on CD4 or CD8 without CD3e -- `Imm-T-cell` is `all_of ['CD3e']`
        and every T subtype in that table repeats it.  A cluster whose winner fails its own
        clause is Uncalled WITH THE MISSING MARKERS NAMED.
    """
    index = {m: j for j, m in enumerate(st["markers"])}
    auc, rel, enr = st["auc"], st["rel"], st["enriched"]
    mag = st["magnitude"]
    primaries = sorted({m for lin in IMMUNE_LINEAGES for m in PRIMARY_IMMUNE[lin]})
    naming = naming_markers(members)
    nk_all, nk_any, nk_none = nk_licence()
    required = subpanel_requirements()

    rows = []
    for i, cid in enumerate(cluster_ids):
        imm_mag, imm_marker, imm_enr = _best(mag[i], enr[i], index, primaries)
        str_mag, str_marker, str_enr = _best(mag[i], enr[i], index, list(STRUCTURAL))

        if not np.isfinite(imm_mag) and not np.isfinite(str_mag):
            is_immune, stage1 = None, "no_channel_carries_a_magnitude_on_this_section"
        elif imm_enr and (not str_enr or not np.isfinite(str_mag) or imm_mag > str_mag):
            is_immune, stage1 = True, "immune_canonical_carries_more_magnitude_than_structural"
        elif str_enr and (not imm_enr or not np.isfinite(imm_mag) or str_mag >= imm_mag):
            is_immune, stage1 = False, "structural_carries_more_magnitude_than_immune"
        elif imm_enr and str_enr:
            is_immune, stage1 = None, "immune_and_structural_evidence_tie"
        else:
            is_immune, stage1 = None, "no_lineage_marker_stands_out_on_this_section"

        rel_at = lambda m: (float(rel[i][index[m]]) if m in index else np.nan)   # noqa: E731
        row = {"donor_id": donor, "cluster_id": cid, "n_cells": int(sizes[i]),
               "panel_offset_auc": float(st["offset"][i]),
               "immune_magnitude": imm_mag, "top_immune_marker": imm_marker,
               "immune_standout": rel_at(imm_marker),
               "top_immune_marker_raw_auc": (float(auc[i][index[imm_marker]])
                                             if imm_marker in index else np.nan),
               "top_structural_marker_raw_auc": (float(auc[i][index[str_marker]])
                                                 if str_marker in index else np.nan),
               "immune_marker_enriched": imm_enr,
               "structural_magnitude": str_mag, "top_structural_marker": str_marker,
               "structural_standout": rel_at(str_marker),
               "structural_marker_enriched": str_enr,
               "magnitude_margin": imm_mag - str_mag,
               "separation_margin": rel_at(imm_marker) - rel_at(str_marker),
               "is_immune": is_immune, "stage1_rule": stage1,
               "immune_class": "", "stage2_reason": "", "immune_subtype": ""}
        for name in SUBPANEL_FILES:
            row[f"standout_{name}"] = np.nan
            row[f"magnitude_{name}"] = np.nan
            row[f"enriched_{name}"] = False

        # --- stage 2 --------------------------------------------------------------------------
        scores: dict[str, float] = {}
        for name, markers in naming.items():
            if name in SUPPORTING_ONLY_SUBPANEL:
                continue
            a, m, e = _best(mag[i], enr[i], index, markers)
            row[f"magnitude_{name}"], row[f"enriched_{name}"] = a, e
            row[f"standout_{name}"] = rel_at(m)
            row[f"marker_{name}"] = m
            if np.isfinite(a):
                scores[name] = a
        # NK_ILC has no naming marker on this panel; its licence is the lane's rule table.
        licensed = (all(enr[i][index[m]] for m in nk_all if m in index)
                    and all(any(enr[i][index[m]] for m in g if m in index) for g in nk_any)
                    and not any(enr[i][index[m]] for m in nk_none if m in index))
        row["nk_licensed_by_rule_table"] = bool(licensed)
        if "CD56" in index and np.isfinite(mag[i][index["CD56"]]):
            row["magnitude_NK_ILC"] = float(mag[i][index["CD56"]])
            row["standout_NK_ILC"] = rel_at("CD56")
            row["enriched_NK_ILC"] = bool(licensed)
            row["marker_NK_ILC"] = "CD56"
            if licensed:
                scores["NK_ILC"] = float(mag[i][index["CD56"]])

        if is_immune is True:
            if not scores:
                row["immune_class"] = "Uncalled"
                row["stage2_reason"] = "no subpanel has a measurable marker on this section"
            else:
                # The winner is the strongest subpanel WHATEVER its enrichment; it is then held to
                # the requirements.  Taking the argmax over enriched subpanels only would quietly
                # reassign the cluster to the runner-up, which is exactly what must not happen.
                won = max(scores, key=lambda k: scores[k])
                missing = [m for m in required.get(won, ())
                           if m not in index or not enr[i][index[m]]]
                if not row[f"enriched_{won}"]:
                    row["immune_class"] = "Uncalled"
                    row["stage2_reason"] = (
                        f"winning subpanel {won} has no enriched canonical of its own lineage "
                        f"-- not reassigned to a smaller-scoring subpanel")
                elif missing:
                    row["immune_class"] = "Uncalled"
                    row["stage2_reason"] = (
                        f"winning subpanel {won} does not satisfy its own rule-table requirement "
                        f"{'+'.join(required[won])} -- {'+'.join(missing)} is not enriched on "
                        f"this cluster; not reassigned to a smaller-scoring subpanel")
                else:
                    row["immune_class"] = won
                    row["stage2_reason"] = "subpanel_argmax_with_an_enriched_canonical"
                    row["immune_subtype"] = members[won].get(row.get(f"marker_{won}", ""), "")
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------
# Per-cluster heterogeneity, measured at the cell level and threshold-free.
# --------------------------------------------------------------------------------------------
def structural_veto_rows(donor: str, calls: pd.DataFrame,
                         dyn: dict[str, float]) -> list[dict]:
    """Every exclusion a structural channel caused on one section, with the numbers behind it.

    Both statistics are on the row, and so is the deciding channel's whole span on that section:
    the arbitration reads the MAGNITUDE, the rank is what admitted the channel as evidence, and a
    reader needs all three to judge an exclusion rather than take it on trust.  This is also the
    table on which the defect corrected on 2026-09-18 was visible for a generation -- EpCAM
    deciding 35,251 cells on donor 115 at a span of 1.30x was printed here and read as a
    curiosity.
    """
    rows = []
    vetoed = calls.loc[calls["is_immune"] == False]                              # noqa: E712
    for marker, group in vetoed.groupby("top_structural_marker"):
        rows.append({
            "donor_id": donor, "structural_marker": str(marker),
            "n_clusters_excluded": len(group),
            "n_cells_excluded": int(group["n_cells"].sum()),
            "median_structural_magnitude": float(group["structural_magnitude"].median()),
            "median_immune_magnitude_of_those_clusters":
                float(group["immune_magnitude"].median()),
            "median_structural_standout": float(group["structural_standout"].median()),
            "median_immune_standout_of_those_clusters":
                float(group["immune_standout"].median()),
            "marker_dynamic_range_on_this_section": dyn.get(str(marker), float("nan")),
        })
    return rows


def cell_evidence(st: dict, members: dict, n: int) -> dict[str, np.ndarray]:
    """Each cell's OWN evidence, on the same statistic the clusters were named on.

    Nothing here is a threshold.  `immune_magnitude` is the cell's largest value-over-section-
    median among the immune canonicals, `structural_magnitude` the same among the structural
    markers, and `argmax_subpanel` the subpanel whose own naming marker the cell reads highest on.
    They exist so that a cluster's call can be checked against the cells inside it: a cluster
    whose cells mostly disagree with its own label is a mixed cluster, and this is how it becomes
    visible instead of hiding behind a majority.

    THE STATISTIC IS THE MAGNITUDE, FOR THE REASON STAGE 1 USES IT.  These are comparisons
    BETWEEN CHANNELS inside one cell, so a within-section rank cannot make them: a cell's place in
    a channel that is background everywhere is still a place.  Measured on the fourteen sections,
    the share of cells called immune that agree with their own cluster's class is 0.747 after this
    correction against 0.662 before it, up on thirteen of the fourteen donors and down on 6505.
    That number is a diagnostic and gates nothing -- it is reported, not aimed at -- and it moves
    for two reasons at once, this statistic and the stage-1 change, which is why it is quoted
    across the correction rather than attributed to either half.
    """
    unit = st["unit_magnitude"]

    def top(markers):
        present = [unit[m] for m in markers if m in unit]
        if not present:
            return np.full(n, np.nan, np.float32), np.full(n, -1, np.int8)
        block = np.vstack(present)
        filled = np.where(np.isnan(block), -np.inf, block)
        idx = filled.argmax(axis=0).astype(np.int8)
        best = filled.max(axis=0)
        return np.where(np.isfinite(best), best, np.nan).astype(np.float32), idx

    primaries = sorted({m for lin in IMMUNE_LINEAGES for m in PRIMARY_IMMUNE[lin]})
    imm, imm_i = top(primaries)
    stc, _ = top(list(STRUCTURAL))
    marker_names = np.array(primaries, dtype=object)
    return {"immune_magnitude": imm, "structural_magnitude": stc,
            "top_immune_marker": np.where(imm_i >= 0, marker_names[np.clip(imm_i, 0, None)], "")}


def per_cell_subpanel(calls: pd.DataFrame, cluster_ids: list[str], codes: np.ndarray,
                      st: dict) -> np.ndarray:
    """Which subpanel each cell reads highest on -- ONE marker per subpanel, never a max of many.

    A max over a subpanel's whole marker set is a max of k, and k differs: the myeloid subpanel
    has five naming markers on this panel and the granulocyte subpanel has one (MPO), so a max-of-
    set comparison would hand myeloid cells for free and make every granulocyte cluster look mixed
    whatever it contains.  Each subpanel is therefore represented by the single marker that spoke
    for it ON THIS CLUSTER -- the one `call_clusters` recorded -- so the per-cell comparison is one
    marker against one marker and the statistic is not a restatement of the panel's shape.
    """
    unit = st["unit_magnitude"]
    names = list(SUBPANEL_FILES)
    by_cluster = calls.set_index("cluster_id")
    out = np.full(len(codes), "", dtype=object)
    for i, cid in enumerate(cluster_ids):
        mask = codes == i
        if not mask.any():
            continue
        rows, labels = [], []
        for name in names:
            marker = str(by_cluster.at[cid, f"marker_{name}"]) if f"marker_{name}" in                 by_cluster.columns else ""
            if marker and marker in unit:
                rows.append(unit[marker][mask])
                labels.append(name)
        if not rows:
            continue
        block = np.vstack(rows)
        filled = np.where(np.isnan(block), -np.inf, block)
        best = filled.max(axis=0)
        won = np.array(labels, dtype=object)[filled.argmax(axis=0)]
        out[mask] = np.where(np.isfinite(best), won, "")
    return out


def heterogeneity(calls: pd.DataFrame, cluster_ids: list[str], codes: np.ndarray,
                  ev: dict, argmax_subpanel: np.ndarray) -> pd.DataFrame:
    """Per cluster, how concentrated its canonical evidence is at the CELL level.

    `cell_frac_immune_lean`  -- cells whose own immune rank beats their own structural rank.  For
                               a cluster called immune this is the share of its cells that agree
                               with that, cell by cell, on a comparison that needs no threshold.
    `cell_frac_agrees_class` -- cells whose own strongest subpanel is the one the cluster was
                               named.  A cluster at 0.9 is one population; a cluster at 0.35 is a
                               mixture wearing a majority label, and the number says which.
    `cell_argmax_entropy`    -- Shannon entropy of the per-cell subpanel argmax within the
                               cluster, normalised by log(5).  0 is one population, 1 is an even
                               mixture of all five.
    """
    lean = np.where(np.isfinite(ev["immune_magnitude"]) & np.isfinite(ev["structural_magnitude"]),
                    ev["immune_magnitude"] > ev["structural_magnitude"], np.nan)
    names = list(SUBPANEL_FILES)
    rows = []
    for i, cid in enumerate(cluster_ids):
        mask = codes == i
        n = int(mask.sum())
        call = calls.loc[calls["cluster_id"] == cid].iloc[0]
        sub = argmax_subpanel[mask]
        counts = pd.Series(sub).value_counts()
        counts = counts[counts.index != ""]
        share = (counts / counts.sum()) if counts.sum() else pd.Series(dtype=float)
        ent = float(-(share * np.log(share)).sum() / np.log(len(names))) if len(share) > 1 else 0.0
        this_lean = lean[mask]
        rows.append({
            "donor_id": call["donor_id"], "cluster_id": cid, "n_cells": n,
            "is_immune": call["is_immune"], "immune_class": call["immune_class"],
            "cell_frac_immune_lean": float(np.nanmean(this_lean)) if n else np.nan,
            "cell_frac_evaluable": float(np.isfinite(this_lean).mean()) if n else np.nan,
            "cell_frac_agrees_class": float((sub == call["immune_class"]).mean()) if n else np.nan,
            "cell_argmax_entropy": ent,
            "cell_argmax_top": str(share.idxmax()) if len(share) else "",
            "cell_argmax_top_frac": float(share.max()) if len(share) else np.nan,
        })
    return pd.DataFrame(rows)


def coherence(donor: str, calls: pd.DataFrame, cluster_ids: list[str], codes: np.ndarray,
              st: dict) -> dict:
    """Internal contradictions, checked WITHOUT reference to what a node should contain.

    None of these compares the answer against an expectation.  Each is a statement that cannot be
    true of the same cells at once, so it is a defect whatever the tissue is.

      cd4_plus_cd8_pct_at_q  CD4-positive % + CD8-positive % among the cells of the T clusters,
                             where positive means "in the top q of this section for this marker".
                             Past 100 % the two together claim more cells than exist -- the
                             contradiction the withdrawn attempt hit at CD4 69 % / CD8 58 %.  The
                             grid is reported so the statement does not rest on one q, and q
                             decides nothing: no call reads it.
      t_cd20_pct             T-cluster cells that are also CD20-positive, and the mirror.
      immune_struct_pct      immune-cluster cells whose own structural rank beats their own immune
                             rank -- permitted by stage 1, which arbitrates on the cluster, so it
                             is reported rather than asserted to be zero.
    """
    unit = st["unit_rank"]
    t_clusters = set(calls.loc[calls["immune_class"] == "T_cell", "cluster_id"])
    b_clusters = set(calls.loc[calls["immune_class"] == "B_Plasma", "cluster_id"])
    imm_clusters = set(calls.loc[calls["is_immune"] == True, "cluster_id"])          # noqa: E712
    ids = np.array(cluster_ids, dtype=object)[codes]
    in_t = np.isin(ids, list(t_clusters))
    in_b = np.isin(ids, list(b_clusters))
    in_imm = np.isin(ids, list(imm_clusters))

    def pct(mask, within):
        return float(100.0 * np.asarray(mask)[within].mean()) if within.any() else float("nan")

    def pos(marker, q):
        u = unit.get(marker)
        if u is None:
            return np.zeros(len(ids), dtype=bool)
        return np.isfinite(u) & (u >= 1.0 - q)

    out = {"donor_id": donor, "n_T_cells": int(in_t.sum()), "n_B_cells": int(in_b.sum()),
           "n_immune_cells": int(in_imm.sum())}
    for q in COHERENCE_Q:
        cd4, cd8 = pos("CD4", q), pos("CD8", q)
        tag = f"{q:g}".replace("0.", "p")
        out[f"cd4_pct_of_T_at_{tag}"] = pct(cd4, in_t)
        out[f"cd8_pct_of_T_at_{tag}"] = pct(cd8, in_t)
        out[f"cd4_and_cd8_pct_of_T_at_{tag}"] = pct(cd4 & cd8, in_t)
        out[f"cd4_plus_cd8_pct_of_T_at_{tag}"] = (out[f"cd4_pct_of_T_at_{tag}"]
                                                  + out[f"cd8_pct_of_T_at_{tag}"])
        # The union.  CD4 % + CD8 % may pass 100 % when many cells are positive for both, which a
        # loose q produces on its own; the union passing 100 % could not happen at all and would
        # be an arithmetic fault rather than a loose definition.  Both are reported so a reader
        # can tell the two apart instead of reading one number.
        out[f"cd4_or_cd8_pct_of_T_at_{tag}"] = (out[f"cd4_plus_cd8_pct_of_T_at_{tag}"]
                                                - out[f"cd4_and_cd8_pct_of_T_at_{tag}"])
    tags = [f"{q:g}".replace("0.", "p") for q in COHERENCE_Q]
    over = [t for t in tags if out[f"cd4_plus_cd8_pct_of_T_at_{t}"] > 100.0]
    out["contradiction_cd4_cd8"] = bool(over)
    out["contradiction_cd4_cd8_at"] = ";".join(over)
    out["union_over_100_at"] = ";".join(
        t for t in tags if out[f"cd4_or_cd8_pct_of_T_at_{t}"] > 100.0 + 1e-9)
    out["t_cd20_pct_of_T_at_p1"] = pct(pos("CD20", 0.10), in_t)
    out["b_cd3e_pct_of_B_at_p1"] = pct(pos("CD3e", 0.10), in_b)
    imm_mag = st["_cell_immune_magnitude"]
    str_mag = st["_cell_structural_magnitude"]
    beaten = np.isfinite(imm_mag) & np.isfinite(str_mag) & (str_mag > imm_mag)
    out["immune_cells_whose_own_structural_magnitude_wins_pct"] = pct(beaten, in_imm)
    return out


# --------------------------------------------------------------------------------------------
# Stage 3 -- the T_NK compartment split by NEGATIVE GATING.
#
# `lineages.LINEAGE_MARKERS["T_NK"] = ("CD3e", "CD8", "CD4", "FOXP3", "CD56")` is ONE COMBINED
# CLASS, and that is a statement about the coarse label, not about the panel's reach.  No channel
# on this panel is NK-specific, and none has to be: NK is read from a COMBINATION, in which the
# absence of a channel carries as much as the presence of one.
#
#     CD3e+ CD4+                   ->  T_CD4
#     CD3e+ CD8+                   ->  T_CD8
#     CD3e+ FOXP3+                 ->  T_reg
#     CD3e- CD8+  or  CD3e- CD56+  ->  NK_ILC
#
# NEGATIVE EVIDENCE IS EVIDENCE.  The load-bearing term is CD3e-NEGATIVE, and it is decided the
# same way this lane decides every other term -- WITHIN THE SECTION, AGAINST THE OTHER CLUSTERS,
# THRESHOLD-FREE.  "Negative" here means AUC <= 0.5 against the rest of this same section: 0.5 is
# not a cutoff, it is the value the statistic takes when a cluster is indistinguishable from the
# rest of the slide, so a cluster at or below it reads no higher on CD3e than the slide at large.
# No intensity, level, quantile or imported constant enters, and nothing is carried between
# sections.
#
# THE RESOLUTION MAY BE SUB-CLUSTER.  A cluster can hold T cells and NK cells at once, so the
# split is attempted INSIDE every cluster of the compartment before the cluster is labelled whole.
# `best_two_group_split` finds the 2-group partition of the cluster's own CD3e ranks that
# separates them best, and the same calculation over EVERY channel of the panel supplies the null:
# any variable splits in two somehow, so CD3e has to split these cells better than the panel's
# typical channel does before the split is read as structure.  A split that is not taken is a real
# answer -- the cluster stays T_NK -- and the numbers behind that answer are reported per cluster.
# --------------------------------------------------------------------------------------------

#: The stage-2 class that carries the combined compartment.  `03_immune_T_cell.csv` is the
#: subpanel that names it, so the compartment arrives here under the name `T_cell`; it is the
#: `T_NK` lineage of `lineages.py` and stage 3 is what tells the two apart inside it.
T_NK_CLASS = "T_cell"

#: The negative term.  It is one marker by design: what makes a cell NK here is that it does NOT
#: read CD3e while it does read CD8 or CD56.
T_NK_NEGATIVE = "CD3e"

#: The three CD3e-POSITIVE terms of the maintainer's table, resolved by argmax so that they are
#: mutually exclusive and so that they SUM to the CD3e-positive cells.  CD4 % + CD8 % among these
#: cells therefore cannot pass 100 % arithmetically -- see `t_nk_summary`.
#:
#: FOXP3 IS `lineages.SUPPORTING_ONLY` AND THAT MATTERS.  It may not admit a cell to the immune
#: compartment on its own and it does not: every cell reaching this argmax is already inside a
#: cluster admitted on CD3e/CD4/CD8.  What it does here is separate cells that are ALREADY T
#: cells, which is the one job `lineages` leaves it -- and `T_reg` is written as the CD3e-positive
#: cells that read highest on FOXP3, which is a measurement of this panel and not a claim that
#: they are a lineage apart from CD4.
T_NK_POSITIVE_SUBSETS = (("T_CD8", "CD8"), ("T_CD4", "CD4"), ("T_reg", "FOXP3"))

#: The positive half of the NK term.  Either is enough; both are reported.
NK_POSITIVE = ("CD8", "CD56")

#: What a cell of the compartment is called when the data does not separate it.  It is the coarse
#: label it already had, kept deliberately -- not a failure and not an empty string.
T_NK_RESIDUAL = "T_NK_unresolved"

T_NK_CLASSES = ("T_CD4", "T_CD8", "T_reg", "NK_ILC", T_NK_RESIDUAL)


def best_two_group_split(x: np.ndarray) -> tuple[float, float, int]:
    """The 2-group split of one variable that separates these cells best, and how cleanly.

    Exact rather than iterative: the values are sorted once and every cut is scored by its
    between-group sum of squares, which is the whole of what a 2-group partition of a single
    variable can be scored on.  NOTHING IS IMPORTED -- the cut is a property of these cells alone,
    and `eta2`, between-group over total sum of squares, is unit-free, so it says how cleanly the
    variable splits THIS cluster without reference to a level, a quantile or another section.

    Returns ``(cut, eta2, n_below)``.  The cut falls between two distinct observed values, so a
    run of ties is never split down the middle.  A variable with no variation here returns eta2 0
    and splits nothing.
    """
    v = np.sort(np.asarray(x, dtype=float))
    n = v.size
    if n < 2:
        return float("nan"), float("nan"), 0
    total = float(((v - v.mean()) ** 2).sum())
    if not (total > 0):
        return float("nan"), 0.0, 0
    cs = np.cumsum(v)
    k = np.arange(1, n, dtype=float)
    between = k * (n - k) / n * (cs[:-1] / k - (cs[-1] - cs[:-1]) / (n - k)) ** 2
    boundary = v[:-1] != v[1:]
    if not boundary.any():
        return float("nan"), 0.0, 0
    scored = np.where(boundary, between, -np.inf)
    j = int(np.argmax(scored))
    return float((v[j] + v[j + 1]) / 2.0), float(between[j] / total), int(j + 1)


def two_group_rank_test(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """AUC of a over b with a tie-corrected rank-sum p, on the pooled ranks of a and b only.

    The same statistic `standout` uses, one level down: there the comparison is a cluster against
    the rest of its section, here it is one side of a split against the other side of the same
    split.  Cells carrying no measurement are dropped rather than ranked at the bottom.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    na, nb = a.size, b.size
    if na < 1 or nb < 1:
        return float("nan"), float("nan")
    pooled = np.concatenate([a, b])
    ranks = rankdata(pooled)
    stat = float(ranks[:na].sum()) - na * (na + 1.0) / 2.0
    auc = stat / (na * nb)
    _, tie_counts = np.unique(pooled, return_counts=True)
    tie_term = float(((tie_counts ** 3) - tie_counts).sum())
    n = na + nb
    var = na * nb / 12.0 * ((n + 1.0) - tie_term / (n * (n - 1.0))) if n > 1 else 0.0
    if not (var > 0):
        return float(auc), float("nan")
    z = (stat - na * nb / 2.0) / np.sqrt(var)
    return float(auc), float(2.0 * norm.sf(abs(z)))


def resolve_t_nk(donor: str, calls: pd.DataFrame, cluster_ids: list[str], codes: np.ndarray,
                 st: dict, channels: list[str]) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Split every cluster of the T_NK compartment into T_CD4 / T_CD8 / T_reg / NK_ILC.

    Returns the per-cluster resolution table and the per-cell arrays.  Every cell of the
    compartment leaves here with exactly one of `T_NK_CLASSES`, so the four subsets plus the
    residual SUM to the compartment by construction and the coherence check is arithmetic.

    EVERY COMPARISON IS OFFSET-CORRECTED, at whichever level it is made.

      cluster level   S = AUC - the cluster's median AUC over the whole panel, as stage 1 and
                      stage 2 already use.  CD3e-NEGATIVE is S <= 0: the cluster reads CD3e no
                      higher than it reads the panel at large.  The raw AUC against the rest of
                      the section is reported beside it, because on a cluster of ordinary
                      brightness the two say the same thing -- but on a bright cluster the raw
                      number says only that the cluster is bright.
      cell level      r = the cell's place on this section for a marker, minus that cell's own
                      median place over the whole panel.  Without it a split on CD3e sorts the
                      cluster by brightness and the dim side reads lower on CD8 and CD56 too,
                      which is measured: on all 32 compartment clusters of these 14 sections the
                      raw comparison puts the CD3e-low side BELOW the CD3e-high side on CD8, AUC
                      0.193 to 0.513.

    Neither correction is a threshold and neither imports anything: both are medians of this
    section's own ranks.  A constant per cell cancels inside one cell, so the CD4/CD8/FOXP3
    argmax is identical with and without it; only comparisons between cells move.

    THE ORDER IS FINEST-FIRST.  A sub-cluster split is attempted before the cluster is labelled
    whole, because a mixed cluster labelled whole would hide exactly the population this is for.
    """
    index = {m: j for j, m in enumerate(st["markers"])}
    unit = st["unit_rank"]
    cell_offset = st["cell_offset"]
    auc, rel, enr = st["auc"], st["rel"], st["enriched"]
    alpha = float(st["alpha"])

    def r_of(m: str) -> np.ndarray | None:
        """The section's rank for marker m with each cell's own panel-wide level taken out."""
        u = unit.get(m)
        return None if u is None else (u.astype(np.float32) - cell_offset)

    labels = np.full(len(codes), "", dtype=object)
    modes = np.full(len(codes), "", dtype=object)
    margins = np.full(len(codes), np.nan, dtype=np.float32)

    compartment = set(calls.loc[calls["immune_class"] == T_NK_CLASS, "cluster_id"])
    subset_markers = [(lab, mk) for lab, mk in T_NK_POSITIVE_SUBSETS if mk in unit]
    r_cd3e = r_of(T_NK_NEGATIVE)
    r_nk = {m: r_of(m) for m in NK_POSITIVE}

    # THE GATE IS ITS OWN FAMILY OF TESTS, AND IT IS CORRECTED OVER THAT FAMILY.  `st["alpha"]` is
    # Bonferroni over every cluster x marker test of stage 1 -- on these sections 952 to 1,344 of
    # them, alpha 3.7e-5 to 5.3e-5.  Stage 3 runs one test per compartment cluster per NK term,
    # which is 2 to 8 tests on a section, so borrowing stage 1's alpha would correct for a family
    # two to three orders of magnitude larger than the one being run.  The count is fixed BEFORE
    # any test is performed, so no result can change it.
    n_gate_tests = max(len(compartment) * sum(1 for m in NK_POSITIVE if m in unit), 1)
    gate_alpha = ALPHA / n_gate_tests

    def a_of(i, m):
        j = index.get(m)
        return float(auc[i][j]) if j is not None and np.isfinite(auc[i][j]) else float("nan")

    def s_of(i, m):
        j = index.get(m)
        return float(rel[i][j]) if j is not None and np.isfinite(rel[i][j]) else float("nan")

    def e_of(i, m):
        j = index.get(m)
        return bool(enr[i][j]) if j is not None else False

    def t_subset(sel: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """argmax over the three CD3e-positive terms, on this section's own ranks.

        The cell's own offset is a constant across the three, so it cancels and the argmax is
        the same whether it is subtracted or not -- which is why it is not subtracted here.
        """
        k = int(sel.sum())
        if not subset_markers:
            return np.full(k, T_NK_RESIDUAL, dtype=object), np.full(k, np.nan)
        block = np.vstack([unit[mk][sel] for _, mk in subset_markers]).astype(float)
        filled = np.where(np.isnan(block), -np.inf, block)
        best = filled.max(axis=0)
        names = np.array([lab for lab, _ in subset_markers], dtype=object)[filled.argmax(axis=0)]
        second = np.sort(filled, axis=0)[-2] if filled.shape[0] > 1 else np.full(best.shape, -np.inf)
        good = np.isfinite(best)
        return (np.where(good, names, T_NK_RESIDUAL),
                np.where(good & np.isfinite(second), best - second, np.nan))

    rows = []
    for i, cid in enumerate(cluster_ids):
        if cid not in compartment:
            continue
        mask = codes == i
        n = int(mask.sum())

        cd3e_auc, cd3e_s = a_of(i, T_NK_NEGATIVE), s_of(i, T_NK_NEGATIVE)
        cd3e_enriched = e_of(i, T_NK_NEGATIVE)
        # THE NEGATIVE TERM AT THE CLUSTER LEVEL.  S <= 0 is "this cluster reads CD3e no higher
        # than it reads the panel at large"; 0 is where S sits for a marker that carries none of
        # the cluster's own level, so this is a comparison and not a cutoff.
        cd3e_at_or_below = bool(np.isfinite(cd3e_s) and cd3e_s <= 0.0)
        cd3e_below_section = bool(np.isfinite(cd3e_auc) and cd3e_auc <= AUC_NULL)
        nk_enr = {m: e_of(i, m) for m in NK_POSITIVE}
        cluster_is_nk = cd3e_at_or_below and any(nk_enr.values())

        row = {"donor_id": donor, "cluster_id": cid, "n_cells": n,
               "compartment_class": T_NK_CLASS,
               "cd3e_S_vs_its_own_panel_level": cd3e_s,
               "cd3e_auc_vs_rest_of_section": cd3e_auc,
               "cd3e_enriched": cd3e_enriched,
               "cd3e_at_or_below_its_own_panel_level": cd3e_at_or_below,
               "cd3e_at_or_below_the_rest_of_its_section": cd3e_below_section,
               "stage1_section_bonferroni_alpha": alpha,
               "gate_bonferroni_alpha": gate_alpha, "n_gate_tests_on_this_section": n_gate_tests}
        for m in NK_POSITIVE:
            low = m.lower()
            row[f"{low}_S_vs_its_own_panel_level"] = s_of(i, m)
            row[f"{low}_auc_vs_rest_of_section"] = a_of(i, m)
            row[f"{low}_enriched"] = nk_enr[m]
            row[f"{low}_measured_on_this_section"] = bool(
                m in unit and np.isfinite(unit[m][mask]).any())
        row["cluster_level_nk_pattern"] = bool(cluster_is_nk)

        # --- the sub-cluster gate, run on every cluster of the compartment ----------------------
        #
        # NOT A 2-GROUP SPLIT.  The best 2-group partition of a cluster maximises between-group
        # sum of squares, which on a largely unimodal cluster cuts it near the middle: measured on
        # 6523_c00 it put 44,466 of 101,138 cells on the low side, and both sides read CD3e ABOVE
        # their own panel level (+0.060 and +0.234).  A minority population cannot be found that
        # way, and NK cells in a T zone are a minority.
        #
        # THE GATE IS THE PRESENCE/ABSENCE PATTERN, PER CELL, on the same offset-corrected ranks:
        #
        #     r_CD3e <= 0                              this cell reads CD3e no higher than it
        #                                              reads its own panel -- CD3e-NEGATIVE
        #     r_CD8 > 0  or  r_CD56 > 0                and it does read one of these above it
        #
        # 0 is not a cutoff: it is where a marker sits for a cell that carries none of that cell's
        # own level, exactly as S > 0 is what "enriched" means for a cluster one level up.  The
        # best 2-group split is still COMPUTED and reported -- `split_eta2_cd3e` beside
        # `eta2_at_the_zero_cut` says whether the zero cut falls anywhere near a natural valley --
        # but it decides nothing.
        #
        # THE GROUP-LEVEL CHECK IS NOT CIRCULAR.  The two sides are formed by CD3e ALONE; the
        # check is then whether the CD3e-negative side reads higher than the CD3e-positive side of
        # the SAME cluster on CD8 or CD56, at the section's own Bonferroni alpha.  Only when it
        # does are that side's CD8/CD56-positive cells named NK_ILC.
        cut = eta2 = null_eta2 = eta2_zero = float("nan")
        neg = np.zeros(n, dtype=bool)
        pos = np.zeros(n, dtype=bool)
        lo_med = hi_med = lo_med_raw = float("nan")
        nk_marker, nk_auc, nk_p = "", float("nan"), float("nan")
        per_nk_marker: dict[str, tuple[float, float]] = {}
        split_separates = nk_side_ok = False
        n_nk_positive = 0
        if r_cd3e is not None:
            cu = r_cd3e[mask].astype(float)
            fin = np.isfinite(cu)
            if int(fin.sum()) >= 2:
                cut, eta2, _ = best_two_group_split(cu[fin])
                nulls = []
                for m in channels:
                    rm = r_of(m)
                    if rm is None:
                        continue
                    v = rm[mask]
                    v = v[np.isfinite(v)]
                    if v.size >= 2:
                        nulls.append(best_two_group_split(v)[1])
                nulls = [v for v in nulls if np.isfinite(v)]
                null_eta2 = float(np.median(nulls)) if nulls else float("nan")
                split_separates = bool(np.isfinite(eta2) and np.isfinite(null_eta2)
                                       and eta2 > null_eta2)
                neg[:] = fin & (cu <= 0.0)
                pos[:] = fin & (cu > 0.0)
                if neg.any() and pos.any():
                    lo_med, hi_med = float(np.median(cu[neg])), float(np.median(cu[pos]))
                    lo_med_raw = float(np.median(unit[T_NK_NEGATIVE][mask][neg]))
                    grand = float(cu[fin].mean())
                    total = float(((cu[fin] - grand) ** 2).sum())
                    a_, b_ = cu[neg], cu[pos]
                    between = (a_.size * (a_.mean() - grand) ** 2
                               + b_.size * (b_.mean() - grand) ** 2)
                    eta2_zero = float(between / total) if total > 0 else float("nan")
                    for m in NK_POSITIVE:
                        rm = r_nk.get(m)
                        if rm is None:
                            continue
                        block = rm[mask]
                        a, pv = two_group_rank_test(block[neg], block[pos])
                        # Both terms are reported whatever they say; the gate reads the stronger,
                        # and 'stronger' is decided on the effect size before any p is looked at.
                        per_nk_marker[m] = (a, pv)
                        if np.isfinite(a) and (not np.isfinite(nk_auc) or a > nk_auc):
                            nk_marker, nk_auc, nk_p = m, a, pv
                    nk_side_ok = bool(np.isfinite(nk_auc) and nk_auc > AUC_NULL
                                      and np.isfinite(nk_p) and nk_p < gate_alpha)
        # HOW ENTANGLED THE TERMS ARE ON THIS CLUSTER.  A negative gate can only work if CD3e and
        # the NK terms are separable within the cluster at the cell level.  Where they are not --
        # where a cell reading high on CD3e reads high on CD56 too, which dense membrane staining
        # in a T zone can produce -- the gate refuses, and this is the number that says so.  It
        # is reported, never applied: nothing is called or refused on it.
        for m in NK_POSITIVE:
            rm = r_nk.get(m)
            row[f"within_cluster_corr_cd3e_{m.lower()}"] = float("nan")
            if rm is None or r_cd3e is None:
                continue
            a_, b_ = r_cd3e[mask].astype(float), rm[mask].astype(float)
            ok = np.isfinite(a_) & np.isfinite(b_)
            if int(ok.sum()) >= 3:
                ra, rb = rankdata(a_[ok]), rankdata(b_[ok])
                sd = ra.std() * rb.std()
                if sd > 0:
                    row[f"within_cluster_corr_cd3e_{m.lower()}"] = float(
                        ((ra - ra.mean()) * (rb - rb.mean())).mean() / sd)

        # the positive half of the NK term, per cell, among the CD3e-negative cells only
        nk_cell = np.zeros(n, dtype=bool)
        for m in NK_POSITIVE:
            rm = r_nk.get(m)
            if rm is None:
                continue
            v = rm[mask]
            nk_cell |= neg & np.isfinite(v) & (v > 0.0)
        n_nk_positive = int(nk_cell.sum())
        gate_ran = bool(neg.any() and pos.any() and nk_marker)
        split_accepted = bool(gate_ran and nk_side_ok and n_nk_positive > 0)

        row.update({
            "n_cd3e_negative": int(neg.sum()), "n_cd3e_positive": int(pos.sum()),
            "n_cd3e_unmeasured": int(n - neg.sum() - pos.sum()),
            "frac_cd3e_negative": float(neg.sum() / n) if n else float("nan"),
            "cd3e_negative_median_offset_corrected_rank": lo_med,
            "cd3e_positive_median_offset_corrected_rank": hi_med,
            "cd3e_negative_median_raw_unit_rank": lo_med_raw,
            "eta2_at_the_zero_cut": eta2_zero,
            "split_eta2_cd3e": eta2,
            "split_cut_on_cd3e_offset_corrected_rank": cut,
            "split_eta2_panel_median": null_eta2,
            "split_separates_better_than_the_panel_typically_does": split_separates,
            "nk_side_marker": nk_marker,
            "nk_side_auc_negative_over_positive": nk_auc, "nk_side_p": nk_p,
            **{f"nk_side_auc_{m.lower()}": per_nk_marker.get(m, (np.nan, np.nan))[0]
               for m in NK_POSITIVE},
            **{f"nk_side_p_{m.lower()}": per_nk_marker.get(m, (np.nan, np.nan))[1]
               for m in NK_POSITIVE},
            "nk_side_stands_out": nk_side_ok,
            "n_cd3e_negative_and_cd8_or_cd56_positive": n_nk_positive,
            "gate_could_run": gate_ran,
            "split_accepted": split_accepted,
        })

        # --- the verdict, finest first ----------------------------------------------------------
        cell_labels = np.full(n, T_NK_RESIDUAL, dtype=object)
        cell_margin = np.full(n, np.nan)
        if split_accepted:
            mode = "sub_cluster_gate_on_CD3e"
            reason = (
                f"{int(neg.sum()):,} of {n:,} cells read CD3e at or below their own panel level "
                f"(median {lo_med:+.3f} against {hi_med:+.3f} on the CD3e-positive side); that "
                f"side reads higher than the CD3e-positive side on {nk_marker} at AUC "
                f"{nk_auc:.3f}, p {nk_p:.3g}, and {n_nk_positive:,} of them also read CD8 or "
                f"CD56 above their own panel level")
            cell_labels[nk_cell] = "NK_ILC"
            # A cell carrying no CD3e measurement cannot be gated on it at all, so it is on
            # neither side: it stays T_NK_unresolved rather than being swept onto the
            # CD3e-positive side.  An unavailable measurement is not a negative one.  So does a
            # CD3e-negative cell that reads neither CD8 nor CD56 above its own level.
            pos_mask = np.zeros(len(codes), dtype=bool)
            pos_mask[np.flatnonzero(mask)[pos]] = True
            names, marg = t_subset(pos_mask)
            cell_labels[pos] = names
            cell_margin[pos] = marg
        elif cluster_is_nk:
            mode = "whole_cluster_CD3e_at_or_below_its_panel_level"
            won = [m for m in NK_POSITIVE if nk_enr[m]]
            reason = (f"the whole cluster reads CD3e at S {cd3e_s:+.3f} -- at or below its own "
                      f"panel level -- while {', '.join(won)} stands out on it")
            cell_labels[:] = "NK_ILC"
        elif cd3e_enriched:
            mode = "whole_cluster_CD3e_enriched"
            names, marg = t_subset(mask)
            cell_labels[:] = names
            cell_margin[:] = marg
            # The CD3e-negative cells of this cluster are NOT named NK: the gate ran and its
            # group-level check failed, so what it would have named is reported, not written.
            reason = (
                f"CD3e stands out for this whole cluster at S {cd3e_s:+.3f} (AUC {cd3e_auc:.3f}). "
                f"{int(neg.sum()):,} of {n:,} cells do read it at or below their own panel level, "
                f"but that side does NOT read higher than the CD3e-positive side on "
                f"{nk_marker or 'CD8 or CD56'} (AUC {nk_auc:.3f}, p {nk_p:.3g} against alpha "
                f"{gate_alpha:.3g}), so the {n_nk_positive:,} cells the positive term would have "
                f"admitted are not named NK.  Every cell is called on the argmax of "
                f"{', '.join(mk for _, mk in subset_markers)}")
        else:
            mode = "did_not_resolve"
            reason = (
                f"CD3e is neither enriched nor at or below this cluster's own panel level "
                f"(S {cd3e_s:+.3f}, AUC {cd3e_auc:.3f}), and the per-cell gate did not separate "
                f"these cells: {int(neg.sum()):,} of {n:,} read CD3e at or below their own level "
                f"and that side reads {nk_marker or 'neither CD8 nor CD56'} at AUC {nk_auc:.3f}, "
                f"p {nk_p:.3g} against alpha {gate_alpha:.3g}.  The compartment stays combined "
                f"here, which is a measurement of this cluster and not a missing one")

        labels[mask] = cell_labels
        modes[mask] = mode
        margins[mask] = cell_margin.astype(np.float32)
        row["resolution_mode"] = mode
        row["resolution_reason"] = reason
        counts = pd.Series(cell_labels).value_counts()
        for name in T_NK_CLASSES:
            row[f"n_{name}"] = int(counts.get(name, 0))
        row["median_subset_argmax_margin"] = (float(np.nanmedian(cell_margin))
                                              if np.isfinite(cell_margin).any() else float("nan"))
        assert sum(row[f"n_{name}"] for name in T_NK_CLASSES) == n, \
            f"{donor}/{cid}: the subsets lost a cell"
        rows.append(row)

    frame = pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=["donor_id", "cluster_id", "n_cells", "resolution_mode"])
    return frame, {"t_nk_class": labels, "t_nk_resolution_mode": modes,
                   "t_nk_subset_argmax_margin": margins}


def t_nk_summary(donor: str, resolution: pd.DataFrame, compartment_cells: int) -> dict:
    """Per donor: the compartment, its four subsets, and whether NK could be decided at all.

    NK_ILC IS NEVER WRITTEN AS 0 WHEN IT WAS NOT DECIDABLE.  A zero here is a measurement only
    when the negative gate could actually run -- CD3e measured, and CD8 or CD56 measured, on at
    least one cluster of this donor's compartment.  When it could not, the count is left EMPTY
    with the reason beside it, because an unavailable measurement is not zero.
    """
    out = {"donor_id": donor, "n_T_NK_compartment": int(compartment_cells),
           "n_clusters_in_compartment": int(len(resolution))}
    if len(resolution):
        evaluable = resolution.loc[resolution["gate_could_run"].astype(bool)]
        modes = resolution["resolution_mode"].value_counts()
    else:
        evaluable, modes = resolution, pd.Series(dtype=int)
    n_evaluable = int(evaluable["n_cells"].sum()) if len(evaluable) else 0
    out["n_cells_where_the_CD3e_gate_could_run"] = n_evaluable
    for name in T_NK_CLASSES:
        col = f"n_{name}"
        out[name] = int(resolution[col].sum()) if col in resolution.columns else 0
    MODES = ("sub_cluster_gate_on_CD3e", "whole_cluster_CD3e_at_or_below_its_panel_level",
             "whole_cluster_CD3e_enriched", "did_not_resolve")
    for mode in MODES:
        out[f"clusters_{mode}"] = int(modes.get(mode, 0))
    n_clusters = int(len(resolution))
    nothing_resolved = bool(n_clusters and int(modes.get("did_not_resolve", 0)) == n_clusters)
    if not compartment_cells:
        out["NK_ILC"] = ""
        out["nk_ilc_state"] = ("unavailable -- this donor has no T_NK compartment on this "
                               "section, so there is nothing here to resolve or to count")
    elif n_evaluable == 0:
        out["NK_ILC"] = ""
        out["nk_ilc_state"] = (
            "unavailable -- no cluster of this donor's T_NK compartment has both a CD3e-negative "
            "and a CD3e-positive side with CD8 or CD56 measured, so the gate could not be run at "
            "all.  This is NOT zero NK cells.")
    elif nothing_resolved:
        out["NK_ILC"] = ""
        out["nk_ilc_state"] = (
            f"unavailable -- every one of this donor's {n_clusters} compartment clusters came "
            f"back did_not_resolve, so the compartment stays combined and NK was never decided "
            f"one way or the other.  This is NOT zero NK cells; t_nk_resolution.csv carries "
            f"what was measured on each of them.")
    else:
        out["nk_ilc_state"] = "measured"
    total = sum(int(resolution[f"n_{n}"].sum()) if f"n_{n}" in resolution.columns else 0
                for n in T_NK_CLASSES)
    out["subsets_sum_to_the_compartment"] = bool(total == int(compartment_cells))
    cd3e_pos = out["T_CD4"] + out["T_CD8"] + out["T_reg"]
    out["n_CD3e_positive"] = cd3e_pos
    out["pct_CD4_of_CD3e_positive"] = 100.0 * out["T_CD4"] / cd3e_pos if cd3e_pos else float("nan")
    out["pct_CD8_of_CD3e_positive"] = 100.0 * out["T_CD8"] / cd3e_pos if cd3e_pos else float("nan")
    out["pct_CD4_plus_CD8_of_CD3e_positive"] = (out["pct_CD4_of_CD3e_positive"]
                                                + out["pct_CD8_of_CD3e_positive"])
    return out


# --------------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--stage", choices=("standout", "all"), default="all")
    args = ap.parse_args()

    warnings.filterwarnings("ignore", r"All-NaN", category=RuntimeWarning)
    warnings.filterwarnings("ignore", r"Mean of empty slice", category=RuntimeWarning)
    warnings.filterwarnings("ignore", r"invalid value encountered", category=RuntimeWarning)

    REF.mkdir(parents=True, exist_ok=True)
    every = declared_written_donors()
    have_clusters = clustered_donors()
    if not have_clusters:
        print("The clustering track has written no assignment under "
              f"{CLUSTERS/'assignments'}.  This lane consumes the clusters and does not make "
              "them, so the run is PENDING that track.  Nothing was written.")
        return 0

    columns = pln_columns()
    channels = panel_channels()
    members = subpanel_membership(columns)
    naming = naming_markers(members)
    primaries = sorted({m for lin in IMMUNE_LINEAGES for m in PRIMARY_IMMUNE[lin]})
    nk_all, nk_any, nk_none = nk_licence()
    decision = sorted(set(primaries) | set(STRUCTURAL) | set(SUPPORTING_IMMUNE)
                      | {m for v in naming.values() for m in v}
                      | set(nk_all) | {m for g in nk_any for m in g} | set(nk_none))
    # Everything the export carries is measured.  The decision markers name the clusters; the rest
    # describe them -- CD21 for the follicular dendritic network, Ki67/PCNA/Bcl_2 for the germinal
    # centre, TCF_1/TOX/PD_1 for T-zone state -- and every channel enters the cluster's panel-wide
    # offset, which is why that offset cannot be steered by the choice of decision markers.
    described = sorted(set(channels) | set(decision))
    absent = [m for m in decision if m not in columns]

    print("NO PANCREAS NUMBER IS READ, AND NO LEVEL IS FITTED AT ALL.  Every statistic below is a"
          "\ncomparison between the clusters of one section, on that section's own ranks.\n")
    print("STAGE 1 -- immune is a union of lineage canonicals, arbitrated against the structural set.")
    for lin in IMMUNE_LINEAGES:
        print(f"   {lin:12s} {', '.join(PRIMARY_IMMUNE[lin])}")
    print(f"   {'supporting':12s} {', '.join(SUPPORTING_IMMUNE)}  (never alone)")
    print(f"   {'structural':12s} {', '.join(STRUCTURAL)}")
    for marker, why in STRUCTURAL_DROPPED.items():
        print(f"   structural dropped: {marker} -- {why}")
    print(f"   wanted by the method but NOT on this export (unavailable, never zero): "
          f"{', '.join(absent) or 'none'}")
    print("\nSTAGE 2 -- subpanel membership from core/panels/data/identity_core:")
    for name, hints in members.items():
        keep = naming[name]
        supp = [m for m in hints if m in LIN.SUPPORTING_ONLY]
        lic = SUPPORTING_ONLY_SUBPANEL.get(name)
        print(f"   {name:12s} {', '.join(keep) or '(no naming marker on this panel)'}"
              + (f"   [supporting: {', '.join(supp)}]" if supp else "")
              + (f"   [licensed by {lic}: {nk_all} + any{nk_any} - {nk_none}]" if lic else ""))
    print(f"\nENRICHED  AUC > {AUC_NULL}, rank-sum p < {ALPHA} Bonferroni-corrected over every "
          f"cluster x marker test on that\n          section, AND above the cluster's own median "
          f"AUC over all {len(channels)} channels.\n")
    print("ARBITRATION  between two channels of ONE cluster, on the MAGNITUDE -- that cluster's "
          "median over\n             the section's median on the channel.  The rank statistic "
          "admits the evidence; it\n             cannot compare two channels, because it is "
          "scale-free and the ordering of background\n             is still an ordering.  See "
          "the module docstring for the measurement.\n")
    print("STAGE-2 REQUIREMENT  from subpanel.immune_rules(), the same table nk_licence() reads:")
    for name, need in subpanel_requirements().items():
        print(f"   {name:12s} {' + '.join(need) or '(the table requires nothing of it)'}")
    print()

    donors = args.donors or have_clusters

    # ------------------------------------------------------------------------------------------
    # D4.  A PARTIAL RUN MUST NEVER OVERWRITE A COHORT-WIDE ARTIFACT.  `--donors X` and
    # `--stage standout` both used to write `cluster_calls.csv` and `per_donor_separation.csv`
    # straight over the cohort's own copies -- measured 2026-09-17, 293 cluster rows down to 19
    # and 14 donor rows down to 1 -- so a file named for the cohort silently held one donor.  A
    # run that does not cover every clustered donor and every stage is now SCOPED: it writes the
    # same table names under `_partial/<scope>/`, and the cohort files are left untouched.
    # ------------------------------------------------------------------------------------------
    covers_every_donor = set(donors) >= set(have_clusters)
    covers_every_stage = args.stage == "all"
    partial = not (covers_every_donor and covers_every_stage)
    if partial:
        scope = ("stage-" + args.stage if covers_every_donor
                 else "donors-" + "_".join(sorted(map(str, donors)))
                 + ("" if covers_every_stage else "-stage-" + args.stage))
        COHORT = OUT / "_partial" / scope[:120]
        COHORT.mkdir(parents=True, exist_ok=True)
        print(f"PARTIAL RUN -- {len(donors)} of {len(have_clusters)} clustered donors, stage "
              f"{args.stage!r}.  The cohort-wide tables under {OUT} are NOT touched; this run "
              f"writes its own copies to {COHORT}.\n")
        (COHORT / "PARTIAL_RUN.json").write_text(json.dumps({
            "donors_run": sorted(map(str, donors)),
            "donors_with_a_clustering": have_clusters,
            "stage": args.stage,
            "why": "a partial run does not produce a cohort-wide table, so it is not allowed to "
                   "write one; these files describe ONLY the donors and stages listed here",
        }, indent=1))
    else:
        COHORT = OUT
    REF_OUT = COHORT / "_reference"
    REF_OUT.mkdir(parents=True, exist_ok=True)

    status: dict[str, dict] = {}
    counts, cluster_rows, het_rows, coh_rows, sep_rows, veto_rows = [], [], [], [], [], []
    tnk_rows, tnk_donor_rows = [], []

    for donor in sorted(set(declared_manifest_donors()) - set(every)):
        status[donor] = {"state": "withheld", "reason": withheld_reason(donor)}
    for donor in sorted(set(every) - set(have_clusters)):
        status[donor] = {"state": "unavailable",
                         "reason": "the clustering track produced no cluster for this donor; "
                                   "its QuPath export carries no intensity measurement, so there "
                                   "is nothing to cluster and nothing to name -- never 0 immune "
                                   "cells"}

    t0 = time.time()
    for donor in donors:
        frame = read_pln(donor, ["object_id", "qc_analysis_eligible"] + described)
        n_rows = len(frame)
        frame["object_id"] = frame["object_id"].astype(str)
        # A cell with no nucleus is excluded from every calculation.  `qc_analysis_eligible` is
        # the ingest's own verdict and `missing_nucleus_geometry` is one of its reasons.
        frame = frame.loc[frame["qc_analysis_eligible"].to_numpy(dtype=bool)].reset_index(drop=True)

        assignment = cluster_assignment(donor)
        merged = frame.merge(assignment, on="object_id", how="inner")
        if len(merged) == 0:
            status[donor] = {"state": "unavailable", "n_rows": n_rows, "n_eligible": len(frame),
                             "reason": "no analysis-eligible cell of this donor carries a cluster"}
            print(f"\n{donor}: UNAVAILABLE -- {status[donor]['reason']}")
            continue
        cluster_ids = sorted(merged["cluster_id"].unique())
        codes = pd.Categorical(merged["cluster_id"], categories=cluster_ids).codes.astype(np.int64)
        sizes = np.bincount(codes, minlength=len(cluster_ids))

        profile = cluster_profile(donor)
        joined = int(profile["n_cells"].sum()) if len(profile) else len(merged)
        print(f"\n{donor}: {n_rows:,} rows, {len(frame):,} analysis-eligible, "
              f"{len(merged):,} carrying a cluster over {len(cluster_ids)} clusters "
              f"(clustering track says {joined:,})   [{time.time()-t0:.0f}s]", flush=True)

        markers = [m for m in described if m in merged.columns]
        st = standout(merged[markers], codes, len(cluster_ids), markers,
                      offset_over=[m for m in channels if m in markers])
        # The magnitude the rank statistic discards, on the same clusters and the same cells.
        # Stage 1 and stage 2 arbitrate on it; `dynamic_range` is its per-channel maximum and is
        # reported beside every exclusion a structural marker causes.
        st["magnitude"] = magnitude(merged[markers], codes, len(cluster_ids), markers)
        dyn = dynamic_range(merged[[m for m in decision if m in markers]], codes,
                            len(cluster_ids), [m for m in decision if m in markers])
        ev = cell_evidence(st, members, len(merged))
        st["_cell_immune_magnitude"] = ev["immune_magnitude"]
        st["_cell_structural_magnitude"] = ev["structural_magnitude"]

        index = {m: j for j, m in enumerate(markers)}
        n_enriched_canonical = int(sum(
            any(st["enriched"][i][index[m]] for m in primaries if m in index)
            for i in range(len(cluster_ids))))
        n_enriched_structural = int(sum(
            any(st["enriched"][i][index[m]] for m in STRUCTURAL if m in index)
            for i in range(len(cluster_ids))))
        if n_enriched_canonical == 0 and n_enriched_structural == 0:
            status[donor] = {"state": "unavailable", "n_rows": n_rows, "n_eligible": len(frame),
                             "reason": f"no cluster on this section is enriched for any lineage "
                                       f"marker over the other clusters "
                                       f"({len(cluster_ids)} clusters, {st['n_tests']:,} tests, "
                                       f"alpha {st['alpha']:.3g}) -- the clusters do not separate "
                                       f"on identity, so nothing is forced"}
            print(f"   UNAVAILABLE -- {status[donor]['reason']}")
            continue

        calls = call_clusters(donor, cluster_ids, sizes, st, members)
        cluster_rows.append(calls)
        argmax_subpanel = per_cell_subpanel(calls, cluster_ids, codes, st)
        het = heterogeneity(calls, cluster_ids, codes, ev, argmax_subpanel)
        het_rows.append(het)

        # --- stage 3: the T_NK compartment, split by negative gating on CD3e -------------------
        resolution, tnk = resolve_t_nk(donor, calls, cluster_ids, codes, st, channels)
        n_compartment = int((calls.loc[calls["immune_class"] == T_NK_CLASS, "n_cells"]).sum())
        if len(resolution):
            tnk_rows.append(resolution)
        tnk_summary = t_nk_summary(donor, resolution, n_compartment)
        tnk_donor_rows.append(tnk_summary)

        sep_rows.append({
            "donor_id": donor, "n_clusters": len(cluster_ids), "n_cells": len(merged),
            "n_tests": st["n_tests"], "bonferroni_alpha": st["alpha"],
            "clusters_with_an_enriched_immune_canonical": n_enriched_canonical,
            "clusters_with_an_enriched_structural_marker": n_enriched_structural,
            "panel_offset_min": float(np.nanmin(st["offset"])),
            "panel_offset_max": float(np.nanmax(st["offset"])),
            "panel_offset_spread": float(np.nanmax(st["offset"]) - np.nanmin(st["offset"])),
            "n_channels_in_the_offset": len(st["offset_over"]),
            "max_immune_standout": float(np.nanmax(calls["immune_standout"])),
            "median_immune_standout": float(np.nanmedian(calls["immune_standout"])),
            "max_structural_standout": float(np.nanmax(calls["structural_standout"])),
            "median_separation_margin": float(np.nanmedian(calls["separation_margin"])),
            "max_separation_margin": float(np.nanmax(calls["separation_margin"])),
            "min_separation_margin": float(np.nanmin(calls["separation_margin"])),
            # The margin the arbitration actually read.  `separation_margin` above is the same two
            # markers on the rank statistic, kept so a reader can see both of them at once.
            "median_magnitude_margin": float(np.nanmedian(calls["magnitude_margin"])),
            "max_magnitude_margin": float(np.nanmax(calls["magnitude_margin"])),
            "min_magnitude_margin": float(np.nanmin(calls["magnitude_margin"])),
        })

        veto_rows.extend(structural_veto_rows(donor, calls, dyn))

        # --- cells inherit their cluster's call -------------------------------------------------
        by_cluster = calls.set_index("cluster_id")
        take = lambda col: by_cluster[col].reindex(cluster_ids).to_numpy()[codes]  # noqa: E731
        is_immune = pd.array(take("is_immune"), dtype="boolean")
        cells = pd.DataFrame({
            "object_id": merged["object_id"].to_numpy(),
            "cluster_id": merged["cluster_id"].astype("category"),
            "is_immune": is_immune,
            "immune_class": pd.Categorical(take("immune_class").astype(str)),
            "immune_subtype": pd.Categorical(take("immune_subtype").astype(str)),
            "call_source": "inherited from this cell's cluster",
            "stage1_rule": pd.Categorical(take("stage1_rule").astype(str)),
            "stage2_reason": pd.Categorical(take("stage2_reason").astype(str)),
            "cluster_immune_standout": take("immune_standout").astype(float),
            "cluster_structural_standout": take("structural_standout").astype(float),
            "cell_immune_magnitude": ev["immune_magnitude"],
            "cell_structural_magnitude": ev["structural_magnitude"],
            "cell_top_immune_marker": pd.Categorical(ev["top_immune_marker"].astype(str)),
            "cell_argmax_subpanel": pd.Categorical(argmax_subpanel.astype(str)),
            # Stage 3.  Empty outside the T_NK compartment -- it is a refinement OF that
            # compartment and says nothing about a myeloid or B cell.
            "t_nk_class": pd.Categorical(tnk["t_nk_class"].astype(str)),
            "t_nk_resolution_mode": pd.Categorical(tnk["t_nk_resolution_mode"].astype(str)),
            "t_nk_subset_argmax_margin": tnk["t_nk_subset_argmax_margin"],
        })
        cells["cell_agrees_with_cluster_class"] = pd.array(
            np.where(argmax_subpanel.astype(str) == "", None,
                     argmax_subpanel.astype(str) == take("immune_class").astype(str)),
            dtype="boolean")
        assert cells["object_id"].is_unique, f"{donor}: object_id is not unique"

        target = OUT / f"donor_id={donor}/roi=pln_1"
        target.mkdir(parents=True, exist_ok=True)
        cells.to_parquet(target / "data.parquet", index=False)
        stale = target / "UNAVAILABLE.json"
        if stale.exists():
            # The only file this script removes, and one it wrote itself: a stale unavailability
            # marker beside a real table is how a tree comes to stack two generations.
            stale.unlink()

        n_imm = int((cells["is_immune"] == True).sum())                         # noqa: E712
        n_non = int((cells["is_immune"] == False).sum())                        # noqa: E712
        n_unc = int(cells["is_immune"].isna().sum())
        breakdown = cells.loc[cells["is_immune"] == True, "immune_class"].value_counts()  # noqa
        row = {"donor_id": donor, "n_rows": n_rows, "n_eligible": len(frame),
               "n_clustered": len(merged), "n_clusters": len(cluster_ids),
               "n_immune": n_imm, "n_non_immune": n_non, "n_uncalled": n_unc,
               "pct_immune": 100.0 * n_imm / len(merged)}
        for name in list(SUBPANEL_FILES) + ["Uncalled"]:
            row[f"class_{name}"] = int(breakdown.get(name, 0))
        # Stage 3.  These five sum to class_T_cell, which is the combined T_NK compartment; the
        # NK count is left EMPTY, never 0, when the CD3e gate could not be run at all.
        for name in T_NK_CLASSES:
            row[f"tnk_{name}"] = tnk_summary[name]
        row["nk_ilc_state"] = tnk_summary["nk_ilc_state"]
        counts.append(row)
        coh_rows.append(coherence(donor, calls, cluster_ids, codes, st))
        status[donor] = {"state": "typed", "n_rows": n_rows, "n_eligible": len(frame),
                         "n_clustered": len(merged), "n_clusters": len(cluster_ids),
                         "n_immune": n_imm, "n_non_immune": n_non, "n_uncalled": n_unc}
        print(f"   immune {n_imm:,}  non-immune {n_non:,}  uncalled {n_unc:,}   "
              + "  ".join(f"{k} {v:,}" for k, v in breakdown.items() if v))

    if cluster_rows:
        allc = pd.concat(cluster_rows, ignore_index=True)
        allc.to_csv(COHORT / "cluster_calls.csv", index=False)
        allc.to_parquet(REF_OUT / "cluster_calls.parquet", index=False)
    if het_rows:
        pd.concat(het_rows, ignore_index=True).to_csv(COHORT / "cluster_heterogeneity.csv",
                                                      index=False)
    if sep_rows:
        pd.DataFrame(sep_rows).to_csv(COHORT / "per_donor_separation.csv", index=False)
    if veto_rows:
        pd.DataFrame(veto_rows).sort_values(
            ["donor_id", "n_cells_excluded"], ascending=[True, False]).to_csv(
            COHORT / "structural_veto_accounting.csv", index=False)
    if coh_rows:
        pd.DataFrame(coh_rows).to_csv(REF_OUT / "coherence.csv", index=False)
    if tnk_rows:
        pd.concat(tnk_rows, ignore_index=True).to_csv(COHORT / "t_nk_resolution.csv", index=False)
    if tnk_donor_rows:
        pd.DataFrame(tnk_donor_rows).to_csv(COHORT / "per_donor_t_nk.csv", index=False)
    if args.stage != "all":
        return 0

    pd.DataFrame(counts, columns=[
        "donor_id", "n_rows", "n_eligible", "n_clustered", "n_clusters", "n_immune",
        "n_non_immune", "n_uncalled", "pct_immune",
        *(f"class_{n}" for n in list(SUBPANEL_FILES) + ["Uncalled"]),
        *(f"tnk_{n}" for n in T_NK_CLASSES), "nk_ilc_state"]).to_csv(
        COHORT / "per_donor_counts.csv", index=False)
    pd.DataFrame([{"donor_id": d, "state": s["state"], "n_eligible": s.get("n_eligible", 0),
                   "n_clusters": s.get("n_clusters", 0), "n_immune": s.get("n_immune", ""),
                   "reason": s.get("reason", "")}
                  for d, s in sorted(status.items())]).to_csv(
        COHORT / "per_donor_status.csv", index=False)
    for donor, state in status.items():
        if state["state"] == "typed" or partial:
            # A partial run has not looked at the donors it did not run, so it may not mark any
            # donor unavailable -- that is a cohort-wide statement.
            continue
        target = OUT / f"donor_id={donor}/roi=pln_1"
        target.mkdir(parents=True, exist_ok=True)
        (target / "UNAVAILABLE.json").write_text(json.dumps({"donor_id": donor, **state}, indent=1))

    (REF_OUT / "SUPERSEDED_2026-09-17.json").write_text(json.dumps({
        "what": "outputs of the WITHDRAWN reference-level generation of this lane, left on disk "
                "because nothing under data/ is deleted here, and named so that this tree does "
                "not stack two generations with nothing marking which is which",
        "superseded_files": ["_reference/pln_reference_levels.parquet",
                             "_reference/pln_reference_levels.csv",
                             "_reference/reference_population.csv"],
        "why": "they are p99 levels fitted on `context == 'tissue'`, which on this export is the "
               "whole slide, so each marker called 1 % of the slide positive by construction. "
               "There is no context column and no reference population on a node.  Do not cite "
               "them; the current generation writes cluster_calls.csv and no level at all.",
        "current_outputs": ["cluster_calls.csv", "cluster_heterogeneity.csv",
                            "structural_veto_accounting.csv",
                            "per_donor_separation.csv", "per_donor_counts.csv",
                            "per_donor_status.csv", "donor_id=*/roi=pln_1/data.parquet"],
    }, indent=1))

    (COHORT / "manifest.json").write_text(json.dumps({
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "method": "clusters from the clustering track are named on their own section: the rank "
                  "statistic admits a marker as evidence (ENRICHED against the rest of that "
                  "section) and the MAGNITUDE arbitrates between the admitted markers (stage 1), "
                  "then the core/panels identity_core subpanels resolve within immune on the same "
                  "magnitude and are held to their own rule-table requirement (stage 2); cells "
                  "inherit their cluster's call",
        "clustering_consumed_from": str(CLUSTERS.relative_to(REPO)),
        "no_pancreas_input": "no pancreas cell, level, threshold, scale or quantile is read by "
                             "this script; the pancreas transfer of 2026-09-17 is WITHDRAWN",
        "no_reference_level": "no level is fitted at all.  `subpanel.fit_levels`, "
                              "`subpanel_io.reference_level` and `context == 'tissue'` are "
                              "WITHDRAWN on the pLN: there is no context column and no reference "
                              "population on a node, and none is synthesised",
        "no_absolute_floor": "IMMUNE_REQUIRE_MIN is not applied to any pLN value; it is defined "
                             "against a reference level that does not exist here",
        "statistic": "Mann-Whitney common-language effect size (AUC) of a cluster against the "
                     "rest of its own section, from that section's own ranks, with a "
                     "tie-corrected rank-sum p",
        "enriched_at": {"auc_above": AUC_NULL, "alpha": ALPHA,
                        "correction": "Bonferroni over every cluster x marker test on the section"},
        "arbitration_statistic": "the cluster's median over the section's median on that channel "
                                 "-- the magnitude the rank statistic discards.  A ratio inside "
                                 "one channel and one section: invariant to that channel's scale, "
                                 "dimensionless so two channels of one cluster can be compared, "
                                 "and unavailable (never zero) where the section median is not "
                                 "positive",
        "why_not_the_rank": "AUC is scale-free by design, so a channel with no target on the "
                            "section competes on equal terms with one that marks a population.  "
                            "Measured on the 14 sections before this correction: 1,669,194 of "
                            "5,122,269 cells excluded as non-immune, 429,278 of them decided by "
                            "an epithelial channel in a tissue with no epithelium",
        "references_tested_and_rejected": {
            "cluster_as_the_unit": "recomputing the statistic against the other CLUSTERS of the "
                                   "section, one unit per cluster, moves 17 of 293 clusters and "
                                   "the non-immune share from 32.6% to 31.1%; on donor 115 it "
                                   "moves nothing",
            "fraction_of_set_standing_out": "the Xenium pLN lane's stage-1 rule (>=2 of a set "
                                            "enriched, larger fraction wins) moves the non-immune "
                                            "share to 47.8%"},
        "subpanel_rule_table_requirement": {k: list(v) for k, v in
                                            subpanel_requirements().items()},
        "primary_immune": {k: list(v) for k, v in PRIMARY_IMMUNE.items()},
        "supporting_immune": list(SUPPORTING_IMMUNE),
        "structural": list(STRUCTURAL), "structural_dropped": STRUCTURAL_DROPPED,
        "markers_absent_from_export": absent,
        "subpanel_members": members,
        "subpanel_naming_markers": naming,
        "nk_licence_from_rule_table": {"rule": SUPPORTING_ONLY_SUBPANEL["NK_ILC"],
                                       "all_of": nk_all, "any_of": nk_any, "none_of": nk_none},
        "nk_ilc_reachability": "class_NK_ILC counts the STAGE-2 route only -- a whole cluster "
                               "reaching NK through the subpanel argmax.  CD56 is the NK "
                               "subpanel's only member on this panel and lineages.SUPPORTING_ONLY "
                               "holds it back from naming a cell alone, so that route is narrow "
                               "and a zero on it is a property of this panel, NOT a measurement "
                               "that the node holds no NK cells.  NK IS NAMED AT STAGE 3 "
                               "INSTEAD, by negative gating -- see t_nk_resolution.",
        "t_nk_resolution": {
            "what": "the combined T_NK compartment (stage-2 class T_cell) split into T_CD4, "
                    "T_CD8, T_reg and NK_ILC by the presence/absence pattern of CD3e with CD8, "
                    "CD4, FOXP3 and CD56.  No channel on this panel is NK-specific and none has "
                    "to be: NK is CD3e-NEGATIVE with CD8 or CD56, and negative evidence is "
                    "evidence.",
            "negative_term": "a cluster, or a side of a within-cluster split, is CD3e-NEGATIVE "
                             "when it reads CD3e at AUC <= 0.5 against the rest of its own "
                             "section.  0.5 is the value the statistic takes for a population "
                             "indistinguishable from the rest of the slide, so this is a "
                             "comparison and not a cutoff; no intensity, level or quantile "
                             "enters and nothing is carried between sections.",
            "sub_cluster": "a split is attempted INSIDE every cluster of the compartment before "
                           "the cluster is labelled whole.  CD3e's best 2-group split of the "
                           "cluster's own ranks is scored by eta2 and compared against the same "
                           "calculation over every channel of the panel on the same cells -- any "
                           "variable splits in two somehow, so the panel is the null.  The split "
                           "is taken only when it beats that null, the CD3e-low side sits at or "
                           "below the rest of the section, and that side reads higher than the "
                           "CD3e-high side on CD8 or CD56 at the section's Bonferroni alpha.",
            "positive_terms": "argmax over the within-section ranks of CD8, CD4 and FOXP3, so "
                              "the three are mutually exclusive and SUM to the CD3e-positive "
                              "cells; CD4 % + CD8 % among them cannot pass 100 % arithmetically.",
            "unresolved": "a cluster that neither splits nor reads CD3e clearly one way stays "
                          "T_NK_unresolved.  That is a measurement of the cluster and it is "
                          "reported with the numbers behind it, never as a missing value.",
            "nk_never_zero_when_undecidable": "a donor's NK_ILC count is written EMPTY, with the "
                                              "reason, when no cluster of its compartment carries "
                                              "both a CD3e measurement and a CD8 or CD56 "
                                              "measurement.  An unavailable measurement is not "
                                              "zero.",
            "classes": list(T_NK_CLASSES),
            "positive_subsets": {lab: mk for lab, mk in T_NK_POSITIVE_SUBSETS},
            "nk_positive": list(NK_POSITIVE),
        },
        "partial_runs_are_scoped": "a run that does not cover every clustered donor and every "
                                   "stage writes its tables under _partial/<scope>/ and never "
                                   "over the cohort's own copies (D4, 2026-09-17)",
        "coherence_q_grid": list(COHERENCE_Q),
        "donors": status}, indent=1))

    pd.set_option("display.width", 250)
    if counts:
        print("\n" + pd.DataFrame(counts).to_string(index=False,
                                                    float_format=lambda v: f"{v:.2f}"))
    if sep_rows:
        print("\nSEPARATION -- can the clusters of each section be told apart on identity at all?")
        print(pd.DataFrame(sep_rows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    if coh_rows:
        keep = ["donor_id", "n_T_cells"] + [c for c in pd.DataFrame(coh_rows).columns
                                            if "plus_cd8" in c or c == "contradiction_cd4_cd8"]
        print("\nCOHERENCE (internal contradictions only -- no expectation is used):")
        print(pd.DataFrame(coh_rows)[keep].to_string(index=False,
                                                     float_format=lambda v: f"{v:.2f}"))
    if tnk_donor_rows:
        t = pd.DataFrame(tnk_donor_rows)
        keep = ["donor_id", "n_T_NK_compartment", "T_CD4", "T_CD8", "T_reg", "NK_ILC",
                T_NK_RESIDUAL, "subsets_sum_to_the_compartment",
                "pct_CD4_plus_CD8_of_CD3e_positive", "nk_ilc_state"]
        print("\nSTAGE 3 -- the T_NK compartment split on CD3e (negative evidence is evidence):")
        print(t[keep].to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    typed = [d for d, s in status.items() if s["state"] == "typed"]
    print(f"\n{len(typed)} typed, "
          f"{sum(1 for s in status.values() if s['state'] == 'unavailable')} unavailable, "
          f"{sum(1 for s in status.values() if s['state'] == 'withheld')} withheld")
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
