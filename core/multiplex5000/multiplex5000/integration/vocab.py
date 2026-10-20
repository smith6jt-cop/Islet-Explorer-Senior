"""
S2 — harmonising two label vocabularies and two molecular alphabets.

Both pipelines call broad classes, which makes them look comparable. They are not the same
partition of the tissue:

    PhenoCycler   Epithelial  Fibroblast  Muscle  Immune              Endocrine  Endothelial  Neural
    Xenium        Exocrine    Stromal ----------  Immune_T/B/Myeloid  Endocrine  Vascular     Neural

(PhenoCycler additionally carries a ``Neutrophil`` class while MPO is gated separately; see
the granulocyte note below for what happens when it is folded into ``Immune``.)

Four specific mismatches, each handled explicitly rather than averaged away:

*Immune.* PhenoCycler gates CD3e/CD20/CD163 but collapses all three into one ``Immune``
class. Xenium resolves T, B and myeloid as separate broad lineages. The split is recoverable
on the PhenoCycler side without re-running anything, because the per-marker ``_pos`` columns
survive in ``restore_gated_redsea/`` — ``export_pheno`` derives ``immune_subclass`` from
them, and this module maps that, not the broad label, onto the common vocabulary.

*Fibroblast vs Muscle.* PhenoCycler separates them by Vimentin/SMA argmax. Xenium folds
stellate cells and pericyte/smooth-muscle into one ``Stromal`` lineage. Worse, ``ACTA2`` is
not in the Xenium 5K reference at all — it arrives only via the custom-100 add-on, and
``12_pericyte_smooth_muscle`` contributes exactly one SMC-hinted gene (``MYLK``). So the
common vocabulary keeps ``Stromal`` as the comparable unit and marks the finer
Fibroblast/Muscle split PhenoCycler-only.

*Granulocytes.* This one tracks the core pipeline rather than asserting a fixed answer.
While a dedicated ``Neutrophil`` lineage exists (MPO gate), PhenoCycler resolves granulocytes
at broad level and Xenium does not — Xenium folds ``07_immune_granulocyte`` into ``Myeloid``,
reaching ``Granulocyte`` only through a fine label. When MPO is instead folded into ``Immune``
(neutrophils are immune cells), PhenoCycler stops resolving them too: ``immune_subclass``
gates CD3e/CD20/CD163 only, so an MPO-only cell has no subclass and lands in ``Myeloid`` —
exactly where Xenium puts it. The two sides then agree, at the cost of the distinction.

The PhenoCycler half is supplied by a declared :class:`PhenotypeSpace`; it is never imported
from core-pipeline Python configuration. A producer-label mapping that outlives its producer is
not harmless: it keeps advertising a capability that no current release has. ``Granulocyte``
never enters the comparable set because Xenium cannot call it at broad level — a lineage only
one side can see is not evidence of a difference between modalities.

Restoring the distinction later is additive on both sides: an ``MPO`` entry in
``export_pheno.IMMUNE_GATES`` for protein, the granulocyte identity-core score for RNA.

*Epithelial.* Xenium deliberately has no epithelial mask — ``16_epithelial_general`` is
``embedding_extra``, not ``broad_mask``, because KRT/EPCAM are shared by ductal, acinar and
endocrine cells and scoring them as one lineage biases the argmax. PhenoCycler's
``Epithelial``, meanwhile, is a *default sink*. Mapping it to Exocrine is defensible but the
asymmetry is real, which is why the contract carries ``evidence_positive``.

The mapping is *derived* from ``core/panels/data/panel_roles.csv`` and the
``identity_core/*.csv`` files rather than hardcoded — the upstream repo ships those two
artifacts specifically as a machine-readable contract for Python consumers, and warns that
unpinned masks drift silently. The taxonomy used to be a pinned submodule; in the monorepo it
is the ``core/panels`` subtree, so the pin is a directory rather than a gitlink.

The molecular crosswalk has its own trap. 46 of 55 PhenoCycler markers have a gene on the
deployed 5,081-gene panel, which sounds like plenty — but **INS, GCG, SST and Vimentin are
off-panel**, and those are precisely the markers that define PhenoCycler's Endocrine and
Fibroblast classes. So for those lineages there is no 1:1 protein<->gene pair and the
comparison has to route through Xenium's curated *surrogate* identity panels (PDX1, ISL1,
NEUROD1, ABCC8... for endocrine). That is orthogonal evidence rather than circular
confirmation, which is scientifically a feature — but it must be stated, not assumed away,
and anything that consumes ``protein_gene_pairs()`` must handle the reduced anchor set.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from ..config import PipelineConfig
from ..marker_taxonomy import PROCESS, TYPE
from . import tissues as _tissues
from .phenospace import PhenotypeSpace, load_phenotype_space

# --------------------------------------------------------------------------- #
# The common vocabulary
# --------------------------------------------------------------------------- #

#: Harmonised broad classes. Deliberately finer than PhenoCycler on the immune axis and
#: coarser on the stromal axis — each level is one that *both* modalities can actually
#: support, given `Stromal` cannot be split by Xenium and `Immune` cannot be left unsplit
#: without discarding Xenium's resolution.
COMMON_LINEAGES: tuple[str, ...] = (
    "Endocrine",
    "Exocrine",       # PhenoCycler Epithelial (Pan-CK) <-> Xenium Acinar/Ductal
    "Stromal",        # PhenoCycler Fibroblast + Muscle <-> Xenium Stellate + pericyte/SMC
    "Vascular",
    "T_NK",
    "B_Plasma",
    "Myeloid",
    "Granulocyte",
    "Neural",
    "Unassigned",
)

#: The former seven-class mapping is retained only as an explicit audit/test declaration.
#: Production callers receive the space from their configured identity release.
LEGACY_SPACE_PATH = Path(__file__).with_name("data") / "legacy_broad_space.json"


def load_legacy_space() -> PhenotypeSpace:
    """Load the retired broad map explicitly; never install it as mutable module state."""

    return load_phenotype_space(LEGACY_SPACE_PATH)

#: Xenium celltype_broad -> common. Mirrors `panel_roles.csv` `lineage_or_group`, which is
#: verified at load time by :func:`check_panel_roles_agreement`.
XENIUM_TO_COMMON: dict[str, str] = {
    "Endocrine": "Endocrine",
    "Exocrine": "Exocrine",
    "Stromal": "Stromal",
    "Vascular": "Vascular",
    "Immune_T": "T_NK",
    "Immune_B": "B_Plasma",
    "Myeloid": "Myeloid",
    "Neural": "Neural",
}

#: Xenium fine labels that should override the broad mapping. Granulocytes have no Xenium
#: broad class (folded into Myeloid), so the only route to a Granulocyte call is a fine
#: label or the granulocyte identity-core score.
XENIUM_FINE_TO_COMMON: dict[str, str] = {
    "Mast": "Granulocyte",
    "Neutrophil": "Granulocyte",
    "Eosinophil": "Granulocyte",
    "Basophil": "Granulocyte",
}

#: Which modality can resolve which common lineage, and why not when it cannot. Consumers
#: check this before reporting a lineage — a class only one side can see is not evidence of
#: a difference between modalities.
RESOLVABLE: dict[str, dict[str, str]] = {
    # Panel-only default: INS/GCG/SST are absent from the 5K panel, so Xenium infers
    # endocrine identity from surrogate cores (PDX1, ISL1, NEUROD1, ABCC8...). Post-Xenium IF
    # staining for INS/GCG makes it a direct protein measurement on the same cells — see
    # `resolvable_for`, which upgrades this when [integration] xenium_hormone_if is set.
    "Endocrine":   {"phenocycler": "yes", "xenium": "surrogate"},
    "Exocrine":    {"phenocycler": "default_sink", "xenium": "yes"},
    "Stromal":     {"phenocycler": "yes", "xenium": "yes"},
    "Vascular":    {"phenocycler": "yes", "xenium": "yes"},
    "T_NK":        {"phenocycler": "subclass", "xenium": "yes"},
    "B_Plasma":    {"phenocycler": "subclass", "xenium": "yes"},
    "Myeloid":     {"phenocycler": "subclass", "xenium": "yes"},
    # The default says no current broad release. ``resolvable_for(space=...)`` upgrades this
    # only when that release explicitly maps one of its labels to Granulocyte.
    "Granulocyte": {"phenocycler": "fine_only", "xenium": "fine_only"},
    "Neural":      {"phenocycler": "yes", "xenium": "yes"},
    "Unassigned":  {"phenocycler": "never", "xenium": "no"},
}

#: How a non-call row's two resolution columns read. Kept apart from `_RESOLVABLE_NOTES`
#: because these are states of a NON-call and those are states of a lineage call.
_NONCALL_NOTES = {
    "noncall": "this modality records this absence, with this label and a stated reason",
    "not_recorded": ("this modality writes no such label; the absence is either handled before "
                     "a class is written or does not arise in this method"),
}

_RESOLVABLE_NOTES = {
    "yes": "resolved by a positive broad-level call",
    "surrogate": "no 1:1 marker; resolved via curated surrogate identity panel",
    "subclass": "recovered from per-marker gates, not the broad label",
    "default_sink": "assigned by fallback when all positive gates fail (epi_default)",
    "fine_only": "no broad class; only reachable as a fine subtype",
    "never": "this modality has no Unassigned bucket by construction",
    "no": "not resolvable",
}

#: PhenoCycler-only distinctions that the common vocabulary deliberately collapses. Kept so
#: the loss is documented in the crosswalk CSV rather than being invisible.
PHENO_ONLY_SPLITS = {
    "Stromal": ["Fibroblast (Vimentin)", "Muscle (SMA)"],
}


# --------------------------------------------------------------------------- #
# Honest non-calls -- recorded, paired, and NOT lineages
# --------------------------------------------------------------------------- #
#
# Both producers write labels that mean "we looked and resolved nothing". They are not
# lineages and they must not become lineages here:
#
#   * `analysis/xenium/scripts/panc_xenium_release.py` declares them in a SEPARATE
#     ``"uncalled"`` array from ``"labels"``, each with a reason, and every such row carries
#     ``common_lineage = ""``. The release refuses to give them a lineage; so does this module.
#   * `scripts/senior/modality_concordance.py::comparable_lineages` keeps any lineage row that
#     names a label on both sides. A non-call added to COMMON_LINEAGES would therefore become a
#     "comparable lineage" and be compared as though it were a cell population -- the precise
#     failure this section exists to prevent.
#
# They must not silently become ``Unassigned`` either. ``Unassigned`` is THIS CROSSWALK failing
# to map a label, which is a defect in the mapping; a non-call is the producer's own recorded
# answer and is not a defect at all. Collapsing the two loses the difference between "we could
# not read your label" and "you told us there was nothing there", and loses the four distinct
# reasons the Xenium side takes the trouble to record.
#
# So: non-calls resolve to the EMPTY common lineage -- the release's own convention -- and the
# pairing is carried as its own block of the crosswalk CSV, queryable through
# :func:`noncall_for_label`.

#: The protein side's one non-call, written by `core/multiplex5000/scripts/typing/refine.py`
#: (``refine.UNRESOLVED``). It does NOT follow the Xenium ``Uncalled_*`` naming convention, so
#: it is matched by its exact spelling and nothing else.
PHENO_UNRESOLVED = "Unresolved (no lineage evidence)"


@dataclass(frozen=True)
class NonCall:
    """One recorded reason for having no identity, on one or both sides."""

    reason_id: str
    phenocycler: str   # "" where this modality records no such label
    xenium: str        # "" where this modality records no such label
    comparable: bool   # may a protein-vs-RNA comparison put these two side by side?
    note: str

    @property
    def paired(self) -> bool:
        return bool(self.phenocycler) and bool(self.xenium)


#: The non-call vocabulary, ordered as the crosswalk CSV writes it. The first four rows keep
#: the positions they have always had; the pLN rows were APPENDED on 2026-09-18 rather than
#: interleaved, so nothing that already read this block moved.
#:
#: ONE PAIR IS COMPARABLE and six are not, and the asymmetry is the finding, not a gap to be
#: filled: the Xenium side records WHY it did not call, in six ways across two tissues, and the
#: protein side records one. Pairing the six against the one protein label would manufacture an
#: agreement.
#:
#: TWO XENIUM PRODUCERS, NOT ONE. `analysis/xenium/scripts/panc_xenium_release.py` writes the
#: pancreas `Uncalled_*` family; `analysis/xenium/scripts/pln_immune_typing_v3.py` writes the
#: pLN one. Their reasons are NOT interchangeable -- the pancreas lane names a cluster against
#: identity cores, the pLN lane first has to decide immune-vs-structural at all -- so each label
#: carries its own row rather than being folded into a similarly-worded one.
#:
#: `Non_immune` IS NOT HERE, deliberately. The pLN release also writes it with
#: `common_lineage = ""`, but it is a measured negative -- stage 1 looked and answered "not
#: immune" -- and an answer must never be collapsed into an absence. An empty common lineage is
#: not by itself evidence of a non-call; the producer's stated reason is.
#:
#: `Uncalled_not_clustered` IS NOT HERE either. Measured 2026-09-18: its only producer is
#: `pln_immune_typing_v2.py`, which its own header declares SUPERSEDED by v3 and which
#: `scripts/retire/dispositions.py` records as a verified orphan (code tier 1). v3 replaced it
#: with `Uncalled_not_partitioned` below. A mapping that outlives its producer keeps advertising
#: a capability no current release has, so the dead spelling is recorded here in prose and not
#: in the table.
NON_CALLS: tuple[NonCall, ...] = (
    NonCall(
        reason_id="no_identity_evidence",
        phenocycler=PHENO_UNRESOLVED,
        xenium="Uncalled_no_naming_evidence",
        comparable=True,
        note=("both sides looked and found no evidence that names an identity. Protein: the "
              "refinement pass demoted the cell out of Exocrine, the complement class, because "
              "none of Exocrine's own markers reached the class bound among its concordant "
              "neighbours -- so there is no complement beneath it. RNA: the winning identity "
              "core has no naming gene standing out above the cluster's own background, and the "
              "call is never reassigned to the runner-up. Same statement, two alphabets"),
    ),
    NonCall(
        reason_id="no_enrichment",
        phenocycler="",
        xenium="Uncalled_no_enrichment",
        comparable=False,
        note=("RNA only: the winning core's score is not above the section mean at all. The "
              "protein classification has no analogue -- it scores a cell against reference "
              "levels drawn from the ordinary tissue, not against the other clusters of its own "
              "section -- so there is nothing on the protein side this could be compared with"),
    ),
    NonCall(
        reason_id="below_cluster_qc",
        phenocycler="",
        xenium="Uncalled_below_cluster_qc",
        comparable=False,
        note=("RNA only: too little was measured on the cell to partition on. The protein side "
              "records its own measurement QC through `qc_analysis_eligible` rather than as a "
              "label, so the two are not the same quantity and are not paired here"),
    ),
    NonCall(
        reason_id="no_nucleus",
        phenocycler="",
        xenium="Unavailable_no_nucleus",
        comparable=False,
        note=("NOT A CELL. `nucleus_count == 0`; excluded from every calculation, never a zero "
              "and never a member of any denominator. The protein side excludes such detections "
              "before a class is written rather than labelling them, so there is no protein "
              "label to pair -- which is why this row exists: a consumer that pairs it with "
              "anything has already broken the rule"),
    ),
    # ---- the pLN lane, `pln_immune_typing_v3.py`. Three labels that fell through as calls --
    # measured 2026-09-18: `xenium_to_common` returned "Unassigned" for all three, which says
    # this crosswalk could not read the label. It could; the producer had answered.
    NonCall(
        reason_id="ambiguous_at_stage1",
        phenocycler="",
        xenium="Uncalled_ambiguous",
        comparable=False,
        note=("RNA only, pLN: stage 1 asks immune-vs-structural on a within-section cluster and "
              "reached NEITHER side. `stage1_decide` records which of three ways it failed -- "
              "`no_set_reaches_two_standouts` (neither gene set reaches 2 standout genes), "
              "`equal_standout_fractions`, or `standout_count_and_fraction_disagree` -- and the "
              "release carries that in its own `reason` column, so this row is the label and not "
              "the whole reason. NOT the same absence as the pancreas lane's: there the identity "
              "cores compete, here it is not yet settled that the cell is immune at all. The "
              "protein side never asks this question -- a pancreas cell gets its context from the "
              "QuPath pixel classifier, not from a partition of its own section -- so there is "
              "nothing on the protein side to compare this with"),
    ),
    NonCall(
        reason_id="immune_without_lineage",
        phenocycler="",
        xenium="Immune_uncalled",
        comparable=False,
        note=("RNA only, pLN: A POSITIVE CALL THAT STOPS SHORT. Stage 1 called the cluster "
              "Immune and that call STANDS; stage 2's winning identity core then had zero "
              "standout naming genes, so the lineage is withheld and the cluster is never "
              "reassigned to the runner-up. So this is an absence at the lineage level over a "
              "measurement at the compartment level, and a consumer that counts it as "
              "'not immune' has inverted the producer. `scripts/senior/build_grid_composition.py` "
              "already reads it that way, into its `partial` bucket. It does not follow the "
              "`Uncalled_*` prefix; it is the `*_uncalled` spelling, in the same release. The "
              "protein immune pass reports its own undecidables per donor rather than on the "
              "cell, so there is no per-cell protein label to pair"),
    ),
    NonCall(
        reason_id="not_partitioned",
        phenocycler="",
        xenium="Uncalled_not_partitioned",
        comparable=False,
        note=("RNA only, pLN: the cell HAS a nucleus and still fell below the partition QC "
              "(`nCount_Xenium < 10` or `nFeature_Xenium < 5`), so it was never inside any "
              "cluster and inherited no cluster call -- carried as uncalled, never as 0. Reads "
              "close to `below_cluster_qc` and is a different measurement: that one is the "
              "pancreas lane's whole-cluster QC floor, this one is a per-cell floor in a "
              "different release, with its own thresholds. They are kept apart rather than "
              "merged, because one reason_id standing for two producers' floors would make the "
              "crosswalk unable to say which release an absence came from. Its v2 predecessor "
              "spelled this `Uncalled_not_clustered`; that producer is superseded and the "
              "spelling is deliberately absent from this table"),
    ),
)

_NONCALL_BY_LABEL: dict[str, NonCall] = {}
for _nc in NON_CALLS:
    for _label in (_nc.phenocycler, _nc.xenium):
        if _label:
            _NONCALL_BY_LABEL[_label] = _nc


def noncall_for_label(label: str) -> Optional[NonCall]:
    """The recorded non-call a producer label names, or ``None`` if it names a call.

    Exact match on either side's spelling. Deliberately NOT a prefix rule: a prefix rule here
    would claim authority over labels no release has declared, and the four reasons above are
    read off the release declarations rather than guessed from their shape.
    """
    return _NONCALL_BY_LABEL.get((label or "").strip())


def is_noncall(label: str) -> bool:
    """True when this producer label is a recorded absence rather than an identity."""
    return noncall_for_label(label) is not None


def _check_noncalls_are_well_formed() -> None:
    """Checked at import, so a malformed addition fails here and not in a released CSV."""
    ids = [nc.reason_id for nc in NON_CALLS]
    duplicated = sorted({i for i in ids if ids.count(i) > 1})
    if duplicated:
        raise AssertionError(f"NON_CALLS reason_id repeated: {duplicated}")
    for nc in NON_CALLS:
        if not nc.phenocycler and not nc.xenium:
            raise AssertionError(f"NON_CALLS {nc.reason_id!r} names no label on either side")
        if not nc.note:
            raise AssertionError(
                f"NON_CALLS {nc.reason_id!r} has no note. An absence without a recorded reason "
                "is indistinguishable from an oversight."
            )
        if nc.comparable and not nc.paired:
            raise AssertionError(
                f"NON_CALLS {nc.reason_id!r} is marked comparable but only one modality names "
                "it; a lineage only one side can see is not evidence of a difference"
            )
    labels = [lab for nc in NON_CALLS for lab in (nc.phenocycler, nc.xenium) if lab]
    repeated = sorted({lab for lab in labels if labels.count(lab) > 1})
    if repeated:
        raise AssertionError(f"one label carries two non-call reasons: {repeated}")
    clash = sorted(set(labels) & (set(COMMON_LINEAGES) | set(XENIUM_TO_COMMON)))
    if clash:
        raise AssertionError(f"non-call label collides with a lineage name: {clash}")


_check_noncalls_are_well_formed()


# --------------------------------------------------------------------------- #
# Is-immune, and where that identity holds
# --------------------------------------------------------------------------- #
#
# A lineage's name is not a machine-readable statement that it is immune. Nothing downstream
# may infer "T_NK is immune" from the string, so the answer is recorded here and written into
# the crosswalk CSV, and a consumer that cannot find it is right to withhold rather than
# guess. `scripts/senior/build_distance_profiles.py` is the first such consumer.
#
# TISSUE SCOPE. Maintainer ruling, 2026-09-17: "All immune cell calls from pancreas can be
# applied to pLN." The immune lineage vocabulary established on pancreas is the same
# vocabulary in the pancreatic lymph node, and the typing method that produced the pancreas
# calls is licensed to run on pLN cells. That is a ruling about IDENTITY. It is not a ruling
# about anatomy or coordinates: the pancreas and pLN sections mounted on one slide are cut
# from different blocks, they share a pixel grid and nothing anatomically, and no distance,
# adjacency, neighbourhood or KDE may cross between them. It also does not license copying
# label rows across -- the two ROIs are segmented separately, so their cell_id namespaces are
# disjoint and the typer must be RUN on the pLN cells.

#: Common lineages that are immune.
IMMUNE_LINEAGES: tuple[str, ...] = ("T_NK", "B_Plasma", "Myeloid", "Granulocyte")

#: Common lineages that are positively not immune -- a recorded "no", not a default.
NON_IMMUNE_LINEAGES: tuple[str, ...] = (
    "Endocrine", "Exocrine", "Stromal", "Vascular", "Neural",
)

#: Common lineages for which is-immune has no answer. ``Unassigned`` is the absence of a
#: lineage call, so it is neither immune nor non-immune; it must never be counted as either.
IMMUNE_UNRESOLVED_LINEAGES: tuple[str, ...] = ("Unassigned",)

#: The tissues in which :data:`IMMUNE_LINEAGES` carry an immune identity, per the ruling.
#: Names match ``PipelineConfig.tissue_list`` / ``tissue_for_roi``.
IMMUNE_TISSUE_SCOPE: tuple[str, ...] = (_tissues.PANCREAS, _tissues.LYMPH_NODE)

_IMMUNE_NOTES = {
    "yes": ("immune lineage; immune identity holds in {scope} "
            "(maintainer ruling 2026-09-17: all immune cell calls from pancreas can be "
            "applied to pLN -- identity only, no spatial relation crosses the tissues)"),
    "no": "not an immune lineage",
    "unresolvable": ("no lineage call; is-immune is not answerable and this bucket is "
                     "neither immune nor non-immune"),
    "not_applicable": "",
}

#: What a marker row carries in the is-immune column. Marker rows are not lineages, and an
#: empty string there would read as "not immune".
IMMUNE_NOT_APPLICABLE = "not_applicable"


def immune_status(lineage: str) -> str:
    """``yes`` / ``no`` / ``unresolvable`` for a common lineage.

    Raises on an unrecognised lineage rather than returning ``no``: a lineage nobody has
    classified is not evidence that it is not immune.
    """
    if lineage in IMMUNE_LINEAGES:
        return "yes"
    if lineage in NON_IMMUNE_LINEAGES:
        return "no"
    if lineage in IMMUNE_UNRESOLVED_LINEAGES:
        return "unresolvable"
    raise KeyError(
        f"{lineage!r} has no recorded is-immune classification. Add it to IMMUNE_LINEAGES, "
        "NON_IMMUNE_LINEAGES or IMMUNE_UNRESOLVED_LINEAGES in phenocycler.integration.vocab; "
        "do not let it default."
    )


def immune_note(status: str) -> str:
    """The reason string written beside the is-immune value, the way every other column
    in the crosswalk carries a note."""
    return _IMMUNE_NOTES.get(status, "").format(scope=" + ".join(IMMUNE_TISSUE_SCOPE))


def immune_tissue_scope(status: str) -> str:
    """Tissues in which this is-immune value applies, ``|``-joined, or '' when it does not
    describe an immune lineage."""
    return "|".join(IMMUNE_TISSUE_SCOPE) if status == "yes" else ""


def _check_immune_classification_is_total() -> None:
    """Every common lineage is classified exactly once, checked at import.

    A lineage added to COMMON_LINEAGES without a decision here fails loudly instead of
    silently arriving in the crosswalk as non-immune.
    """
    buckets = (IMMUNE_LINEAGES, NON_IMMUNE_LINEAGES, IMMUNE_UNRESOLVED_LINEAGES)
    classified = [lin for bucket in buckets for lin in bucket]
    duplicated = sorted({lin for lin in classified if classified.count(lin) > 1})
    if duplicated:
        raise AssertionError(f"lineages classified more than once as immune/not: {duplicated}")
    missing = [lin for lin in COMMON_LINEAGES if lin not in set(classified)]
    unknown = sorted(set(classified) - set(COMMON_LINEAGES))
    if missing or unknown:
        raise AssertionError(
            f"is-immune classification does not cover COMMON_LINEAGES: missing {missing}, "
            f"not a common lineage {unknown}"
        )


_check_immune_classification_is_total()

# --------------------------------------------------------------------------- #
# Protein <-> gene crosswalk
# --------------------------------------------------------------------------- #

#: PhenoCycler marker -> candidate gene symbol(s), best first. Candidates are *proposals*;
#: :func:`protein_gene_pairs` keeps only those actually present in the run's panel, so a
#: panel change silently shrinks the anchor set instead of producing phantom features.
PROTEIN_GENE_CANDIDATES: dict[str, list[str]] = {
    # epithelial / exocrine
    "Pan_Cytokeratin": ["KRT19", "KRT8", "KRT18", "KRT7"],
    "EpCAM": ["EPCAM"],
    "E_cadherin": ["CDH1"],
    "TP63": ["TP63"],
    "CD66": ["CEACAM1", "CEACAM5", "CEACAM6"],
    # mesenchyme / muscle
    "Vimentin": ["VIM"],                    # off-panel -> surrogate
    "Collagen_IV": ["COL4A1", "COL4A2"],
    "SMA": ["ACTA2"],                       # custom-100 only
    # neural
    "B3TUBB": ["TUBB3"],
    "CD56": ["NCAM1"],
    # vascular
    "CD31": ["PECAM1"],
    "CD34": ["CD34"],
    "Podoplanin": ["PDPN"],
    "Caveolin": ["CAV1"],                   # off-panel
    # endocrine
    "INS": ["INS"],                         # off-panel -> surrogate
    "GCG": ["GCG"],                         # off-panel -> surrogate
    "SST": ["SST"],                         # off-panel -> surrogate
    "CD99": ["CD99"],
    # immune
    "CD3e": ["CD3E", "CD3D", "CD3G"],
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["CD4"],
    "FOXP3": ["FOXP3"],
    "CD20": ["MS4A1"],
    "CD79a": ["CD79A"],
    "CD68": ["CD68"],
    "CD163": ["CD163"],
    "Iba1": ["AIF1"],                       # off-panel
    "CD206": ["MRC1"],
    "M2Gal3": ["LGALS3"],                   # off-panel
    "CD11c": ["ITGAX"],
    "CD209": ["CD209"],
    "CD11b": ["ITGAM"],
    "MPO": ["MPO"],
    "b_Catenin1": ["CTNNB1"],
    # process / state
    "Ki67": ["MKI67"],
    "PCNA": ["PCNA"],
    "Granzyme_B": ["GZMB"],
    "CD107a": ["LAMP1"],
    "PD_1": ["PDCD1"],
    "PD_L1": ["CD274"],
    "LAG3": ["LAG3"],
    "TOX": ["TOX"],
    "VISTA": ["VSIR"],
    "ICOS": ["ICOS"],
    "TCF_1": ["TCF7"],
    "CD39": ["ENTPD1"],
    "HLA_DR": ["HLA-DRA", "HLA-DRB1"],      # off-panel
    "HLA_A": ["HLA-A"],                     # off-panel
    "iNOS": ["NOS2"],
    "IDO1": ["IDO1"],
    "SOX2": ["SOX2"],
    "Bcl_2": ["BCL2"],
    "CD44": ["CD44"],
    "CD57": ["B3GAT1"],
    "CD38": ["CD38"],
}

#: Where a protein has no gene on the panel, the comparison routes through the Xenium
#: identity-core panel for the lineage that protein defines. `identity_core/<stem>.csv`
#: supplies the actual gene list; this only names which panel to use.
SURROGATE_PANEL: dict[str, str] = {
    "INS": "01_pancreas_endocrine",
    "GCG": "01_pancreas_endocrine",
    "SST": "01_pancreas_endocrine",
    "Vimentin": "13_fibroblast_stellate",
    "Caveolin": "11_endothelial_vascular",
    "Iba1": "06_immune_myeloid",
    "M2Gal3": "06_immune_myeloid",
    "HLA_DR": "06_immune_myeloid",
    "HLA_A": "",                            # ubiquitous; no lineage surrogate is meaningful
}

#: Markers that gate PhenoCycler's broad lineage call. Anchors drawn from this set carry the
#: most weight for cross-modal identity comparison.
GATE_MARKERS = ("Pan_Cytokeratin", "Vimentin", "SMA", "CD31", "INS", "GCG", "SST",
                "CD3e", "CD20", "CD163", "CD99", "B3TUBB", "MPO")


class VocabError(RuntimeError):
    """The panel taxonomy could not be loaded or is internally inconsistent."""


@dataclass
class MarkerPair:
    """One resolved protein <-> gene(s) correspondence."""

    protein: str
    genes: list[str]                  # present on the panel; empty when off-panel
    on_panel: bool
    is_gate: bool
    role: str                         # 'TYPE' | 'PROCESS' | 'OTHER'
    surrogate_panel: str = ""         # identity_core stem used when off-panel
    surrogate_genes: list[str] = field(default_factory=list)

    @property
    def usable_anchor(self) -> bool:
        """Whether this pair can anchor a shared feature space.

        Requires a real gene on the panel — a surrogate panel supports *lineage-level*
        comparison but is not a per-cell feature correspondence, so it must not be counted
        toward the pseudo-cell linking anchor budget.
        """
        return self.on_panel and bool(self.genes)

    def anchor_in(self, mode: str) -> bool:
        """Whether this pair is a per-cell anchor under a given acquisition ``mode``.

        Protein-or-gene substitution is licensed by the slide, not by convenience. On one
        physical section the protein was measured on *this very cell*, so a marker with no
        on-panel gene is still a real per-cell feature — it simply arrives from the protein
        side. Across serial sections the identical substitution is invalid, because the
        protein was measured on a different cell; there it stays a lineage-level surrogate.

        This is what lets the same-slide pair carry identity the RNA panel cannot express:
        INS/GCG/SST are absent from both Xenium panels but are PhenoCycler gates, and
        GHRL/AMY1A are on the +100 T1D add-on but not the stock Prime 5k.
        """
        if self.usable_anchor:
            return True
        return mode == "same_slide" and bool(self.protein)

    def anchor_source(self, mode: str) -> str:
        """Which modality supplies this anchor: 'gene', 'protein', or '' if neither."""
        if self.usable_anchor:
            return "gene"
        if mode == "same_slide" and self.protein:
            return "protein"
        return ""


@dataclass
class PanelTaxonomy:
    """The XeniumPanelExplorer identity layer, loaded for one tissue."""

    tissue: str
    root: Path
    #: identity_core stem -> Xenium broad lineage (from panel_roles.csv broad_mask rows)
    panel_to_lineage: dict[str, str] = field(default_factory=dict)
    #: identity_core stem -> [gene]
    core_genes: dict[str, list[str]] = field(default_factory=dict)
    #: identity_core stem -> {gene: subtype_hint}
    core_hints: dict[str, dict[str, str]] = field(default_factory=dict)
    #: stems flagged embedding_extra rather than broad_mask (Epithelial)
    embedding_extra: list[str] = field(default_factory=list)

    def lineage_genes(self, lineage: str) -> list[str]:
        """Union of identity-core genes mapped to a Xenium broad lineage."""
        out: list[str] = []
        for stem, lin in self.panel_to_lineage.items():
            if lin == lineage:
                out.extend(self.core_genes.get(stem, []))
        return sorted(set(out))

    def genes_by_hint(self, stem: str, hint: str) -> list[str]:
        """Genes in one identity core carrying a given ``subtype_hint``.

        This is how a Granulocyte score is recovered from ``07_immune_granulocyte``, whose
        neutrophil-hinted genes are the only route to a class Xenium has no broad label for.
        """
        return sorted(g for g, h in self.core_hints.get(stem, {}).items() if h == hint)


def _read_csv_rows(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def load_panel_taxonomy(cfg: PipelineConfig, tissue: Optional[str] = None) -> PanelTaxonomy:
    """Load ``panel_roles.csv`` + the ``identity_core`` CSVs from the pinned submodule.

    ``tissue`` selects which cores apply. XeniumPanelExplorer ships ``pancreas`` and
    ``thymus``; a tissue with no directory of its own (pancreatic lymph node) resolves to the
    **shared** cores — T cell, B/plasma, NK/ILC, myeloid, granulocyte, endothelial, pericyte,
    neural, epithelial — and drops the three pancreas-specific ones. That is the correct
    answer for a lymph node rather than a degraded fallback: scoring lymph-node cells against
    ``01_pancreas_endocrine`` would put weight on an endocrine axis that cannot exist there,
    and the broad call is an argmax over those axes.

    Raises with an actionable message when the taxonomy is missing — this is the single most
    likely first-run failure, and a bare 'FileNotFoundError: panel_roles.csv' would not tell
    anyone where it is meant to come from.
    """
    from .tissues import applies as tissue_applies
    from .tissues import panel_tissue_dir

    tissue = tissue or cfg.tissue
    root = Path(cfg.panel_explorer)
    roles_csv = root / "data" / "panel_roles.csv"
    if not roles_csv.exists():
        raise VocabError(
            f"Panel taxonomy not found at {roles_csv}.\n"
            f"It ships as the core/panels subtree; point [integration] panel_explorer\n"
            f"at a checkout if you have moved it."
        )

    tax = PanelTaxonomy(tissue=tissue, root=root)
    # Which row set in panel_roles.csv describes this tissue's cores. A tissue with no
    # directory borrows the pancreas rows and then filters out the pancreas-specific stems,
    # which is how the shared pool is reached without duplicating the role table upstream.
    panel_dir = panel_tissue_dir(tissue)
    roles_tissue = panel_dir or "pancreas"

    for row in _read_csv_rows(roles_csv):
        if row.get("tissue") != roles_tissue:
            continue
        stem = (row.get("subpanel") or "").strip()
        use = (row.get("phenotyping_use") or "").strip()
        lineage = (row.get("lineage_or_group") or "").strip()
        if not stem or not tissue_applies(tissue, stem):
            continue
        if use == "broad_mask":
            tax.panel_to_lineage[stem] = lineage
        elif use == "embedding_extra":
            tax.embedding_extra.append(stem)

    # identity_core lives per-tissue for tissue-specific panels and in a shared pool for the
    # rest; panel_roles.phenotyping_source records which, but checking both is simpler and
    # equivalent (the tissue copy wins).
    search_dirs = [
        *( [root / "data" / "tissues" / panel_dir / "identity_core"] if panel_dir else [] ),
        root / "data" / "identity_core",
    ]
    for stem in list(tax.panel_to_lineage) + tax.embedding_extra:
        for d in search_dirs:
            path = d / f"{stem}.csv"
            if not path.exists():
                continue
            rows = _read_csv_rows(path)
            genes, hints = [], {}
            for r in rows:
                g = (r.get("gene") or "").strip()
                if not g:
                    continue
                genes.append(g)
                hints[g] = (r.get("subtype_hint") or "").strip()
            tax.core_genes[stem] = genes
            tax.core_hints[stem] = hints
            break

    missing = [s for s in tax.panel_to_lineage if s not in tax.core_genes]
    if missing:
        raise VocabError(
            f"identity_core CSVs missing for {missing} under {root}/data — the pinned "
            f"submodule looks incomplete"
        )
    return tax


def check_panel_roles_agreement(tax: PanelTaxonomy) -> list[str]:
    """Verify every upstream lineage has a place in the common vocabulary.

    A drift guard, not a formality: XeniumPanelExplorer is a separate repo on its own release
    cycle, and a renamed or added lineage would otherwise land silently in `Unassigned`.
    """
    problems = []
    for stem, lineage in sorted(tax.panel_to_lineage.items()):
        if lineage not in XENIUM_TO_COMMON:
            problems.append(
                f"panel_roles.csv maps {stem} -> lineage {lineage!r}, which XENIUM_TO_COMMON "
                f"does not cover (known: {sorted(XENIUM_TO_COMMON)})"
            )
    return problems


def _declared_pheno_label(
    lineage_native: str,
    immune_subclass: str = "",
) -> str:
    """Build the one producer label the declaration must resolve.

    The historical broad lane stored the immune branch in a second column.  Its JSON space
    declares those combinations as ``Immune/T`` etc.  A non-empty, undeclared subclass is a
    hard error instead of falling back to the broad label and erasing producer information.
    """

    native = (lineage_native or "").strip()
    subclass = (immune_subclass or "").strip()
    if native == "Immune" and subclass:
        return f"{native}/{subclass}"
    return native


def pheno_to_common(
    space: PhenotypeSpace,
    lineage_native: str,
    immune_subclass: str = "",
) -> str:
    """Declared PhenoCycler producer label -> common lineage.

    ``space`` is positional and has no default deliberately.  An undeclared producer label
    raises in :meth:`PhenotypeSpace.common_for`; it never becomes a plausible-looking
    ``Unassigned`` call through a stale module-global dictionary.

    A RECORDED NON-CALL RETURNS ``""``, checked before the space is consulted.  That is the
    Xenium releases' own convention for these rows and it is the honest answer: the producer
    resolved no identity, so there is no lineage, and there is also no mapping failure to
    report.  ``noncall_for_label()`` answers which absence it was.  The non-call is NOT
    declared in the label space, so leaving it to ``common_for`` would raise ``undeclared`` --
    which would be the crosswalk calling the producer's honest answer a defect.
    """

    if is_noncall(lineage_native):
        return ""
    common, resolution = space.common_for(
        _declared_pheno_label(lineage_native, immune_subclass)
    )
    return common if resolution == "mapped" else "Unassigned"


def xenium_to_common(celltype_broad: str, celltype_fine: str = "") -> str:
    """Xenium broad (+ fine) -> common lineage. Fine wins where it is more specific.

    A recorded non-call returns ``""`` rather than ``Unassigned``, on either column.  The
    release writes ``common_lineage = ""`` on exactly these rows; ``Unassigned`` would say
    instead that this crosswalk could not read the label, and would put an absence the
    producer measured into the same bucket as a mapping defect.
    """
    fine = (celltype_fine or "").strip()
    broad = (celltype_broad or "").strip()
    if is_noncall(broad) or is_noncall(fine):
        return ""
    if fine in XENIUM_FINE_TO_COMMON:
        return XENIUM_FINE_TO_COMMON[fine]
    return XENIUM_TO_COMMON.get(broad, "Unassigned")


def map_pheno_series(
    space: PhenotypeSpace,
    lineage_native: pd.Series,
    immune_subclass: Optional[pd.Series] = None,
) -> pd.Series:
    """Vectorised :func:`pheno_to_common` over a column pair."""
    if immune_subclass is None:
        immune_subclass = pd.Series([""] * len(lineage_native), index=lineage_native.index)
    return pd.Series(
        [pheno_to_common(space, a, b) for a, b in zip(lineage_native.fillna(""),
                                                      immune_subclass.fillna(""))],
        index=lineage_native.index, dtype="object",
    )


def map_xenium_series(celltype_broad: pd.Series, celltype_fine: Optional[pd.Series] = None
                      ) -> pd.Series:
    """Vectorised :func:`xenium_to_common` over a column pair."""
    if celltype_fine is None:
        celltype_fine = pd.Series([""] * len(celltype_broad), index=celltype_broad.index)
    return pd.Series(
        [xenium_to_common(a, b) for a, b in zip(celltype_broad.fillna(""),
                                                celltype_fine.fillna(""))],
        index=celltype_broad.index, dtype="object",
    )


def protein_gene_pairs(
    panel_genes: Iterable[str],
    tax: Optional[PanelTaxonomy] = None,
    markers: Optional[Iterable[str]] = None,
) -> list[MarkerPair]:
    """Resolve the protein <-> gene crosswalk against a run's actual panel.

    Parameters
    ----------
    panel_genes
        The gene symbols actually present in the Xenium run (``adata.var_names``). Passing
        the *nominal* panel here rather than the run's own genes is the mistake this
        signature exists to prevent: the deployed panel is the 5K base plus a custom add-on,
        and several textbook markers live only in the add-on.
    tax
        Panel taxonomy, used to attach surrogate gene lists for off-panel proteins.
    markers
        Restrict to these PhenoCycler markers; default is every marker with a candidate.
    """
    present = {str(g).strip() for g in panel_genes}
    wanted = list(markers) if markers is not None else list(PROTEIN_GENE_CANDIDATES)
    type_set, process_set = set(TYPE), set(PROCESS)

    out: list[MarkerPair] = []
    for protein in wanted:
        candidates = PROTEIN_GENE_CANDIDATES.get(protein, [])
        hits = [g for g in candidates if g in present]
        role = "TYPE" if protein in type_set else ("PROCESS" if protein in process_set
                                                   else "OTHER")
        stem = SURROGATE_PANEL.get(protein, "") if not hits else ""
        surrogate_genes: list[str] = []
        if stem and tax is not None:
            surrogate_genes = [g for g in tax.core_genes.get(stem, []) if g in present]
        out.append(MarkerPair(
            protein=protein,
            genes=hits,
            on_panel=bool(hits),
            is_gate=protein in GATE_MARKERS,
            role=role,
            surrogate_panel=stem,
            surrogate_genes=surrogate_genes,
        ))
    return out


def usable_anchors(pairs: Iterable[MarkerPair], *, type_only: bool = False,
                   mode: str = "sequential") -> list[MarkerPair]:
    """Pairs that can anchor a shared feature space.

    ``type_only`` keeps identity markers and drops state/process ones — appropriate when the
    anchors must discriminate cell *type* rather than activation state.

    ``mode`` is the pair's acquisition design (see :mod:`phenocycler.integration.modes`).
    Under ``same_slide`` a protein-only marker counts as a per-cell anchor; under
    ``sequential`` it does not. The default stays ``sequential`` so existing callers keep
    their current, stricter behaviour until they opt in.
    """
    out = [p for p in pairs if p.anchor_in(mode)]
    if type_only:
        out = [p for p in out if p.role == "TYPE"]
    return out


def resolvable_for(
    cfg: Optional[PipelineConfig] = None,
    *,
    space: Optional[PhenotypeSpace] = None,
) -> dict[str, dict[str, str]]:
    """`RESOLVABLE`, adjusted for what this run's Xenium sections actually carry.

    Post-Xenium immunofluorescence for INS/GCG turns Xenium's endocrine call from a surrogate
    inference into a direct protein measurement *on the same cells as the RNA* — the single
    biggest gap in the panel crosswalk, since the hormones are the markers that define the
    compartment a T1D study is about. Reported rather than assumed: with the flag off, the
    crosswalk keeps saying `surrogate`, which is the truthful answer for panel-only data.
    """
    out = {k: dict(v) for k, v in RESOLVABLE.items()}
    if space is not None:
        mapped = {
            space.common_for(label)[0]
            for label in space.rules
            if space.common_for(label)[1] == "mapped"
        }
        for lineage in COMMON_LINEAGES:
            if lineage == "Unassigned":
                continue
            if lineage not in mapped:
                out[lineage]["phenocycler"] = "no"
        if "Granulocyte" in mapped:
            out["Granulocyte"]["phenocycler"] = "yes"
    if cfg is not None and getattr(cfg, "xenium_hormone_if", False):
        out["Endocrine"]["xenium"] = "yes"
    return out


def crosswalk_frame(
    space: PhenotypeSpace,
    pairs: Iterable[MarkerPair],
    tax: Optional[PanelTaxonomy] = None,
) -> pd.DataFrame:
    """The full crosswalk as a table — this is what gets written to disk and reviewed.

    Three blocks: one row per common lineage (how each modality resolves it), one row per
    recorded non-call (an absence each producer measured), then one row per PhenoCycler marker
    (its gene, or the surrogate panel standing in for it).

    Lineage rows also carry ``is_immune`` (``yes``/``no``/``unresolvable``),
    ``immune_tissue_scope`` (the tissues that value holds in) and ``immune_note`` (why), so a
    consumer asking "is T_NK an immune lineage in the pancreatic lymph node" reads the answer
    instead of inferring it from the label. Marker and non-call rows carry ``not_applicable``:
    neither is a lineage, and an empty string there would read as "not immune".

    ``kind == "noncall"`` rows sit OUTSIDE the lineage block on purpose. Both current readers
    of this CSV -- ``build_distance_profiles.load_immune_vocabulary`` and
    ``modality_concordance.comparable_lineages`` -- select ``kind == "lineage"``, so the
    non-calls are visible to anyone reading the file and can never be folded into a lineage
    comparison by a consumer that was written before they existed.
    """
    rows: list[dict] = []
    resolutions = resolvable_for(space=space)
    pheno_by_common: dict[str, list[str]] = {lineage: [] for lineage in COMMON_LINEAGES}
    for label in space.rules:
        common, resolution = space.common_for(label)
        if resolution == "mapped" and common in pheno_by_common:
            pheno_by_common[common].append(label)

    for lineage in COMMON_LINEAGES:
        res = resolutions.get(lineage, {})
        p, x = res.get("phenocycler", "no"), res.get("xenium", "no")
        xen_native = sorted({k for k, v in XENIUM_TO_COMMON.items() if v == lineage})
        pheno_native = sorted(pheno_by_common[lineage])
        immune = immune_status(lineage)
        rows.append({
            "kind": "lineage",
            "common": lineage,
            "phenocycler": "|".join(pheno_native),
            "xenium": "|".join(xen_native),
            "pheno_resolution": p,
            "xenium_resolution": x,
            "pheno_note": _RESOLVABLE_NOTES.get(p, ""),
            "xenium_note": _RESOLVABLE_NOTES.get(x, ""),
            "collapsed": "|".join(PHENO_ONLY_SPLITS.get(lineage, [])),
            "identity_core_genes": (len(tax.lineage_genes(
                next((k for k, v in XENIUM_TO_COMMON.items() if v == lineage), "")))
                if tax is not None else ""),
            "is_immune": immune,
            "immune_tissue_scope": immune_tissue_scope(immune),
            "immune_note": immune_note(immune),
        })

    for nc in NON_CALLS:
        rows.append({
            "kind": "noncall",
            # NOT a lineage, so it has no common name. The pairing lives in the two label
            # columns and in `collapsed`, which names the reason both sides share.
            "common": "",
            "phenocycler": nc.phenocycler,
            "xenium": nc.xenium,
            "pheno_resolution": "noncall" if nc.phenocycler else "not_recorded",
            "xenium_resolution": "noncall" if nc.xenium else "not_recorded",
            "pheno_note": _NONCALL_NOTES["noncall" if nc.phenocycler else "not_recorded"],
            "xenium_note": _NONCALL_NOTES["noncall" if nc.xenium else "not_recorded"],
            "collapsed": nc.reason_id,
            "identity_core_genes": "",
            # A non-call is not a lineage. `unresolvable` is reserved for the `Unassigned`
            # LINEAGE, so using it here would put a row that is not a lineage into a column
            # readers scan for lineages.
            "is_immune": IMMUNE_NOT_APPLICABLE,
            "immune_tissue_scope": "",
            "immune_note": ("comparable with the other side's label on this row"
                            if nc.comparable
                            else "one modality only -- NOT comparable across modalities")
                           + ". " + nc.note,
        })

    for mp in pairs:
        rows.append({
            "kind": "marker",
            "common": "",
            "phenocycler": mp.protein,
            "xenium": "|".join(mp.genes),
            "pheno_resolution": mp.role,
            "xenium_resolution": "on_panel" if mp.on_panel else "off_panel",
            "pheno_note": "lineage gate" if mp.is_gate else "",
            "xenium_note": (f"surrogate: {mp.surrogate_panel} "
                            f"({len(mp.surrogate_genes)} genes)" if mp.surrogate_panel else ""),
            "collapsed": "",
            "identity_core_genes": len(mp.surrogate_genes) if mp.surrogate_panel else "",
            # A marker is not a lineage. An empty string here would read as "not immune".
            "is_immune": IMMUNE_NOT_APPLICABLE,
            "immune_tissue_scope": "",
            "immune_note": "",
        })

    return pd.DataFrame(rows)


def coverage_report(
    space: PhenotypeSpace,
    pairs: Iterable[MarkerPair],
    tax: Optional[PanelTaxonomy] = None,
) -> str:
    """Human-readable summary — printed by the stage and pasted into notebooks."""
    pairs = list(pairs)
    on = [p for p in pairs if p.on_panel]
    off = [p for p in pairs if not p.on_panel]
    gates_on = [p for p in on if p.is_gate]
    gates_off = [p for p in off if p.is_gate]
    anchors = usable_anchors(pairs, type_only=True)
    resolutions = resolvable_for(space=space)

    lines = [
        "PhenoCycler <-> Xenium vocabulary crosswalk",
        "=" * 60,
        f"common lineages      : {len(COMMON_LINEAGES) - 1} (+Unassigned)",
        f"recorded non-calls   : {len(NON_CALLS)} (not lineages; see below)",
        f"markers with a gene  : {len(on)}/{len(pairs)}",
        f"  lineage gates      : {len(gates_on)}/{len(GATE_MARKERS)}",
        f"identity anchors     : {len(anchors)} (TYPE markers with a panel gene)",
        "",
    ]
    if gates_off:
        lines.append("OFF-PANEL LINEAGE GATES — these route through surrogate panels, so the")
        lines.append("corresponding lineage comparison is orthogonal evidence, not a 1:1 check:")
        for p in gates_off:
            sg = f"{len(p.surrogate_genes)} surrogate genes" if p.surrogate_genes else "NO surrogate"
            lines.append(f"    {p.protein:<18} -> {p.surrogate_panel or '(none)'}  [{sg}]")
        lines.append("")
    other_off = [p for p in off if not p.is_gate]
    if other_off:
        lines.append(f"other off-panel markers ({len(other_off)}): "
                     f"{', '.join(p.protein for p in other_off)}")
        lines.append("")

    lines.append("lineage resolution by modality:")
    lines.append(f"  {'common':<13} {'phenocycler':<15} {'xenium':<15} {'immune':<13} note")
    for lin in COMMON_LINEAGES:
        r = resolutions.get(lin, {})
        p, x = r.get("phenocycler", "no"), r.get("xenium", "no")
        note = ""
        if lin in PHENO_ONLY_SPLITS:
            note = "collapses " + " + ".join(PHENO_ONLY_SPLITS[lin])
        elif x == "fine_only":
            note = "no Xenium broad class"
        elif p == "default_sink":
            note = "PhenoCycler fallback bucket — compare evidence-positive only"
        lines.append(f"  {lin:<13} {p:<15} {x:<15} {immune_status(lin):<13} {note}")

    lines.append("")
    lines.append("recorded non-calls (NOT lineages -- an absence each producer measured):")
    lines.append(f"  {'reason':<22} {'phenocycler':<34} {'xenium':<28} comparable")
    for nc in NON_CALLS:
        lines.append(f"  {nc.reason_id:<22} {nc.phenocycler or '(none)':<34} "
                     f"{nc.xenium or '(none)':<28} {'yes' if nc.comparable else 'no'}")
    paired = [nc.reason_id for nc in NON_CALLS if nc.comparable]
    lines.append(f"  {len(paired)} of {len(NON_CALLS)} pair across modalities: "
                 f"{', '.join(paired) or '(none)'}. The rest are one-sided and a comparison "
                 "must not read a difference from them.")

    if tax is not None:
        probs = check_panel_roles_agreement(tax)
        if probs:
            lines.append("")
            lines.append("DRIFT WARNINGS (XeniumPanelExplorer changed under us):")
            lines.extend(f"    {p}" for p in probs)
    return "\n".join(lines)


def run_vocab(cfg: PipelineConfig, panel_genes: Optional[Iterable[str]] = None,
              *, tissue: Optional[str] = None, write: bool = True) -> pd.DataFrame:
    """Stage entry point: build the crosswalk, print coverage, write the CSV.

    ``panel_genes`` defaults to the union of the tissue's identity cores, which is enough to
    build the lineage half of the crosswalk without any Xenium run mounted. Pass the run's
    real ``var_names`` for an accurate marker-level report.

    The crosswalk is tissue-specific because the identity cores are: a lymph node is scored
    against the shared immune/vascular/stromal cores only, so its lineage rows carry
    different core-gene counts than the pancreas rows. Each tissue therefore gets its own
    ``vocab_crosswalk_{tissue}.csv``, and the canonical ``vocab_crosswalk.csv`` is the
    concatenation of them, tagged by tissue.

    ``is_immune`` is not tissue-derived: it is the recorded classification from
    :data:`IMMUNE_LINEAGES` and friends, and ``immune_tissue_scope`` names the tissues it
    holds in, so the same answer is readable from any tissue's rows.
    """
    from .export_pheno import require_current_pheno

    release = require_current_pheno("integration phenotype vocabulary", cfg=cfg)
    space = release.space
    tissue = tissue or cfg.tissue
    tax = load_panel_taxonomy(cfg, tissue)
    if panel_genes is None:
        panel_genes = sorted({g for genes in tax.core_genes.values() for g in genes})
        print(f"[vocab] {tissue}: no run panel supplied — marker coverage is computed against "
              "the identity-core union only and will understate on-panel markers", flush=True)

    pairs = protein_gene_pairs(panel_genes, tax=tax)
    print(coverage_report(space, pairs, tax), flush=True)

    df = crosswalk_frame(space, pairs, tax=tax)
    df.insert(0, "tissue", tissue)
    if write:
        out = cfg.vocab_crosswalk_csv
        out.parent.mkdir(parents=True, exist_ok=True)
        per_tissue = out.with_name(f"{out.stem}_{tissue}{out.suffix}")
        df.to_csv(per_tissue, index=False)
        # Merge into the canonical file rather than overwrite it, so a per-tissue re-run does
        # not silently delete the other tissue's rows.
        combined = df
        if out.exists():
            try:
                prev = pd.read_csv(out)
                if "tissue" in prev.columns:
                    prev = prev[prev["tissue"] != tissue]
                    combined = pd.concat([prev, df], ignore_index=True)
            except Exception:  # noqa: BLE001 - a corrupt/legacy file is replaced, not fatal
                combined = df
        combined.to_csv(out, index=False)
        print(f"\n[vocab] wrote {per_tissue} ({len(df)} rows) and {out} "
              f"({len(combined)} rows)", flush=True)
    return df


def main(argv: Optional[list[str]] = None) -> int:
    import argparse

    from ..config import load_config

    ap = argparse.ArgumentParser(description="Build the PhenoCycler<->Xenium vocab crosswalk")
    ap.add_argument("--config", default=None)
    ap.add_argument("--panel", default=None,
                    help="text file of gene symbols (one per line) from the actual run")
    ap.add_argument("--tissue", action="append", dest="tissues",
                    help="build for this tissue (repeatable); default is every configured one")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    genes = None
    if args.panel:
        genes = [ln.strip() for ln in Path(args.panel).read_text().splitlines() if ln.strip()]
    for tissue in (args.tissues or cfg.tissue_list):
        run_vocab(cfg, genes, tissue=tissue, write=not args.no_write)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
