"""Which markers speak for which broad lineage.

This is a statement of what each antibody *means*, made before any result is seen.  It is not a
set of thresholds and it never decides a call on its own: clusters are found without reference to
it, and it is used only to name them.

Every marker on the panel that can speak to identity appears here.  The panel's own design is the
authority for that -- a marker was put on the panel because it marks something -- so the only
channels held back are the two that cannot be read at all: ``IAPP`` (the registry's own
``failed_marker``) and ``DAPI`` (the nuclear counterstain, kept as the segmentation and brightness
reference).  In particular B3TUBB, CD34, CD56, CD66 and CD11b are used here; an earlier lane
dropped them, and with B3TUBB it dropped the Neural lineage entirely and folded granulocytes into
the myeloid label.

State markers (proliferation, exhaustion, activation) are deliberately absent: they describe what
a cell is *doing*, not what it is.  They are carried through the substrate and reported per
cluster, so a Ki67-high T cell is still a T cell.
"""

from __future__ import annotations

#: lineage -> the markers that mark it.  Ordered most to least specific within each lineage.
LINEAGE_MARKERS: dict[str, tuple[str, ...]] = {
    # Islet hormones.  Granule-localised, read from the cytoplasm.
    "Endocrine":    ("INS", "GCG", "SST"),
    # Exocrine pancreas.  The panel carries no acinar-specific and no ductal-specific marker, so
    # this stays one label; Keratin_5/TP63 are basal-associated and inform it but cannot split it.
    "Exocrine":     ("Pan_Cytokeratin", "Ker8_18", "E_cadherin", "EpCAM", "Keratin_5", "TP63"),
    "Vascular":     ("CD31", "Caveolin", "Podoplanin", "CD34"),
    "Stromal":      ("SMA", "Vimentin", "Collagen_IV"),
    "T_NK":         ("CD3e", "CD8", "CD4", "FOXP3", "CD56"),
    "B_Plasma":     ("CD20", "CD79a"),
    "Myeloid":      ("CD68", "Iba1", "CD163", "CD206", "CD11c", "CD209", "M2Gal3", "CD11b"),
    # Neutrophils.  MPO is the panel's cleanest immune channel and is a granulocyte marker, not a
    # general myeloid one; CD66 marks the same compartment.
    "Granulocyte":  ("MPO", "CD66"),
    "Neural":       ("B3TUBB",),
}

#: Markers that may NOT carry a lineage call on their own.  They still enter the z statistic
#: and still describe the cluster -- they simply cannot be the marker whose lift and positive
#: fraction admit the label, because each is measured in this donor set to be too broad or too
#: shared to name a cell by itself:
#:
#:   M2Gal3, CD11b  -- declared last in Myeloid's own ordering; CD11b is carried by
#:                     granulocytes and NK cells as well, and on donor 6442 clusters whose
#:                     CD68 and Iba1 sit BELOW background were named Myeloid on CD11b alone.
#:   FOXP3          -- above the positivity model's own boundary in 90.4% of cells.
#:   CD56           -- shared with the Neural lineage; on 6436 a cluster with CD3e at +0.06
#:                     was named T_NK on CD56 alone.
#:   TP63           -- above the boundary in 98.6% of cells.
#:
#: A lineage whose evidence rests only on these is Ambiguous, not named.
SUPPORTING_ONLY: frozenset = frozenset({"M2Gal3", "CD11b", "FOXP3", "CD56", "TP63"})


#: Subtypes resolved *within* a parent lineage, never across the whole section.  A marker that is
#: flat in its marginal distribution can still be strongly informative inside its parent -- CD4 is
#: 34x over-represented among CD3e-high cells on 6442 -- which is how these are read.
#:
#: THESE MUST BE MUTUALLY EXCLUSIVE.  They are resolved by argmax, so a subtype that is really a
#: SUBSET of another one wins cells that belong to its parent.  Treg was listed here as a peer of
#: CD4_T until it produced 11,463 Tregs, a third of every T cell on donor 6442 -- because a Treg
#: *is* a CD4 T cell carrying FOXP3, and because Nimbus scores FOXP3 above its own boundary in
#: 90.4% of cells.  Anything that refines a subtype rather than competing with it belongs in
#: :data:`SUBTYPE_MODIFIERS`.
SUBTYPE_MARKERS: dict[str, dict[str, tuple[str, ...]]] = {
    "Endocrine":  {"Beta": ("INS",), "Alpha": ("GCG",), "Delta": ("SST",)},
    "T_NK":       {"CD8_T": ("CD8",), "CD4_T": ("CD4",), "NK": ("CD56",)},
    "B_Plasma":   {"B": ("CD20",), "Plasma": ("CD79a",)},
    "Myeloid":    {"Macrophage": ("CD68", "Iba1"), "Dendritic": ("CD11c", "CD209")},
    "Vascular":   {"Endothelial": ("CD31",), "Lymphatic": ("Podoplanin",)},
}

#: Refinements *of* a subtype, reported in their own column rather than competing for the label.
#: Scored within the cells that already carry the parent subtype, so the comparison is against
#: that population and nothing else.
#:
#: NAMED FOR WHAT IS MEASURED, NOT FOR WHAT IT MIGHT BE.  These are the upper quartile of the
#: parent population on the listed markers -- so the count is 25% of the parent BY CONSTRUCTION
#: and carries no information beyond the ranking.  Calling that quartile "Treg" or "M2" would
#: dress a percentile as a biological call: on donor 6442 it would have reported 2,909 Tregs, and
#: the only thing that number measures is that 11,634 cells were called CD4_T.  A real Treg call
#: needs a FOXP3 threshold this panel cannot support -- Nimbus scores FOXP3 above its own
#: boundary in 90.4% of cells -- so the flag says which cells rank highest and the reviewer
#: judges from the image what that means.
SUBTYPE_MODIFIERS: dict[str, tuple[str, tuple[str, ...]]] = {
    "CD4_T":      ("FOXP3_high", ("FOXP3",)),
    "Macrophage": ("CD163_CD206_high", ("CD163", "CD206")),
}

#: Cells that satisfy no lineage, or satisfy more than one with no separation between them.
AMBIGUOUS = "Ambiguous"
#: Cells with no measurement (geometry QC).  Carries the QC reason.
UNAVAILABLE = "Unavailable"


def identity_markers() -> list[str]:
    """Every marker used to name a lineage, de-duplicated, in declaration order."""
    seen: dict[str, None] = {}
    for markers in LINEAGE_MARKERS.values():
        for m in markers:
            seen.setdefault(m, None)
    return list(seen)
