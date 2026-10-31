"""Which REDSEA compartment each marker is measured in, and why.

The run stores four compartment means per marker: ``Nucleus``, ``Cytoplasm``, ``Membrane``
(the boundary band) and ``Cell`` (the whole polygon).  Everything downstream of the pipeline
has historically read ``Cell``, and that is the wrong basis for most of the panel: the whole-cell
mean averages a surface protein over nuclear pixels that cannot contain it, which both dilutes
the signal and imports nuclear brightness into it.

Measured on donors 6436 and 6442 (all analysis-eligible cells), correlation with nuclear DAPI
drops sharply when the marker is read in its own compartment rather than over the whole cell --
CD163 0.53 -> 0.15, Keratin_5 0.72 -> 0.38, CD8 0.37 -> 0.14, CD79a 0.45 -> 0.15 -- and the
bright tail separates correspondingly better.  The assignment below is made from where each
protein actually sits, not by picking the compartment that scored best; the measurement is the
check on the assignment, not its source.

Coverage matters for the choice between ``Cytoplasm`` and ``Membrane``.  7.3% (6436) and 9.0%
(6442) of analysis-eligible cells have *zero* cytoplasm pixels, so a cytoplasm mean is NaN for
them, while the boundary band is never empty (1st percentile 39 px, median ~71 px).  Surface
proteins therefore read from ``Membrane``: it is where they are, and it is complete.
"""

from __future__ import annotations

#: Nuclear proteins.  These are also the markers the registry declares ``compartment: Nucleus``
#: and ``spillover_policy: uncompensated`` -- REDSEA does not touch them, so ``Cell`` for these
#: is a whole-polygon average of a signal confined to the nucleus.
NUCLEUS = (
    "FOXP3", "TP63", "Ki67", "PCNA", "SOX2", "TCF_1", "TOX",
)

#: Cytoplasmic, granule and secreted proteins, plus the cytoskeletal filaments.  Hormones live in
#: secretory granules; cytokeratins, vimentin and smooth-muscle actin are cytoplasmic filaments.
CYTOPLASM = (
    "INS", "GCG", "SST",
    "MPO", "Granzyme_B", "IDO1", "iNOS", "Iba1", "M2Gal3",
    "Ker8_18", "Keratin_5", "Pan_Cytokeratin", "Vimentin", "SMA", "B3TUBB",
    "Collagen_IV", "Beta_actin", "Bcl_2",
    # Spleen PhenoCycler panel (2026-10-30): a secreted cytokine read where it is made, and a
    # cytokeratin filament, as the marker_compartment_ruling.csv rows for them say.
    "IFNG", "Keratin_14",
)

#: Everything else on the panel is a surface or membrane-associated protein and is read from the
#: boundary band.  Listed explicitly rather than left as a default so a new marker has to be
#: placed deliberately.
MEMBRANE = (
    "CD3e", "CD4", "CD8", "CD20", "CD79a", "CD11b", "CD11c", "CD31", "CD34", "CD38",
    "CD39", "CD44", "CD56", "CD57", "CD66", "CD68", "CD99", "CD107a", "CD163", "CD206",
    "CD209", "Caveolin", "E_cadherin", "EpCAM", "Podoplanin", "HLA_A", "HLA_DR",
    "ICOS", "LAG3", "PD_1", "PD_L1", "VISTA", "b_Catenin1",
    # Spleen PhenoCycler panel (2026-10-30): surface CD14, and CD45 with its RO isoform.
    "CD14", "CD45", "CD45RO",
)

#: Not usable for cell typing.  ``IAPP`` is the registry's own ``failed_marker``; ``DAPI`` is the
#: nuclear counterstain and is kept for segmentation QC and as the brightness reference against
#: which every other marker's nuclear coupling is measured.  These two, and only these two.
NOT_FOR_TYPING = ("IAPP", "DAPI")

_ASSIGNMENT: dict[str, str] = {}
for _m in NUCLEUS:
    _ASSIGNMENT[_m] = "Nucleus"
for _m in CYTOPLASM:
    _ASSIGNMENT[_m] = "Cytoplasm"
for _m in MEMBRANE:
    _ASSIGNMENT[_m] = "Membrane"


def compartment_for(marker: str) -> str:
    """The compartment ``marker`` is measured in.

    Raises for a marker with no declared placement, so adding a channel to the panel forces a
    deliberate decision rather than silently falling back to the whole-cell mean.
    """
    try:
        return _ASSIGNMENT[marker]
    except KeyError:
        raise KeyError(
            f"{marker!r} has no declared compartment; add it to NUCLEUS, CYTOPLASM or MEMBRANE "
            f"in scripts/typing/compartments.py with a reason"
        ) from None


def column_for(marker: str) -> str:
    """The REDSEA column to read for ``marker`` (e.g. ``Membrane__CD3e``)."""
    return f"{compartment_for(marker)}__{marker}"


def typing_markers(available: list[str]) -> list[str]:
    """``available`` restricted to markers usable for typing, in the given order."""
    skip = set(NOT_FOR_TYPING)
    return [m for m in available if m not in skip and m in _ASSIGNMENT]
