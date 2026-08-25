"""What this lane treats as the panel, and the groups a marker's level may step with.

Everything here is DERIVED from files that predate any result in this lane -- the marker
registry, the QuPath manifest and the cohort acquisition table -- or DECLARED from
expression biology in writing.  Nothing is chosen after seeing a statistic.

The one judgement call, ``LINEAGE_PARTNERS``, is declared biology: which two markers label
the *same cell*.  It is checkable by someone who has never seen this cohort, and each entry
carries two ``NEGATIVE_CONTROLS`` -- markers that label a *different* cell -- so every
enrichment is reported against something that must not be enriched.
"""
from __future__ import annotations

import json
import pathlib

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parents[1]                      # core/phenocycler
REPO = PKG.parents[1]                      # repo root
RUN_ID = "3bf4fc8e391d007766d1"
EXPRESSION = PKG / "data" / "runs" / RUN_ID / "expression"
REDSEA = PKG / "data" / "runs" / RUN_ID / "redsea"
OUT = REPO / "data" / "development" / "celltype"

#: Markers the registry refuses for identity work, for reasons recorded in the registry
#: itself: CD11b is ``assay_qc`` and IAPP is ``failed_marker``.
REGISTRY_EXCLUDED = ("CD11b", "IAPP")

#: Stored in ``expression`` as the RAW UNCOMPENSATED QuPath *nucleus* mean
#: (``qupath:Nucleus__<m>``), not the REDSEA spillover-corrected whole-cell mean that the
#: other 52 markers carry.  A different measurement, and it must not be compared to the
#: others as though it were the same one.
NUCLEUS_MEASURED = ("DAPI", "FOXP3", "Ki67", "PCNA", "SOX2", "TCF_1", "TOX", "TP63")

#: Acquired on only some donors, so a per-donor statistic exists but a cohort one does not.
PANEL_VARIANT_MARKERS = ("CD38", "b_Catenin1")

#: Declared biology: (a, b) label the same cell, so their bright tails should overlap.
#: Written before any tail statistic was computed on this cohort.
LINEAGE_PARTNERS = [
    ("CD20", "CD79a", "B cell"),
    ("CD3e", "CD8", "T cell"),
    ("CD3e", "CD4", "CD4 T cell"),
    ("CD3e", "FOXP3", "regulatory T cell"),
    ("CD3e", "TOX", "T cell, nuclear"),
    ("CD3e", "TCF_1", "T cell, nuclear"),
    ("CD8", "Granzyme_B", "cytotoxic CD8 T cell"),
    ("CD3e", "Granzyme_B", "cytotoxic T cell"),
    ("CD68", "Iba1", "macrophage"),
    ("CD68", "CD163", "M2 macrophage"),
    ("Iba1", "CD206", "M2 macrophage"),
    ("CD11c", "CD209", "dendritic cell"),
    ("MPO", "CD66", "neutrophil"),
    ("E_cadherin", "Ker8_18", "epithelium"),
    ("Ker8_18", "Pan_Cytokeratin", "epithelium"),
    ("E_cadherin", "EpCAM", "epithelium"),
    ("Ker8_18", "EpCAM", "epithelium"),
    ("Ker8_18", "Keratin_5", "epithelium"),
    ("Keratin_5", "TP63", "basal epithelium"),
    ("CD31", "CD34", "blood endothelium"),
    ("CD31", "Caveolin", "blood endothelium"),
    ("B3TUBB", "CD56", "nerve"),
    ("SMA", "Vimentin", "mesenchyme"),
    ("INS", "GCG", "islet endocrine, adjacent cells"),
    ("GCG", "SST", "islet endocrine, adjacent cells"),
    ("Ki67", "PCNA", "proliferating cell"),
]

#: One well-populated marker per broad compartment, in a fixed order.  Negative controls are
#: DERIVED from this rather than hand-listed -- a hand-picked control set could drift toward
#: pairs that happen to behave; a rule cannot.  A control must come from a DIFFERENT broad
#: compartment than the marker: cells of one compartment sit together in tissue, and
#: PhenoCycler spillover crosses touching cell boundaries, so a same-compartment marker is a
#: poor negative control by construction.  Measured on the first cohort run: non-partner
#: pairs where both markers are immune have median overlap 1.61x (75th pct 5.82x) against
#: 0.00x (75th 0.88x) for everything else.
CONTROL_REPRESENTATIVES = {
    "Epithelial": "Ker8_18",
    "Immune": "CD3e",
    "Mesenchymal": "SMA",
    "Endothelial": "CD31",
    "Neural": "B3TUBB",
}

#: Extra declared non-partner pairs, kept because they are the ones reported in the plan and
#: because two of them are informative in their own right: CD31/Podoplanin should be
#: ANTI-correlated (blood vs lymphatic endothelium are different cells), and DAPI/SMA is the
#: "both present on many cells, but not the same cells" case.
DECLARED_CONTROLS = [
    ("CD20", "SMA"), ("CD20", "INS"), ("TOX", "E_cadherin"), ("FOXP3", "SMA"),
    ("TP63", "CD3e"), ("CD31", "Podoplanin"), ("DAPI", "SMA"),
]


def cohort() -> pd.DataFrame:
    """The 20 eligible donors, with the acquisition group their level may step with.

    ``acquisition_software`` is a proxy for the filter-responsivity calibration epoch and
    reproduces it exactly on this cohort (2.3.1 -> 2025-09-10, 2.4.0 -> 2026-03-12,
    2.4.1 -> 2026-06-11), so the group is derived from the cohort table alone.  Deriving it
    from ``data/pilots/acquisition_by_biomarker.parquet`` instead would lose donor 6523,
    whose slide is filed under another name.
    """
    df = pd.read_csv(REPO / "data" / "donor_acquisition_batch.csv", dtype={"donor_id": str})
    df = df[df["eligible"]].copy()
    df["epoch"] = df["acquisition_software"].str.replace("Fusion ", "", regex=False)
    df["acq_group"] = df["epoch"] + "/" + df["panel_variant"]
    return df[["donor_id", "batch", "epoch", "panel_variant", "acq_group",
               "disease.status", "n_channels"]].rename(columns={"disease.status": "status"})


def markers(include_variant: bool = False) -> list:
    """Registry markers measured on every eligible donor, minus the registry's exclusions.

    ``include_variant`` adds CD38 / b_Catenin1, which exist on some donors only and so have
    a per-donor statistic but no cohort-wide one.
    """
    reg = json.loads((PKG / "phenocycler" / "marker_registry.json").read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    names = [m["name"] for m in reg if m["name"] not in REGISTRY_EXCLUDED]
    if not include_variant:
        names = [n for n in names if n not in PANEL_VARIANT_MARKERS]
    return sorted(names)


def compartments() -> dict:
    """marker -> broad compartment, for the negative-control rule.

    The registry assigns one to 30 markers and declines for the rest (every state marker is
    ``uses=['deferred_state'], broad_compartment=null``).  A marker it declines takes the
    compartment of its declared partners -- Granzyme_B has none in the registry but partners
    CD8 and CD3e, so it is Immune.  Derived, so it cannot be steered per marker.
    """
    reg = json.loads((PKG / "phenocycler" / "marker_registry.json").read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    comp = {m["name"]: m.get("broad_compartment") for m in reg}
    partners = partner_map()
    for m, c in list(comp.items()):
        if c:
            continue
        near = [comp.get(p) for p, _ in partners.get(m, []) if comp.get(p)]
        comp[m] = near[0] if near else None
    return comp


def partner_map() -> dict:
    """marker -> [(partner, what they jointly label), ...], both directions."""
    out: dict = {}
    for a, b, what in LINEAGE_PARTNERS:
        out.setdefault(a, []).append((b, what))
        out.setdefault(b, []).append((a, what))
    return out


def control_map(n: int = 2) -> dict:
    """marker -> [non-partner marker, ...], derived so every marker gets ``n`` controls.

    A partner enrichment with no control beside it is unreadable, so no marker may go
    without one.  ``DECLARED_CONTROLS`` are merged in first, then representatives fill up to
    ``n``; if a marker somehow cannot reach ``n``, it is returned short and the caller
    records that rather than silently reporting a bare enrichment.
    """
    partners = {m: {p for p, _ in v} for m, v in partner_map().items()}
    comp = compartments()
    out: dict = {}
    for a, b in DECLARED_CONTROLS:
        out.setdefault(a, []).append(b)
        out.setdefault(b, []).append(a)
    for m in set(partners) | set(out):
        got = out.setdefault(m, [])
        for c, rep in CONTROL_REPRESENTATIVES.items():
            if len(got) >= n:
                break
            if c == comp.get(m) or rep == m or rep in partners.get(m, ()) or rep in got:
                continue
            got.append(rep)
    return out
