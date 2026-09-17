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

import functools
import json
import pathlib
import warnings

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parents[1]                      # core/phenocycler
REPO = PKG.parents[1]                      # repo root
#: The retired run used when a donor has no expression in the stable root.
#:
#: New exports use ``data/expression`` since 2026-09-16. The original 20 donors still
#: live in this retired tree; ``expression_path`` selects the stable root independently
#: for each donor and warns when it has to read the retired copy.
#:
#: **It is regenerable, which the deleted run was not.**  That run could not be rebuilt: the
#: registry demanded ``measurement_source = qupath_norm`` on CD68, Iba1 and M2Gal3, which needs
#: ``QuPathNorm__<marker>`` columns that ``cells/`` has never carried, so ``expression`` failed
#: closed.  Those three now take the ``redsea`` default their ``redsea_corrected`` spillover
#: policy already implied.  **Measured value-neutral**: the rebuilt partitions reproduced the
#: old ones exactly on the donors that were checked -- all 62 numeric columns ``allclose``,
#: ``object_id`` order identical -- because the demand was never satisfiable in the first
#: place, so removing it restored what the data always was.  That is why ``method_version``
#: moved from ``selected-expression-v1`` to ``selected-expression-v2-qupath-norm-inputs``
#: while no value changed.
#:
#: Provenance: built 2026-08-28 from commit ``0f5ea04`` with a **dirty** worktree, so its
#: source is identified by ``source_sha256 0542e887...`` rather than by the commit alone.
#: Its manifests record a canonical ``output.root``, so unlike the two deleted runs it is
#: adoptable (``adoption.py`` refused those for the pre-monorepo path).
RUN_ID = "30f6689e4d27d295f890"

#: Assertions for the retired stage only; these do not describe the mixed-root cohort.
EXPRESSION_CONTENT_SHA256 = "35fa90fa66887329a3fb90ce23c0c944ede2fa73594fb1caac744f1bc8931ea0"
EXPRESSION_TOTAL_ROWS = 31_454_041
# Keep the historical root available to older consumers; new reads resolve per donor.
EXPRESSION = PKG / "data" / "runs" / RUN_ID / "expression"
EXPRESSION_FLAT = PKG / "data" / "expression"
REDSEA = PKG / "data" / "runs" / RUN_ID / "redsea"
OUT = REPO / "data" / "development" / "celltype"


def expression_path(donor: str) -> pathlib.Path:
    """Resolve one donor's parquet, stable root first, then the warned retired fallback.

    Numeric partition names are compared without leading zeroes: ``0112`` and ``112``
    name one donor. Multiple populated aliases in the same root are ambiguous.
    """
    def identity(value: str) -> str:
        return str(int(value)) if value.isdecimal() else value

    donor = identity(str(donor))
    for root in (EXPRESSION_FLAT, EXPRESSION):
        matches = sorted(
            part / "data.parquet" for part in root.glob("donor_id=*")
            if identity(part.name.removeprefix("donor_id=")) == donor
            and (part / "data.parquet").is_file()
        )
        if len(matches) > 1:
            raise ValueError(f"Ambiguous expression partitions for donor {donor}: {matches}")
        if matches:
            if root == EXPRESSION:
                warnings.warn(
                    f"Donor {donor} has no expression in {EXPRESSION_FLAT}; "
                    f"reading the retired run tree {matches[0]}",
                    RuntimeWarning, stacklevel=2,
                )
            return matches[0]
    raise FileNotFoundError(
        f"Expression unavailable for donor {donor}: no data.parquet in "
        f"{EXPRESSION_FLAT} or retired run {EXPRESSION}"
    )

#: IAPP only.  The registry records it ``kind: excluded, uses: failed_marker``.
#:
#: CD11b carries the same registry flag (``assay_qc``) and is NOT excluded here: per the
#: maintainer, it was set aside only for the retired RESTORE lane, and nothing in this lane
#: depends on that decision.  It matters because the rule table's ``Myeloid`` gate is
#: ``CD68 | CD11c | CD11b``, so dropping it would have silently removed a third of that
#: gate's evidence.
REGISTRY_EXCLUDED = ("IAPP",)

#: Stored in ``expression`` as the RAW UNCOMPENSATED QuPath *nucleus* mean
#: (``qupath:Nucleus__<m>``), not the REDSEA spillover-corrected whole-cell mean that the
#: other 52 markers carry.  A different measurement, and it must not be compared to the
#: others as though it were the same one.
NUCLEUS_MEASURED = ("DAPI", "FOXP3", "Ki67", "PCNA", "SOX2", "TCF_1", "TOX", "TP63")

#: Acquired on only some donors.  NOT a reason to drop them: the panel is four variants, so
#: marker availability is a per-donor property and a comparison that uses one of these simply
#: runs on fewer donors.  CD38 is on 15 of 22 (2 of them also carry b_Catenin1), b_Catenin1 is
#: on 6, and 3 donors carry neither.  Availability of one is not exclusive of the other, so
#: ``donors_with(["CD38"])`` and ``donors_with(["b_Catenin1"])`` overlap by those 2 donors.
#: :func:`donor_markers` is the per-donor authority; ``markers()`` is the cohort intersection,
#: for the places that genuinely need one common space.
#:
#: The 12 metabolic-panel markers added on 2026-09-16 belong here for exactly the same reason:
#: they are carried by 6533 and 112 ONLY -- 2 of 22 -- because those two slides were stained
#: with a 70-channel panel no other donor has. Leaving them out of this tuple put them in
#: ``markers()``, whose contract is "measured on every eligible donor", and the first caller to
#: build a donor frame from that list asked pyarrow for ASCT2 on a donor that never had it.
METABOLIC_PANEL_MARKERS = (
    "ASCT2", "ATP5A", "CD141", "CD21", "CS", "G6PD",
    "GLUT1", "HK1", "Histone_H3_(p_Ser28)", "IDH2", "LDHA", "SDHA",
)

PANEL_VARIANT_MARKERS = ("CD38", "b_Catenin1") + METABOLIC_PANEL_MARKERS

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
    # Alpha, beta and delta are DIFFERENT cells, so the hormones must not co-occur.  An
    # earlier draft declared these partners "islet endocrine, adjacent cells", which admits
    # they do not satisfy the definition of a partner -- labelling the same cell.  Their
    # measured 1.2-1.4x is the right answer, and it belongs on this side of the ledger.
    ("INS", "GCG"), ("GCG", "SST"), ("INS", "SST"),
]


def cohort() -> pd.DataFrame:
    """The 22 eligible donors, with the acquisition group their level may step with.

    ``acquisition_software`` is a proxy for the filter-responsivity calibration epoch and
    reproduces it exactly on this cohort (2.3.1 -> 2025-09-10, 2.4.0 -> 2026-03-12,
    2.4.1 -> 2026-06-11), so the group is derived from the cohort table alone.  Deriving it
    from ``data/pilots/acquisition_by_biomarker.parquet`` instead would lose donor 6523,
    whose slide is filed under another name.

    Five acquisition groups since 6533 and 112 were added 2026-09-16: 2.3.1/CD38,
    2.4.0/CD38, 2.4.1/b_Catenin1, 2.4.1/both (the new panel variant, carrying both CD38 and
    b_Catenin1) and 2.4.1/neither.
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


@functools.lru_cache(maxsize=None)
def _acquired_cached(donor: str, wanted: tuple) -> frozenset:
    """Which of ``wanted`` this donor actually ACQUIRED, from parquet statistics only.

    A column being present is not evidence the marker was measured. The expression schema is
    the union across the cohort, so a donor that never carried a marker still gets a column --
    filled with NaN, which is correct (a missing measurement is unavailable, not zero) but
    which made ``donor_markers`` claim the marker.

    Measured 2026-09-17: ``Ker8_18`` is entirely non-finite for 6533 and 112 and finite for all
    20 other donors, because those two slides were stained with a 70-plex carrying Keratin_5
    instead -- an acquisition fact, confirmed against their channel_map. Left unhandled it made
    every cell of both donors "incomplete", so ``donor_frame`` returned 0 complete cells.

    Uses per-row-group null counts rather than reading the column: exact for the all-null case,
    which is the one that matters, and costs no data read. Cached because donor_markers() is
    called per donor inside loops.
    """
    import pyarrow.parquet as _pq
    try:
        md = _pq.ParquetFile(expression_path(donor)).metadata
    except (FileNotFoundError, OSError):
        # Absence of the file is not evidence a marker is absent.
        return frozenset(wanted)
    names = {md.schema.column(i).name: i for i in range(md.num_columns)}
    out = set()
    for m in wanted:
        i = names.get(m)
        if i is None:
            continue
        nulls = sum(md.row_group(g).column(i).statistics.null_count
                    for g in range(md.num_row_groups)
                    if md.row_group(g).column(i).statistics is not None)
        if nulls < md.num_rows:
            out.add(m)
    return frozenset(out)


def _acquired(donor: str, wanted) -> frozenset:
    return _acquired_cached(str(donor), tuple(sorted(wanted)))


def donor_markers(donor: str | None = None) -> dict | list:
    """Which markers each donor actually carries.

    The cohort is four panel variants, so this is a property of the donor and not of the
    project.  A rule that names CD38 is evaluable on the 15 donors that have it and
    ``unavailable`` on the other 7 -- which is a smaller n, not an unencodable rule.  CD38
    and b_Catenin1 are not mutually exclusive: 2 donors (the "both" variant) carry both.
    """
    coh = pd.read_csv(REPO / "data" / "donor_acquisition_batch.csv", dtype={"donor_id": str})
    coh = coh[coh["eligible"]]
    base = [m for m in markers(include_variant=True) if m not in PANEL_VARIANT_MARKERS]
    variant = set(PANEL_VARIANT_MARKERS)
    out = {}
    for _, r in coh.iterrows():
        got = list(base)
        if r["has_CD38"]:
            got.append("CD38")
        if r["has_b_Catenin1"]:
            got.append("b_Catenin1")
        # The two has_* columns can only express two variants, and as of 2026-09-16 there is a
        # third: 6533 and 112 carry 12 metabolic markers no other donor has. Rather than add a
        # has_* column per panel change, ask the donor's own expression table which of the
        # remaining variant markers it actually carries. The table is the authority on what was
        # acquired; the CSV is the authority on the cohort. Absence of the file is not evidence
        # that a marker is absent, so a donor with no expression yet keeps the CSV-derived set.
        got = [m for m in got if m in _acquired(r["donor_id"], got)]
        rest = sorted(variant - {"CD38", "b_Catenin1"})
        if rest:
            got.extend(m for m in _acquired(r["donor_id"], rest) if m in set(rest))
        out[r["donor_id"]] = sorted(set(got))
    return out[donor] if donor is not None else out


def donors_with(markers_needed) -> list:
    """The donors on which every named marker was acquired -- the n for a comparison."""
    dm = donor_markers()
    need = set(markers_needed)
    return sorted(d for d, got in dm.items() if need <= set(got))


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
