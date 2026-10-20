"""
What each tissue supports — the single place that knows pancreas from lymph node.

Since D2 (2026-10-20) the values come from the tissue profiles
(``core/multiplex5000/profiles/tissues/*.yaml``); this module reads them and holds the logic.

The package supports pancreas and pancreatic lymph node, but support is not a claim that both
modalities measured every section. The current PhenoCycler run has 20 pancreas sections and
no pLN sections; the manifest keeps Xenium-only pLN rows explicit. When both tissues are
present, they diverge in two specific biological ways:

**Structures.** Islets exist in pancreas. Lymph nodes have no islets, so hormone-seeded DBSCAN
has nothing to seed on and islet matching (S6) has nothing to match. Rather than invent a
substitute landmark, the structure level is simply **not run** for lymph node — the
registration-free levels (donor, niche) and the registered grid still apply, and those are
the levels that were always going to carry the lymph-node analysis anyway.

That is a deliberate choice, not an omission. B-cell follicles are the obvious candidate
landmark and the code would support them (``call_landmarks`` is already generic over a seed
lineage), but a follicle called from a 5K RNA panel on one side and a 55-plex protein panel on
the other is a much weaker correspondence than an islet, and shipping it would invite
structure-level claims the data cannot support. ``STRUCTURE_SEEDS`` below is where it would
go if that changes.

**Identity panels.** XeniumPanelExplorer ships curated identity cores per tissue, and has
``pancreas`` and ``thymus`` — no lymph node. But the cores a lymph node needs are exactly the
*shared* ones: T cell, B/plasma, NK/ILC, myeloid, granulocyte, endothelial, pericyte, neural,
epithelial. The three pancreas-specific cores (endocrine, exocrine, fibroblast/stellate) do
not apply. So a tissue with no directory of its own resolves to the shared pool, which is the
correct answer rather than a fallback.
"""

from __future__ import annotations

from typing import Optional

from .. import profiles as _profiles


class TissueError(ValueError):
    """An unknown or unsupported tissue."""


def _load_tissue_profiles() -> dict[str, _profiles.TissueProfile]:
    """Every tissue profile, by id. Since D2 (2026-10-20) they are the authority, not a mirror."""
    loaded = _profiles.load_all()
    found = {pid: prof for (kind, pid), prof in loaded.items() if kind == "tissue"}
    if not found:
        raise TissueError(f"no tissue profiles under {_profiles.PROFILES_DIR}; the tissue names "
                          "are read from there, and MultiPlex5000 runs from a checkout that has them")
    return dict(sorted(found.items()))


#: The tissue profiles (``profiles/tissues/*.yaml``), by id, in id order.
PROFILES: dict[str, _profiles.TissueProfile] = _load_tissue_profiles()

#: Canonical tissue identifiers. Anything else is rejected rather than silently creating a third
#: tissue nobody configured. These two are named because the code holds logic that differs by
#: tissue (the plugins); each must have a profile.
PANCREAS = "pancreas"
LYMPH_NODE = "pancreatic_lymph_node"
for _named in (PANCREAS, LYMPH_NODE):
    if _named not in PROFILES:
        raise TissueError(f"tissue {_named!r} has no profile under {_profiles.PROFILES_DIR}")
TISSUES: tuple[str, ...] = tuple(PROFILES)

#: Human-readable, for figures and reports.
DISPLAY_NAME = {t: prof.display_name for t, prof in PROFILES.items()}

#: Every accepted spelling -> its tissue. Lowercased, with spaces and hyphens as underscores.
_SPELLINGS: dict[str, str] = {}
for _t, _prof in PROFILES.items():
    for _spelling in (_t, *_prof.aliases):
        if _SPELLINGS.setdefault(_spelling, _t) != _t:
            raise TissueError(f"{_spelling!r} is an alias of two tissues: "
                              f"{_SPELLINGS[_spelling]} and {_t}")

#: ROI prefix -> tissue. One donor has one pancreas ROI and up to three lymph-node ROIs.
ROI_PREFIX_TO_TISSUE = {prefix: t for t, prof in PROFILES.items() for prefix in prof.roi_prefixes}

#: The XeniumPanelExplorer tissue directory backing each tissue's identity cores.
#: ``None`` means "no tissue-specific cores — use the shared pool", which is the right answer
#: for lymph node rather than a degraded one.
PANEL_TISSUE_DIR: dict[str, Optional[str]] = {t: prof.panel_tissue_dir
                                              for t, prof in PROFILES.items()}

#: Identity cores each tissue owns, which must NOT be applied to any other tissue.
TISSUE_ONLY_PANELS: dict[str, frozenset[str]] = {t: frozenset(prof.tissue_only_panels)
                                                 for t, prof in PROFILES.items()}
#: The pancreas-specific cores, under the name the integration layer has always used.
PANCREAS_ONLY_PANELS = TISSUE_ONLY_PANELS[PANCREAS]

#: Structure kinds each tissue supports. Empty tuple = the structure level does not run, and
#: everything gated on it (S6 islet matching, the islet QC metric) is skipped for that tissue.
STRUCTURE_KINDS: dict[str, tuple[str, ...]] = {t: tuple(prof.structures)
                                               for t, prof in PROFILES.items()}

#: Seed definition per structure kind, as (common lineage, description). The registry exists
#: so adding follicles later is a one-line change here plus a config flag, not a redesign.
STRUCTURE_SEEDS: dict[str, tuple[str, str]] = {
    "islet": ("Endocrine", "hormone-positive endocrine cells"),
    "duct": ("Exocrine", "ductal/acinar epithelium"),
    "vessel": ("Vascular", "endothelium"),
    # "follicle": ("B_Plasma", "B-cell follicles") — deliberately not enabled; see module docstring.
}

#: Common lineages worth comparing per tissue (the profile's ``comparable_lineages``). A lymph
#: node has no endocrine or exocrine compartment, so reporting a near-zero fraction for those
#: would be noise dressed as a finding; a pancreas has all of them. Granulocyte is absent from
#: both: neither modality calls it at broad level, so it is `fine_only` on both sides and
#: `donor.comparable_for` would drop it anyway.
COMPARABLE_LINEAGES: dict[str, tuple[str, ...]] = {t: tuple(prof.comparable_lineages)
                                                   for t, prof in PROFILES.items()}


def normalize(tissue: str) -> str:
    """Canonicalise a tissue name, accepting the obvious spellings.

    Lenient on input (``pLN``, ``lymph node``, ``pancreatic_lymph_node`` all work) and strict
    on output, so config files and CLI flags stay forgiving without letting a typo become a
    third tissue that quietly matches nothing.
    """
    t = str(tissue).strip().lower().replace(" ", "_").replace("-", "_")
    if not t:
        raise TissueError("tissue is empty; it must be declared when the dataset is imported")
    if t in _SPELLINGS:
        return _SPELLINGS[t]
    raise TissueError(f"unknown tissue {tissue!r}; expected one of {TISSUES}")


def for_roi(roi: str) -> str:
    """ROI name -> tissue. Returns '' for an unrecognised ROI rather than guessing."""
    r = str(roi).strip().lower()
    for prefix, tissue in ROI_PREFIX_TO_TISSUE.items():
        if r.startswith(prefix):
            return tissue
    return ""


#: Lymph nodes per donor when the caller does not say. Every donor in this cohort has a
#: pancreas and 1-4 nodes; `PipelineConfig.max_pln` is the configurable authority and passes
#: its value in. Kept here as well so this module stands alone for callers without a config.
DEFAULT_MAX_PLN = 4


def rois(tissue: str, *, max_pln: int = DEFAULT_MAX_PLN) -> tuple[str, ...]:
    """Known ROI names for a tissue: ``panc``, or ``pln_1..max_pln``.

    The node count is a parameter because it describes the cohort, not the code. It bounds
    only the *enumeration* used when nothing is on disk yet — `sections.parse` accepts any
    index, so a `pln_5` partition is readable; it is the stages that must be told to look
    for it, and they refuse it loudly rather than skipping it.
    """
    t = normalize(tissue)
    if t == PANCREAS:
        return ("panc",)
    n = max(1, int(max_pln))
    return tuple(f"pln_{i}" for i in range(1, n + 1))


def structure_kinds(tissue: str) -> tuple[str, ...]:
    """Structure kinds to extract. Empty means the structure level does not apply."""
    return STRUCTURE_KINDS.get(normalize(tissue), ())


def has_structures(tissue: str) -> bool:
    """Whether structure extraction and islet matching are meaningful for this tissue."""
    return bool(structure_kinds(tissue))


def comparable_lineages(tissue: str) -> tuple[str, ...]:
    """Lineages worth reporting a cross-modal comparison for, in this tissue."""
    return COMPARABLE_LINEAGES.get(normalize(tissue), ())


def panel_tissue_dir(tissue: str) -> Optional[str]:
    """XeniumPanelExplorer tissue directory, or ``None`` for shared-cores-only."""
    return PANEL_TISSUE_DIR.get(normalize(tissue))


def applies(tissue: str, panel_stem: str) -> bool:
    """Whether an identity-core panel applies to this tissue.

    Lymph node gets the shared immune / vascular / stromal / neural cores and none of the
    pancreas-specific ones. Scoring a lymph-node cell against ``01_pancreas_endocrine`` would
    not merely be uninformative — it would put weight on an endocrine axis that cannot exist
    there, and the broad call is an argmax over those axes.
    """
    t = normalize(tissue)
    return not any(panel_stem in owned for other, owned in TISSUE_ONLY_PANELS.items() if other != t)


def describe(tissue: str) -> str:
    """One-line summary of what this tissue supports — used in stage logs."""
    t = normalize(tissue)
    kinds = structure_kinds(t)
    struct = ", ".join(kinds) if kinds else "none (registration-free + grid levels only)"
    return (f"{DISPLAY_NAME[t]}: rois={'/'.join(rois(t))} | structures={struct} | "
            f"{len(comparable_lineages(t))} comparable lineages")
