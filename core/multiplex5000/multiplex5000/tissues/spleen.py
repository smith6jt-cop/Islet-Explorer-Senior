"""Spleen protein reference derived from the selected 6551 hierarchical sheet.

This loads marker knowledge and reports unavailable rules. It does not assign
cell types before normalization, controls and the Stage 06 gate are reviewed.
"""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

import yaml

from .lymph_node import Rule, RuleSource, RuleSourceError
from .lymph_node import load_rule_source as load_pln_rule_source


def normalized_markers(manifest_path: Path, *, table_path: Path | None = None,
                       expected_entries: set[str] | None = None) -> frozenset[str]:
    """Audit the recorded HPG scales; return markers scaled for every PC67 section.

    The 67 measured export columns must not be mistaken for 67 normalized
    markers. This table currently covers only the 33 classifier channels.
    ``expected_entries`` may bind the audit to Stage 00's exact paired sections.
    """
    manifest = yaml.safe_load(Path(manifest_path).read_text())
    if not isinstance(manifest, dict) or manifest.get("tissue_id") != "spleen":
        raise RuleSourceError("expected the spleen panel manifest")
    spec = manifest.get("protein_normalization")
    if not isinstance(spec, dict):
        raise RuleSourceError("spleen normalization source is not declared")
    path = Path(table_path) if table_path is not None else Path(spec["path"])
    if hashlib.sha256(path.read_bytes()).hexdigest() != spec["sha256"]:
        raise RuleSourceError("spleen normalization table differs from the recorded source")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["image", "group", "marker", "reference", "scale",
                                 "shift_log", "provenance", "reagent_changed"]:
            raise RuleSourceError("unexpected spleen normalization columns")
        rows = list(reader)
    if len(rows) != spec["rows"]:
        raise RuleSourceError("spleen normalization row count changed")
    entries: dict[str, set[str]] = {}
    seen: set[tuple[str, str]] = set()
    for row in rows:
        image, marker = row["image"], row["marker"]
        key = (image, marker)
        if key in seen or not image or not marker:
            raise RuleSourceError("repeated or blank spleen image-marker scale")
        seen.add(key)
        try:
            scale = float(row["scale"])
        except (TypeError, ValueError) as err:
            raise RuleSourceError(f"{image}/{marker}: invalid UniFORM scale") from err
        if not math.isfinite(scale) or scale <= 0:
            raise RuleSourceError(f"{image}/{marker}: invalid UniFORM scale")
        if "_PC67_" in image:
            if row["group"] != "PC67":
                raise RuleSourceError(f"{image}: wrong normalization group")
            entries.setdefault(image, set()).add(marker)
    if len(entries) != spec["pc67_entries"]:
        raise RuleSourceError("spleen normalization does not cover the eight PC67 sections")
    if expected_entries is not None and set(entries) != expected_entries:
        raise RuleSourceError("spleen normalization entries differ from Stage 00 pairs")
    markers = set.intersection(*entries.values())
    if any(entry != markers for entry in entries.values()):
        raise RuleSourceError("PC67 sections do not share the same normalized marker set")
    if len(markers) != spec["markers_per_pc67_entry"]:
        raise RuleSourceError("PC67 normalized marker count changed")
    if not markers <= observed_markers(manifest_path):
        raise RuleSourceError("normalized marker is absent from the PC67 export")
    return frozenset(markers)


def observed_markers(manifest_path: Path) -> frozenset[str]:
    """Read the versioned PC67 export-column inventory, including its aliases."""
    manifest_path = Path(manifest_path)
    manifest = yaml.safe_load(manifest_path.read_text())
    if not isinstance(manifest, dict) or manifest.get("tissue_id") != "spleen":
        raise RuleSourceError("expected the spleen panel manifest")
    panel = manifest.get("protein_panel")
    if not isinstance(panel, dict) or panel.get("checked_exports") != 8:
        raise RuleSourceError("the eight PC67 export headers have not been declared")
    panel_path = manifest_path.parent / panel["file"]
    header = panel_path.read_bytes()
    if hashlib.sha256(header).hexdigest() != panel.get("source_header_sha256"):
        raise RuleSourceError("PC67 export header differs from the recorded source")
    if header.count(b"\n") != 1:
        raise RuleSourceError("PC67 inventory must contain exactly one export header")
    columns = header.decode("utf-8").strip().split("\t")
    if columns[:3] != ["object_id", "x_um", "y_um"]:
        raise RuleSourceError("PC67 header does not name the expected cell identity")
    names = [column.removeprefix("Cell: ").removesuffix(": Mean") for column in columns
             if column.startswith("Cell: ") and column.endswith(": Mean")]
    if len(names) != 67 or len(set(names)) != len(names):
        raise RuleSourceError("expected 67 unique PC67 marker columns")
    return frozenset(names)


def load_rule_source(manifest_path: Path) -> RuleSource:
    """Apply the recorded spleen adaptation to the one selected source sheet."""
    manifest_path = Path(manifest_path)
    manifest = yaml.safe_load(manifest_path.read_text())
    if not isinstance(manifest, dict) or manifest.get("tissue_id") != "spleen":
        raise RuleSourceError("expected the spleen panel manifest")
    spec = manifest.get("rule_table")
    if not isinstance(spec, dict) or spec.get("assay") != "phenocycler":
        raise RuleSourceError("spleen rules require the phenocycler assay")
    source_manifest = (manifest_path.parent / spec["source_manifest"]).resolve()
    source = load_pln_rule_source(source_manifest)
    if source.sha256 != spec.get("source_sha256"):
        raise RuleSourceError("spleen source digest differs from the selected 6551 sheet")
    adaptation = spec.get("adaptation")
    if not isinstance(adaptation, dict):
        raise RuleSourceError("spleen adaptation is missing")
    panel = observed_markers(manifest_path)
    names = source.by_phenotype
    branches = set(adaptation.get("drop_branches", []))
    individual = set(adaptation.get("drop_rules", []))
    if not branches <= names.keys() or not individual <= names.keys():
        raise RuleSourceError("spleen adaptation drops an unknown source phenotype")
    if branches & individual:
        raise RuleSourceError("a dropped branch is also named as an individual rule")
    dropped = set(branches)
    for rule in source.rules:
        current = rule.parent
        while current in names:
            if current in branches:
                dropped.add(rule.phenotype)
                break
            current = names[current].parent
    dropped.update(individual)
    kept = [rule for rule in source.rules if rule.phenotype not in dropped]
    available_names = {rule.phenotype for rule in kept}
    for rule in kept:
        if rule.parent not in {"all", "any"} and rule.parent not in available_names:
            raise RuleSourceError(f"{rule.phenotype}: parent dropped from spleen rules")
    aliases = {**source.panel_aliases, **manifest["protein_panel"].get("aliases", {})}
    kept, omitted = _omit_absent_anypos(kept, adaptation.get("omit_absent_anypos", {}), panel, aliases)
    additions = adaptation.get("add_rules", [])
    if not isinstance(additions, list):
        raise RuleSourceError("spleen add_rules must be a list")
    for item in additions:
        if not isinstance(item, dict) or set(item) != {"parent", "phenotype", "constraints"}:
            raise RuleSourceError("invalid spleen rule addition")
        parent, phenotype, constraints = item["parent"], item["phenotype"], item["constraints"]
        if parent not in available_names or phenotype in names or phenotype in available_names:
            raise RuleSourceError(f"{phenotype}: invalid new rule or parent")
        if not isinstance(constraints, dict) or not constraints:
            raise RuleSourceError(f"{phenotype}: new rule needs marker constraints")
        if set(constraints) - panel or set(constraints.values()) - {"pos", "allpos", "anypos", "allneg"}:
            raise RuleSourceError(f"{phenotype}: new rule uses absent markers or invalid tokens")
        kept.append(Rule(parent, phenotype, tuple(constraints.items())))
        available_names.add(phenotype)
    declared_roots = adaptation.get("candidate_lineage_roots")
    if declared_roots is None:
        candidate_roots = tuple(name for name in source.candidate_lineage_roots
                                if name in available_names)
    else:
        lineage_roots = {rule.phenotype for rule in kept if rule.parent == source.lineage_root_parent}
        if (not isinstance(declared_roots, list) or len(declared_roots) != len(set(declared_roots))
                or not set(declared_roots) <= lineage_roots):
            raise RuleSourceError("spleen candidate lineage roots must be distinct kept lineage roots")
        candidate_roots = tuple(declared_roots)
    if not candidate_roots:
        raise RuleSourceError("spleen adaptation removed every candidate lineage root")
    result = RuleSource(source.sha256, tuple(kept), source.omitted + omitted, aliases,
                        source.lineage_root_parent, source.state_root_parent,
                        candidate_roots)
    for rule in result.rules:
        result.inherited_markers(rule.phenotype)
    _structure_only(adaptation, result)
    return result


def _omit_absent_anypos(rules: list[Rule], spec, panel: frozenset[str],
                        aliases: dict[str, str]) -> tuple[list[Rule], tuple[tuple[str, str], ...]]:
    """Leave the listed unacquired markers out of their rules' OR groups (maintainer, 2026-11-01).

    Only an ``anypos`` clause on a marker absent from the panel (and its alias) may be left out, and at
    least one ``anypos`` clause must remain, so the rule is read on its measured members. Every other
    missing marker still makes its rule unavailable."""
    if not isinstance(spec, dict):
        raise RuleSourceError("spleen omit_absent_anypos must map rules to markers")
    by_name = {rule.phenotype: rule for rule in rules}
    omitted: list[tuple[str, str]] = []
    for phenotype, markers in spec.items():
        rule = by_name.get(phenotype)
        if rule is None or not isinstance(markers, list) or not markers:
            raise RuleSourceError(f"{phenotype}: omit_absent_anypos names no kept rule or no markers")
        tokens = dict(rule.constraints)
        for marker in markers:
            if tokens.get(marker) != "anypos":
                raise RuleSourceError(f"{phenotype}/{marker}: only an anypos clause may be left out")
            if marker in panel or aliases.get(marker) in panel:
                raise RuleSourceError(f"{phenotype}/{marker}: the panel measures it; it cannot be left out")
        remaining = tuple((m, t) for m, t in rule.constraints if m not in markers)
        if not any(t == "anypos" for _, t in remaining):
            raise RuleSourceError(f"{phenotype}: leaving out {markers} would empty its OR group")
        by_name[phenotype] = Rule(rule.parent, phenotype, remaining)
        omitted += [(phenotype, marker) for marker in markers]
    return [by_name[rule.phenotype] for rule in rules], tuple(omitted)


def _structure_only(adaptation: dict, source: RuleSource) -> frozenset[str]:
    names = adaptation.get("structure_only", [])
    roots = {rule.phenotype for rule in source.rules if rule.parent == source.lineage_root_parent}
    if not isinstance(names, list) or not set(names) <= roots:
        raise RuleSourceError("spleen structure_only must name kept lineage roots")
    return frozenset(names)


def structure_only_branches(manifest_path: Path) -> frozenset[str]:
    """The lineage roots that name a cell only inside a structure object of their class (maintainer,
    2026-11-01: the 10-30 structure-only ruling applied to spleen). Without such an object a match in the
    branch is unavailable, never negative (``positivity.hierarchy``)."""
    manifest = yaml.safe_load(Path(manifest_path).read_text())
    return _structure_only(manifest["rule_table"]["adaptation"], load_rule_source(manifest_path))


# --------------------------------------------------------------------------- per-image typing hooks
# The hooks multiplex5000.positivity reads (see tissues/pancreas.py for the CODEX originals). They place no cell:
# the spleen PC67 registry, which decides what the method may use, is a separate step reviewed by the maintainer.
#: populations a spleen registry may name as references: no endocrine or exocrine population, since those are
#: pancreas-only (downstream plan, 2026-10-30); epithelial cells away from leukocytes stay (epithelium is present)
TYPING_POPULATIONS = frozenset({"leukocyte", "epithelial_away", "structure_cells", "t_cells", "leukocyte_clean"})


def start_class(cells):
    """The cell-derived region the cell's centroid lies in (``start_object``: Follicle, Perifollicular, PALS, RedPulp;
    maintainer 2026-11-01), else Tissue. A region names no cell; it is where the cell starts."""
    import pandas as pd

    column = cells["start_object"] if "start_object" in cells else pd.Series([None] * len(cells), index=cells.index)
    return column.fillna("Tissue").to_numpy(dtype=object)


def reference_populations(registry) -> frozenset[str]:
    return TYPING_POPULATIONS


def selection_sections(roster) -> list[str]:
    """Every section of the run, blind to genotype: a method that fails on any section is revised for all of them."""
    return list(roster)


def structure_objects() -> tuple[str, ...]:
    """None yet: no working vessel object exists on PC67 (SmallVessel calls 29-78 % of tissue; SmallVessel_PC67 is
    being retrained). Endothelial, structure-only (manifest), stays unavailable until one does."""
    return ()


#: the check columns' name for the known-negative cells (``{marker}_in_reference``)
CHECK_POPULATION = "reference"


def check_population(cells, registry, roles: dict, pos: dict, background):
    """The known-negative cells of the leukocyte-marker checks (P7, maintainer 2026-11-01): spleen has no
    hormone-positive endocrine cells, so the checks read epithelial and stromal cells instead -- cells positive for an
    Epithelial majority marker or a defining marker of the Mesenchymal compartment. Each marker is read at its chosen
    level where one is chosen, otherwise at its own background peak (the immune stage runs before those levels are
    chosen). ``levels`` then keeps only those with no leukocyte within ISOLATION_UM, as it does hormone cells.

    No leukocyte-marker criterion enters the definition: a CD45-negative reference would pass the CD45 check by
    construction. Vimentin is not a stromal marker here: it is shared with leukocytes (lymphocytes carry it, and B
    cells are not leukocytes in the method's two-marker definition), so it would put B cells in the reference."""
    import numpy as np

    measured = set(registry.measured_markers())
    markers = [m for m in roles["epithelial"] if m in measured] + \
              [m for m in registry.defining("Mesenchymal") if m in measured]
    out = np.zeros(len(cells), bool)
    for m in markers:
        out |= pos[m] if m in pos else background(m)
    return out

