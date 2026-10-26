"""Spleen protein reference derived from the selected 6551 hierarchical sheet.

This loads marker knowledge and reports unavailable rules. It does not assign
cell types before normalization, controls and the Stage 06 gate are reviewed.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

from .lymph_node import Rule, RuleSource, RuleSourceError
from .lymph_node import load_rule_source as load_pln_rule_source


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
    aliases = {**source.panel_aliases, **manifest["protein_panel"].get("aliases", {})}
    result = RuleSource(source.sha256, tuple(kept), source.omitted, aliases)
    for rule in result.rules:
        result.inherited_markers(rule.phenotype)
    return result
