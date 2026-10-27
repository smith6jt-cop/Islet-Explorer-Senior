"""Load and audit the selected pLN protein rules without making cell calls.

The full sheet is reference knowledge for one hierarchical method. Xenium-only
cells need RNA evidence; same-slide cells need both assays and recorded links.
An absent protein marker is never a negative observation.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import yaml


class RuleSourceError(ValueError):
    """The recorded rules or adaptation cannot be used as declared."""


@dataclass(frozen=True)
class Rule:
    parent: str
    phenotype: str
    constraints: tuple[tuple[str, str], ...]

    @property
    def markers(self) -> frozenset[str]:
        return frozenset(marker for marker, _ in self.constraints)


@dataclass(frozen=True)
class RuleSource:
    sha256: str
    rules: tuple[Rule, ...]
    omitted: tuple[tuple[str, str], ...]
    panel_aliases: dict[str, str]
    lineage_root_parent: str = "all"
    state_root_parent: str = "any"
    candidate_lineage_roots: tuple[str, ...] = ()

    @property
    def by_phenotype(self) -> dict[str, Rule]:
        return {rule.phenotype: rule for rule in self.rules}

    def inherited_markers(self, phenotype: str) -> frozenset[str]:
        rules = self.by_phenotype
        if phenotype not in rules:
            raise KeyError(phenotype)
        current, seen, markers = phenotype, set(), set()
        while current in rules:
            if current in seen:
                raise RuleSourceError(f"cycle in rule parents at {current}")
            seen.add(current)
            rule = rules[current]
            markers.update(rule.markers)
            current = rule.parent
        return frozenset(markers)

    def missing_by_phenotype(self, measured_markers: Iterable[str]) -> dict[str, tuple[str, ...]]:
        """Return missing inherited markers; gaps are availability, not cell calls."""
        measured = set(measured_markers)
        return {
            rule.phenotype: tuple(sorted(
                marker for marker in self.inherited_markers(rule.phenotype)
                if marker not in measured and self.panel_aliases.get(marker) not in measured
            ))
            for rule in self.rules
        }


def load_rule_source(manifest_path: Path, *, table_path: Path | None = None) -> RuleSource:
    """Verify the source and apply only the maintainer-approved pLN adaptation."""
    manifest_path = Path(manifest_path)
    manifest = yaml.safe_load(manifest_path.read_text())
    if not isinstance(manifest, dict) or manifest.get("tissue_id") != "pancreatic_lymph_node":
        raise RuleSourceError("expected the pancreatic_lymph_node panel manifest")
    spec = manifest.get("rule_table")
    if not isinstance(spec, dict):
        raise RuleSourceError("manifest has no rule_table")
    if spec.get("assay") != "phenocycler":
        raise RuleSourceError("the pLN protein rules require the phenocycler assay")
    table_path = Path(table_path) if table_path else manifest_path.parent / spec["file"]
    digest = hashlib.sha256(table_path.read_bytes()).hexdigest()
    if digest != spec["source_sha256"]:
        raise RuleSourceError(f"{table_path}: digest differs from the selected source")
    with table_path.open(newline="") as source:
        reader = csv.DictReader(source)
        headers = reader.fieldnames or []
        if headers[:2] != ["CellType", "Phenotype"] or len(headers) != len(set(headers)):
            raise RuleSourceError("expected unique CellType, Phenotype, then marker columns")
        rows = list(reader)
    if len(rows) != spec["rows"]:
        raise RuleSourceError(f"expected {spec['rows']} rule rows, found {len(rows)}")
    allowed = set(spec["tokens"])
    root_roles = spec.get("root_roles")
    if (not isinstance(root_roles, dict) or set(root_roles) != {"lineage", "state"}
            or set(root_roles.values()) != {"all", "any"}):
        raise RuleSourceError("manifest must distinguish lineage and state root parents")
    omit = set(spec["pln_adaptation"]["omit_from_root_allneg"])
    if not omit <= set(headers[2:]):
        raise RuleSourceError(f"adaptation names unknown marker(s): {sorted(omit-set(headers))}")
    aliases = dict(spec.get("panel_aliases") or {})
    if not set(aliases) <= set(headers[2:]):
        raise RuleSourceError("panel_aliases names a marker absent from the sheet")
    rules, omitted, names = [], [], set()
    for row in rows:
        parent, phenotype = row["CellType"].strip(), row["Phenotype"].strip()
        if not phenotype or phenotype in names:
            raise RuleSourceError(f"blank or duplicate phenotype: {phenotype!r}")
        names.add(phenotype)
        clauses = []
        for marker in headers[2:]:
            token = (row[marker] or "").strip()
            if token and token not in allowed:
                raise RuleSourceError(f"{phenotype}/{marker}: unsupported token {token!r}")
            if parent == "all" and marker in omit and token == "allneg":
                omitted.append((phenotype, marker))
            elif token:
                clauses.append((marker, token))
        rules.append(Rule(parent, phenotype, tuple(clauses)))
    for rule in rules:
        if rule.parent not in {"all", "any"} and rule.parent not in names:
            raise RuleSourceError(f"{rule.phenotype}: missing parent {rule.parent}")
    candidate_roots = spec.get("candidate_lineage_roots")
    if (not isinstance(candidate_roots, list) or not candidate_roots
            or not all(isinstance(name, str) for name in candidate_roots)
            or len(candidate_roots) != len(set(candidate_roots))):
        raise RuleSourceError("manifest needs distinct candidate lineage roots")
    available_roots = {rule.phenotype for rule in rules
                       if rule.parent == root_roles["lineage"]}
    if not set(candidate_roots) <= available_roots:
        raise RuleSourceError("candidate lineage root is absent from the selected sheet")
    result = RuleSource(digest, tuple(rules), tuple(omitted), aliases,
                        root_roles["lineage"], root_roles["state"],
                        tuple(candidate_roots))
    for rule in rules:
        result.inherited_markers(rule.phenotype)
    if {marker for _, marker in omitted} != omit:
        raise RuleSourceError("an approved omitted marker was not a root allneg constraint")
    return result
