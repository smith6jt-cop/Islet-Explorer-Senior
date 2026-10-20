#!/usr/bin/env python3
"""Resolve explicit protein semantics without modifying a typing generation.

``validate`` reads only the compact overlay contracts. ``counts`` explicitly
streams one selected base scope. Neither command creates labels or activates a
release. See protein_semantic_view_contract.md for the JSON contract.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import pyarrow as pa

try:
    from . import typing_release_view as base_view
except ImportError:  # Explicit CLI: python analysis/shared/protein_semantic_view.py
    import typing_release_view as base_view


PROTEIN_SCOPES = {"pancreas_protein": 25, "dedicated_pLN_protein": 19,
                  "embedded_pLN_protein": 12}
HELD_SCOPE = "pancreas_protein|6523|ordinary_pancreas"
SCIENTIFIC_PURPOSES = ("programme", "broad_lineage", "structure")
PURPOSES = (*SCIENTIFIC_PURPOSES, "native_provenance")
EVIDENCE_LEVELS = {
    "hormone_informed_native_allocation", "affirmative_broad_immune_allocation",
    "qualified_broad_lineage_programme", "structure_associated_allocation",
    "structural_or_contextual_pattern", "identity_unresolved", "analysis_status",
}
NON_LINEAGE_LEVELS = {"structure_associated_allocation", "structural_or_contextual_pattern",
                      "identity_unresolved", "analysis_status"}
KEY_FIELDS = ("scope_id", "cluster_id", "prior_candidate_lineage")
DERIVED_COLUMNS = ("annotation_name", "evidence_level", "annotation_reason",
                   "semantic_decision_status", "source_use_status", "source_use_reason",
                   *(f"eligible_for_{p}" for p in PURPOSES),
                   *(f"eligibility_reason_for_{p}" for p in PURPOSES))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def nonempty(value, name):
    require(isinstance(value, str) and bool(value.strip()), name + " must be a nonempty string")
    return value


def checked_binding(item, directory):
    require(isinstance(item, dict), "A file binding must be an object")
    path = Path(nonempty(item.get("path"), "binding path"))
    if not path.is_absolute():
        path = directory / path
    path = path.resolve()
    digest = nonempty(item.get("sha256"), "binding sha256")
    require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
            "Binding SHA256 must be lowercase hexadecimal")
    require(base_view.sha256(path) == digest, "Changed semantic input: " + str(path))
    return {"path": str(path), "sha256": digest}


def records(binding):
    path = Path(binding["path"])
    if path.suffix == ".csv":
        with path.open(newline="") as stream:
            return list(csv.DictReader(stream))
    value = json.loads(path.read_text())
    require(isinstance(value, list) and all(isinstance(r, dict) for r in value),
            "Semantic record file must contain a JSON list of objects or CSV records")
    return value


def part_key(record):
    return tuple(nonempty(record.get(k), k) for k in KEY_FIELDS)


@dataclass(frozen=True)
class Decision:
    annotation_name: str
    evidence_level: str
    reason: str
    eligibility: tuple[bool | None, ...]


def decision(value):
    if value is None:
        return None
    require(isinstance(value, dict), "A root decision must be an object or null")
    level = nonempty(value.get("evidence_level"), "evidence_level")
    require(level in EVIDENCE_LEVELS, "Unknown evidence_level: " + level)
    eligibility = value.get("eligibility")
    require(isinstance(eligibility, dict) and set(eligibility) == set(SCIENTIFIC_PURPOSES),
            "Root eligibility must explicitly specify programme, broad_lineage and structure")
    values = tuple(eligibility[p] for p in SCIENTIFIC_PURPOSES)
    require(all(v is None or type(v) is bool for v in values),
            "Root eligibility must use JSON true, false or null; numbers/strings are invalid")
    require(not (level in NON_LINEAGE_LEVELS and eligibility["broad_lineage"] is True),
            "Pattern, structure, unresolved and analysis-status decisions cannot admit a named lineage")
    return Decision(nonempty(value.get("annotation_name"), "annotation_name"), level,
                    nonempty(value.get("reason"), "decision reason"), values)


@dataclass(frozen=True)
class ProteinSemanticOverlay:
    scopes: Mapping
    pan_routes: Mapping[str, Decision | None]
    pln_inventory: Mapping
    pln_decisions: Mapping
    source_hold_reason: str
    bindings: Mapping

    def summary(self):
        return {
            "status": "compact_semantic_contract_validated_not_activated",
            "n_base_scopes": len(self.scopes), "scope_counts": dict(PROTEIN_SCOPES),
            "n_pLN_inventory_parts": len(self.pln_inventory),
            "n_pLN_inventory_units": len({k[:2] for k in self.pln_inventory}),
            "n_pLN_inventory_scopes": len({k[0] for k in self.pln_inventory}),
            "n_pLN_root_decisions": sum(d is not None for d in self.pln_decisions.values()),
            "n_pLN_root_decisions_unavailable": sum(
                self.pln_decisions.get(k) is None for k in self.pln_inventory),
            "n_PAN_dictionary_routes": len(self.pan_routes),
            "n_PAN_dictionary_decisions_unavailable": sum(d is None for d in self.pan_routes.values()),
            "held_scientific_source_scope": HELD_SCOPE,
            "bindings": dict(self.bindings), "base_cell_files_read_or_hashed": False,
            "across_scope_counts_are_additive": False, "selection_activated": False,
            "biological_identity_accuracy_validated": False,
        }


def load_overlay(path) -> ProteinSemanticOverlay:
    """Validate compact pins/keys only. Base cell files are not opened or hashed."""
    path = Path(path).resolve()
    contract = json.loads(path.read_text())
    require(type(contract.get("schema_version")) is int and contract["schema_version"] == 1,
            "Unsupported protein semantic contract schema")
    require(contract.get("kind") == "protein_semantic_overlay", "Wrong semantic contract kind")
    require(contract.get("across_scope_counts_are_additive") is False,
            "Protein semantic views must explicitly remain nonadditive")
    directory = path.parent
    bindings = {k: checked_binding(contract.get(k), directory) for k in
                ("registry", "PAN_route_dictionary", "pLN_unit_claim_inventory", "pLN_root_decisions")}
    registry = base_view.load_registry(bindings["registry"]["path"], verify_bindings=False)
    scopes = {s["scope_id"]: s for s in registry["scopes"] if s["lane"] in PROTEIN_SCOPES}
    require(Counter(s["lane"] for s in scopes.values()) == Counter(PROTEIN_SCOPES),
            "Overlay requires exactly the selected 25 PAN,19 dedicated and12 embedded protein scopes")
    bases = contract.get("base_files")
    require(isinstance(bases, list) and len(bases) == 56, "All56 base scope pins are required")
    pinned = {}
    for item in bases:
        require(isinstance(item, dict), "A base scope pin must be an object")
        sid = nonempty(item.get("scope_id"), "base scope_id")
        require(sid not in pinned and sid in scopes, "Duplicate or foreign base scope: " + sid)
        source = Path(nonempty(item.get("path"), "base path"))
        if not source.is_absolute():
            source = directory/source
        s = scopes[sid]
        require(str(source.resolve()) == str(Path(s["path"]).resolve()) and item.get("sha256") == s["sha256"],
                "Base pin differs from exact selected scope: " + sid)
        require(s["QC_kind"] == "boolean", "Protein scopes require an explicit boolean native QC")
        expected_claim = "candidate_cell_class" if s["lane"] == "pancreas_protein" else "candidate_lineage"
        require(s["claim_column"] == expected_claim, "Unexpected selected protein claim column")
        pinned[sid] = s
    require(set(pinned) == set(scopes) and HELD_SCOPE in scopes, "Incomplete selected protein base binding")
    hold = contract.get("source_hold")
    require(isinstance(hold, dict) and hold.get("scope_id") == HELD_SCOPE,
            "Only the exact ordinary-PAN6523 acquisition may receive this source hold")
    hold_reason = nonempty(hold.get("reason"), "source hold reason")
    pan = {}
    for r in records(bindings["PAN_route_dictionary"]):
        label = nonempty(r.get("prior_candidate_cell_class"), "prior_candidate_cell_class")
        require(label not in pan, "Duplicate PAN route dictionary key")
        pan[label] = decision(r.get("decision"))
    require(bool(pan), "An explicit PAN route dictionary is required")
    inventory = {}
    for r in records(bindings["pLN_unit_claim_inventory"]):
        key = part_key(r)
        require(key not in inventory, "Duplicate pLN unit/prior-claim inventory key")
        require(key[0] in scopes and scopes[key[0]]["lane"] != "pancreas_protein",
                "Foreign or PAN scope in pLN inventory")
        n = r.get("current_native_QC_claim_rows")
        require(type(n) is int or (isinstance(n, str) and n.isdecimal()),
                "Inventory requires an explicit native-QC part count")
        require(int(n) > 0, "Inventory part count must be positive")
        inventory[key] = int(n)
    require(len(inventory) == 649 and len({k[0] for k in inventory}) == 31,
            "Current pLN inventory requires exactly649 parts in31 selected scopes")
    require(len({k[:2] for k in inventory}) == 609,
            "Current pLN inventory requires exactly609 unique scope/cluster units")
    inventory_qc = Counter()
    for key, n in inventory.items():
        inventory_qc[key[0]] += n
    require(all(inventory_qc[sid] == s["expected_QC_rows"] for sid, s in scopes.items()
                if s["lane"] != "pancreas_protein"),
            "pLN inventory native-QC part counts differ from selected per-scope denominators")
    pln = {}
    for r in records(bindings["pLN_root_decisions"]):
        key = part_key(r)
        require(key in inventory and key not in pln, "Duplicate or foreign pLN root decision key")
        pln[key] = decision(r.get("decision"))
    has_decisions = any(d is not None for d in (*pan.values(), *pln.values()))
    if has_decisions:
        bindings["root_acceptance"] = checked_binding(contract.get("root_acceptance"), directory)
    elif contract.get("root_acceptance") is not None:
        bindings["root_acceptance"] = checked_binding(contract["root_acceptance"], directory)
    if "root_acceptance" in bindings:
        acceptance = json.loads(Path(bindings["root_acceptance"]["path"]).read_text())
        require(isinstance(acceptance, dict) and
                acceptance.get("status") == "ROOT_ACCEPTED_PROTEIN_SEMANTIC_DECISIONS",
                "Root acceptance must explicitly accept the protein semantic decisions")
        accepted = acceptance.get("bindings")
        require(isinstance(accepted, dict) and all(accepted.get(k) == bindings[k] for k in
                ("registry", "PAN_route_dictionary", "pLN_unit_claim_inventory", "pLN_root_decisions")),
                "Root acceptance does not endorse the exact checked semantic input bindings")
        require(acceptance.get("source_hold") == hold,
                "Root acceptance source hold differs from the exact overlay source hold")
    bindings["overlay_contract"] = {"path": str(path), "sha256": base_view.sha256(path)}
    return ProteinSemanticOverlay(MappingProxyType(scopes), MappingProxyType(pan),
        MappingProxyType(inventory), MappingProxyType(pln), hold_reason, MappingProxyType(bindings))


def effective_eligibility(d, purpose, qc, admission, held):
    if purpose == "native_provenance":
        return True, "Source provenance includes every retained source row; no identity admission implied"
    if d is None:
        return None, "Root semantic decision unavailable"
    explicit = d.eligibility[SCIENTIFIC_PURPOSES.index(purpose)]
    if explicit is None:
        return None, "Root purpose eligibility unavailable"
    if held:
        return False, "Exact ordinary-PAN6523 source identity held; native QC/admission preserved"
    if qc is False or admission is False:
        return False, "Native QC or native additional admission is ineligible"
    if qc is None or admission is None:
        return None, "Native QC or native additional admission is unavailable"
    return explicit, d.reason


def resolve_batch(overlay: ProteinSemanticOverlay, scope_id: str, table: pa.Table) -> pa.Table:
    """Append nullable semantic fields in original row order; never overwrite native fields."""
    require(scope_id in overlay.scopes, "Unknown selected protein scope")
    s = overlay.scopes[scope_id]
    require(not set(DERIVED_COLUMNS).intersection(table.column_names),
            "Refusing to overwrite existing semantic columns")
    required = {s["cell_key_column"], s["claim_column"], s["QC_column"]}
    additional = s.get("additional_admission_column")
    if additional:
        required.add(additional)
    is_pan = s["lane"] == "pancreas_protein"
    if not is_pan:
        required.add("cluster_id")
    require(required.issubset(table.column_names), "Base batch lacks contracted semantic lookup columns")
    require(pa.types.is_boolean(table[s["QC_column"]].type), "Native QC must remain boolean")
    if additional:
        require(pa.types.is_boolean(table[additional].type), "Native additional admission must remain boolean")
    qc = table[s["QC_column"]].to_pylist()
    admitted = table[additional].to_pylist() if additional else [True]*len(table)
    require(not any(a is True and q is not True for q, a in zip(qc, admitted)) if additional else True,
            "Native additional admission cannot bypass native QC")
    claims = table[s["claim_column"]].to_pylist()
    clusters = [None]*len(table) if is_pan else table["cluster_id"].to_pylist()
    held = scope_id == HELD_SCOPE
    values = {name: [] for name in DERIVED_COLUMNS}
    for q, a, label, cluster in zip(qc, admitted, claims, clusters, strict=True):
        d = overlay.pan_routes.get(label) if is_pan else overlay.pln_decisions.get((scope_id, cluster, label))
        values["annotation_name"].append(d.annotation_name if d else None)
        values["evidence_level"].append(d.evidence_level if d else None)
        values["annotation_reason"].append(d.reason if d else None)
        values["semantic_decision_status"].append("root_decision_present" if d else "root_decision_unavailable")
        values["source_use_status"].append("held_source_identity" if held else "source_bound")
        values["source_use_reason"].append(overlay.source_hold_reason if held else
            "Exact selected source binding; independent biological identity truth is not implied")
        for purpose in PURPOSES:
            eligible, reason = effective_eligibility(d, purpose, q, a, held)
            values[f"eligible_for_{purpose}"].append(eligible)
            values[f"eligibility_reason_for_{purpose}"].append(reason)
    for name, data in values.items():
        table = table.append_column(name, pa.array(data, type=pa.bool_() if name.startswith("eligible_for_") else pa.string()))
    return table


def iter_resolved_batches(overlay, scope_id, *, columns=None, verify_file=True, batch_size=131072):
    """Explicitly read one immutable selected base; verify its bytes by default."""
    require(scope_id in overlay.scopes, "Unknown selected protein scope")
    s = overlay.scopes[scope_id]
    if columns is not None:
        columns = set(columns)
        if s["lane"] != "pancreas_protein":
            columns.add("cluster_id")
    for table in base_view.iter_scope_batches(s, columns=columns, verify_file=verify_file, batch_size=batch_size):
        yield resolve_batch(overlay, scope_id, table)


def scope_counts(overlay, scope_id, *, purpose, verify_file=True, batch_size=131072):
    """Count only explicit eligibility for a requested purpose; retain false and unavailable."""
    require(purpose in PURPOSES, "Unknown analytical purpose")
    require(scope_id in overlay.scopes, "Unknown selected protein scope")
    s = overlay.scopes[scope_id]
    counts = Counter()
    annotations, joint = Counter(), Counter()
    extra = s.get("additional_admission_column")
    cols = [s["claim_column"]]
    for table in iter_resolved_batches(overlay, scope_id, columns=cols, verify_file=verify_file, batch_size=batch_size):
        n = len(table)
        qc = table[s["QC_column"]].to_pylist()
        admission = table[extra].to_pylist() if extra else [None]*n
        eligible = table[f"eligible_for_{purpose}"].to_pylist()
        names = table["annotation_name"].to_pylist()
        levels = table["evidence_level"].to_pylist()
        status = table["semantic_decision_status"].to_pylist()
        sources = table["source_use_status"].to_pylist()
        for q, a, e, name, level, state, source in zip(qc, admission, eligible, names, levels, status, sources, strict=True):
            counts["source_rows"] += 1
            counts["native_QC_rows"] += q is True
            counts["native_QC_unavailable_rows"] += q is None
            counts["native_additional_admission_rows"] += a is True
            counts["native_additional_admission_unavailable_rows"] += a is None and extra is not None
            counts["held_source_rows"] += source == "held_source_identity"
            counts["root_decision_unavailable_rows"] += state == "root_decision_unavailable"
            counts[{True: "purpose_eligible_rows", False: "purpose_ineligible_rows", None: "purpose_unavailable_rows"}[e]] += 1
            joint[(q, a, source, state, e)] += 1
            if e is True:
                annotations[(name, level)] += 1
    require(counts["source_rows"] == s["expected_source_rows"] and
            counts["native_QC_rows"] == s["expected_QC_rows"], "Native scope denominator mismatch")
    if extra:
        require(counts["native_additional_admission_rows"] == s["expected_additional_admission_rows"],
                "Native additional-admission denominator mismatch")
    result = {k: counts[k] for k in ("source_rows", "native_QC_rows", "native_QC_unavailable_rows",
        "held_source_rows", "root_decision_unavailable_rows", "purpose_eligible_rows",
        "purpose_ineligible_rows", "purpose_unavailable_rows")}
    result.update(scope_id=scope_id, purpose=purpose,
        native_additional_admission_rows=counts["native_additional_admission_rows"] if extra else None,
        native_additional_admission_unavailable_rows=counts["native_additional_admission_unavailable_rows"] if extra else None,
        purpose_annotation_counts=[{"annotation_name": name, "evidence_level": level, "n": n}
            for (name, level), n in sorted(annotations.items(), key=lambda r: str(r[0]))],
        joint_status_counts=[{"native_QC": q, "native_additional_admission": a,
            "source_use_status": source, "semantic_decision_status": state, "purpose_eligible": e, "n": n}
            for (q, a, source, state, e), n in sorted(joint.items(), key=lambda r: str(r[0]))],
        separate_native_additional_admission_declared=bool(extra),
        across_scope_counts_are_additive=False, biological_identity_accuracy_validated=False,
        source_or_selection_mutations=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overlay", type=Path, required=True)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate", help="Check only compact contracts; never read base cell files")
    counts = sub.add_parser("counts", help="Explicitly stream and hash one selected protein base scope")
    counts.add_argument("--scope-id", required=True)
    counts.add_argument("--purpose", choices=PURPOSES, required=True)
    counts.add_argument("--output", type=Path, help="Optional NEW count-receipt file; never overwrites")
    args = parser.parse_args()
    overlay = load_overlay(args.overlay)
    result = overlay.summary() if args.command == "validate" else scope_counts(
        overlay, args.scope_id, purpose=args.purpose)
    text = json.dumps(result, indent=2, allow_nan=False)+"\n"
    if args.command == "counts" and args.output:
        with args.output.open("x") as stream:
            stream.write(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
