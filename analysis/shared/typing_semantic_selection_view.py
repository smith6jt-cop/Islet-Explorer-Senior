#!/usr/bin/env python3
"""Read an explicit schema2 semantic selection; never create or activate one.

Compact validation opens contracts and receipts only. Counts and iterators read
one explicitly requested scope. Protein scientific use requires the reviewed
semantic overlay and an explicit analytical purpose; raw QC claims are never
the fallback. See typing_semantic_selection_view_contract.md.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, field
import json
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import pyarrow.compute as pc

try:
    from . import protein_semantic_view as protein
    from . import typing_release_view as base
except ImportError:  # Explicit file CLI
    import protein_semantic_view as protein
    import typing_release_view as base


CAPABILITIES = frozenset({"protein_semantic_overlay_v1", "purpose_qualified_protein_consumption_v1"})
LANES = {**protein.PROTEIN_SCOPES, "pancreas_RNA": 25, "pLN_RNA": 27}
READBACK_STATUS = "TERMINAL_INDEPENDENT_FULL56_PROTEIN_SEMANTIC_OVERLAY_ALLROW_ALLFIELD_READBACK_PASS_UNSELECTED"
ACCEPTANCE_STATUS = "ROOT_ACCEPTED_INDEPENDENT_FULL56_PROTEIN_SEMANTIC_READBACK"
BATCH_SIZE = 8192
READBACK_GATES = frozenset({
    "all56_source_bytes_freshly_hashed_before_and_after",
    "all_original_columns_order_values_and_Arrow_metadata_exact",
    "all_derived_values_independently_reconstructed_and_exact",
    "all_native_QC_admissions_and_exclusions_preserved",
    "all649_QC_claim_parts_and609_ordered_UUID_hashes_exact",
    "all_QC_rows_have_exact_root_decision_and_inventory",
    "all40_T_NK_uncertainty_parts_preserved",
    "exact_PAN6523_scientific_hold_and_dedicated6523_distinguished",
    "all_purpose_eligible_ineligible_unavailable_denominators_reconciled",
    "no_source_candidate_or_selected_mutations",
})
require = protein.require


def freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({k: freeze(v) for k, v in value.items()})
    if isinstance(value, list):
        return tuple(freeze(v) for v in value)
    return value


def count(value, name):
    require(type(value) is int and value >= 0, name + " must be a nonnegative integer")
    return value


def checked(item, directory):
    return protein.checked_binding(item, directory)


def require_executed_api(binding, module_path, name):
    require(binding["path"] == str(Path(module_path).resolve()),
            "Bound " + name + " differs from the executing API file")


def validate_full56_readback(receipt, expected_bindings, scopes):
    """Validate accepted proof structure, never reread source cells here."""
    require(isinstance(receipt, dict) and type(receipt.get("schema_version")) is int and receipt["schema_version"] == 1 and
            receipt.get("status") == READBACK_STATUS, "A successful independent full56 semantic readback is required")
    readback_bindings = receipt.get("bindings")
    require(isinstance(readback_bindings, dict) and
            all(isinstance(readback_bindings.get(k), dict) and
                readback_bindings[k].get("path") == v["path"] and
                readback_bindings[k].get("sha256") == v["sha256"]
                for k, v in expected_bindings.items()),
            "Independent readback bindings differ from the semantic selection")
    require(receipt.get("across_scope_counts_are_additive") is False and
            receipt.get("biological_identity_accuracy_validated") is False,
            "Readback must retain nonadditive scope counts and no identity-accuracy claim")
    gates = receipt.get("selection_gate")
    require(isinstance(gates, dict) and READBACK_GATES.issubset(gates) and
            all(v is True for v in gates.values()),
            "Every required independent full56 operational selection gate must be present and true")
    require(receipt.get("expected_semantics_from_consumer_helpers") is False and
            type(receipt.get("n_derived_columns")) is int and receipt["n_derived_columns"] == 14,
            "Readback must independently reconstruct and verify all14 derived columns")
    coverage = receipt.get("coverage", {})
    require(coverage.get("n_pLN_units") == 609 and coverage.get("n_pLN_prior_claim_parts") == 649 and
            coverage.get("n_pLN_ordered_UUID_membership_hashes") == 609 and
            coverage.get("ordered_UUID_membership_hashes_verified") is True and
            coverage.get("exact_inventory_parts_verified") is True,
            "Independent readback must verify609 ordered UUID-membership hashes/units and649 prior-claim parts")
    rows = receipt.get("scopes")
    require(isinstance(rows, list) and len(rows) == 56 and receipt.get("n_scopes") == 56,
            "Independent readback must contain all56 protein scopes")
    seen = set()
    source_total = qc_total = 0
    for row in rows:
        require(isinstance(row, dict), "Readback scope must be an object")
        sid = row.get("scope_id")
        require(sid in scopes and sid not in seen, "Duplicate or foreign readback protein scope")
        seen.add(sid)
        s = scopes[sid]
        require(row.get("base_file") == {"path": str(Path(s["path"]).resolve()), "sha256": s["sha256"]},
                "Readback source file differs from selected scope: " + sid)
        n = count(row.get("source_rows"), "readback source_rows")
        require(n == s["expected_source_rows"] and
                count(row.get("native_QC_rows"), "readback native_QC_rows") == s["expected_QC_rows"],
                "Readback source/native QC denominator mismatch: " + sid)
        expected_admission = s.get("expected_additional_admission_rows") if s.get("additional_admission_column") else None
        observed_admission = row.get("native_additional_admission_rows")
        if expected_admission is not None:
            count(observed_admission, "readback native_additional_admission_rows")
        require(observed_admission == expected_admission,
                "Readback native additional-admission denominator mismatch: " + sid)
        require(row.get("original_row_order_preserved") is True and
                row.get("original_columns_values_preserved") is True and
                row.get("Arrow_schema_metadata_preserved") is True,
                "Independent readback must preserve original row order, every original column value and Arrow metadata")
        totals = row.get("purpose_totals")
        require(isinstance(totals, dict) and set(totals) == set(protein.PURPOSES),
                "Independent readback requires all four explicit purpose totals")
        for purpose, values in totals.items():
            require(isinstance(values, dict), "Purpose totals must be an object")
            e, f, u = [count(values.get(k), purpose + " " + k) for k in
                       ("eligible_rows", "ineligible_rows", "unavailable_rows")]
            require(e + f + u == n, "Purpose totals do not partition retained source rows")
            if purpose == "native_provenance":
                require((e, f, u) == (n, 0, 0), "Native provenance must retain every source row")
            else:
                require(e <= s["expected_QC_rows"] and
                        (expected_admission is None or e <= expected_admission),
                        "Purpose eligibility bypasses native QC or native additional admission")
                require(sid != protein.HELD_SCOPE or e == 0,
                        "Held ordinary-PAN6523 source cannot admit a scientific purpose")
        source_total += n
        qc_total += s["expected_QC_rows"]
    require(seen == set(scopes), "Independent readback omits a protein scope")
    require(count(receipt.get("n_source_rows_nonadditive"), "n_source_rows_nonadditive") == source_total,
            "Independent readback nonadditive source row total mismatch")
    require(count(receipt.get("n_native_QC_observations_nonadditive"), "n_native_QC_observations_nonadditive") == qc_total,
            "Independent readback nonadditive native-QC observation total mismatch")


@dataclass(frozen=True)
class ScopeMetadata:
    """Metadata snapshot; deliberately not a legacy scope Mapping."""
    scope_id: str
    lane: str
    fields: Mapping


@dataclass(frozen=True)
class SemanticSelection:
    _scopes: Mapping = field(repr=False)
    _overlay: protein.ProteinSemanticOverlay = field(repr=False)
    bindings: Mapping

    def _scope(self, scope_id):
        require(scope_id in self._scopes, "Unknown semantic-selection scope")
        return self._scopes[scope_id]

    def scope_metadata(self, scope_id):
        s = self._scope(scope_id)
        return ScopeMetadata(scope_id, s["lane"], freeze(s))

    def summary(self):
        return {"status": "SCHEMA2_SEMANTIC_SELECTION_CONTRACT_VALIDATED_NOT_ACTIVATED",
                "n_scopes": len(self._scopes), "scope_counts": dict(LANES),
                "bindings": dict(self.bindings), "base_cell_files_read_or_hashed": False,
                "protein_default_known_lineage_counts_available": False,
                "protein_requires_explicit_purpose": True, "RNA_scope_contracts_unchanged": True,
                "across_scope_counts_are_additive": False, "selection_activated": False,
                "biological_identity_accuracy_validated": False}

    def counts(self, scope_id, *, purpose=None, verify_file=True, batch_size=BATCH_SIZE):
        s = self._scope(scope_id)
        if s["lane"] in protein.PROTEIN_SCOPES:
            require(purpose in protein.PURPOSES, "Protein counts require an explicit supported purpose")
            return protein.scope_counts(self._overlay, scope_id, purpose=purpose,
                                        verify_file=verify_file, batch_size=batch_size)
        require(purpose is None, "RNA counts retain the unchanged selected RNA contract; omit protein purpose")
        return base.scope_counts(s, verify_file=verify_file)

    def iter_annotations(self, scope_id, *, columns=None, verify_file=True, batch_size=BATCH_SIZE):
        s = self._scope(scope_id)
        if s["lane"] in protein.PROTEIN_SCOPES:
            yield from protein.iter_resolved_batches(self._overlay, scope_id, columns=columns,
                                                      verify_file=verify_file, batch_size=batch_size)
        else:
            yield from base.iter_scope_batches(s, columns=columns, verify_file=verify_file, batch_size=batch_size)

    def iter_for_purpose(self, scope_id, *, purpose, columns=None, verify_file=True, batch_size=BATCH_SIZE):
        require(self._scope(scope_id)["lane"] in protein.PROTEIN_SCOPES,
                "Purpose filtering is protein-only; use unchanged RNA annotations")
        require(purpose in protein.PURPOSES, "Protein filtering requires an explicit supported purpose")
        for table in self.iter_annotations(scope_id, columns=columns, verify_file=verify_file, batch_size=batch_size):
            # Missing remains missing in annotation/count views. This explicit use
            # iterator selects true only and never turns missing into a negative call.
            selected = table.filter(pc.fill_null(table[f"eligible_for_{purpose}"], False))
            if len(selected):
                yield selected


def load_semantic_selection(path):
    """Validate the explicit schema2 entry and compact independent acceptance proof."""
    path = Path(path).resolve()
    entry = json.loads(path.read_text())
    require(type(entry.get("schema_version")) is int and entry["schema_version"] == 2 and
            entry.get("kind") == "typing_semantic_selection", "An explicit schema2 semantic selection is required")
    capabilities = entry.get("required_capabilities")
    require(isinstance(capabilities, list) and len(capabilities) == len(CAPABILITIES) and
            set(capabilities) == CAPABILITIES, "Missing or unsupported required semantic consumer capability")
    require(entry.get("across_scope_counts_are_additive") is False and
            entry.get("biological_identity_accuracy_validated") is False,
            "Semantic selection must retain nonadditivity and no identity-accuracy claim")
    directory = path.parent
    bindings = {k: checked(entry.get(k), directory) for k in
                ("base_registry", "previous_selection_inputs", "consumer", "root_readback_acceptance")}
    require_executed_api(bindings["consumer"], __file__, "semantic-selection consumer")
    registry = base.load_registry(bindings["base_registry"]["path"], verify_bindings=False)
    scopes = entry.get("scopes")
    require(isinstance(scopes, list) and scopes == registry["scopes"],
            "All108 RNA/protein scope dictionaries and ordering must equal the unchanged base registry")
    require(Counter(s["lane"] for s in scopes) == Counter(LANES), "Expected52 RNA and56 protein scopes")
    previous = json.loads(Path(bindings["previous_selection_inputs"]["path"]).read_text())
    require(previous.get("registry") == bindings["base_registry"], "Previous selected inputs bind another registry")
    generations = registry.get("current_RNA_lane_generations")
    require(isinstance(generations, dict) and entry.get("current_RNA_lane_generations") == generations ==
            previous.get("current_RNA_lane_generations"), "Selected RNA generation metadata changed")
    ps, rs = entry.get("protein_semantics"), entry.get("RNA_semantics")
    require(isinstance(ps, dict) and isinstance(rs, dict), "Explicit protein and unchanged RNA interfaces are required")
    require(ps.get("required_scientific_purposes") == list(protein.SCIENTIFIC_PURPOSES) and
            ps.get("explicit_native_provenance_purpose") == "native_provenance",
            "Protein scientific purposes and separate native provenance must remain explicit")
    require(rs.get("policy") == "unchanged_selected_RNA_scope_contracts", "RNA semantics must remain unchanged")
    bindings.update({k: checked(ps.get(k), directory) for k in
                    ("overlay", "semantic_api", "root_acceptance", "independent_full_row_readback")})
    bindings["base_api"] = checked(rs.get("base_api"), directory)
    require_executed_api(bindings["semantic_api"], protein.__file__, "protein semantic API")
    require_executed_api(bindings["base_api"], base.__file__, "native base API")
    overlay = protein.load_overlay(bindings["overlay"]["path"])
    require(overlay.bindings["registry"] == bindings["base_registry"] and
            overlay.bindings.get("root_acceptance") == bindings["root_acceptance"],
            "Overlay registry or scientific acceptance differs from the semantic selection")
    require(all(overlay.pln_decisions.get(k) is not None for k in overlay.pln_inventory) and
            all(d is not None for d in overlay.pan_routes.values()),
            "A final selection requires all649 pLN root decisions and explicit PAN dictionary decisions")
    expected = {"registry": bindings["base_registry"], "overlay_contract": bindings["overlay"],
                "actual_consumer": bindings["semantic_api"],
                "actual_consumer_base_helper": bindings["base_api"],
                "root_acceptance": bindings["root_acceptance"]}
    receipt = json.loads(Path(bindings["independent_full_row_readback"]["path"]).read_text())
    protein_scopes = {s["scope_id"]: s for s in scopes if s["lane"] in protein.PROTEIN_SCOPES}
    validate_full56_readback(receipt, expected, protein_scopes)
    acceptance = json.loads(Path(bindings["root_readback_acceptance"]["path"]).read_text())
    require(isinstance(acceptance, dict) and acceptance.get("status") == ACCEPTANCE_STATUS and
            acceptance.get("full56_semantics_accepted") is True,
            "Root must explicitly accept independent full56 semantics, not mere completion")
    require(acceptance.get("bindings") == {**expected, "independent_full_row_readback": bindings["independent_full_row_readback"]},
            "Root readback acceptance does not endorse these exact semantic/readback bindings")
    bindings["semantic_selection"] = {"path": str(path), "sha256": base.sha256(path)}
    return SemanticSelection(MappingProxyType({s["scope_id"]: s for s in scopes}), overlay, freeze(bindings))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, required=True)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate", help="Compact contracts/receipts only; no source cell reads")
    sub.add_parser("inventory", help="Scope metadata only")
    counts = sub.add_parser("counts", help="Explicit one-scope source readback")
    counts.add_argument("--scope-id", required=True)
    counts.add_argument("--purpose", choices=protein.PURPOSES)
    counts.add_argument("--output", type=Path, help="Optional NEW receipt; refuses overwrite")
    args = parser.parse_args()
    selection = load_semantic_selection(args.selection)
    if args.command == "counts":
        result = selection.counts(args.scope_id, purpose=args.purpose)
    else:
        result = selection.summary()
        if args.command == "inventory":
            result["scopes"] = [{k: s[k] for k in ("scope_id", "lane", "expected_source_rows", "expected_QC_rows")}
                                for s in selection._scopes.values()]
    # Nested immutable bindings are MappingProxyType; serialize their JSON snapshot.
    text = json.dumps(result, indent=2, allow_nan=False, default=lambda x: dict(x) if isinstance(x, Mapping) else x)+"\n"
    if args.command == "counts" and args.output:
        with args.output.open("x") as stream:
            stream.write(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
