#!/usr/bin/env python3
"""Consume explicitly reviewed exact-cell annotations in a reversible candidate view.

The accepted schema2 selection and every original column remain unchanged. This
opt-in view appends effective candidate annotations and broad-lineage eligibility;
it cannot establish truth, train/tune labels, activate an alias, or lift a hold.
Only the two reviewed source scopes are supported; use the unchanged base API
elsewhere. Exhaust either iterator completely to validate all reviewed key
membership, including rows in batches that produce no purpose-filtered output.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass, field
import json
from pathlib import Path
from types import MappingProxyType

import pyarrow as pa
import pyarrow.compute as pc

try:
    from . import typing_semantic_selection_view as semantic
except ImportError:
    import typing_semantic_selection_view as semantic

require = semantic.require
DERIVED = {
    'effective_candidate_annotation_name': pa.string(),
    'effective_candidate_broad_lineage_eligible': pa.bool_(),
    'effective_candidate_identity_status': pa.string(),
    'effective_candidate_identity_uncertainty': pa.bool_(),
    'exact_cell_review_applied': pa.bool_(),
    'exact_cell_review_reason': pa.string(),
    'exact_cell_review_evidence_sha256': pa.string(),
    'exact_cell_review_subtype': pa.string(),
    'exact_cell_review_function': pa.string(),
    'reviewed_annotation_truth_training_tuning_eligible': pa.bool_(),
}
PURPOSE = 'conditional_candidate_broad_lineage'
POLICY = {
    'identity_status': 'conditional_source_bound_exact_cell_review',
    'independent_truth': False, 'training_or_tuning': False,
    'production_alias_activation': False, 'release_latch_change': False,
    'subtype_or_function_assignment': False, 'across_scope_counts_are_additive': False,
    'review_provenance': 'separate_agent_computational_case_review_not_blinded_expert_truth',
    'RNA_nucleus_requirement': 'native_count_exactly_one',
    'protein_nucleus_requirement': 'positive_explicit_linked_native_geometry_count_not_established',
}
REVIEWED_SCOPES = frozenset({
    'pLN_RNA|0037502|6380',
    'embedded_pLN_protein|6580|ee6522ed-8576-443a-aaf8-357463cb9ea4',
})


def _binding_in_receipt(receipt, binding):
    return any(r.get('path') == binding['path'] and r.get('sha256') == binding['sha256']
               for r in receipt if isinstance(r, dict))


@dataclass(frozen=True)
class ExactCellDecision:
    scope_id: str
    cell_key: str
    cluster_id: str
    prior_claim: str
    prior_annotation: str
    annotation: str
    reason: str
    evidence_sha256: str
    positive_linked_nucleus: bool


@dataclass(frozen=True)
class ConditionalCellSelection:
    base: semantic.SemanticSelection = field(repr=False)
    decisions: dict = field(repr=False)
    bindings: dict

    def summary(self):
        return {'status': 'CONDITIONAL_EXACT_CELL_CANDIDATE_VALIDATED_NOT_ACTIVATED',
                'reviewed_exact_cells': len(self.decisions), 'affected_scopes': sorted({k[0] for k in self.decisions}),
                'explicit_purpose': PURPOSE, 'policy': POLICY, 'bindings': self.bindings,
                'supported_scopes': sorted(REVIEWED_SCOPES),
                'all_other_scopes_require_unchanged_base_API': True,
                'complete_iterator_exhaustion_required_for_membership_validation': True,
                'native_and_base_semantic_columns_preserved': True,
                'base_cell_files_read_or_hashed': False, 'production_alias_or_latch_changed': False,
                'biological_identity_accuracy_validated': False}

    def scope_metadata(self, scope_id):
        require(scope_id in REVIEWED_SCOPES,
                'This bounded candidate supports only the two reviewed scopes; use the unchanged base API elsewhere')
        return self.base.scope_metadata(scope_id)

    def _resolve_batch(self, scope_id, table):
        """Append operative candidate fields; never replace any source/base value."""
        scope = self.scope_metadata(scope_id).fields
        require(not set(DERIVED).intersection(table.column_names), 'Refusing to overwrite existing candidate-view fields')
        protein = scope['lane'] in semantic.protein.PROTEIN_SCOPES
        key_column = scope['cell_key_column']
        required = {key_column, scope['claim_column'], scope['QC_column']}
        if protein:
            required.update({'annotation_name', 'eligible_for_broad_lineage', 'eligible_for_programme',
                             'source_use_status', 'cluster_id'})
        else:
            required.update({'candidate_identity_uncertainty', 'candidate_identity_withheld_after_exact_unit_review'})
        additional = scope.get('additional_admission_column')
        if additional:
            required.add(additional)
        require(required.issubset(table.column_names), 'Candidate batch lacks required preserved source fields')
        if scope['QC_kind'] == 'boolean':
            require(pa.types.is_boolean(table[scope['QC_column']].type), 'Native QC must remain boolean')
        else:
            require(scope['QC_kind'] == 'nonnull', 'Unsupported native QC contract')
        if additional:
            require(pa.types.is_boolean(table[additional].type), 'Additional admission must remain boolean')
        values = {name: [] for name in DERIVED}
        for row in table.to_pylist():
            key = row[key_column]
            d = self.decisions.get((scope_id, key))
            qc_value = row[scope['QC_column']]
            qc = qc_value if scope['QC_kind'] == 'boolean' else qc_value is not None
            held = scope_id == semantic.protein.HELD_SCOPE or row.get('source_use_status') == 'held_source_identity'
            admission = row[additional] if additional else True
            annotation = row['annotation_name'] if protein else row[scope['claim_column']]
            eligible = row['eligible_for_broad_lineage'] if protein else admission
            uncertainty = None if protein else row['candidate_identity_uncertainty']
            status = row.get('evidence_level', 'preserved_base_candidate_identity')
            if d:
                require(d.scope_id == scope_id and d.cell_key == key, 'Exact reviewed key mismatch')
                require(row[scope['claim_column']] == d.prior_claim and row['cluster_id'] == d.cluster_id,
                        'Reviewed exact cell no longer matches its pinned prior claim/cluster')
                require(annotation == d.prior_annotation, 'Reviewed exact cell base annotation changed')
                require(d.positive_linked_nucleus is True, 'Reviewed identity requires positive linked-native-nucleus evidence')
                require(qc is True, 'Reviewed exact-cell annotation cannot bypass native QC or missing QC')
                require(not held, 'Reviewed exact-cell annotation cannot bypass a physical-source hold')
                if protein:
                    require(admission is True, 'Reviewed protein annotation cannot bypass independent additional admission')
                    require(row['eligible_for_programme'] is True,
                            'These reviewed protein cases require preserved base scientific program admission')
                    require('mandatory_analysis_exclusion' not in row or row['mandatory_analysis_exclusion'] is False,
                            'Reviewed protein annotation cannot bypass a mandatory analysis exclusion')
                    require('analysis_eligible' not in row or row['analysis_eligible'] is True,
                            'Reviewed protein annotation cannot bypass independent analysis ineligibility')
                else:
                    # This exact column is a named-lineage use decision, not native QC.
                    # Preserve its old value and expose the separately reviewed candidate use.
                    require(additional == 'candidate_common_lineage_analysis_eligible',
                            'Cannot reinterpret an unknown RNA additional-admission policy')
                    require(row['candidate_identity_withheld_after_exact_unit_review'] is False,
                            'Reviewed RNA annotation cannot bypass an explicit prior exact-unit withholding')
                annotation, eligible, uncertainty, status = d.annotation, True, True, POLICY['identity_status']
            # The accepted base semantic view owns unreviewed eligibility,
            # including unavailable decisions on excluded rows. Preserve its
            # missingness exactly; apply the above hard checks only to reviews.
            values['effective_candidate_annotation_name'].append(annotation)
            values['effective_candidate_broad_lineage_eligible'].append(eligible)
            values['effective_candidate_identity_status'].append(status)
            values['effective_candidate_identity_uncertainty'].append(uncertainty)
            values['exact_cell_review_applied'].append(d is not None)
            values['exact_cell_review_reason'].append(d.reason if d else None)
            values['exact_cell_review_evidence_sha256'].append(d.evidence_sha256 if d else None)
            values['exact_cell_review_subtype'].append(None)
            values['exact_cell_review_function'].append(None)
            values['reviewed_annotation_truth_training_tuning_eligible'].append(False)
        for name, typ in DERIVED.items():
            table = table.append_column(name, pa.array(values[name], type=typ))
        return table

    def iter_annotations(self, scope_id, *, verify_file=True, batch_size=semantic.BATCH_SIZE):
        """Append candidate fields; exhaust the stream to verify complete membership."""
        require(verify_file is True, 'Conditional candidate iteration requires source-byte hash verification')
        seen = set()
        expected = {key for sid, key in self.decisions if sid == scope_id}
        key_column = self.scope_metadata(scope_id).fields['cell_key_column']
        for table in self.base.iter_annotations(scope_id, verify_file=verify_file, batch_size=batch_size):
            keys = table[key_column].to_pylist()
            hits = [key for key in keys if key in expected]
            require(len(hits) == len(set(hits)) and not seen.intersection(hits), 'Duplicate reviewed exact cell in source stream')
            seen.update(hits)
            yield self._resolve_batch(scope_id, table)
        require(seen == expected, 'A reviewed exact cell is absent from the complete selected scope stream')

    def iter_for_purpose(self, scope_id, *, purpose, verify_file=True, batch_size=semantic.BATCH_SIZE):
        """Select explicit candidate use; exhaust even when some batches yield no rows."""
        require(purpose == PURPOSE, 'Only the explicit conditional candidate broad-lineage purpose is supported')
        for table in self.iter_annotations(scope_id, verify_file=verify_file, batch_size=batch_size):
            selected = table.filter(pc.fill_null(table['effective_candidate_broad_lineage_eligible'], False))
            if len(selected):
                yield selected

    def counts(self, scope_id, *, purpose, verify_file=True, batch_size=semantic.BATCH_SIZE):
        require(purpose == PURPOSE, 'Only the explicit conditional candidate broad-lineage purpose is supported')
        status, annotations, total, reviewed = Counter(), Counter(), 0, 0
        for table in self.iter_annotations(scope_id, verify_file=verify_file, batch_size=batch_size):
            for row in table.select(['effective_candidate_annotation_name', 'effective_candidate_broad_lineage_eligible', 'exact_cell_review_applied']).to_pylist():
                total += 1
                reviewed += row['exact_cell_review_applied']
                status[row['effective_candidate_broad_lineage_eligible']] += 1
                if row['effective_candidate_broad_lineage_eligible'] is True:
                    annotations[row['effective_candidate_annotation_name']] += 1
        return {'scope_id': scope_id, 'purpose': purpose, 'source_rows': total, 'reviewed_exact_cells': reviewed,
                'eligible_rows': status[True], 'ineligible_rows': status[False], 'unavailable_rows': status[None],
                'eligible_annotation_counts': dict(annotations), 'counts_across_scopes_additive': False,
                'conditional_annotations_are_independent_truth': False}


def load_candidate_selection(path):
    path = Path(path).resolve()
    contract = json.loads(path.read_text())
    require(type(contract.get('schema_version')) is int and contract.get('schema_version') == 1 and
            contract.get('kind') == 'conditional_exact_cell_typing_candidate',
            'An explicit conditional exact-cell candidate contract is required')
    require(contract.get('policy') == POLICY, 'Conditional identity/use policy changed or incomplete')
    names = ['base_selection', 'candidate_api', 'root_scientific_acceptance', 'RNA_scientific_decisions', 'protein_scientific_decisions']
    bindings = {k: semantic.checked(contract.get(k), path.parent) for k in names}
    semantic.require_executed_api(bindings['candidate_api'], __file__, 'conditional exact-cell candidate API')
    base = semantic.load_semantic_selection(bindings['base_selection']['path'])
    root = json.loads(Path(bindings['root_scientific_acceptance']['path']).read_text())
    rna = json.loads(Path(bindings['RNA_scientific_decisions']['path']).read_text())
    protein = json.loads(Path(bindings['protein_scientific_decisions']['path']).read_text())
    require(root.get('typing_accuracy_established') is False and root.get('goal_complete') is False,
            'Root acceptance must not promote conditional cases to accuracy or completed goal')
    require(_binding_in_receipt(root.get('completed_evidence', []), bindings['base_selection']) and
            _binding_in_receipt(root.get('completed_evidence', []), bindings['RNA_scientific_decisions']) and
            _binding_in_receipt(root.get('completed_evidence', []), bindings['protein_scientific_decisions']),
            'Root scientific record does not bind the exact base and both independent decisions')
    membership_inputs = [r for r in rna.get('source_inputs', []) if isinstance(r, dict) and
                         r.get('path', '').endswith('/independent_design_and_membership_review_direct_validation_r1.json')]
    require(len(membership_inputs) == 1, 'RNA case review requires its exact saved membership-review binding')
    bindings['RNA_membership_review'] = semantic.checked(membership_inputs[0], path.parent)
    membership = json.loads(Path(bindings['RNA_membership_review']['path']).read_text())
    checked_inputs = membership.get('checked_saved_inputs', {})
    require(checked_inputs.get('selection') == bindings['base_selection'] and
            _binding_in_receipt(protein.get('reviewed_source_inputs', []), bindings['base_selection']),
            'Separate agent case reviews do not bind this exact base selection')
    rna_scope = base._scopes['pLN_RNA|0037502|6380']
    require(checked_inputs.get('current_sidecar') ==
            {'path': str(Path(rna_scope['path']).resolve()), 'sha256': rna_scope['sha256']},
            'RNA membership review does not bind the exact selected source sidecar')
    require(rna.get('independent_accuracy_established') is False and protein.get('independent_accuracy_established') is False,
            'Independent case reviews must retain no cohort-accuracy claim')
    decisions = {}
    accepted_rna = [r for r in rna['cases'] if r.get('accepted_conditional_broad_annotation')]
    root_rna = set(root['new_RNA_evidence']['source_qualified_accepted_keys'])
    require({r['serial_qualified_cell_key'] for r in accepted_rna} == root_rna,
            'Root RNA accepted exact keys differ from independent case decisions')
    for r in accepted_rna:
        key = r['serial_qualified_cell_key']
        matches = [s for s in base._scopes.values() if s['lane'] == 'pLN_RNA' and
                   key == 'pLN|' + s['section_id'] + '|' + s['source_serial'] + '|' + r['cell_id']]
        require(len(matches) == 1, 'RNA reviewed key must identify exactly one native section and serial')
        s = matches[0]
        require(type(r['native_nucleus_count']) is int and r['native_nucleus_count'] == 1 and
                r['accepted_subtype'] is None and r['accepted_function'] is None,
                'RNA conditional case must remain single-nucleus broad-only')
        require(r['review_decision'] in {'accept', 'accept_qualified_sparse'} and r['accepted_conditional_broad_annotation'] == 'T-lineage cell',
                'Unsupported RNA exact-cell scientific decision')
        d = ExactCellDecision(s['scope_id'], key, '6380_c02', 'T_NK_unresolved', 'T_NK_unresolved',
                             r['accepted_conditional_broad_annotation'], r['literal_reason'], bindings['RNA_scientific_decisions']['sha256'], True)
        require((d.scope_id, key) not in decisions, 'Duplicate reviewed exact-cell key')
        decisions[(d.scope_id, key)] = d
    accepted_protein = protein['cases']
    root_protein = {(r['scope_id'], r['object_id'], r['annotation'], r['decision']) for r in root['new_protein_evidence']['accepted_conditional_case_assignments']}
    require({(r['scope_id'], r['object_id'], r['accepted_annotation'], r['review_decision']) for r in accepted_protein} == root_protein,
            'Root protein accepted exact keys differ from independent case decisions')
    for r in accepted_protein:
        sid, key = r['scope_id'], r['object_id']
        require(sid in base._scopes and base._scopes[sid]['lane'] in semantic.protein.PROTEIN_SCOPES and sid != semantic.protein.HELD_SCOPE,
                'Reviewed protein source unknown or held')
        require(r['accepted_annotation'] in {'B-lineage cell', 'T-lineage cell'} and r['accepted_subtype'] is None and r['accepted_function'] is None,
                'Protein conditional case must remain broad-only')
        require(r['numeric_checks']['nucleus_pixels'] > 0 and r['numeric_checks']['independent_pixel_center_raster_matches_saved_masks'] is True,
                'Protein conditional case lacks linked-native-nucleus proof')
        d = ExactCellDecision(sid, key, r['cluster_id'], r['exact_current_prior_claim'], r['current_annotation_name'],
                             r['accepted_annotation'], r['literal_visual_evidence'], bindings['protein_scientific_decisions']['sha256'], True)
        require((sid, key) not in decisions, 'Duplicate reviewed exact-cell key')
        decisions[(sid, key)] = d
    require(len(decisions) == contract.get('expected_exact_reviewed_cells') == 14 and len(accepted_rna) == 12 and len(accepted_protein) == 2,
            'This bounded candidate requires exactly the accepted 12 RNA and two protein cells')
    require({sid for sid, key in decisions} == REVIEWED_SCOPES,
            'Accepted case reviews must remain bounded to the two explicitly supported scopes')
    rows = contract.get('reviewed_cells')
    require(isinstance(rows, list) and len(rows) == 14, 'The candidate contract must contain the exact14 reviewed cells')
    seen = set()
    for row in rows:
        key = (row.get('scope_id'), row.get('cell_key'))
        require(key in decisions and key not in seen, 'Duplicate or foreign exact-cell contract key')
        seen.add(key)
        d = decisions[key]
        scope = base._scopes[key[0]]
        require(row.get('source_scope') == {'path': str(Path(scope['path']).resolve()), 'sha256': scope['sha256']},
                'Exact-cell source path/hash differs from the selected scope')
        require(row.get('cluster_id') == d.cluster_id and row.get('prior_claim') == d.prior_claim and
                row.get('annotation') == d.annotation and row.get('evidence_sha256') == d.evidence_sha256,
                'Exact-cell contract differs from its independent scientific decision')
        require(row.get('subtype') is None and row.get('function') is None,
                'These14 conditional cases require null subtype/function')
        proof = semantic.checked(row.get('nucleus_proof'), path.parent)
        evidence_inputs = rna['source_inputs'] if scope['lane'] == 'pLN_RNA' else protein['reviewed_source_inputs']
        require(_binding_in_receipt(evidence_inputs, proof), 'Native nucleus proof is not bound by the independent review')
        if scope['lane'] == 'pLN_RNA':
            require(type(row.get('linked_native_nucleus_count')) is int and row.get('linked_native_nucleus_count') == 1 and
                    row.get('positive_linked_native_nucleus_geometry') is None,
                    'RNA reviewed cases require a source-counted single native nucleus')
            require(row.get('nucleus_proof_kind') == 'native_count_csv', 'Unsupported RNA nucleus proof')
            with Path(proof['path']).open(newline='') as handle:
                matching = [r for r in csv.DictReader(handle) if r['serial_qualified_cell_key'] == key[1]]
            require(len(matching) == 1 and matching[0]['native_nucleus_count'] == '1' and matching[0]['source_native_QC'] == 'True',
                    'RNA native proof does not establish this exact QC-positive one-nucleus segment')
        else:
            require(row.get('linked_native_nucleus_count') is None and row.get('positive_linked_native_nucleus_geometry') is True,
                    'Protein proof must state positive linked geometry without inventing a nucleus count')
            require(row.get('nucleus_proof_kind') == 'linked_nucleusGeometry_geojson', 'Unsupported protein nucleus proof')
            geometry = json.loads(Path(proof['path']).read_text())
            matching = [r for r in geometry['features'] if r.get('id') == key[1]]
            require(len(matching) == 1 and matching[0].get('nucleusGeometry', {}).get('type') == 'Polygon' and
                    bool(matching[0]['nucleusGeometry'].get('coordinates')),
                    'Protein native proof requires exactly one explicit linked nucleus geometry')
        bindings['nucleus_proof:' + key[1]] = proof
    require(seen == set(decisions), 'The candidate contract omits a reviewed exact cell')
    bindings['candidate_contract'] = {'path': str(path), 'sha256': semantic.base.sha256(path)}
    return ConditionalCellSelection(base, MappingProxyType(decisions), bindings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('validate')
    counts = sub.add_parser('counts')
    counts.add_argument('--scope-id', required=True)
    counts.add_argument('--purpose', choices=[PURPOSE], required=True)
    args = parser.parse_args()
    selection = load_candidate_selection(args.candidate)
    result = selection.summary() if args.command == 'validate' else selection.counts(args.scope_id, purpose=args.purpose)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
