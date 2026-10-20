#!/usr/bin/env python3
"""Opt into one human-attested capture-name correction without changing typing.

All other scopes delegate to the immutable schema-2 semantic selection. Compact
validation reads no source cells or raw image. Exhaust annotation/purpose streams
to complete their denominator checks; an early stop is not a full readback.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Mapping

import pyarrow as pa
import pyarrow.compute as pc

try:
    from . import typing_semantic_selection_view as semantic
except ImportError:
    import typing_semantic_selection_view as semantic

protein = semantic.protein
require = protein.require
ROOT = Path(__file__).resolve().parents[2]
EXPECTED_ROOT = Path('/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior-v3')
PEER = ROOT / 'results/v3/islet_immune_reassessment_20261006/cell_typing_validation/peer_review_typing_20261008'
SCOPE_ID = 'pancreas_protein|6523|ordinary_pancreas'
STATUS = 'source_bound_user_attested_name_correction'
SOURCE_REASON = ('Direct user attestation corrects only the scan setup name for donor 6523; '
                 'native SlideID 6623panc and capture provenance are preserved. '
                 'Biological identity accuracy is not established.')
ENTRY = PEER / 'selection_inputs/pancreas_r13_pLN_r10_protein_semantic_r1/typing_semantic_selection.json'
ENTRY_SHA = '103bb8ffd909f2640c86fb5cf221333c393559737979a33b8029c067c004b333'
ATTESTATION_DIR = PEER / 'accuracy_and_all_outcomes_20261010_r1/PAN6523_independent_primary_attestation_and_amendment_review_B_r1'
ATTESTATION = ATTESTATION_DIR / 'primary_human_capture_identity_attestation_r2.json'
ATTESTATION_SHA = '7042d883f819e51ee0f56a4766a205dc83b30624b66d35340a5f95b97116ddbe'
HEADER_READBACK = ATTESTATION_DIR / 'root_current_raw_header_metadata_readback.json'
HEADER_READBACK_SHA = 'fba0b8d6e806bca25a0f821d6aa2bd98efb697ec941d0b445ac9130a3fbc3b6a'
SOURCE_SHA = 'cd7920aa6d54015a75fd143ad993a222930129f7d92938933a0095aba056086c'
CAPTURE = {
    'scope_id': SCOPE_ID, 'donor_id': '6523',
    'acquisition_identifier': '6523panc_Scan1.er.qptiff - resolution #1',
    'raw_image_path': '/home/smith6jt/IO60panc2nd/Images/6523panc_Scan1.er.qptiff',
    'raw_full_TIFF_SHA256': 'dc42efc2624637886530aa4f95c69459ff4ca3001b8e1f91d383843c7f2c273b',
    'native_DAPI_Identifier': '3f5a1d0b-38e5-4108-9aca-a6e649bd2c60',
    'native_SlideID': '6623panc', 'native_TIFF_DateTime': '2026:06:21 14:17:16',
    'raw_full_resolution_shape_yx': [30240, 24000], 'n_acquired_channels': 58,
    'pixel_type': 'UINT16', 'native_Barcode': None,
}
POLICY = {
    'physical_source_identity_only': True,
    'preserve_native_header_and_all_typing_annotations': True,
    'preserve_QC_exclusions_and_nullable_decisions': True,
    'across_scope_counts_are_additive': False,
    'production_activation': False,
    'production_alias_changes': False,
    'biological_identity_accuracy_validated': False,
    'independent_cell_identity_truth_established': False,
    'training_or_tuning_eligibility_created': False,
    'subtype_or_function_upgrade': False,
    'conditional14_candidate_changes': False,
}
CONTRACT_KEYS = {'schema_version', 'kind', 'scope_id', 'capture', 'policy',
                 'base_selection', 'base_source', 'base_selection_api',
                 'base_semantic_api', 'base_source_api', 'amendment_api',
                 'human_attestation', 'current_raw_header_metadata_readback', 'builder'}


def pin(path):
    path = Path(path).resolve()
    return {'path': str(path), 'sha256': semantic.base.sha256(path)}


def exact_pin(item, path, sha, directory):
    checked = protein.checked_binding(item, directory)
    require(checked == {'path': str(Path(path).resolve()), 'sha256': sha},
            'Binding is not the exact authorized source: ' + str(path))
    return checked


def amend_batch(base, scope_id, table):
    """Replace only source status/reason and existing scientific eligibility pairs."""
    require(scope_id == SCOPE_ID, 'Only the exact ordinary-PAN6523 scope may be amended')
    s = base._scope(scope_id)
    require(s['lane'] == 'pancreas_protein' and s['donor_id'] == '6523', 'Foreign source scope')
    require(set(protein.DERIVED_COLUMNS).issubset(table.column_names), 'Missing base semantic fields')
    require(pa.types.is_boolean(table[s['QC_column']].type), 'Native QC must remain boolean')
    additional = s.get('additional_admission_column')
    require(additional and pa.types.is_boolean(table[additional].type),
            'Exact scope requires preserved boolean additional admission')
    qc, admission = table[s['QC_column']].to_pylist(), table[additional].to_pylist()
    require(not any(a is True and q is not True for q, a in zip(qc, admission, strict=True)),
            'Native additional admission cannot bypass native QC')
    require(all(v == 'held_source_identity' for v in table['source_use_status'].to_pylist()),
            'Amendment requires the unchanged held base, never a previously amended alias')
    values = {'source_use_status': [STATUS] * len(table),
              'source_use_reason': [SOURCE_REASON] * len(table)}
    labels = table[s['claim_column']].to_pylist()
    for purpose in protein.SCIENTIFIC_PURPOSES:
        pairs = [protein.effective_eligibility(base._overlay.pan_routes.get(label), purpose, q, a, False)
                 for label, q, a in zip(labels, qc, admission, strict=True)]
        values['eligible_for_' + purpose] = [p[0] for p in pairs]
        values['eligibility_reason_for_' + purpose] = [p[1] for p in pairs]
    for name, value in values.items():
        index = table.schema.get_field_index(name)
        field = table.schema.field(index)
        table = table.set_column(index, field, pa.array(value, type=field.type))
    return table


@dataclass(frozen=True)
class PhysicalSourceIdentityAmendment:
    _base: semantic.SemanticSelection
    bindings: Mapping

    def scope_metadata(self, scope_id):
        old = self._base.scope_metadata(scope_id)
        if scope_id != SCOPE_ID:
            return old
        fields = dict(old.fields)
        fields['physical_source_identity_amendment'] = {
            'effective_source_use_status': STATUS,
            'attested_corrected_donor_id': '6523',
            'native_SlideID_preserved': CAPTURE['native_SlideID'],
            'native_capture_UUID_preserved': CAPTURE['native_DAPI_Identifier'],
            'native_TIFF_DateTime_preserved': CAPTURE['native_TIFF_DateTime'],
            'human_attestation': dict(self.bindings['human_attestation']),
            'base_source_identity_status_preserved': fields.get('source_identity_status'),
            'biological_identity_accuracy_validated': False,
        }
        return semantic.ScopeMetadata(old.scope_id, old.lane, semantic.freeze(fields))

    def summary(self):
        return {'status': 'source_identity_amendment_validated_not_activated',
                'n_scopes': len(self._base._scopes), 'n_amended_scopes': 1,
                'amended_scope_id': SCOPE_ID, 'n_unchanged_delegated_scopes': 107,
                'source_use_status': STATUS, 'capture': CAPTURE,
                'bindings': self.bindings, 'policy': POLICY,
                'base_cell_files_read_or_hashed': False,
                'raw_images_read_or_hashed': False,
                'stream_verification_requires_full_exhaustion': True}

    def iter_annotations(self, scope_id, *, columns=None, verify_file=True, batch_size=8192):
        if scope_id != SCOPE_ID:
            yield from self._base.iter_annotations(scope_id, columns=columns,
                verify_file=verify_file, batch_size=batch_size)
            return
        s = self._base._scope(scope_id)
        n = q = a = 0
        for table in self._base.iter_annotations(scope_id, columns=columns,
                verify_file=verify_file, batch_size=batch_size):
            n += len(table)
            q += sum(v is True for v in table[s['QC_column']].to_pylist())
            a += sum(v is True for v in table[s['additional_admission_column']].to_pylist())
            yield amend_batch(self._base, scope_id, table)
        require((n, q, a) == (s['expected_source_rows'], s['expected_QC_rows'],
                s['expected_additional_admission_rows']), 'Amended native denominator mismatch')

    def iter_for_purpose(self, scope_id, *, purpose, columns=None, verify_file=True, batch_size=8192):
        if scope_id != SCOPE_ID:
            yield from self._base.iter_for_purpose(scope_id, purpose=purpose, columns=columns,
                verify_file=verify_file, batch_size=batch_size)
            return
        require(purpose in protein.PURPOSES, 'Unknown analytical purpose')
        for table in self.iter_annotations(scope_id, columns=columns,
                verify_file=verify_file, batch_size=batch_size):
            selected = table.filter(pc.fill_null(table['eligible_for_' + purpose], False))
            if len(selected):
                yield selected

    def counts(self, scope_id, *, purpose=None, verify_file=True, batch_size=8192):
        if scope_id != SCOPE_ID:
            return self._base.counts(scope_id, purpose=purpose, verify_file=verify_file,
                                     batch_size=batch_size)
        require(purpose in protein.PURPOSES, 'Amended protein counts require an explicit purpose')
        s = self._base._scope(scope_id)
        counts, annotations, joint = Counter(), Counter(), Counter()
        for table in self.iter_annotations(scope_id, columns=[s['claim_column']],
                verify_file=verify_file, batch_size=batch_size):
            for r in table.to_pylist():
                q, a, e = r[s['QC_column']], r[s['additional_admission_column']], r['eligible_for_' + purpose]
                counts['source_rows'] += 1
                counts['native_QC_rows'] += q is True
                counts['native_QC_unavailable_rows'] += q is None
                counts['native_additional_admission_rows'] += a is True
                counts['native_additional_admission_unavailable_rows'] += a is None
                counts['root_decision_unavailable_rows'] += r['semantic_decision_status'] == 'root_decision_unavailable'
                counts[{True: 'purpose_eligible_rows', False: 'purpose_ineligible_rows', None: 'purpose_unavailable_rows'}[e]] += 1
                joint[(q, a, r['source_use_status'], r['semantic_decision_status'], e)] += 1
                if e is True:
                    annotations[(r['annotation_name'], r['evidence_level'])] += 1
        result = {k: counts[k] for k in ('source_rows', 'native_QC_rows', 'native_QC_unavailable_rows',
            'native_additional_admission_rows', 'native_additional_admission_unavailable_rows',
            'held_source_rows', 'root_decision_unavailable_rows', 'purpose_eligible_rows',
            'purpose_ineligible_rows', 'purpose_unavailable_rows')}
        result.update(scope_id=scope_id, purpose=purpose,
            source_identity_amendment_applied=True,
            purpose_annotation_counts=[{'annotation_name': name, 'evidence_level': level, 'n': n}
                for (name, level), n in sorted(annotations.items(), key=lambda r: str(r[0]))],
            joint_status_counts=[{'native_QC': q, 'native_additional_admission': a,
                'source_use_status': source, 'semantic_decision_status': state, 'purpose_eligible': e, 'n': n}
                for (q, a, source, state, e), n in sorted(joint.items(), key=lambda r: str(r[0]))],
            separate_native_additional_admission_declared=True,
            across_scope_counts_are_additive=False, biological_identity_accuracy_validated=False,
            source_or_selection_mutations=False, production_activation=False)
        return result


def load_source_identity_amendment(path):
    """Load only the authorized additive contract and compact proof; no cell reads."""
    path = Path(path).resolve()
    c = json.loads(path.read_text())
    require(ROOT == EXPECTED_ROOT, 'Amendment API must execute in the exact authorized V3 repository')
    require(type(c) is dict and set(c) == CONTRACT_KEYS, 'Unknown/missing amendment contract fields or alias')
    require(type(c.get('schema_version')) is int and c['schema_version'] == 1,
            'Unsupported source amendment schema')
    require(c.get('kind') == 'physical_source_identity_amendment', 'Wrong amendment kind')
    require(c.get('scope_id') == SCOPE_ID and c.get('capture') == CAPTURE,
            'Foreign or incomplete capture/scope tuple')
    require(type(c.get('capture', {}).get('n_acquired_channels')) is int,
            'Capture channel count must be an integer')
    shape = c['capture']['raw_full_resolution_shape_yx']
    require(type(shape) is list and len(shape) == 2 and all(type(v) is int for v in shape),
            'Capture shape must preserve two actual integer dimensions')
    require(type(c.get('policy')) is dict and set(c['policy']) == set(POLICY) and
            all(type(c['policy'][k]) is bool and c['policy'][k] is v for k, v in POLICY.items()),
            'Source-only policy requires exact JSON booleans; truth, tuning and activation are prohibited')
    directory = path.parent
    for name in CONTRACT_KEYS - {'schema_version', 'kind', 'scope_id', 'capture', 'policy'}:
        require(type(c[name]) is dict and set(c[name]) == {'path', 'sha256'},
                'File binding must contain only an exact path and SHA256: ' + name)
    bindings = {
        'base_selection': exact_pin(c['base_selection'], ENTRY, ENTRY_SHA, directory),
        'human_attestation': exact_pin(c['human_attestation'], ATTESTATION, ATTESTATION_SHA, directory),
        'current_raw_header_metadata_readback': exact_pin(c['current_raw_header_metadata_readback'],
            HEADER_READBACK, HEADER_READBACK_SHA, directory),
    }
    for name, executing in [('base_selection_api', semantic.__file__),
                             ('base_semantic_api', protein.__file__),
                             ('base_source_api', semantic.base.__file__),
                             ('amendment_api', __file__)]:
        bindings[name] = protein.checked_binding(c[name], directory)
        semantic.require_executed_api(bindings[name], executing, name)
    bindings['builder'] = protein.checked_binding(c['builder'], directory)
    require(Path(bindings['builder']['path']) == ROOT / 'analysis/shared/build_physical_source_identity_amendment.py',
            'Foreign amendment builder')
    base = semantic.load_semantic_selection(bindings['base_selection']['path'])
    s = base._scope(SCOPE_ID)
    require(s['lane'] == 'pancreas_protein' and s['donor_id'] == '6523' and
            s['acquisition_identifier'] == CAPTURE['acquisition_identifier'] and
            s['QC_kind'] == 'boolean' and s.get('additional_admission_column') == 'analysis_eligible',
            'Exact base source semantics differ')
    require(c['base_source'] == {'path': s['path'], 'sha256': SOURCE_SHA},
            'Base source tuple differs from exact unchanged selected source')
    bindings['base_source'] = dict(c['base_source'])
    require(s['sha256'] == SOURCE_SHA, 'Base source hash differs from authorized selected source')
    attestation = json.loads(Path(bindings['human_attestation']['path']).read_text())
    require(attestation.get('schema_version') == 1 and
            attestation.get('kind') == 'human_attested_capture_source_identity' and
            attestation.get('status') == 'ACCEPTED_DIRECT_USER_SLIDE_NAME_CORRECTION' and
            all(attestation.get(k) == v for k, v in CAPTURE.items()), 'Invalid direct human capture attestation')
    message = attestation.get('root_verified_active_user_message', {})
    require(message.get('thread_id') == '01a11c26-4123-72a3-a371-4696e65ca8f5' and
            message.get('turn_id') == '01a1268f-8647-7152-b36f-fd27b776656a' and
            message.get('user_message_id') == '01a12695-5f6c-7151-924d-08b9c49c79a7' and
            message.get('exact_content') == 'I confirm the slide-name error\n',
            'Missing exact direct human confirmation')
    require({k: attestation['unchanged_selected_semantic_registry'][k] for k in ('path', 'sha256')} ==
            bindings['base_selection'], 'Attestation does not bind unchanged base selection')
    for k in ('native_header_XML', 'existing_raw_source_identity_audit', 'exact_capture_binding'):
        bindings[k] = protein.checked_binding(attestation[k], directory)
    bindings['preserved_original_human_record'] = protein.checked_binding(
        attestation['original_user_attestation']['exact_bounded_record'], directory)
    current = json.loads(Path(bindings['current_raw_header_metadata_readback']['path']).read_text())
    require(current.get('status') == 'PASS_CURRENT_RAW_DAPI_HEADER_BYTES_EQUAL_PINNED_CAPTURE_METADATA' and
            current.get('scope_id') == SCOPE_ID and current.get('current_raw_header_bytes_equal_pinned_XML_bytes') is True and
            current.get('native_XML_SHA256') == bindings['native_header_XML']['sha256'] and
            current.get('attestation_sha256') == ATTESTATION_SHA and
            current.get('biological_accuracy_established') is False,
            'Current capture metadata does not bind the accepted attestation')
    bindings['amendment_contract'] = pin(path)
    return PhysicalSourceIdentityAmendment(base, semantic.freeze(bindings))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--amendment', type=Path, required=True)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('validate', help='Compact contract and saved proof only')
    count = sub.add_parser('counts', help='Stream one exact scope with native denominator checks')
    count.add_argument('--scope-id', required=True)
    count.add_argument('--purpose', choices=protein.PURPOSES)
    count.add_argument('--output', type=Path)
    args = parser.parse_args()
    view = load_source_identity_amendment(args.amendment)
    result = view.summary() if args.command == 'validate' else view.counts(args.scope_id, purpose=args.purpose)
    text = json.dumps(result, indent=2, allow_nan=False,
                      default=lambda x: dict(x) if isinstance(x, Mapping) else x) + '\n'
    if args.command == 'counts' and args.output:
        with args.output.open('x') as stream:
            stream.write(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
