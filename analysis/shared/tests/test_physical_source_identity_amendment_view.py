"""Meaningful source-only eligibility, missingness, and exact-binding regressions.

Synthetic software fixtures are not cell-identity truth or accuracy evaluations.
"""
from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pyarrow as pa
import pytest

from analysis.shared import physical_source_identity_amendment_view as view

OTHER_PROTEIN = 'dedicated_pLN_protein|6523|unrelated-capture'
OTHER_RNA = 'pancreas_RNA|6623|unrelated-capture'


class FakeBase:
    def __init__(self, scope, table, routes):
        self.scope, self.table = scope, table
        self._scopes = {scope['scope_id']: scope}
        self._overlay = SimpleNamespace(pan_routes=routes)
        self.calls = []

    def _scope(self, sid):
        assert sid == self.scope['scope_id']
        return self.scope

    def scope_metadata(self, sid):
        self.calls.append(('metadata', sid))
        return view.semantic.ScopeMetadata(sid, self.scope['lane'], view.semantic.freeze(self.scope))

    def counts(self, sid, **kwargs):
        self.calls.append(('counts', sid, kwargs))
        return {'same_base_result': True, 'scope_id': sid}

    def iter_annotations(self, sid, **kwargs):
        self.calls.append(('annotations', sid, kwargs))
        for offset in range(0, len(self.table), kwargs['batch_size']):
            yield self.table.slice(offset, kwargs['batch_size'])

    def iter_for_purpose(self, sid, **kwargs):
        self.calls.append(('purpose', sid, kwargs))
        yield self.table


def synthetic(qc=True, admission=True, explicit=True, missing_decision=False):
    s = {'scope_id': view.SCOPE_ID, 'lane': 'pancreas_protein', 'donor_id': '6523',
         'QC_column': 'qc_analysis_eligible', 'QC_kind': 'boolean',
         'additional_admission_column': 'analysis_eligible',
         'cell_key_column': 'object_id', 'claim_column': 'candidate_cell_class',
         'expected_source_rows': 1, 'expected_QC_rows': int(qc is True),
         'expected_additional_admission_rows': int(admission is True),
         'source_identity_status': 'TIFF_SlideID_discrepancy_under_review'}
    d = None if missing_decision else view.protein.Decision('Same conditional annotation',
        'qualified_broad_lineage_programme', 'Same root scientific reason', (explicit,) * 3)
    table = pa.table({'object_id': ['same-native-id'], 'candidate_cell_class': ['same-route'],
        'native_opaque_missing': pa.array([None], type=pa.string()),
        'qc_analysis_eligible': pa.array([qc], type=pa.bool_()),
        'analysis_eligible': pa.array([admission], type=pa.bool_())})
    field = table.schema.field('object_id').with_metadata({b'native': b'original ID'})
    table = table.set_column(0, field, table['object_id']).replace_schema_metadata({b'source': b'unchanged'})
    overlay = SimpleNamespace(scopes={view.SCOPE_ID: s}, pan_routes={'same-route': d},
                              source_hold_reason='Same exact physical source hold')
    table = view.protein.resolve_batch(overlay, view.SCOPE_ID, table)
    base = FakeBase(s, table, {'same-route': d})
    return base, table


@pytest.mark.parametrize('qc,admission,explicit,missing,expected', [
    (True, True, True, False, True), (True, True, False, False, False),
    (False, False, True, False, False), (True, False, True, False, False),
    (None, None, True, False, None), (True, None, True, False, None),
    (False, None, True, False, False), (False, False, None, False, None),
    (True, True, None, False, None), (False, False, True, True, None),
    (True, True, True, True, None),
])
def test_existing_nullable_decision_QC_and_admission_semantics(qc, admission, explicit, missing, expected):
    base, original = synthetic(qc, admission, explicit, missing)
    got = view.amend_batch(base, view.SCOPE_ID, original)
    for p in view.protein.SCIENTIFIC_PURPOSES:
        assert got['eligible_for_' + p].to_pylist() == [expected]
        assert got['eligibility_reason_for_' + p].to_pylist() == [
            view.protein.effective_eligibility(base._overlay.pan_routes['same-route'], p, qc, admission, False)[1]]
    changed = {'source_use_status', 'source_use_reason',
        *(f'eligible_for_{p}' for p in view.protein.SCIENTIFIC_PURPOSES),
        *(f'eligibility_reason_for_{p}' for p in view.protein.SCIENTIFIC_PURPOSES)}
    unchanged = [c for c in original.column_names if c not in changed]
    assert got.select(unchanged).equals(original.select(unchanged), check_metadata=True)
    assert got.schema.equals(original.schema, check_metadata=True)
    assert got['source_use_status'].to_pylist() == [view.STATUS]
    assert got['eligible_for_native_provenance'].to_pylist() == [True]
    assert len(got.column_names) - 5 == 14


@pytest.mark.parametrize('sid', [OTHER_PROTEIN, OTHER_RNA])
def test_unrelated_scopes_delegate_values_missingness_metadata_and_purpose(sid):
    base, original = synthetic(None, None, None)
    amendment = view.PhysicalSourceIdentityAmendment(base, {})
    assert list(amendment.iter_annotations(sid, verify_file=False, batch_size=1))[0].equals(original, check_metadata=True)
    assert list(amendment.iter_for_purpose(sid, purpose='native_provenance', verify_file=False, batch_size=1))[0] is original
    assert amendment.counts(sid, purpose=None, verify_file=False) == {'same_base_result': True, 'scope_id': sid}
    assert amendment.scope_metadata(sid).fields == base.scope_metadata(sid).fields
    assert all(call[1] == sid for call in base.calls)


def test_direct_amendment_rejects_foreign_scope_and_previously_amended_alias():
    base, original = synthetic()
    with pytest.raises(ValueError, match='exact ordinary'):
        view.amend_batch(base, OTHER_PROTEIN, original)
    amended = view.amend_batch(base, view.SCOPE_ID, original)
    with pytest.raises(ValueError, match='held base'):
        view.amend_batch(base, view.SCOPE_ID, amended)


def test_full_stream_denominator_failure_is_not_hidden_by_empty_purpose_batches():
    base, _ = synthetic(False, False, False)
    base.scope['expected_source_rows'] = 2
    amendment = view.PhysicalSourceIdentityAmendment(base, {})
    with pytest.raises(ValueError, match='denominator'):
        list(amendment.iter_for_purpose(view.SCOPE_ID, purpose='broad_lineage'))


def contract():
    base = view.semantic.load_semantic_selection(view.ENTRY)
    s = base._scope(view.SCOPE_ID)
    return {'schema_version': 1, 'kind': 'physical_source_identity_amendment',
        'scope_id': view.SCOPE_ID, 'capture': deepcopy(view.CAPTURE), 'policy': deepcopy(view.POLICY),
        'base_selection': view.pin(view.ENTRY), 'base_source': {'path': s['path'], 'sha256': s['sha256']},
        'base_selection_api': view.pin(view.semantic.__file__), 'base_semantic_api': view.pin(view.protein.__file__),
        'base_source_api': view.pin(view.semantic.base.__file__), 'amendment_api': view.pin(view.__file__),
        'builder': view.pin(view.ROOT / 'analysis/shared/build_physical_source_identity_amendment.py'),
        'human_attestation': view.pin(view.ATTESTATION),
        'current_raw_header_metadata_readback': view.pin(view.HEADER_READBACK)}


def write_contract(tmp_path, value):
    path = tmp_path / 'amendment.json'
    path.write_text(json.dumps(value))
    return path


def test_real_compact_binding_load_does_not_read_or_hash_base_cells(tmp_path, monkeypatch):
    c = contract()
    real_hash = view.semantic.base.sha256
    def guarded_hash(path):
        assert not str(path).endswith('.parquet')
        return real_hash(path)
    monkeypatch.setattr(view.semantic.base, 'sha256', guarded_hash)
    loaded = view.load_source_identity_amendment(write_contract(tmp_path, c))
    assert loaded.summary()['n_unchanged_delegated_scopes'] == 107
    for sid in loaded._base._scopes:
        if sid != view.SCOPE_ID:
            assert loaded.scope_metadata(sid) == loaded._base.scope_metadata(sid)
    metadata = loaded.scope_metadata(view.SCOPE_ID).fields
    assert metadata['source_identity_status'] == 'TIFF_SlideID_discrepancy_under_review'
    assert metadata['physical_source_identity_amendment']['attested_corrected_donor_id'] == '6523'
    assert metadata['physical_source_identity_amendment']['native_SlideID_preserved'] == '6623panc'


@pytest.mark.parametrize('mutation', [
    'scope', 'donor', 'UUID', 'timestamp', 'header', 'barcode_missing', 'policy_missing',
    'accuracy', 'numeric_bool', 'string_bool', 'nullable_bool', 'alias', 'source_hash',
    'source_path', 'selection_hash', 'foreign_repo', 'attestation_hash', 'api_alias', 'float_shape',
])
def test_real_contract_rejects_foreign_or_incomplete_authority(tmp_path, mutation):
    c = contract()
    if mutation == 'scope': c['scope_id'] = OTHER_PROTEIN
    elif mutation == 'donor': c['capture']['donor_id'] = '6623'
    elif mutation == 'UUID': c['capture']['native_DAPI_Identifier'] = 'foreign-native-id'
    elif mutation == 'timestamp': c['capture']['native_TIFF_DateTime'] = '2026:06:16 16:56:37'
    elif mutation == 'header': c['capture']['native_SlideID'] = '6523panc'
    elif mutation == 'barcode_missing': del c['capture']['native_Barcode']
    elif mutation == 'policy_missing': del c['policy']['physical_source_identity_only']
    elif mutation == 'accuracy': c['policy']['biological_identity_accuracy_validated'] = True
    elif mutation == 'numeric_bool': c['policy']['production_activation'] = 0
    elif mutation == 'string_bool': c['policy']['production_activation'] = 'false'
    elif mutation == 'nullable_bool': c['policy']['production_activation'] = None
    elif mutation == 'alias': c['production_alias'] = 'current'
    elif mutation == 'source_hash': c['base_source']['sha256'] = '0' * 64
    elif mutation == 'source_path': c['base_source']['path'] = str(tmp_path / 'foreign.parquet')
    elif mutation == 'selection_hash': c['base_selection']['sha256'] = '0' * 64
    elif mutation == 'foreign_repo':
        copy = tmp_path / 'foreign_selection.json'
        copy.write_bytes(view.ENTRY.read_bytes())
        c['base_selection'] = view.pin(copy)
    elif mutation == 'attestation_hash': c['human_attestation']['sha256'] = '0' * 64
    elif mutation == 'api_alias':
        copy = tmp_path / 'alias.py'
        copy.write_bytes(Path(view.__file__).read_bytes())
        c['amendment_api'] = view.pin(copy)
    elif mutation == 'float_shape': c['capture']['raw_full_resolution_shape_yx'][0] = 30240.0
    with pytest.raises(ValueError):
        view.load_source_identity_amendment(write_contract(tmp_path, c))
