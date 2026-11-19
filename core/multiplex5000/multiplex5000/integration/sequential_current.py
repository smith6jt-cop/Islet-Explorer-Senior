"""Read-only contracts for an existing, reviewed Islet-only result package.

This verifies saved results and their inputs; it never reruns assignments or
creates scientific judgments. Hashing figures does not independently validate
their anatomy or repeat raw image decoding.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


CORE_FILES = (
    'all_islet_relationships.parquet', 'comparable_islet_relationships.parquet',
    'not_comparable_islet_relationships.parquet', 'all_image_reviews.parquet',
    'all_donor_comparability.csv', 'donor_inventory.csv',
    'comparable_islet_geometries.parquet', 'comparability_summary.json',
    'analysis_plan.json', 'workflow_binding.json', 'channel_sources.json',
    'REPORT.html', 'display_revision.json', 'final_verification.json',
    'numeric/independent_numeric_checks.json',
)


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(path):
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    _require((before.st_ino, before.st_size, before.st_mtime_ns) ==
             (after.st_ino, after.st_size, after.st_mtime_ns),
             f'File changed during verification: {path}')
    return digest.hexdigest()


def _local(root, relative):
    path = Path(relative)
    _require(not path.is_absolute() and '..' not in path.parts,
             f'Unsafe package-relative path: {relative}')
    path = (root / path).resolve()
    _require(path.is_relative_to(root), f'Path escapes result package: {relative}')
    return path


def _keys(frame):
    return set(frame[['donor_id', 'moving_id', 'fixed_id']].astype(str)
               .itertuples(index=False, name=None))


def verify_current_islets(root, *, expected_manifest_sha256=None,
                          verify_images=True):
    """Return a receipt, or raise ValueError on saved-result contract failure.

    Requires the pinned native object tables and transforms on this host.
    Native tables are read only for native IDs, categories and row counts; their
    complete byte content and original multigigabyte images are not rehashed.
    Figure hashes are checked only when ``verify_images`` is True. The artifact
    manifest must itself be pinned by the caller for authenticated reuse.
    """
    import pandas as pd
    import pyarrow.parquet as pq

    root = Path(root).expanduser().resolve()
    manifest_path = root / 'artifact_manifest.json'
    manifest_hash = _sha(manifest_path)
    if expected_manifest_sha256 is not None:
        _require(manifest_hash == expected_manifest_sha256,
                 'Current Islet manifest differs from the selected recipe')
    rows = json.loads(manifest_path.read_text())
    _require(isinstance(rows, list), 'Artifact manifest must be a list')
    manifest = {}
    for row in rows:
        relative = row['path']
        _local(root, relative)
        _require(relative not in manifest, f'Duplicate manifest path: {relative}')
        manifest[relative] = row
    verified_files = {}

    def verify(relative):
        if relative in verified_files:
            return verified_files[relative]
        _require(relative in manifest, f'Unmanifested required artifact: {relative}')
        path = _local(root, relative)
        entry = manifest[relative]
        _require(path.is_file(), f'Missing required artifact: {relative}')
        _require(path.stat().st_size == entry['bytes'], f'Artifact size mismatch: {relative}')
        digest = _sha(path)
        _require(digest == entry['sha256'], f'Artifact hash mismatch: {relative}')
        verified_files[relative] = {'path': relative, 'bytes': entry['bytes'], 'sha256': digest}
        return verified_files[relative]

    for relative in CORE_FILES:
        verify(relative)
    summary = json.loads((root / 'comparability_summary.json').read_text())
    display = json.loads((root / 'display_revision.json').read_text())
    display_directory = display['directory']
    _local(root, display_directory)
    pairs = pd.read_parquet(root / 'all_islet_relationships.parquet')
    comparable = pd.read_parquet(root / 'comparable_islet_relationships.parquet')
    excluded = pd.read_parquet(root / 'not_comparable_islet_relationships.parquet')
    reviews = pd.read_parquet(root / 'all_image_reviews.parquet')
    geometry = pd.read_parquet(root / 'comparable_islet_geometries.parquet')
    donors = pd.read_csv(root / 'all_donor_comparability.csv', dtype={'donor_id': str, 'source_serial': str})
    inventory = pd.read_csv(root / 'donor_inventory.csv', dtype={'donor_id': str, 'source_serial': str})
    _require(len(pairs) == summary['candidate_islet_pairs'], 'Candidate count differs from summary')
    _require(len(comparable) == summary['comparable'], 'Comparable count differs from summary')
    _require(len(excluded) == summary['not_comparable'], 'Not-comparable count differs from summary')
    _require(len(reviews) == summary['image_reviewed_pairs'], 'Image-review count differs from summary')
    _require(len(donors) == summary['inventory_donors'], 'Donor count differs from summary')
    _require(not pairs.pair_id.duplicated().any() and len(_keys(pairs)) == len(pairs),
             'Duplicate current pair identifiers')
    _require(not pairs.duplicated(['donor_id', 'moving_id']).any() and
             not pairs.duplicated(['donor_id', 'fixed_id']).any(), 'Assignment is not one to one')
    _require(pairs.RNA_object_class.eq('Islet').all() and pairs.protein_object_class.eq('Islet').all(),
             'Current correspondence includes a non-Islet object')
    _require(pairs.correspondence_policy.eq('native_islet_only_before_assignment_no_cluster_or_single_v1').all(),
             'Unexpected correspondence policy')
    _require(not pairs.same_cell_identity_claim.any() and not pairs.comparability_is_expert_ground_truth.any(),
             'Unsupported identity or ground-truth claim')
    accepted_name = 'Comparable as the same islet'
    excluded_name = 'Not comparable as the same islet'
    _require(pairs.comparability_designation.isin([accepted_name, excluded_name]).all(), 'Unknown designation')
    try:
        pd.testing.assert_frame_equal(comparable.reset_index(drop=True),
            pairs[pairs.comparability_designation.eq(accepted_name)].reset_index(drop=True))
        pd.testing.assert_frame_equal(excluded.reset_index(drop=True),
            pairs[pairs.comparability_designation.eq(excluded_name)].reset_index(drop=True))
    except AssertionError as exc:
        raise ValueError('Comparable/not-comparable exports differ from the full table') from exc
    _require(not reviews.pair_id.duplicated().any() and
             set(reviews.pair_id) == set(pairs.loc[pairs.eligible_for_visual_review, 'pair_id']),
             'Missing or duplicate recorded image reviews')
    _require(set(comparable.pair_id) <= set(reviews.pair_id), 'Comparable pair lacks recorded image evidence')
    _require(not donors.donor_id.duplicated().any() and not inventory.donor_id.duplicated().any() and
             set(donors.donor_id) == set(inventory.donor_id), 'Donor inventory mismatch')
    unavailable = donors.assessment_status.eq('Not assessable')
    _require(int(unavailable.sum()) == summary['not_assessable_donors'], 'Not-assessable count mismatch')
    _require(donors.loc[unavailable, ['candidate_pairs', 'comparable_as_same_islet',
                 'not_comparable_as_same_islet']].isna().all().all(),
             'Not-assessable donors have misleading zero or negative pair counts')
    _require(set(pairs.donor_id.astype(str)) == set(donors.loc[~unavailable, 'donor_id']),
             'Paired donors differ from roster')
    sources = {str(row['donor_id']): row for row in json.loads((root / 'channel_sources.json').read_text())}
    native_sources = []
    transforms = []
    for donor, group in pairs.groupby('donor_id'):
        donor = str(donor)
        record = inventory.set_index('donor_id').loc[donor]
        source = sources[donor]
        transform = Path(source['transform_path'])
        _require(_sha(transform) == source['transform_sha256'], f'Transform changed for donor {donor}')
        transforms.append({'donor_id': donor, 'path': str(transform), 'sha256': source['transform_sha256']})
        for role, id_column, path_column in [('RNA', 'moving_id', 'RNA_native_objects_path'),
                                             ('protein', 'fixed_id', 'protein_native_objects_path')]:
            path = Path(record[path_column])
            objects = pd.read_parquet(path, columns=['structure_id', 'object_class'])
            _require(not objects.structure_id.duplicated().any(), f'Duplicate native IDs: {path}')
            islet_ids = set(objects.loc[objects.object_class.eq('Islet'), 'structure_id'])
            _require(set(group[id_column]) <= islet_ids, f'Current IDs absent from native Islets: {donor} {role}')
            native_sources.append({'donor_id': donor, 'role': role, 'path': str(path),
                                   'rows': pq.read_metadata(path).num_rows, 'native_islets': len(islet_ids),
                                   'full_content_hash_checked': False})
    _require(len(geometry) == 2 * len(comparable) == summary['comparable_geometry_rows'], 'Geometry count mismatch')
    _require(geometry.object_class.eq('Islet').all() and not geometry.same_cell_identity_claim.any(),
             'Unexpected geometry class or identity claim')
    _require(not geometry.duplicated(['pair_id', 'role']).any(), 'Duplicate pair/role geometry')
    for pair in comparable.itertuples(index=False):
        rows = geometry[geometry.pair_id.eq(pair.pair_id)].set_index('role')
        _require(set(rows.index) == {'RNA', 'protein'}, f'Missing pair geometry: {pair.pair_id}')
        _require(rows.loc['RNA', 'structure_id'] == pair.moving_id and
                 rows.loc['protein', 'structure_id'] == pair.fixed_id, f'Geometry/native ID mismatch: {pair.pair_id}')
        _require(rows.coordinate_frame.eq(f'phenocycler:{pair.donor_id}:um').all(),
                 f'Unexpected coordinate frame: {pair.pair_id}')
        _require(rows.loc['RNA', 'transform_sha256'] == sources[str(pair.donor_id)]['transform_sha256'],
                 f'Geometry transform pin mismatch: {pair.pair_id}')
    image_records = []
    pair_index = pairs.set_index('pair_id')
    for review in reviews.itertuples(index=False):
        pair_id = str(review.pair_id)
        _require(re.fullmatch(r'[0-9]+-I[0-9]+', pair_id) is not None, f'Unsafe pair ID: {pair_id}')
        reviewed = f'pair_figures/{pair_id}.jpg'
        current = f'{display_directory}/figures/{pair_id}.jpg'
        current_metadata = f'{display_directory}/figures/{pair_id}.json'
        _require(reviewed in manifest and manifest[reviewed]['sha256'] == review.image_sha256,
                 f'Recorded review image is not pinned: {pair_id}')
        _require(current in manifest and current_metadata in manifest,
                 f'Current display image is not pinned: {pair_id}')
        verify(current_metadata)
        metadata = json.loads((root / current_metadata).read_text())
        pair = pair_index.loc[pair_id]
        _require(metadata['pair_id'] == pair_id and metadata['moving_id'] == pair.moving_id and
                 metadata['fixed_id'] == pair.fixed_id and metadata['designation'] == pair.comparability_designation,
                 f'Current display metadata differs from result: {pair_id}')
        if verify_images:
            verify(reviewed)
            verify(current)
        image_records.append({'pair_id': pair_id, 'reviewed_image': reviewed,
                              'reviewed_sha256': manifest[reviewed]['sha256'],
                              'current_image': current, 'current_sha256': manifest[current]['sha256'],
                              'current_metadata': current_metadata})
    return {'status': 'PASS', 'root': str(root), 'manifest_sha256': manifest_hash,
            'counts': {'candidates': len(pairs), 'comparable': len(comparable), 'not_comparable': len(excluded),
                       'image_reviews': len(reviews), 'inventory_donors': len(donors),
                       'paired_donors': int((~unavailable).sum()), 'not_assessable_donors': int(unavailable.sum()),
                       'geometry_rows': len(geometry)},
            'display_directory': display_directory, 'images_verified': bool(verify_images),
            'image_records': image_records, 'verified_files': list(verified_files.values()),
            'transforms': transforms, 'native_sources': native_sources,
            'limitations': ['Saved-result contract verification; no numerical or anatomical reassessment',
                            'No full native-table or original TIFF content hashing',
                            'No same-cell identity or expert ground-truth inference']}
