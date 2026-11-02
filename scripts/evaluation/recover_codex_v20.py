#!/usr/bin/env python3
"""Recover release-bound v20 files by content without modifying the source run."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def recover(root: Path, out: Path):
    release_path = root / 'data/development/uniform-20donor-20261030-v20/uniform_release.json'
    release = json.loads(release_path.read_text())
    alternate = root / 'output/v19_v20_component_comparison_20261102/endocrine/authoritative_run'
    binding = release['structure_generation']
    items = [('typing/levels.csv', Path(release['run']) / 'typing/levels.csv',
              alternate / 'typing/levels.csv', release['levels_sha256']),
             ('inputs/donor_metadata.csv', Path(release['donor_metadata']),
              Path(release['donor_metadata']), release['donor_metadata_sha256']),
             ('structure_generation.json', Path(binding['binding_path']),
              Path(binding['manifest_path']).with_name('binding_candidate.json'), binding['binding_sha256']),
             ('inputs/structure_selection_manifest.json', Path(binding['manifest_path']),
              Path(binding['manifest_path']), binding['manifest_sha256'])]
    items += [(f'typing/identity/{tag}.parquet', Path(release['run']) / f'typing/identity/{tag}.parquet',
               alternate / f'typing/identity/{tag}.parquet', expected)
              for tag, expected in release['identity_sha256'].items()]
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for relative, declared, candidate, expected in items:
        original = digest(declared) if declared.is_file() else None
        candidates = list(dict.fromkeys([declared, candidate]))
        selected = next((p for p in candidates if p.is_file() and digest(p) == expected), None)
        row = {'declared_path': str(declared), 'declared_current_sha256': original,
               'expected_sha256': expected, 'candidates_checked': [str(p) for p in candidates],
               'recovered_path': None, 'status': 'unavailable'}
        if selected is not None:
            target = out / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(selected, target)
            if digest(target) != expected:
                raise ValueError(f'Recovered content changed during copy: {selected}')
            row.update(recovered_path=str(target.resolve()), selected_source=str(selected),
                       recovered_sha256=expected, status='recovered')
        rows.append(row)
    shutil.copy2(release_path, out / 'uniform_release.json')
    summary = {'schema': 'codex.v20_content_recovery/1', 'release': str(release_path),
               'release_sha256': digest(release_path), 'producer': release['git_commit'],
               'release_status': release['status'], 'bindings': rows,
               'expected': len(items), 'recovered': sum(r['status'] == 'recovered' for r in rows),
               'source_run_modified': False,
               'scope': 'Exact recovery of the 24 audited release bindings. No biological acceptance, activation, or claim that the full v20 execution input closure has been recovered.'}
    (out / 'recovery.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f"v20 bindings: {summary['recovered']}/{summary['expected']} recovered", flush=True)
    return int(summary['recovered'] != summary['expected'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--codex-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    raise SystemExit(recover(args.codex_root, args.out))
