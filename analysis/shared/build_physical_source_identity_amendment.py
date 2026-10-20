#!/usr/bin/env python3
"""Create a separate, source-only contract for the confirmed PAN6523 scan-name error."""
import argparse
import json
from pathlib import Path
from typing import Mapping

try:
    from . import physical_source_identity_amendment_view as view
except ImportError:
    import physical_source_identity_amendment_view as view


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    base = view.semantic.load_semantic_selection(view.ENTRY)
    s = base._scope(view.SCOPE_ID)
    contract = {'schema_version': 1, 'kind': 'physical_source_identity_amendment',
        'scope_id': view.SCOPE_ID, 'capture': view.CAPTURE, 'policy': view.POLICY,
        'base_selection': view.pin(view.ENTRY), 'base_source': {'path': s['path'], 'sha256': s['sha256']},
        'base_selection_api': view.pin(view.semantic.__file__),
        'base_semantic_api': view.pin(view.protein.__file__),
        'base_source_api': view.pin(view.semantic.base.__file__),
        'amendment_api': view.pin(view.__file__), 'builder': view.pin(__file__),
        'human_attestation': view.pin(view.ATTESTATION),
        'current_raw_header_metadata_readback': view.pin(view.HEADER_READBACK)}
    args.output_dir.mkdir(parents=True, exist_ok=False)
    path = args.output_dir / 'physical_source_identity_amendment.json'
    with path.open('x') as stream:
        json.dump(contract, stream, indent=2, allow_nan=False)
        stream.write('\n')
    selected = view.load_source_identity_amendment(path)
    with (args.output_dir / 'terminal_compact_source_identity_amendment_validation.json').open('x') as stream:
        json.dump(selected.summary(), stream, indent=2, allow_nan=False,
                  default=lambda x: dict(x) if isinstance(x, Mapping) else x)
        stream.write('\n')
    print(json.dumps(view.pin(path)))


if __name__ == '__main__':
    main()
