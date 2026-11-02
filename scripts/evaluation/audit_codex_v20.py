#!/usr/bin/env python3
"""Check v20 bindings independently of the older CODEX-port parity prerequisite."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--codex-root', type=Path, default=Path('/data/panc_CODEX/Islet-Explorer'))
    parser.add_argument('--out', type=Path, required=True, help='New output root; historical audits are never replaced')
    parser.add_argument('--recovery-manifest', type=Path, help='Explicit content-bound relocation map')
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    root = args.codex_root
    source = root / "data/development/uniform-20donor-20261030-v20/uniform_release.json"
    release = json.loads(source.read_text())
    recovered = {}
    if args.recovery_manifest:
        recovery = json.loads(args.recovery_manifest.read_text())
        actual_release = hashlib.sha256(source.read_bytes()).hexdigest()
        if recovery['release_sha256'] != actual_release:
            raise ValueError('Recovery map describes a different release manifest')
        if len({r['declared_path'] for r in recovery['bindings']}) != len(recovery['bindings']):
            raise ValueError('Recovery map has duplicate source paths')
        recovered = {r['declared_path']: r for r in recovery['bindings']}
    reconciliation = root / "output/typing_20261030/final_reconciliation_v20/reconciliation.json"
    historical = json.loads(reconciliation.read_text())
    for path in (source, reconciliation):
        shutil.copyfile(path, out / path.name)
    paths = [(Path(release["run"]) / "typing/levels.csv", release["levels_sha256"]),
             (Path(release["donor_metadata"]), release["donor_metadata_sha256"])]
    binding = release["structure_generation"]
    paths += [(Path(binding[k + "_path"]), binding[k + "_sha256"]) for k in ("binding", "manifest")]
    paths += [(Path(release["run"]) / f"typing/identity/{name}.parquet", sha)
              for name, sha in release["identity_sha256"].items()]
    checks = []
    for path, expected in paths:
        declared = path
        remapping = recovered.get(str(path))
        if args.recovery_manifest:
            if remapping is None or remapping['status'] != 'recovered' or remapping['expected_sha256'] != expected:
                raise ValueError(f'Recovery map lacks the exact release binding: {path}')
            path = Path(remapping['recovered_path'])
        if path.is_file():
            with path.open("rb") as handle:
                actual = hashlib.file_digest(handle, "sha256").hexdigest()
        else:
            actual = None
        checks.append({"path": str(path), "declared_path": str(declared), "recovered": bool(remapping),
                       "expected": expected, "actual": actual,
                       "status": "pass" if actual == expected else "fail"})
    summary = {"producer": release["git_commit"], "release_status": release["status"],
               "source_release": str(source), "checks": checks,
               "declared_reconciliation_status": historical["status"],
               "declared_reconciled_images": len(historical["images"]),
               "declared_cells_analysed": release["cells_analysed"],
               "recovery_manifest": str(args.recovery_manifest) if args.recovery_manifest else None,
               "source_release_modified": False,
               "validation_scope": "Checked release levels, 20 identity tables, structure binding/selection manifest, and donor metadata hashes, using an explicit content-bound relocation map when supplied. Historical reconciliation preserved, not reexecuted. This does not establish full v20 execution-input closure, biological acceptance, or release activation."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"v20: {len(checks)} hash checks; {sum(r['status'] == 'fail' for r in checks)} mismatches")
    return int(any(row['status'] == 'fail' for row in checks))


if __name__ == "__main__":
    raise SystemExit(main())
