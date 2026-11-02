#!/usr/bin/env python3
"""Check v20 bindings independently of the older CODEX-port parity prerequisite."""
import hashlib
import json
import shutil
from pathlib import Path


def main():
    base = Path(__file__).resolve().parents[2]
    out = base / "outputs/dynamic_portability_20261102/codex_v20"
    out.mkdir(parents=True, exist_ok=True)
    root = Path("/data/panc_CODEX/Islet-Explorer")
    source = root / "data/development/uniform-20donor-20261030-v20/uniform_release.json"
    release = json.loads(source.read_text())
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
        if path.is_file():
            with path.open("rb") as handle:
                actual = hashlib.file_digest(handle, "sha256").hexdigest()
        else:
            actual = None
        checks.append({"path": str(path), "expected": expected, "actual": actual,
                       "status": "pass" if actual == expected else "fail"})
    summary = {"producer": release["git_commit"], "release_status": release["status"],
               "source_release": str(source), "checks": checks,
               "declared_reconciliation_status": historical["status"],
               "declared_reconciled_images": len(historical["images"]),
               "declared_cells_analysed": release["cells_analysed"],
               "validation_scope": "Checked release levels, 20 identity tables, structure binding/selection manifest, and donor metadata hashes against current files. Historical reconciliation preserved, not reexecuted. Port compatibility additionally requires older baseline parity."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"v20: {len(checks)} hash checks; {sum(r['status'] == 'fail' for r in checks)} mismatches")


if __name__ == "__main__":
    main()
