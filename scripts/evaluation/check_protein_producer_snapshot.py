#!/usr/bin/env python3
"""Compare the protein producers and their local import closure with a launch snapshot.

The snapshot is deliberately broader than this route. This audit does not pretend that
every archived file was executed, or that a static import closure is a runtime trace.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path, help="source snapshot for the exact attempt")
    parser.add_argument("--stage", type=Path, help="explicit complete or incomplete Stage 06 directory")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    snapshot_path = args.snapshot or args.run / "producer_snapshot.json"
    snapshot = json.loads(snapshot_path.read_text())
    frozen = {row["path"]: row["sha256"] for row in snapshot["sources"]}
    package = repo / "core/multiplex5000"
    script_roots = [package / "scripts/celltype", package / "scripts/typing"]
    stage = args.stage or args.run / "releases/protein-availability-candidate/06_typing"
    if not stage.exists():
        stages = list((stage.parent).glob(".06_typing.incomplete-*"))
        assert len(stages) == 1, "specify a run with one complete or incomplete candidate"
        stage = stages[0]
    run_path = stage / "pancreas/run.json"
    runner = json.loads(run_path.read_text())
    seeds = [Path(path) for path in runner["source_scripts"]]
    seeds += [package / "multiplex5000/stages/io60_mixed_protein_pilot.py",
              package / "multiplex5000/stages/typing.py"]
    seen = set()
    queue = list(seeds)

    def resolve(module, origin):
        roots = [package, origin.parent, *script_roots]
        for root in roots:
            target = root.joinpath(*module.split("."))
            for candidate in (target.with_suffix(".py"), target / "__init__.py"):
                if candidate.is_file() and candidate.is_relative_to(package):
                    return candidate.resolve()
        return None

    while queue:
        path = queue.pop()
        if path in seen:
            continue
        seen.add(path)
        for parent in path.parents:
            if parent == package:
                break
            init = parent / "__init__.py"
            if init.is_file():
                queue.append(init)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                modules = [item.name for item in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    relative = path.parent.relative_to(package).parts
                    prefix = ".".join(relative[:len(relative) - node.level + 1])
                    base = ".".join(part for part in (prefix, node.module) if part)
                else:
                    base = node.module or ""
                modules = [base, *[f"{base}.{item.name}" for item in node.names]]
            else:
                continue
            # The repository's documented compatibility package resolves these
            # modules to multiplex5000 through a meta-path import hook.
            modules += ["multiplex5000" + module[len("phenocycler"):]
                        for module in modules if module.startswith("phenocycler.")]
            for module in modules:
                target = resolve(module, path)
                # The Stage 06 dispatcher selects this route explicitly. Imports in
                # other route branches are not part of this candidate's execution.
                if (path.name == "typing.py" and path.parent.name == "stages"
                        and target is not None and target.parent == path.parent
                        and target.name not in {"__init__.py", "io60_mixed_protein_pilot.py"}):
                    continue
                if target is not None and target not in seen:
                    queue.append(target)
    rows = []
    for path in sorted(seen):
        relative = str(path.relative_to(repo))
        current = sha(path)
        rows.append({"path": relative, "launch_sha256": frozen.get(relative),
                     "current_sha256": current,
                     "unchanged": frozen[relative] == current if relative in frozen else None,
                     "status": ("not_snapshotted" if relative not in frozen else
                                "unchanged" if frozen[relative] == current else "changed")})
    closure_names = {row["path"] for row in rows}
    unrelated_changes = []
    for path, digest in frozen.items():
        if path not in closure_names and (repo / path).is_file() and sha(repo / path) != digest:
            unrelated_changes.append(path)
    runner_checks = [{"path": path, "launch_sha256": digest,
                      "current_sha256": sha(Path(path)), "unchanged": sha(Path(path)) == digest}
                     for path, digest in runner["source_scripts"].items()]
    result = {
        "status": ("fail" if any(row["unchanged"] is False for row in rows + runner_checks)
                   else "pass_with_unavailable_static_candidates"
                   if any(row["unchanged"] is None for row in rows) else "pass"),
        "scope": "declared protein producers and conservative static local-import closure; not a runtime trace",
        "limitations": "Conditional imports may overinclude files. File comparison proves equality at launch and audit; it does not continuously monitor intermediate writes. Dataset and registry inputs are checked by stage receipts separately.",
        "snapshot": {"path": str(snapshot_path), "sha256": sha(snapshot_path), "files": len(frozen)},
        "runner_manifest": {"path": str(run_path), "sha256": sha(run_path)},
        "closure": rows, "runner_source_checks": runner_checks,
        "snapshot_changes_outside_protein_closure": unrelated_changes,
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(result["status"], len(rows), "protein closure files;", len(unrelated_changes), "changes outside closure")
    if result["status"] == "fail":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
