#!/usr/bin/env python3
"""Static parameter/default and registry inventory; does not claim dynamic execution coverage."""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path


def catalogue(path, category, producer):
    source = path.read_text()
    tree = ast.parse(source)
    sha = hashlib.sha256(source.encode()).hexdigest()
    rows = []
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for name in targets:
                if isinstance(name, ast.Name) and name.id.isupper() and node.value is not None:
                    uses = [n.lineno for n in ast.walk(tree) if isinstance(n, ast.Name)
                            and isinstance(n.ctx, ast.Load) and n.id == name.id]
                    rows.append({"name": name.id, "kind": "module_setting", "line": node.lineno,
                                 "expression": ast.unparse(node.value), "read_lines": sorted(set(uses))})
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            args = node.args.posonlyargs + node.args.args
            defaults = list(zip(args[-len(node.args.defaults):], node.args.defaults)) if node.args.defaults else []
            defaults += [(a, d) for a, d in zip(node.args.kwonlyargs, node.args.kw_defaults) if d is not None]
            for arg, default in defaults:
                rows.append({"name": f"{node.name}.{arg.arg}", "kind": "function_default",
                             "line": default.lineno, "expression": ast.unparse(default),
                             "read_lines": sorted({n.lineno for n in ast.walk(node) if isinstance(n, ast.Name)
                                                   and isinstance(n.ctx, ast.Load) and n.id == arg.arg})})
    return [{**row, "source": str(path), "sha256": sha, "producer": producer,
             "dependency_class": category, "execution_evidence": "static; see route receipts and dependency matrix"}
            for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    pkg = Path("core/multiplex5000/multiplex5000")
    groups = {
        "tissue_biology": ["tissues/pancreas_xenium_worker.py", "tissues/lymph_node_xenium_worker.py",
                           "tissues/pancreas_rna_naming.py", "tissues/pln_clustering.py", "tissues/protein_section.py",
                           "tissues/t_cell.py", "tissues/macrophage.py", "tissues/pancreas_mixed.py",
                           "../scripts/typing/subpanel.py", "../scripts/typing/refine.py",
                           "../scripts/typing/build_cells.py", "../scripts/typing/immune_subtype.py"],
        "assay_specific_measurement": ["protein_evidence.py", "xenium_typing_input.py", "maldi.py", "cohort_scale.py"],
        "registration_geometry": ["integration/register.py", "integration/sameslide.py", "integration/paired_typing.py",
                                  "stages/maldi_serial_registration.py", "stages/maldi_serial_integration.py",
                                  "stages/same_slide_registration.py", "stages/same_slide_integration.py"],
        "shared_computation": ["identity.py", "relationships.py", "tissues/evaluate.py"],
        "unresolved_assumption": ["stages/io60_mixed_protein_pilot.py", "stages/io60_rna_typing_pilot.py",
                                  "stages/spleen_typing_pilot.py", "stages/pln_typing_pilot.py"],
    }
    rows = []
    for category, paths in groups.items():
        for relative in paths:
            rows.extend(catalogue(pkg / relative, category, "maldi/serial-platform@00ec032"))
    registries = []
    for snapshot, producer in [("typing_port", "port@df3fe1d"), ("remote/spleen-code", "spleen@a9549cc")]:
        p = args.evaluation / snapshot / pkg
        for folder in (p / "positivity", p / "tissues"):
            for source in sorted(folder.glob("*.py")):
                if folder.name == "tissues" and source.name not in ("pancreas.py", "spleen.py"):
                    continue
                rows.extend(catalogue(source, "unresolved_assumption", producer))
        for table in sorted((p / "positivity/registry").glob("*.csv")):
            registries.append({"source": str(table), "producer": producer,
                               "sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
                               "dependency_class": "panel_configuration",
                               "rows": list(csv.DictReader(table.open()))})
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "parameter_catalog.json").write_text(json.dumps({"scope": "static constants, function defaults and registry rows of the audited route modules; not a runtime trace", "parameters": rows, "registries": registries}, indent=2) + "\n")
    print(len(rows), "settings/defaults;", len(registries), "registry tables")


if __name__ == "__main__":
    main()
