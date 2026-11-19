#!/usr/bin/env python3
"""Evidence-bounded static parameter traces; never execute or change scientific routes.

This is a source inventory, not a profiler. A receipt authenticates a producer and
its declared settings; it does not establish execution of every function/branch.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import shlex
import subprocess
import tarfile
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

PKG = Path("core/multiplex5000/multiplex5000")
RUN = Path("outputs/dynamic_portability_20261102")
DOC = Path("docs/evaluations/dynamic_portability_20261102")
FUNCTIONS = (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)
SCOPES = (*FUNCTIONS, ast.ClassDef, ast.Module)
WRITERS = {"to_csv", "to_parquet", "write_table", "write_text", "write_bytes", "dump",
           "dumps", "save", "savez", "savetxt", "writerow", "writerows", "write_h5ad"}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def normalized(path):
    return str(Path(path).resolve().relative_to(Path.cwd()))


def repo_path(path):
    parts = Path(normalized(path)).parts
    return str(Path(*parts[parts.index("core"):]))


def namespace(path):
    path = normalized(path)
    if "/typing_port/" in path:
        return "port@df3fe1d"
    if "/remote/spleen-code/" in path:
        return "spleen@a9549cc"
    return "catalog_working_tree"


def module_name(path):
    path = repo_path(path)
    if path.startswith(str(PKG) + "/"):
        path = "multiplex5000/" + path[len(str(PKG)) + 1:]
    else:
        path = path.removeprefix("core/multiplex5000/")
    return path.removesuffix(".py").removesuffix("/__init__").replace("/", ".")


class Source:
    def __init__(self, path, ns=None):
        self.path = normalized(path)
        self.ns = ns or namespace(path)
        self.module = module_name(path)
        self.bytes = Path(path).read_bytes()
        self.text = self.bytes.decode()
        self.lines = self.text.splitlines(keepends=True)
        self.hash = sha(self.bytes)
        self.tree = ast.parse(self.text)
        self.nodes = list(ast.walk(self.tree))
        self.parent = {child: node for node in ast.walk(self.tree)
                       for child in ast.iter_child_nodes(node)}
        self.bindings = defaultdict(dict)
        self.globals = defaultdict(set)
        for node in ast.walk(self.tree):
            scope = self.scope(node)
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        self.bindings[scope][alias.asname or alias.name.split(".")[0]] = (
                            alias.name if alias.asname else alias.name.split(".")[0])
                else:
                    base = node.module or ""
                    if node.level:
                        package = self.module.split(".")
                        if not self.path.endswith("/__init__.py"):
                            package = package[:-1]
                        base = ".".join(package[:len(package) - node.level + 1] + ([base] if base else []))
                    for alias in node.names:
                        self.bindings[scope][alias.asname or alias.name] = base + "." + alias.name
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                self.bindings[scope][node.name] = self.qualified(node)
            elif isinstance(node, ast.Global):
                self.globals[scope].update(node.names)
            elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
                self.bindings[scope].setdefault(node.id, self.scope_name(scope) + "." + node.id)
        for scope in ast.walk(self.tree):
            if isinstance(scope, FUNCTIONS):
                for arg in (*scope.args.posonlyargs, *scope.args.args, *scope.args.kwonlyargs,
                            *([scope.args.vararg] if scope.args.vararg else []),
                            *([scope.args.kwarg] if scope.args.kwarg else [])):
                    self.bindings[scope][arg.arg] = self.qualified(scope) + "." + arg.arg
        self.read_index = defaultdict(list)
        self.call_index = defaultdict(list)
        for node in self.nodes:
            if isinstance(node, (ast.Name, ast.Attribute)) and (not isinstance(node, ast.Name)
                                                               or isinstance(node.ctx, ast.Load)):
                self.read_index[self.resolve_expr(node)].append(node)
            elif isinstance(node, ast.Call):
                self.call_index[self.resolve_expr(node.func)].append(node)

    def ancestors(self, node):
        while node in self.parent:
            node = self.parent[node]
            yield node

    def scope(self, node):
        return next((a for a in self.ancestors(node) if isinstance(a, SCOPES)), self.tree)

    def scope_name(self, node):
        return self.module if isinstance(node, ast.Module) else self.qualified(node)

    def qualified(self, node):
        scopes = [a.name for a in self.ancestors(node)
                  if isinstance(a, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))]
        return ".".join([self.module, *reversed(scopes), getattr(node, "name", "<lambda>")])

    def resolve_name(self, node, name):
        scope = self.scope(node)
        if isinstance(scope, FUNCTIONS) and any(node in ast.walk(d) for d in
                [*scope.args.defaults, *[d for d in scope.args.kw_defaults if d is not None]]):
            scope = self.scope(scope)  # Defaults capture the defining scope, not arguments.
        while True:
            if name in self.globals[scope]:
                scope = self.tree
            if name in self.bindings[scope]:
                return self.bindings[scope][name]
            if isinstance(scope, ast.Module):
                return self.module + "." + name
            scope = self.scope(scope)

    def resolve_expr(self, node):
        if isinstance(node, ast.Name):
            return self.resolve_name(node, node.id)
        if isinstance(node, ast.Attribute):
            base = self.resolve_expr(node.value)
            return base + "." + node.attr if base else None
        return None

    def ref(self, node):
        start, end = node.lineno - 1, node.end_lineno - 1
        # AST column offsets are UTF-8 byte offsets.
        lines = [line.encode() for line in self.lines[start:end + 1]]
        if len(lines) == 1:
            text = lines[0][node.col_offset:node.end_col_offset].decode()
        else:
            text = b"".join([lines[0][node.col_offset:], *lines[1:-1],
                              lines[-1][:node.end_col_offset]]).decode()
        return {"source": self.path, "sha256": self.hash, "line": node.lineno,
                "end_line": getattr(node, "end_lineno", node.lineno),
                "scope": self.scope_name(self.scope(node)), "expression": text[:700],
                **({"expression_truncated": True, "expression_sha256": sha(text.encode())}
                   if len(text) > 700 else {})}

    def declaration(self, row):
        if row["kind"] == "module_setting":
            for node in self.tree.body:
                if isinstance(node, (ast.Assign, ast.AnnAssign)) and node.lineno == row["line"]:
                    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                    if any(isinstance(n, ast.Name) and n.id == row["name"] for n in targets):
                        return node, node.value, self.module + "." + row["name"]
        else:
            name, argname = row["name"].rsplit(".", 1)
            for node in ast.walk(self.tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name:
                    args = node.args.posonlyargs + node.args.args
                    pairs = list(zip(args[-len(node.args.defaults):], node.args.defaults)) if node.args.defaults else []
                    pairs += list(zip(node.args.kwonlyargs, node.args.kw_defaults))
                    for arg, default in pairs:
                        if default is not None and arg.arg == argname and default.lineno == row["line"]:
                            return node, default, self.qualified(node) + "." + argname
        raise ValueError(f"exact declaration not found: {row}")


def direct_consumers(source, symbol, sources):
    found = []
    for src in sources:
        if src.ns != source.ns:
            continue
        for node in src.read_index.get(symbol, []):
            # Attribute bases are part of the same read; retain the exact parameter reference.
            stmt = next((p for p in src.ancestors(node) if isinstance(p, ast.stmt)), node)
            enclosing = src.scope(node)
            default_of = next((p for p in src.ancestors(node) if isinstance(p, FUNCTIONS)), None)
            in_default = bool(default_of and any(node in ast.walk(d)
                              for d in [*default_of.args.defaults,
                                        *[x for x in default_of.args.kw_defaults if x is not None]]))
            found.append({**src.ref(node), "statement": src.ref(stmt),
                          "read_phase": "function_definition_default_capture" if in_default else "body_or_module_read",
                          "consumer": src.scope_name(enclosing),
                          "resolution": "lexical_name_or_explicit_import; static_only"})
    return found


def output_candidates(source, symbol, sources, consumers):
    """Conservative syntactic assignment closure, separately labelled from runtime dataflow."""
    by_scope = defaultdict(list)
    lookup = {s.path: s for s in sources}
    for read in consumers:
        by_scope[(read["source"], read["consumer"])].append(read)
    outputs = []
    for path, scope_name in by_scope:
        src = lookup[path]
        scope = next((s for s in ast.walk(src.tree) if isinstance(s, SCOPES)
                      and src.scope_name(s) == scope_name), None)
        if scope is None:
            continue
        nodes = [n for n in ast.walk(scope) if src.scope(n) is scope]
        tainted = {symbol}

        def uses(node, src=src, tainted=tainted):
            return node is not None and any(src.resolve_expr(n) in tainted
                                           for n in ast.walk(node) if isinstance(n, (ast.Name, ast.Attribute)))

        assignments = []
        for _ in range(len(nodes) + 1):
            before = len(tainted)
            for node in nodes:
                targets, value = [], None
                if isinstance(node, ast.Assign):
                    targets, value = node.targets, node.value
                elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
                    targets, value = [node.target], node.value
                elif isinstance(node, (ast.For, ast.AsyncFor)):
                    targets, value = [node.target], node.iter
                if uses(value):
                    for target in targets:
                        for name in ast.walk(target):
                            if isinstance(name, ast.Name) and isinstance(name.ctx, ast.Store):
                                resolved = src.resolve_name(name, name.id)
                                if resolved not in tainted:
                                    tainted.add(resolved)
                                    assignments.append(src.ref(node))
            if len(tainted) == before:
                break
        boundaries = []
        for node in nodes:
            if isinstance(node, ast.Return) and uses(node.value):
                boundaries.append({"kind": "return_to_caller", **src.ref(node)})
            elif isinstance(node, ast.Call) and (any(uses(a) for a in node.args)
                                                or any(uses(k.value) for k in node.keywords)):
                callee = ast.unparse(node.func)
                boundaries.append({"kind": "write_call" if callee.split(".")[-1] in WRITERS
                                   else "call_argument_boundary", "callee": callee, **src.ref(node)})
        outputs.append({"consumer_source": path, "consumer": scope_name,
                        "method": "flow_insensitive_syntactic_assignment_closure; possible_not_proven",
                        "assignment_edges": assignments, "boundaries": boundaries,
                        "disposition": "static_boundary_candidates_only" if boundaries else "output_unresolved",
                        "limitations": "No branch execution, alias/mutation analysis, return-value interprocedural propagation, or artifact attribution."})
    return outputs


def callsites(source, function, argument, sources):
    target = source.qualified(function)
    positional = [a.arg for a in function.args.posonlyargs + function.args.args]
    found = []
    for src in sources:
        if src.ns != source.ns:
            continue
        for node in src.call_index.get(target, []):
            supplied = next((k.value for k in node.keywords if k.arg == argument), None)
            if supplied is not None:
                disposition = "explicit_keyword_override"
            elif argument in positional and positional.index(argument) < len(node.args) and not any(
                    isinstance(n, ast.Starred) for n in node.args[:positional.index(argument) + 1]):
                supplied = node.args[positional.index(argument)]
                disposition = "explicit_positional_override"
            elif any(isinstance(n, ast.Starred) for n in node.args) or any(k.arg is None for k in node.keywords):
                disposition = "unresolved_star_arguments"
            else:
                disposition = "default_selected_if_call_executes"
            found.append({**src.ref(node), "argument": argument, "disposition": disposition,
                          "supplied_expression": ast.unparse(supplied) if supplied is not None else None})
    return found


def classification(row):
    name, source = row["name"], row["source"]
    # Existing catalog categories are source-scoped, not evidence of portability.
    category = row["dependency_class"]
    basis = "catalog source-scope classification; no biological transfer approval"
    if name in {"TARGET_SUM", "DETECTION_MIN_LOGNORM", "DETECTION_PREDICATE", "CONTROL_PREFIXES",
                "QC_MIN_COUNTS", "QC_MIN_FEATURES", "CLUSTER_MIN_COUNTS", "CLUSTER_MIN_FEATURES",
                "NUCLEUS_RULE_MEASURED", "NUCLEUS_RULE_UNAVAILABLE", "NUCLEUS_RULE_STATES"}:
        category, basis = "assay_specific_measurement", "RNA measurement, detection, or admission convention"
    elif name in {"IDENTITY_CORE", "IDENTITY_CORE_FLAT", "IDENTITY_CORE_PANC", "LINEAGE_CORES",
                  "IMMUNE_CORES", "NON_IMMUNE_CORES", "EMBEDDING_EXTRA", "PANEL", "PANELS", "RULES"}:
        category, basis = "panel_configuration", "explicit panel/core/rule source selection"
    elif name in {"N_HVG", "N_PCS", "PCA_FIT_MAX", "LEIDEN_SUBSAMPLE", "LEIDEN_RESOLUTION", "KNN_K",
                  "REFINE_MIN_CELLS", "REFINE_RESOLUTION", "CLUSTER_SUBSAMPLE", "CLUSTER_RESOLUTION",
                  "CLUSTER_SEED_BASE", "Z_SENSITIVITY", "NK_SENSITIVITY"}:
        category, basis = "shared_computation", "algorithm or sensitivity setting; its numerical transfer remains unvalidated"
    elif "/positivity/" in source:
        if source.endswith(("/tissue.py", "/sources.py")) or name in {
            "REFERENCE_COMPARTMENTS", "STRUCTURE_COLUMNS", "STRUCTURE_CONTEXT", "OVERLAP",
            "HORMONES_AND_LEUKOCYTES", "POPULATION_STAGE", "POPULATIONS", "NON_HAEM_COMPARTMENTS"}:
            category, basis = "tissue_biology", "reference/context/hook definition; execution must be checked separately"
        elif source.endswith("/coverage.py"):
            category, basis = "assay_specific_measurement", "coverage, compartment, and pixel measurement contract"
        elif source.endswith("/registry.py"):
            category, basis = "panel_configuration", "registry loading, measurement selection, or review vocabulary"
        elif source.endswith("/levels.py") and name in {"EDGES", "CEN", "SMOOTH_BINS", "HWHM_TO_SIGMA", "MAX_WIDTH", "AGAINST_ZERO", "PEAK_SIDE"}:
            category, basis = "shared_computation", "numerical fitting setting; reference and transfer validity separate"
    return {"class": category, "basis": basis, "catalog_class": row["dependency_class"]}


def evidence(path, pointer=None):
    path = Path(path)
    return {"path": str(path), "sha256": sha(path.read_bytes()),
            **({"json_pointer": pointer} if pointer is not None else {})}


def archive_producers(out):
    routes = {
        "rna": (RUN / "remote/rna/receipt.json", [
            "stages/io60_rna_typing_pilot.py", "tissues/pancreas_xenium_worker.py",
            "tissues/lymph_node_xenium_worker.py", "tissues/pancreas_rna_naming.py",
            "xenium_typing_input.py", "xenium_region.py", "xenium_sources.py"]),
        "pln": (RUN / "pln/receipt.json", ["stages/pln_typing_pilot.py", "protein_evidence.py",
                                                "tissues/protein_section.py", "tissues/evaluate.py",
                                                "tissues/lymph_node.py"]),
        "spleen": (RUN / "remote/spleen/receipt.json", ["stages/spleen_typing.py"]),
    }
    records, snapshots = [], []
    for route, (receiptpath, paths) in routes.items():
        receipt = json.loads(receiptpath.read_text())
        commit = receipt["code"]["commit"]
        for relative in paths:
            source = PKG / relative
            result = subprocess.run(["git", "show", f"{commit}:{source}"], capture_output=True, check=False)
            if result.returncode:
                records.append({"route": route, "commit": commit, "source": str(source),
                                "disposition": "git_object_unavailable", "stderr": result.stderr.decode()[:300]})
                continue
            target = out / "producer_sources" / commit / source
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(result.stdout)
            snapshots.append(Source(target, route + "@" + commit))
            current = sha(source.read_bytes()) if source.exists() else None
            records.append({"route": route, "commit": commit, "source": str(source),
                            "snapshot": str(target), "sha256": sha(result.stdout),
                            "current_sha256": current, "matches_current": current == sha(result.stdout),
                            "receipt": evidence(receiptpath, "/code")})
    return records, snapshots


def route_evidence(row, source, producers, primary):
    path = repo_path(source.path)
    found = []
    for record in producers:
        if record["source"] == path and "sha256" in record:
            found.append({"route": record["route"], "producer_source": record,
                          "source_matches_declaration": record["sha256"] == source.hash,
                          "execution_disposition": "route_receipt_exists; individual_consumer_execution_unrecorded",
                          "runtime_value": None})
    if source.ns == "spleen@a9549cc":
        receipt = RUN / "remote/spleen/receipt.json"
        commit = json.loads(receipt.read_text())["code"]["commit"]
        result = subprocess.run(["git", "show", f"{commit}:{path}"], capture_output=True, check=False)
        found.append({"route": "spleen", "receipt": evidence(receipt, "/code"),
                      "producer_source_commit": commit,
                      "producer_source_sha256": sha(result.stdout) if result.returncode == 0 else None,
                      "source_matches_declaration": result.returncode == 0 and sha(result.stdout) == source.hash,
                      "execution_disposition": "route_receipt_exists; module/function execution unrecorded (legacy identity and CODEX registries need not run)",
                      "runtime_value": None})
    for declared_path, expected in primary.get("source_scripts", {}).items():
        if Path(declared_path).resolve() == Path(source.path).resolve():
            found.append({"route": "protein_primary", "receipt": evidence(RUN / "protein_primary/pancreas/run.json", "/source_scripts"),
                          "producer_source_sha256": expected, "source_matches_declaration": expected == source.hash,
                          "execution_disposition": "partial_run: subpanel/build_cells complete; refine failed before final cells; subtype outputs not proven",
                          "runtime_value": None})
    if source.ns == "port@df3fe1d":
        completion = Path("outputs/handoff_completion_20261102/codex/completion_summary.json")
        found.append({"route": "CODEX_port", "evidence": evidence(completion),
                      "execution_disposition": "unchanged originating-source-to-port computation: all 60 exact comparisons pass across 20 images; all 20 historical identity assertions pass for 3,199,606 eligible cells. Original levels/canary default-CSV assertions retain 40 failures and the original mtime fixture still rejects 20 files. Controlled original-interpreter, single-thread, round-trip parsing diagnostic reproduces all 20 tables exactly. Individual parameter values are not thereby runtime-instrumented; v20 compatibility remains separate and unbound.",
                      "runtime_value": None})
    if path.endswith("/maldi.py"):
        found.append({"route": "MALDI_ingest_audit", "evidence": evidence(RUN / "maldi/summary.json"),
                      "execution_disposition": "evaluation audit; no per-call parameter trace or Stage07/08 receipt", "runtime_value": None})
    if "/integration/" in path or path.endswith(("/same_slide_registration.py", "/same_slide_integration.py",
                                                "/maldi_serial_registration.py", "/maldi_serial_integration.py")):
        found.append({"route": "spatial_registration_or_matching", "execution_disposition": "no accepted cohort registration/matching execution; synthetic probes separate",
                      "runtime_value": None})
    return found


def registry_traces(catalog, sources):
    result = []
    for table in catalog["registries"]:
        path = Path(table["source"])
        actual = list(csv.DictReader(path.open()))
        if actual != table["rows"] or sha(path.read_bytes()) != table["sha256"]:
            raise ValueError(f"registry changed: {path}")
        pool = [s for s in sources if s.ns == namespace(path) and
                ("/positivity/" in s.path or "/tissues/" in s.path)]
        fields = sorted({k for r in actual for k in r})
        field_reads = {field: [] for field in fields}
        const_reads = defaultdict(list)
        loaders = []
        for src in pool:
            for node in ast.walk(src.tree):
                if isinstance(node, ast.Call):
                    callee = ast.unparse(node.func)
                    if callee.endswith(".read_csv"):
                        loaders.append(src.ref(node))
                    if isinstance(node.func, ast.Attribute) and node.func.attr == "const" and node.args and isinstance(node.args[0], ast.Constant):
                        const_reads[str(node.args[0].value)].append(src.ref(node))
                    if isinstance(node.func, ast.Attribute) and node.func.attr == "get" and node.args and isinstance(node.args[0], ast.Constant):
                        key = node.args[0].value
                        if key in field_reads:
                            field_reads[key].append({"kind": "literal_get_key; receiver_type_not_proven", **src.ref(node)})
                elif isinstance(node, ast.Subscript):
                    for key in ast.walk(node.slice):
                        if isinstance(key, ast.Constant) and isinstance(key.value, str) and key.value in field_reads:
                            field_reads[key.value].append({"kind": "literal_subscript_key; receiver_type_not_proven", **src.ref(node)})
        rows = []
        for index, row in enumerate(actual, start=2):
            rowname = row.get("name") or row.get("marker") or row.get("id") or row.get("label")
            constant = "__constants.csv" in str(path)
            rows.append({"csv_line": index, "row_identifier": rowname, "values": row,
                         "class": "panel_configuration", "field_consumer_keys": sorted(row),
                         "constant_call_sites": const_reads.get(rowname, []) if constant else [],
                         "runtime_disposition": "row_value_execution_unrecorded",
                         "review_state": row.get("status") if "__proposals.csv" in str(path) else None,
                         "review_interpretation": "recorded source state only; no approval added by evaluation",
                         "output_disposition": "exact row-to-output lineage uninstrumented; table-level route evidence only"})
        result.append({"source": str(path), "sha256": table["sha256"], "producer": table["producer"],
                       "source_verified": True, "rows": rows, "field_read_candidates": field_reads,
                       "csv_read_calls_in_version_namespace": loaders,
                       "runtime_disposition": "PC67 tables may be bound by spleen receipt; no row-level execution trace. CODEX tables not established by PC67 run.",
                       "limitations": "Field-name matches are lookup candidates with receiver expressions, not proof that every registry row reaches that consumer."})
    return result


def literal_candidates(sources, declarations):
    result = []
    spans = defaultdict(list)
    for source, value in declarations:
        spans[source.path].append((value.lineno, value.col_offset,
                                   value.end_lineno, value.end_col_offset))
    for src in sources:
        if src.path not in spans:
            continue
        candidates = []
        for node in ast.walk(src.tree):
            if not isinstance(node, ast.Constant) or not isinstance(node.value, (str, int, float)) or isinstance(node.value, bool):
                continue
            parent = src.parent.get(node)
            if isinstance(parent, ast.Expr) and isinstance(src.parent.get(parent), (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                body = src.parent[parent].body
                if body and body[0] is parent and isinstance(node.value, str):
                    continue
            if any((a, b) <= (node.lineno, node.col_offset) and (node.end_lineno, node.end_col_offset) <= (c, d)
                   for a, b, c, d in spans[src.path]):
                continue
            context = next((p for p in src.ancestors(node) if isinstance(p, (ast.Call, ast.Compare, ast.BinOp, ast.Assign, ast.Return))), parent)
            disposition = literal_disposition(src, node)
            candidates.append({**src.ref(node), "literal_type": type(node.value).__name__,
                               "context": src.ref(context) if context is not None and hasattr(context, "lineno") else None,
                               "evidence_level": "candidate_not_dynamic_proof", **disposition})
        result.append({"source": src.path, "sha256": src.hash, "count": len(candidates),
                       "literals": candidates})
    return result


def literal_disposition(src, node):
    """Explicitly heuristic triage of remaining literals, never scientific approval."""
    ancestors = list(src.ancestors(node))
    parent = src.parent.get(node)
    call = next((a for a in ancestors if isinstance(a, ast.Call)), None)
    keyword = next((a.arg for a in ancestors if isinstance(a, ast.keyword)), None)
    field = None
    for ancestor in ancestors:
        if isinstance(ancestor, ast.Dict):
            for key, value in zip(ancestor.keys, ancestor.values):
                if isinstance(key, ast.Constant) and isinstance(key.value, str) and node in ast.walk(value):
                    field = key.value
                    break
            if field:
                break
    if "/integration/" in src.path or "serial_registration" in src.path or "serial_integration" in src.path:
        domain = "registration_geometry"
    elif src.path.endswith(("/cohort_scale.py", "/protein_evidence.py", "/maldi.py")):
        domain = "assay_specific_measurement"
    elif "/tissues/" in src.path or "/scripts/typing/" in src.path:
        domain = "tissue_biology"
    else:
        domain = "unresolved_assumption"
    kind, category, basis = "unresolved_literal_context", "unresolved_assumption", "role cannot be established from immediate AST context"
    method_keys = {"n_neighbors", "n_components", "n_top_genes", "resolution", "random_state", "seed",
                   "n_jobs", "n_iterations", "max_iter", "bins", "q", "quantile", "percentile", "min_dist",
                   "max_dist_um", "target_sum", "method", "metric", "flavor", "n_per_donor", "k", "tol"}
    criteria_words = {"threshold", "cutoff", "min_", "max_", "floor", "resolution", "fdr", "quantile", "percentile", "pixel", "pitch", "tolerance", "seed"}
    if isinstance(node.value, str) and any(isinstance(a, ast.Raise) for a in ancestors):
        kind, category, basis = "nonparameter_diagnostic_text_candidate", "shared_computation", "literal occurs in raised diagnostic text"
    elif keyword in method_keys:
        kind, category, basis = "scientific_method_argument_candidate", domain, f"explicit method keyword {keyword}"
    elif field and any(word in field.lower() for word in criteria_words) and isinstance(node.value, (int, float)):
        kind, category, basis = "scientific_criterion_or_reported_criterion_candidate", domain, f"numeric value under criterion-like key {field}; metadata equality is not consumer control"
    elif isinstance(node.value, (int, float)) and call and ast.unparse(call.func).split(".")[-1] in {
            "percentile", "quantile", "histogram", "log1p", "arcsinh", "clip", "threshold", "erode", "buffer"}:
        kind, category, basis = "scientific_numeric_transform_candidate", domain, "numeric argument to explicit fitting/transform operation"
    elif isinstance(node.value, (int, float)) and any(isinstance(a, ast.Compare) for a in ancestors):
        kind, category, basis = "unresolved_threshold_or_validation_bound", domain, "comparison bound: scientific criterion versus structural validity requires consumer context"
    elif isinstance(node.value, str) and call and ast.unparse(call.func).split(".")[-1] in {"astype", "to_numpy", "asarray", "array"}:
        kind, category, basis = "measurement_or_numeric_representation_candidate", "shared_computation", "numeric representation/dtype literal"
    elif isinstance(node.value, str) and ("/" in node.value or node.value.endswith((".csv", ".parquet", ".json", ".yml", ".h5ad"))):
        kind, category, basis = "unresolved_input_or_output_path_configuration", "unresolved_assumption", "path-like literal; input selection versus output naming must retain routing context"
    elif isinstance(node.value, str) and node.value.isupper() and node.value.isalnum() and 2 <= len(node.value) <= 12 and "/tissues/" in src.path:
        kind, category, basis = "tissue_gene_or_marker_literal_candidate", "tissue_biology", "uppercase tissue-worker token may select biological evidence; requires panel coverage"
    elif isinstance(node.value, str) and isinstance(parent, ast.Dict) and node in parent.keys:
        kind, category, basis = "nonparameter_schema_key_candidate", "shared_computation", "dictionary key; possible biological label semantics remain visible in context"
    elif isinstance(node.value, (int, float)) and isinstance(parent, ast.Subscript) and node is parent.slice:
        kind, category, basis = "nonparameter_index_candidate", "shared_computation", "literal positional index"
    return {"disposition": kind, "classification": category,
            "triage_basis": basis, "classification_method": "heuristic_AST_context; not scientific applicability or runtime proof"}


def verify_static_resolver(out):
    """Focused checks for shadowing, imports, version isolation, and overridden defaults."""
    with tempfile.TemporaryDirectory(dir=out) as temporary:
        directory = Path(temporary) / PKG
        directory.mkdir(parents=True)
        a = directory / "trace_fixture_a.py"
        b = directory / "trace_fixture_b.py"
        a.write_text("X = 7\ndef foo(x=3):\n    y = x * X\n    return y\ndef shadow(X):\n    return X\ndef capture(X=X):\n    return X\ndef caller():\n    return foo(x=9)\ndef defaultcaller():\n    return foo()\ndef starcaller(args):\n    return foo(**args)\n")
        b.write_text("from .trace_fixture_a import X as Y, foo as f\nZ = Y\ndef caller():\n    return f(11)\n")
        sa, sb = Source(a, "fixture"), Source(b, "fixture")
        reads = direct_consumers(sa, "multiplex5000.trace_fixture_a.X", [sa, sb])
        assert {(Path(r["source"]).name, r["line"]) for r in reads} == {
            (a.name, 3), (a.name, 7), (b.name, 2)}
        fn = next(n for n in sa.nodes if isinstance(n, ast.FunctionDef) and n.name == "foo")
        calls = callsites(sa, fn, "x", [sa, sb])
        assert Counter(c["disposition"] for c in calls) == {
            "explicit_keyword_override": 1, "explicit_positional_override": 1,
            "default_selected_if_call_executes": 1, "unresolved_star_arguments": 1}
        assert len(direct_consumers(sa, "multiplex5000.trace_fixture_a.X", [sa, Source(b, "different_version")])) == 2
        own_reads = direct_consumers(sa, "multiplex5000.trace_fixture_a.foo.x", [sa, sb])
        flows = output_candidates(sa, "multiplex5000.trace_fixture_a.foo.x", [sa, sb], own_reads)
        assert any(b["kind"] == "return_to_caller" for flow in flows for b in flow["boundaries"])
    return {"status": "pass", "checks": ["shadowed global excluded", "definition-time default capture", "explicit import aliases",
            "keyword/positional/default/star call distinctions", "producer namespace isolation", "assignment-to-return candidate"]}


def receipt_settings(snapshots):
    result = []
    for route, receiptpath in [("rna", RUN / "remote/rna/receipt.json"),
                               ("pln", RUN / "pln/receipt.json"),
                               ("spleen", RUN / "remote/spleen/receipt.json")]:
        receipt = json.loads(receiptpath.read_text())
        commit = receipt["code"]["commit"]
        pool = [s for s in snapshots if s.ns == route + "@" + commit]
        for key, value in receipt["settings"].items():
            declarations, reads = [], []
            for src in pool:
                for node in src.nodes:
                    if isinstance(node, ast.Dict) and src.scope_name(src.scope(node)).endswith(".settings"):
                        for k, val in zip(node.keys, node.values):
                            if isinstance(k, ast.Constant) and k.value == key:
                                declarations.append({"key": src.ref(k), "value_expression": src.ref(val)})
                    if isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Constant) and node.slice.value == key:
                        if ast.unparse(node.value).endswith(".settings"):
                            reads.append(src.ref(node))
                    elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                          and node.func.attr == "get" and ast.unparse(node.func.value).endswith(".settings")
                          and node.args and isinstance(node.args[0], ast.Constant) and node.args[0].value == key):
                        reads.append(src.ref(node))
            result.append({"route": route, "producer_commit": commit, "setting": key,
                           "resolved_value": value, "role": "receipt_declared_setting",
                           "evidence": evidence(receiptpath, "/settings/" + key),
                           "source_declarations": declarations, "explicit_settings_reads": reads,
                           "output_receipts": evidence(receiptpath, "/outputs"),
                           "execution_disposition": "receipt_value_verified; individual_consumer_branch_execution_unrecorded",
                           "consumer_disposition": "explicit_static_settings_read" if reads else "no_explicit_settings_read_in_archived_route_sources; may_be_metadata_or_separately_hardcoded",
                           "scientific_acceptance": "not_established_by_receipt"})
    return result


def contextual_findings(sources, snapshots):
    def function_refs(suffix, names, ns=None):
        return [s.ref(n) for s in [*sources, *snapshots] if s.path.endswith(suffix) and
                (ns is None or s.ns == ns) for n in s.nodes
                if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name in names]

    port = "port@df3fe1d"
    hook_functions = {"multiplex5000.positivity.tissue.hooks", "multiplex5000.positivity.tissue.check_registry"}
    explicit_calls = [s.ref(n) for s in sources if s.ns == port and s.path.endswith(("/levels.py", "/identity.py"))
                      for n in s.nodes if isinstance(n, ast.Call) and s.resolve_expr(n.func) in hook_functions]
    normalization = RUN / "remote/rna/pln_embedding_normalization.json"
    normalization_data = json.loads(normalization.read_text())
    candidates = [r["normalization_target"] for r in normalization_data["acquisitions"]
                  if r["normalization_target"] is not None]
    return [
        {"id": "ignored_tissue_hooks", "class": "unresolved_assumption",
         "source": function_refs("/positivity/tissue.py", {"hooks", "check_registry"}, port),
         "consumer_source": function_refs("/positivity/levels.py", {"choose", "image_pass"}, port) +
                            function_refs("/positivity/identity.py", {"classify_image", "default_start_class"}, port),
         "exact_explicit_hook_calls_in_levels_and_identity": explicit_calls,
         "finding": "The standalone hook API is present. No lexically resolved hooks/check_registry call exists in port levels.py or identity.py. classify_image permits caller-supplied start_class and otherwise takes the pancreas default. Absence of an explicit call is not proof about dynamic external wrappers.",
         "disposition": "hook_enforcement_not_established", "remediation": "Specify and verify tissue dispatch at the route boundary, reject unsupported reference/context requirements, and record hook version and selection before fitting."},
        {"id": "protein_availability_and_refinement", "class": "panel_configuration",
         "source": function_refs("/scripts/typing/refine.py", {"donor_values", "refine_donor", "neighbours"}),
         "evidence": [evidence(RUN / "protein_primary_expected_universe.json"),
                      evidence(RUN / "protein_primary/per_acquisition.json"), evidence(RUN / "protein_primary/pancreas/run.json")],
         "finding": "The schema-presence have list includes an all-NaN unacquired marker. The all-finite refinement admission mask therefore empties the graph for the four sections lacking Ker8_18. All 25 primary partitions reconcile to 31,385,695 IDs; refinement fails at 0112, so complete final refinement and downstream output lineage are unavailable.",
         "unavailable_marker_sections": ["0112", "6347", "6533", "6580"],
         "disposition": "failed_existing_route_no_scientific_fix_applied",
         "remediation": "Define assay/panel marker availability independently from schema presence; specify allowed representation marker set and explicit insufficient-signal outcome. Reconcile every primary ID and preserve unavailable evidence before rerunning all 25 sections."},
        {"id": "protein_normalization_contexts", "class": "assay_specific_measurement",
         "source": function_refs("/cohort_scale.py", {"build_scale", "background_mode_floor", "marker_scale", "apply_scale"}) +
                   function_refs("/scripts/typing/refine.py", {"representation_params", "donor_values", "to_Z"}) +
                   function_refs("/stages/io60_mixed_protein_pilot.py", {"settings", "run"}),
         "evidence": [evidence(RUN / "protein_primary/pancreas/run.json", "/representation_donors"),
                      evidence(RUN / "protein_primary/pancreas/refinement/representation_provenance.json")],
         "finding": "The primary cohort scale filters QC plus in_alignment_block, sorts native IDs, then samples up to 200,000/donor. Refinement representation joins geometry QC and samples without the alignment-block predicate; it fits cofactor/centre/spread across 25 donors although wrapper metadata says single_donor. Missing mode floor falls back to 0 in marker_scale; invalid scale/representation spreads use 1. These are separate operations and contexts, not transferable tissue defaults.",
         "disposition": "fit_artifact_available; scope_label_inconsistent; final_refinement_unavailable",
         "remediation": "Record separately the exact fit population, available markers, fitted values, selected fixed settings and units for each normalization operation; correct metadata only in a separately authorized scientific change."},
        {"id": "rna_normalization_contexts", "class": "assay_specific_measurement",
         "source": function_refs("/tissues/lymph_node_xenium_worker.py", {"partition_section", "process_exact_section"}) +
                   function_refs("/tissues/pancreas_xenium_worker.py", {"cluster_section", "process_section"}),
         "evidence": evidence(normalization),
         "finding": "RNA detection uses a fixed 10,000-count normalization target. pLN embedding calls normalize_total without target_sum, which uses the post-QC median in the pinned Scanpy implementation. The embedding targets were reconstructed from admitted cells; they are not original receipt settings.",
         "fixed_detection_target": 10000.0,
         "reconstructed_embedding_target_range": [min(candidates), max(candidates)],
         "reconstructed_values": normalization_data["acquisitions"],
         "disposition": "reconstructed_fitted_values_separate_from_fixed_detection_settings",
         "remediation": "Emit normalization target and admitted-cell count from each operation and keep fitted values distinct from fixed detection thresholds."},
        {"id": "pln_cutoff_metadata_vs_consumer", "class": "unresolved_assumption",
         "source": function_refs("/stages/pln_typing_pilot.py", {"settings", "run"}) +
                   function_refs("/protein_evidence.py", {"candidate_calls"}),
         "evidence": evidence(RUN / "pln/receipt.json", "/settings/protein_cutoff_fraction"),
         "finding": "The receipt records protein_cutoff_fraction=0.1. The wrapper does not pass that settings field to evaluate_protein_section; protein_evidence uses its REQUIRE_MIN constant. Equal declarations do not establish that configuration controls execution.",
         "disposition": "recorded_setting_not_governing_consumer", "remediation": "Specify a single resolved source of the cutoff, enforce agreement if metadata remains separate, and record the actual selected value at the consumer."},
    ]


def fitted_artifact_records():
    import pyarrow.parquet as pq

    records = []
    paths = [(RUN / "protein_primary/pancreas/refinement/representation_params.parquet", "protein_representation", "fitted_before_route_failure"),
             (RUN / "pln_recheck/protein_thresholds.csv", "pln_thresholds", "recomputed_evaluation_values")]
    for path, operation, status in paths:
        if not path.exists():
            records.append({"operation": operation, "path": str(path), "status": "artifact_unavailable"})
            continue
        values = pq.read_table(path).to_pandas().reset_index().to_json(orient="records") if path.suffix == ".parquet" else None
        rows = json.loads(values) if values is not None else list(csv.DictReader(path.open()))
        records.append({"operation": operation, "role": "fitted_values", "status": status,
                        "artifact": evidence(path), "rows": rows,
                        "consumer_execution_limit": "Artifact values exist; they do not prove downstream consumer completion or biological acceptance."})
    for name in ["remote/spleen/tables", "remote/rna/reports"]:
        directory = RUN / name
        records.append({"operation": name, "role": "fitted_or_diagnostic_artifact_index",
                        "status": "producer_artifacts; see receipt output hashes",
                        "artifacts": [evidence(p) for p in sorted(directory.rglob("*")) if p.is_file()]})
    return records


def original_rna_reports(out, fetch=False):
    receipt = json.loads((RUN / "remote/rna/receipt.json").read_text())
    declared = [r for r in receipt["outputs"] if r["path"].endswith(".report.json")]
    target = out / "producer_rna_reports"
    if fetch:
        base = "/blue/maigan/smith6jt/mplx5k-io60-rna-20261027/releases/io60-candidate-1/06_typing"
        command = shlex.join(["tar", "-C", base, "-cf", "-", "--", *[r["path"] for r in declared]])
        archive = out / "producer_rna_reports.tar"
        archive.parent.mkdir(parents=True, exist_ok=True)
        with archive.open("wb") as handle:
            subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "hpg", command],
                           stdout=handle, check=True)
        with tarfile.open(archive) as handle:
            handle.extractall(target, filter="data")
    records = []
    for declared_output in declared:
        path = target / declared_output["path"]
        if not path.is_file():
            raise ValueError(f"missing exact producer report: {path}; rerun with --fetch-rna-reports for read-only HPG recovery")
        digest = sha(path.read_bytes())
        if digest != declared_output["sha256"]:
            raise ValueError(f"original RNA report hash mismatch: {path}")
        local = RUN / "remote/rna/reports" / (path.parent.name + ".json")
        equal = json.loads(local.read_text()) == json.loads(path.read_text())
        if not equal:
            raise ValueError(f"evaluation report semantics differ from original: {local}")
        records.append({"producer_report": evidence(path), "receipt_output": declared_output,
                        "evaluation_report": evidence(local), "json_semantically_equal": equal,
                        "byte_identity_required_only_for_original": True})
    dump(out / "producer_rna_reports_manifest.json", {"receipt": evidence(RUN / "remote/rna/receipt.json"),
         "reports": records, "interpretation": "Evaluation copies are reserialized JSON. Original reports match the scientific receipt exactly; all evaluation copies are semantically identical."})
    return target


def bind_reported_values(rows, snapshots, report_root):
    """Bind a small explicit set of direct source-expression -> report-field edges.

    No numeric matching: each field must directly read the parameter symbol at
    the pinned producer, and each complete report must match the receipt hash.
    """
    receiptpath = RUN / "remote/rna/receipt.json"
    receipt = json.loads(receiptpath.read_text())
    output_hashes = {r["sha256"]: r["path"] for r in receipt["outputs"]}
    bindings = {
        ("pancreas_xenium_worker.py", "process_section.score_on"): ("score_on", "score_on"),
        ("pancreas_xenium_worker.py", "process_section.naming_policy"): ("naming_policy", "naming_policy"),
        ("pancreas_xenium_worker.py", "process_section.source_serial"): ("source_serial", "source_serial"),
        ("pancreas_xenium_worker.py", "Y_FLOOR"): ("y_floor_group", "y_floor_group"),
        ("pancreas_xenium_worker.py", "KNN_K"): ("knn_assignment_k", "cluster_meta/knn_assignment_k"),
        ("pancreas_xenium_worker.py", "cluster_section.resolution"): ("leiden_resolution", "cluster_meta/leiden_resolution"),
        ("lymph_node_xenium_worker.py", "Y_FLOOR"): ("y_floor_group", "y_floor_group"),
        ("lymph_node_xenium_worker.py", "DETECTION_PREDICATE"): ("detection_predicate", "detection_predicate"),
    }
    for row in rows:
        declaration = row["declaration"]
        key = Path(declaration["source"]).name, declaration["name"]
        if key not in bindings:
            continue
        field, pointer = bindings[key]
        producer = next(s for s in snapshots if s.ns.startswith("rna@") and
                        Path(s.path).name == key[0] and s.hash == declaration["sha256"])
        edges = []
        for node in producer.nodes:
            if isinstance(node, ast.Dict):
                for k, value in zip(node.keys, node.values):
                    if isinstance(k, ast.Constant) and k.value == field and producer.resolve_expr(value) == row["resolved_symbol"]:
                        edges.append({"key": producer.ref(k), "value": producer.ref(value)})
        if not edges:
            raise ValueError(f"no exact producer source-to-report edge: {key}")
        observed = []
        pattern = "*__panc.report.json" if key[0].startswith("pancreas_") else "*__pln_1.report.json"
        for path in sorted(report_root.rglob(pattern)):
            data = json.loads(path.read_text())
            value = data
            for component in pointer.split("/"):
                if not isinstance(value, dict) or component not in value:
                    raise ValueError(f"missing expected report field {path}: {pointer}")
                value = value[component]
            digest = sha(path.read_bytes())
            if digest not in output_hashes:
                raise ValueError(f"report does not match receipt output hash: {path}")
            observed.append({"acquisition_id": data["acquisition_id"], "section_key": data["section_key"],
                             "value": value, "report": evidence(path, "/" + pointer),
                             "receipt_output_path": output_hashes[digest]})
            if key[1] in {"KNN_K", "cluster_section.resolution"}:
                for index, sub in enumerate(data.get("refine_meta", [])):
                    subfield = "meta_" + field
                    if subfield in sub:
                        observed.append({"acquisition_id": data["acquisition_id"], "parent_cluster": sub["parent"],
                                         "value": sub[subfield],
                                         "report": evidence(path, f"/refine_meta/{index}/{subfield}"),
                                         "receipt_output_path": output_hashes[digest]})
        row["resolved_runtime_value"] = {
            "status": "producer_report_field_verified_against_receipt",
            "source_to_field_edges": edges, "observations": observed,
            "receipt": evidence(receiptpath, "/outputs"),
            "scope": "These reported invocations only. A function default may have been overridden; observed values are not blanket branch coverage or scientific acceptance."}
        row["disposition"] = "reported_parameter_values_verified_for_recorded_invocations"
        row["runtime_limit"] = "Direct source-to-report value binding verified for reported sections/subclusters; other invocations and unreported branches remain unresolved."


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=DOC / "parameter_catalog.json")
    parser.add_argument("--out", type=Path, default=Path("outputs/handoff_completion_20261102"))
    parser.add_argument("--fetch-rna-reports", action="store_true",
                        help="Recover only receipt-declared report JSONs read-only over the existing hpg SSH connection.")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    resolver_checks = verify_static_resolver(args.out)
    catalog = json.loads(args.catalog.read_text())
    source_paths = {normalized(r["source"]) for r in catalog["parameters"]}
    # Resolve import references only inside explicit, separate producer namespaces.
    source_paths.update(str(p) for p in PKG.rglob("*.py"))
    source_paths.update(str(p) for p in Path("core/multiplex5000/scripts/typing").glob("*.py"))
    sources = [Source(path) for path in sorted(source_paths)]
    lookup = {s.path: s for s in sources}
    producers, snapshots = archive_producers(args.out)
    primary = json.loads((RUN / "protein_primary/pancreas/run.json").read_text())
    rows, declarations = [], []
    for index, row in enumerate(catalog["parameters"]):
        src = lookup[normalized(row["source"])]
        if src.hash != row["sha256"]:
            raise ValueError(f"catalog source changed: {src.path}")
        declaration, value, symbol = src.declaration(row)
        if ast.unparse(value) != row["expression"]:
            raise ValueError(f"catalog expression changed: {row}")
        declarations.append((src, value))
        consumers = direct_consumers(src, symbol, sources)
        defaults = callsites(src, declaration, row["name"].rsplit(".", 1)[1], sources) if row["kind"] == "function_default" else []
        boundaries = output_candidates(src, symbol, sources, consumers)
        runtime = route_evidence(row, src, producers, primary)
        producer_defaults = []
        if row["kind"] == "function_default":
            for actual in snapshots:
                if actual.module == src.module and actual.hash == src.hash:
                    actual_fn, _, _ = actual.declaration(row)
                    producer_defaults.extend(callsites(actual, actual_fn, row["name"].rsplit(".", 1)[1], snapshots))
        rows.append({"catalog_index": index, "parameter_id": sha(f"{row['source']}:{row['line']}:{row['name']}:{row['kind']}".encode())[:20],
                     "declaration": row, "source_verified": True, "namespace": src.ns, "resolved_symbol": symbol,
                     "classification": classification(row), "consumers": consumers, "default_callsites": defaults,
                     "output_trace": boundaries, "route_evidence": runtime,
                     "actual_producer_default_callsites": producer_defaults,
                     "resolved_runtime_value": None,
                     "disposition": "static_trace_runtime_unresolved" if consumers else "no_resolved_read_in_inspected_sources",
                     "runtime_limit": "No per-parameter executed-value record; declaration/route receipt alone cannot prove branch execution or selection of a default.",
                     "remediation": "Record resolved settings at consumer entry, including overrides, units, source/registry hashes and output receipt. Preserve unavailable fitted values; do not infer execution from static reads."})
    registries = registry_traces(catalog, sources)
    report_root = original_rna_reports(args.out, args.fetch_rna_reports)
    bind_reported_values(rows, snapshots, report_root)
    literals = literal_candidates(sources, declarations)
    dump(args.out / "parameter_literal_candidates.json", {"schema": 1, "scope": "all non-docstring numeric/string literals outside catalog declaration expressions in catalog source files",
         "not_exhaustive": "Not all literals are scientific settings. Dynamic config, imported library defaults, environment and data-fitted values require separate traces.",
         "sources": literals})
    source_manifest = [{"path": s.path, "namespace": s.ns, "module": s.module, "sha256": s.hash} for s in sources]
    dump(args.out / "parameter_source_manifest.json", {"sources": source_manifest, "producer_snapshots": producers})
    dependencies = list(csv.DictReader((DOC / "dependencies.csv").open()))
    resolved_settings = receipt_settings(snapshots)
    findings = contextual_findings(sources, snapshots)
    fits = fitted_artifact_records()
    result = {"schema": 1, "catalog": evidence(args.catalog), "generator": evidence(Path(__file__)),
              "scope": "Every catalog declaration and registry row; exact static consumer references, conservative output candidates, explicit unproven runtime values.",
              "method_limits": ["AST is static evidence, not execution coverage.",
                                "Lexical/import resolution does not prove runtime monkeypatching, dynamic import, method receiver type or branches.",
                                "Output assignment closure is flow-insensitive and may overapproximate dependencies; no artifact-level causality is claimed.",
                                "Default callsite analysis covers explicit lexical/import bindings only; unresolvable methods, star arguments and dynamic calls remain unresolved.",
                                "Complete route receipts establish route outputs; individual source declarations can remain unused or overridden."],
              "parameters": rows, "registries": registries, "dependencies": dependencies,
              "resolver_verification": resolver_checks,
              "resolved_receipt_settings": resolved_settings,
              "fitted_values_and_diagnostics": fits,
              "contextual_findings": findings,
              "producer_snapshots": producers,
              "literal_inventory": evidence(args.out / "parameter_literal_candidates.json"),
              "source_manifest": evidence(args.out / "parameter_source_manifest.json"),
              "original_rna_reports_manifest": evidence(args.out / "producer_rna_reports_manifest.json"),
              "summary": {"catalog_declarations": len(catalog["parameters"]), "traced_declarations": len(rows),
                          "unique_parameter_ids": len({r["parameter_id"] for r in rows}),
                          "registry_tables": len(registries), "registry_rows": sum(len(r["rows"]) for r in registries),
                          "catalog_source_hashes_verified": len({r["declaration"]["source"] for r in rows}),
                          "declaration_dispositions": dict(Counter(r["disposition"] for r in rows)),
                          "classifications": dict(Counter(r["classification"]["class"] for r in rows)),
                          "declarations_with_receipt_bound_reported_values": sum(r["resolved_runtime_value"] is not None for r in rows),
                          "receipt_declared_settings": len(resolved_settings),
                          "fitted_artifact_groups": len(fits),
                          "literal_dispositions": dict(Counter(v["disposition"] for r in literals for v in r["literals"])),
                          "remaining_literal_candidates": sum(r["count"] for r in literals)},
              "remaining_work": [
                  {"id": "CORE-1-runtime", "requirement": "Add opt-in consumer-entry parameter receipts with explicit fixed/selected/fitted/review roles, qualified producer identities, unavailable values, input/output receipts and branch scope. Evaluation does not authorize this scientific instrumentation change."},
                  {"id": "CORE-1-literals", "requirement": "The literal inventory explicitly disposes every candidate using labelled AST heuristics. Review the unresolved contexts and confirm scientific candidates when specifying production resolved-parameter emission; preserve nonparameter roles. Do not call the 572-entry declaration catalog exhaustive."},
                  {"id": "CORE-1-library-defaults", "requirement": "Pin and report external library defaults (including pLN normalize_total target_sum=None), CLI/env-derived values, fit populations and actual overrides; retain reconstructed values separately."},
                  {"id": "CORE-2-hooks", "requirement": "Specify which tissue hooks must be invoked by each route, test dispatch and absent context separately; current hook declarations are not evidence of enforcement."}],
              "complete": False,
              "completion_interpretation": "Static accounting is complete for the catalog. Runtime source-to-output trace and scientific transfer acceptance remain incomplete."}
    assert len(rows) == len(catalog["parameters"]) == len({r["parameter_id"] for r in rows})
    assert len(registries) == len(catalog["registries"])
    dump(args.out / "parameter_traces.json", result)
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
