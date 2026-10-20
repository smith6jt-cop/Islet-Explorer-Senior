"""Superseded work must be UNREACHABLE, not merely documented as retired.

WHY THIS FILE EXISTS.  The refinement pass -- ``refine.py``, ``REFINEMENT.md`` and every test
guarding the typing -- sat UNTRACKED in a working tree for days while the maintainer had said it
was of paramount importance that it be set in stone so that no older, incorrect work could be
used in its place.  Nothing failed.  A working tree is not a record: an editor that reverts, a
stash, a checkout, a second clone, and the method is gone while the superseded generation beside
it is still perfectly readable.  That is the whole episode, and it has two halves:

1.  the current method was not in git, and
2.  the superseded generation was still readable by anything that asked for it.

This file is the guard for both.  It asserts four properties, all DERIVED FROM THE REPOSITORY --
never from a list of results kept here, never from a cell count, and never by reproducing a
biological finding, which would measure the expectation rather than the build:

``test_every_typing_module_and_guard_is_tracked_in_git``
    Every ``.py`` and ``.md`` under ``scripts/typing/`` and every test under ``tests/current/``
    is known to git.  ``git ls-files --error-unmatch``, one file at a time, so the failure names
    the file rather than a count.

``test_no_typing_file_lives_only_as_a_working_tree_edit``
    The weaker half of the same failure: a file that IS tracked but whose current content exists
    only as an uncommitted edit is still one ``git checkout`` from gone.

``test_no_live_reader_of_a_superseded_artifact``
    No producer under ``scripts/senior/`` or ``core/multiplex5000/scripts/`` READS a superseded
    artifact.  A refusal is fine -- that is the fix.  Prose is fine -- history belongs in
    comments, and a comment is not in the AST at all, so it cannot be confused for a read.
    WHAT IS SUPERSEDED, AT WHICH SEVERITY, AND WHO MAY READ IT ANYWAY is not decided here: it is
    read from the one table, ``core/multiplex5000/scripts/typing/superseded.py``.  A read the
    registry names in a record's ``allowed_readers`` -- ``<file>::<function>`` -- is cleared by
    that entry and by nothing typed into this file.  A read inside a module the repository
    itself declares retired (``scripts/retire/dispositions.py``, code tier 1) is not a LIVE
    reader; it is reported by ``test_no_retired_producer_still_reads_a_superseded_artifact``,
    which is a strict xfail until the retire lane removes those modules.

``test_every_artifact_read_as_an_input_has_a_producer``
    An orphan is never a live input.  A file with no producer in the repository cannot be
    re-made, so whatever is on disk is frozen at whatever generation happened to write it, and
    nothing can tell you which.  ``scripts/retire/dispositions.py`` names this failure at its
    code-tier banner: deleting a producer "makes a live artifact producerless, which is the
    exact failure ``typing/umap/refinement_log.json`` is already in".  The property asserted
    here is the general one, not that one case.  An artifact the registry names as superseded is
    that guard's, not this one's; an artifact the registry declares in ``PRODUCERLESS`` -- a
    ruling waiting to be made -- is reported by its own strict xfail rather than here, so the
    ruling stays in every test report until it is made.

``test_reconstruction_names_every_concept_the_refinement_rests_on``
``test_reconstruction_names_every_artifact_the_refinement_touches``
    ``RECONSTRUCTION.md`` is the document the maintainer named -- "read it rather than
    re-deriving it".  A reader who does that and finds a concept missing re-derives it, and a
    re-derivation is how a second generation gets made.  The concepts are DERIVED from the
    intersection of ``refine.py``'s module-level definitions and the identifiers
    ``REFINEMENT.md`` marks in code spans, so the test keeps working as the method grows.

TWO THINGS A READER OF THE FIRST FAILURE SHOULD EXPECT.  This file names ITSELF until it is
committed -- that is the property working, not a bug, and the same is true of any new module a
sibling lane has just written beside the typing code.  And a file that is tracked but modified
is named by the second test rather than the first, because the two failures want different
fixes: ``git add`` for one, a commit for the other.

WHAT THIS FILE IS NOT.  It makes no method decision, wires nothing, reads no result, asserts no
count, and touches no file on any ``data/`` root.  It runs read-only git subcommands only --
``ls-files``, ``status``, ``rev-parse`` -- and :func:`_git` refuses any other.

Run from the repository root:

    .venv/bin/python -m pytest \
      core/multiplex5000/tests/current/test_no_superseded_work_is_reachable.py -v
"""

from __future__ import annotations

import ast
import fnmatch
import pathlib
import re
import subprocess
import sys

import pytest

# tests/current/ is two below core/multiplex5000/, which is two below the repo root.
REPO = pathlib.Path(__file__).resolve().parents[4]
TYPING = REPO / "core/multiplex5000/scripts/typing"
CURRENT = REPO / "core/multiplex5000/tests/current"

sys.path.insert(0, str(TYPING))
import superseded as SUP  # noqa: E402

#: The trees the lane calls producers.  ``tests/`` subtrees inside them are excluded: a test that
#: names a superseded path in order to assert it is refused is doing this file's job, not
#: defeating it.
PRODUCER_TREES = (REPO / "scripts/senior", REPO / "core/multiplex5000/scripts")

#: The superseded artifacts, as patterns matched against the text of a path expression -- DERIVED
#: from the one table, never kept here.  Each fragment and alias of every registry record maps to
#: its record, so the severity, the reason and the allowed readers all come from
#: ``core/multiplex5000/scripts/typing/superseded.py``.  This file carries no policy of its own:
#: an artifact is superseded because the registry records it, a reader is cleared because the
#: registry names it, and ``test_the_superseded_table_is_itself_declared_in_the_repository``
#: still cross-checks EVERY fragment and alias of every record against the repository's prose,
#: so the registry cannot invent a retirement the maintainer never declared either.
SUPERSEDED: dict[str, SUP.Record] = {
    fragment: record
    for record in SUP.REGISTRY
    for fragment in (*record.fragments, *record.aliases)
}

#: Calls that read.  ``exists`` / ``is_dir`` / ``is_file`` are probes rather than reads, but a
#: probe on a superseded root is how a silent fallback begins -- the caller finds nothing and
#: carries on with a default -- so they are reported, under their own label.
READ_FUNCS = {
    "read_parquet", "read_table", "read_csv", "read_text", "read_bytes", "read_json",
    "read_feather", "read_excel", "read_h5ad", "read_zarr", "ParquetFile", "ParquetDataset",
    "load", "loads", "safe_load", "glob", "iglob", "open", "listdir", "scandir",
}
READ_METHODS = {"glob", "rglob", "iterdir", "read_text", "read_bytes", "open", "joinpath"}
PROBE_METHODS = {"exists", "is_dir", "is_file", "stat", "lstat"}
WRITE_FUNCS = {
    "to_parquet", "to_csv", "to_json", "to_feather", "write_table", "write_dataset", "write_text",
    "write_bytes", "write_h5ad", "write_zarr", "dump", "savefig", "mkdir", "makedirs", "to_excel",
    "write_parquet", "write_csv", "write_feather", "savez", "save",
}
#: Names whose function is to refuse.  A literal inside one of these is the fix, not the defect.
REFUSAL_PREFIXES = ("refuse", "reject", "forbid", "_refuse", "_reject", "assert_not")
READ_ONLY_GIT = {"ls-files", "status", "rev-parse", "diff"}


# --------------------------------------------------------------------- plumbing ----
def _git(*args: str) -> subprocess.CompletedProcess:
    """Run a READ-ONLY git subcommand.  Anything that could change state is refused here."""
    assert args and args[0] in READ_ONLY_GIT, (
        f"this test may only run read-only git subcommands {sorted(READ_ONLY_GIT)}; "
        f"refusing {args[0]!r}")
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)


def _py_files(tree: pathlib.Path) -> list[pathlib.Path]:
    return sorted(p for p in tree.rglob("*.py")
                  if "__pycache__" not in p.parts and "tests" not in p.parts)


def _rel(path: pathlib.Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:  # pragma: no cover - every path here is inside the repo
        return str(path)


def _parents(tree: ast.AST) -> dict[int, ast.AST]:
    out: dict[int, ast.AST] = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            out[id(child)] = node
    return out


def _ancestors(node: ast.AST, parents: dict[int, ast.AST]) -> list[ast.AST]:
    chain: list[ast.AST] = []
    while id(node) in parents:
        node = parents[id(node)]
        chain.append(node)
    return chain


def _docstring_ids(tree: ast.AST) -> set[int]:
    """Every Constant that is a docstring.  Comments never reach the AST, so they need no rule."""
    out: set[int] = set()
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) \
                and isinstance(body, list) and body \
                and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
            out.add(id(body[0].value))
    return out


def _callee(call: ast.Call) -> str | None:
    func = call.func
    if isinstance(func, ast.Attribute):
        return func.attr
    if isinstance(func, ast.Name):
        return func.id
    return None


def _path_text(node: ast.AST, env: dict[str, str] | None = None) -> str:
    """The literal text of a path expression, with each f-string hole collapsed to ``*``.

    Only the path-building spine is followed -- ``/`` joins, f-strings, ``Path(...)`` and
    ``str(...)``.  Walking the whole subtree instead would splice unrelated literals from the same
    call into one nonsense key (``regions.json`` + a ``frame`` column name), and a key that is not
    a real path cannot be matched against a real producer.
    """
    if isinstance(node, ast.Constant):
        return node.value if isinstance(node.value, str) else ""
    if isinstance(node, ast.JoinedStr):
        out = []
        for part in node.values:
            if isinstance(part, ast.Constant) and isinstance(part.value, str):
                out.append(part.value)
            elif env is not None and isinstance(part, ast.FormattedValue) \
                    and isinstance(part.value, ast.Name) and part.value.id in env:
                out.append(env[part.value.id])
            else:
                out.append("*")
        return "".join(out)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left, right = _path_text(node.left, env), _path_text(node.right, env)
        if left and right:
            return f"{left.rstrip('/')}/{right.lstrip('/')}"
        return left or right
    if isinstance(node, ast.Call) and _callee(node) in {"Path", "str", "PurePath", "PosixPath"} \
            and node.args:
        return _path_text(node.args[0], env)
    if isinstance(node, ast.Attribute) and node.attr in {"parent", "resolve"}:
        return _path_text(node.value, env)
    if env is not None and isinstance(node, ast.Name) and node.id in env:
        return env[node.id]
    return ""


def _names_in(node: ast.AST) -> set[str]:
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def _enclosing_def(node: ast.AST, parents: dict[int, ast.AST]) -> str | None:
    for anc in _ancestors(node, parents):
        if isinstance(anc, (ast.FunctionDef, ast.AsyncFunctionDef)):
            return anc.name
    return None


# -------------------------------------------------- 1. nothing important is untracked ----
def _method_and_guard_files() -> list[pathlib.Path]:
    """The set that must be in git, derived from the tree rather than listed here."""
    files = [p for p in TYPING.rglob("*") if p.is_file()
             and p.suffix in {".py", ".md", ".R"} and "__pycache__" not in p.parts]
    files += [p for p in CURRENT.glob("*.py") if p.is_file()]
    return sorted(set(files))


def test_every_typing_module_and_guard_is_tracked_in_git():
    """The exact failure that lost the refinement for days: important work, not in git.

    ``refine.py``, ``REFINEMENT.md`` and every test guarding the typing sat untracked in a
    working tree.  Nothing said so, and the superseded generation beside them stayed readable.
    """
    checked = _method_and_guard_files()
    assert checked, f"derived no files to check under {_rel(TYPING)} and {_rel(CURRENT)}"
    untracked = []
    for path in checked:
        done = _git("ls-files", "--error-unmatch", "--", _rel(path))
        if done.returncode != 0:
            untracked.append(path)
    assert not untracked, (
        f"{len(untracked)} of {len(checked)} files that carry the typing method or guard it are "
        f"NOT IN GIT.  A working tree is not a record: a checkout, a stash or a second clone "
        f"loses them while the superseded generation beside them stays readable.  That is the "
        f"episode this file exists for.\n  "
        + "\n  ".join(f"{_rel(p)}  -- `git add` it and commit it with its tests"
                      for p in untracked))


def test_no_typing_file_lives_only_as_a_working_tree_edit():
    """Tracked is not enough: content that exists only as an uncommitted edit is still one
    ``git checkout`` from gone, which is how the refinement was lost the first time."""
    checked = _method_and_guard_files()
    done = _git("status", "--porcelain", "--", *[_rel(p) for p in checked])
    assert done.returncode == 0, f"git status failed: {done.stderr.strip()}"
    dirty = [line for line in done.stdout.splitlines() if line.strip()]
    assert not dirty, (
        f"{len(dirty)} file(s) carrying the typing method or its guards have content that exists "
        f"ONLY in the working tree.  Commit them -- an uncommitted edit is not set in stone.\n  "
        + "\n  ".join(dirty))


# ------------------------------------------ 2. no live reader of a superseded artifact ----
def _superseded_hits(path: pathlib.Path) -> list[dict]:
    """Classify every non-docstring string literal that names a superseded artifact.

    Returns one record per hit with ``kind`` in ``read`` / ``probe`` / ``declared_input`` /
    ``refusal`` / ``inert``.  A comment cannot appear here at all -- comments are not in the AST
    -- and a docstring is dropped explicitly, so prose about history can never be mistaken for a
    read.
    """
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError:
        return []
    parents = _parents(tree)
    docstrings = _docstring_ids(tree)

    # Names a superseded path literal flows into: name -> (line bound, which artifact).
    tainted: dict[str, tuple[int, str]] = {}
    # Names that exist only to be raised: the refusal message pattern.
    raised_names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Raise):
            raised_names |= _names_in(node)

    def _is_refusal(node: ast.AST) -> bool:
        chain = _ancestors(node, parents)
        if any(isinstance(a, ast.Raise) for a in chain):
            return True
        enclosing = _enclosing_def(node, parents)
        if enclosing and enclosing.lower().startswith(REFUSAL_PREFIXES):
            return True
        # The literal is handed to the guard inline -- ``read_parquet(refuse_if_superseded(p))``
        # or ``ROOT = refuse_if_superseded(...)`` -- which is the loud form of the read.
        for anc in chain:
            if isinstance(anc, ast.Call) \
                    and (_callee(anc) or "").lower().startswith(REFUSAL_PREFIXES):
                return True
        for anc in chain:
            if isinstance(anc, (ast.Assign, ast.AnnAssign)):
                targets = anc.targets if isinstance(anc, ast.Assign) else [anc.target]
                bound = {t.id for t in targets if isinstance(t, ast.Name)}
                if bound & raised_names:
                    return True
        return False

    hits: list[dict] = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Constant) and isinstance(node.value, str)):
            continue
        if id(node) in docstrings:
            continue
        matched = [p for p in SUPERSEDED if p in node.value]
        if not matched:
            continue
        pattern = matched[0]
        enclosing = _enclosing_def(node, parents) or "<module>"
        if _is_refusal(node):
            hits.append(dict(line=node.lineno, pattern=pattern, kind="refusal", detail="refusal",
                             enclosing=enclosing))
            continue
        chain = _ancestors(node, parents)
        # A bare string statement is prose that happens not to be the first in its block.
        if chain and isinstance(chain[0], ast.Expr):
            hits.append(dict(line=node.lineno, pattern=pattern, kind="inert", detail="bare string",
                             enclosing=enclosing))
            continue
        kind, detail = "inert", "named but not read here"
        for anc in chain:
            if isinstance(anc, ast.Call):
                name = _callee(anc)
                if name in READ_FUNCS or name in READ_METHODS:
                    kind, detail = "read", f"{name}()"
                    break
                if name in PROBE_METHODS:
                    kind, detail = "probe", f"{name}()"
                    break
                if name in WRITE_FUNCS:
                    kind, detail = "write", f"{name}()"
                    break
        if kind == "inert":
            for anc in chain:
                if isinstance(anc, (ast.Assign, ast.AnnAssign)):
                    targets = anc.targets if isinstance(anc, ast.Assign) else [anc.target]
                    bound = [t.id for t in targets if isinstance(t, ast.Name)]
                    for b in bound:
                        tainted[b] = (node.lineno, pattern)
                    if bound:
                        detail = f"bound to {', '.join(bound)}"
                    break
                if isinstance(anc, ast.keyword) and anc.arg in {"inputs", "outputs", "sources"}:
                    kind, detail = "declared_input", f"declared as {anc.arg}="
                    break
        hits.append(dict(line=node.lineno, pattern=pattern, kind=kind, detail=detail,
                         enclosing=enclosing))

    # Second pass: a read THROUGH a name the literal was bound to.
    if tainted:
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = _callee(node)
            if name not in READ_FUNCS and name not in READ_METHODS and name not in PROBE_METHODS:
                continue
            used = _names_in(node) & set(tainted)
            # A method call on a name bound from a tainted one, e.g. ``gdir.glob(...)``.
            if isinstance(node.func, ast.Attribute):
                used |= _names_in(node.func) & set(tainted)
            for via in sorted(used):
                bound_at, artifact = tainted[via]
                kind = "probe" if name in PROBE_METHODS else "read"
                hits.append(dict(line=node.lineno, pattern=artifact, kind=kind,
                                 detail=f"{name}() on {via}, bound at line {bound_at}",
                                 enclosing=_enclosing_def(node, parents) or "<module>"))
        # Names bound FROM a tainted name inherit the taint (one level).
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                src = sorted(_names_in(node.value) & set(tainted))
                if src:
                    for t in node.targets:
                        if isinstance(t, ast.Name) and t.id not in tainted:
                            tainted[t.id] = (node.lineno, tainted[src[0]][1])
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = _callee(node)
            if name not in READ_FUNCS and name not in READ_METHODS and name not in PROBE_METHODS:
                continue
            used = set()
            if isinstance(node.func, ast.Attribute):
                used |= _names_in(node.func) & set(tainted)
            used |= _names_in(node) & set(tainted)
            for via in sorted(used):
                bound_at, artifact = tainted[via]
                rec = dict(line=node.lineno, pattern=artifact,
                           kind="probe" if name in PROBE_METHODS else "read",
                           detail=f"{name}() on {via}, bound at line {bound_at}",
                           enclosing=_enclosing_def(node, parents) or "<module>")
                if rec not in hits:
                    hits.append(rec)
    return hits


def _retired_code_paths() -> set[str]:
    """What the repository itself declares retired code, read from ``scripts/retire``.

    Resolution is BY PATH, never by basename: ``dispositions.py`` records at its own head that
    two files here are called ``aggregate_islets.py`` and three ``pins.py``, and that matching on
    the name takes the wrong one.  So the module's simple string constants are folded first and
    each entry's path is resolved against them.
    """
    disp = REPO / "scripts/retire/dispositions.py"
    if not disp.exists():
        return set()
    try:
        tree = ast.parse(disp.read_text())
    except SyntaxError:  # pragma: no cover
        return set()
    env: dict[str, str] = {}
    for node in tree.body:  # module level, in order, so later constants see earlier ones
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and isinstance(node.targets[0], ast.Name):
            text = _path_text(node.value, env)
            if text and "*" not in text:
                env[node.targets[0].id] = text
    # dispositions.py names its paths under ITS OWN root constant, which is an absolute path to
    # one checkout.  Strip that, not this test's REPO: in a worktree the two differ, and stripping
    # the wrong one derived an empty set -- a retired module then looked live.
    root = env.get("REPO", str(REPO))
    out: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        names = {t.id for t in node.targets if isinstance(t, ast.Name)}
        if not any(n.startswith("CODE_TIER") for n in names):
            continue
        for call in ast.walk(node.value):
            if not isinstance(call, ast.Call):
                continue
            expr = next((kw.value for kw in call.keywords if kw.arg == "path"), None)
            if expr is None and call.args:
                expr = call.args[0]
            if expr is None:
                continue
            text = _path_text(expr, env)
            if text.startswith(root + "/"):
                out.add(text[len(root) + 1:])
    return out


def _superseded_reads() -> tuple[int, list[str], list[str]]:
    """``(producers scanned, live findings, findings inside retired producers)``.

    A hit is a read, probe or declared input of a registry artifact.  It is cleared -- by the
    registry, never by this file -- when the record names ``<file>::<function>`` in its
    ``allowed_readers``.  A hit inside a module that ``scripts/retire/dispositions.py`` declares
    retired code (tier 1: verified orphan, no surviving output) is not a live reader and goes in
    the second list, for the strict xfail that keeps it in view until the retire lane acts.
    """
    retired = _retired_code_paths()
    live: list[str] = []
    in_retired: list[str] = []
    scanned = 0
    for tree_root in PRODUCER_TREES:
        if not tree_root.exists():
            continue
        for path in _py_files(tree_root):
            scanned += 1
            rel = _rel(path)
            for hit in _superseded_hits(path):
                if hit["kind"] not in {"read", "probe", "declared_input"}:
                    continue
                record = SUPERSEDED[hit["pattern"]]
                if f"{rel}::{hit['enclosing']}" in record.allowed_readers:
                    continue  # the registry names this reader, and says why
                line = (f"{rel}:{hit['line']}  {hit['kind'].upper()} of {hit['pattern']} "
                        f"[{record.key}, {record.severity}] via {hit['detail']}")
                (in_retired if rel in retired else live).append(line)
    return scanned, live, in_retired


def _superseded_why() -> str:
    return "\n".join(
        f"    [{r.severity}] {r.key}: {', '.join(r.fragments)} -- {r.why}" for r in SUP.REGISTRY)


def test_no_live_reader_of_a_superseded_artifact():
    """A superseded artifact may be refused or described.  It may not be read.

    Every producer under ``scripts/senior/`` and ``core/multiplex5000/scripts/`` is parsed, and
    each string literal naming a registry artifact is classified.  A literal inside a ``raise``,
    inside a ``refuse_*`` function, handed to ``refuse_if_superseded`` inline, or bound to a name
    that is only raised, is the FIX and passes.  A docstring passes.  A comment never reaches the
    AST.  A read the registry names as allowed, with its reason, passes.  What fails is an
    undeclared literal that reaches a read -- at either severity, because a literal read of an
    ANNOUNCE tree is exactly the silent read the announcement exists to prevent.
    """
    scanned, findings, _ = _superseded_reads()
    assert scanned, f"scanned no producers under {[_rel(t) for t in PRODUCER_TREES]}"
    assert not findings, (
        f"{len(findings)} live reference(s) to superseded work, across {scanned} producers.\n"
        f"Superseded, and why (from core/multiplex5000/scripts/typing/superseded.py):\n"
        f"{_superseded_why()}\n"
        f"WHAT TO DO: make each one REFUSE, naming the file, why it is superseded and what to "
        f"read instead -- the pattern is `refuse_superseded` in "
        f"scripts/senior/build_distance_profiles.py, which raises with the counts in the "
        f"message -- or, for a tree the registry ANNOUNCES, hand the path to "
        f"`refuse_if_superseded` inline so the read is loud.  A silent fallback is the defect; a "
        f"loud failure is the fix.  Do not delete any data root, and do not add an exclusion "
        f"here: a reader that must read anyway is recorded in the registry record's "
        f"allowed_readers, with why.\n  "
        + "\n  ".join(sorted(findings)))


def test_no_retired_producer_still_reads_a_superseded_artifact():
    """A retired producer that still opens superseded work can regenerate it, silently."""
    scanned, _, findings = _superseded_reads()
    assert scanned
    assert not findings, (
        f"{len(findings)} reference(s) to superseded work inside modules the repository declares "
        f"RETIRED, still on disk and runnable.  WHAT TO DO: retire them -- delete, or make the "
        f"read refuse by name -- per their disposition in scripts/retire/dispositions.py.\n  "
        + "\n  ".join(sorted(findings)))


#: The prose a retirement must be declared in before the registry may enforce it.
DECLARATION_SOURCES = ("CLAUDE.md", "scripts/CLAUDE.md", "core/multiplex5000/README.md",
                       "scripts/retire/dispositions.py", "README.md")

#: A glob the prose documents, as a code span made only of path characters and ``*`` -- the root
#: ``CLAUDE.md`` retires "flat `cells_redsea*` / `restore*` roots" and ``scripts/CLAUDE.md``
#: repeats them as ``data/cells_redsea*`` and ``data/restore*``.  The span is the whole match,
#: so a sentence that happens to sit between two code spans on one line is never read as one.
DOCUMENTED_GLOB = re.compile(r"`([A-Za-z0-9_./-]*\*[A-Za-z0-9_./*-]*)`")


def _documented_globs(text: str) -> set[str]:
    """The globs the prose documents, keeping only those whose basename begins with a literal.

    ``*.csv``, ``cells/*.csv`` and a hypothetical ``data/*`` constrain nothing about the name of
    a root, so honouring them would let the registry name any root at all and call it declared.
    A glob counts only when the text before its first ``*``, after the last ``/``, is non-empty.
    """
    globs: set[str] = set()
    for glob in DOCUMENTED_GLOB.findall(text):
        if glob.rsplit("/", 1)[-1].split("*", 1)[0]:
            globs.add(glob)
    return globs


def _declared_as(pattern: str, text: str, globs: set[str]) -> str | None:
    """How the prose declares ``pattern`` -- its literal basename, or a documented glob -- or None.

    A glob with a ``/`` is matched against the whole repo-relative pattern (``data/restore*``
    against ``data/restore_pair_validation``); one without is matched against the basename.
    """
    rel = pattern.strip("/")
    base = rel.rsplit("/", 1)[-1]
    if base in text:
        return f"basename {base!r}"
    for glob in sorted(globs):
        if fnmatch.fnmatchcase(rel if "/" in glob else base, glob):
            return f"documented glob `{glob}`"
    return None


@pytest.mark.parametrize("record", SUP.REGISTRY, ids=[r.key for r in SUP.REGISTRY])
def test_the_superseded_table_is_itself_declared_in_the_repository(record):
    """The registry may not invent a retirement.  EVERY fragment and alias is declared.

    Declared means the basename appears in one of :data:`DECLARATION_SOURCES`, or the path
    matches a glob those sources document.  A record's key declares nothing, and one declared
    fragment does not clear its siblings: the any-of reading introduced in 5137741 passed
    ``restore_roots`` on ``restore_gated_redsea`` alone while ``data/restore_gated_redsea_extra``
    and ``data/restore_pair_validation`` were spelled out nowhere -- they are declared by the
    maintainer's ``restore*`` glob, which is what this test honours instead.
    """
    sources = [REPO / s for s in DECLARATION_SOURCES]
    text = "\n".join(p.read_text(errors="replace") for p in sources if p.exists())
    globs = _documented_globs(text)
    patterns = (*record.fragments, *record.aliases)
    undeclared = [p for p in patterns if _declared_as(p, text, globs) is None]
    assert not undeclared, (
        f"registry record {record.key!r} names {len(undeclared)} of {len(patterns)} pattern(s) "
        f"that are declared nowhere in {list(DECLARATION_SOURCES)} -- neither by basename nor "
        f"by a documented glob {sorted(globs)}: {undeclared}.  Either the declaration was "
        f"removed or the registry invented a retirement.  Declare it in one of those sources; "
        f"do not weaken this check to any-of.")


# ----------------------------------------- 3. an orphan is never a live input ----
def _artifact_calls(path: pathlib.Path) -> tuple[list[tuple[str, int]], set[str]]:
    """``(reads_under_an_artifact_root, everything_this_module_can_write)``.

    Writes are followed the same way reads are, because a producer that builds its destination in
    one statement and writes it in the next is still a producer:

    * a literal path inside the write call itself;
    * a literal path bound to a name that a write call then uses
      (``shifts_out = args.shifts_out or (panel.OUT / "shifts.json")`` ... ``.write_text``);
    * an ``argparse`` ``default=`` whose ``args.<dest>`` a write call then uses
      (``--out`` defaulting to ``donor_metadata.csv`` ... ``args.out.to_csv``).

    Without those two, a module that parameterises its destination reads as producing nothing,
    and every artifact it makes looks like an orphan.  That would be noise, and a guard that
    cries wolf gets muted.
    """
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError:
        return [], set()
    reads: list[tuple[str, int]] = []
    writes: set[str] = set()

    def _keys(expr: ast.AST) -> str:
        key = _path_text(expr).strip("/")
        return key if key and "." in key.rsplit("/", 1)[-1] else ""

    # Names bound to a literal path, and argparse dests whose default is a literal path.
    bound: dict[str, str] = {}
    arg_dest: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            key = _keys(node.value)
            if key:
                for target in node.targets:
                    if isinstance(target, ast.Name):
                        bound[target.id] = key
            else:  # `a or b` / `a if c else b`: take whichever side carries the literal
                for sub in ast.walk(node.value):
                    key = _keys(sub)
                    if key:
                        for target in node.targets:
                            if isinstance(target, ast.Name):
                                bound.setdefault(target.id, key)
                        break
        if isinstance(node, ast.Call) and _callee(node) == "add_argument":
            flags = [a.value for a in node.args
                     if isinstance(a, ast.Constant) and isinstance(a.value, str)]
            default = next((kw.value for kw in node.keywords if kw.arg == "default"), None)
            dest = next((kw.value.value for kw in node.keywords
                         if kw.arg == "dest" and isinstance(kw.value, ast.Constant)), None)
            if default is not None:
                key = _keys(default)
                if key:
                    for flag in flags:
                        name = dest or flag.lstrip("-").replace("-", "_")
                        arg_dest[name] = key

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = _callee(node)
        if name is None:
            continue
        exprs: list[ast.AST] = list(node.args)
        if isinstance(node.func, ast.Attribute):
            exprs.append(node.func.value)
        for expr in exprs:
            if name in WRITE_FUNCS:
                key = _keys(expr)
                if key:
                    writes.add(key)
                for sub in ast.walk(expr):
                    if isinstance(sub, ast.Name) and sub.id in bound:
                        writes.add(bound[sub.id])
                    if isinstance(sub, ast.Attribute) and isinstance(sub.value, ast.Name) \
                            and sub.value.id == "args" and sub.attr in arg_dest:
                        writes.add(arg_dest[sub.attr])
            elif name in READ_FUNCS or name in READ_METHODS:
                key = _keys(expr)
                if not key:
                    continue
                # Only artifacts rooted at a module's declared output root are in scope here.
                roots = {n.attr for n in ast.walk(expr) if isinstance(n, ast.Attribute)}
                roots |= _names_in(expr)
                if "OUT" in roots or "data/development/typing" in key:
                    reads.append((key, node.lineno))
    return reads, writes


def _declared_producerless(key: str) -> str | None:
    """The registry's PRODUCERLESS entry naming this read key, if any."""
    for path in SUP.PRODUCERLESS:
        if path == key or path.endswith("/" + key):
            return path
    return None


#: What the orphan guard found among the registry's declared producerless inputs, for the strict
#: xfail below to report.  Filled by test_every_artifact_read_as_an_input_has_a_producer.
_ORPHANS_DECLARED_IN_THE_REGISTRY: list[str] = []


def _tracked_basenames() -> set[str]:
    """Files git already carries.  A checked-in input HAS a producer: source control."""
    done = _git("ls-files")
    if done.returncode != 0:  # pragma: no cover - a repo without git is not this test's business
        return set()
    return {line.rsplit("/", 1)[-1] for line in done.stdout.splitlines() if line.strip()}


def test_every_artifact_read_as_an_input_has_a_producer():
    """A file read as an input must have something in the repository that can re-make it.

    An orphan cannot be re-made, so whatever sits on disk is frozen at whatever generation wrote
    it and nothing can say which -- the same shape as the superseded calls table, whose producer
    was lost before ``subpanel.py`` was written back from ``RECONSTRUCTION.md``.  The property
    asserted is the general one, not any single case.

    Three things count as a producer, all derived: a write call anywhere in the repository that
    names the same file, a write into a directory on the read path, and git itself -- a file that
    is checked in beside its reader is not an orphan.
    """
    producers: set[str] = set()
    for tree_root in (REPO / "scripts", REPO / "core/multiplex5000", REPO / "apps",
                      REPO / "analysis"):
        if not tree_root.exists():
            continue
        for path in sorted(p for p in tree_root.rglob("*.py")
                           if "__pycache__" not in p.parts and ".venv" not in p.parts):
            producers |= _artifact_calls(path)[1]
    produced_basenames = {k.rsplit("/", 1)[-1] for k in producers} | _tracked_basenames()
    assert produced_basenames, "found no producers at all; the derivation is broken"

    orphans: set[str] = set()
    declared: set[str] = set()
    for tree_root in PRODUCER_TREES:
        if not tree_root.exists():
            continue
        for path in _py_files(tree_root):
            for key, line in _artifact_calls(path)[0]:
                base = key.rsplit("/", 1)[-1]
                if base in produced_basenames:
                    continue
                if any(key.startswith(w.rstrip("/") + "/") or w.startswith(key + "/")
                       for w in producers):
                    continue
                # A superseded artifact is recorded, refused and guarded by the registry and by
                # test_no_live_reader_of_a_superseded_artifact; policing it here again would make
                # the same finding under two names.
                if SUP.match(REPO / key) is not None:
                    continue
                finding = (f"{_rel(path)}:{line}  reads {key} -- NOTHING in the repository "
                           f"writes {base}, and it is not checked in")
                (declared if _declared_producerless(key) else orphans).add(finding)
    _ORPHANS_DECLARED_IN_THE_REGISTRY.clear()
    _ORPHANS_DECLARED_IN_THE_REGISTRY.extend(sorted(declared))
    assert not orphans, (
        f"{len(orphans)} artifact(s) are read as a live input with no producer in the "
        f"repository.  An orphan cannot be re-made, so it is frozen at whatever generation "
        f"happened to write it and nothing downstream can say which -- that is the failure "
        f"scripts/retire/dispositions.py names at its code-tier banner, and the one "
        f"umap/refinement_log.json was in.  WHAT TO DO: either commit the producer, or make the "
        f"read REFUSE by name and say what to read instead.  Do not delete the file.\n  "
        + "\n  ".join(sorted(orphans)))


@pytest.mark.parametrize("path", sorted(SUP.PRODUCERLESS), ids=lambda p: p.rsplit("/", 1)[-1])
def test_a_declared_producerless_input_is_still_producerless(path):
    """The declaration may not outlive its truth: once the file is tracked or written by
    something in the repository, the PRODUCERLESS entry must go."""
    base = path.rsplit("/", 1)[-1]
    assert base not in _tracked_basenames(), (
        f"{path} is now tracked by git, which is a producer; remove it from "
        f"superseded.PRODUCERLESS")
    writers: set[str] = set()
    for tree_root in (REPO / "scripts", REPO / "core/multiplex5000", REPO / "apps",
                      REPO / "analysis"):
        if not tree_root.exists():
            continue
        for module in sorted(p for p in tree_root.rglob("*.py")
                             if "__pycache__" not in p.parts and ".venv" not in p.parts):
            if any(w.rsplit("/", 1)[-1] == base for w in _artifact_calls(module)[1]):
                writers.add(_rel(module))
    assert not writers, (
        f"{path} now has a producer ({sorted(writers)}); remove it from superseded.PRODUCERLESS")


@pytest.mark.xfail(strict=True, reason=(
    "superseded.PRODUCERLESS declares data/development/typing/subpanel/subpanel_log.json: read "
    "by subpanel.py::_log and trajectory_singlecell_features.py::feature_markers, 688 bytes, "
    "mtime 2026-09-11 14:10, written by nothing in the repository, ignored by .gitignore's "
    "anchored /data/.  A ruling is needed -- commit it as a tracked sidecar, or give it a "
    "producer -- and it is not made in a test.  This xfail fails the moment the declaration is "
    "emptied, so the marker is removed with the ruling."))
def test_no_live_input_is_declared_producerless():
    """Every live input has a producer, and the registry declares no exception to that."""
    test_every_artifact_read_as_an_input_has_a_producer()
    assert not SUP.PRODUCERLESS and not _ORPHANS_DECLARED_IN_THE_REGISTRY, (
        f"declared producerless, awaiting a ruling: {sorted(SUP.PRODUCERLESS)}\n  "
        + "\n  ".join(_ORPHANS_DECLARED_IN_THE_REGISTRY))


# -------------------------------- 4. the document the maintainer named is complete ----
def _refinement_concepts() -> list[str]:
    """Derived: identifiers ``refine.py`` defines at module level AND ``REFINEMENT.md`` names in
    a code span.  Two committed artifacts have to agree that a name matters before this file
    calls it a concept, so plumbing drops out and the set grows with the method."""
    refine = TYPING / "refine.py"
    doc = TYPING / "REFINEMENT.md"
    if not (refine.exists() and doc.exists()):
        pytest.skip(f"the refinement is not on disk: {_rel(refine)} / {_rel(doc)}")
    tree = ast.parse(refine.read_text())
    defined: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            defined.add(node.name)
        elif isinstance(node, ast.Assign):
            defined |= {t.id for t in node.targets if isinstance(t, ast.Name)}
    spans = set(re.findall(r"`([^`\n]{2,60})`", doc.read_text()))
    idents = {s for s in spans if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", s)}
    return sorted(defined & idents)


def _names_concept(doc: str, concept: str) -> bool:
    """The identifier itself, or its plain-English spelling with the underscores opened out."""
    if re.search(r"\b" + re.escape(concept) + r"\b", doc):
        return True
    if "_" in concept:
        spaced = concept.replace("_", " ")
        return bool(re.search(r"\b" + re.escape(spaced) + r"\b", doc, re.IGNORECASE))
    return False


def test_reconstruction_names_every_concept_the_refinement_rests_on():
    """``RECONSTRUCTION.md`` is the document the maintainer said to read rather than re-derive.

    A reader who does that and finds a concept absent re-derives it, and a re-derivation is
    exactly how a second generation gets made beside the first.
    """
    recon = TYPING / "RECONSTRUCTION.md"
    assert recon.exists(), f"{_rel(recon)} is missing -- it is the document the maintainer named"
    concepts = _refinement_concepts()
    assert concepts, (
        "derived no concepts from refine.py x REFINEMENT.md; the derivation is broken, which "
        "would make this test pass by measuring nothing")
    text = recon.read_text()
    missing = [c for c in concepts if not _names_concept(text, c)]
    assert not missing, (
        f"{len(missing)} of {len(concepts)} concepts the committed refinement rests on are NOT "
        f"mentioned in {_rel(recon)}: {missing}.\n"
        f"Each is defined at module level in scripts/typing/refine.py and named in a code span "
        f"in REFINEMENT.md, so both committed artifacts already treat it as load-bearing.  WHAT "
        f"TO DO: add each to the refinement section of RECONSTRUCTION.md, marked [CHAIN NOW], "
        f"so a reader who follows the maintainer's instruction to read that file instead of "
        f"re-deriving the method is not sent back to the source.")


def test_reconstruction_names_every_artifact_the_refinement_touches():
    """Every file ``refine.py`` writes, reads or refuses has to be named where the method is
    documented -- otherwise the next reader cannot tell the refusal from the artifact."""
    refine = TYPING / "refine.py"
    recon = TYPING / "RECONSTRUCTION.md"
    if not refine.exists():
        pytest.skip(f"{_rel(refine)} is not on disk")
    reads, writes = _artifact_calls(refine)
    tree = ast.parse(refine.read_text())
    docstrings = _docstring_ids(tree)
    named: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                and id(node) not in docstrings:
            for token in re.findall(r"[\w./-]+\.(?:parquet|json|csv|geojson)", node.value):
                named.add(token.rsplit("/", 1)[-1])
    named |= {k.rsplit("/", 1)[-1] for k, _ in reads} | {k.rsplit("/", 1)[-1] for k in writes}
    text = recon.read_text()
    missing = sorted(n for n in named if n not in text)
    assert not missing, (
        f"{len(missing)} artifact(s) that scripts/typing/refine.py names in code are absent from "
        f"{_rel(recon)}: {missing}.  WHAT TO DO: name each one where the refinement is "
        f"documented, and say whether it is written, read or REFUSED -- a reader who cannot tell "
        f"a refused file from a live input is one step from reading the superseded generation.")
