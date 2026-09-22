#!/usr/bin/env python3
"""Re-verify every declared disposition and write the retirement manifest.

Nothing is retired on the strength of the table in dispositions.py. Each entry names a CHECK; this
script runs it against the filesystem now and records pass/fail with the detail. `retire.py` acts
only on entries whose check PASSED, so a table that has drifted from disk fails closed.

    python scripts/retire/build_manifest.py            # writes data/ARCHIVE_2026-09-15/MANIFEST.csv
    python scripts/retire/build_manifest.py --quiet
    python scripts/retire/build_manifest.py --kind code --out /tmp/x.csv  # verify elsewhere

Code rows carry two tiers that are not tiers: "LIBRARY" (declared reached by live code; the
`library` check asserts presence) and "CANDIDATE" (no live reference measured; the `code_orphan`
check is run so the row shows today's measurement). `retire.py` selects tiers 1 and 2 only.
"""
from __future__ import annotations

import argparse
import ast
import csv
import functools
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dispositions as D  # noqa: E402

ARCHIVE = Path(D.REPO) / "data" / "ARCHIVE_2026-09-15"
SHA_LIMIT = 100 * 1024 * 1024  # don't hash multi-GB artifacts; record the skip instead


def walk(p: Path):
    if p.is_file():
        yield p
        return
    for root, _dirs, files in os.walk(p):
        for f in files:
            yield Path(root) / f


def measure(p: Path) -> tuple[int, int, str]:
    """(bytes, n_files, iso mtime) -- mtime is the newest file for a directory."""
    if not p.exists():
        return 0, 0, ""
    total = n = 0
    newest = 0.0
    for f in walk(p):
        try:
            st = f.lstat()
        except OSError:
            continue
        total += st.st_size
        n += 1
        newest = max(newest, st.st_mtime)
    if p.is_file():
        newest = p.lstat().st_mtime
    return total, n, datetime.fromtimestamp(newest, timezone.utc).isoformat(timespec="seconds")


def sha256(p: Path) -> str:
    if not p.is_file():
        return "dir"
    if p.stat().st_size > SHA_LIMIT:
        return "skipped_large"
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rel_files(root: Path) -> set[str]:
    return {str(f.relative_to(root)) for f in walk(root)}


# ---------------------------------------------------------------- checks ----
def check_plain(p: Path, arg):
    return (True, "exists") if p.exists() else (False, "MISSING")


def check_orphan_temp(p: Path, arg):
    ev = Path(arg)
    if not p.exists():
        return False, "MISSING"
    if not ev.exists():
        return False, f"completed product absent: {ev}"
    return True, f"completed product present: {ev} ({ev.stat().st_size:,} B)"


def check_subset_dir(p: Path, arg):
    other = Path(arg)
    if not p.exists():
        return False, "MISSING"
    if not other.exists():
        return False, f"replacement absent: {other}"
    mine, theirs = rel_files(p), rel_files(other)
    missing = mine - theirs
    if missing:
        return False, f"{len(missing)} path(s) not in replacement, e.g. {sorted(missing)[:3]}"
    return True, f"all {len(mine)} paths present in replacement"


def check_names_exist_under(p: Path, arg):
    other = Path(arg)
    if not p.exists():
        return False, "MISSING"
    if not other.exists():
        return False, f"replacement absent: {other}"
    mine = {f.name for f in walk(p)}
    theirs = {f.name for f in walk(other)} - mine if other == p.parent else {f.name for f in walk(other)}
    if other == p.parent:  # parent contains p; exclude p's own files from the comparison
        theirs = {f.name for f in walk(other) if p not in f.parents and f.parent != p}
    missing = mine - theirs
    if missing:
        return False, f"{len(missing)} name(s) absent from replacement, e.g. {sorted(missing)[:3]}"
    return True, f"all {len(mine)} basenames present in replacement"


def check_identical_to(p: Path, arg):
    other = Path(arg)
    if not (p.exists() and other.exists()):
        return False, "MISSING one side"
    a, b = sha256(p), sha256(other)
    return (a == b, f"sha256 {'match' if a == b else 'DIFFER'} {a[:16]} vs {b[:16]}")


def _records(p: Path):
    j = json.loads(p.read_text())
    return j if isinstance(j, list) else [j]


def check_donor_local_cover(p: Path, arg):
    if not p.exists():
        return False, "MISSING"
    valis = Path(D.REPO) / "data/development/restain/valis"
    donors = [str(r.get("donor_id")) for r in _records(p) if r.get("donor_id")]
    if not donors:
        return False, "no donor_id in records"
    missing = [d for d in donors if not (valis / f"donor_id={d}" / arg).exists()]
    if missing:
        return False, f"donor(s) without {arg}: {missing}"
    return True, f"all {len(donors)} donor(s) {donors} have donor-local {arg}"


def check_identical_keys_subset(p: Path, arg):
    other = Path(arg)
    if not (p.exists() and other.exists()):
        return False, "MISSING one side"
    mine = {str(r.get("donor_id")): r for r in _records(p)}
    theirs = {str(r.get("donor_id")): r for r in _records(other)}
    for d, rec in mine.items():
        if d not in theirs:
            return False, f"donor {d} absent from replacement"
        extra = set(rec.get("hormone", rec)) - set(theirs[d].get("hormone", theirs[d]))
        if extra:
            return False, f"donor {d} has key(s) the replacement lacks: {sorted(extra)}"
    return True, f"all {len(mine)} donor(s) covered; replacement keys are a superset"


def check_sibling_name_under(p: Path, arg):
    """A single file whose basename exists under `arg` (a sibling directory)."""
    other = Path(arg)
    if not p.exists():
        return False, "MISSING"
    if not other.exists():
        return False, f"replacement dir absent: {other}"
    hit = [f for f in walk(other) if f.name == p.name]
    if not hit:
        return False, f"{p.name} not found under {other}"
    return True, f"{p.name} present at {hit[0]} ({hit[0].stat().st_size:,} B)"


# ------------------------------------------------------- code-tier checks ----
#
# A code entry is retired only when NOTHING outside the retiring set reaches it. "Reaches" is
# split three ways, deliberately:
#
#   HARD   a call path or a path literal: a Python `import`/`from` of the module (read with
#          `ast`, so a mention inside a string or a comment cannot masquerade as one), an R
#          `source()`/`sys.source()`, an interpreter invocation, or a path-qualified occurrence
#          in a non-Markdown file -- a test that does readLines("apps/R/ieCore/R/ie_holds.R")
#          breaks when the file moves just as surely as an import does. One hard hit FAILS the
#          check and `retire.py` refuses the entry.
#   DOCS   an occurrence in a .md file. Reported with file:line so the retiring commit fixes the
#          instruction, never blocking: prose is not a caller.
#   SOFT   a bare basename with no path in front of it, in a file where that basename is not
#          unique in the repo, or an occurrence on a comment line. Reported, never blocking.
#
# Resolution is by PATH, not by basename. `scripts/senior/aggregate_islets.py` and
# `core/phenocycler/scripts/typing/aggregate_islets.py` are different files, and the second is a
# live node in `scripts/senior/build_layers.py`'s DAG; a basename match would have retired the
# wrong one.
#
# `data/` and `scripts/retire/` are not scanned: a provenance note under data/app_data records how
# a file was made, which is history, and this table naming a path is not a caller either.
SCAN_EXTS = {".py", ".R", ".r", ".sh", ".ini", ".toml", ".yml", ".yaml", ".json",
             ".cfg", ".Rmd", ".groovy", ".md"}
SCAN_SKIP = {".git", ".venv", "external", "node_modules", "__pycache__", "data",
             "renv", ".Rproj.user", "site-packages", ".pytest_cache"}
#: This lane's own table names every candidate by path; it is a declaration, not a call site.
SCAN_SKIP_PATHS = (f"{D.REPO}/scripts/retire",)
EXEC_MARK = ("subprocess", "Popen", "check_call", "check_output", "os.system",
             "Rscript ", "python ", "python3 ", "bash ", "sh ")


@functools.lru_cache(maxsize=1)
def scan_files() -> tuple[Path, ...]:
    """Every file in the repo that could hold a reference, walked once."""
    root = Path(D.REPO)
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SCAN_SKIP and not d.endswith(".egg-info")]
        if any(dirpath == s or dirpath.startswith(s + "/") for s in SCAN_SKIP_PATHS):
            continue
        for name in filenames:
            if Path(name).suffix in SCAN_EXTS:
                out.append(Path(dirpath) / name)
    return tuple(sorted(out))


@functools.lru_cache(maxsize=1)
def basename_index() -> dict:
    """basename -> how many files in the repo carry it. A bare mention of a name that is not
    unique cannot be attributed, so it is reported rather than acted on."""
    counts: dict[str, int] = {}
    root = Path(D.REPO)
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SCAN_SKIP and not d.endswith(".egg-info")]
        for name in filenames:
            counts[name] = counts.get(name, 0) + 1
    return counts


def _relrepo(p: Path) -> str:
    return str(p).replace(D.REPO + "/", "")


def _targets(p: Path) -> tuple[set[str], dict]:
    """(python module stems this path removes, {token -> repo-relative path that owns it})."""
    stems: set[str] = set()
    owner: dict[str, str] = {}
    files = [p] if p.is_file() else [f for f in walk(p) if f.suffix in SCAN_EXTS]
    for f in files:
        owner.setdefault(f.name, _relrepo(f))
        if f.suffix == ".py" and f.name != "__init__.py":
            stems.add(f.stem)
    if p.is_dir():
        owner[p.name] = _relrepo(p)
    return stems, owner


@functools.lru_cache(maxsize=None)
def _py_doc_lines(f: Path) -> frozenset:
    """Line numbers covered by module/class/function docstrings.

    A docstring that names a path is documentation, not a call site. Without this, this very
    lane's own `rebuild_adata_ins_root.py` -- whose opening paragraph explains WHY the two legacy
    producers are not its ancestors -- would read as a live caller of both of them.
    """
    try:
        tree = ast.parse(f.read_text(errors="replace"))
    except (SyntaxError, ValueError, OSError):
        return frozenset()
    lines: set[int] = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        doc = node.body[0] if node.body else None
        if (isinstance(doc, ast.Expr) and isinstance(doc.value, ast.Constant)
                and isinstance(doc.value.value, str)):
            lines.update(range(doc.lineno, (doc.end_lineno or doc.lineno) + 1))
    return frozenset(lines)


@functools.lru_cache(maxsize=None)
def _sys_path_text(f: Path) -> str:
    """The raw text of every `sys.path` line in `f` -- how the file says where its modules live."""
    try:
        return "\n".join(ln for ln in f.read_text(errors="replace").splitlines()
                          if "sys.path" in ln)
    except OSError:
        return ""


def _import_resolves_here(f: Path, stem: str, module_dir: Path) -> bool:
    """Does `f`'s `import <stem>` reach THIS module, or another file of the same name?

    Three lanes carry a `pins.py` (broadtype, report, xenium_typing) and two carry an
    `aggregate_islets.py`. An import matched on the stem alone retires the wrong one, so the
    import is resolved the way the file itself resolves it: same directory, or a `sys.path` line
    that names the module's directory, or a stem that is unique in the repo.
    """
    if f.parent == module_dir:
        return True
    text = _sys_path_text(f)
    name = module_dir.name
    if name and (f'"{name}"' in text or f"'{name}'" in text or f"/{name}" in text):
        return True
    if _relrepo(module_dir) in text:
        return True
    return basename_index().get(stem + ".py", 0) == 1


@functools.lru_cache(maxsize=None)
def _py_imports(f: Path) -> frozenset:
    """Top-level module names `f` imports, from the AST -- never from the raw text."""
    try:
        tree = ast.parse(f.read_text(errors="replace"))
    except (SyntaxError, ValueError, OSError):
        return frozenset()
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            found.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            found.add(node.module.split(".")[0])
    return frozenset(found)


def _inside(f: Path, roots) -> bool:
    s = str(f)
    return any(s == r or s.startswith(r.rstrip("/") + "/") for r in roots)


@functools.lru_cache(maxsize=None)
def _token_re(name: str):
    return re.compile(r"(?P<pre>[\w./{}$-]*/)?" + re.escape(name) + r"(?![\w])")


def _code_orphan(p: Path, internal, scope: str):
    """The reference walk. ``internal`` is the set of roots whose references do not count;
    ``scope`` names it in the verdict."""
    if not p.exists():
        return False, "MISSING"
    stems, owner = _targets(p)
    if not owner:
        return False, "no scannable file under it"
    counts = basename_index()
    hard, docs, soft = [], [], []
    for f in scan_files():
        if _inside(f, [str(p)]) or _inside(f, internal):
            continue
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        shared = stems & _py_imports(f)
        if f.suffix == ".py" and shared:
            resolved = sorted(s for s in shared
                              if _import_resolves_here(f, s, Path(D.REPO) / Path(owner[s + ".py"]).parent))
            if resolved:
                hard.append(f"{_relrepo(f)}: imports {resolved}")
                continue
            soft.append(f"{_relrepo(f)}: imports {sorted(shared)} -- resolves to another lane's "
                        f"module of the same name")
            continue
        if not any(tok in text for tok in owner):
            continue
        for i, line in enumerate(text.splitlines(), 1):
            for tok, tok_owner in owner.items():
                for m in _token_re(tok).finditer(line):
                    pre = m.group("pre") or ""
                    if pre:
                        if not tok_owner.endswith(pre + tok):
                            continue          # a different file that shares the basename
                        qualified = True
                    else:
                        if counts.get(tok, 0) != 1:
                            soft.append(f"{_relrepo(f)}:{i} (bare '{tok}', not unique)")
                            continue
                        qualified = False
                    where = f"{_relrepo(f)}:{i}: {line.strip()[:100]}"
                    if f.suffix == ".md" or (f.suffix == ".py" and i in _py_doc_lines(f)):
                        docs.append(where)
                    elif line.lstrip().startswith("#"):
                        soft.append(f"{_relrepo(f)}:{i} (comment)")
                    elif qualified or any(k in line for k in EXEC_MARK) or "source(" in line:
                        hard.append(where)
                    else:
                        soft.append(f"{_relrepo(f)}:{i}")
    tail = (f"; {len(docs)} doc instruction(s) to update: {docs[:3]}" if docs else "") + \
           (f"; {len(soft)} prose mention(s)" if soft else "")
    if hard:
        return False, f"{len(hard)} LIVE reference(s) {scope}: {hard[:3]}{tail}"
    return True, f"no live import, source(), exec or path literal {scope}{tail}"


def check_code_orphan(p: Path, arg):
    return _code_orphan(p, D.CODE_RETIRING, "outside the retiring set")


def check_candidate_orphan(p: Path, arg):
    """A CODE_CANDIDATE row: the walk the maintainer's decision rests on. A reference from
    another candidate is reported as internal, because the decision is about the set -- keep
    eval_detection.py and context.py stays with it -- and the row's own evidence names it."""
    return _code_orphan(p, D.CODE_RETIRING + D.CODE_CANDIDATE_PATHS,
                        "outside the retiring and candidate sets")


CONFIG_INI = Path(D.REPO) / "core/phenocycler/config.ini"


def _config_binding(ini: Path, key: str) -> tuple[int, str] | None:
    """``(line, value)`` of ``key = value`` in ``ini``, or None when the key is absent."""
    for i, line in enumerate(ini.read_text(errors="replace").splitlines(), 1):
        m = re.match(r"\s*" + re.escape(key) + r"\s*=(.*)$", line)
        if m:
            return i, m.group(1).strip()
    return None


def check_code_orphan_and_config_blank(p: Path, arg):
    """The orphan walk, AND the config.ini key ``arg`` must be blank.

    For the Xenium pancreas comparison apparatus the walk is not enough: its importers are all
    inside the retiring set, but core/phenocycler/config.ini binds the RELEASE DIRECTORY it
    writes, and that binding is what carries its class column into cells_xen and two build_layers
    nodes. A directory binding is invisible to a reference walk over .py names, so the
    prerequisite is checked by name here and the entry fails closed until Step 5 blanks it.
    """
    ok, detail = check_code_orphan(p, None)
    if not ok:
        return ok, detail
    found = _config_binding(CONFIG_INI, arg)
    if found is None:
        return True, f"{detail}; {arg} is not bound in {_relrepo(CONFIG_INI)}"
    line, value = found
    if value:
        return False, (f"prerequisite not met: {_relrepo(CONFIG_INI)}:{line} still binds "
                       f"{arg} = {value}; {detail}")
    return True, f"{detail}; {_relrepo(CONFIG_INI)}:{line} has {arg} blank"


def check_no_reference_anywhere(p: Path, arg):
    if not p.exists():
        return False, "MISSING"
    hits = []
    for f in scan_files():
        try:
            text = f.read_text(errors="replace")
        except OSError:
            continue
        if p.name in text or p.stem in text:
            hits.append(_relrepo(f))
    if hits:
        return False, f"referenced by {len(hits)} file(s): {hits[:4]}"
    return True, f"'{p.stem}' appears in no scanned file -- no producer and no consumer"


def check_library(p: Path, arg):
    """A declared-library entry. Asserted so the guard is real, never acted on.

    Absence IS a failure: every library was measured present on main on 2026-09-22 (the
    merge/hpg-sameslide additions this check once tolerated as "not landed yet" have landed).
    A library that is gone means the table drifted, and a drifted table is corrected, not acted
    on.
    """
    if not p.exists():
        return False, "MISSING -- a declared library is gone; correct the table"
    n = 1 if p.is_file() else sum(1 for _ in walk(p))
    return True, f"present, {n} file(s)"


CHECKS = {
    "plain": check_plain,
    "code_orphan": check_code_orphan,
    "candidate_orphan": check_candidate_orphan,
    "code_orphan_and_config_blank": check_code_orphan_and_config_blank,
    "no_reference_anywhere": check_no_reference_anywhere,
    "library": check_library,
    "sibling_name_under": check_sibling_name_under,
    "orphan_temp": check_orphan_temp,
    "subset_dir": check_subset_dir,
    "names_exist_under": check_names_exist_under,
    "identical_to": check_identical_to,
    "donor_local_cover": check_donor_local_cover,
    "identical_keys_subset": check_identical_keys_subset,
}

FIELDS = ["kind", "tier", "path", "exists", "bytes", "n_files", "mtime_utc", "sha256", "check",
          "check_arg", "verified", "verify_detail", "superseded_by", "evidence", "produces",
          "prereq", "blocked_by"]


def row(tier, path, check=None, arg=None, superseded_by="", evidence="", kind="data",
        produces=(), prereq="", blocked_by=()):
    p = Path(path)
    b, n, m = measure(p)
    verified, detail = ("", "not checked (untouched tier)")
    if check:
        fn = CHECKS[check]
        ok, detail = fn(p, arg)
        verified = "PASS" if ok else "FAIL"
    return dict(kind=kind, tier=tier, path=str(p), exists=p.exists(), bytes=b, n_files=n,
                mtime_utc=m, sha256=sha256(p) if p.is_file() else "dir", check=check or "",
                check_arg=arg or "", verified=verified, verify_detail=detail,
                superseded_by=superseded_by, evidence=evidence,
                produces=" ".join(produces), prereq=prereq,
                blocked_by=" ".join(blocked_by))


def data_rows() -> list[dict]:
    rows = []
    for e in D.TIER1:
        rows.append(row(1, e["path"], e["check"], e.get("arg"), e["superseded_by"], e["evidence"]))
    for p, why in D.KEEP_BESIDE_TIER1:
        rows.append(row("1-KEEP", p, "plain", None, "-", why))
    for e in D.TIER2:
        r = row(2, e["path"], "plain", None, e["superseded_by"], e["evidence"])
        r["check"] = f"carve:{e['carve']}"
        rows.append(r)
    for p, why in D.TIER3:
        rows.append(row(3, p, None, None, "-", why))
    for p, why in D.TIER4:
        rows.append(row(4, p, None, None, "-", why))
    return rows


def code_rows() -> list[dict]:
    """The 2026-09-18 code tier. `kind` keeps it out of the data totals."""
    rows = []
    for e in D.CODE_TIER1 + D.CODE_TIER2:
        rows.append(row(e["tier"], e["path"], e["check"], e["arg"], e["superseded_by"],
                        e["evidence"], kind="code", produces=e["produces"], prereq=e["prereq"],
                        blocked_by=e["blocked_by"]))
    for p, why in D.CODE_LIBRARY:
        rows.append(row("LIBRARY", p, "library", None, "-", why, kind="code"))
    for p, why in D.CODE_CANDIDATE:
        # The orphan walk runs on a candidate too, so the manifest carries today's measurement
        # beside the recorded evidence. A candidate is on no tier: retire.py never selects it.
        rows.append(row("CANDIDATE", p, "candidate_orphan", None, "-", why, kind="code"))
    for p, why in D.CODE_TIER3:
        rows.append(row(3, p, None, None, "-", why, kind="code"))
    for p, why in D.CODE_TIER4:
        rows.append(row(4, p, None, None, "-", why, kind="code"))
    for e in D.ARTIFACT_TIER2:
        rows.append(row(e["tier"], e["path"], e["check"], None, e["superseded_by"],
                        e["evidence"], kind="artifact", produces=e["produces"],
                        prereq=e["prereq"], blocked_by=e["blocked_by"]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--kind", choices=("data", "code", "all"), default="all",
                    help="which half of the table to verify (default: all)")
    ap.add_argument("--out", type=Path, default=None,
                    help="where to write MANIFEST.csv. Defaults to "
                         "data/ARCHIVE_2026-09-15/MANIFEST.csv; point it elsewhere to verify "
                         "without writing under data/.")
    a = ap.parse_args()

    rows = []
    if a.kind in ("data", "all"):
        rows += data_rows()
    if a.kind in ("code", "all"):
        rows += code_rows()

    out = a.out if a.out is not None else ARCHIVE / "MANIFEST.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    if not a.quiet:
        fails = [r for r in rows if r["verified"] == "FAIL"]
        for kind in ("data", "code", "artifact"):
            sel_kind = [r for r in rows if r["kind"] == kind]
            if not sel_kind:
                continue
            print(f"--- {kind} ---")
            for t in (1, "1-KEEP", "LIBRARY", "CANDIDATE", 2, 3, 4):
                sel = [r for r in sel_kind if r["tier"] == t]
                if not sel:
                    continue
                gb = sum(r["bytes"] for r in sel) / 1e9
                print(f"  tier {str(t):<7} {len(sel):>3} entries  {gb:>8.3f} GB")
        print()
        for r in rows:
            if r["tier"] not in (1, "1-KEEP", "LIBRARY", "CANDIDATE", 2):
                continue
            if r["kind"] == "data" and r["tier"] == 2:
                continue          # data tier 2 is checked by carve.py, not here
            flag = {"PASS": "ok  ", "FAIL": "FAIL"}.get(r["verified"], "    ")
            print(f"  [{flag}] {r['kind']:<8} t{str(r['tier']):<6} "
                  f"{r['path'].replace(D.REPO + '/', '')}")
            print(f"            {r['verify_detail']}")
            if r["prereq"]:
                print(f"            PREREQ: {r['prereq']}")
        if fails:
            print(f"\n{len(fails)} CHECK(S) FAILED -- retire.py will refuse these.")
        print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
