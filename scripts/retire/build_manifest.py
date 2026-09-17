#!/usr/bin/env python3
"""Re-verify every declared disposition and write the retirement manifest.

Nothing is retired on the strength of the table in dispositions.py. Each entry names a CHECK; this
script runs it against the filesystem now and records pass/fail with the detail. `retire.py` acts
only on entries whose check PASSED, so a table that has drifted from disk fails closed.

    python scripts/retire/build_manifest.py            # writes data/ARCHIVE_2026-09-15/MANIFEST.csv
    python scripts/retire/build_manifest.py --quiet
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
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


CHECKS = {
    "plain": check_plain,
    "sibling_name_under": check_sibling_name_under,
    "orphan_temp": check_orphan_temp,
    "subset_dir": check_subset_dir,
    "names_exist_under": check_names_exist_under,
    "identical_to": check_identical_to,
    "donor_local_cover": check_donor_local_cover,
    "identical_keys_subset": check_identical_keys_subset,
}

FIELDS = ["tier", "path", "exists", "bytes", "n_files", "mtime_utc", "sha256", "check",
          "check_arg", "verified", "verify_detail", "superseded_by", "evidence"]


def row(tier, path, check=None, arg=None, superseded_by="", evidence=""):
    p = Path(path)
    b, n, m = measure(p)
    verified, detail = ("", "not checked (untouched tier)")
    if check:
        fn = CHECKS[check]
        ok, detail = fn(p, arg)
        verified = "PASS" if ok else "FAIL"
    return dict(tier=tier, path=str(p), exists=p.exists(), bytes=b, n_files=n, mtime_utc=m,
                sha256=sha256(p) if p.is_file() else "dir", check=check or "",
                check_arg=arg or "", verified=verified, verify_detail=detail,
                superseded_by=superseded_by, evidence=evidence)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    ARCHIVE.mkdir(parents=True, exist_ok=True)
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

    out = ARCHIVE / "MANIFEST.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    if not a.quiet:
        fails = [r for r in rows if r["verified"] == "FAIL"]
        for t in (1, "1-KEEP", 2, 3, 4):
            sel = [r for r in rows if r["tier"] == t]
            gb = sum(r["bytes"] for r in sel) / 1e9
            print(f"tier {str(t):<7} {len(sel):>3} entries  {gb:>8.2f} GB")
        print()
        for r in rows:
            if r["tier"] in (1, "1-KEEP"):
                flag = {"PASS": "ok  ", "FAIL": "FAIL"}.get(r["verified"], "    ")
                print(f"  [{flag}] {r['bytes']/1e9:8.3f} GB  {r['path'].replace(D.REPO + '/', '')}")
                print(f"            {r['verify_detail']}")
        if fails:
            print(f"\n{len(fails)} CHECK(S) FAILED -- retire.py will refuse these.")
        print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
