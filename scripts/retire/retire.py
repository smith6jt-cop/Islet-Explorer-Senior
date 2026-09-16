#!/usr/bin/env python3
"""Execute the retirement: delete Tier 1, quarantine Tier 2.

Fails closed. An entry is acted on only if the manifest says its check PASSED. Tier 3 (overturned
by the still-read / narrower adversarial passes) and Tier 4 (the maintainer's call) are refused
outright, as is anything on the keep-beside list.

    python scripts/retire/retire.py --tier 1            # dry run
    python scripts/retire/retire.py --tier 1 --execute
    python scripts/retire/retire.py --tier 2 --execute  # requires carve.py to have run
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dispositions as D  # noqa: E402

REPO = Path(D.REPO)
ARCHIVE = REPO / "data" / "ARCHIVE_2026-09-15"
MANIFEST = ARCHIVE / "MANIFEST.csv"
CARVE = ARCHIVE / "carve" / "CARVE.json"
LOG = ARCHIVE / "RETIRED.csv"

REFUSE = {p for p, _ in D.KEEP_BESIDE_TIER1} | {p for p, _ in D.TIER3} | {p for p, _ in D.TIER4}

# Quarantine must land on the SAME filesystem as the source, or shutil.move copies 80 GB instead of
# renaming it. /data and /home are separate devices here.
def quarantine_for(p: Path) -> Path:
    s = str(p)
    if s.startswith(str(REPO)):
        return ARCHIVE / "quarantine" / p.relative_to(REPO)
    if s.startswith("/data/"):
        return Path("/data/ARCHIVE_2026-09-15_quarantine") / p.relative_to("/data")
    raise SystemExit(f"no quarantine root for {p} -- refusing rather than copying across devices")


def load() -> list[dict]:
    if not MANIFEST.exists():
        sys.exit(f"no manifest at {MANIFEST} -- run build_manifest.py first")
    return list(csv.DictReader(MANIFEST.open()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", type=int, choices=(1, 2), required=True)
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--mode", choices=("move", "delete"), default="move",
                    help="move (default) quarantines on the SAME filesystem, so it is a rename: "
                         "costs no space and is reversible. delete is irreversible and is the "
                         "maintainer's call.")
    a = ap.parse_args()

    rows = load()
    sel = [r for r in rows if r["tier"] == str(a.tier)]
    if not sel:
        sys.exit(f"no tier-{a.tier} entries in the manifest")

    if a.tier == 2 and not CARVE.exists():
        sys.exit(f"refusing tier 2: no carve record at {CARVE}. Run carve.py --execute first.")

    act, skip = [], []
    for r in sel:
        p = Path(r["path"])
        if r["path"] in REFUSE:
            skip.append((r, "on the keep/overturned list"))
        elif not p.exists():
            skip.append((r, "already gone"))
        elif a.tier == 1 and r["verified"] != "PASS":
            skip.append((r, f"check did not pass: {r['verify_detail']}"))
        else:
            act.append(r)

    verb = "DELETE (irreversible)" if a.mode == "delete" else "MOVE to quarantine (reversible)"
    total = sum(int(r["bytes"]) for r in act)
    print(f"tier {a.tier}: {verb}\n")
    for r in act:
        print(f"  {int(r['bytes'])/1e9:8.3f} GB  {r['path'].replace(str(REPO) + '/', '')}")
    for r, why in skip:
        print(f"  SKIP           {r['path'].replace(str(REPO) + '/', '')}  -- {why}")
    print(f"\n  {len(act)} entries, {total/1e9:.2f} GB")

    if not a.execute:
        print("\ndry run; pass --execute")
        return 0

    done = []
    for r in act:
        p = Path(r["path"])
        if a.mode == "delete":
            shutil.rmtree(p) if p.is_dir() else p.unlink()
            where = "deleted"
        else:
            dest = quarantine_for(p)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), str(dest))  # same filesystem -> rename, no copy
            where = str(dest)
        done.append(dict(retired_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                         tier=a.tier, path=r["path"], bytes=r["bytes"], action=where,
                         superseded_by=r["superseded_by"], evidence=r["evidence"]))
        print(f"  {'rm' if a.mode == 'delete' else 'mv'} {p}")

    new = not LOG.exists()
    with LOG.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(done[0]))
        if new:
            w.writeheader()
        w.writerows(done)
    print(f"\n{len(done)} entries, {total/1e9:.2f} GB -> {LOG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
