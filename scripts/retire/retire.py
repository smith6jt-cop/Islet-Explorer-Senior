#!/usr/bin/env python3
"""Execute the retirement: delete Tier 1, quarantine Tier 2.

Fails closed. An entry is acted on only if the manifest says its check PASSED. Tier 3 (overturned
by the still-read / narrower adversarial passes) and Tier 4 (the maintainer's call) are refused
outright, as is anything on the keep-beside list.

    python scripts/retire/retire.py --tier 1                        # dry run (data)
    python scripts/retire/retire.py --tier 1 --execute
    python scripts/retire/retire.py --tier 2 --execute              # requires carve.py to have run
    python scripts/retire/retire.py --kind code --tier 1            # dry run (code)
    python scripts/retire/retire.py --kind code --tier 2 --execute  # quarantine only

Code entries add two refusals of their own:

* everything on ``KEEP_BESIDE_CODE`` is refused outright. That list is the ``merge/hpg-sameslide``
  exclusion, encoded so a future run cannot retire a file the HiPerGator session is developing on
  even if someone puts it on a tier by mistake.
* an entry is held while anything in its ``blocked_by`` is not moving in the same invocation.
  A within-set import (broadtype -> dct/regions.py, aggregate_islets -> fixed_islet_aggregation)
  only stays safe if both ends move together.

``--mode delete`` is refused for code. Quarantine is a rename that a ``git mv`` back undoes;
deleting an entrypoint is a second, reviewable commit, and this script is not where that decision
gets made.
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

REFUSE = ({p for p, _ in D.KEEP_BESIDE_TIER1} | {p for p, _ in D.TIER3} | {p for p, _ in D.TIER4}
          | {p for p, _ in D.KEEP_BESIDE_CODE} | {p for p, _ in D.CODE_TIER3}
          | {p for p, _ in D.CODE_TIER4})


def under_protected_tree(path: str) -> str | None:
    """The keep-beside TREE a path sits under, if any. A file added inside
    core/phenocycler/scripts/sameslide/ after this list was written is still protected."""
    for tree in D.KEEP_BESIDE_CODE_TREES:
        if path == tree or path.startswith(tree.rstrip("/") + "/"):
            return tree
    return None

# Quarantine must land on the SAME filesystem as the source, or shutil.move copies 80 GB instead of
# renaming it. /data and /home are separate devices here.
def quarantine_for(p: Path) -> Path:
    s = str(p)
    if s.startswith(str(REPO)):
        return ARCHIVE / "quarantine" / p.relative_to(REPO)
    if s.startswith("/data/"):
        return Path("/data/ARCHIVE_2026-09-15_quarantine") / p.relative_to("/data")
    raise SystemExit(f"no quarantine root for {p} -- refusing rather than copying across devices")


def load(manifest: Path) -> list[dict]:
    if not manifest.exists():
        sys.exit(f"no manifest at {manifest} -- run build_manifest.py first")
    return list(csv.DictReader(manifest.open()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", type=int, choices=(1, 2), required=True)
    ap.add_argument("--kind", choices=("data", "code", "artifact"), default="data",
                    help="which half of the table to act on (default: data)")
    ap.add_argument("--manifest", type=Path, default=MANIFEST,
                    help="the manifest to read (default: the one under data/ARCHIVE_2026-09-15)")
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--mode", choices=("move", "delete"), default="move",
                    help="move (default) quarantines on the SAME filesystem, so it is a rename: "
                         "costs no space and is reversible. delete is irreversible and is the "
                         "maintainer's call. Refused for --kind code.")
    a = ap.parse_args()

    if a.kind in ("code", "artifact") and a.mode == "delete":
        sys.exit(f"refusing --mode delete for --kind {a.kind}: a retired entrypoint is quarantined "
                 "first and deleted in a second, reviewable commit")

    rows = load(a.manifest)
    sel = [r for r in rows if r["tier"] == str(a.tier) and r.get("kind", "data") == a.kind]
    if not sel:
        sys.exit(f"no tier-{a.tier} {a.kind} entries in {a.manifest}")

    if a.tier == 2 and a.kind == "data" and not CARVE.exists():
        sys.exit(f"refusing tier 2: no carve record at {CARVE}. Run carve.py --execute first.")

    act, skip = [], []
    for r in sel:
        p = Path(r["path"])
        tree = under_protected_tree(r["path"])
        if r["path"] in REFUSE:
            skip.append((r, "on the keep/overturned list"))
        elif tree:
            skip.append((r, f"under the protected tree {tree} ({D.HPG_BRANCH} is developing on it)"))
        elif not p.exists():
            skip.append((r, "already gone"))
        elif a.kind == "data" and a.tier == 2:
            act.append(r)          # data tier 2 is gated by carve.py, not by a check here
        elif r["verified"] != "PASS":
            skip.append((r, f"check did not pass: {r['verify_detail']}"))
        else:
            act.append(r)

    # A blocker that is not itself moving turns a safe group move into a broken import. A blocker
    # already off the tree counts as satisfied, so a tier-2 pass that took the importer away
    # unblocks the importee on the next tier-1 pass. Iterate: holding one entry can drop another
    # entry's blocker out of the acting set.
    def unsatisfied(r, moving):
        return [b for b in (r.get("blocked_by") or "").split()
                if b and b not in moving and Path(b).exists()]

    while True:
        moving = {r["path"] for r in act}
        held = [r for r in act if unsatisfied(r, moving)]
        if not held:
            break
        for r in held:
            missing = unsatisfied(r, moving)
            skip.append((r, "held: blocked_by not moving in this run: "
                            + ", ".join(m.replace(str(REPO) + "/", "") for m in missing)))
        act = [r for r in act if r not in held]

    verb = "DELETE (irreversible)" if a.mode == "delete" else "MOVE to quarantine (reversible)"
    total = sum(int(r["bytes"]) for r in act)
    unit, scale = ("GB", 1e9) if a.kind == "data" else ("KB", 1e3)
    print(f"tier {a.tier} / {a.kind}: {verb}\n")
    for r in act:
        print(f"  {int(r['bytes'])/scale:9.3f} {unit}  {r['path'].replace(str(REPO) + '/', '')}")
        if r.get("produces"):
            print(f"               produces (kept, only the producer moves): {r['produces']}")
    for r, why in skip:
        print(f"  SKIP           {r['path'].replace(str(REPO) + '/', '')}  -- {why}")
        if r.get("prereq"):
            print(f"               PREREQ: {r['prereq']}")
    print(f"\n  {len(act)} entries, {total/scale:.3f} {unit}; {len(skip)} skipped")

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
                         kind=a.kind, tier=a.tier, path=r["path"], bytes=r["bytes"], action=where,
                         superseded_by=r["superseded_by"], evidence=r["evidence"]))
        print(f"  {'rm' if a.mode == 'delete' else 'mv'} {p}")

    new = not LOG.exists()
    with LOG.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(done[0]))
        if new:
            w.writeheader()
        w.writerows(done)
    print(f"\n{len(done)} entries, {total/scale:.3f} {unit} -> {LOG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
