#!/usr/bin/env python3
"""Drop a SUPERSEDED.md into every directory that stacks more than one generation.

This is the part of the retirement that stops the error recurring. Citing
typing/subpanel/subpanel_calls.parquet as the current endocrine typing happened twice in one
session because nothing on disk said otherwise. Each note names the CURRENT artifact first, then
every superseded sibling with its replacement and the measurement that settles it.

Generated from scripts/retire/dispositions.py, so the notes cannot drift from the table.

    python scripts/retire/write_superseded_notes.py [--execute]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dispositions as D  # noqa: E402

REPO = Path(D.REPO)

# The current artifact in each stacked directory, stated once, in the maintainer's own numbers.
CURRENT = {
    "data/development/typing": (
        "cells/  -- 20 donors, 25,230,109 cells, 24 columns\n"
        "    Beta 226,183 | Alpha 279,686 | Delta 90,658 | is_endocrine 596,527\n"
        "    resolved_by: fallthrough 20,258,524 | definition 2,779,374 | structure 1,440,160 |\n"
        "                 hormone 596,527 | supporting 155,524\n"
        "    Built Fri 2026-09-11 -> 09-14 -> 09-15. Eight classes plus Excluded. It carries\n"
        "    NO CD99 question and no Unresolved class: the 2026-09-18 islet chain (structure,\n"
        "    exocrine-inside-an-islet excluded, CD99 fitted PER DONOR AND PER SIZE STRATUM with\n"
        "    a floor guard, nucleus_elongation as a SECOND ROUTE that can only ADD and only in\n"
        "    the >10-cell stratum, then hormone argmax over INS/GCG/SST) is the decided method\n"
        "    and NO DONOR HAS BEEN RE-RUN with it. Do not read this build as carrying it.\n"
        "    refined_class does not exist here either: scripts/typing/refine.py is written and\n"
        "    not wired, and its 5th-percentile demotion is an open question, not a decision.\n"
        "    Supporting, measured 2026-09-18: subpanel_rebuilt/ qupath_fundamentals/\n"
        "    islet_membership/ crossed/ immune/ geometry/ trajectory/ all 20 donors |\n"
        "    revised_cells/ (7, batch-2). immune/ was rewritten 2026-09-17: 20 donors,\n"
        "    2,934,898 Immune cells, a second pass that never modifies cell_class."),
    "data/development/restain": (
        "valis/donor_id=<d>/  -- the donor-local files are the resumable authority.\n"
        "    Per the lane README: 'The resumable authority is the donor directory\n"
        "    data/development/restain/valis/donor_id=<donor>/, not a batch summary that a later\n"
        "    subset can replace.' 13 donors PASS at 0.041-0.162 um median offset."),
    "core/phenocycler/data/integration": (
        "registration/  -- 20 donors, live. The .backup_* siblings are pre-change snapshots."),
    "data/development/celltype": (
        "islet_residents_all.parquet  -- the variant with live readers.\n"
        "    uniform/ is a LIVE input to the current typing lane; do not retire it."),
    "data/development/typing/subpanel": (
        "NOTHING HERE IS CURRENT for cell calls. See ../cells/.\n"
        "    subpanel_log.json IS current: it is one of the two parameter sidecars that\n"
        "    typed-orbit Phase 0 rebuilds core/phenocycler/scripts/typing/subpanel.py from."),
}

EXTRA = {
    "data/development/typing/subpanel": [
        ("subpanel_calls.parquet",
         "SUPERSEDED by ../cells/ -- Beta 191,485 / Alpha 246,055 / Delta 62,797 here against "
         "226,183 / 279,686 / 90,658 there. The row deficit is exactly the 7 batch-2 donors. "
         "KEPT ON DISK ANYWAY: typed-orbit Phase 0 names it as the row-for-row reconstruction "
         "gate, so Phase 0 cannot run once it is gone. Do not cite it; do not delete it."),
    ],
    "data/development/broadtype": [
        ("release/, arms/, frame/, scratch/",
         "The broadtype PRODUCER was WITHDRAWN on 2026-09-08 "
         "(docs/HANDOFF_2026-09-08_TYPING_REBUILD.md: 'Withdrawn. The producer choice. "
         "data/development/broadtype/release/ still exists on disk and must not be read as a "
         "verified typing.'). The directories stay because config.ini binds release/, the "
         "--verify manifest declares frame/ files as inputs and _pair/bakeoff as outputs, and "
         "five live modules use scratch/ as their cache root. Withdrawn is not the same as dead."),
    ],
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--execute", action="store_true")
    a = ap.parse_args()

    by_dir: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    for e in D.TIER1:
        p = Path(e["path"])
        by_dir[str(p.parent)].append((p.name, "retired (tier 1)", e["evidence"]))
    for e in D.TIER2:
        p = Path(e["path"])
        sb = e["superseded_by"].replace(str(REPO) + "/", "")
        by_dir[str(p.parent)].append((p.name, f"superseded by `{sb}`", e["evidence"]))

    for d, extras in EXTRA.items():
        for name, why in extras:
            by_dir[str(REPO / d)].append((name, "see below", why))
    for d in CURRENT:
        by_dir.setdefault(str(REPO / d), [])

    written = 0
    for d, entries in sorted(by_dir.items()):
        dd = Path(d)
        if not dd.is_dir():
            continue
        # Don't drop notes into directories a server reads as a data source.
        if "/data/tissuumaps/" in str(dd):
            continue
        rel = str(dd.relative_to(REPO)) if str(dd).startswith(str(REPO)) else str(dd)
        cur = CURRENT.get(rel)
        lines = [f"# Which artifact here is current",
                 "",
                 f"Generated {date.today().isoformat()} by `scripts/retire/write_superseded_notes.py`",
                 "from `scripts/retire/dispositions.py`. Regenerate rather than hand-edit.", ""]
        if cur:
            lines += ["## Current", "", f"    {cur}", ""]
        if entries:
            lines += ["## Superseded or retired", ""]
            for name, disp, why in sorted(entries):
                lines += [f"- **`{name}`** — {disp}", f"  {why}", ""]
        lines += ["---", "",
                  "If you are about to cite a number out of this directory, check the mtime of the",
                  "file you are reading against its siblings first. Several generations of the same",
                  "measurement sit here and only one of them is current.", ""]
        out = dd / "SUPERSEDED.md"
        text = "\n".join(lines)
        if a.execute:
            out.write_text(text)
        written += 1
        print(f"{'wrote' if a.execute else 'would write'} {out.relative_to(REPO) if str(out).startswith(str(REPO)) else out}"
              f"  ({len(entries)} entr{'y' if len(entries) == 1 else 'ies'}{', current named' if cur else ''})")
    print(f"\n{written} note(s){'' if a.execute else '; dry run, pass --execute'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
