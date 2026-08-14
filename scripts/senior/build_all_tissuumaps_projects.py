#!/usr/bin/env python3
"""Build a TissUUmaps project for every donor in the cohort.

A thin, resumable loop over ``build_tissuumaps_project.py``. It adds no
extraction logic of its own — per-channel caching, stale-settings detection and
the region/marker falsification all live in the builder — but it decides WHICH
donors to build, keeps a log per donor, refuses to fill the disk, and prints the
coverage table that says what each donor actually got.

    python scripts/senior/build_all_tissuumaps_projects.py --dry-run
    python scripts/senior/build_all_tissuumaps_projects.py
    python scripts/senior/build_all_tissuumaps_projects.py --cases 6539 6450

Donors come from the ``MAP`` in ``convert_viewer_ometiff.sh``, which
``build_acquisition_batch_table.py`` names as the authoritative donor->filename
map, minus the central exclusions in ``phenocycler.cohort``. Neither list is
restated here: a second copy of a cohort definition is a second thing to get
wrong, and ``scripts/CLAUDE.md`` says scripts must not recreate or override the
exclusions.

Batch-1 donors are built first. They are the ones with islet outlines and cell
data, so the app becomes fully useful earliest; the batch-2 slides that have
channels only come last.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BUILDER = REPO / "scripts" / "senior" / "build_tissuumaps_project.py"
CONVERTER = REPO / "scripts" / "senior" / "convert_viewer_ometiff.sh"


def log(msg: str) -> None:
    print(f"[all] {msg}", flush=True)


def read_cohort_map() -> dict[str, str]:
    """case_id -> qptiff filename, parsed from the converter's bash MAP.

    Parsed rather than duplicated. The alternative is a second copy of the donor
    list that silently drifts from the one the slides were actually converted
    with.
    """
    text = CONVERTER.read_text()
    block = re.search(r"declare -A MAP=\((.*?)\n\)", text, re.S)
    if not block:
        raise RuntimeError(f"could not find the MAP array in {CONVERTER}")
    return dict(re.findall(r"\[(\w+)\]=\"([^\"]+)\"", block.group(1)))


def read_batch1(cohort: dict[str, str]) -> set[str]:
    """The batch-1 case ids, from the converter's ALL15 array."""
    m = re.search(r"ALL15=\(([^)]*)\)", CONVERTER.read_text())
    return set(m.group(1).split()) if m else set(cohort)


def read_exclusions() -> tuple[set[str], str]:
    """The central donor exclusions, imported — never restated.

    Returns the set and where it came from. If the submodule is not installed we
    say so and return nothing rather than guessing: silently building an excluded
    donor is exactly the outcome the central list exists to prevent.
    """
    sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
    try:
        from phenocycler.cohort import DONOR_EXCLUSIONS  # noqa: PLC0415

        return {str(d) for d in DONOR_EXCLUSIONS}, "phenocycler.cohort"
    except Exception as exc:  # noqa: BLE001 - submodule may be uninitialised
        return set(), f"UNAVAILABLE ({type(exc).__name__})"


def free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 1e9


def project_summary(out_dir: Path, case_id: str) -> dict:
    """What a built project actually contains, read back from its .tmap."""
    tmap = out_dir / f"{case_id}.tmap"
    if not tmap.is_file():
        return {"channels": 0, "regions": {}, "markers": 0}
    try:
        p = json.loads(tmap.read_text())
    except ValueError:
        return {"channels": 0, "regions": {}, "markers": 0}
    regions = {}
    for r in p.get("regionFiles", []):
        f = out_dir / r["path"]
        try:
            regions[r["title"].split(" (")[0]] = len(
                json.loads(f.read_text()).get("features", []))
        except (OSError, ValueError):
            regions[r["title"].split(" (")[0]] = -1
    markers = 0
    for mf in p.get("markerFiles", []):
        f = out_dir / mf["path"]
        try:
            with open(f) as fh:
                markers = sum(1 for _ in fh) - 1
        except OSError:
            markers = -1
    return {"channels": len(p.get("layers", [])), "regions": regions, "markers": markers}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--images", type=Path, default=Path("/data/islet_ome_tiff"),
                    help="Directory of <case>.ome.tiff slides")
    ap.add_argument("--out-dir", type=Path, default=Path("/data/tissuumaps"),
                    help="TissUUmaps project root")
    ap.add_argument("--regions", type=Path, default=REPO / "data" / "Islets",
                    help="Directory of QuPath islet exports")
    ap.add_argument("--structures", type=Path,
                    default=REPO / "data" / "app_data" / "json" / "structures",
                    help="Directory of <case>_structures.geojson (build_viewer_assets.py "
                         "--structures)")
    ap.add_argument("--markers", type=Path,
                    default=REPO / "data" / "tissuumaps_inputs" / "markers",
                    help="Directory of <case>.csv (export_viewer_cell_markers.py)")
    ap.add_argument("--log-dir", type=Path, default=REPO / "logs" / "tissuumaps")
    ap.add_argument("--cases", nargs="*", default=None,
                    help="Restrict to these case ids (exclusions still apply)")
    ap.add_argument("--min-free-gb", type=float, default=300.0,
                    help="Skip a donor when free space is below this (default 300, matching "
                         "convert_viewer_ometiff.sh)")
    ap.add_argument("--channels", default=None,
                    help="Pass through to the builder: extract only these markers. Use for a "
                         "fast first pass over the whole cohort; a later full run backfills "
                         "the rest, since channel caching makes the second pass additive")
    ap.add_argument("--dry-run", action="store_true", help="Print the plan and exit")
    ap.add_argument("--report", action="store_true",
                    help="Print the coverage table for what is already built, and exit")
    args = ap.parse_args()

    cohort = read_cohort_map()
    batch1 = read_batch1(cohort)
    excluded, source = read_exclusions()
    if not excluded:
        log(f"ERROR: could not import the central donor exclusions ({source}). "
            f"Refusing to guess — run `git submodule update --init Phenocycler_Analysis`, "
            f"or pass --cases explicitly if you know what you are doing.")
        if not args.cases:
            return 1

    # Batch-1 first: those donors have outlines and cells, so the app becomes
    # fully useful earliest.
    ordered = sorted(cohort, key=lambda c: (c not in batch1, c))
    cases = [c for c in ordered if c not in excluded]
    skipped_excluded = [c for c in ordered if c in excluded]
    if args.cases:
        wanted = set(args.cases)
        unknown = wanted - set(cohort)
        if unknown:
            log(f"WARNING: not in the cohort map, ignoring: {', '.join(sorted(unknown))}")
        blocked = wanted & excluded
        if blocked:
            log(f"WARNING: centrally excluded, not building: {', '.join(sorted(blocked))}")
        cases = [c for c in cases if c in wanted]

    out_dir = args.out_dir.resolve()
    log(f"cohort {len(cohort)}, excluded {len(skipped_excluded)} "
        f"({', '.join(skipped_excluded) or 'none'}) via {source}, building {len(cases)}")

    if args.report:
        return report(out_dir, cases, args)

    missing = [c for c in cases if not (args.images / f"{c}.ome.tiff").is_file()]
    if missing:
        log(f"WARNING: no slide for {len(missing)} donor(s): {', '.join(missing)}")
    cases = [c for c in cases if c not in missing]

    log(f"plan: {len(cases)} donor(s) -> {out_dir}")
    for c in cases:
        gb = (args.images / f"{c}.ome.tiff").stat().st_size / 1e9
        log(f"   {c:<8} {gb:>6.1f} GB slide{'  [batch-1]' if c in batch1 else '  [batch-2]'}")
    log(f"free space at {out_dir}: {free_gb(out_dir if out_dir.exists() else Path('/')):.0f} GB")
    if args.dry_run:
        log("DRY RUN — nothing built")
        return 0

    args.log_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict] = {}
    for n, case_id in enumerate(cases, 1):
        avail = free_gb(out_dir)
        if avail < args.min_free_gb:
            log(f"[{n}/{len(cases)}] {case_id}: SKIPPED — {avail:.0f} GB free is below "
                f"--min-free-gb {args.min_free_gb:.0f}")
            results[case_id] = {"status": "disk", "free_gb": round(avail, 1)}
            continue

        cmd = [sys.executable, str(BUILDER),
               "--image", str(args.images / f"{case_id}.ome.tiff"),
               "--out-dir", str(out_dir)]
        # Attach whatever exists. A donor without outlines or cells still gets a
        # perfectly good channel viewer, and the coverage table records the gap
        # rather than the run failing on it.
        if args.regions.exists():
            cmd += ["--regions", str(args.regions)]
        if args.structures.exists():
            cmd += ["--structures", str(args.structures)]
        marker_csv = args.markers / f"{case_id}.csv"
        if marker_csv.is_file():
            cmd += ["--markers", str(marker_csv)]
        if args.channels:
            cmd += ["--channels", args.channels]

        log_path = args.log_dir / f"{case_id}.log"
        t0 = time.time()
        log(f"[{n}/{len(cases)}] {case_id}: building -> {log_path}")
        with open(log_path, "w") as fh:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
        secs = time.time() - t0
        # Exit 3 is the builder's "project written, but a requested overlay was
        # refused" — worth surfacing distinctly from a failure to build at all.
        status = {0: "ok", 3: "ok_no_markers"}.get(proc.returncode, "failed")
        results[case_id] = {"status": status, "returncode": proc.returncode,
                            "seconds": round(secs)}
        log(f"[{n}/{len(cases)}] {case_id}: {status} in {secs / 60:.1f} min")
        if status == "failed":
            tail = log_path.read_text().splitlines()[-3:]
            for line in tail:
                log(f"      {line}")

    (args.log_dir / "build_results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    return report(out_dir, cases, args, results)


def report(out_dir: Path, cases: list[str], args, results: dict | None = None) -> int:
    """The coverage table — the artifact that says what actually got built.

    Silence about a missing overlay is the failure mode worth engineering
    against: a donor that quietly has no cell layer looks identical to one that
    was never meant to have one.
    """
    log("")
    log(f"{'case':<8}{'channels':>9}{'islets':>8}{'structs':>9}{'cells':>10}   status")
    gaps = []
    for case_id in cases:
        s = project_summary(out_dir, case_id)
        islets = s["regions"].get("Islets", 0)
        structs = s["regions"].get("Structures", 0)
        status = (results or {}).get(case_id, {}).get("status", "—")
        log(f"{case_id:<8}{s['channels']:>9}{islets:>8}{structs:>9}{s['markers']:>10}   {status}")
        if s["channels"] == 0:
            gaps.append(f"{case_id}: no project")
        elif not islets:
            gaps.append(f"{case_id}: no islet outlines")
    if gaps:
        log("")
        log(f"{len(gaps)} gap(s):")
        for g in gaps:
            log(f"   {g}")
    failed = [c for c, r in (results or {}).items() if r["status"] == "failed"]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
