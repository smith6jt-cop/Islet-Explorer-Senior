#!/usr/bin/env python
"""Postconditions for the analysis tree. Run after a full notebook re-run.

Each check here exists because the corresponding mistake was actually made and
was invisible at the time — a passing run, a plausible-looking figure, a table
that simply did not move. None of them raise on their own.

    python scripts/verify_analysis_state.py [--since 'YYYY-MM-DD HH:MM']

Exits non-zero if any check fails.
"""
import argparse
import datetime as dt
import glob
import json
import os
import re
import sys
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT))
from paths import OUTPUT_ROOT  # noqa: E402

NOTEBOOKS = sorted(glob.glob(str(CODE_ROOT / "analysis" / "H*.ipynb")))

# data_utils names whose absence from an explicit import list is a NameError
# waiting to happen. Notebooks using `from data_utils import *` are exempt.
DU_NAMES = ["load_region_layer", "assign_regions_by_image", "signed_distance_by_image",
            "signed_distance_to_regions", "get_vessels", "compute_density",
            "REGION_ORDER", "load_all_data", "load_single_cells",
            "assign_vessels_to_follicles", "get_regions"]

failures = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        failures.append(name)


def code_of(path):
    nb = json.load(open(path))
    return "\n".join("".join(c["source"]) for c in nb["cells"]
                     if c["cell_type"] == "code")


def uncommented(code):
    return [l for l in code.split("\n") if not l.strip().startswith("#")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--since", default=None,
                    help="figures/tables older than this are stale "
                         "(default: skip the freshness checks)")
    ap.add_argument("--static-only", action="store_true",
                    help="check notebook source only; skip generated outputs")
    args = ap.parse_args()

    print("\n== region assignment ==")
    # The Parent hierarchy strands 1,014 vessels at the image root, 21.5% of
    # them in one donor. Nothing may resolve a region through it.
    # Only `.str.extract` on the Parent string is region assignment. H12 and H13
    # also `.str.replace` the same pattern, but that harmonises CODEX class
    # names inside Parent and is not a region lookup — flagging it would train
    # the reader to ignore this check.
    for pat, label in [(r"Parent.*str\.extract", "Parent regex for region assignment"),
                       (r"assign_region_by_distance\(", "assign_region_by_distance call"),
                       (r"['\"]geojson['\"]", "reads superseded analysis/geojson/"),
                       (r"_PIXEL_SIZES|_DEFAULT_PIXEL_SIZE", "local pixel-size table")]:
        bad = []
        for f in NOTEBOOKS:
            for l in uncommented(code_of(f)):
                if re.search(pat, l) and "geojson_regions" not in l:
                    bad.append(os.path.basename(f))
                    break
        check(f"no notebook uses {label}", not bad, ", ".join(sorted(set(bad))))

    print("\n== cropping ==")
    # Cropping to a CODEX-sized window discarded 75-88% of the data to address
    # only the unequal-effort mechanism. Thresholds must stay unreachable.
    bad = []
    for f in NOTEBOOKS:
        for m in re.finditer(r"CROP_THRESHOLD[_A-Z]*\s*(?:,\s*CROP_THRESHOLD_Y\s*)?=\s*([^\n#]+)",
                             code_of(f)):
            if "inf" not in m.group(1):
                bad.append(f"{os.path.basename(f)}: {m.group(1).strip()}")
    check("no notebook can crop a Phenocycler section", not bad, "; ".join(bad))

    print("\n== imports ==")
    # A missing name in an explicit import list fails only when the cell runs,
    # which for the single-cell notebooks is after an 11.8 GB read.
    bad = []
    for f in NOTEBOOKS:
        code = code_of(f)
        if "from data_utils import *" in code:
            continue
        m = re.search(r"from data_utils import \(([^)]*)\)", code)
        imported = set(re.findall(r"\w+", m.group(1))) if m else set()
        used = {n for n in DU_NAMES if re.search(rf"(?<![\w.]){n}\s*[(\[]", code)}
        if used - imported:
            bad.append(f"{os.path.basename(f)}: {sorted(used - imported)}")
    check("every data_utils name used is imported", not bad, "; ".join(bad))

    print("\n== numpy 2.0 ==")
    bad = [os.path.basename(f) for f in NOTEBOOKS
           if any(re.search(r"[\]\)]\.ptp\(\)", l) for l in uncommented(code_of(f)))]
    check("no removed ndarray.ptp()", not bad, ", ".join(bad))

    if args.static_only:
        print("\n== outputs ==\n  SKIP  generated outputs (--static-only)")
    else:
        print("\n== outputs ==")
        figs = glob.glob(str(OUTPUT_ROOT / "analysis" / "figures" / "*.png"))
        check("figures exist", len(figs) > 0, f"{len(figs)} PNGs")
    if not args.static_only and args.since:
        cutoff = dt.datetime.strptime(args.since, "%Y-%m-%d %H:%M").timestamp()
        # A figure older than the run is one whose producing cell was removed;
        # it sits there looking current, and analysis/figures is gitignored.
        stale = [os.path.basename(f) for f in figs if os.path.getmtime(f) < cutoff]
        check("no figure predates the run", not stale,
              f"{len(stale)} stale: {', '.join(sorted(stale)[:5])}")
        # Only tables a notebook still asks for are judged against the notebook
        # cutoff. The region layer's own tables come from scripts that run
        # before the notebooks and legitimately carry an earlier timestamp;
        # holding them to this cutoff reports seven false positives and buries
        # the real ones. A notebook-produced table older than the run is an
        # orphan — its save_table() call is gone, and analysis/tables is
        # gitignored, so nothing else would ever tell you.
        wanted = set()
        for f in NOTEBOOKS:
            wanted |= set(re.findall(r'save_table\([^,]+,\s*["\']([^"\']+)["\']',
                                     code_of(f)))
        tables = glob.glob(str(OUTPUT_ROOT / "analysis" / "tables" / "*.csv"))
        stale_t = [os.path.basename(t) for t in tables
                   if os.path.getmtime(t) < cutoff
                   and os.path.basename(t)[:-4] in wanted]
        check("no notebook table predates the run", not stale_t,
              f"{len(stale_t)} stale: {', '.join(sorted(stale_t)[:5])}")

        orphans = [os.path.basename(t) for t in tables
                   if os.path.getmtime(t) < cutoff
                   and os.path.basename(t)[:-4] not in wanted
                   and re.match(r"H\d+_", os.path.basename(t))]
        check("no orphaned H* table", not orphans,
              f"{len(orphans)}: {', '.join(sorted(orphans)[:5])}")

    print("\n== region layer ==")
    if args.static_only:
        print("  SKIP  generated region layer (--static-only)")
    else:
        import pandas as pd
        qc = OUTPUT_ROOT / "analysis" / "tables" / "region_layer_qc.csv"
        if qc.exists():
            q = pd.read_csv(qc)
            check("14 donors in the QC table", len(q) == 14, f"{len(q)} rows")
            check("Trabeculae present for every donor",
                  (q["trabeculae_source"] == "qupath").all(),
                  f"{(q['trabeculae_source'] != 'qupath').sum()} missing")
            flagged = q.loc[q["flagged"].notna(), "Sample"].tolist()
            check("only HDL073 flagged", flagged == ["HDL073"], f"flagged: {flagged}")
        else:
            check("region_layer_qc.csv exists", False)

    print()
    if failures:
        print(f"FAILED {len(failures)}: {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
