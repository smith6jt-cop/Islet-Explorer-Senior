#!/usr/bin/env python3
"""Unit tests for export_viewer_cell_markers.py.

No pandas, no Parquet, no cohort data — the legacy-CSV path, the skip
behaviour, coordinate hygiene and the coverage record are all exercised on
stdlib alone, which is also what makes them runnable anywhere.

    python scripts/check_export_viewer_cell_markers.py
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "senior"))
import export_viewer_cell_markers as m  # noqa: E402

SCRIPT = Path(__file__).resolve().parent / "senior" / "export_viewer_cell_markers.py"
fails = []


def check(label, got, want):
    ok = got == want
    print(f"{'ok  ' if ok else 'FAIL'}  {label}: {got!r}" + ("" if ok else f" != {want!r}"))
    if not ok:
        fails.append(label)


_tmp = tempfile.TemporaryDirectory()
DATA = Path(_tmp.name)
donors = DATA / "app_data" / "donors"
donors.mkdir(parents=True)
(DATA / "parquet" / "tissue").mkdir(parents=True)

# A donor with legacy cell data, including two rows the viewer must not be given:
# an empty coordinate and a literal NaN.
rows = [
    "X_centroid,Y_centroid,phenotype,cell_region,islet_name",
    "100.123456,200.987654,Beta cell,core,Islet_1",
    "300.5,400.5,Alpha cell,peri,Islet_1",
    ",500.0,Delta cell,tissue,",
    "600.0,NaN,Macrophage,tissue,",
    "700.25,800.75,Endothelial,tissue,",
]
(donors / "6539.csv").write_text("\n".join(rows) + "\n")

# --- discovery ------------------------------------------------------------
check("discovers donors with legacy data",
      m.discover_cases(DATA / "parquet" / "tissue", donors), ["6539"])

# --- export ---------------------------------------------------------------
out = DATA / "tissuumaps_inputs" / "markers"
proc = subprocess.run([sys.executable, str(SCRIPT), "--data-dir", str(DATA)],
                      capture_output=True, text=True)
check("exits 0", proc.returncode, 0)

written = list(csv.DictReader(open(out / "6539.csv")))
check("header is the builder's expected columns",
      list(written[0].keys()) if written else [], m.MARKER_COLS)
check("non-finite coordinates dropped", len(written), 3)
check("coordinates rounded", written[0]["X_centroid"], "100.12")
check("phenotype preserved", written[0]["phenotype"], "Beta cell")
check("cell_region carried through", written[1]["cell_region"], "peri")

cov = json.loads((out / "coverage.json").read_text())
check("coverage records the source", cov["6539"]["source"], "legacy_csv")
check("coverage counts kept cells", cov["6539"]["n_cells"], 3)
check("coverage counts dropped cells", cov["6539"]["n_dropped"], 2)
check("coverage records extent for the builder's scale check",
      (cov["6539"]["max_x_um"], cov["6539"]["max_y_um"]), (700.25, 800.75))

# --- caching --------------------------------------------------------------
proc = subprocess.run([sys.executable, str(SCRIPT), "--data-dir", str(DATA)],
                      capture_output=True, text=True)
check("second run reuses the CSV", "cached" in proc.stdout, True)
cov = json.loads((out / "coverage.json").read_text())
check("cached run still records coverage", cov["6539"]["status"], "cached")

# --- a donor with a slide but no cells ------------------------------------
proc = subprocess.run([sys.executable, str(SCRIPT), "--data-dir", str(DATA),
                       "--cases", "6450"], capture_output=True, text=True)
check("missing donor is not an error", proc.returncode, 0)
check("missing donor is named, not silent", "6450" in proc.stdout and "SKIPPED" in proc.stdout, True)
cov = json.loads((out / "coverage.json").read_text())
check("missing donor recorded as a data fact", cov["6450"]["status"], "no_cell_data")
check("no CSV written for it", (out / "6450.csv").exists(), False)

# --- nothing at all -------------------------------------------------------
empty = Path(tempfile.mkdtemp())
proc = subprocess.run([sys.executable, str(SCRIPT), "--data-dir", str(empty)],
                      capture_output=True, text=True)
check("empty data dir fails loudly", proc.returncode, 1)

print()
print(f"{len(fails)} failed" if fails else "all checks passed")
sys.exit(1 if fails else 0)
