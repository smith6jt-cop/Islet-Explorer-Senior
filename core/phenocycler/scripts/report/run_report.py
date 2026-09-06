"""The two-donor figure set driver.

    cd core/phenocycler && ../../.venv/bin/python scripts/report/run_report.py --section pheno xen [--donors 6436 6442] [--only P01 X03]
    ../../.venv/bin/python scripts/report/run_report.py --preflight
    ../../.venv/bin/python scripts/report/run_report.py --verify

Sections: ``pheno`` (P01-P05, immutable run + Step 2 only), ``xen`` (X01-X03, bundle + base only).
Identity and multi-platform sections are added once the releases exist.  Every figure is recorded in
``figure_manifest.json`` with the sha256 of every input it read; ``--verify`` re-hashes them.
"""
from __future__ import annotations

import argparse
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402
from manifest import Manifest, verify  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=list(pins.DONORS))
    ap.add_argument("--section", nargs="+", default=["pheno", "xen"], choices=["pheno", "xen"])
    ap.add_argument("--only", nargs="+", default=None)
    ap.add_argument("--preflight", action="store_true")
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    if a.verify:
        ok, rows = verify()
        return 0 if ok else 1
    pins.assert_pins()
    if a.preflight:
        print("[report] pins ok:", pins.RUN_ID, pins.EXPRESSION_CONTENT_ID[:12])
        return 0
    prov = pins.provenance(a.donors)
    man = Manifest()
    if "pheno" in a.section:
        import pheno_qc
        for d in a.donors:
            for rec in pheno_qc.run(d, prov, a.only):
                man.add_figure(rec)
    if "xen" in a.section:
        import xen_qc
        tables = {}
        for d in a.donors:
            import xen_io
            tables[d] = xen_io.load_tables(d)
            for fid, fn in xen_qc.FIGURES.items():
                if a.only and fid not in a.only:
                    continue
                rec = fn(d, prov, tables[d]); man.add_figure(rec)
                print(f"[report] {fid} donor {d}: {rec.files[0]}", flush=True)
        if not a.only or "X03" in a.only:
            rec = xen_qc.x03_negative_controls(a.donors, prov, tables); man.add_figure(rec)
            print(f"[report] X03 pair: {rec.files[0]}", flush=True)
    p = man.write()
    print(f"[report] manifest {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
