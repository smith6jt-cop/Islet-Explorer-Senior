"""The two-donor figure set driver.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/report/run_report.py --section pheno xen [--donors 6436 6442] [--only P01 X03]
    ../../.venv/bin/python scripts/report/run_report.py --preflight
    ../../.venv/bin/python scripts/report/run_report.py --verify

Sections: ``pheno`` (P01-P05, immutable run + Step 2 only), ``xen`` (X01-X03, bundle + base only),
``identity`` (P06-P11 from the PhenoCycler release, X04-X09 from the Xenium release), ``multi`` (M01-M06,
serial-section integration + same-section IF), ``bakeoff`` (B01-B04: every arm beside its nulls, from the
screen verdicts and the Xenium gate records), ``all``.  A figure whose release or stage output does not
exist yet renders a stated "not available" panel.  Every figure is recorded in ``figure_manifest.json``
with the sha256 of every input it read; ``--verify`` re-hashes them.
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
    ap.add_argument("--section", nargs="+", default=["pheno", "xen"], choices=["pheno", "xen", "identity", "multi", "bakeoff", "kronos2doc", "all"])
    ap.add_argument("--only", nargs="+", default=None)
    ap.add_argument("--preflight", action="store_true")
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    if "all" in a.section:
        a.section = ["pheno", "xen", "identity", "multi", "bakeoff", "kronos2doc"]
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
    def _run_section(mod):
        for d in a.donors:
            for fid, fn in mod.PER_DONOR.items():
                if a.only and fid not in a.only:
                    continue
                rec = fn(d, prov); man.add_figure(rec)
                print(f"[report] {fid} donor {d}: {rec.files[0]}", flush=True)
        for fid, fn in mod.PAIR.items():
            if a.only and fid not in a.only:
                continue
            rec = fn(a.donors, prov); man.add_figure(rec)
            print(f"[report] {fid} pair: {rec.files[0]}", flush=True)
    if "identity" in a.section:
        import pheno_identity, xen_identity
        _run_section(pheno_identity); _run_section(xen_identity)
    if "multi" in a.section:
        import multi
        _run_section(multi)
    if "kronos2doc" in a.section:
        import kronos2_doc_fig
        _run_section(kronos2_doc_fig)
    if "bakeoff" in a.section:
        import bakeoff
        for rec in bakeoff.run(prov, a.only):
            man.add_figure(rec)
    p = man.write()
    print(f"[report] manifest {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
