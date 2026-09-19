#!/usr/bin/env python
"""SUPERSEDED 2026-09-18. The current pLN Xenium typer is
``analysis/xenium/scripts/pln_immune_typing_v3.py``, which writes the release
``data/development/xenium_typing/release_pln_immune_v3``; that release declares it supersedes
``release_pln_immune_v2``, and ``scripts/senior/xenium_release.py`` honours the declaration, so
no builder resolves anything this file produced. Recorded as a verified orphan in
``scripts/retire/dispositions.py`` (code tier 1). Kept as the record of the ROI-rebinding and
the two guards below; do not read it as the current method.

Run the pancreas Xenium typing lane on the pLN sections (ruling 2026-09-17).

The maintainer's ruling: *all immune cell calls from pancreas can be applied to pLN*.  The
immune lineage vocabulary and the typing METHOD carry over; the labels themselves do not, and
cannot -- panc and pLN are separately segmented bundles from the same instrument run, so their
``cell_id`` namespaces are disjoint (measured on 6442: 44 of 177,799 pLN ids collide with panc
ids, and the collisions are coincidental, not a mapping).  So the typer is RUN on the pLN cells.

This driver owns nothing of the lane's method.  It reuses, unmodified:

  core/phenocycler/scripts/xenium_celltype/build_phenotyped_base.py   (bases; --roi already exists)
  analysis/xenium/scripts/finalize_panel.py                           (identity/subpanel masks)
  analysis/xenium/scripts/hierarchical_phenotype.py                   (the two-tier typer)
  core/phenocycler/scripts/xenium_typing/{gate,nullgate_typer,release}.py

What it adds is the ROI rebinding those modules do not have, and two guards.

ROI REBINDING.  ``pins.py`` pins ``ROI = "panc"`` and binds it into the *default arguments* of
``paths_table``/``source_path``/``serial``/..., so those defaults were frozen at import and
setting ``pins.ROI`` alone would not move them.  ``patch_pins()`` rebinds every ROI-dependent
accessor to ``pln_1`` and redirects the lane's arm/gate output to a ``pln/`` subtree so the
existing pancreas evidence is never touched.

GUARD 1 -- no pancreas IF on a pLN cell.  Every one of the 11 cohort donors with a local pLN
bundle carries pLN and panc on the SAME slide serial (measured: 112, 115, 6374, 6442, 6461,
6482, 6521, 6533, 6551, 6566, 6593 -- all 11).  ``pins.if_values_path`` keys the re-stain array
on ``<donor>_<serial>_<compartment>_<marker>.npy``, i.e. on donor and serial only, so under an
unpatched ROI it resolves a pLN request to the PANCREAS section's per-cell IF array (6442 has
one on disk).  That is a protein measurement of a different block.  ``patch_pins`` routes it to
a path that cannot exist, so ``release.build_donor`` records IF as unavailable -- which it is.

GUARD 2 -- the h5ad write defect.  ``build_phenotyped_base.build_one`` stamps
``uns['prep']['gdna']['area_table']`` as a list of dicts, which anndata cannot write; the build
raises ``TypeError: Can't implicitly convert non-string objects to strings`` before writing a
byte.  This is ROI-independent and reproduces in six lines (see the write-up).  Rather than edit
a file this lane does not own, ``_json_safe_write`` JSON-encodes that one provenance field at
write time.  Nothing downstream reads it.

    .venv/bin/python analysis/xenium/scripts/pln_typing.py bases --donors 6533
    .venv/bin/python analysis/xenium/scripts/pln_typing.py panel --donors 6533
    <scvi-env>/bin/python analysis/xenium/scripts/pln_typing.py type  --donors 6533
    .venv/bin/python analysis/xenium/scripts/pln_typing.py calls --donors 6533 --judge
    .venv/bin/python analysis/xenium/scripts/pln_typing.py release --donors 6533
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[3]
XT = REPO / "core" / "phenocycler" / "scripts" / "xenium_typing"
XC = REPO / "core" / "phenocycler" / "scripts" / "xenium_celltype"
for _p in (str(XT), str(XC)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

ROI = "pln_1"

#: The eleven donors that are BOTH in the lane's acquired cohort (pins.donors()) and have a
#: pLN bundle on this workstation.  The other eleven cohort donors declare a pln_1 bundle that
#: has not been transferred here; 6591 declares none at all.  Not typed is not zero.
DEFAULT_DONORS = ("112", "115", "6374", "6442", "6461", "6482",
                  "6521", "6533", "6551", "6566", "6593")


def patch_pins():
    """Rebind every ROI-dependent accessor in ``pins`` to pLN, and block the panc IF arrays."""
    import pins

    pins.ROI = ROI
    _orig_paths_table = pins.paths_table

    def paths_table(csv=pins.XENIUM_PATHS, roi=ROI):
        return _orig_paths_table(csv, roi)

    def source_path(donor, csv=pins.XENIUM_PATHS, roi=ROI):
        t = paths_table(csv, roi)
        hit = t.loc[t["donor_id"] == str(donor), "source_path"]
        if len(hit) != 1:
            raise pins.PinError(f"donor {donor} roi={roi}: {len(hit)} bundle rows in {csv}")
        return str(hit.iloc[0])

    def serial(donor):
        return pins.serial_of(source_path(donor))

    def xenium_run(donor):
        return pins.xenium_run_of(source_path(donor))

    def bundle_dir(donor):
        return pins.XENIUM_ROOT / source_path(donor)

    def if_values_path(donor, marker, compartment="cell"):
        # Deliberately unresolvable: the re-stain arrays are pancreas-section protein and the
        # filename carries no ROI, so the real accessor would hand a pLN release the pancreas
        # block's values.  There is no pLN IF; unavailable is the honest answer.
        return (pins.RESTAIN_PROTEIN_RNA
                / f"__no_pln_if__/{donor}_{serial(donor)}_{compartment}_{marker}.npy")

    pins.paths_table = paths_table
    pins.source_path = source_path
    pins.serial = serial
    pins.xenium_run = xenium_run
    pins.bundle_dir = bundle_dir
    pins.serial_map = lambda csv=pins.XENIUM_PATHS, roi=ROI: {
        d: pins.serial_of(p) for d, p in zip(paths_table(csv, roi)["donor_id"],
                                             paths_table(csv, roi)["source_path"])}
    pins.run_map = lambda csv=pins.XENIUM_PATHS, roi=ROI: {
        d: pins.xenium_run_of(p) for d, p in zip(paths_table(csv, roi)["donor_id"],
                                                 paths_table(csv, roi)["source_path"])}
    pins.analysis_dir = lambda donor: bundle_dir(donor) / "analysis"
    pins.cells_parquet = lambda donor: bundle_dir(donor) / "cells.parquet"
    pins.base_path = lambda donor: (pins.PROCESSED_ROOT / serial(donor) / ROI
                                    / f"{serial(donor)}_phenotyped.h5ad")
    pins.sidecar_path = lambda donor: (pins.PROCESSED_ROOT / serial(donor) / ROI
                                       / f"{serial(donor)}_hierarchical.h5ad")
    pins.if_values_path = if_values_path
    # Lane evidence for pLN lives in its own subtree; the pancreas arms/gate are untouched.
    pins.arm_dir = lambda arm, donor: pins.OUT / "pln" / "arms" / arm / f"donor_id={donor}"
    pins.gate_dir = lambda arm: pins.OUT / "pln" / "gate" / arm
    pins.gate_json = lambda arm, donor: pins.gate_dir(arm) / f"donor_id={donor}.json"
    return pins


# --------------------------------------------------------------------------- stage: bases
def _json_safe_write():
    """See GUARD 2.  Coerce the one un-writable provenance field; touch nothing else."""
    import anndata as ad
    if getattr(ad.AnnData.write_h5ad, "_pln_patch", False):
        return
    orig = ad.AnnData.write_h5ad

    def write_h5ad(self, *a, **k):
        g = (self.uns or {}).get("prep", {}).get("gdna", {})
        for key in ("area_table", "vendor_check"):
            if isinstance(g.get(key), (list, dict)):
                g[key] = json.dumps(g[key], default=str)
        return orig(self, *a, **k)

    write_h5ad._pln_patch = True
    ad.AnnData.write_h5ad = write_h5ad


def stage_bases(donors, overwrite=False):
    _json_safe_write()
    import build_phenotyped_base as bpb
    rows = []
    for d in donors:
        r = bpb.build_one(d, roi=ROI, out_root=bpb.DEFAULT_OUT,
                          paths_csv=bpb.DEFAULT_PATHS_CSV, overwrite=overwrite)
        print(f"  {r}", flush=True)
        rows.append(r)
    return rows


# --------------------------------------------------------------------------- stage: panel
def stage_panel(donors, python=sys.executable):
    pins = patch_pins()
    paths = [str(pins.base_path(d)) for d in donors]
    missing = [p for p in paths if not pathlib.Path(p).is_file()]
    if missing:
        raise SystemExit("bases absent; run `bases` first:\n  " + "\n  ".join(missing))
    cmd = [python, str(REPO / "analysis/xenium/scripts/finalize_panel.py"),
           "--tissue", "pancreas",
           "--audit-csv", str(REPO / "analysis/xenium/panel_audit.csv"),
           "--subpanel-dir", str(REPO / "data/development/xenium_panels/subpanels"),
           "--paths", *paths]
    print(" ".join(cmd), flush=True)
    return subprocess.run(cmd, cwd=REPO).returncode


# --------------------------------------------------------------------------- stage: type
def stage_type(donors, python=sys.executable):
    pins = patch_pins()
    rc = 0
    for d in donors:
        base = pins.base_path(d)
        if not base.is_file():
            raise SystemExit(f"donor {d}: no pLN base at {base}")
        cmd = [python, str(REPO / "analysis/xenium/scripts/hierarchical_phenotype.py"),
               "--tissue", "pancreas", "--output", "slim", "--no-fig",
               "--subpanel-dir", str(REPO / "data/development/xenium_panels/subpanels"),
               "--paths", str(base)]
        print(f"\n=== donor {d} roi={ROI}\n" + " ".join(cmd), flush=True)
        rc |= subprocess.run(cmd, cwd=REPO).returncode
    return rc


# --------------------------------------------------------------------------- stage: calls
def stage_calls(donors, judge=False):
    pins = patch_pins()
    import gate
    import nullgate_typer
    for d in donors:
        p = nullgate_typer.run(d)
        print(f"{nullgate_typer.ARM} {d} roi={ROI}: wrote {p}", flush=True)
        if judge:
            gate.judge(nullgate_typer.ARM, d)
    return 0


# --------------------------------------------------------------------------- stage: release
def stage_release(donors, out_dir):
    pins = patch_pins()
    import release as rel_mod
    rel_mod.ROI = ROI                      # module-level copy, bound at its import
    out = rel_mod.write("scanvi-per-donor", list(donors), pathlib.Path(out_dir))
    print(f"[release] wrote {out}", flush=True)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=("bases", "panel", "type", "calls", "release"))
    ap.add_argument("--donors", nargs="+", default=list(DEFAULT_DONORS))
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--judge", action="store_true")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--out", default=str(REPO / "data/development/xenium_typing/release_pln"))
    a = ap.parse_args(argv)
    print(f"pln_typing: stage={a.stage} roi={ROI} donors={a.donors}", flush=True)
    if a.stage == "bases":
        stage_bases(a.donors, overwrite=a.overwrite)
        return 0
    if a.stage == "panel":
        return stage_panel(a.donors, python=a.python)
    if a.stage == "type":
        return stage_type(a.donors, python=a.python)
    if a.stage == "calls":
        return stage_calls(a.donors, judge=a.judge)
    stage_release(a.donors, a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
