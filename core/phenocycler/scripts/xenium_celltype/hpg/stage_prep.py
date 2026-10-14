#!/usr/bin/env python
"""Prep one donor for the production run: counts -> lognorm -> PCA -> KODAMA's three inputs.

Reuses :mod:`run_donor`'s ``prepare`` unchanged, so the production path and the workstation
path compute the same thing; what this adds is the site resolution, the run-manifest
provenance, and the refusal to subsample.

NO SUBSAMPLE, AND THE REASON IS SCIENTIFIC. ``run_donor.DEFAULT_MAX_CELLS`` is 50,000, which
was the right call on a 16-core workstation and is the wrong call here. The deliverable is
not only donor-level composition; the study question is whether immune cells sit near islets,
and the peri-islet compartment is a small fraction of a section. A uniform 50k-of-1.1M sample
keeps ~4.5% of the peri-islet cells, and ``nullgate.MIN_GROUP_CELLS`` is 50, so a population
at 0.05% frequency contributes ~25 cells and cannot be annotated at all -- it is not measured
as rare, it is invisible. At full n the same population is ~550 cells and the cluster
granularity can go fine enough to isolate it. The subsample was a resource decision with a
computable cost; running all 14,893,162 cells removes the cost rather than reporting it.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))                 # xenium_celltype/
sys.path.insert(0, str(HERE.parents[2]))             # core/phenocycler

import numpy as np                                    # noqa: E402
import pandas as pd                                   # noqa: E402
import scipy.sparse as sp                             # noqa: E402

import sites                                          # noqa: E402
import xenio                                          # noqa: E402
import run_donor                                      # noqa: E402


def load_at_site(donor: str, roi: str, xpaths: Path,
                 serial: str | None = None) -> "xenio.DonorData":
    """xenio.load_donor, but with the paths this site actually has.

    xenio hard-codes ``/data/xenium``; sites.py knows the two /orange trees and which of
    them is the only copy of MM046/MM047. Rather than fork the loader, the resolved
    locations are pushed into its module globals for the duration of the call -- one place,
    named, and reversed afterwards, so nothing else in the process sees a different root.

    ``serial`` names WHICH section for a donor+roi carrying several bundles (6533, 112 at
    roi=panc). Omitting it there raises out of ``sites.resolve`` rather than picking one.
    """
    r = sites.resolve(donor, xpaths, roi=roi, serial=serial)
    if r.expression_source == "h5ad":
        import anndata as ad
        a = ad.read_h5ad(r.h5ad)
        ids = np.array(list(map(str, a.obs_names)))
        genes = list(map(str, a.var_names))
        counts = sp.csr_matrix(a.X)
    else:
        ids, genes, counts = xenio._read_bundle_matrix(r.bundle)

    cells = pd.read_parquet(r.cells_parquet, columns=["cell_id", "x_centroid", "y_centroid"])
    cells["cell_id"] = cells.cell_id.astype(str)
    pos = pd.Series(np.arange(len(cells)), index=cells.cell_id.values)
    hit = pos.reindex(ids)
    ok = hit.notna().to_numpy()
    if not ok.any():
        raise xenio.XeniumIOError(
            f"donor {donor}: no cell_id from {r.expression_source} matched cells.parquet")
    rows = hit[ok].to_numpy().astype(int)
    xy = cells.loc[:, ["x_centroid", "y_centroid"]].to_numpy()[rows]
    return xenio.DonorData(donor_id=donor, cell_id=ids[ok], genes=genes,
                           counts=counts[ok].tocsr(), xy=xy, source=r.expression_source,
                           n_dropped_no_spatial=int((~ok).sum()),
                           serial=r.serial, bundle=r.bundle,
                           h5ad_refused_reason=r.h5ad_refused_reason)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--serial", default=None,
                    help="Xenium slide serial naming WHICH section (e.g. 0080870). Required "
                         "for 6533 / 112 at roi=panc, which each carry a 2024 nPOD run and a "
                         "2026 Pheno-Xeno run; taken from the run manifest's donor entry "
                         "when it records one.")
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)

    man = json.loads((a.run / "work_manifest.json").read_text())
    n_comps = int(man["pca"]["n_comps"])
    seed = int(man["kodama"]["seed"])
    repo = HERE.parents[2]

    dest = a.run / f"donor_id={a.donor}" / f"roi={a.roi}" / "input"
    done = dest / "prep.json"
    if done.is_file() and not a.force:
        prev = json.loads(done.read_text())
        if prev.get("complete") and prev.get("n_comps") == n_comps:
            print(f"[prep] {a.donor} already complete ({prev['n_cells']:,} cells) — skipping")
            return 0

    xpaths = repo / "data" / "integration" / "xenium_paths.csv"
    expected = next(d for d in man["donors"] if d["donor_id"] == a.donor)
    # CLI wins, then whatever the manifest recorded for this donor. Neither invents one.
    serial = a.serial or expected.get("serial") or None
    dd = load_at_site(a.donor, a.roi, xpaths, serial=serial)
    if len(dd.cell_id) != expected["n_cells_expression"] - dd.n_dropped_no_spatial:
        raise SystemExit(
            f"[prep] {a.donor}: loaded {len(dd.cell_id):,} cells; the manifest recorded "
            f"{expected['n_cells_expression']:,} in the expression source and this load "
            f"dropped {dd.n_dropped_no_spatial:,} for want of coordinates. The counts do "
            f"not reconcile, so the source has changed under the manifest.")

    meta = run_donor.prepare(dd, dest, n_comps=n_comps, seed=seed)
    meta.update({"complete": True, "subsampled": False,
                 "n_available": int(len(dd.cell_id)), "n_used": int(len(dd.cell_id)),
                 "site": sites.resolve(a.donor, xpaths, roi=a.roi, serial=serial).site,
                 "serial": dd.serial, "bundle": str(dd.bundle),
                 "h5ad_refused_reason": dd.h5ad_refused_reason,
                 "run_manifest_created": man["created_at_utc"]})
    done.write_text(json.dumps(meta, indent=2))
    print(f"[prep] {a.donor}: {meta['n_cells']:,} cells from {meta['source']}; "
          f"{n_comps} PCs explain {meta['explained_variance_ratio_sum']:.1%}; "
          f"{dd.n_dropped_no_spatial:,} cells had no coordinates and were dropped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
