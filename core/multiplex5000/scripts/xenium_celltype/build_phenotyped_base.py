#!/usr/bin/env python
"""Turn a Xenium bundle into the h5ad the established typer expects.

`analysis/xenium/scripts/hierarchical_phenotype.py` is the reference Xenium typer -- two-tier
identity-core anchoring with scANVI propagation, canonical under that lane's do-not-regress
decisions #16/#17 -- and `phenocycler/integration/import_xenium.py` is already written against
its output. It has never run on this cohort, because what it consumes does not exist here: the
46 h5ads under `/data/xenium/h5ad_out/` carry raw counts and donor metadata and **nothing
else** -- no `obsm['spatial']`, no layers, no cell types. Four cohort donors have no h5ad at
all.

This builds that input. It is deliberately NOT a port of `notebooks/01_preprocessing_v2.ipynb`:
the typer needs exactly three things -- `layers['counts']`, `layers['lognorm']` and
`var['panel_for_embedding']` (the last stamped afterwards by `finalize_panel.py`) -- and it
fails loud on each. The notebook's remaining outputs (a flat scVI latent, a leiden, a UMAP)
are consumed by nothing downstream of the hierarchical rewire.

**Why `prefer="bundle"` for every donor, not just the four without an h5ad.** The dual source
is a fact about which files were transferred, not a methodological fork, and `xenio`'s own
cross-validation on donor 6539 settles that they agree: identical gene order, scores to 6e-7,
100% argmax agreement over 20,000 cells. What differs is the *population* -- the h5ad is an
upstream Seurat QC subset (317,610 of 6539's 319,005 cells). Reading the bundle everywhere
means all 20 donors start from the same cells and take the same QC recorded here, instead of
16 inheriting an undocumented upstream filter and 4 taking none. That asymmetry would sit
underneath every cross-donor comparison the study makes.

QC follows the Xenium lane's decision #10: absolute LOWER bounds (`min_genes`, `min_counts`)
because those are sanity floors for noise and segmentation artifacts, and a DENSITY-based
upper bound (`counts / cell_area` above a quantile) because in FFPE Xenium `cell_area` varies
5-10x and an absolute counts cap predominantly removes large *healthy* cells -- the large
exocrine and epithelial ones -- rather than doublets.  (Wording corrected 2026-09-20: this
read "acinar, ductal, large epithelial".  EXOCRINE IS ONE LABEL on this panel and no acinar
or ductal claim is made; the QC rationale is unchanged, only the label is.)

    cd core/multiplex5000
    ../../.venv/bin/python scripts/xenium_celltype/build_phenotyped_base.py --donors 6539
    ../../.venv/bin/python scripts/xenium_celltype/build_phenotyped_base.py --workers 8

Writes `<out>/<serial>/<roi>/<serial>_phenotyped.h5ad`. The nesting is not cosmetic:
`integration/manifest.py` only falls back to the flat `<serial>/<serial>_phenotyped.h5ad` when
a serial carries exactly one ROI, and 25 of this cohort's 26 serials carry more than one
(every pancreas slide also carries pLN regions). The flat name would never be found.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))             # core/multiplex5000, for `phenocycler`

import xenio                                          # noqa: E402

DEFAULT_OUT = Path("/data/xenium/processed")
DEFAULT_PATHS_CSV = HERE.parents[1] / "data" / "integration" / "xenium_paths.csv"

#: Decision #10. Lower bounds absolute (noise floors, not size-correlated); upper bound on
#: DENSITY, because Spearman(total_counts, cell_area) is ~0.78 in FFPE Xenium and an absolute
#: cap would drop large healthy cells rather than doublets.
MIN_GENES = 20
MIN_COUNTS = 50
DENSITY_QUANTILE = 0.98

#: The bundle directory name carries the slide serial, which is the sample id every downstream
#: consumer keys on -- `manifest.py`'s h5ad discovery, and the existing integration fixture.
_SERIAL = re.compile(r"output-XETG\d+__(\d+)__")


class BuildError(RuntimeError):
    """A donor cannot be built -- never silently skipped."""


def serial_of(source_path: str) -> str:
    m = _SERIAL.search(str(source_path))
    if not m:
        raise BuildError(f"no XETG serial in bundle path {source_path!r}")
    return m.group(1)


def bundle_extras(bundle: Path, cell_ids: np.ndarray) -> pd.DataFrame:
    """`cell_area` and friends for the cells we kept, in their order.

    `xenio.load_donor` reads only the three columns it needs. `cell_area` is read here rather
    than added there because it is a QC input, not a coordinate, and the density filter is
    this script's decision -- not something every consumer of `load_donor` should inherit.
    """
    cols = ["cell_id", "cell_area", "nucleus_area", "transcript_counts", "segmentation_method",
            # The vendor's own totals. `total_counts` INCLUDES control and genomic counts and
            # is the only record of them once the matrix is loaded gene-expression-only;
            # `sc.pp.calculate_qc_metrics` overwrites `obs['total_counts']` with a gene-only
            # row sum wherever it is called, so it is preserved here under its own name.
            "total_counts", "genomic_control_counts"]
    df = pd.read_parquet(bundle / "cells.parquet", columns=cols)
    df["cell_id"] = df.cell_id.astype(str)
    df = df.set_index("cell_id").reindex(pd.Index(cell_ids, name="cell_id"))
    if df.cell_area.isna().any():
        raise BuildError(f"{bundle}: cell_area missing for "
                         f"{int(df.cell_area.isna().sum())} of {len(df)} cells")
    return df


def qc_mask(counts: sp.csr_matrix, cell_area: np.ndarray, *,
            gdna_lambda: np.ndarray | None = None,
            n_panel_genes: int | None = None) -> tuple[np.ndarray, dict]:
    """Decision #10's filter, returning the mask AND what it removed.

    ``MIN_COUNTS`` is applied to **gDNA-corrected** counts when ``gdna_lambda`` is supplied.
    On these runs the genomic-DNA background is ~126 transcripts in a median 433-transcript
    cell, so a floor of 50 on raw counts admits cells that are majority background: measured,
    14.8% (0080870) and 18.6% (0081366) of the cells a raw floor of 50 keeps have more
    expected background than real signal. Correcting the quantity rather than raising the
    number costs 7.0% / 9.4% of cells; reaching the same place by raising the RAW floor to
    300 costs 31% / 34% and still admits them.

    The correction subtracts an EXPECTATION, not identified reads -- nothing says which counts
    are genomic. That is sound for a threshold and for aggregates, and must not be mistaken
    for a cleaned expression vector.
    """
    n_counts = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    density = np.divide(n_counts, cell_area, out=np.zeros_like(n_counts, dtype=float),
                        where=cell_area > 0)
    cap = float(np.quantile(density[np.isfinite(density)], DENSITY_QUANTILE))

    info = {"min_genes": MIN_GENES, "min_counts": MIN_COUNTS,
            "density_quantile": DENSITY_QUANTILE, "density_cap_counts_per_um2": cap,
            "dropped_min_genes": int((n_genes < MIN_GENES).sum()),
            "dropped_density": int((density > cap).sum())}

    if gdna_lambda is None or n_panel_genes is None:
        counts_for_floor = n_counts
        info["min_counts_applied_to"] = "raw"
    else:
        expected_bg = np.asarray(gdna_lambda, dtype=float) * float(n_panel_genes)
        counts_for_floor = np.maximum(n_counts - expected_bg, 0.0)
        info.update({
            "min_counts_applied_to": "gdna_corrected",
            "expected_background_median": float(np.median(expected_bg)),
            "frac_cells_majority_background": float((expected_bg > 0.5 * n_counts).mean()),
            "frac_cells_background_exceeds_counts": float((expected_bg > n_counts).mean()),
        })
    info["dropped_min_counts"] = int((counts_for_floor < MIN_COUNTS).sum())

    keep = (n_genes >= MIN_GENES) & (counts_for_floor >= MIN_COUNTS) & (density <= cap)
    return keep, info


def _object_index(values) -> pd.Index:
    """An unnamed index of plain Python strings -- see the note in `build_one`."""
    return pd.Index(np.asarray(values, dtype=object), dtype=object)


def _object_categorical(values) -> pd.Categorical:
    """A categorical whose CATEGORIES are `object`, not pandas 3's inferred StringDtype."""
    arr = np.asarray(values, dtype=object)
    return pd.Categorical(arr, categories=pd.Index(np.unique(arr), dtype=object))


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_one(donor: str, *, roi: str, out_root: Path, paths_csv: Path,
              overwrite: bool = False) -> dict:
    import anndata as ad

    xen = pd.read_csv(paths_csv, dtype=str)
    hit = xen[(xen.donor_id == donor) & (xen.roi == roi)]
    hit = hit[hit.source_path.notna() & (hit.source_path != "NA")]
    if hit.empty:
        raise BuildError(f"donor {donor} roi={roi} has no bundle path in {paths_csv}")
    source_path = str(hit.iloc[0].source_path).rstrip("/")
    serial = serial_of(source_path)
    out = out_root / serial / roi / f"{serial}_phenotyped.h5ad"
    if out.exists() and not overwrite:
        return {"donor_id": donor, "serial": serial, "status": "exists", "path": str(out)}

    dd = xenio.load_donor(donor, paths_csv, roi=roi, prefer="bundle")
    bundle = xenio.XENIUM_ROOT / source_path
    extras = bundle_extras(bundle, dd.cell_id)

    # Per-cell genomic-DNA background, from the bundle's own genomic control probes.
    # `summarise` re-derives the vendor's headline false-positive number and refuses if it
    # does not reproduce, so the correction below rests on a checked identity.
    import gdna as gdna_mod
    comp = gdna_mod.panel_composition(bundle)
    lam, lam_raw, area_tbl, shrink_a = gdna_mod.background_rates(
        extras["genomic_control_counts"].to_numpy(), extras.cell_area.to_numpy(),
        comp["genomic_control"])
    gdna_check = gdna_mod.verify_against_vendor(
        float(lam.mean()), comp["gene"], gdna_mod.vendor_metrics(bundle))

    keep, qc = qc_mask(dd.counts, extras.cell_area.to_numpy(),
                       gdna_lambda=lam, n_panel_genes=comp["gene"])
    counts = dd.counts[keep].tocsr().astype(np.float32)
    lognorm = xenio.lognorm(dd.counts[keep].tocsr()).astype(np.float32)

    # Every string here is written as `object`, not pandas 3's default `StringDtype`.
    # anndata < 0.13 refuses to write a nullable StringArray unless
    # `settings.allow_write_nullable_strings` is toggled on, and this file has to be readable
    # by the typer's env (anndata 0.12.4) as well as by this one (0.13.2).  Categorical for
    # the low-cardinality columns, which is also how they should be stored.
    # Two traps here, both silent, both caught only because the round-trip was checked.
    #
    # 1. The indices are DELIBERATELY UNNAMED.  A named index is written as an ordinary
    #    column plus an `_index` attribute, and round-tripping that through a second writer
    #    (finalize_panel.py, under a different anndata) reset BOTH indices to positional
    #    strings -- var_names became '0','1','2'...  Every gene lookup then missed: the
    #    identity panel came back 0/5101 and donor_sex 'unknown', and nothing raised.  A
    #    typing run on that file would have been meaningless and looked fine.
    #
    # 2. Every string is forced to `object`, not pandas 3's inferred `StringDtype`.
    #    `pd.Index(np.array([...], dtype=object))` returns an ArrowStringArray under pandas
    #    3.0, and anndata < 0.13 REFUSES to write one unless
    #    `settings.allow_write_nullable_strings` is toggled.  Toggling it would make the file
    #    unreadable by the very env the typer runs in (scvi-env, anndata 0.12.4), so the
    #    dtype is pinned instead of the setting.
    #
    # Both are pinned by `tests/test_build_phenotyped_base.py`.
    n_keep = int(keep.sum())
    obs = pd.DataFrame(index=_object_index(dd.cell_id[keep]))
    for c in ("cell_area", "nucleus_area", "transcript_counts",
              "total_counts", "genomic_control_counts"):
        obs[c] = extras[c].to_numpy()[keep]
    # The background every downstream detection decision has to know about. Carried on the
    # cell so no later stage has to re-open the bundle to find it.
    obs["gdna_lambda"] = lam[keep].astype(np.float32)
    obs["gdna_raw"] = lam_raw[keep].astype(np.float32)
    obs["gdna_expected_counts"] = (lam[keep] * comp["gene"]).astype(np.float32)
    obs["segmentation_method"] = _object_categorical(
        extras["segmentation_method"].to_numpy()[keep])
    obs["donor_id"] = _object_categorical([donor] * n_keep)
    obs["roi"] = _object_categorical([roi] * n_keep)
    obs["sample"] = _object_categorical([serial] * n_keep)

    adata = ad.AnnData(
        X=lognorm, obs=obs,
        var=pd.DataFrame(index=_object_index(dd.genes)))
    adata.layers["counts"] = counts
    adata.layers["lognorm"] = lognorm
    adata.obsm["spatial"] = dd.xy[keep].astype(np.float32)

    # Provenance. This is NOT a StageManifest and must not pretend to be one -- nothing here
    # writes under data/runs/. But it must be attributable, because it is the input every
    # Xenium cell type in this cohort descends from.
    adata.uns["prep"] = {
        "built_by": "core/multiplex5000/scripts/xenium_celltype/build_phenotyped_base.py",
        "source": dd.source,
        "bundle": source_path,
        "xenium_paths_csv_sha256": _sha256(paths_csv),
        "donor_id": donor, "roi": roi, "serial": serial,
        "n_cells_in": int(dd.counts.shape[0]),
        "n_cells_out": int(keep.sum()),
        "n_dropped_no_spatial": int(dd.n_dropped_no_spatial),
        "gdna": {"n_genomic_probes": int(comp["genomic_control"]),
                 "n_panel_genes": int(comp["gene"]),
                 "shrink_a": float(shrink_a),
                 "lambda_mean": float(lam.mean()),
                 "lambda_raw_mean": float(lam_raw.mean()),
                 "vendor_check": gdna_check,
                 "area_table": area_tbl.to_dict("records")},
        "n_genes": int(len(dd.genes)),
        "qc": qc,
        "lognorm": "median-depth normalize_total then log1p (xenio.lognorm)",
    }

    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="gzip")
    return {"donor_id": donor, "serial": serial, "status": "written", "path": str(out),
            "n_cells_in": int(dd.counts.shape[0]), "n_cells_out": int(keep.sum()), **qc}


def xenium_eligible(donors) -> list[str]:
    """The donors, or a refusal naming any whose Xenium is centrally excluded.

    Asked of ``multiplex5000.cohort.excluded_on(donor, "xenium")``: 6457's exclusion covers its
    PhenoCycler data only (maintainer, 2026-10-22), so its Xenium bundle is built here.
    """
    from multiplex5000.cohort import DONOR_EXCLUSIONS, donor_of, excluded_on
    donors = [str(d) for d in donors]
    blocked = [d for d in donors if excluded_on(d, "xenium")]
    if blocked:
        details = "; ".join(f"{d}: {DONOR_EXCLUSIONS[donor_of(d)]}" for d in blocked)
        raise ValueError(f"xenium phenotyped base includes excluded donor(s) on xenium: {details}")
    return donors


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None,
                    help="default: the whole eligible cohort")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--xenium-paths", type=Path, default=DEFAULT_PATHS_CSV)
    ap.add_argument("--workers", type=int, default=1,
                    help="gzip-write bound, so parallelism helps; each worker peaks ~3 GB")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)

    if args.donors:
        donors = list(args.donors)
    else:
        import json as _json
        spec = _json.loads((HERE.parents[1] / "cohort-spec.json").read_text())
        donors = list(spec["expected_donors"])
    donors = xenium_eligible(donors)

    print(f"building {len(donors)} donors  roi={args.roi}  -> {args.out}", flush=True)
    rows = []
    if args.workers <= 1:
        for d in donors:
            rows.append(build_one(d, roi=args.roi, out_root=args.out,
                                  paths_csv=args.xenium_paths, overwrite=args.overwrite))
            print(f"  {rows[-1]}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            fut = {ex.submit(build_one, d, roi=args.roi, out_root=args.out,
                             paths_csv=args.xenium_paths, overwrite=args.overwrite): d
                   for d in donors}
            for f in as_completed(fut):
                rows.append(f.result())
                print(f"  {rows[-1]}", flush=True)

    summary = args.out / f"build_summary_{args.roi}.json"
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
